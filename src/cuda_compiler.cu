#include <vector>
#include <string>
#include <unistd.h>

#ifdef ENABLE_KERNEL_LAUNCHER
#include "kernel_launcher.h"
#include "cuda_compiler.h"

namespace kl = kernel_launcher;

const std::string& home_directory();

GridKernel::GridKernel(
        GridFunctor meta,
        kernel_launcher::TypeInfo functor_type,
        std::vector<kl::TypeInfo> param_types,
        Grid_layout grid):
    meta(std::move(meta)),
    functor_type(functor_type),
    param_types(std::move(param_types)),
    grid(grid) {}

std::string GridKernel::tuning_key() const {
    std::string name = meta.name;

    for (auto p: param_types) {
        auto base = p.remove_pointer().remove_const();

        if (base == kl::type_of<float>()) {
            return name + "_float";
        }

        if (base == kl::type_of<double>()) {
            return name + "_double";
        }
    }

    return name;
}

bool GridKernel::equals(const KernelDescriptor& that) const {
    if (const GridKernel* g = dynamic_cast<const GridKernel*>(&that)) {
        return g->meta.name == meta.name &&
                g->functor_type == functor_type &&
                g->param_types == param_types &&
                g->grid == grid;
    }

    return false;
}

size_t GridKernel::hash() const {
    return kl::hash_fields(meta.name, functor_type);
}

kl::KernelBuilder GridKernel::build() const {
    std::string source = "#include \"" + std::string(meta.file) + "\"\n";
    source += R"(
        #include "cuda_tiling.h"

        template <typename F, typename... Args>
        __global__
        __launch_bounds__(BLOCK_SIZE_X * BLOCK_SIZE_Y * BLOCK_SIZE_Z, BLOCKS_PER_SM)
        void kernel(Args... args) {
            dim3 num_blocks = {NUM_BLOCKS_X, NUM_BLOCKS_Y, NUM_BLOCKS_Z};
            dim3 block_index = unravel_dim3(blockIdx.x, num_blocks, AXES_PERMUTATION);

            Grid_layout gd = {
                GRID_START_I,
                GRID_START_J,
                GRID_START_K,
                GRID_END_I,
                GRID_END_J,
                GRID_END_K,
                GRID_STRIDE_I,
                GRID_STRIDE_J,
                GRID_STRIDE_K
            };

            using Tiling = StaticTilingStrategy<
                BLOCK_SIZE_X,
                BLOCK_SIZE_Y,
                BLOCK_SIZE_Z,
                TILE_FACTOR_X,
                TILE_FACTOR_Y,
                TILE_FACTOR_Z,
                UNROLL_FACTOR_X,
                UNROLL_FACTOR_Y,
                UNROLL_FACTOR_Z,
                TILE_CONTIGUOUS_X,
                TILE_CONTIGUOUS_Y,
                TILE_CONTIGUOUS_Z,
            >;

            cta_execute_tiling_with_edges<EDGE_LEVELS, Tiling>(gd, block_index, F{}, args...);
        }
    )";

    std::vector<uint32_t> block_size_x = {16, 32, 64, 128, 256};
    std::vector<uint32_t> block_size_y = {1, 2, 4, 8, 16};
    std::vector<uint32_t> block_size_z = {1, 2, 4, 8, 16};

    kl::KernelBuilder builder("kernel", kl::KernelSource("kernel.cu", source));
    auto bx = builder.tune("BLOCK_SIZE_X", block_size_x, meta.block_size.x);
    auto by = builder.tune("BLOCK_SIZE_Y", block_size_y, meta.block_size.y);
    auto bz = builder.tune("BLOCK_SIZE_Z", block_size_z, meta.block_size.z);
    auto blocks_per_sm = builder.tune_define("BLOCKS_PER_SM", {1, 2, 3, 4, 5, 6});

    auto tx = builder.tune_define("TILE_FACTOR_X", {1, 2, 3, 4});
    auto ty = builder.tune_define("TILE_FACTOR_Y", {1, 2, 3, 4});
    auto tz = builder.tune_define("TILE_FACTOR_Z", {1, 2, 3, 4, 8, 16, 32});

    // How many loops to unroll
    // - 0: no unroll
    // - 1: only inner loop
    // - 2: two inner loops
    // - 3: all loops
    auto unroll_depth = builder.tune_define("LOOP_UNROLL_DEPTH", {3, 2, 1, 0});

    // What order to unravel the block index.
    builder.tune_define("AXES_PERMUTATION", {0, 1, 2, 3, 4, 5});

    // Tiling is contiguous or block strided
    auto tile_cont = builder.tune_define("TILE_CONTIGUOUS_XY", {0, 1});

    // Number of thread blocks
    auto nx = kl::div_ceil(grid.iend - grid.istart, bx * tx);
    auto ny = kl::div_ceil(grid.jend - grid.jstart, by * ty);
    auto nz = kl::div_ceil(grid.kend - grid.kstart, bz * tz);

    builder
        .block_size(bx, by, bz)
        .grid_size(nx * ny * nz);

    builder
        .define(bx)
        .define(by)
        .define(bz)
        .define("GRID_START_I", std::to_string(grid.istart))
        .define("GRID_START_J", std::to_string(grid.jstart))
        .define("GRID_START_K", std::to_string(grid.kstart))
        .define("GRID_END_I", std::to_string(grid.iend))
        .define("GRID_END_J", std::to_string(grid.jend))
        .define("GRID_END_K", std::to_string(grid.kend))
        .define("GRID_STRIDE_I", std::to_string(grid.ii))
        .define("GRID_STRIDE_J", std::to_string(grid.jj))
        .define("GRID_STRIDE_K", std::to_string(grid.kk))
        .define("NUM_BLOCKS_X", nx)
        .define("NUM_BLOCKS_Y", ny)
        .define("NUM_BLOCKS_Z", nz)
        .define("UNROLL_FACTOR_X", kl::ifelse(unroll_depth >= 1, tx, 1))
        .define("UNROLL_FACTOR_Y", kl::ifelse(unroll_depth >= 2, ty, 1))
        .define("UNROLL_FACTOR_Z", kl::ifelse(unroll_depth >= 3, tz, 1))
        .define("TILE_CONTIGUOUS_X", "0")
        .define("TILE_CONTIGUOUS_Y", tile_cont)
        .define("TILE_CONTIGUOUS_Z", tile_cont)
        .define("EDGE_LEVELS", std::to_string(meta.edge_levels));

    builder.compiler_flags(
            "--restrict",
            "-std=c++17",
            "-I" + home_directory() + "/include");

    // restrictions:
    // - Threads per block should not too small (>=64) or too big (<=1024)
    // - Threads per SM should not be too big (<=2048)
    // - Items per thread should not be too many (<= 32)
    auto threads_per_block = bx * by * bz;
    builder.restriction(threads_per_block >= 64 && threads_per_block <= 1024);
    builder.restriction(threads_per_block * blocks_per_sm <= 2048);
    builder.restriction(tx * ty * tz <= 32);

    builder.template_arg(functor_type);
    for (const auto& p: param_types) {
        builder.template_arg(p);
    }

    return builder;
}

std::string find_home_directory() {
    std::string result;
    const char *env = getenv("MICROHH_HOME");

    // Check environment KEY
    if (env) {
        result = env;
    }

    // No success, try from __FILE__
    if (result.empty()) {
        std::string file = __FILE__;
        size_t index = file.rfind("/src/");

        if (index != std::string::npos) {
            result = file.substr(0, index);
        }
    }

    // No success, try "cwd"?
    if (result.empty()) {
        char buffer[1024] = {0};

       if (getcwd(buffer, sizeof(buffer)) != nullptr) {
           result = buffer;
       }
    }

    // No success, try "."
    if (result.empty()) {
        result = ".";
    }

    // Add trailing slash
    if (result.back() != '/') {
        result += "/";
    }

    if (env == nullptr) {
        std::cerr << "WARNING: environment variable MICROHH_HOME is not set, best guess: " << result << "\n";
    }

    return result;
}

const std::string& home_directory() {
    static std::string dir = find_home_directory();
    return dir;
}

bool launch_kernel(
        cudaStream_t stream,
        dim3 problem_size,
        GridKernel grid_kernel,
        const std::vector<kl::KernelArg>& args
) {
    static bool initialized = false;
    static bool error_state = false;

    // If an exception was thrown at some point, just return false immediately.
    if (error_state) {
        return false;
    }

    kl::AnyKernelDescriptor kernel = std::move(grid_kernel);

    try {
        if (!initialized) {
            initialized = true;
            kernel_launcher::set_global_wisdom_directory(home_directory() + "wisdom");
            kernel_launcher::set_global_tuning_directory(home_directory() + "tuning");
        }

        kernel_launcher::default_registry()
                .lookup(kernel)
                .launch(stream, problem_size, args);
        return true;
    } catch (const std::exception &e) {
        kl::log_warning() << "error occurred while compiling the following kernel: " <<
                kernel.descriptor().tuning_key() << ":" << "\n"  << e.what() << std::endl;
        kl::log_warning() << "CUDA dynamic kernel compilation is now disabled and the application is in FALLBACK mode. "
                             "No more new kernels will be compiled using Kernel Launcher." << std::endl;

        error_state = true;
        return false;
    }
}
#endif