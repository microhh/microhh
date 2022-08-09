#include <vector>
#include <string>

#include "kernel_launcher.h"
#include "cuda_compiler.h"

namespace kl = kernel_launcher;

std::string& home_directory() {
    static std::string dir = "";
    if (dir.empty()) {
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

        // No success, try "."?
        if (result.empty()) {
            result = ".";
        }

        if (result.back() != '/') {
            result += "/";
        }

        if (env == nullptr) {
            std::cerr << "WARNING: environment variable MICROHH_HOME is not set, best guess: " << result << "\n";
        }

        dir = result;
    }

    return dir;
}

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
                TILE_CONTIGUOUS_Z
            >;

            cta_execute_tiling_border<BORDER_SIZE, Tiling>(gd, F{}, args...);
        }
    )";

    std::vector<uint32_t> block_size_x = {32, 64, 128, 256};
    std::vector<uint32_t> block_size_y = {1, 2, 4, 8};
    std::vector<uint32_t> block_size_z = {1, 2, 4, 8};

    kl::KernelBuilder builder("kernel", kl::KernelSource("kernel.cu", source));
    auto bx = builder.tune("BLOCK_SIZE_X", block_size_x, meta.block_size.x);
    auto by = builder.tune("BLOCK_SIZE_Y", block_size_y, meta.block_size.y);
    auto bz = builder.tune("BLOCK_SIZE_Z", block_size_z, meta.block_size.z);

    auto tx = builder.tune_define("TILE_FACTOR_X", {1, 2, 4});
    auto ty = builder.tune_define("TILE_FACTOR_Y", {1, 2, 4});
    auto tz = builder.tune_define("TILE_FACTOR_Z", {1, 2, 4, 8});

    auto blocks_per_sm = builder.tune_define("BLOCKS_PER_SM", {1, 2, 3, 4, 5, 6});

    builder
        .block_size(bx, by, bz)
        .grid_divisors(bx * tx, by * ty, bz * tz);

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
        .define("UNROLL_FACTOR_X", tx)
        .define("UNROLL_FACTOR_Y", ty)
        .define("UNROLL_FACTOR_Z", tz)
        .define("TILE_CONTIGUOUS_X", "0")
        .define("TILE_CONTIGUOUS_Y", "1")
        .define("TILE_CONTIGUOUS_Z", "1")
        .define("BORDER_SIZE", std::to_string(meta.border_size));

    builder.compiler_flags(
            "--restrict",
            "-std=c++17",
            "-I" + home_directory() + "/include"
    );

    auto threads_per_block = bx * by * bz;
    builder.restriction(threads_per_block >= 64 && threads_per_block <= 1024);
    builder.restriction(threads_per_block * blocks_per_sm <= 2048);

    builder.template_arg(functor_type);
    for (const auto& p: param_types) {
        builder.template_arg(p);
    }

    return builder;
}

void launch_kernel(
        cudaStream_t stream,
        dim3 problem_size,
        GridKernel kernel,
        const std::vector<kl::KernelArg>& args
) {
    static bool initialized = false;
    if (!initialized) {
        initialized = true;
        kernel_launcher::set_global_wisdom_directory(home_directory() + "wisdom");
        kernel_launcher::set_global_tuning_directory(home_directory() + "tuning");
    }

    kernel_launcher::default_registry()
        .lookup(std::move(kernel))
        .launch(stream, problem_size, args);
}