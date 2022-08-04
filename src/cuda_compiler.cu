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
        std::string kernel_name,
        std::string kernel_file,
        kernel_launcher::TypeInfo functor_type,
        std::vector<kl::TypeInfo> param_types,
        Grid_layout grid):
    kernel_name(std::move(kernel_name)),
    kernel_file(std::move(kernel_file)),
    functor_type(functor_type),
    param_types(std::move(param_types)),
    grid(grid) {}

std::string GridKernel::tuning_key() const {
    for (auto p: param_types) {
        auto base = p.remove_pointer().remove_const();

        if (base == kl::type_of<float>()) {
            return kernel_name + "_float";
        }

        if (base == kl::type_of<double>()) {
            return kernel_name + "_double";
        }
    }

    return kernel_name;
}

bool GridKernel::equals(const KernelDescriptor& that) const {
    if (const GridKernel* g = dynamic_cast<const GridKernel*>(&that)) {
        return g->kernel_name == kernel_name &&
                g->functor_type == functor_type &&
                g->param_types == param_types &&
                g->grid == grid;
    }

    return false;
}

size_t GridKernel::hash() const {
    return kl::hash_fields(kernel_name, functor_type);
}

kl::KernelBuilder GridKernel::build() const {
    std::string source = R"(
        #include "cuda_tiling.h"
        #include "advec_2i5_kernels.cuh"

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

            cta_execute_tiling_border<3, Tiling>(gd, F{}, args...);
        }
    )";

    kl::KernelBuilder builder("kernel", kl::KernelSource("kernel.cu", source));
    auto bx = builder.tune_define("BLOCK_SIZE_X", {32, 64, 128, 256});
    auto by = builder.tune_define("BLOCK_SIZE_Y", {2, 1, 4, 8});
    auto bz = builder.tune_define("BLOCK_SIZE_Z", {2, 1, 4, 8});

    auto tx = builder.tune_define("TILE_FACTOR_X", {1, 2, 4});
    auto ty = builder.tune_define("TILE_FACTOR_Y", {1, 2, 4});
    auto tz = builder.tune_define("TILE_FACTOR_Z", {1, 2, 4, 8});

    builder.tune_define("BLOCKS_PER_SM", {1, 2, 3, 4, 5, 6});

    builder
        .block_size(bx, by, bz)
        .grid_divisors(bx * tx, by * ty, bz * tz);

    builder.define("GRID_START_I", std::to_string(grid.istart));
    builder.define("GRID_START_J", std::to_string(grid.jstart));
    builder.define("GRID_START_K", std::to_string(grid.kstart));
    builder.define("GRID_END_I", std::to_string(grid.iend));
    builder.define("GRID_END_J", std::to_string(grid.jend));
    builder.define("GRID_END_K", std::to_string(grid.kend));
    builder.define("GRID_STRIDE_I", std::to_string(grid.ii));
    builder.define("GRID_STRIDE_J", std::to_string(grid.jj));
    builder.define("GRID_STRIDE_K", std::to_string(grid.kk));
    builder.define("UNROLL_FACTOR_X", tx);
    builder.define("UNROLL_FACTOR_Y", ty);
    builder.define("UNROLL_FACTOR_Z", tz);
    builder.define("TILE_CONTIGUOUS_X", "0");
    builder.define("TILE_CONTIGUOUS_Y", "1");
    builder.define("TILE_CONTIGUOUS_Z", "1");

    builder.compiler_flag(std::string("--restrict"));
    builder.compiler_flag(std::string("--std=c++17"));
    builder.compiler_flag("-I" + home_directory() + "/include");

    auto bxyz = bx * by * bz;
    builder.restriction(bxyz >= 64 && bxyz <= 1024);

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