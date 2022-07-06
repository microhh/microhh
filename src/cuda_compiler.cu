#include <vector>
#include <string>
#include <utility>
#include <sstream>
#include <iostream>

#include "cuda_tiling.h"
#include "kernel_launcher.h"
#include "cuda_compiler.h"

namespace kl = kernel_launcher;

std::string GridKernelDescriptor::tuning_key() const {
    dim3 problem_size(gd.iend - gd.istart, gd.jend - gd.jstart, gd.kend - gd.kstart);

    std::stringstream ss;
    ss << kernel_info.name() << "_" << problem_size.x << "x" << problem_size.y << "x" << problem_size.z;
    return ss.str();
}

std::string GridKernelDescriptor::unique_key() const {
    std::stringstream ss;
    ss << "grid" << kernel_info.name() << "|"
        << gd.istart << "|" << gd.jstart << "|" << gd.kstart << "|"
        << gd.iend << "|" << gd.jend << "|" << gd.kend << "|"
        << gd.ii << "|" << gd.jj << "|" << gd.kk;
    return ss.str();
}

kl::KernelBuilder GridKernelDescriptor::build_kernel(std::vector<kl::TypeInfo> param_types) const {
    std::string source = R"(
        #include "/var/scratch/sheldens/microhh/include/cuda_tiling.h"
        #include "/var/scratch/sheldens/microhh/include/advec_2i5_kernels.cuh"

        template <typename F, typename... Args>
        __global__
        __launch_bounds__(BLOCK_SIZE_X * BLOCK_SIZE_Y * BLOCK_SIZE_Z, BLOCKS_PER_SM)
        void custom_kernel_wrapper(Args... args) {
            using CustomTilingStrategy = StaticTilingStrategy<
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
                BLOCKS_PER_SM
            >;

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

            TilingKernel<F, CustomTilingStrategy>{}(gd, args...);
        }
    )";

    kl::KernelBuilder builder("custom_kernel_wrapper", kl::KernelSource("kernel.cu", source));
    builder.template_arg(kernel_info);
    for (const auto& ty: param_types) {
        builder.template_arg(ty);
    }

    builder.define("GRID_START_I", gd.istart);
    builder.define("GRID_START_J", gd.jstart);
    builder.define("GRID_START_K", gd.kstart);
    builder.define("GRID_END_I", gd.iend);
    builder.define("GRID_END_J", gd.jend);
    builder.define("GRID_END_K", gd.kend);
    builder.define("GRID_STRIDE_I", gd.ii);
    builder.define("GRID_STRIDE_J", gd.jj);
    builder.define("GRID_STRIDE_K", gd.kk);

    auto bx = builder.tune_define("BLOCK_SIZE_X", {256, 32, 64, 128});
    auto by = builder.tune_define("BLOCK_SIZE_Y", {1, 2, 4, 8, 16, 32});
    auto bz = builder.tune_define("BLOCK_SIZE_Z", {1, 2, 4, 8, 16, 32});

    auto tx = builder.tune_define("TILE_FACTOR_X", {1, 2, 4});
    auto ty = builder.tune_define("TILE_FACTOR_Y", {1, 2, 4});
    auto tz = builder.tune_define("TILE_FACTOR_Z", {1, 2});

    auto ux = builder.tune_define("UNROLL_FACTOR_X", {1, 2, 4});
    auto uy = builder.tune_define("UNROLL_FACTOR_Y", {1, 2, 4});
    auto uz = builder.tune_define("UNROLL_FACTOR_Z", {1, 2});

    builder.define("TILE_CONTIGUOUS_X", 0);
    builder.define("TILE_CONTIGUOUS_Y", 1);
    builder.define("TILE_CONTIGUOUS_Z", 1);

    builder.tune_define("BLOCKS_PER_SM", {1, 2, 3, 4, 5, 6});

    builder.compiler_flag(std::string("--restrict"));
    builder.compiler_flag(std::string("--std=c++17"));

    builder.restriction(bx * by * bz <= 1024);
    builder.restriction(bx * by * bz >= 64);

    builder.restriction(ux <= tx);
    builder.restriction(uy <= ty);
    builder.restriction(uz <= tz);

    return builder;
}

std::string append_parameter_types_to_key(std::string key, const std::vector<const kl::KernelArg*>& args)
{
    for (const auto& arg: args)
    {
        key += "|";
        key += arg->type_info().name();
    }

    return key;
}

struct KernelCache {
    KernelCache() {
        wisdom_dir_ = "/var/scratch/sheldens/microhh/wisdom";
        compiler_ = kl::NvrtcCompiler();
    }

    const kl::KernelInstance* find(
            dim3 problem_size,
            const KernelDescriptor& kernel,
            const std::vector<const kl::KernelArg*>& args
    ) {
        std::string cache_key = append_parameter_types_to_key(kernel.unique_key(), args);

        auto it = cache_.find(cache_key);
        if (it != cache_.end()) {
            return it->second.get();
        }

        std::vector<kl::TypeInfo> param_types;
        for (const auto& arg: args) {
            param_types.push_back(arg->type_info());
        }

        kl::KernelBuilder builder = kernel.build_kernel(param_types);
        const std::string tuning_key = kernel.tuning_key();

        kl::Config config;
        kl::WisdomResult result = kl::read_wisdom_file(
                tuning_key,
                builder,
                wisdom_dir_,
                config);

        if (result != kl::WisdomResult::Success) {
            if (result == kl::WisdomResult::NotFound) {
                kl::write_wisdom_file(
                        tuning_key,
                        builder,
                        wisdom_dir_,
                        wisdom_dir_,
                        problem_size,
                        args
                );
            }

            config = builder.default_config();
        }

        std::unique_ptr<kl::KernelInstance> value;

        try {
            kl::KernelInstance instance = builder.compile(config, param_types, compiler_);
            value = std::make_unique<kl::KernelInstance>(std::move(instance));
        } catch (const std::exception& e) {
            std::cerr << "kernel compilation failed: " << e.what() << std::endl;
        }

        auto p = cache_.emplace(std::move(cache_key), std::move(value));
        return p.first->second.get();
    }

private:
    std::string wisdom_dir_;
    kl::Compiler compiler_;
    std::unordered_map<std::string, std::unique_ptr<kl::KernelInstance>> cache_;
};

void launch_kernel(
        cudaStream_t stream,
        dim3 problem_size,
        const KernelDescriptor& kernel,
        const std::vector<const kl::KernelArg*>& args,
        const std::function<void()>& fallback
) {
    static KernelCache cache;
    const kl::KernelInstance* instance = cache.find(problem_size, kernel, args);

    if (instance) {
        std::vector<const void*> raw_args(args.size());
        for (size_t i = 0; i < args.size(); i++) {
            raw_args[i] = args[i]->as_ptr();
        }

        instance->launch(stream, problem_size, const_cast<void**>(raw_args.data()));
    } else {
        fallback();
    }
}