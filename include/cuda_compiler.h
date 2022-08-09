/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CUDA_COMPILER_H
#define CUDA_COMPILER_H

#include <type_traits>

#include "kernel_launcher.h"
#include "kernel_launcher/registry.h"
#include "cuda_buffer.h"
#include "cuda_tiling.h"

struct GridKernel: kernel_launcher::KernelDescriptor {
    GridKernel(
        GridFunctor meta,
        kernel_launcher::TypeInfo functor_type,
        std::vector<kernel_launcher::TypeInfo> param_types,
        Grid_layout grid);
    std::string tuning_key() const override;
    kernel_launcher::KernelBuilder build() const override;
    bool equals(const KernelDescriptor& that) const override;
    size_t hash() const override;

private:
    GridFunctor meta;
    kernel_launcher::TypeInfo functor_type;
    std::vector<kernel_launcher::TypeInfo> param_types;
    Grid_layout grid;
};

void launch_kernel(
        cudaStream_t stream,
        dim3 problem_size,
        GridKernel kernel,
        const std::vector<kernel_launcher::KernelArg>& args);

template <typename T>
struct KernelArgConvert {
    using type = T;

    static T call(T input) {
        return std::move(input);
    }
};

template <typename T>
struct KernelArgConvert<T&>: KernelArgConvert<T> {};

template <typename T>
struct KernelArgConvert<T&&>: KernelArgConvert<T> {};

template <typename T>
struct KernelArgConvert<cuda_span<T>> {
    using type = T*;

    static kernel_launcher::CudaSpan<T> call(cuda_span<T> input) {
        return {input.data(), input.size()};
    }
};

template <typename T>
struct KernelArgConvert<const cuda_span<T>>: KernelArgConvert<cuda_span<T>> {};

template <typename T>
struct KernelArgConvert<cuda_vector<T>>: KernelArgConvert<cuda_span<const T>> {};

template <typename T>
struct KernelArgConvert<const cuda_vector<T>>: KernelArgConvert<cuda_span<const T>> {};

template <typename F, typename... Args>
void launch_grid_kernel(
        const Grid_layout &gd,
        Args&&... args
)
{
    GridFunctor meta = F::meta;
    dim3 problem_size = {
            uint(gd.iend - gd.istart),
            uint(gd.jend - gd.jstart),
            uint(gd.kend - gd.kstart)};

    auto fallback = [&] {
        dim3 block_size = meta.block_size;
        dim3 grid_size = {
                (problem_size.x / block_size.x) + bool(problem_size.x % block_size.x != 0),
                (problem_size.y / block_size.y) + bool(problem_size.y % block_size.y != 0),
                (problem_size.z / block_size.z) + bool(problem_size.z % block_size.z != 0)
        };

        grid_tiling_kernel<F><<<grid_size, block_size>>>(
                gd,
                typename KernelArgConvert<Args>::type(args)...);
    };

    std::vector<kernel_launcher::TypeInfo> param_types = {
            kernel_launcher::TypeInfo::of<typename KernelArgConvert<Args>::type>()...
    };

    std::vector<kernel_launcher::KernelArg> kernel_args = {
            kernel_launcher::into_kernel_arg(KernelArgConvert<Args>::call(args))...
    };

    GridKernel kernel(
        meta,
        kernel_launcher::TypeInfo::of<F>(),
        std::move(param_types),
        gd
    );

    try {
        launch_kernel(
                nullptr,
                problem_size,
                std::move(kernel),
                kernel_args);
    } catch (const std::exception& e) {
        std::cerr << "WARNING: error occurred while compiling kernel: " << e.what() << "\n";
        fallback();
    }
}

template <typename F, typename TF, typename... Args>
void launch_grid_kernel(
        const Grid_data<TF> &gd,
        Args&&... args
)
{
    launch_grid_kernel<F>(
            Grid_layout::from_grid_data(gd),
            std::forward<Args>(args)...
    );
}

template <typename F, typename TF, typename... Args>
void launch_grid_kernel(
        const Grid<TF> &grid,
        Args&&... args
)
{
    launch_grid_kernel<F>(
            Grid_layout::from_grid_data(grid.get_grid_data()),
            std::forward<Args>(args)...
    );
}

#endif //CUDA_COMPILER_H
