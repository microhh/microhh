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
#include <functional>
#include <iostream>

#include "tools.h"
#include "cuda_buffer.h"
#include "cuda_tiling.h"

#ifdef ENABLE_KERNEL_LAUNCHER
#include "kernel_launcher.h"
#include "kernel_launcher/registry.h"

struct GridKernel: kernel_launcher::IKernelDescriptor {
    GridKernel(
        GridFunctor meta,
        kernel_launcher::TypeInfo functor_type,
        std::vector<kernel_launcher::TypeInfo> param_types,
        Grid_layout grid);
    kernel_launcher::KernelBuilder build() const override;
    bool equals(const IKernelDescriptor& that) const override;
    size_t hash() const override;

    GridFunctor meta;
    kernel_launcher::TypeInfo functor_type;
    std::vector<kernel_launcher::TypeInfo> param_types;
    Grid_layout grid;
};

namespace kernel_launcher {
    template <typename T>
    struct IntoKernelArg<::cuda_span<T>> {
        static KernelArg convert(::cuda_span<T> s) {
            return KernelArg::from_array(s.data(), s.size());
        }
    };
}

bool launch_kernel(
        cudaStream_t stream,
        GridKernel kernel,
        const std::vector<kernel_launcher::KernelArg>& args);
#endif

template <typename T>
struct convert_kernel_arg {
    using type = T;

    static T call(T input) {
        return std::move(input);
    }
};

template <typename T>
struct convert_kernel_arg<T&>: convert_kernel_arg<T> {};

template <typename T>
struct convert_kernel_arg<T&&>: convert_kernel_arg<T> {};

template <typename T>
struct convert_kernel_arg<cuda_span<T>> {
    using type = T*;

    static cuda_span<T> call(cuda_span<T> input) {
        return {input.data(), input.size()};
    }
};

template <typename T>
struct convert_kernel_arg<const cuda_span<T>>: convert_kernel_arg<cuda_span<T>> {};

template <typename T>
struct convert_kernel_arg<cuda_vector<T>>: convert_kernel_arg<cuda_span<const T>> {};

template <typename T>
struct convert_kernel_arg<const cuda_vector<T>>: convert_kernel_arg<cuda_span<const T>> {};

template <typename F, typename... Args>
void launch_grid_kernel(
        const Grid_layout &gd,
        Args&&... args
)
{
    GridFunctor meta = F::meta;

#ifdef ENABLE_KERNEL_LAUNCHER
    std::vector<kernel_launcher::TypeInfo> param_types = {
            kernel_launcher::TypeInfo::of<typename convert_kernel_arg<Args>::type>()...
    };

    std::vector<kernel_launcher::KernelArg> kernel_args = {
            kernel_launcher::into_kernel_arg(convert_kernel_arg<Args>::call(args))...
    };

    GridKernel kernel(
        meta,
        kernel_launcher::TypeInfo::of<F>(),
        std::move(param_types),
        gd
    );

    bool success = launch_kernel(
            nullptr,
            std::move(kernel),
            kernel_args);

    // If the kernel was launched successfully, exit the function now.
    if (success) {
        return;
    }
#endif

    // Fallback function. This kernel is called either if kernel_launcher is disabled or if `launch_kernel` failed.
    dim3 problem_size = {
            uint(gd.iend - gd.istart),
            uint(gd.jend - gd.jstart),
            uint(gd.kend - gd.kstart)};

    dim3 block_size = meta.block_size;
    dim3 grid_size = {
            (problem_size.x / block_size.x) + bool(problem_size.x % block_size.x != 0),
            (problem_size.y / block_size.y) + bool(problem_size.y % block_size.y != 0),
            (problem_size.z / block_size.z) + bool(problem_size.z % block_size.z != 0)
    };

    grid_tiling_kernel<F><<<grid_size, block_size>>>(
            gd,
            typename convert_kernel_arg<Args>::type(args)...);
    cuda_check_error();
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

template <typename F, typename TF, typename... Args>
void launch_grid_kernel(
        const Grid_layout &grid_layout,
        Args&&... args
)
{
    launch_grid_kernel<F>(
            grid_layout,
            std::forward<Args>(args)...
    );
}

#endif //CUDA_COMPILER_H
