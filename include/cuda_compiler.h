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
#include "cuda_buffer.h"
#include "cuda_tiling.h"

template <typename T, typename Enabled = void>
struct KernelArgConvert;

template <typename T>
struct KernelArgConvert<T&>: KernelArgConvert<T> {};

template <typename T>
struct KernelArgConvert<T&&>: KernelArgConvert<T> {};

template <typename T>
struct KernelArgConvert<T, typename std::enable_if<
        std::is_trivial<T>::value &&
        !std::is_pointer<T>::value &&
        !std::is_reference<T>::value>::type>: kernel_launcher::KernelArgScalar<T> {
    using type = T;
    using base_type = kernel_launcher::KernelArgScalar<T>;

    KernelArgConvert(T v): base_type(v) {}
};

template <typename T>
struct KernelArgConvert<cuda_span<T>>: kernel_launcher::KernelArgArray<T> {
    using type = T*;
    using base_type = kernel_launcher::KernelArgArray<T>;
    KernelArgConvert(cuda_span<T> v): base_type(v.data(), v.size()) {}
};

template <typename T>
struct KernelArgConvert<const cuda_span<T>>: kernel_launcher::KernelArgArray<T> {
    using type = T*;
    using base_type = kernel_launcher::KernelArgArray<T>;
    KernelArgConvert(cuda_span<T> v): base_type(v.data(), v.size()) {}
};

template <typename T>
struct KernelArgConvert<cuda_vector<T>>: kernel_launcher::KernelArgArray<T> {
    using type = T*;
    using base_type = kernel_launcher::KernelArgArray<T>;
    KernelArgConvert(cuda_span<T> v): base_type(v.data(), v.size()) {}
};

template <typename T>
struct KernelArgConvert<const cuda_vector<T>>: kernel_launcher::KernelArgArray<const T> {
    using type = const T*;
    using base_type = kernel_launcher::KernelArgArray<const T>;
    KernelArgConvert(cuda_span<const T> v): base_type(v.data(), v.size()) {}
};

struct KernelDescriptor {
    virtual std::string tuning_key() const = 0;
    virtual std::string unique_key() const = 0;
    virtual kernel_launcher::KernelBuilder build_kernel(std::vector<kernel_launcher::TypeInfo> param_types) const = 0;
};

struct GridKernelDescriptor: KernelDescriptor {
    GridKernelDescriptor(Grid_layout gd, kernel_launcher::TypeInfo kernel_info): gd(gd), kernel_info(kernel_info) {}

    std::string tuning_key() const override;
    std::string unique_key() const override;
    kernel_launcher::KernelBuilder build_kernel(std::vector<kernel_launcher::TypeInfo> param_types) const override;

private:
    Grid_layout gd;
    kernel_launcher::TypeInfo kernel_info;
};

void launch_kernel(
        cudaStream_t stream,
        dim3 problem_size,
        const KernelDescriptor& kernel,
        const std::vector<const kernel_launcher::KernelArg*>& args,
        const std::function<void()>& fallback
);

template <typename Tuple, size_t... Is>
std::vector<const kernel_launcher::KernelArg*> arg_tuple_to_vector_helper(
        const Tuple& tuple,
        std::index_sequence<Is...>
) {
    return {&std::get<Is>(tuple)...};
}

template <typename... Args>
void launch_kernel(
        cudaStream_t stream,
        dim3 problem_size,
        const KernelDescriptor& kernel,
        const std::function<void()>& fallback,
        Args&&... args
) {
    std::tuple<KernelArgConvert<typename std::remove_reference<Args>::type>...> tuple_args{
            KernelArgConvert<typename std::remove_reference<Args>::type>(args)...
    };

    std::vector<const kernel_launcher::KernelArg*> kernel_args = arg_tuple_to_vector_helper(
            tuple_args, std::index_sequence_for<Args...>{}
    );

    launch_kernel(stream, problem_size, kernel, kernel_args, fallback);
}

template <typename F, typename... Args>
void launch_grid_kernel(
        Grid_layout gd,
        Args&&... args
) {
    dim3 problem_size(gd.iend - gd.istart, gd.jend - gd.jstart, gd.kend - gd.kstart);

    auto fallback = [&]() {
        dim3 block_size(32, 8, 1);
        dim3 grid_size(
                problem_size.x / block_size.x + (problem_size.x % block_size.x != 0),
                problem_size.y / block_size.y + (problem_size.y % block_size.y != 0),
                problem_size.z / block_size.z + (problem_size.z % block_size.z != 0)
        );

        kernel_wrapper<TilingKernel<F>><<<grid_size, block_size>>>(
                gd, typename KernelArgConvert<Args>::type(args)...
        );
    };

    launch_kernel(
            nullptr,
            problem_size,
            GridKernelDescriptor(gd, kernel_launcher::TypeInfo::of<F>()),
            fallback,
            std::forward<Args>(args)...
    );
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

#endif //CUDA_COMPILER_H
