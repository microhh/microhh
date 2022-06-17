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

#include <utility>
#include "cuda_buffer.h"


#ifdef USECUDA
#include <cuda_runtime_api.h>
#include "tools.h"
#endif //USECUDA

cuda_raw_buffer::cuda_raw_buffer(size_t size_in_bytes)
{
    resize(size_in_bytes);
}

cuda_raw_buffer::cuda_raw_buffer(cuda_raw_buffer&& that) noexcept
{
    *this = std::move(that);
}

cuda_raw_buffer& cuda_raw_buffer::operator=(cuda_raw_buffer&& that)
{
    std::swap(ptr_, that.ptr_);
    return *this;
}

void cuda_raw_buffer::resize(size_t size_in_bytes)
{
#ifdef USECUDA
    if (ptr_)
    {
        cuda_safe_call(cudaFree(ptr_));
        ptr_ = nullptr;
    }

    if (size_in_bytes > 0)
    {
        cuda_safe_call(cudaMalloc(&ptr_, size_in_bytes));
    }
#else
    if (size_in_bytes > 0)
    {
        throw std::runtime_exception("CUDA is not enabled, allocating GPU memory is not possible");
    }
#endif //USECUDA
}
