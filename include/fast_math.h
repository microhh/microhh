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

#ifndef FAST_MATH_H
#define FAST_MATH_H

// In case the code is compiled with NVCC, add the macros for CUDA
#ifdef __CUDACC__
#  define CUDA_MACRO __host__ __device__
#else
#  define CUDA_MACRO
#endif

namespace Fast_math
{
    template<typename TF>
    CUDA_MACRO inline TF pow2(const TF a)
    {
        return a*a;
    }

    template<typename TF>
    CUDA_MACRO inline TF pow3(const TF a)
    {
        return a*a*a;
    }

    template<typename TF>
    CUDA_MACRO inline TF pow4(const TF a)
    {
        return a*a*a*a;
    }

    template<typename TF>
    CUDA_MACRO inline TF pow7(const TF a)
    {
        return a*a*a*a*a*a*a;
    }

    template<typename TF>
    CUDA_MACRO inline TF pow9(const TF a)
    {
        return a*a*a*a*a*a*a*a*a;
    }
}
#endif
