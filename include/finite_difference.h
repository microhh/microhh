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

#ifndef FINITE_DIFFERENCE_H
#define FINITE_DIFFERENCE_H

// In case the code is compiled with NVCC, add the macros for CUDA
#ifdef __CUDACC__
#  define CUDA_MACRO __host__ __device__
#else
#  define CUDA_MACRO
#endif

namespace Finite_difference
{
    namespace O2
    {
        template<typename TF>
        CUDA_MACRO inline TF interp2(const TF a, const TF b)
        {
            return TF(0.5) * (a + b);
        }

        template<typename TF>
        CUDA_MACRO inline TF interp22(const TF a, const TF b, const TF c, const TF d)
        {
            return TF(0.25) * (a + b + c + d);
        }

        template<typename TF>
        CUDA_MACRO inline TF grad2(const TF a, const TF b)
        {
            return (b - a);
        }
    }

    namespace O4
    {
        template <typename TF> constexpr TF ci0  = -1./16.;
        template <typename TF> constexpr TF ci1  =  9./16.;
        template <typename TF> constexpr TF ci2  =  9./16.;
        template <typename TF> constexpr TF ci3  = -1./16.;

        template <typename TF> constexpr TF bi0  =  5./16.;
        template <typename TF> constexpr TF bi1  = 15./16.;
        template <typename TF> constexpr TF bi2  = -5./16.;
        template <typename TF> constexpr TF bi3  =  1./16.;

        template <typename TF> constexpr TF ti0  =  1./16.;
        template <typename TF> constexpr TF ti1  = -5./16.;
        template <typename TF> constexpr TF ti2  = 15./16.;
        template <typename TF> constexpr TF ti3  =  5./16.;

        template <typename TF> constexpr TF cg0  =   1./24.;
        template <typename TF> constexpr TF cg1  = -27./24.;
        template <typename TF> constexpr TF cg2  =  27./24.;
        template <typename TF> constexpr TF cg3  =  -1./24.;

        template <typename TF> constexpr TF bg0  = -23./24.;
        template <typename TF> constexpr TF bg1  =  21./24.;
        template <typename TF> constexpr TF bg2  =   3./24.;
        template <typename TF> constexpr TF bg3  =  -1./24.;

        template <typename TF> constexpr TF tg0  =   1./24.;
        template <typename TF> constexpr TF tg1  =  -3./24.;
        template <typename TF> constexpr TF tg2  = -21./24.;
        template <typename TF> constexpr TF tg3  =  23./24.;

        template <typename TF> constexpr TF cdg0 = -1460./576.;
        template <typename TF> constexpr TF cdg1 =   783./576.;
        template <typename TF> constexpr TF cdg2 =   -54./576.;
        template <typename TF> constexpr TF cdg3 =     1./576.;

        template<typename TF>
        CUDA_MACRO inline TF interp4c(const TF a, const TF b, const TF c, const TF d)
        {
            return ci0<TF>*(a+d) + ci1<TF>*(b+c);
        }

        template<typename TF>
        CUDA_MACRO inline TF interp4b(const TF a, const TF b, const TF c, const TF d)
        {
            return bi0<TF>*a + bi1<TF>*b - bi2<TF>*c + bi3<TF>*d;
        }

        template<typename TF>
        CUDA_MACRO inline TF interp4t(const TF a, const TF b, const TF c, const TF d)
        {
            return ti0<TF>*a + ti1<TF>*b + ti2<TF>*c + ti3<TF>*d;
        }

        template<typename TF>
        CUDA_MACRO inline TF interp4_ws(const TF a, const TF b, const TF c, const TF d)
        {
            constexpr TF c0 = TF(7./12.);
            constexpr TF c1 = TF(1./12.);
            return c0*(b+c) - c1*(a+d);
        }

        template<typename TF>
        CUDA_MACRO inline TF interp3_ws(const TF a, const TF b, const TF c, const TF d)
        {
            constexpr TF c0 = TF(3./12.);
            constexpr TF c1 = TF(1./12.);
            return c0*(c-b) - c1*(d-a);
        }

        template<typename TF>
        CUDA_MACRO inline TF grad4(const TF a, const TF b, const TF c, const TF d)
        {
            return ( - cg0<TF>*(d-a) - cg1<TF>*(c-b) );
        }
    }

    namespace O6
    {
        template<typename TF>
        CUDA_MACRO inline TF interp6_ws(
                const TF a, const TF b, const TF c, const TF d, const TF e, const TF f)
        {
            constexpr TF c0 = TF(37./60.);
            constexpr TF c1 = TF(8./60.);
            constexpr TF c2 = TF(1./60.);

            return c0*(c+d) - c1*(b+e) + c2*(a+f);
        }

        template<typename TF>
        CUDA_MACRO inline TF interp5_ws(
                const TF a, const TF b, const TF c, const TF d, const TF e, const TF f)
        {
            constexpr TF c0 = TF(10./60.);
            constexpr TF c1 = TF(5./60.);
            constexpr TF c2 = TF(1./60.);

            return c0*(d-c) - c1*(e-b) + c2*(f-a);
        }
    }
}
#endif
