/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

#ifndef MONIN_OBUKHOV_H
#define MONIN_OBUKHOV_H

#include "fast_math.h"
#include "constants.h"
namespace fm = Fast_math;

// In case the code is compiled with NVCC, add the macros for CUDA
#ifdef __CUDACC__
#  define CUDA_MACRO __host__ __device__
#else
#  define CUDA_MACRO
#endif

namespace Monin_obukhov
{
    //
    // GRADIENT FUNCTIONS
    //
    template<typename TF>
    CUDA_MACRO inline TF phim_unstable(const TF zeta)
    {
        // Wilson, 2001 functions, see Wyngaard, page 222.
        return std::pow(TF(1.) + TF(3.6)*std::pow(std::abs(zeta), TF(2./3.)), TF(-1./2.));
    }

    template<typename TF>
    CUDA_MACRO inline TF phim_stable(const TF zeta)
    {
        // Hogstrom, 1988
        //return TF(1.) + TF(4.8)*zeta;

        // IFS
        return TF(1) + TF(5)*zeta;
    }

    template<typename TF>
    CUDA_MACRO inline TF phim(const TF zeta)
    {
        return (zeta <= TF(0.)) ? phim_unstable(zeta) : phim_stable(zeta);
    }

    template<typename TF>
    CUDA_MACRO inline TF phih_unstable(const TF zeta)
    {
        // Wilson, 2001 functions, see Wyngaard, page 222.
        return std::pow(TF(1.) + TF(7.9)*std::pow(std::abs(zeta), TF(2./3.)), TF(-1./2.));
    }

    template<typename TF>
    CUDA_MACRO inline TF phih_stable(const TF zeta)
    {
        // Hogstrom, 1988
        //return TF(1.) + TF(7.8)*zeta;

        // IFS
        return fm::pow2(TF(1) + TF(4)*zeta);
    }

    template<typename TF>
    CUDA_MACRO inline TF phih(const TF zeta)
    {
        return (zeta <= TF(0.)) ? phih_unstable(zeta) : phih_stable(zeta);
    }

    //
    // INTEGRATED FUNCTIONS
    //
    template<typename TF>
    CUDA_MACRO inline TF psim_unstable(const TF zeta)
    {
        // Wilson, 2001 functions, see Wyngaard, page 222.
        return TF(3.)*std::log( ( TF(1.) + TF(1.)/phim_unstable(zeta) ) / TF(2.));
    }

    template<typename TF>
    CUDA_MACRO inline TF psim_stable(const TF zeta)
    {
        // Hogstrom, 1988
        //return TF(-4.8)*zeta;

        // IFS
        constexpr TF a = TF(1);
        constexpr TF b = TF(2)/TF(3);
        constexpr TF c = TF(5);
        constexpr TF d = TF(0.35);

        return -b * (zeta - (c/d)) * std::exp(-d * zeta) - a*zeta - (b*c)/d;
    }

    template<typename TF>
    CUDA_MACRO inline TF psih_unstable(const TF zeta)
    {
        // Wilson, 2001 functions, see Wyngaard, page 222.
        return TF(3.) * std::log( ( TF(1.) + TF(1.) / phih_unstable(zeta) ) / TF(2.));
    }

    template<typename TF>
    CUDA_MACRO inline TF psih_stable(const TF zeta)
    {
        // Hogstrom, 1988
        //return TF(-7.8)*zeta;

        // IFS
        constexpr TF a = TF(1);
        constexpr TF b = TF(2)/TF(3);
        constexpr TF c = TF(5);
        constexpr TF d = TF(0.35);

        return -b * (zeta - (c/d)) * std::exp(-d * zeta) - std::pow(TF(1)+ b*a*zeta, TF(1.5)) -(b*c)/d + TF(1);

    }

    template<typename TF>
    CUDA_MACRO inline TF fm(const TF zsl, const TF z0m, const TF L)
    {
        return (L <= TF(0.))
            ? Constants::kappa<TF> / (std::log(zsl/z0m) - psim_unstable(zsl/L) + psim_unstable(z0m/L))
            : Constants::kappa<TF> / (std::log(zsl/z0m) - psim_stable  (zsl/L) + psim_stable  (z0m/L));
    }

    template<typename TF>
    CUDA_MACRO inline TF fh(const TF zsl, const TF z0h, const TF L)
    {
        return (L <= TF(0.))
            ? Constants::kappa<TF> / (std::log(zsl/z0h) - psih_unstable(zsl/L) + psih_unstable(z0h/L))
            : Constants::kappa<TF> / (std::log(zsl/z0h) - psih_stable  (zsl/L) + psih_stable  (z0h/L));
    }
}
#endif
