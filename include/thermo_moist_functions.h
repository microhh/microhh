/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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

#ifndef THERMO_MOIST_FUNCTIONS
#define THERMO_MOIST_FUNCTIONS

// In case the code is compiled with NVCC, add the macros for CUDA
#ifdef __CUDACC__
#  define CUDA_MACRO __host__ __device__
#else
#  define CUDA_MACRO
#endif

#include "constants.h"

namespace Thermo_moist_functions
{
    using namespace Constants;

    // INLINE FUNCTIONS
    template<typename TF>
    CUDA_MACRO inline TF buoyancy(const TF exn, const TF thl, const TF qt, const TF ql, const TF thvref)
    {
        return grav * ((thl + Lv<TF>*ql/(cp<TF>*exn)) * (1. - (1. - Rv<TF>/Rd<TF>)*qt - Rv<TF>/Rd<TF>*ql) - thvref) / thvref;
    }

    template<typename TF>
    CUDA_MACRO inline TF virtual_temperature(const TF exn, const TF thl, const TF qt, const TF ql)
    {
        return (thl + Lv<TF>*ql/(cp<TF>*exn)) * (1. - (1. - Rv<TF>/Rd<TF>)*qt - Rv<TF>/Rd<TF>*ql);
    }

    template<typename TF>
    CUDA_MACRO inline TF virtual_temperature_no_ql(const TF exn, const TF thl, const TF qt)
    {
        return thl * (1. - (1. - Rv<TF>/Rd<TF>)*qt);
    }

    template<typename TF>
    CUDA_MACRO inline TF buoyancy_no_ql(const TF thl, const TF qt, const TF thvref)
    {
        return grav * (thl * (1. - (1. - Rv<TF>/Rd<TF>)*qt) - thvref) / thvref;
    }

    template<typename TF>
    CUDA_MACRO inline TF buoyancy_flux_no_ql(const TF thl, const TF thlflux, const TF qt, const TF qtflux, const TF thvref)
    {
        return grav/thvref * (thlflux * (1. - (1.-Rv<TF>/Rd<TF>)*qt) - (1.-Rv<TF>/Rd<TF>)*thl*qtflux);
    }

    //CUDA_MACRO inline TF esat(const TF T)
    //{
    //    #ifdef __CUDACC__
    //    const TF x=fmax(-80.,T-T0);
    //    #else
    //    const TF x=std::max(-80.,T-T0);
    //    #endif

    //    return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))));
    //}

    // Saturation vapor pressure, using Taylor expansion at T=T0 around the Arden Buck (1981) equation:
    // es = 611.21 * exp(17.502 * Tc / (240.97 + Tc)), with Tc=T-T0
    template<typename TF>
    CUDA_MACRO inline TF esat(const TF T)
    {
        #ifdef __CUDACC__
        const TF x=fmax(-75.,T-T0<TF>);
        #else
        const TF x=std::max(-75.,T-T0<TF>);
        #endif

        return c00+x*(c10+x*(c20+x*(c30+x*(c40+x*(c50+x*(c60+x*(c70+x*(c80+x*(c90+x*c100)))))))));
    }

    template<typename TF>
    CUDA_MACRO inline TF qsat(const TF p, const TF T)
    {
        return ep<TF>*esat(T)/(p-(1-ep<TF>)*esat(T));
    }

    template<typename TF>
    CUDA_MACRO inline TF exner(const TF p)
    {
        return pow((p/p0<TF>),(Rd<TF>/cp<TF>));
    }
}
#endif
