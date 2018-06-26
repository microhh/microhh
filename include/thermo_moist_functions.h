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
        return grav<TF> * ((thl + Lv<TF>*ql/(cp<TF>*exn)) * (TF(1.) - (TF(1.) - Rv<TF>/Rd<TF>)*qt - Rv<TF>/Rd<TF>*ql) - thvref) / thvref;
    }

    template<typename TF>
    CUDA_MACRO inline TF virtual_temperature(const TF exn, const TF thl, const TF qt, const TF ql)
    {
        return (thl + Lv<TF>*ql/(cp<TF>*exn)) * (TF(1.) - (TF(1.) - Rv<TF>/Rd<TF>)*qt - Rv<TF>/Rd<TF>*ql);
    }

    template<typename TF>
    CUDA_MACRO inline TF virtual_temperature_no_ql(const TF thl, const TF qt)
    {
        return thl * (TF(1.) - (TF(1.) - Rv<TF>/Rd<TF>)*qt);
    }

    template<typename TF>
    CUDA_MACRO inline TF buoyancy_no_ql(const TF thl, const TF qt, const TF thvref)
    {
        return grav<TF> * (thl * (TF(1.) - (TF(1.) - Rv<TF>/Rd<TF>)*qt) - thvref) / thvref;
    }

    template<typename TF>
    CUDA_MACRO inline TF buoyancy_flux_no_ql(const TF thl, const TF thlflux, const TF qt, const TF qtflux, const TF thvref)
    {
        return grav<TF>/thvref * (thlflux * (TF(1.) - (TF(1.)-Rv<TF>/Rd<TF>)*qt) - (TF(1.)-Rv<TF>/Rd<TF>)*thl*qtflux);
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
        const TF x=fmax(TF(-75.), T-T0<TF>);
        #else
        const TF x=std::max(TF(-75.), T-T0<TF>);
        #endif

        return c00<TF>+x*(c10<TF>+x*(c20<TF>+x*(c30<TF>+x*(c40<TF>+x*(c50<TF>+x*(c60<TF>+x*(c70<TF>+x*(c80<TF>+x*(c90<TF>+x*c100<TF>)))))))));
    }

    template<typename TF>
    CUDA_MACRO inline TF qsat(const TF p, const TF T)
    {
        return ep<TF>*esat(T)/(p-(TF(1.)-ep<TF>)*esat(T));
    }

    template<typename TF>
    CUDA_MACRO inline TF exner(const TF p)
    {
        return pow((p/p0<TF>),(Rd<TF>/cp<TF>));
    }


    template<typename TF>
    struct Struct_sat_adjust
    {
        TF ql;
        TF t;
        TF qs;
    };

    template<typename TF>
    inline Struct_sat_adjust<TF> sat_adjust(const TF thl, const TF qt, const TF p, const TF exn)
    {
        int niter = 0;
        int nitermax = 30;
        TF ql, tl, tnr_old = 1.e9, tnr, qs=0;
        tl = thl * exn;
        Struct_sat_adjust<TF> ans;

        // Calculate if q-qs(Tl) <= 0. If so, return 0. Else continue with saturation adjustment
        ans.ql = 0;
        ans.t = tl;
        ans.qs = qsat(p, tl);
        if(qt-ans.qs <= 0)
            return ans;

        tnr = tl;
        while (std::fabs(tnr-tnr_old)/tnr_old> 1e-5 && niter < nitermax)
        {
            ++niter;
            tnr_old = tnr;
            qs = qsat(p, tnr);
            tnr = tnr - (tnr+(Lv<TF>/cp<TF>)*qs-tl-(Lv<TF>/cp<TF>)*qt)/(1+(std::pow(Lv<TF>, 2)*qs)/ (Rv<TF>*cp<TF>*std::pow(tnr, 2)));
        }

        if (niter == nitermax)
            throw std::runtime_error("Non-converging saturation adjustment.");

        ql = std::max(TF(0.), qt - qs);

        ans.ql = ql;
        ans.t  = tnr;
        ans.qs = qs;
        return ans;
    }

    template<typename TF>
    void  calc_base_state(TF* restrict pref,    TF* restrict prefh,
                          TF* restrict rho,     TF* restrict rhoh,
                          TF* restrict thv,     TF* restrict thvh,
                          TF* restrict ex,      TF* restrict exh,
                          TF* restrict thlmean, TF* restrict qtmean, const TF pbot,
                          const int kstart, const int kend,
                          const TF* restrict z, const TF* restrict dz, const TF* const dzh)
    {
        const TF thlsurf = TF(0.5)*(thlmean[kstart-1 ]+ thlmean[kstart]);
        const TF qtsurf  = TF(0.5)*(qtmean[kstart-1] +  qtmean[kstart]);

        TF ql;

        // Calculate the values at the surface (half level == kstart)
        prefh[kstart] = pbot;
        exh[kstart]   = exner(prefh[kstart]);
        ql            = sat_adjust(thlsurf, qtsurf, prefh[kstart], exh[kstart]).ql;
        thvh[kstart]  = virtual_temperature(exh[kstart], thlsurf, qtsurf, ql);
        rhoh[kstart]  = pbot / (Rd<TF> * exh[kstart] * thvh[kstart]);

        // Calculate the first full level pressure
        pref[kstart]  = prefh[kstart] * std::exp(-grav<TF> * z[kstart] / (Rd<TF> * exh[kstart] * thvh[kstart]));

        for (int k=kstart+1; k<kend+1; ++k)
        {
            // 1. Calculate remaining values (thv and rho) at full-level[k-1]
            ex[k-1]  = exner(pref[k-1]);
            ql       = sat_adjust(thlmean[k-1], qtmean[k-1], pref[k-1], ex[k-1]).ql;
            thv[k-1] = virtual_temperature(ex[k-1], thlmean[k-1], qtmean[k-1], ql);
            rho[k-1] = pref[k-1] / (Rd<TF> * ex[k-1] * thv[k-1]);

            // 2. Calculate pressure at half-level[k]
            prefh[k] = prefh[k-1] * std::exp(-grav<TF> * dz[k-1] / (Rd<TF> * ex[k-1] * thv[k-1]));
            exh[k]   = exner(prefh[k]);

            // 3. Use interpolated conserved quantities to calculate half-level[k] values
            const TF thli = TF(0.5)*(thlmean[k-1] + thlmean[k]);
            const TF qti  = TF(0.5)*(qtmean [k-1] + qtmean [k]);
            const TF qli  = sat_adjust(thli, qti, prefh[k], exh[k]).ql;

            thvh[k]  = virtual_temperature(exh[k], thli, qti, qli);
            rhoh[k]  = prefh[k] / (Rd<TF> * exh[k] * thvh[k]);

            // 4. Calculate pressure at full-level[k]
            pref[k] = pref[k-1] * std::exp(-grav<TF> * dzh[k] / (Rd<TF> * exh[k] * thvh[k]));
        }

        pref[kstart-1] = TF(2.)*prefh[kstart] - pref[kstart];
    }
    
    template<typename TF>
    void  calc_base_state_no_ql(TF* restrict pref,    TF* restrict prefh,
                          TF* restrict rho,     TF* restrict rhoh,
                          TF* restrict thv,     TF* restrict thvh,
                          TF* restrict ex,      TF* restrict exh,
                          TF* restrict thlmean, TF* restrict qtmean, const TF pbot,
                          const int kstart, const int kend,
                          const TF* restrict z, const TF* restrict dz, const TF* const dzh)
    {
        const TF thlsurf = TF(0.5)*(thlmean[kstart-1 ]+ thlmean[kstart]);
        const TF qtsurf  = TF(0.5)*(qtmean[kstart-1] +  qtmean[kstart]);

        // Calculate the values at the surface (half level == kstart)
        prefh[kstart] = pbot;
        exh[kstart]   = exner(prefh[kstart]);
        thvh[kstart]  = virtual_temperature_no_ql(thlsurf, qtsurf);
        rhoh[kstart]  = pbot / (Rd<TF> * exh[kstart] * thvh[kstart]);

        // Calculate the first full level pressure
        pref[kstart]  = prefh[kstart] * std::exp(-grav<TF> * z[kstart] / (Rd<TF> * exh[kstart] * thvh[kstart]));

        for (int k=kstart+1; k<kend+1; ++k)
        {
            // 1. Calculate remaining values (thv and rho) at full-level[k-1]
            ex[k-1]  = exner(pref[k-1]);
            thv[k-1] = virtual_temperature_no_ql(thlmean[k-1], qtmean[k-1]);
            rho[k-1] = pref[k-1] / (Rd<TF> * ex[k-1] * thv[k-1]);

            // 2. Calculate pressure at half-level[k]
            prefh[k] = prefh[k-1] * std::exp(-grav<TF> * dz[k-1] / (Rd<TF> * ex[k-1] * thv[k-1]));
            exh[k]   = exner(prefh[k]);

            // 3. Use interpolated conserved quantities to calculate half-level[k] values
            const TF thli = TF(0.5)*(thlmean[k-1] + thlmean[k]);
            const TF qti  = TF(0.5)*(qtmean [k-1] + qtmean [k]);

            thvh[k]  = virtual_temperature_no_ql(thli, qti);
            rhoh[k]  = prefh[k] / (Rd<TF> * exh[k] * thvh[k]);

            // 4. Calculate pressure at full-level[k]
            pref[k] = pref[k-1] * std::exp(-grav<TF> * dzh[k] / (Rd<TF> * exh[k] * thvh[k]));
        }

        pref[kstart-1] = TF(2.)*prefh[kstart] - pref[kstart];
    }
}
#endif
