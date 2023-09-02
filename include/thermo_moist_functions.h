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

#ifndef THERMO_MOIST_FUNCTIONS_H
#define THERMO_MOIST_FUNCTIONS_H

// In case the code is compiled with NVCC, add the macros for CUDA
#ifdef __CUDACC__
#  define CUDA_MACRO __host__ __device__
#else
#  define CUDA_MACRO
#endif


#include <iostream>
#include <iomanip>

#include "constants.h"
#include "fast_math.h"

namespace Thermo_moist_functions
{
    using namespace Constants;
    using Fast_math::pow2;

    // INLINE FUNCTIONS
    template<typename TF>
    CUDA_MACRO inline TF virtual_temperature(const TF exn, const TF thl, const TF qt, const TF ql, const TF qi)
    {
        const TF th = thl + Lv<TF>*ql/(cp<TF>*exn) + Ls<TF>*qi/(cp<TF>*exn);
        return th * (TF(1.) - (TF(1.) - Rv<TF>/Rd<TF>)*qt - Rv<TF>/Rd<TF>*(ql+qi));
    }

    template<typename TF>
    CUDA_MACRO inline TF buoyancy(const TF exn, const TF thl, const TF qt, const TF ql, const TF qi, const TF thvref)
    {
        return grav<TF> * (virtual_temperature(exn, thl, qt, ql, qi) - thvref) / thvref;
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

    template<typename TF>
    CUDA_MACRO inline TF virtual_temperature_flux_no_ql(const TF thl, const TF thlflux, const TF qt, const TF qtflux)
    {
        return (thlflux * (TF(1.) - (TF(1.)-Rv<TF>/Rd<TF>)*qt) - (TF(1.)-Rv<TF>/Rd<TF>)*thl*qtflux);
    }

    // Saturation vapor pressure, using Taylor expansion at T=T0 around the Arden Buck (1981) equation:
    // es = 611.21 * exp(17.502 * Tc / (240.97 + Tc)), with Tc=T-T0
    template<typename TF>
    CUDA_MACRO inline TF esat_liq(const TF T)
    {
        #ifdef __CUDACC__
        const TF x = fmax(TF(-75.), T-T0<TF>);
        #else
        const TF x = std::max(TF(-75.), T-T0<TF>);
        #endif

        // return TF(611.21)*std::exp(TF(17.502)*x / (TF(240.97)+x));
        return c00<TF>+x*(c10<TF>+x*(c20<TF>+x*(c30<TF>+x*(c40<TF>+x*(c50<TF>+x*(c60<TF>+x*(c70<TF>+x*(c80<TF>+x*(c90<TF>+x*c100<TF>)))))))));
    }

    template<typename TF>
    CUDA_MACRO inline TF qsat_liq(const TF p, const TF T)
    {
        return ep<TF>*esat_liq(T)/(p-(TF(1.)-ep<TF>)*esat_liq(T));
    }

    // Saturation vapor pressure over ice, Arden Buck (1981) equation:
    // es = 611.15 * exp(22.452 * Tc / (272.55 + Tc)), with Tc=T-T0
    template<typename TF>
    CUDA_MACRO inline TF esat_ice(const TF T)
    {
        #ifdef __CUDACC__
        const TF x = fmax(TF(-100.), T-T0<TF>);
        #else
        const TF x = std::max(TF(-100.), T-T0<TF>);
        #endif

        return TF(611.15)*std::exp(TF(22.452)*x / (TF(272.55)+x));
    }

    template<typename TF>
    CUDA_MACRO inline TF qsat_ice(const TF p, const TF T)
    {
        return ep<TF>*esat_ice(T)/(p-(TF(1.)-ep<TF>)*esat_ice(T));
    }

    // Compute water fraction of condensate following Tomita, 2008.
    template<typename TF>
    CUDA_MACRO inline TF water_fraction(const TF T)
    {
        #ifdef __CUDACC__
        return fmax(TF(0.), fmin((T - TF(233.15)) / (T0<TF> - TF(233.15)), TF(1.)));
        #else
        return std::max(TF(0.), std::min((T - TF(233.15)) / (T0<TF> - TF(233.15)), TF(1.)));
        #endif
    }

    // Combine the ice and water saturated specific humidities following Tomita, 2008.
    template<typename TF>
    CUDA_MACRO inline TF qsat(const TF p, const TF T)
    {
        const TF alpha = water_fraction(T);
        return alpha*qsat_liq(p, T) + (TF(1.)-alpha)*qsat_ice(p, T);
    }

    template<typename TF>
    CUDA_MACRO inline TF esat(const TF T)
    {
        const TF alpha = water_fraction(T);
        return alpha*esat_liq(T) + (TF(1.)-alpha)*esat_ice(T);
    }

    template<typename TF>
    CUDA_MACRO inline TF dqsatdT_liq(const TF p, const TF T)
    {
        const TF den = p - esat_liq(T)*(TF(1.) - ep<TF>);
        return (ep<TF>/den - (TF(1.) + ep<TF>)*ep<TF>*esat_liq(T)/pow2(den)) * Lv<TF>*esat_liq(T) / (Rv<TF>*pow2(T));
    }

    template<typename TF>
    CUDA_MACRO inline TF dqsatdT_ice(const TF p, const TF T)
    {
        const TF den = p - esat_ice(T)*(TF(1.) - ep<TF>);
        return (ep<TF>/den + (TF(1.) - ep<TF>)*ep<TF>*esat_ice(T)/pow2(den)) * Ls<TF>*esat_ice(T) / (Rv<TF>*pow2(T));
    }

    template<typename TF>
    CUDA_MACRO inline TF dqsatdT(const TF p, const TF T)
    {
        const TF alpha = water_fraction(T);
        return alpha*dqsatdT_liq(p,T) + (TF(1.)-alpha)*dqsatdT_ice(p,T);
    }

    template<typename TF>
    CUDA_MACRO inline TF exner(const TF p)
    {
        return pow((p/p0<TF>), (Rd<TF>/cp<TF>));
    }

    template<typename TF>
    struct Struct_sat_adjust
    {
        TF ql;
        TF qi;
        TF t;
        TF qs;
    };

    template<typename TF>
    inline Struct_sat_adjust<TF> sat_adjust(const TF thl, const TF qt, const TF p, const TF exn)
    {
        using Fast_math::pow2;

        int niter = 0;
        int nitermax = 10;
        TF tnr_old = TF(1.e9);

        const TF tl = thl * exn;
        TF qs = qsat_liq(p, tl);

        Struct_sat_adjust<TF> ans =
        {
            TF(0.), // ql
            TF(0.), // qi
            tl, // t
            qs, // qs
        };

        // Calculate if q-qs(Tl) <= 0. If so, return 0. Else continue with saturation adjustment.
        if (qt-ans.qs <= TF(0.))
            return ans;

        /* Saturation adjustment solver.
         * Root finding function is f(T) = T - tnr - Lv/cp*qt + alpha_w * Lv/cp*qs(T) + alpha_i*Ls/cp*qs(T)
         * dq_sat/dT derivatives can be rewritten using Claussius-Clapeyron (desat/dT = L{v,s}*esat / (Rv*T^2)).
         */

        TF tnr = tl;

        // Warm adjustment.
        if (tl >= T0<TF>)
        {
            while (std::fabs(tnr-tnr_old)/tnr_old > TF(1.e-5) && niter < nitermax)
            {
                ++niter;
                tnr_old = tnr;
                qs = qsat_liq(p, tnr);
                const TF f =
                    tnr - tl - Lv<TF>/cp<TF>*(qt - qs);

                const TF f_prime = TF(1.) + Lv<TF>/cp<TF>*dqsatdT_liq(p, tnr);

                tnr -= f / f_prime;
            }

            qs = qsat_liq(p, tnr);
            ans.ql = std::max(TF(0.), qt - qs);
            ans.t  = tnr;
            ans.qs = qs;
        }
        // Cold adjustment.
        else
        {
            while (std::fabs(tnr-tnr_old)/tnr_old > TF(1.e-5) && niter < nitermax)
            {
                ++niter;
                tnr_old = tnr;
                qs = qsat(p, tnr);
                const TF alpha_w = water_fraction(tnr);
                const TF alpha_i = TF(1.) - alpha_w;
                const TF dalphadT = (alpha_w > TF(0.) && alpha_w < TF(1.)) ? TF(0.025) : TF(0.);
                const TF dqsatdT_w = dqsatdT_liq(p, tnr);
                const TF dqsatdT_i = dqsatdT_ice(p, tnr);

                const TF f =
                    tnr - tl - alpha_w*Lv<TF>/cp<TF>*qt - alpha_i*Ls<TF>/cp<TF>*qt
                             + alpha_w*Lv<TF>/cp<TF>*qs + alpha_i*Ls<TF>/cp<TF>*qs;

                const TF f_prime = TF(1.)
                    - dalphadT*Lv<TF>/cp<TF>*qt + dalphadT*Ls<TF>/cp<TF>*qt
                    + dalphadT*Lv<TF>/cp<TF>*qs - dalphadT*Ls<TF>/cp<TF>*qs
                    + alpha_w*Lv<TF>/cp<TF>*dqsatdT_w
                    + alpha_i*Ls<TF>/cp<TF>*dqsatdT_i;

                tnr -= f / f_prime;
            }

            const TF alpha_w = water_fraction(tnr);
            const TF alpha_i = TF(1.) - alpha_w;

            qs = qsat(p, tnr);
            const TF ql_qi = std::max(TF(0.), qt - qs);

            ans.ql = alpha_w*ql_qi;
            ans.qi = alpha_i*ql_qi;
            ans.t  = tnr;
            ans.qs = qs;
        }

        if (niter == nitermax)
        {
            std::string error = "Non-converging saturation adjustment: thl, qt, p = "
                + std::to_string(thl) + ", " + std::to_string(qt) + ", " + std::to_string(p);

            #ifdef USEMPI
            std::cout << "SINGLE PROCESS EXCEPTION: " << error << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            #else
            throw std::runtime_error(error);
            #endif
        }

        return ans;
    }

    template<typename TF>
    void calc_base_state(TF* restrict pref, TF* restrict prefh,
                         TF* restrict rho, TF* restrict rhoh,
                         TF* restrict thv, TF* restrict thvh,
                         TF* restrict ex, TF* restrict exh,
                         const TF* restrict thlmean, const TF* restrict qtmean, const TF pbot,
                         const int kstart, const int kend,
                         const TF* restrict z, const TF* restrict dz, const TF* const dzh)
    {
        const TF thlsurf = TF(0.5)*(thlmean[kstart-1] + thlmean[kstart]);
        const TF qtsurf  = TF(0.5)*(qtmean [kstart-1] + qtmean[kstart]);

        // Calculate the values at the surface (half level == kstart)
        prefh[kstart] = pbot;
        exh[kstart]   = exner(prefh[kstart]);

        Struct_sat_adjust<TF> ssa = sat_adjust(thlsurf, qtsurf, prefh[kstart], exh[kstart]);
        TF ql = ssa.ql;
        TF qi = ssa.qi;

        thvh[kstart] = virtual_temperature(exh[kstart], thlsurf, qtsurf, ql, qi);
        rhoh[kstart] = pbot / (Rd<TF> * exh[kstart] * thvh[kstart]);

        // Calculate the first full level pressure
        pref[kstart]  = prefh[kstart] * std::exp(-grav<TF> * z[kstart] / (Rd<TF> * exh[kstart] * thvh[kstart]));

        for (int k=kstart+1; k<kend+1; ++k)
        {
            // 1. Calculate remaining values (thv and rho) at full-level[k-1]
            ex[k-1]  = exner(pref[k-1]);
            ssa      = sat_adjust(thlmean[k-1], qtmean[k-1], pref[k-1], ex[k-1]);
            ql       = ssa.ql;
            qi       = ssa.qi;
            thv[k-1] = virtual_temperature(ex[k-1], thlmean[k-1], qtmean[k-1], ql, qi);
            rho[k-1] = pref[k-1] / (Rd<TF> * ex[k-1] * thv[k-1]);

            // 2. Calculate pressure at half-level[k]
            prefh[k] = prefh[k-1] * std::exp(-grav<TF> * dz[k-1] / (Rd<TF> * ex[k-1] * thv[k-1]));
            exh[k]   = exner(prefh[k]);

            // 3. Use interpolated conserved quantities to calculate half-level[k] values
            const TF thli = TF(0.5)*(thlmean[k-1] + thlmean[k]);
            const TF qti  = TF(0.5)*(qtmean [k-1] + qtmean [k]);

            ssa = sat_adjust(thli, qti, prefh[k], exh[k]);
            const TF qli = ssa.ql;
            const TF qii = ssa.qi;

            thvh[k]  = virtual_temperature(exh[k], thli, qti, qli, qii);
            rhoh[k]  = prefh[k] / (Rd<TF> * exh[k] * thvh[k]);

            // 4. Calculate pressure at full-level[k]
            pref[k] = pref[k-1] * std::exp(-grav<TF> * dzh[k] / (Rd<TF> * exh[k] * thvh[k]));
        }

        pref[kstart-1] = TF(2.)*prefh[kstart] - pref[kstart];
    }

    template<typename TF>
    void calc_base_state_no_ql(TF* restrict pref,    TF* restrict prefh,
                               TF* restrict rho,     TF* restrict rhoh,
                               TF* restrict thv,     TF* restrict thvh,
                               TF* restrict ex,      TF* restrict exh,
                               TF* restrict thlmean, TF* restrict qtmean, const TF pbot,
                               const int kstart, const int kend,
                               const TF* restrict z, const TF* restrict dz, const TF* const dzh)
    {
        const TF thlsurf = TF(0.5)*(thlmean[kstart-1] + thlmean[kstart]);
        const TF qtsurf  = TF(0.5)*(qtmean[kstart-1] + qtmean[kstart]);

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
