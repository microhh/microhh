/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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
#include "thermo_moist.h"

namespace Thermo_moist_functions
{
    using namespace Constants;
    using Fast_math::pow2;

    // INLINE FUNCTIONS
    template<typename TF, Satadjust_type sw_satadjust>
    CUDA_MACRO inline TF virtual_temperature(const TF exn, const TF thl, const TF qt, const TF ql, const TF qi,
                                             const TF T, const TF qhm)
    {
        if (sw_satadjust == Satadjust_type::Disabled)
        {
            return thl * (TF(1.) - (TF(1.) - Rv<TF>/Rd<TF>)*qt);
        }
        else if (sw_satadjust == Satadjust_type::Liquid_shallow || sw_satadjust == Satadjust_type::Liquid_deep)
        {
            const TF th = T / exn;
            return th * (TF(1.) - (TF(1.) - Rv<TF>/Rd<TF>)*qt - Rv<TF>/Rd<TF>*(ql) - qhm);
        }
        else if (sw_satadjust == Satadjust_type::Liquid_ice || sw_satadjust == Satadjust_type::Liquid_ice_deep)
        {
            const TF th = T / exn;
            return th * (TF(1.) - (TF(1.) - Rv<TF>/Rd<TF>)*qt - Rv<TF>/Rd<TF>*(ql+qi) - qhm);
        }
    }

    template<typename TF>
    CUDA_MACRO inline TF virtual_temperature_no_ql(const TF thl, const TF qt)
    {
        return thl * (TF(1.) - (TF(1.) - Rv<TF>/Rd<TF>)*qt);
    }

    template<typename TF, Satadjust_type sw_satadjust>
    CUDA_MACRO inline TF buoyancy(const TF exn, const TF thl, const TF qt, const TF ql, const TF qi, const TF thvref, const TF T, const TF qhm)
    {
        return grav<TF> * (virtual_temperature<TF, sw_satadjust>(exn, thl, qt, ql, qi, T, qhm) - thvref) / thvref;
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
        // const TF x = fmax(TF(-75.), T-T0<TF>);
        const TF x = fmin(fmax(TF(-75.), T-T0<TF>), TF(50.));       // Limit the temperature range to avoid numerical errors
        #else
        // const TF x = std::max(TF(-75.), T-T0<TF>);
        // const TF x = std::min(std::max(TF(-75.), T-T0<TF>), TF(50.));     // Limit the temperature range to avoid numerical errors
        const TF x = std::min(std::max(TF(-100.), T-T0<TF>), TF(50.));
        #endif

        return TF(611.21)*std::exp(TF(17.502)*x / (TF(240.97)+x));
        //return c00<TF>+x*(c10<TF>+x*(c20<TF>+x*(c30<TF>+x*(c40<TF>+x*(c50<TF>+x*(c60<TF>+x*(c70<TF>+x*(c80<TF>+x*(c90<TF>+x*c100<TF>)))))))));
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
        // const TF x = fmax(TF(-100.), T-T0<TF>);
        const TF x = fmin(fmax(TF(-100.), T-T0<TF>), TF(50.));     // Limit the temperature range to avoid numerical errors
        #else
        // const TF x = std::max(TF(-100.), T-T0<TF>);
        const TF x = std::min(std::max(TF(-100.), T-T0<TF>), TF(50.));     // Limit the temperature range to avoid numerical errors
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
        return (ep<TF>/den + (TF(1.) - ep<TF>)*ep<TF>*esat_liq(T)/pow2(den)) * Lv<TF>*esat_liq(T) / (Rv<TF>*pow2(T));
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
    CUDA_MACRO inline TF f_D(const TF p, const TF T, const TF qt, const TF tl)
    {
        const TF qs = qsat_liq(p, T);
        const TF f = T - tl - Lv<TF>/cp<TF>*(qt - qs);

        return f;
    }

    template<typename TF>
    CUDA_MACRO inline TF f_D_deep(const TF p, const TF T, const TF qt, const TF tl)
    {
        const TF alpha_w = water_fraction(T);
        const TF alpha_i = TF(1.) - alpha_w;
        const TF qs = qsat(p, T);
        const TF f = T - tl - alpha_w*Lv<TF>/cp<TF>*qt - alpha_i*Ls<TF>/cp<TF>*qt
                + alpha_w*Lv<TF>/cp<TF>*qs + alpha_i*Ls<TF>/cp<TF>*qs;

        return f;
    }

    template<typename TF>
    CUDA_MACRO inline TF f_E(const TF p, const TF T, const TF qt, const TF tl)
    {

        const TF ql = qt - qsat_liq(p, T);
        const TF f = - tl + T * pow(1 + (Lv<TF> * ql) / (cp<TF> * T), -1) ;

        return f;
    }

    template<typename TF>
    CUDA_MACRO inline TF f_E_deep(const TF p, const TF T, const TF qt, const TF tl)
    {

        const TF alpha_w = water_fraction(T);
        const TF alpha_i = TF(1.) - alpha_w;
        const TF qs = qsat(p, T);
        const TF ql = alpha_w * (qt - qs);
        const TF qi = alpha_i * (qt - qs);

        const TF f = - tl + T * pow(1 + (Lv<TF> * ql) / (cp<TF> * T) + (Ls<TF> * qi)/(cp<TF> * T), -1) ;

        return f;
    }

    template<typename TF>
    CUDA_MACRO inline TF f_F(const TF p, const TF T, const TF qt, const TF tl)
    {

        const TF ql = qt - qsat_liq(p, T);
        const TF f = - tl + T * pow(1 + (Lv<TF> * ql) / (cp<TF> * std::max(T, TF(253))), -1) ;

        return f;
    }

    template<typename TF>
    CUDA_MACRO inline TF f_F_deep(const TF p, const TF T, const TF qt, const TF tl)
    {

        const TF alpha_w = water_fraction(T);
        const TF alpha_i = TF(1.) - alpha_w;
        const TF qs = qsat(p, T);
        const TF ql = alpha_w * (qt - qs);
        const TF qi = alpha_i * (qt - qs);

        const TF f = - tl + T * pow(1 + (Lv<TF> * ql) / (cp<TF> * std::max(T, TF(253))) + (Ls<TF> * qi)/(cp<TF> * std::max(T, TF(253))), -1) ;

        return f;
    }

    template<typename TF>
    CUDA_MACRO inline TF f_G(const TF p, const TF T, const TF qt, const TF thl)
    {
        const TF chi = (Rd<TF> + Rv<TF> * qt) / (cp<TF> + cpv<TF> * qt);
        const TF gamma = (Rv<TF> * qt) / (cp<TF> + cpv<TF> * qt);
        const TF epsilon = Rd<TF>  / Rv<TF>;
        const TF ql = qt - qsat_liq(p, T);

        const TF cl = TF(4186);
        const TF lv1 = Lv<TF> + (cl - cpv<TF>) * T0<TF>;
        const TF lv2 = cl - cpv<TF>;
        const TF Lv_T = lv1 - lv2 * T;

        const TF f = -thl + T * pow((p0<TF>/p), chi) * pow((1 - ql / (epsilon + qt)), chi) * pow((1 - ql / qt), -gamma)
                            * std::exp(((-Lv_T * ql) / ((cp<TF> + cpv<TF> * qt) * T)));
        return f;
    }

    template<typename TF>
    CUDA_MACRO inline TF f_G_deep(const TF p, const TF T, const TF qt, const TF thl)
    {
        const TF alpha_w = water_fraction(T);
        const TF alpha_i = TF(1.) - alpha_w;
        const TF qs = qsat(p, T);
        const TF ql = alpha_w * (qt - qs);
        const TF qi = alpha_i * (qt - qs);

        const TF chi = (Rd<TF> + Rv<TF> * qt) / (cp<TF> + cpv<TF> * qt);
        const TF gamma = (Rv<TF> * qt) / (cp<TF> + cpv<TF> * qt);
        const TF epsilon = Rd<TF>  / Rv<TF>;

        const TF cl = TF(4186);
        const TF lv1 = Lv<TF> + (cl - cpv<TF>) * T0<TF>;
        const TF lv2 = cl - cpv<TF>;
        const TF Lv_T = lv1 - lv2 * T;

        const TF ci = TF(2106);
        const TF ls1 = Ls<TF> + (ci - cpv<TF>) * T0<TF>;
        const TF ls2 = ci - cpv<TF>;
        const TF Ls_T = ls1 - ls2 * T;

        const TF f = -thl + T * pow((p0<TF>/p), chi)
                    * pow((1 - (ql + qi) / (epsilon + qt)), chi)
                    * pow((1 - (ql + qi) / qt), -gamma)
                    * std::exp(((-Lv_T * ql - Ls_T * qi) / ((cp<TF> + cpv<TF> * qt) * T)));
        return f;
    }

    template<typename TF>
    struct Struct_sat_adjust
    {
        TF ql;
        TF qi;
        TF t;
        TF qs;
    };

    template<typename TF, Satadjust_type sw_satadjust>
    inline Struct_sat_adjust<TF> sat_adjust(
            const TF thl, const TF qt, const TF p, const TF exn)
    {
        // saturation adjustment for different formulations of thl following BF04 (Bryan, G. H., & Fritsch, J. M. (2004).
        // A reevaluation of ice–liquid water potential temperature. Monthly weather review, 132(10), 2421-2431.)
        // BF04 D is the formulation from Betts, A. K. (1973). Non‐precipitating cumulus convection and its parameterization.
        // Quarterly Journal of the Royal Meteorological Society, 99(419), 178-196.

        int niter = 0;
        int nitermax = 10;
        TF tnr_old = TF(1.e9);

        TF tl;

        if (sw_satadjust == Satadjust_type::Liquid_deep || sw_satadjust == Satadjust_type::Liquid_ice_deep)
        {
            // BF04 E, F
            tl = thl * exn;

            // BF04 G
            // const TF chi = (Rd<TF> + Rv<TF> * qt) / (cp<TF> + cpv<TF> * qt);
            // tl = thl / std::pow(p0<TF>/p, chi);
        }
        else
        {
            // BF04 D
            tl = thl * exn;
        }

        TF qs = qsat_liq(p, tl);

        Struct_sat_adjust<TF> ans =
        {
            TF(0.), // ql
            TF(0.), // qi
            tl,     // t
            qs,     // qs
        };

        // Calculate if q-qs(Tl) <= 0. If so, return 0. Else continue with saturation adjustment.
        if (sw_satadjust == Satadjust_type::Disabled || qt-ans.qs <= TF(0.))
            return ans;

        /* Saturation adjustment solver.
         * Root finding function is f(T) = T - tnr - Lv/cp*qt + alpha_w * Lv/cp*qs(T) + alpha_i*Ls/cp*qs(T)
         * dq_sat/dT derivatives can be rewritten using Claussius-Clapeyron (desat/dT = L{v,s}*esat / (Rv*T^2)).
         */

        TF tnr = tl;

        if (sw_satadjust == Satadjust_type::Liquid_shallow)
        {
            // Warm adjustment.
            while (std::fabs(tnr-tnr_old)/tnr_old > TF(1.e-5) && niter < nitermax)
            {
                ++niter;
                tnr_old = tnr;
                // const TF epsilon = 0.01;

                // BF04 D
                qs = qsat_liq(p, tnr);
                const TF f = tnr - tl - Lv<TF>/cp<TF>*(qt - qs);
                const TF f_prime = TF(1.) + Lv<TF>/cp<TF>*dqsatdT_liq(p, tnr);
                // const TF f = f_D(p, tnr, qt, tl);
                // const TF f_prime = (f_D(p, tnr+epsilon, qt, tl) - f)/epsilon;

                tnr -= f / f_prime;
            }

            qs = qsat_liq(p, tnr);
            ans.ql = std::max(TF(0.), qt - qs);
            ans.t  = tnr;
            ans.qs = qs;
        }

        else if (sw_satadjust == Satadjust_type::Liquid_deep)
        {
            // Warm adjustment.
            while (std::fabs(tnr-tnr_old)/tnr_old > TF(1.e-5) && niter < nitermax)
            {
                ++niter;
                tnr_old = tnr;
                const TF epsilon = 0.1;

                // BF04 E
                const TF f = f_E(p, tnr, qt, tl);
                const TF f_prime = (f_E(p, tnr + epsilon, qt, tl) - f)/epsilon;

                // BF04 F
                // const TF f = f_F(p, tnr, qt, tl);
                // const TF f_prime = (f_F(p, tnr + epsilon, qt, tl) - f)/epsilon;

                //BF04 G
                // const TF f = f_G(p, tnr, qt, thl);
                // const TF f_prime = (f_G(p, tnr+epsilon, qt, thl) - f)/epsilon;

                tnr -= f / f_prime;
            }

            qs = qsat_liq(p, tnr);
            ans.ql = std::max(TF(0.), qt - qs);
            ans.t  = tnr;
            ans.qs = qs;
        }

        else if (sw_satadjust == Satadjust_type::Liquid_ice)
        {
            if (tl >= T0<TF>)
            {
                // Warm adjustment.
                while (std::fabs(tnr-tnr_old)/tnr_old > TF(1.e-5) && niter < nitermax)
                {
                    ++niter;
                    tnr_old = tnr;

                    qs = qsat_liq(p, tnr);
                    const TF f = tnr - tl - Lv<TF>/cp<TF>*(qt - qs);
                    const TF f_prime = TF(1.) + Lv<TF>/cp<TF>*dqsatdT_liq(p, tnr);

                    tnr -= f / f_prime;
                }
            }
            else
            {
                // Cold adjustment.
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
            }

            const TF alpha_w = water_fraction(tnr);
            const TF alpha_i = TF(1.) - alpha_w;

            qs = qsat(p, tnr);
            const TF qlqi = std::max(TF(0.), qt - qs);

            ans.ql = alpha_w*qlqi;
            ans.qi = alpha_i*qlqi;
            ans.t  = tnr;
            ans.qs = qs;
        }
        else    // sw_satadjust == Satadjust_type::Liquid_ice_deep
        {
            if (tl >= T0<TF>)
            {
                // Warm adjustment.
                while (std::fabs(tnr-tnr_old)/tnr_old > TF(1.e-5) && niter < nitermax)
                {
                    ++niter;
                    tnr_old = tnr;

                    const TF epsilon = 0.1;

                    // BF04 D
                    // const TF f = f_D(p, tnr, qt, tl);
                    // const TF f_prime = (f_D(p, tnr+epsilon, qt, tl) - f)/epsilon;

                    // BF04 E
                    const TF f = f_E(p, tnr, qt, tl);
                    const TF f_prime = (f_E(p, tnr + epsilon, qt, tl) - f)/epsilon;

                    // BF04 F
                    // const TF f = f_F(p, tnr, qt, tl);
                    // const TF f_prime = (f_F(p, tnr + epsilon, qt, tl) - f)/epsilon;

                    //BF04 G
                    // const TF f = f_G(p, tnr, qt, thl);
                    // const TF f_prime = (f_G(p, tnr+epsilon, qt, thl) - f)/epsilon;

                    tnr -= f / f_prime;
                }
            }
            else
            {
                // Cold adjustment.
                while (std::fabs(tnr-tnr_old)/tnr_old > TF(1.e-5) && niter < nitermax)
                {
                    ++niter;
                    tnr_old = tnr;

                    const TF epsilon = 0.1;

                    // const TF f = f_D_deep(p, tnr, qt, tl);
                    // const TF f_prime = (f_D_deep(p, tnr + epsilon, qt, tl) - f)/epsilon;

                    // BF04 E
                    const TF f = f_E_deep(p, tnr, qt, tl);
                    const TF f_prime = (f_E_deep(p, tnr + epsilon, qt, tl) - f)/epsilon;

                    // BF04 F
                    // const TF f = f_F_deep(p, tnr, qt, tl);
                    // const TF f_prime = (f_F_deep(p, tnr + epsilon, qt, tl) - f)/epsilon;

                    //BF04 G
                    // const TF f = f_G_deep(p, tnr, qt, thl);
                    // const TF f_prime = (f_G_deep(p, tnr+epsilon, qt, thl) - f)/epsilon;

                    tnr -= f / f_prime;
                }
            }

            const TF alpha_w = water_fraction(tnr);
            const TF alpha_i = TF(1.) - alpha_w;

            qs = qsat(p, tnr);
            const TF qlqi = std::max(TF(0.), qt - qs);

            ans.ql = alpha_w*qlqi;
            ans.qi = alpha_i*qlqi;
            ans.t  = tnr;
            ans.qs = qs;
        }

        if (niter == nitermax)
        {
            std::string error = "Non-converging saturation-adjustment. Input: thl="
                    + std::to_string(thl) + " K, qt="
                    + std::to_string(qt) + " kg/kg, p="
                    + std::to_string(p) + " Pa";

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
    inline Struct_sat_adjust<TF> sat_adjust_absolute_T(
            const TF T, const TF qt, const TF p, const TF qc, const TF qv, const TF Lv, const TF cp)
    {
        int niter = 0;
        int nitermax = 10;
        TF tnr_old = TF(1.e9);

        // the essential difference with the satadjust above is in the tl here:
        TF tl = T - (Lv / cp) * qc;
        TF qs = qsat_liq(p, tl);

        Struct_sat_adjust<TF> ans =
                {
                        TF(0.), // ql
                        TF(0.), // qi
                        tl,     // t
                        qs,     // qs
                };

        // Calculate if q-qs(Tl) <= 0. If so, return 0. Else continue with saturation adjustment.
        if (qt-ans.qs <= TF(0.))
            return ans;

        else {
            /* Saturation adjustment solver.
             * Root finding function is f(T) = T - tnr - Lv/cp*qt + alpha_w * Lv/cp*qs(T) + alpha_i*Ls/cp*qs(T)
             * dq_sat/dT derivatives can be rewritten using Claussius-Clapeyron (desat/dT = L{v,s}*esat / (Rv*T^2)).
             */

            TF tnr = tl;

            // Warm adjustment.
            while (std::fabs(tnr - tnr_old) / tnr_old > TF(1.e-5) && niter < nitermax) {
                ++niter;
                tnr_old = tnr;

                qs = qsat_liq(p, tnr);
                const TF f = tnr - tl - Lv / cp * (qt - qs);
                const TF f_prime = TF(1.) + Lv / cp * dqsatdT_liq(p, tnr);

                tnr -= f / f_prime;
            }

            qs = qsat_liq(p, tnr);

            ans.ql = std::max(TF(0.), qt - qs);
            ans.qi = TF(0.);
            ans.t = tnr;
            ans.qs = qs;

            if (niter == nitermax) {
                std::string error = "Non-converging saturation-adjustment. Input: T="
                                    + std::to_string(T) + " K, qt="
                                    + std::to_string(qt) + " kg/kg, p="
                                    + std::to_string(p) + " Pa";

                #ifdef USEMPI
                std::cout << "SINGLE PROCESS EXCEPTION: " << error << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
                #else
                throw std::runtime_error(error);
                #endif
            }

            return ans;
        }
    }


    template<typename TF>
    inline Struct_sat_adjust<TF> sat_adjust_absolute_T_ice(
            const TF T, const TF qt, const TF p, const TF qc, const TF qv, const TF qi, const TF Lv, const TF Ls, const TF cp)
    {
        int niter = 0;
        int nitermax = 10;
        TF tnr_old = TF(1.e9);

        // the essential difference with the satadjust above is in the tl here:
        TF tl = T - (Lv / cp) * qc -  Ls / cp * qi;
        TF qs = qsat_liq(p, tl);

        Struct_sat_adjust<TF> ans =
                {
                        TF(0.), // ql
                        TF(0.), // qi
                        tl,     // t
                        qs,     // qs
                };

        // Calculate if q-qs(Tl) <= 0. If so, return 0. Else continue with saturation adjustment.
        if (qt-ans.qs <= TF(0.))
            return ans;

        /* Saturation adjustment solver.
         * Root finding function is f(T) = T - tnr - Lv/cp*qt + alpha_w * Lv/cp*qs(T) + alpha_i*Ls/cp*qs(T)
         * dq_sat/dT derivatives can be rewritten using Claussius-Clapeyron (desat/dT = L{v,s}*esat / (Rv*T^2)).
         */

        TF tnr = tl;

        if (tl >= T0<TF>)
        {
            // Warm adjustment.
            while (std::fabs(tnr-tnr_old)/tnr_old > TF(1.e-5) && niter < nitermax)
            {
                ++niter;
                tnr_old = tnr;

                qs = qsat_liq(p, tnr);
                const TF f = tnr - tl - Lv/cp*(qt - qs);
                const TF f_prime = TF(1.) + Lv/cp*dqsatdT_liq(p, tnr);

                tnr -= f / f_prime;
            }
        }
        else
        {
            // Cold adjustment.
            while (std::fabs(tnr - tnr_old) / tnr_old > TF(1.e-5) && niter < nitermax) {
                ++niter;
                tnr_old = tnr;
                qs = qsat(p, tnr);
                const TF alpha_w = water_fraction(tnr);
                const TF alpha_i = TF(1.) - alpha_w;
                const TF dalphadT = (alpha_w > TF(0.) && alpha_w < TF(1.)) ? TF(0.025) : TF(0.);
                const TF dqsatdT_w = dqsatdT_liq(p, tnr);
                const TF dqsatdT_i = dqsatdT_ice(p, tnr);

                const TF f =
                        tnr - tl - alpha_w * Lv / cp * qt - alpha_i * Ls / cp * qt
                        + alpha_w * Lv / cp * qs + alpha_i * Ls / cp * qs;

                const TF f_prime = TF(1.)
                                   - dalphadT * Lv / cp * qt + dalphadT * Ls / cp * qt
                                   + dalphadT * Lv / cp * qs - dalphadT * Ls / cp * qs
                                   + alpha_w * Lv / cp * dqsatdT_w
                                   + alpha_i * Ls / cp * dqsatdT_i;

                tnr -= f / f_prime;
            }
        }

        const TF alpha_w = water_fraction(tnr);
        const TF alpha_i = TF(1.) - alpha_w;

        qs = qsat(p, tnr);
        const TF qlqi = std::max(TF(0.), qt - qs);

        ans.ql = alpha_w*qlqi;
        ans.qi = alpha_i*qlqi;
        ans.t  = tnr;
        ans.qs = qs;


        if (niter == nitermax)
        {
            std::string error = "Non-converging saturation-adjustment. Input: T="
                                + std::to_string(T) + " K, qt="
                                + std::to_string(qt) + " kg/kg, p="
                                + std::to_string(p) + " Pa";

            #ifdef USEMPI
            std::cout << "SINGLE PROCESS EXCEPTION: " << error << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            #else
            throw std::runtime_error(error);
            #endif
        }

        return ans;
    }


    template<typename TF, Satadjust_type sw_satadjust>
    void calc_base_state(
            TF* restrict pref,
            TF* restrict prefh,
            TF* restrict rho,
            TF* restrict rhoh,
            TF* restrict thv,
            TF* restrict thvh,
            TF* restrict ex,
            TF* restrict exh,
            const TF* restrict thlmean,
            const TF* restrict qtmean,
            const TF* restrict qhmmean,
            const TF pbot,
            const int kstart,
            const int kend,
            const TF* restrict z,
            const TF* restrict dz,
            const TF* restrict dzh)
    {
        const TF thlsurf = TF(0.5)*(thlmean[kstart-1] + thlmean[kstart]);
        const TF qtsurf = TF(0.5)*(qtmean [kstart-1] + qtmean[kstart]);
        const TF qhmsurf = TF(0.5)*(qhmmean [kstart-1] + qhmmean[kstart]);

        // Calculate the values at the surface (half level == kstart)
        prefh[kstart] = pbot;
        exh[kstart] = exner(prefh[kstart]);

        Struct_sat_adjust<TF> ssa =
            sat_adjust<TF, sw_satadjust>(thlsurf, qtsurf, prefh[kstart], exh[kstart]);

        thvh[kstart] = virtual_temperature<TF, sw_satadjust>(exh[kstart], thlsurf, qtsurf, ssa.ql, ssa.qi, ssa.t, qhmsurf);
        rhoh[kstart] = pbot / (Rd<TF> * exh[kstart] * thvh[kstart]);

        // Calculate the first full level pressure
        pref[kstart] = prefh[kstart] * std::exp(-grav<TF> * z[kstart] / (Rd<TF> * exh[kstart] * thvh[kstart]));

        for (int k=kstart+1; k<kend+1; ++k)
        {
            // 1. Calculate remaining values (thv and rho) at full-level[k-1]
            ex[k-1]  = exner(pref[k-1]);
            ssa = sat_adjust<TF, sw_satadjust>(thlmean[k-1], qtmean[k-1], pref[k-1], ex[k-1]);
            thv[k-1] = virtual_temperature<TF, sw_satadjust>(ex[k-1], thlmean[k-1], qtmean[k-1], ssa.ql, ssa.qi, ssa.t, qhmmean[k-1]);
            rho[k-1] = pref[k-1] / (Rd<TF> * ex[k-1] * thv[k-1]);

            // 2. Calculate pressure at half-level[k]
            prefh[k] = prefh[k-1] * std::exp(-grav<TF> * dz[k-1] / (Rd<TF> * ex[k-1] * thv[k-1]));
            exh[k] = exner(prefh[k]);

            // 3. Use interpolated conserved quantities to calculate half-level[k] values
            const TF thli = TF(0.5)*(thlmean[k-1] + thlmean[k]);
            const TF qti = TF(0.5)*(qtmean [k-1] + qtmean [k]);
            const TF qhmi = TF(0.5)*(qhmmean [k-1] + qhmmean [k]);

            ssa = sat_adjust<TF, sw_satadjust>(thli, qti, prefh[k], exh[k]);

            thvh[k] = virtual_temperature<TF, sw_satadjust>(exh[k], thli, qti, ssa.ql, ssa.qi, ssa.t, qhmi);
            rhoh[k] = prefh[k] / (Rd<TF> * exh[k] * thvh[k]);

            // 4. Calculate pressure at full-level[k]
            pref[k] = pref[k-1] * std::exp(-grav<TF> * dzh[k] / (Rd<TF> * exh[k] * thvh[k]));
        }

        pref[kstart-1] = TF(2.)*prefh[kstart] - pref[kstart];
    }

    template<typename TF>
    void calc_base_state_no_ql(
            TF* restrict pref,
            TF* restrict prefh,
            TF* restrict rho,
            TF* restrict rhoh,
            TF* restrict thv,
            TF* restrict thvh,
            TF* restrict ex,
            TF* restrict exh,
            TF* restrict thlmean,
            TF* restrict qtmean,
            const TF pbot,
            const int kstart,
            const int kend,
            const TF* restrict z,
            const TF* restrict dz,
            const TF* const dzh)
    {
        const TF thlsurf = TF(0.5)*(thlmean[kstart-1] + thlmean[kstart]);
        const TF qtsurf  = TF(0.5)*(qtmean[kstart-1] + qtmean[kstart]);

        // Calculate the values at the surface (half level == kstart)
        prefh[kstart] = pbot;
        exh[kstart]   = exner(prefh[kstart]);
        thvh[kstart]  = virtual_temperature<TF, Satadjust_type::Disabled>(thlsurf, qtsurf);
        rhoh[kstart]  = pbot / (Rd<TF> * exh[kstart] * thvh[kstart]);

        // Calculate the first full level pressure
        pref[kstart]  = prefh[kstart] * std::exp(-grav<TF> * z[kstart] / (Rd<TF> * exh[kstart] * thvh[kstart]));

        for (int k=kstart+1; k<kend+1; ++k)
        {
            // 1. Calculate remaining values (thv and rho) at full-level[k-1]
            ex[k-1]  = exner(pref[k-1]);
            thv[k-1] = virtual_temperature<TF, Satadjust_type::Disabled>(thlmean[k-1], qtmean[k-1]);
            rho[k-1] = pref[k-1] / (Rd<TF> * ex[k-1] * thv[k-1]);

            // 2. Calculate pressure at half-level[k]
            prefh[k] = prefh[k-1] * std::exp(-grav<TF> * dz[k-1] / (Rd<TF> * ex[k-1] * thv[k-1]));
            exh[k]   = exner(prefh[k]);

            // 3. Use interpolated conserved quantities to calculate half-level[k] values
            const TF thli = TF(0.5)*(thlmean[k-1] + thlmean[k]);
            const TF qti  = TF(0.5)*(qtmean [k-1] + qtmean [k]);

            thvh[k]  = virtual_temperature<TF, Satadjust_type::Disabled>(thli, qti);
            rhoh[k]  = prefh[k] / (Rd<TF> * exh[k] * thvh[k]);

            // 4. Calculate pressure at full-level[k]
            pref[k] = pref[k-1] * std::exp(-grav<TF> * dzh[k] / (Rd<TF> * exh[k] * thvh[k]));
        }

        pref[kstart-1] = TF(2.)*prefh[kstart] - pref[kstart];
    }

    template<typename TF>
    void add_hydrometeor(TF* restrict qhm_total,
                         TF* restrict qhm,
                         const int istart, const int iend,
                         const int jstart, const int jend,
                         const int kstart, const int kend,
                         const int jj, const int kk)
    {
        for (int k = kstart; k < kend; k++)
            for (int j = jstart; j < jend; j++)
                for (int i = istart; i < iend; i++)
                {
                    const int ijk = i + j * jj + k * kk;
                    qhm_total[ijk] += qhm[ijk];
                }
    }
}
#endif
