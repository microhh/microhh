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

#include <iostream>
#include <cmath>
#include <vector>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "diff.h"
#include "stats.h"
#include "cross.h"
#include "thermo.h"
#include "thermo_moist_functions.h"

//  #include "microphys_2mom_warm.h"

#include "tools.h"
#include "field3d_operators.h"

#include "constants.h"
#include "microphys.h"
#include "microphys_nsw6.h"
// Constants, move out later.
using namespace Constants;
using namespace Thermo_moist_functions;
//  using namespace Micro_2mom_warm_constants;
//  using namespace Micro_2mom_warm_functions;

namespace
{
    template<typename TF> constexpr TF qv_min = 1.e-7; // Threshold qv for calculating microphysical terms.
    template<typename TF> constexpr TF ql_min = 1.e-7; // Threshold ql for calculating microphysical terms.
    template<typename TF> constexpr TF qi_min = 1.e-7; // Threshold qi for calculating microphysical terms.
    template<typename TF> constexpr TF qr_min = 1.e-12; // Threshold qr for calculating microphysical terms.
    template<typename TF> constexpr TF qs_min = 1.e-12; // Threshold qs for calculating microphysical terms.
    template<typename TF> constexpr TF qg_min = 1.e-12; // Threshold qg for calculating microphysical terms.
    template<typename TF> constexpr TF q_tiny = 1.e-15; // Threshold qg for calculating microphysical terms.

    template<typename TF> constexpr TF pi = M_PI; // Pi constant.
    template<typename TF> constexpr TF pi_2 = M_PI*M_PI; // Pi constant squared.

    template<typename TF> constexpr TF rho_w = 1.e3; // Density of water.
    template<typename TF> constexpr TF rho_s = 1.e2; // Density of snow.
    template<typename TF> constexpr TF rho_g = 4.e2; // Density of graupel.

    template<typename TF> constexpr TF N_0r = 8.e6; // Intercept parameter rain (m-4).
    template<typename TF> constexpr TF N_0s = 3.e6; // Intercept parameter snow (m-4).
    template<typename TF> constexpr TF N_0g = 4.e6; // Intercept parameter graupel (m-4).

    template<typename TF> constexpr TF a_r = M_PI*rho_w<TF>/6.; // Empirical constant for m_r.
    template<typename TF> constexpr TF a_s = M_PI*rho_s<TF>/6.; // Empirical constant for m_s.
    template<typename TF> constexpr TF a_g = M_PI*rho_g<TF>/6.; // Empirical constant for m_g.

    template<typename TF> constexpr TF b_r = 3.; // Empirical constant for m_r.
    template<typename TF> constexpr TF b_s = 3.; // Empirical constant for m_s.
    template<typename TF> constexpr TF b_g = 3.; // Empirical constant for m_g.

    template<typename TF> constexpr TF c_r = 130.; // Empirical constant for v_r (wrong value in Tomita's paper).
    template<typename TF> constexpr TF c_s = 4.84; // Empirical constant for v_s.
    template<typename TF> constexpr TF c_g = 82.5; // Empirical constant for v_g.

    template<typename TF> constexpr TF d_r = 0.5;  // Empirical constant for v_r.
    template<typename TF> constexpr TF d_s = 0.25; // Empirical constant for v_s.
    template<typename TF> constexpr TF d_g = 0.5;  // Empirical constant for v_g.

    template<typename TF> constexpr TF C_i = 2006.; // Specific heat of solid water.
    template<typename TF> constexpr TF C_l = 4218.; // Specific heat of liquid water.

    template<typename TF> constexpr TF f_1r = 0.78; // First coefficient of ventilation factor for rain.
    template<typename TF> constexpr TF f_1s = 0.65; // First coefficient of ventilation factor for snow.
    template<typename TF> constexpr TF f_1g = 0.78; // First coefficient of ventilation factor for graupel.

    template<typename TF> constexpr TF f_2r = 0.27; // First coefficient of ventilation factor for rain.
    template<typename TF> constexpr TF f_2s = 0.39; // First coefficient of ventilation factor for snow.
    template<typename TF> constexpr TF f_2g = 0.27; // First coefficient of ventilation factor for graupel.

    template<typename TF> constexpr TF E_ri = 1.;  // Collection efficiency of ice for rain.
    template<typename TF> constexpr TF E_rw = 1.;  // Collection efficiency of rain for cloud water.
    template<typename TF> constexpr TF E_sw = 1.;  // Collection efficiency of snow for cloud water.
    template<typename TF> constexpr TF E_gw = 1.;  // Collection efficiency of graupel for cloud water.
    template<typename TF> constexpr TF E_gi = 0.1; // Collection efficiency of graupel for cloud ice.
    template<typename TF> constexpr TF E_sr = 1.;  // Collection efficiency of snow for rain.
    // template<typename TF> constexpr TF E_gr = 1.;  // Collection efficiency of graupel for rain.
    template<typename TF> constexpr TF E_gr = 0.1; // CvH: I reduced the value as the current convection is too active in cold pools, leading to conversion from rain to graupel that extends too far to the surface.

    template<typename TF> constexpr TF K_a = 2.43e-2;  // Thermal diffusion coefficient of air.
    template<typename TF> constexpr TF K_d = 2.26e-5;  // Diffusion coefficient of water vapor in air.

    template<typename TF> constexpr TF M_i = 4.19e-13; // Mass of one cloud ice particle.

    template<typename TF> constexpr TF beta_saut = 6.e-3; // Values from SCALE-LES model.
    template<typename TF> constexpr TF beta_gaut = 0.e-3; // Values from SCALE-LES model.

    template<typename TF> constexpr TF gamma_sacr = 25.e-3;
    template<typename TF> constexpr TF gamma_saut = 60.e-3; // Tomita's code in SCALE is different than paper (0.025).
    template<typename TF> constexpr TF gamma_gacs = 90.e-3;
    template<typename TF> constexpr TF gamma_gaut = 90.e-3;

    template<typename TF> constexpr TF nu = 1.5e-5; // Kinematic viscosity of air

    template<typename TF> CUDA_MACRO
    inline TF minmod(const TF a, const TF b)
    {
        #ifdef __CUDACC__
        return copysign(TF(1.), a) * fmax(TF(0.), fmin(fabs(a), TF(copysign(TF(1.), a))*b));
        #else
        return copysign(TF(1.), a) * std::max(TF(0.), std::min(std::abs(a), TF(copysign(TF(1.), a))*b));
        #endif
    }
    // Taken from include/microphys_2mom_warm.h , can likely delete the cxx version
}
// Above are all constants imported from microphys.cxx

namespace
{
    using namespace Constants;
    using namespace Thermo_moist_functions;
    using namespace Fast_math;

    // Compute all microphysical tendencies.
    template<typename TF> __global__
    void conversion_g(
            TF* const __restrict__ qrt, TF* const __restrict__ qst, TF* const __restrict__ qgt,
            TF* const __restrict__ qtt, TF* const __restrict__ thlt,
            const TF* const __restrict__ qr, const TF* const __restrict__ qs, const TF* const __restrict__ qg,
            const TF* const __restrict__ qt, const TF* const __restrict__ thl,
            const TF* const __restrict__ ql, const TF* const __restrict__ qi,
            const TF* const __restrict__ rho, const TF* const __restrict__ exner, const TF* const __restrict__ p,
            const TF* const __restrict__ dzi, const TF* const __restrict__ dzhi,
            const TF N_d, const TF dt,
            const int istart, const int jstart, const int kstart,
            const int iend, const int jend, const int kend,
            const int jj, const int kk, TF* const __restrict__ temp_array)
    {
        // Tomita Eq. 51. N_d is converted from SI units (m-3 instead of cm-3).
        const TF D_d = TF(0.146) - TF(5.964e-2)* log((N_d*TF(1.e-6)) / TF(2.e3));
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        //const int k = blockIdx.z + kstart;

        for (int k=kstart; k<kend; ++k)
        {
            const TF rho0_rho_sqrt = pow(rho[kstart]/rho[k], TF(0.5));
            // Setup the relevant variables for each z slice

            // Part of Tomita Eq. 29
            const TF fac_iacr =
            pi_2<TF> * E_ri<TF> * N_0r<TF> * c_r<TF> * rho_w<TF> * tgamma(TF(6.) + d_r<TF>)
            / (TF(24.) * M_i<TF>)
            * rho0_rho_sqrt;

            // Part of Tomita Eq. 32
            const TF fac_raci =
            pi<TF> * E_ri<TF> * N_0r<TF> * c_r<TF> * tgamma(TF(3.) + d_r<TF>)
            / TF(4.)
            * rho0_rho_sqrt;

            // Part of Tomita Eq. 34
            const TF fac_racw =
            pi<TF> * E_rw<TF> * N_0r<TF> * c_r<TF> * tgamma(TF(3.) + d_r<TF>)
            / TF(4.)
            * rho0_rho_sqrt;

            // Part of Tomita Eq. 35
            const TF fac_sacw =
            pi<TF> * E_sw<TF> * N_0s<TF> * c_s<TF> * tgamma(TF(3.) + d_s<TF>)
            / TF(4.)
            * rho0_rho_sqrt;

            // Part of Tomita Eq. 36 (E_si is temperature dependent and missing therefore here).
            const TF fac_saci =
            pi<TF> * N_0s<TF> * c_s<TF> * tgamma(TF(3.) + d_s<TF>)
            / TF(4.)
            * rho0_rho_sqrt;

            // Part of Tomita Eq. 37
            const TF fac_gacw =
            pi<TF> * E_gw<TF> * N_0g<TF> * c_g<TF> * tgamma(TF(3.) + d_g<TF>)
            / TF(4.)
            * rho0_rho_sqrt;

            // Part of Tomita Eq. 38
            const TF fac_gaci =
            pi<TF> * E_gi<TF> * N_0g<TF> * c_g<TF> * tgamma(TF(3.) + d_g<TF>)
            / TF(4.)
            * rho0_rho_sqrt;

            // end z slice var setup

            if (i < iend && j < jend)
            {
                const int ijk = i + j*jj + k*kk;

                // Compute the T out of the known values of ql and qi, this saves memory and sat_adjust.
                const TF T = exner[k]*thl[ijk] + Lv<TF>/cp<TF>*ql[ijk] + Ls<TF>/cp<TF>*qi[ijk];
                const TF qv = qt[ijk] - ql[ijk] - qi[ijk];

                // Flag the sign of the absolute temperature.
                const TF T_pos = TF(T >= T0<TF>);
                const TF T_neg = TF(1.) - T_pos;

                // Check which species are present.
                const bool has_vapor   = (qv      > qv_min<TF>);
                const bool has_liq     = (ql[ijk] > ql_min<TF>);
                const bool has_ice     = (qi[ijk] > qi_min<TF>);
                const bool has_rain    = (qr[ijk] > qr_min<TF>);
                const bool has_snow    = (qs[ijk] > qs_min<TF>);
                const bool has_graupel = (qg[ijk] > qg_min<TF>);

                if (! (has_liq || has_ice || has_rain || has_snow || has_graupel) )
                    continue;

                // Tomita Eq. 27
                const TF lambda_r = pow(
                        a_r<TF> * N_0r<TF> * tgamma(b_r<TF> + TF(1.))
                        / (rho[k] * (qr[ijk] + q_tiny<TF>)),
                        TF(1.) / (b_r<TF> + TF(1.)) );

                const TF lambda_s = pow(
                        a_s<TF> * N_0s<TF> * tgamma(b_s<TF> + TF(1.))
                        / (rho[k] * (qs[ijk] + q_tiny<TF>)),
                        TF(1.) / (b_s<TF> + TF(1.)) );

                const TF lambda_g = pow(
                        a_g<TF> * N_0g<TF> * tgamma(b_g<TF> + TF(1.))
                        / (rho[k] * (qg[ijk] + q_tiny<TF>)),
                        TF(1.) / (b_g<TF> + TF(1.)) );

                // Tomita Eq. 28
                const TF V_Tr = !(has_rain) ? TF(0.) :
                    c_r<TF> * rho0_rho_sqrt
                    * tgamma(b_r<TF> + d_r<TF> + TF(1.)) / tgamma(b_r<TF> + TF(1.))
                    * pow(lambda_r, -d_r<TF>);

                const TF V_Ts = !(has_snow) ? TF(0.) :
                    c_s<TF> * rho0_rho_sqrt
                    * tgamma(b_s<TF> + d_s<TF> + TF(1.)) / tgamma(b_s<TF> + TF(1.))
                    * pow(lambda_s, -d_s<TF>);

                const TF V_Tg = !(has_graupel) ? TF(0.) :
                    c_g<TF> * rho0_rho_sqrt
                    * tgamma(b_g<TF> + d_g<TF> + TF(1.)) / tgamma(b_g<TF> + TF(1.))
                    * pow(lambda_g, -d_g<TF>);

                // ACCRETION
                // Tomita Eq. 29
                const TF P_iacr = !(has_rain && has_ice) ? TF(0.) :
                    fac_iacr / pow(lambda_r, TF(6.) + d_r<TF>) * qi[ijk];

                // Tomita Eq. 30
                const TF delta_1 = TF(qr[ijk] >= TF(1.e-4));

                // Tomita Eq. 31
                TF P_iacr_s = (TF(1.) - delta_1) * P_iacr;
                TF P_iacr_g = delta_1 * P_iacr;

                // Tomita Eq. 32
                const TF P_raci = !(has_rain && has_ice) ? TF(0.) :
                    fac_raci / pow(lambda_r, TF(3.) + d_r<TF>) * qi[ijk];

                // Tomita Eq. 33
                TF P_raci_s = (TF(1.) - delta_1) * P_raci;
                TF P_raci_g = delta_1 * P_raci;

                // Tomita Eq. 34, 35
                TF P_racw = !(has_liq && has_rain) ? TF(0.) :
                    fac_racw / pow(lambda_r, TF(3.) + d_r<TF>) * ql[ijk];
                TF P_sacw = !(has_liq && has_snow) ? TF(0.) :
                    fac_sacw / pow(lambda_s, TF(3.) + d_s<TF>) * ql[ijk];

                // Tomita Eq. 39
                const TF E_si = exp(gamma_sacr<TF> * (T - T0<TF>));

                // Tomita Eq. 36 - 38
                TF P_saci = !(has_snow && has_ice) ? TF(0.) :
                    fac_saci * E_si / pow(lambda_s, TF(3.) + d_s<TF>) * qi[ijk];
                TF P_gacw = !(has_graupel && has_liq) ? TF(0.) :
                    fac_gacw / pow(lambda_g, TF(3.) + d_g<TF>) * ql[ijk];
                TF P_gaci = !(has_graupel && has_ice) ? TF(0.) :
                    fac_gaci / pow(lambda_g, TF(3.) + d_g<TF>) * qi[ijk];

                // Accretion of falling hydrometeors.
                // Tomita Eq. 42
                const TF delta_2 = TF(1.) - TF( (qr[ijk] >= TF(1.e-4)) || (qs[ijk] >= TF(1.e-4)) );

                // Tomita Eq. 41
                TF P_racs =
                !(has_rain && has_snow) ? TF(0.) :
                    (TF(1.) - delta_2)
                    * pi<TF> * a_s<TF> * abs(V_Tr - V_Ts) * E_sr<TF> * N_0s<TF> * N_0r<TF> / (TF(4.)*rho[k])
                    * (          tgamma(b_s<TF> + TF(3.)) * tgamma(TF(1.)) / ( pow(lambda_s, b_s<TF> + TF(3.)) * lambda_r )
                    + TF(2.) * tgamma(b_s<TF> + TF(2.)) * tgamma(TF(2.)) / ( pow(lambda_s, b_s<TF> + TF(2.)) * pow2(lambda_r) )
                    +          tgamma(b_s<TF> + TF(1.)) * tgamma(TF(3.)) / ( pow(lambda_s, b_s<TF> + TF(1.)) * pow3(lambda_r) ) );

                // Tomita Eq. 44
                const TF P_sacr =
                !(has_snow && has_rain) ? TF(0.) :
                    pi<TF> * a_r<TF> * abs(V_Ts - V_Tr) * E_sr<TF> * N_0r<TF> * N_0s<TF> / (TF(4.)*rho[k])
                    * (          tgamma(b_r<TF> + TF(1.)) * tgamma(TF(3.)) / ( pow(lambda_r, b_r<TF> + TF(1.)) * pow3(lambda_s) )
                    + TF(2.) * tgamma(b_r<TF> + TF(2.)) * tgamma(TF(2.)) / ( pow(lambda_r, b_r<TF> + TF(2.)) * pow2(lambda_s) )
                    +          tgamma(b_r<TF> + TF(3.)) * tgamma(TF(1.)) / ( pow(lambda_r, b_r<TF> + TF(3.)) * lambda_s ) );

                // Tomita Eq. 43
                TF P_sacr_g = (TF(1.) - delta_2) * P_sacr;
                TF P_sacr_s = delta_2 * P_sacr;

                // Tomita Eq. 49
                const TF E_gs = fmin( TF(1.), exp(gamma_gacs<TF> * (T - T0<TF>)) );

                // Tomita Eq. 47
                TF P_gacr =
                !(has_graupel && has_rain) ? TF(0.) :
                    pi<TF> * a_r<TF> * abs(V_Tg - V_Tr) * E_gr<TF> * N_0g<TF> * N_0r<TF> / (TF(4.)*rho[k])
                    * (          tgamma(b_r<TF> + TF(1.)) * tgamma(TF(3.)) / ( pow(lambda_r, b_r<TF> + TF(1.)) * pow3(lambda_g) )
                    + TF(2.) * tgamma(b_r<TF> + TF(2.)) * tgamma(TF(2.)) / ( pow(lambda_r, b_r<TF> + TF(2.)) * pow2(lambda_g) )
                    +          tgamma(b_r<TF> + TF(3.)) * tgamma(TF(1.)) / ( pow(lambda_r, b_r<TF> + TF(3.)) * lambda_g ) );

                // Tomita Eq. 48
                TF P_gacs =
                 !(has_graupel && has_snow) ? TF(0.) :
                    pi<TF> * a_s<TF> * abs(V_Tg - V_Ts) * E_gs * N_0g<TF> * N_0s<TF> / (TF(4.)*rho[k])
                    * (          tgamma(b_s<TF> + TF(1.)) * tgamma(TF(3.)) / ( pow(lambda_s, b_s<TF> + TF(1.)) * pow3(lambda_g) )
                    + TF(2.) * tgamma(b_s<TF> + TF(2.)) * tgamma(TF(2.)) / ( pow(lambda_s, b_s<TF> + TF(2.)) * pow2(lambda_g) )
                    +          tgamma(b_s<TF> + TF(3.)) * tgamma(TF(1.)) / ( pow(lambda_s, b_s<TF> + TF(3.)) * lambda_g ) );

                // AUTOCONVERSION.
                constexpr TF q_icrt = TF(0.);
                constexpr TF q_scrt = TF(6.e-4);

                // Tomita Eq. 53
                const TF beta_1 = fmin( beta_saut<TF>, beta_saut<TF>*exp(gamma_saut<TF> * (T - T0<TF>)) );

                // Tomita Eq. 54
                const TF beta_2 = fmin( beta_gaut<TF>, beta_gaut<TF>*exp(gamma_gaut<TF> * (T - T0<TF>)) );

                // Tomita Eq. 50. Our N_d is SI units, so conversion is applied.
                TF P_raut = !(has_liq) ? TF(0.) :
                    TF(16.7)/rho[k] * pow2(rho[k]*ql[ijk]) / (TF(5.) + TF(3.66e-2) * TF(1.e-6)*N_d / (D_d*rho[k]*ql[ijk]));

                // Tomita Eq. 52
                TF P_saut = !(has_ice) ? TF(0.) :
                    fmax(beta_1*(qi[ijk] - q_icrt), TF(0.));

                // Tomita Eq. 54
                TF P_gaut = !(has_snow) ? TF(0.) :
                    fmax(beta_2*(qs[ijk] - q_scrt), TF(0.));

                // PHASE CHANGES.
                // Tomita Eq. 57
                const TF G_w = TF(1.) / (
                    Lv<TF> / (K_a<TF> * T) * (Lv<TF> / (Rv<TF> * T) - TF(1.))
                    + Rv<TF>*T / (K_d<TF> * esat_liq(T)) );

                    /// CHANGE THIS BACK


                // Tomita Eq. 62
                const TF G_i = TF(1.) / (
                    Ls<TF> / (K_a<TF> * T) * (Ls<TF> / (Rv<TF> * T) - TF(1.))
                    + Rv<TF>*T / (K_d<TF> * esat_ice(T)) );

                const TF S_w = (qt[ijk] - ql[ijk] - qi[ijk]) / qsat_liq(p[k], T);
                const TF S_i = (qt[ijk] - ql[ijk] - qi[ijk]) / qsat_ice(p[k], T);

                // Tomita Eq. 63
                const TF delta_3 = TF(S_i <= TF(1.)); // Subsaturated, then delta_3 = 1.

                // Tomita Eq. 59
                TF P_revp = !(has_rain) ? TF(0.) : // still requires bounding ?
                - TF(2.)*pi<TF> * N_0r<TF> * (fmin(S_w, TF(1.)) - TF(1.)) * G_w / rho[k]
                * ( f_1r<TF> * tgamma(TF(2.)) / pow(lambda_r,TF(2.))
                + f_2r<TF> * sqrt(c_r<TF> * rho0_rho_sqrt / nu<TF>)
                * tgamma( TF(0.5) * (TF(5.) + d_r<TF>) )
                / pow(lambda_r, TF(0.5) * (TF(5.) + d_r<TF>)) );

                //  fmin( TF(1.e-10), !(has_rain) ? TF(0.) : // still requires bounding by zero!
                //      - TF(2.)*pi<TF> * N_0r<TF> * (fmin(S_w, TF(1.)) - TF(1.)) * G_w / rho[k]
                //      * ( f_1r<TF> * tgamma(TF(2.)) / pow(lambda_r,TF(2.))
                //      + f_2r<TF> * sqrt(c_r<TF> * rho0_rho_sqrt / nu<TF>)
                //      * tgamma( TF(0.5) * (TF(5.) + d_r<TF>) )
                //      / pow(lambda_r, TF(0.5) * (TF(5.) + d_r<TF>)) ));

                temp_array[ijk] = P_revp; //Rv<TF>; // G_w * rho[k]; //

                //******************************************************************************

                // Tomita Eq. 60. Negative for sublimation, positive for deposition.
                const TF P_sdep_ssub =
                    TF(2.)*pi<TF> * N_0s<TF> * (S_i - TF(1.)) * G_i / rho[k]
                    * ( f_1s<TF> * tgamma(TF(2.)) / pow2(lambda_s)
                      + f_2s<TF> * sqrt(c_s<TF> * rho0_rho_sqrt / nu<TF>)
                      * tgamma( TF(0.5) * (TF(5.) + d_s<TF>) )
                      / pow(lambda_s, TF(0.5) * (TF(5.) + d_s<TF>)) );

                // Tomita Eq. 51
                const TF P_gdep_gsub =
                    TF(2.)*pi<TF> * N_0g<TF> * (S_i - TF(1.)) * G_i / rho[k]
                    * ( f_1g<TF> * tgamma(TF(2.)) / pow2(lambda_g)
                      + f_2g<TF> * sqrt(c_g<TF> * rho0_rho_sqrt / nu<TF>)
                      * tgamma( TF(0.5) * (TF(5.) + d_g<TF>) )
                      / pow(lambda_g, TF(0.5) * (TF(5.) + d_g<TF>)) );

                // Tomita Eq. 64
                TF P_sdep = !(has_vapor) ? TF(0.) :
                    (TF(1.) - delta_3) * P_sdep_ssub;
                TF P_gdep = !(has_vapor) ? TF(0.) :
                    (TF(1.) - delta_3) * P_gdep_gsub;

                // Tomita Eq. 65
                // CvH: I swapped the sign with respect to Tomita, the term should be positive.
                TF P_ssub = !(has_snow) ? TF(0.) :
                    - delta_3 * P_sdep_ssub;
                TF P_gsub = !(has_graupel) ? TF(0.) :
                    - delta_3 * P_gdep_gsub;

                // Freezing and melting
                // Tomita Eq. 67, 68 combined.
                TF P_smlt = !(has_snow) ? TF(0.) :
                    TF(2.)*pi<TF> * K_a<TF> * (T - T0<TF>) * N_0s<TF> / (rho[k]*Lf<TF>)
                    * ( f_1s<TF> * tgamma(TF(2.)) / pow2(lambda_s)
                      + f_2s<TF> * sqrt(c_s<TF> * rho0_rho_sqrt / nu<TF>)
                      * tgamma( TF(0.5) * (TF(5.) + d_s<TF>) )
                      / pow(lambda_s, TF(0.5) * (TF(5.) + d_s<TF>)) )
                    + C_l<TF> * (T - T0<TF>) / Lf<TF> * (P_sacw + P_sacr);

                // Tomita Eq. 69
                TF P_gmlt = !(has_graupel) ? TF(0.) :
                    TF(2.)*pi<TF> * K_a<TF> * (T - T0<TF>) * N_0g<TF> / (rho[k]*Lf<TF>)
                    * ( f_1g<TF> * tgamma(TF(2.)) / pow2(lambda_g)
                      + f_2g<TF> * sqrt(c_g<TF> * rho0_rho_sqrt / nu<TF>)
                      * tgamma( TF(0.5) * (TF(5.) + d_g<TF>) )
                      / pow(lambda_g, TF(0.5) * (TF(5.) + d_g<TF>)) )
                    + C_l<TF> * (T - T0<TF>) / Lf<TF> * (P_gacw + P_gacr);

                // Tomita Eq. 70
                constexpr TF A_prime = TF(0.66);
                constexpr TF B_prime = TF(100.);

                TF P_gfrz = !(has_rain) ? TF(0.) :
                    TF(20.) * pi_2<TF> * B_prime * N_0r<TF> * rho_w<TF> / rho[k]
                    * (exp(A_prime * (T0<TF> - T)) - TF(1.)) / pow7(lambda_r);

                // COMPUTE THE TENDENCIES.
                // Limit the production terms to avoid instability.
                auto limit_tend = [&](TF& tend, const TF tend_limit)
                {
                    tend = fmax(TF(0.), fmin(tend, tend_limit));
                };

                const TF dqv_dt_max = qv      / dt;
                const TF dqi_dt_max = qi[ijk] / dt;
                const TF dql_dt_max = ql[ijk] / dt;
                const TF dqr_dt_max = qr[ijk] / dt;
                const TF dqs_dt_max = qs[ijk] / dt;
                const TF dqg_dt_max = qg[ijk] / dt;

                // Limit on the availability of the source.
                // Limit accretion terms.
                limit_tend(P_iacr_s, dqr_dt_max);
                limit_tend(P_iacr_g, dqr_dt_max);
                limit_tend(P_raci_s, dqi_dt_max);
                limit_tend(P_raci_g, dqi_dt_max);
                limit_tend(P_racw  , dql_dt_max);
                limit_tend(P_sacw  , dql_dt_max);
                limit_tend(P_saci  , dqi_dt_max);
                limit_tend(P_gacw  , dql_dt_max);
                limit_tend(P_gaci  , dqi_dt_max);
                limit_tend(P_racs  , dqs_dt_max);
                limit_tend(P_sacr_s, dqr_dt_max);
                limit_tend(P_sacr_g, dqr_dt_max);
                limit_tend(P_gacr  , dqr_dt_max);
                limit_tend(P_gacs  , dqs_dt_max);

                // Limit autoconversion terms.
                limit_tend(P_raut, dql_dt_max);
                limit_tend(P_saut, dqi_dt_max);
                limit_tend(P_gaut, dqs_dt_max);

                // Limit phase changes.
                limit_tend(P_revp, dqr_dt_max);
                limit_tend(P_sdep, dqv_dt_max);
                limit_tend(P_ssub, dqs_dt_max);
                limit_tend(P_gdep, dqv_dt_max);
                limit_tend(P_gsub, dqg_dt_max);
                limit_tend(P_smlt, dqs_dt_max);
                limit_tend(P_gmlt, dqg_dt_max);
                limit_tend(P_gfrz, dqr_dt_max);

                TF vapor_to_snow = P_sdep;
                TF vapor_to_graupel = P_gdep;

                TF cloud_to_rain = P_racw + P_sacw * T_pos + P_raut;
                TF cloud_to_graupel = P_gacw;
                TF cloud_to_snow = P_sacw * T_neg;

                TF rain_to_vapor = P_revp;
                TF rain_to_graupel = P_gacr + P_iacr_g + P_sacr_g * T_neg + P_gfrz * T_neg;
                TF rain_to_snow = P_sacr_s * T_neg + P_iacr_s;

                TF ice_to_snow = P_raci_s + P_saci + P_saut;
                TF ice_to_graupel = P_raci_g + P_gaci;

                TF snow_to_graupel = P_gacs + P_racs + P_gaut;
                TF snow_to_rain = P_smlt;
                TF snow_to_vapor = P_ssub;

                TF graupel_to_rain = P_gmlt * T_pos;
                TF graupel_to_vapor = P_gsub;

                const TF dqv_dt =
                    - vapor_to_snow - vapor_to_graupel;

                const TF dql_dt =
                    - cloud_to_rain - cloud_to_graupel - cloud_to_snow;

                const TF dqi_dt =
                    - ice_to_snow - ice_to_graupel;

                const TF dqr_dt =
                    + cloud_to_rain + snow_to_rain + graupel_to_rain
                    - rain_to_vapor - rain_to_graupel - rain_to_snow;

                const TF dqs_dt =
                    + cloud_to_snow + ice_to_snow + vapor_to_snow
                    - snow_to_graupel - snow_to_vapor - snow_to_rain;

                const TF dqg_dt =
                    + cloud_to_graupel + rain_to_graupel + ice_to_graupel
                    + vapor_to_graupel + snow_to_graupel
                    - graupel_to_rain - graupel_to_vapor;

                // Limit the production terms to avoid instability.
                auto limit_factor = [](const TF tend, const TF tend_limit)
                {
                    return (tend < TF(0.)) ? fmin(-tend_limit/tend, TF(1.)) : TF(1.);
                };

                const TF dqv_dt_fac = limit_factor(dqv_dt, dqv_dt_max);
                const TF dql_dt_fac = limit_factor(dql_dt, dql_dt_max);
                const TF dqi_dt_fac = limit_factor(dqi_dt, dqi_dt_max);
                const TF dqr_dt_fac = limit_factor(dqr_dt, dqr_dt_max);
                const TF dqs_dt_fac = limit_factor(dqs_dt, dqs_dt_max);
                const TF dqg_dt_fac = limit_factor(dqg_dt, dqg_dt_max);

                vapor_to_snow    *= dqv_dt_fac * dqs_dt_fac;
                vapor_to_graupel *= dqv_dt_fac * dqg_dt_fac;

                cloud_to_rain    *= dql_dt_fac * dqr_dt_fac;
                cloud_to_graupel *= dql_dt_fac * dqg_dt_fac;
                cloud_to_snow    *= dql_dt_fac * dqs_dt_fac;

                rain_to_vapor    *= dqr_dt_fac * dqv_dt_fac;
                rain_to_graupel  *= dqr_dt_fac * dqg_dt_fac;
                rain_to_snow     *= dqr_dt_fac * dqs_dt_fac;

                ice_to_snow      *= dqi_dt_fac * dqs_dt_fac;
                ice_to_graupel   *= dqi_dt_fac * dqg_dt_fac;

                snow_to_graupel  *= dqs_dt_fac * dqg_dt_fac;
                snow_to_vapor    *= dqs_dt_fac * dqv_dt_fac;
                snow_to_rain     *= dqs_dt_fac * dqr_dt_fac;

                graupel_to_rain  *= dqg_dt_fac * dqr_dt_fac;
                graupel_to_vapor *= dqg_dt_fac * dqv_dt_fac;

                // Loss from cloud.

                qtt[ijk] -= cloud_to_rain;
                qrt[ijk] += cloud_to_rain;
                thlt[ijk] += Lv<TF> / (cp<TF> * exner[k]) * cloud_to_rain;

                qtt[ijk] -= cloud_to_graupel;
                qgt[ijk] += cloud_to_graupel;
                thlt[ijk] += Ls<TF> / (cp<TF> * exner[k]) * cloud_to_graupel;

                qtt[ijk] -= cloud_to_snow;
                qst[ijk] += cloud_to_snow;
                thlt[ijk] += Ls<TF> / (cp<TF> * exner[k]) * cloud_to_snow;

                // Loss from rain.
                qrt[ijk] -= rain_to_vapor;
                qtt[ijk] += rain_to_vapor;
                thlt[ijk] -= Lv<TF> / (cp<TF> * exner[k]) * rain_to_vapor;

                qrt[ijk] -= rain_to_graupel;
                qgt[ijk] += rain_to_graupel;
                thlt[ijk] += Lf<TF> / (cp<TF> * exner[k]) * rain_to_graupel;

                qrt[ijk] -= rain_to_snow;
                qst[ijk] += rain_to_snow;
                thlt[ijk] += Lf<TF> / (cp<TF> * exner[k]) * rain_to_snow;

                   //  qrt[ijk] = fmax(TF(0.0), qrt[ijk]);
            }
        } // end for loop over k
    } // end conversion function

    template<typename TF> __global__
    void sedimentation_ss08_g(
            TF* const __restrict__ qct, // i have deleted rc_bot !
            TF* const __restrict__ w_qc, TF* const __restrict__ c_qc,
            TF* const __restrict__ slope_qc, TF* const __restrict__ flux_qc,
            const TF* const __restrict__ qc,
            const TF* const __restrict__ rho,
            const TF* const __restrict__ dzi, const TF* const __restrict__ dz,
            const double dt,
            const TF a_c, const TF b_c, const TF c_c, const TF d_c, const TF N_0c,
            const TF qc_min,
            const int istart, const int jstart, const int kstart,
            const int iend, const int jend, const int kend,
            const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        //const int k = blockIdx.z + kstart;

        constexpr TF V_Tmin = TF(0.1);
        constexpr TF V_Tmax = TF(10.);

        // 1. Calculate sedimentation velocity at cell center
        for (int k=kstart; k<kend; ++k)
        {
            const TF rho0_rho_sqrt = sqrt(rho[kstart]/rho[k]);

            if (j<jend)
                if (i<iend)
                {
                    const int ijk = i + j*jj + k*kk;

                    if (qc[ijk] > qc_min)
                    {
                        const TF lambda_c = pow(
                            a_c * N_0c * tgamma(b_c + TF(1.))
                            / (rho[k] * qc[ijk]),
                            TF(1.) / (b_c + TF(1.)) );

                        TF V_T =
                            c_c * rho0_rho_sqrt
                            * tgamma(b_c + d_c + TF(1.)) / tgamma(b_c + TF(1.))
                            * pow(lambda_c, -d_c);

                        // Constrain the terminal velocity between 0.1 and 10.
                        V_T = max(V_Tmin, min(V_T, V_Tmax));

                        w_qc[ijk] = V_T;
                    }
                    else
                        w_qc[ijk] = TF(0.);
                }
        }

        // 1.1 Set one ghost cell to zero
        if (j<jend)
            if (i<iend)
            {
                const int ijk_bot = i + j*jj + (kstart-1)*kk;
                const int ijk_top = i + j*jj + (kend    )*kk;
                w_qc[ijk_bot] = w_qc[ijk_bot+kk];
                w_qc[ijk_top] = TF(0.);
            }

        // 2. Calculate CFL number using interpolated sedimentation velocity
        for (int k=kstart; k<kend; ++k)
            if (j<jend)
                if (i<iend)
                {
                    const int ijk = i + j*jj + k*kk;
                    c_qc[ijk] = TF(0.25) * (w_qc[ijk-kk] + TF(2.)*w_qc[ijk] + w_qc[ijk+kk]) * dzi[k] * dt;
                }

        // 3. Calculate slopes
        for (int k=kstart; k<kend; ++k)
            if (j<jend)
                if (i<iend)
                {
                    const int ijk = i + j*jj + k*kk;
                    slope_qc[ijk] = minmod(qc[ijk]-qc[ijk-kk], qc[ijk+kk]-qc[ijk]);
                }

        // Calculate flux
        // Set the fluxes at the top of the domain (kend) to zero
        if (j<jend)
            if (i<iend)
            {
                const int ijk = i + j*jj + kend*kk;
                flux_qc[ijk] = TF(0.);
            }

        for (int k=kend-1; k>kstart-1; --k)
            if (j<jend)
                if (i<iend)
                {
                    const int ijk = i + j*jj + k*kk;

                    // q_rain
                    int kc = k; // current grid level
                    TF ftot = TF(0.); // cumulative 'flux' (kg m-2)
                    TF dzz = TF(0.); // distance from zh[k]
                    TF cc = fmin(TF(1.), c_qc[ijk]);
                    while (cc > 0 && kc < kend)
                    {
                        const int ijkc = i + j*jj+ kc*kk;

                        ftot += rho[kc] * (qc[ijkc] + TF(0.5) * slope_qc[ijkc] * (TF(1.)-cc)) * cc * dz[kc];

                        dzz += dz[kc];
                        kc += 1;
                        cc = fmin(TF(1.), c_qc[ijkc] - dzz*dzi[kc]);
                    }

                    // Given flux at top, limit bottom flux such that the total rain content stays >= 0.
                    ftot = fmin(ftot, rho[k] * dz[k] * qc[ijk] - flux_qc[ijk+kk] * TF(dt));
                    flux_qc[ijk] = -ftot / dt;
                }

        // Calculate tendency
        for (int k=kstart; k<kend; ++k)
            if (j<jend)
                if (i<iend)
                {
                    const int ijk = i + j*jj + k*kk;
                    qct[ijk] += -(flux_qc[ijk+kk] - flux_qc[ijk]) / rho[k] * dzi[k];
                }

        // Store surface sedimentation flux
        // Sedimentation flux is already multiplied with density (see flux div. calculation), so
        // the resulting flux is in kg m-2 s-1, with rho_water = 1000 kg/m3 this equals a rain rate in mm s-1

        // THIS PART DOES NOT MATTER RIGHT NOW...JUST CALCULATES ACCUMULATION ON GROUND (BOTTOM LAYER):
        // for (int j=jstart; j<jend; ++j)
        //     for (int i=istart; i<iend; ++i)
        //     {
        //         const int ij  = i + j*jj;
        //         const int ijk = i + j*jj + kstart*kk;

        //         rc_bot[ij] = -flux_qc[ijk];
        //     }
   }  // end sedimentation function
} // end namespace

#ifdef USECUDA
template<typename TF>
void Microphys_nsw6<TF>::exec(Thermo<TF>& thermo, const double dt, Stats<TF>& stats){

    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax+1);
    dim3 blockGPU(blocki, blockj, 1);

    dim3 grid2dGPU (gridi, gridj);
    dim3 block2dGPU(blocki, blockj);

    auto ql = fields.get_tmp_g();
    auto qi = fields.get_tmp_g();

    thermo.get_thermo_field_g(*ql, "ql", false);
    thermo.get_thermo_field_g(*qi, "qi", false);

    TF* p     = thermo.get_basestate_fld_g("pref");
    TF* exner = thermo.get_basestate_fld_g("exner");

    auto temp_g = fields.get_tmp_g();
    auto temp = fields.get_tmp();

    conversion_g<<<grid2dGPU, block2dGPU>>>(
           fields.st.at("qr")->fld_g, fields.st.at("qs")->fld_g, fields.st.at("qg")->fld_g,
           fields.st.at("qt")->fld_g, fields.st.at("thl")->fld_g,
           fields.sp.at("qr")->fld_g, fields.sp.at("qs")->fld_g, fields.sp.at("qg")->fld_g,
           fields.sp.at("qt")->fld_g, fields.sp.at("thl")->fld_g,
           ql->fld_g, qi->fld_g,
           fields.rhoref_g, exner, p,
           gd.dzi_g, gd.dzhi_g,
           this->N_d, TF(dt),
           gd.istart, gd.jstart, gd.kstart,
           gd.iend, gd.jend, gd.kend,
           gd.icells, gd.ijcells, temp_g->fld_g);
    cuda_check_error();

    cuda_safe_call(cudaMemcpy(temp->fld.data(), temp_g->fld_g, gd.ncells*sizeof(TF), cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(fields.st.at("qr")->fld.data(), fields.st.at("qr")->fld_g, gd.ncells * sizeof(TF), cudaMemcpyDeviceToHost));

    fields.release_tmp_g(ql);
    fields.release_tmp_g(qi);

    auto tmp1 = fields.get_tmp_g();
    auto tmp2 = fields.get_tmp_g();
    auto tmp3 = fields.get_tmp_g();
    auto tmp4 = fields.get_tmp_g();

    fields.release_tmp_g(temp_g);
    fields.release_tmp(temp);

    // Falling rain.
    sedimentation_ss08_g<<<grid2dGPU, block2dGPU>>>(
            fields.st.at("qr")->fld_g,
            tmp1->fld_g, tmp2->fld_g,
            tmp3->fld_g, tmp4->fld_g,
            fields.sp.at("qr")->fld_g,
            fields.rhoref_g,
            gd.dzi_g, gd.dz_g,
            TF(dt),
            a_r<TF>, b_r<TF>, c_r<TF>, d_r<TF>, N_0r<TF>,
            qr_min<TF>,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend, gd.jend, gd.kend,
            gd.icells, gd.ijcells);
   cuda_check_error();

   fields.release_tmp_g(tmp1);
   fields.release_tmp_g(tmp2);
   fields.release_tmp_g(tmp3);
   fields.release_tmp_g(tmp4);

   // cleanup
   cudaDeviceSynchronize();

   stats.calc_tend(*fields.st.at("thl"), tend_name);
   stats.calc_tend(*fields.st.at("qt" ), tend_name);
   stats.calc_tend(*fields.st.at("qr" ), tend_name);
   stats.calc_tend(*fields.st.at("qs" ), tend_name);
   stats.calc_tend(*fields.st.at("qg" ), tend_name);
}

template<typename TF>
unsigned long Microphys_nsw6<TF>::get_time_limit(unsigned long idt, const double dt)
{
   return 1000000000000; // DO NOT CHANGE
}
#endif

template class Microphys_nsw6<double>;
template class Microphys_nsw6<float>;
