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

#ifndef LAND_SURFACE_KERNELS_H
#define LAND_SURFACE_KERNELS_H

#include <iomanip>

#include "constants.h"
#include "thermo_moist_functions.h"
#include "boundary_surface_kernels.h"
#include "monin_obukhov.h"
#include "fast_math.h"

using namespace Constants;

namespace Land_surface_kernels
{
    namespace tmf = Thermo_moist_functions;
    namespace most = Monin_obukhov;
    namespace bsk = Boundary_surface_kernels;
    namespace fm = Fast_math;

    template<typename TF>
    void init_tile(Surface_tile<TF>& tile, const int ijcells)
    {
        tile.fraction.resize(ijcells);
        tile.thl_bot.resize(ijcells);
        tile.qt_bot.resize(ijcells);

        tile.obuk.resize(ijcells);
        tile.ustar.resize(ijcells);
        tile.bfluxbot.resize(ijcells);
        tile.nobuk.resize(ijcells);

        tile.rs.resize(ijcells);
        tile.H.resize(ijcells);
        tile.LE.resize(ijcells);
        tile.G.resize(ijcells);
        tile.S.resize(ijcells);
    }

    template<typename TF>
    void calc_tile_fractions(
            TF* const restrict tile_frac_veg,
            TF* const restrict tile_frac_soil,
            TF* const restrict tile_frac_wet,
            const TF* const restrict wl,
            const TF* const restrict c_veg,
            const TF* const restrict lai,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;

                const TF wlm = Constants::wlmax<TF> * (TF(1)-c_veg[ij] + c_veg[ij]*lai[ij]);

                tile_frac_wet[ij]  = std::min(TF(1), wl[ij]/wlm);
                tile_frac_veg[ij]  = (TF(1)-tile_frac_wet[ij]) * c_veg[ij];
                tile_frac_soil[ij] = (TF(1)-tile_frac_wet[ij]) * (TF(1)-c_veg[ij]);
            }
    }

    template<typename TF>
    void calc_liquid_water_reservoir(
            TF* const restrict wl_tend,
            TF* const restrict interception,
            TF* const restrict throughfall,
            const TF* const restrict wl,
            const TF* const restrict LE_veg,
            const TF* const restrict LE_soil,
            const TF* const restrict LE_wet,
            const TF* const restrict tile_frac_veg,
            const TF* const restrict tile_frac_soil,
            const TF* const restrict tile_frac_wet,
            const TF* const restrict rain_rate,
            const TF* const restrict c_veg,
            const TF* const restrict lai,
            const double subdt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        const TF intercept_eff = TF(0.5);
        const TF to_ms  = TF(1) / (Constants::rho_w<TF> * Constants::Lv<TF>);
        const TF subdti = TF(1) / subdt;

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;

                // Convert rain rate from kg m-2 s-1 to m s-1
                const TF rr_ms = rain_rate[ij]/Constants::rho_w<TF>;

                // Max `wl` accounting for vegetation fraction and LAI
                const TF wlm = Constants::wlmax<TF> * (TF(1) - c_veg[ij] + c_veg[ij] * lai[ij]);

                // Max and min possible tendencies
                const TF wl_tend_max = (wlm - wl[ij]) * subdti - wl_tend[ij];
                const TF wl_tend_min = (    - wl[ij]) * subdti - wl_tend[ij];

                // Tendency due to evaporation from liquid water reservoir/tile.
                const TF wl_tend_liq = -std::max(TF(0), tile_frac_wet[ij] * LE_wet[ij] * to_ms);

                // Tendency due to dewfall into vegetation/soil/liquid water tiles
                const TF wl_tend_dew = -( std::min(TF(0), tile_frac_wet[ij]  * LE_wet[ij]  * to_ms)
                                        + std::min(TF(0), tile_frac_veg[ij]  * LE_veg[ij]  * to_ms)
                                        + std::min(TF(0), tile_frac_soil[ij] * LE_soil[ij] * to_ms) );

                // Tendency due to interception of precipitation by vegetation
                // Rain rate is positive downwards, so minus is excluded.
                const TF wl_tend_precip = intercept_eff * c_veg[ij] * rr_ms;

                // Total and limited tendencies
                const TF wl_tend_sum = wl_tend_liq + wl_tend_dew + wl_tend_precip;
                const TF wl_tend_lim = std::min(wl_tend_max, std::max(wl_tend_min,  wl_tend_sum));

                // Diagnose throughfall and interception
                throughfall[ij] =
                    -(TF(1)-c_veg[ij]) * rr_ms
                    -(TF(1)-intercept_eff) * c_veg[ij] * rr_ms +
                    std::min(TF(0), wl_tend_lim - wl_tend_sum);

                interception[ij] = std::max(TF(0), wl_tend_lim);

                wl_tend[ij] += wl_tend_lim;
            }
    }

    template<typename TF>
    void calc_resistance_functions(
            TF* const restrict f1,
            TF* const restrict f2,
            TF* const restrict f2b,
            TF* const restrict f3,
            const TF* const restrict sw_dn,
            const TF* const restrict theta,
            const TF* const restrict theta_mean_n,
            const TF* const restrict vpd,
            const TF* const restrict gD,
            const TF* const restrict c_veg,
            const TF* const restrict theta_wp,
            const TF* const restrict theta_fc,
            const TF* const restrict theta_res,
            const int* const restrict soil_index,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kend,
            const int icells, const int ijcells)
    {
        // Constants f1 calculation:
        const TF a_f1 = 0.81;
        const TF b_f1 = 0.004;
        const TF c_f1 = 0.05;

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                const int ijk = i + j*icells + (kend-1)*ijcells;    // Top soil layer
                const int si  = soil_index[ij];

                // f1: reduction vegetation resistance as f(sw_in):
                const TF sw_dn_lim = std::max(TF(0), sw_dn[ij]);
                f1[ij] = TF(1)/std::min( TF(1), (b_f1*sw_dn_lim + c_f1) / (a_f1 * (b_f1*sw_dn_lim + TF(1))) );

                // f2: reduction vegetation resistance as f(theta):
                f2[ij] = TF(1)/std::min( TF(1), std::max(TF(1e-9), theta_mean_n[ij]) );

                // f3: reduction vegetation resistance as f(VPD):
                f3[ij] = TF(1)/exp(-gD[ij] * vpd[ij]);

                // f2b: reduction soil resistance as f(theta)
                const TF theta_min = c_veg[ij] * theta_wp[si] + (TF(1)-c_veg[ij]) * theta_res[si];
                const TF theta_rel = (theta[ijk] - theta_min) / (theta_fc[si] - theta_min);
                f2b[ij] = TF(1)/std::min(TF(1), std::max(TF(1e-9), theta_rel));
            }
    }

    template<typename TF>
    void calc_canopy_resistance(
            TF* const restrict rs,
            const TF* const restrict rs_min,
            const TF* const restrict lai,
            const TF* const restrict f1,
            const TF* const restrict f2,
            const TF* const restrict f3,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                rs[ij] = rs_min[ij] / lai[ij] * f1[ij] * f2[ij] * f3[ij];
            }
    }

    template<typename TF>
    void calc_soil_resistance(
            TF* const restrict rs,
            const TF* const restrict rs_min,
            const TF* const restrict f2b,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                rs[ij] = rs_min[ij] * f2b[ij];
            }
    }

    template<typename TF, bool sw_constant_z0>
    void calc_stability(
            TF* const restrict ustar,
            TF* const restrict obuk,
            TF* const restrict bfluxbot,
            TF* const restrict ra,
            int* const restrict nobuk,
            const TF* const restrict dutot,
            const TF* const restrict b,
            const TF* const restrict bbot,
            const TF* const restrict z0m,
            const TF* const restrict z0h,
            const float* const restrict zL_sl,
            const float* const restrict f_sl,
            const TF db_ref,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int jcells,
            const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;
                const TF db = b[ijk] - bbot[ij] + db_ref;

                if (sw_constant_z0)
                    obuk[ij] = bsk::calc_obuk_noslip_dirichlet_lookup(
                            zL_sl, f_sl, nobuk[ij], dutot[ij], db, zsl);
                else
                    obuk[ij] = bsk::calc_obuk_noslip_dirichlet_iterative(
                            obuk[ij], dutot[ij], db, zsl, z0m[ij], z0h[ij]);

                ustar[ij] = dutot[ij] * most::fm(zsl, z0m[ij], obuk[ij]);
                bfluxbot[ij] = - ustar[ij] * db * most::fh(zsl, z0h[ij], obuk[ij]);
                ra[ij]  = TF(1) / (ustar[ij] * most::fh(zsl, z0h[ij], obuk[ij]));
            }
    }

    template<typename TF>
    void calc_fluxes(
            TF* const restrict H,
            TF* const restrict LE,
            TF* const restrict G,
            TF* const restrict S,
            TF* const restrict thl_bot,
            TF* const restrict qt_bot,
            const TF* const restrict T,
            const TF* const restrict qt,
            const TF* const restrict T_soil,
            const TF* const restrict qsat_bot,
            const TF* const restrict dqsatdT_bot,
            const TF* const restrict ra,
            const TF* const restrict rs,
            const TF* const restrict lambda_stable,
            const TF* const restrict lambda_unstable,
            const TF* const restrict cs_veg,
            const TF* const restrict sw_dn,
            const TF* const restrict sw_up,
            const TF* const restrict lw_dn,
            const TF* const restrict lw_up,
            const TF* const restrict b,
            const TF* const restrict b_bot,
            const TF* const restrict rhorefh,
            const TF* const restrict exnerh,
            const TF db_ref,
            const TF emis_sfc,
            const TF dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend_soil,
            const int icells, const int ijcells,
            bool use_cs_veg, std::string name)
    {
        const TF exner_bot = exnerh[kstart];
        const TF rho_bot = rhorefh[kstart];

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij    = i + j*icells;
                const int ijk   = ij + kstart*ijcells;
                const int ijk_s = ij + (kend_soil-1)*ijcells;

                const TF T_bot = thl_bot[ij] * exner_bot;

                // Disable canopy resistance in case of dew fall
                const TF rs_lim = qsat_bot[ij] < qt[ijk] ? TF(0) : rs[ij];

                // Switch between skin heat capacity or not
                const TF cs_veg_lim = use_cs_veg ? cs_veg[ij] : TF(0);

                // Switch conductivity skin layer stable/unstable conditions
                const TF db = b[ijk] - b_bot[ij] + db_ref;
                const TF lambda = db > 0 ? lambda_stable[ij] : lambda_unstable[ij];

                // Recuring factors
                const TF fH  = rho_bot * cp<TF> / ra[ij];
                const TF fLE = rho_bot * Lv<TF> / (ra[ij] + rs_lim);

                // Net radiation; negative sign = net input of energy at surface
                const TF Qnet = sw_dn[ij] - sw_up[ij] + lw_dn[ij] - lw_up[ij];

                // Solve for the new surface temperature
                const TF num =
                    Qnet + lw_up[ij] + fH*T[ij] +
                    fLE*(qt[ijk] + dqsatdT_bot[ij]*T_bot - qsat_bot[ij]) +
                    lambda*T_soil[ijk_s] + TF(3)*emis_sfc*sigma_b<TF>*fm::pow4(T_bot) - (TF(1)-emis_sfc) * lw_dn[ij];
                const TF denom = fH + fLE*dqsatdT_bot[ij] + lambda + TF(4)*emis_sfc*sigma_b<TF>*fm::pow3(T_bot);
                const TF T_bot_new = (num + cs_veg_lim/dt*T_bot) / (denom + cs_veg_lim/dt);

                // Update qsat with linearised relation, to make sure that the SEB closes
                const TF dT_bot = T_bot_new - T_bot;
                const TF qsat_new  = qsat_bot[ij] + dqsatdT_bot[ij] * dT_bot;

                // Calculate surface fluxes
                H [ij] = fH  * (T_bot_new - T[ij]);
                LE[ij] = fLE * (qsat_new - qt[ijk]);
                G [ij] = lambda * (T_bot_new - T_soil[ijk_s]);
                S [ij] = cs_veg_lim * (T_bot_new - T_bot)/dt;

                // Update skin values
                thl_bot[ij] = T_bot_new / exner_bot;
                qt_bot[ij]  = qt[ijk] + LE[ij] * ra[ij] / (rho_bot * Lv<TF>);
            }
    }

    template<typename TF>
    TF smooth(const TF v1, const TF v2,
              const TF x1, const TF x2,
              const TF x)
    {
        const TF f = std::max(TF(0), std::min(TF(1), (x2-x)/(x2-x1)));
        return f*v1 + (TF(1)-f)*v2;
    }

    //template<typename TF, bool sw_constant_z0>
    template<typename TF>
    void calc_stability_and_fluxes(
            TF* const restrict H,
            TF* const restrict LE,
            TF* const restrict G,
            TF* const restrict S,
            TF* const restrict thl_bot,
            TF* const restrict qt_bot,
            TF* const restrict ustar,
            TF* const restrict obuk,
            int* const restrict nobuk,
            TF* const restrict bfluxbot,
            const TF* const restrict sw_dn,
            const TF* const restrict sw_up,
            const TF* const restrict lw_dn,
            const TF* const restrict lw_up,
            const TF* const restrict du_tot,
            const TF* const restrict T,
            const TF* const restrict qt,
            const TF* const restrict b,
            const TF* const restrict T_soil,
            const TF* const restrict rs,
            const TF* const restrict z0m,
            const TF* const restrict z0h,
            const TF* const restrict lambda_stable,
            const TF* const restrict lambda_unstable,
            const TF* const restrict cs_veg,
            const TF* const restrict rhorefh,
            const TF* const restrict exnerh,
            const TF* const restrict thvrefh,
            const TF* const restrict prefh,
            const float* const restrict f_sl,
            const float* const restrict zL_sl,
            const TF db_ref,
            const TF emis_sfc,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend_soil,
            const int icells, const int ijcells,
            const bool use_cs_veg,
            const bool sw_constant_z0,
            const std::string name,
            const int iter,
            const int subs,
            const int mpiid)
    {
        const TF thvref_bot = thvrefh[kstart];
        const TF exner_bot = exnerh[kstart];
        const TF rho_bot = rhorefh[kstart];
        const TF p_bot = prefh[kstart];

        // Lambda function for calculating SEB closure (SEB = Qn - H - LE - G).
        auto calc_seb_closure = [&](
                const TF thl_bot, const TF qt_bot, TF& obuk, bool sw_constant_z0,
                const int ij, const int ijk, const int ijk_s)
        {
            const TF T_bot = thl_bot * exner_bot;
            const TF qsat_bot = tmf::qsat(p_bot, T_bot);

            // Disable canopy resistance in case of dew fall
            //const TF rs_lim = qsat_bot < qt[ijk] ? TF(0) : rs[ij];

            // Smooth rs around qsat_bot == qt:
            const TF rs_lim = smooth(
                    TF(0), rs[ij], qt[ijk], qt[ijk]+TF(1e-4), qsat_bot);

            // Surface layer calculations
            const TF bbot = tmf::buoyancy_no_ql(thl_bot, qt_bot, thvref_bot);
            const TF db = b[ijk] - bbot + db_ref;

            if (sw_constant_z0)
                obuk = bsk::calc_obuk_noslip_dirichlet_lookup(
                        zL_sl, f_sl, nobuk[ij], du_tot[ij], db, zsl);
            else
                obuk = bsk::calc_obuk_noslip_dirichlet_iterative(
                            obuk, du_tot[ij], db, zsl, z0m[ij], z0h[ij]);

            const TF ustar = du_tot[ij] * most::fm(zsl, z0m[ij], obuk);
            const TF ra = TF(1) / (ustar * most::fh(zsl, z0h[ij], obuk));

            // Switch conductivity skin layer stable/unstable conditions
            //const TF lambda = db > TF(0) ? lambda_stable[ij] : lambda_unstable[ij];

            // Smooth Lambda_skin around db=0:
            const TF db_eps = TF(0.005);
            const TF lambda = smooth(
                lambda_stable[ij], lambda_unstable[ij], db_eps, -db_eps, db);

            // Calculate fluxes
            const TF H  = rho_bot * cp<TF> / ra * (T_bot - T[ij]);
            const TF LE = rho_bot * Lv<TF> / (ra + rs_lim) * (qsat_bot - qt[ijk]);
            const TF G  = lambda * (T_bot - T_soil[ijk_s]);
            const TF lw_up_n = emis_sfc * sigma_b<TF> * fm::pow4(T_bot) + (TF(1) - emis_sfc) * lw_dn[ij];

            const TF Qn = sw_dn[ij] - sw_up[ij] + lw_dn[ij] - lw_up_n;
            const TF seb = Qn - H - LE - G;

            return seb;
        };

        int max_iters = 0;

        const TF eps_thl = TF(1e-6);
        const TF eps_seb = TF(1e-1);

        const TF max_step = TF(1);
        const TF min_step = TF(0.001);
        const int max_it = 1000;

        const TF max_slope = TF(100);

        int did_not_converge = 0;

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij    = i + j*icells;
                const int ijk   = ij + kstart*ijcells;
                const int ijk_s = ij + (kend_soil-1)*ijcells;

                TF qt_bot_p = qt_bot[ij];
                TF obuk_p = obuk[ij];

                // Store data for debugging output
                const TF thl_bot_in = thl_bot[ij];
                const TF qt_bot_in = qt_bot[ij];
                const TF obuk_in = obuk[ij];

                int it;
                for (it=0; it<max_it; ++it)
                {
                    const TF max_dtheta = max_step + it * (min_step-max_step)/(max_it-1);

                    // Solve for actual thl_bot
                    const TF seb_0 = calc_seb_closure(
                            thl_bot[ij], qt_bot[ij], obuk[ij], sw_constant_z0, ij, ijk, ijk_s);

                    // Exit if SEB closes (Qn = H+LE+G)
                    if (std::abs(seb_0) < eps_seb)
                        break;

                    // Solve for perturbed thl_bot
                    const TF seb_1 = calc_seb_closure(
                            thl_bot[ij]+eps_thl, qt_bot_p, obuk_p, sw_constant_z0, ij, ijk, ijk_s);

                    // Calculate increment thl_bot
                    const TF slope = (seb_1 - seb_0) / eps_thl;
                    const TF dtheta = std::max(-max_dtheta, std::min(-seb_0 / slope, max_dtheta));

                    // Increment thl_bot
                    thl_bot[ij] += dtheta;
                }

                // Calculate/set final values
                const TF T_bot = thl_bot[ij] * exner_bot;
                const TF qsat_bot = tmf::qsat(p_bot, T_bot);

                // Disable canopy resistance in case of dew fall
                const TF rs_lim = qsat_bot < qt[ijk] ? TF(0) : rs[ij];

                // Surface layer calculations
                const TF bbot = tmf::buoyancy_no_ql(thl_bot[ij], qt_bot[ij], thvref_bot);
                const TF db = b[ijk] - bbot + db_ref;

                ustar[ij] = du_tot[ij] * most::fm(zsl, z0m[ij], obuk[ij]);
                const TF ra = TF(1) / (ustar[ij] * most::fh(zsl, z0h[ij], obuk[ij]));

                // Switch conductivity skin layer stable/unstable conditions
                const TF lambda = db > TF(0) ? lambda_stable[ij] : lambda_unstable[ij];

                // Calculate fluxes
                H[ij]  = rho_bot * cp<TF> / ra * (T_bot - T[ij]);
                LE[ij] = rho_bot * Lv<TF> / (ra + rs_lim) * (qsat_bot - qt[ijk]);
                G[ij]  = lambda * (T_bot - T_soil[ijk_s]);
                bfluxbot[ij] = -ustar[ij] * db * most::fh(zsl, z0h[ij], obuk[ij]);

                // Update qt_bot
                qt_bot[ij] = qt[ijk] + LE[ij] * ra / (rho_bot * Lv<TF>);

                if (it == max_it)
                {
                    did_not_converge += 1;

                    std::cout << i << " " << j << " " <<  iter << " " << subs << " " << name << std::endl;
                    std::cout << "Input:" << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " thl_bot = " << thl_bot_in << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " qt_bot = " << qt_bot_in << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " T_a = " << T[ij] << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " qt = " << qt[ijk] << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " du_tot = " << du_tot[ij] << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " T_soil = " << T_soil[ijk_s] << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " sw_dn = " << sw_dn[ij] << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " sw_up = " << sw_up[ij] << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " lw_dn = " << lw_dn[ij] << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " lw_up = " << lw_up[ij] << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " exner_bot = " << exner_bot << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " p_bot = " << p_bot << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " rho_bot = " << rho_bot << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " rs = " << rs[ij] << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " b = " << b[ijk] << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " db_ref = " << db_ref << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " thvrefh = " << thvrefh[kstart] << std::endl;
                    std::cout << mpiid << " " << i << " " << j << std::setprecision(20) << " obuk = " << obuk_in << std::endl;

                    // NOTE BvS: I'm letting the non-converged iteration pass for now.
                    // In all situations that I have seen, the iteration does not converge
                    // when `slope` (dSEB/dthl_bot) becomes very large. In those case,
                    // the solver is more or less converged, and just jumping between
                    // two almost similar `thl_bot` values.

                    std::string error = "SEB solver did not converge for tile: \"" + name + "\"";
                    #ifdef USEMPI
                        std::cout << "SINGLE PROCESS EXCEPTION: " << error << std::endl;
                        MPI_Abort(MPI_COMM_WORLD, 1);
                    #else
                        throw std::runtime_error(error);
                    #endif
                }

                // Just for statistics...
                max_iters = std::max(max_iters, it);
            }

            if (did_not_converge > 0)
                std::cout << "Tile \"" + name + "\", SEB solver did not converge for: " << did_not_converge << " grid points!" << std::endl;
            else if (max_iters >= 20)
                std::cout << "Tile \"" + name + "\" required: " << max_iters << " SEB iterations..." << std::endl;
    }

    template<typename TF>
    void calc_bcs(
            TF* const restrict thl_bot,
            TF* const restrict qt_bot,
            const TF* const restrict thl,
            const TF* const restrict qt,
            const TF* const restrict H_veg,
            const TF* const restrict H_soil,
            const TF* const restrict H_wet,
            const TF* const restrict LE_veg,
            const TF* const restrict LE_soil,
            const TF* const restrict LE_wet,
            const TF* const restrict tile_frac_veg,
            const TF* const restrict tile_frac_soil,
            const TF* const restrict tile_frac_wet,
            const TF* const restrict ra,
            const TF* const restrict rhorefh,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int ijcells)
    {
        const TF rhocp_i = TF(1) / (rhorefh[kstart] * cp<TF>);
        const TF rholv_i = TF(1) / (rhorefh[kstart] * Lv<TF>);

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                const int ijk = ij + kstart*ijcells;

                // Tile averaged surface fluxes
                const TF wthl = (
                    tile_frac_veg [ij] * H_veg [ij] +
                    tile_frac_soil[ij] * H_soil[ij] +
                    tile_frac_wet [ij] * H_wet [ij] ) * rhocp_i;

                const TF wqt = (
                    tile_frac_veg [ij] * LE_veg [ij] +
                    tile_frac_soil[ij] * LE_soil[ij] +
                    tile_frac_wet [ij] * LE_wet [ij] ) * rholv_i;

                // Calculate surface values
                thl_bot[ij] = thl[ijk] + wthl * ra[ij];
                qt_bot [ij] = qt[ijk]  + wqt  * ra[ij];
            }
    }

    template<typename TF>
    void set_water_bcs(
            TF* const restrict thl_bot,
            TF* const restrict qt_bot,
            const int* const restrict water_mask,
            const TF* const restrict exnerh,
            const TF* const restrict prefh,
            const TF tskin_water,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int icells)
    {
        // Water temperature is fixed (for now..), so thl_bot and qt_bot are identical everywhere.
        const TF thl_bot_val = tskin_water * TF(1.)/exnerh[kstart];
        const TF qt_bot_val  = Thermo_moist_functions::qsat(prefh[kstart], tskin_water);

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;

                if (water_mask[ij])
                {
                    // Set surface values
                    thl_bot[ij] = thl_bot_val;
                    qt_bot [ij] = qt_bot_val;
                }
            }
    }

    template<typename TF>
    void set_water_tiles(
            TF* const restrict c_veg,
            TF* const restrict c_soil,
            TF* const restrict c_wet,
            TF* const restrict H_veg,
            TF* const restrict H_soil,
            TF* const restrict H_wet,
            TF* const restrict LE_veg,
            TF* const restrict LE_soil,
            TF* const restrict LE_wet,
            TF* const restrict G_veg,
            TF* const restrict G_soil,
            TF* const restrict G_wet,
            TF* const restrict rs_veg,
            TF* const restrict rs_soil,
            TF* const restrict rs_wet,
            const int* const restrict water_mask,
            const TF* const restrict thl,
            const TF* const restrict qt,
            const TF* const restrict thl_bot,
            const TF* const restrict qt_bot,
            const TF* const restrict ra,
            const TF* const restrict rhoh,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int ijcells)
    {
        const TF rhocp = rhoh[kstart] * Constants::cp<TF>;
        const TF rholv = rhoh[kstart] * Constants::Lv<TF>;

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                const int ijk = ij + kstart*ijcells;

                if (water_mask[ij])
                {
                    c_veg[ij]  = TF(0);
                    c_soil[ij] = TF(0);
                    c_wet[ij]  = TF(1);

                    H_veg[ij]  = TF(0);
                    H_soil[ij] = TF(0);
                    H_wet[ij]  = rhocp / ra[ij] * (thl_bot[ij] - thl[ijk]);

                    LE_veg[ij]  = TF(0);
                    LE_soil[ij] = TF(0);
                    LE_wet[ij]  = rholv / ra[ij] * (qt_bot[ij] - qt[ijk]);

                    G_veg[ij]  = TF(0);
                    G_soil[ij] = TF(0);
                    G_wet[ij]  = TF(0);

                    rs_veg[ij]  = TF(0);
                    rs_soil[ij] = TF(0);
                    rs_wet[ij]  = TF(0);
                }
            }
    }

    template<typename TF>
    void calc_tiled_mean(
            TF* const restrict fld_mean,
            const TF* const restrict fld_veg,
            const TF* const restrict fld_soil,
            const TF* const restrict fld_wet,
            const TF* const restrict tile_frac_veg,
            const TF* const restrict tile_frac_soil,
            const TF* const restrict tile_frac_wet,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;

                fld_mean[ij] =
                    tile_frac_veg [ij] * fld_veg [ij] +
                    tile_frac_soil[ij] * fld_soil[ij] +
                    tile_frac_wet [ij] * fld_wet [ij];
            }
    }

    template<typename TF>
    void scale_tile_with_fraction(
            TF* const restrict fld_scaled,
            const TF* const restrict fld,
            const TF* const restrict tile_frac,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                fld_scaled[ij] = fld[ij] * tile_frac[ij];
            }
    }
}
#endif
