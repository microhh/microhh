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
            const TF db_ref, const TF dt,
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
                    lambda*T_soil[ijk_s] + TF(3)*sigma_b<TF>*fm::pow4(T_bot);
                const TF denom = fH + fLE*dqsatdT_bot[ij] + lambda + TF(4)*sigma_b<TF>*fm::pow3(T_bot);
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

    template<typename TF, bool sw_constant_z0>
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
            const TF db_ref, const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend_soil,
            const int icells, const int ijcells,
            const bool use_cs_veg,
            const std::string name)
    {
        const TF thvref_bot = thvrefh[kstart];
        const TF exner_bot = exnerh[kstart];
        const TF rho_bot = rhorefh[kstart];
        const TF p_bot = prefh[kstart];

        int max_iters = 0;

        const TF eps_thl = TF(1e-6);
        const TF max_step = TF(1);
        const TF eps_seb = TF(1e-1);

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij    = i + j*icells;
                const int ijk   = ij + kstart*ijcells;
                const int ijk_s = ij + (kend_soil-1)*ijcells;

                int it;
                const int max_it = 100; // usually 2-4 is enough...
                for (it=0; it<max_it; ++it)
                {
                    const TF T_bot = thl_bot[ij] * exner_bot;
                    const TF qsat_bot = tmf::qsat(p_bot, T_bot);

                    // Disable canopy resistance in case of dew fall
                    const TF rs_lim = qsat_bot < qt[ijk] ? TF(0) : rs[ij];

                    // Switch between skin heat capacity or not
                    const TF cs_veg_lim = use_cs_veg ? cs_veg[ij] : TF(0);

                    // Surface layer calculations
                    const TF bbot = tmf::buoyancy_no_ql(thl_bot[ij], qt_bot[ij], thvrefh[kstart]);
                    const TF db = b[ijk] - bbot + db_ref;

                    if (sw_constant_z0)
                        obuk[ij] = bsk::calc_obuk_noslip_dirichlet_lookup(
                                zL_sl, f_sl, nobuk[ij], du_tot[ij], db, zsl);
                    else
                        obuk[ij] = bsk::calc_obuk_noslip_dirichlet_iterative(
                                obuk[ij], du_tot[ij], db, zsl, z0m[ij], z0h[ij]);

                    ustar[ij] = du_tot[ij] * most::fm(zsl, z0m[ij], obuk[ij]);

                    const TF ra = TF(1) / (ustar[ij] * most::fh(zsl, z0h[ij], obuk[ij]));

                    // Switch conductivity skin layer stable/unstable conditions
                    const TF lambda = db > 0 ? lambda_stable[ij] : lambda_unstable[ij];

                    // Calculate fluxes
                    H[ij]  = rho_bot * cp<TF> / ra * (T_bot - T[ij]);
                    LE[ij] = rho_bot * Lv<TF> / (ra + rs_lim) * (qsat_bot - qt[ij]);
                    G[ij]  = lambda * (T_bot - T_soil[ijk_s]);
                    const TF lw_up_n = sigma_b<TF> * fm::pow4(T_bot);

                    const TF Qn = sw_dn[ij] - sw_up[ij] + lw_dn[ij] - lw_up_n;
                    const TF seb = Qn - H[ij] - LE[ij] - G[ij];

                    if (it > max_it-10)
                    {
                        std::cout << "------- " << name << " , it="  << it << " -------" << std::endl;
                        std::cout << "thl_bot=" << thl_bot[ij] << " , qsat_bot=" << qsat_bot << " , db=" << db << std::endl;
                        std::cout << "obuk=" << obuk[ij] << " , ustar=" << ustar[ij] << " , ra=" << ra << std::endl;
                        std::cout << "H=" << H[ij] << " , LE=" << LE[ij] << " , G=" << G[ij] << " , L_up=" << lw_up_n << std::endl;
                        std::cout << "Qn=" << Qn << " , seb=" << seb << std::endl;
                    }

                    if (std::abs(seb) < eps_seb)
                    {
                        bfluxbot[ij] = -ustar[ij] * db * most::fh(zsl, z0h[ij], obuk[ij]);
                        qt_bot[ij] = qt[ijk] + LE[ij] * ra / (rho_bot * Lv<TF>);
                        break;
                    }

                    // Solve for perturbed thl_bot
                    const TF T_bot_p = (thl_bot[ij]+eps_thl) * exner_bot;
                    const TF qsat_bot_p = tmf::qsat(p_bot, T_bot_p);

                    // Disable canopy resistance in case of dew fall
                    const TF rs_lim_p = qsat_bot_p < qt[ijk] ? TF(0) : rs[ij];

                    // Surface layer calculations
                    const TF bbot_p = tmf::buoyancy_no_ql((thl_bot[ij]+eps_thl), qt_bot[ij], thvrefh[kstart]);
                    const TF db_p = b[ijk] - bbot_p + db_ref;

                    // AARGH, for now limited to iterative solver....
                    // NOTE: I left obuk[ij] as a starting point here......
                    TF obuk_p;
                    if (sw_constant_z0)
                        obuk_p = bsk::calc_obuk_noslip_dirichlet_lookup(
                                zL_sl, f_sl, nobuk[ij], du_tot[ij], db_p, zsl);
                    else
                        obuk_p = bsk::calc_obuk_noslip_dirichlet_iterative(
                                obuk[ij], du_tot[ij], db, zsl, z0m[ij], z0h[ij]);

                    const TF ustar_p = du_tot[ij] * most::fm(zsl, z0m[ij], obuk_p);

                    const TF ra_p = TF(1) / (ustar_p * most::fh(zsl, z0h[ij], obuk_p));

                    // Switch conductivity skin layer stable/unstable conditions
                    const TF lambda_p = db_p > 0 ? lambda_stable[ij] : lambda_unstable[ij];

                    // Calculate fluxes
                    const TF H_p  = rho_bot * cp<TF> / ra_p * (T_bot_p - T[ij]);
                    const TF LE_p = rho_bot * Lv<TF> / (ra_p + rs_lim) * (qsat_bot_p - qt[ij]);
                    const TF G_p  = lambda * (T_bot_p - T_soil[ijk_s]);
                    const TF lw_up_n_p = sigma_b<TF> * fm::pow4(T_bot_p);

                    const TF Qn_p = sw_dn[ij] - sw_up[ij] + lw_dn[ij] - lw_up_n_p;
                    const TF seb_p = Qn_p - H_p - LE_p - G_p;

                    const TF slope = (seb_p - seb) / eps_thl;
                    const TF dtheta = -seb / slope;

                    if (it > max_it-10)
                    {
                        std::cout << "--------- Perturbed: ----------" << std::endl;
                        std::cout << "obuk=" << obuk_p << " , ustar=" << ustar_p << " , ra=" << ra_p << std::endl;
                        std::cout << "H=" << H_p << " , LE=" << LE_p << " , G=" << G_p << " , L_up=" << lw_up_n_p << std::endl;
                        std::cout << "Qn=" << Qn_p << " , seb=" << seb_p << std::endl;
                        std::cout << "slope=" << slope << " , dTs=" << dtheta << std::endl;
                    }

                    thl_bot[ij] += std::max(-max_step, std::min(dtheta, max_step));
                }

                max_iters = std::max(max_iters, it);

                if (it == max_it)
                {
                    std::string error = "SEB solver did not converge for tile: \"" + name + "\"";
                    #ifdef USEMPI
                        std::cout << "SINGLE PROCESS EXCEPTION: " << error << std::endl;
                        MPI_Abort(MPI_COMM_WORLD, 1);
                    #else
                        throw std::runtime_error(error);
                    #endif
                }
            }

            std::cout << name << " max iters=" << max_iters << std::endl;
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
