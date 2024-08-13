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

#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "diff.h"
#include "stats.h"
#include "cross.h"
#include "column.h"
#include "thermo.h"
#include "thermo_moist_functions.h"
#include "fast_math.h"

#include "constants.h"
#include "microphys.h"
#include "microphys_2mom_warm.h"

namespace
{
    using namespace Constants;
    using namespace Thermo_moist_functions;
    using namespace Micro_2mom_warm_constants;
    using namespace Micro_2mom_warm_functions;

    template<typename TF>
    void remove_negative_values(TF* const restrict field,
                                const int istart, const int jstart, const int kstart,
                                const int iend,   const int jend,   const int kend,
                                const int jj,     const int kk)
    {
        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    field[ijk] = std::max(TF(0.), field[ijk]);
                }
    }

    template<typename TF>
    void zero_field(TF* const restrict field, const int ncells)
    {
        for (int n=0; n<ncells; n++)
            field[n] = TF(0);
    }

    template<typename TF>
    TF* get_tmp_slice(std::vector<std::shared_ptr<Field3d<TF>>> &tmp_fields, int &slice_counter,
                      const int jcells, const int ikcells)
    {
        const int tmp_index   = slice_counter / jcells;     // Which tmp field in tmp_fields vector?
        const int fld_index   = slice_counter % jcells;     // Which slice in tmp field?
        const int slice_start = fld_index * ikcells;        // Start index of slice

        slice_counter++;

        return &(tmp_fields[tmp_index]->fld[slice_start]);
    }

}

// Microphysics calculated over entire 3D field
namespace mp3d
{
    namespace fm = Fast_math;

    // Autoconversion: formation of rain drop by coagulating cloud droplets
    template<typename TF>
    void autoconversion(TF* const restrict qrt, TF* const restrict nrt,
                        TF* const restrict qtt, TF* const restrict thlt,
                        const TF* const restrict qr,  const TF* const restrict ql,
                        const TF* const restrict rho, const TF* const restrict exner, const TF nc,
                        const int istart, const int jstart, const int kstart,
                        const int iend,   const int jend,   const int kend,
                        const int jj, const int kk)
    {
        const TF x_star = 2.6e-10;       // SB06, list of symbols, same as UCLA-LES
        const TF k_cc   = 9.44e9;        // UCLA-LES (Long, 1974), 4.44e9 in SB06, p48
        const TF nu_c   = 1;             // SB06, Table 1., same as UCLA-LES
        const TF kccxs  = k_cc / (TF(20.) * x_star) * (nu_c+2)*(nu_c+4) / fm::pow2(nu_c+1);

        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    if (ql[ijk] > ql_min<TF>)
                    {
                        const TF xc      = rho[k] * ql[ijk] / nc;    // Mean mass of cloud drops [kg]
                        const TF tau     = TF(1.) - ql[ijk] / (ql[ijk] + qr[ijk] + dsmall);    // SB06, Eq 5
                        const TF phi_au  = TF(600.) * pow(tau, TF(0.68)) * fm::pow3(TF(1.) - pow(tau, TF(0.68)));    // UCLA-LES
                        const TF au_tend = rho_0<TF> * kccxs * fm::pow2(ql[ijk]) * fm::pow2(xc) *
                                               (TF(1.) + phi_au / fm::pow2(TF(1.)-tau)); // SB06, eq 4

                        qrt[ijk]  += au_tend;
                        nrt[ijk]  += au_tend * rho[k] / x_star;
                        qtt[ijk]  -= au_tend;
                        thlt[ijk] += Lv<TF> / (cp<TF> * exner[k]) * au_tend;
                    }
                }
    }

    // Accreation: growth of raindrops collecting cloud droplets
    template<typename TF>
    void accretion(TF* const restrict qrt, TF* const restrict qtt, TF* const restrict thlt,
                   const TF* const restrict qr,  const TF* const restrict ql,
                   const TF* const restrict rho, const TF* const restrict exner,
                   const int istart, const int jstart, const int kstart,
                   const int iend,   const int jend,   const int kend,
                   const int jj, const int kk)
    {
        const TF k_cr = 5.25; // SB06, p49

        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    if (ql[ijk] > ql_min<TF> && qr[ijk] > qr_min<TF>)
                    {
                        const TF tau     = TF(1.) - ql[ijk] / (ql[ijk] + qr[ijk]); // SB06, Eq 5
                        const TF phi_ac  = fm::pow4(tau / (tau + TF(5e-5))); // SB06, Eq 8
                        const TF ac_tend = k_cr * ql[ijk] *  qr[ijk] * phi_ac * sqrt(rho_0<TF> / rho[k]); // SB06, Eq 7

                        qrt[ijk]  += ac_tend;
                        qtt[ijk]  -= ac_tend;
                        thlt[ijk] += Lv<TF> / (cp<TF> * exner[k]) * ac_tend;
                    }
                }
    }


    // Calculate maximum sedimentation velocity
    template<typename TF>
    TF calc_max_sedimentation_cfl(TF* const restrict w_qr,
                                  const TF* const restrict qr, const TF* const restrict nr,
                                  const TF* const restrict rho, const TF* const restrict dzi,
                                  const double dt,
                                  const int istart, const int jstart, const int kstart,
                                  const int iend,   const int jend,   const int kend,
                                  const int icells, const int ijcells)
    {
        const TF w_max = 9.65; // 9.65=UCLA, 20=SS08, appendix A
        const TF a_R = 9.65;   // SB06, p51
        const TF c_R = 600;    // SB06, p51
        const TF Dv  = 25.0e-6;
        const TF b_R = a_R * exp(c_R*Dv); // UCLA-LES

        // Calculate sedimentation velocity at cell centre
        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;

                    if (qr[ijk] > qr_min<TF>)
                    {
                        // Calculate mean rain drop mass and diameter
                        const TF mr      = calc_rain_mass(qr[ijk], nr[ijk], rho[k]);
                        const TF dr      = calc_rain_diameter(mr);
                        const TF mur     = calc_mu_r(dr);
                        const TF lambdar = calc_lambda_r(mur, dr);

                        w_qr[ijk] = std::min(w_max, std::max(TF(0.1), a_R - b_R * TF(pow(TF(1.) + c_R/lambdar, TF(-1.)*(mur+TF(4.))))));
                    }
                    else
                    {
                        w_qr[ijk] = TF(0.);
                    }
                }

        // Mirror the values of w_qr over the ghost cells.
        // Bottom boundary.
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*icells + kstart*ijcells;
                w_qr[ijk-ijcells] = w_qr[ijk];
            }

        // Top boundary.
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*icells + (kend-1)*ijcells;
                w_qr[ijk+ijcells] = w_qr[ijk];
            }

        // Calculate maximum CFL based on interpolated velocity
        TF cfl_max = 1e-5;
        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;

                    const TF cfl_qr = TF(0.25) * (w_qr[ijk-ijcells] + TF(2.)*w_qr[ijk] + w_qr[ijk+ijcells]) * dzi[k] * TF(dt);
                    cfl_max = std::max(cfl_max, cfl_qr);
                }

        return cfl_max;
    }
}

// Microphysics calculated over 2D (xz) slices
namespace mp2d
{
    namespace fm = Fast_math;

    // Calculate microphysics properties which are used in multiple routines
    template<typename TF>
    void prepare_microphysics_slice(TF* const restrict rain_mass, TF* const restrict rain_diameter,
                                    TF* const restrict mu_r, TF* const restrict lambda_r,
                                    const TF* const restrict qr, const TF* const restrict nr,
                                    const TF* const restrict rho,
                                    const int istart, const int iend,
                                    const int kstart, const int kend,
                                    const int icells, const int ijcells, const int j)
    {

        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*icells + k*ijcells;
                const int ik  = i + k*icells;

                if (qr[ijk] > qr_min<TF>)
                {
                    rain_mass[ik]     = calc_rain_mass(qr[ijk], nr[ijk], rho[k]);
                    rain_diameter[ik] = calc_rain_diameter(rain_mass[ik]);
                    mu_r[ik]          = calc_mu_r(rain_diameter[ik]);
                    lambda_r[ik]      = calc_lambda_r(mu_r[ik], rain_diameter[ik]);
                }
                else
                {
                    rain_mass[ik]     = TF(0);
                    rain_diameter[ik] = TF(0);
                    mu_r[ik]          = TF(0);
                    lambda_r[ik]      = TF(0);
                }
            }
    }

    // Evaporation: evaporation of rain drops in unsaturated environment
    template<typename TF>
    void evaporation(TF* const restrict qrt, TF* const restrict nrt,
                     TF* const restrict qtt, TF* const restrict thlt,
                     const TF* const restrict qr, const TF* const restrict nr,
                     const TF* const restrict ql, const TF* const restrict qt, const TF* const restrict thl,
                     const TF* const restrict rho, const TF* const restrict exner, const TF* const restrict p,
                     const TF* const restrict rain_mass, const TF* const restrict rain_diameter,
                     const int istart, const int jstart, const int kstart,
                     const int iend,   const int jend,   const int kend,
                     const int jj, const int kk, const int j)
    {
        const TF lambda_evap = 1.; // 1.0 in UCLA, 0.7 in DALES

        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ik  = i + k*jj;
                const int ijk = i + j*jj + k*kk;

                if (qr[ijk] > qr_min<TF>)
                {
                    const TF mr  = rain_mass[ik];
                    const TF dr  = rain_diameter[ik];

                    const TF T   = thl[ijk] * exner[k] + (Lv<TF> * ql[ijk]) / (cp<TF> * exner[k]); // Absolute temperature [K]
                    const TF Glv = TF(1.) / (Rv<TF> * T / (esat_liq(T) * D_v<TF>) +
                                       (Lv<TF> / (K_t<TF> * T)) * (Lv<TF> / (Rv<TF> * T) - TF(1.))); // Cond/evap rate (kg m-1 s-1)?

                    const TF S   = (qt[ijk] - ql[ijk]) / qsat_liq(p[k], T) - TF(1.); // Saturation
                    const TF F   = 1.; // Evaporation excludes ventilation term from SB06 (like UCLA, unimportant term? TODO: test)

                    const TF ev_tend = TF(2.) * pi<TF> * dr * Glv * S * F * nr[ijk] / rho[k];

                    qrt[ijk]  += ev_tend;
                    nrt[ijk]  += lambda_evap * ev_tend * rho[k] / mr;
                    qtt[ijk]  -= ev_tend;
                    thlt[ijk] += Lv<TF> / (cp<TF> * exner[k]) * ev_tend;
                }
            }
    }

    // Selfcollection & breakup: growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
    template<typename TF>
    void selfcollection_breakup(TF* const restrict nrt, const TF* const restrict qr, const TF* const restrict nr, const TF* const restrict rho,
                                const TF* const restrict rain_mass, const TF* const restrict rain_diameter,
                                const TF* const restrict lambda_r,
                                const int istart, const int jstart, const int kstart,
                                const int iend,   const int jend,   const int kend,
                                const int jj, const int kk, const int j)
    {
        const TF k_rr     = 7.12;   // SB06, p49
        const TF kappa_rr = 60.7;   // SB06, p49
        const TF D_eq     = 0.9e-3; // SB06, list of symbols
        const TF k_br1    = 1.0e3;  // SB06, p50, for 0.35e-3 <= Dr <= D_eq
        const TF k_br2    = 2.3e3;  // SB06, p50, for Dr > D_eq

        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ik  = i + k*jj;
                const int ijk = i + j*jj + k*kk;

                if (qr[ijk] > qr_min<TF>)
                {
                    // Calculate mean rain drop mass and diameter
                    const TF dr      = rain_diameter[ik];
                    const TF lambdar = lambda_r[ik];

                    // Selfcollection
                    const TF sc_tend =
                        -k_rr * nr[ijk] * qr[ijk]*rho[k] / fm::pow9(TF(1.) + kappa_rr /
                        lambdar * pow(pirhow<TF>, TF(1.)/TF(3.))) * sqrt(rho_0<TF> / rho[k]);

                    nrt[ijk] += sc_tend;

                    // Breakup
                    const TF dDr = dr - D_eq;
                    if (dr > TF(0.35e-3))
                    {
                        TF phi_br;
                        if (dr <= D_eq)
                            phi_br = k_br1 * dDr;
                        else
                            phi_br = TF(2.) * exp(k_br2 * dDr) - TF(1.);

                        const TF br_tend = -(phi_br + TF(1.)) * sc_tend;
                        nrt[ijk] += br_tend;
                    }
                }
            }
    }

    // Sedimentation from Stevens and Seifert (2008)
    template<typename TF>
    void sedimentation_ss08(TF* const restrict qrt, TF* const restrict nrt, TF* const restrict rr_bot,
                            TF* const restrict w_qr, TF* const restrict w_nr,
                            TF* const restrict c_qr, TF* const restrict c_nr,
                            TF* const restrict slope_qr, TF* const restrict slope_nr,
                            TF* const restrict flux_qr, TF* const restrict flux_nr,
                            const TF* const restrict mu_r, const TF* const restrict lambda_r,
                            const TF* const restrict qr, const TF* const restrict nr,
                            const TF* const restrict rho, const TF* const restrict rhoh,
                            const TF* const restrict dzi,
                            const TF* const restrict dz, const double dt,
                            const int istart, const int jstart, const int kstart,
                            const int iend,   const int jend,   const int kend,
                            const int icells, const int kcells, const int ijcells, const int j)
    {
        const TF w_max = 9.65; // 9.65=UCLA, 20=SS08, appendix A
        const TF a_R = 9.65;   // SB06, p51
        const TF c_R = 600;    // SB06, p51
        const TF Dv  = 25.0e-6;
        const TF b_R = a_R * exp(c_R*Dv); // UCLA-LES

        const int kk3d = ijcells;
        const int kk2d = icells;

        // 1. Calculate sedimentation velocity at cell center
        for (int k=kstart; k<kend; k++)
        {
            const TF rho_n = sqrt(TF(1.2) / rho[k]);
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*icells + k*ijcells;
                const int ik  = i + k*icells;

                if (qr[ijk] > qr_min<TF>)
                {
                    // SS08:
                    w_qr[ik] = std::min(w_max, std::max(TF(0.1), rho_n * a_R - b_R * TF(pow(TF(1.) + c_R/lambda_r[ik], TF(-1.)*(mu_r[ik]+TF(4.))))));
                    w_nr[ik] = std::min(w_max, std::max(TF(0.1), rho_n * a_R - b_R * TF(pow(TF(1.) + c_R/lambda_r[ik], TF(-1.)*(mu_r[ik]+TF(1.))))));
                }
                else
                {
                    w_qr[ik] = 0.;
                    w_nr[ik] = 0.;
                }
            }
        }

        // 1.1 Set one ghost cell to zero
        for (int i=istart; i<iend; i++)
        {
            const int ik1  = i + (kstart-1)*icells;
            const int ik2  = i + (kend    )*icells;
            w_qr[ik1] = w_qr[ik1+kk2d];
            w_nr[ik1] = w_nr[ik1+kk2d];
            w_qr[ik2] = TF(0.);
            w_nr[ik2] = TF(0.);
        }

        // 2. Calculate CFL number using interpolated sedimentation velocity
        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ik  = i + k*icells;
                c_qr[ik] = TF(0.25) * (w_qr[ik-kk2d] + TF(2.)*w_qr[ik] + w_qr[ik+kk2d]) * dzi[k] * dt;
                c_nr[ik] = TF(0.25) * (w_nr[ik-kk2d] + TF(2.)*w_nr[ik] + w_nr[ik+kk2d]) * dzi[k] * dt;
            }

        // 3. Calculate slopes
        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*icells + k*ijcells;
                const int ik  = i + k*icells;

                slope_qr[ik] = minmod(qr[ijk]-qr[ijk-kk3d], qr[ijk+kk3d]-qr[ijk]);
                slope_nr[ik] = minmod(nr[ijk]-nr[ijk-kk3d], nr[ijk+kk3d]-nr[ijk]);
            }

        // Calculate flux
        // Set the fluxes at the top of the domain (kend) to zero
        for (int i=istart; i<iend; i++)
        {
            const int ik  = i + kend*icells;
            flux_qr[ik] = TF(0.);
            flux_nr[ik] = TF(0.);
        }

        for (int k=kend-1; k>kstart-1; k--)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*icells + k*ijcells;
                const int ik  = i + k*icells;

                int kk;
                TF ftot, dzz, cc;

                // q_rain
                kk    = k;  // current grid level
                ftot  = TF(0);  // cumulative 'flux' (kg m-2)
                dzz   = TF(0);  // distance from zh[k]
                cc    = std::min(TF(1), c_qr[ik]);
                while (cc > 0 && kk < kend)
                {
                    const int ikk  = i + kk*icells;
                    const int ijkk = i + j*icells + kk*ijcells;

                    ftot  += rho[kk] * (qr[ijkk] + TF(0.5) * slope_qr[ikk] * (TF(1.)-cc)) * cc * dz[kk];

                    dzz   += dz[kk];
                    kk    += 1;
                    cc     = std::min(TF(1.), c_qr[ikk] - dzz*dzi[kk]);
                }

                // Given flux at top, limit bottom flux such that the total rain content stays >= 0.
                ftot = std::min(ftot, rho[k] * dz[k] * qr[ijk] - flux_qr[ik+icells] * TF(dt));
                flux_qr[ik] = -ftot / dt;

                // number density
                kk    = k;  // current grid level
                ftot  = TF(0);  // cumulative 'flux'
                dzz   = TF(0);  // distance from zh[k]
                cc    = std::min(TF(1), c_nr[ik]);
                while (cc > 0 && kk < kend)
                {
                    const int ikk  = i + kk*icells;
                    const int ijkk = i + j*icells + kk*ijcells;

                    ftot += rho[kk] * (nr[ijkk] + TF(0.5) * slope_nr[ikk] * (TF(1.)-cc)) * cc * dz[kk];

                    dzz   += dz[kk];
                    kk    += 1;
                    cc     = std::min(TF(1.), c_nr[ikk] - dzz*dzi[k]);
                }

                // Given flux at top, limit bottom flux such that the number density stays >= 0.
                ftot = std::min(ftot, rho[k] * dz[k] * nr[ijk] - flux_nr[ik+icells] * TF(dt));
                flux_nr[ik] = -ftot / dt;
            }

        // Calculate tendency
        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*icells + k*ijcells;
                const int ik  = i + k*icells;

                qrt[ijk] += -(flux_qr[ik+kk2d] - flux_qr[ik]) / rho[k] * dzi[k];
                nrt[ijk] += -(flux_nr[ik+kk2d] - flux_nr[ik]) / rho[k] * dzi[k];
            }


        // Store surface sedimentation flux
        // Sedimentation flux is already multiplied with density (see flux div. calculation), so
        // the resulting flux is in kg m-2 s-1, with rho_water = 1000 kg/m3 this equals a rain rate in mm s-1
        for (int i=istart; i<iend; i++)
        {
            const int ij  = i + j*icells;
            const int ik  = i + kstart*icells;

            rr_bot[ij] = -flux_qr[ik];
        }
    }
}

template<typename TF>
Microphys_2mom_warm<TF>::Microphys_2mom_warm(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Microphys<TF>(masterin, gridin, fieldsin, inputin)
{
    auto& gd = grid.get_grid_data();
    swmicrophys = Microphys_type::Warm_2mom;

    // Read microphysics switches and settings
    swmicrobudget = inputin.get_item<bool>("micro", "swmicrobudget", "", false);
    cflmax = inputin.get_item<TF>("micro", "cflmax", "", 2.);
    Nc0 = inputin.get_item<TF>("micro", "Nc0", "");

    // Initialize the qr (rain water specific humidity) and nr (droplot number concentration) fields
    const std::string group_name = "thermo";

    fields.init_prognostic_field("qr", "Rain water specific humidity", "kg kg-1", group_name, gd.sloc, false);
    fields.init_prognostic_field("nr", "Number density rain", "m-3", group_name, gd.sloc, false);

    // Load the viscosity for both fields.
    fields.sp.at("qr")->visc = inputin.get_item<TF>("fields", "svisc", "qr");
    fields.sp.at("nr")->visc = inputin.get_item<TF>("fields", "svisc", "nr");
}

template<typename TF>
Microphys_2mom_warm<TF>::~Microphys_2mom_warm()
{
}

template<typename TF>
void Microphys_2mom_warm<TF>::init()
{
    auto& gd = grid.get_grid_data();

    rr_bot.resize(gd.ijcells);     // 2D surface sedimentation flux (rain rate)
}

template<typename TF>
void Microphys_2mom_warm<TF>::create(
        Input& inputin, Netcdf_handle& input_nc,
        Stats<TF>& stats, Cross<TF>& cross, Dump<TF>& dump, Column<TF>& column)
{
    const std::string group_name = "thermo";

    // BvS: for now I have left the init of statistics and cross-sections here
    // If this gets out of hand, move initialisation to separate function like in e.g. thermo_moist

    // Add variables to the statistics
    if (stats.get_switch())
    {
        // Time series
        stats.add_time_series("rr", "Mean surface rain rate", "kg m-2 s-1", group_name);
        stats.add_profs(*fields.sp.at("qr"), "z", {"frac", "cover"}, group_name);

        if (swmicrobudget)
        {
            // Microphysics tendencies for qr, nr, thl and qt
            stats.add_prof("sed_qrt", "Sedimentation tendency of qr", "kg kg-1 s-1", "z", group_name);
            stats.add_prof("sed_nrt", "Sedimentation tendency of nr", "m3 s-1", "z", group_name);

            stats.add_prof("auto_qrt" , "Autoconversion tendency qr",  "kg kg-1 s-1", "z", group_name);
            stats.add_prof("auto_nrt" , "Autoconversion tendency nr",  "m-3 s-1", "z", group_name);
            stats.add_prof("auto_thlt", "Autoconversion tendency thl", "K s-1", "z", group_name);
            stats.add_prof("auto_qtt" , "Autoconversion tendency qt",  "kg kg-1 s-1", "z", group_name);

            stats.add_prof("evap_qrt" , "Evaporation tendency qr",  "kg kg-1 s-1", "z", group_name);
            stats.add_prof("evap_nrt" , "Evaporation tendency nr",  "m-3 s-1", "z", group_name);
            stats.add_prof("evap_thlt", "Evaporation tendency thl", "K s-1", "z", group_name);
            stats.add_prof("evap_qtt" , "Evaporation tendency qt",  "kg kg-1 s-1", "z", group_name);

            stats.add_prof("scbr_nrt" , "Selfcollection and breakup tendency nr", "m-3 s-1", "z", group_name);

            stats.add_prof("accr_qrt" , "Accretion tendency qr",  "kg kg-1 s-1", "z", group_name);
            stats.add_prof("accr_thlt", "Accretion tendency thl", "K s-1", "z", group_name);
            stats.add_prof("accr_qtt" , "Accretion tendency qt",  "kg kg-1 s-1", "z", group_name);
        }

        stats.add_tendency(*fields.st.at("thl"), "z", tend_name, tend_longname);
        stats.add_tendency(*fields.st.at("qt") , "z", tend_name, tend_longname);
        stats.add_tendency(*fields.st.at("qr") , "z", tend_name, tend_longname);
        stats.add_tendency(*fields.st.at("nr") , "z", tend_name, tend_longname);
    }

    // Add variables to column statistics
    if (column.get_switch())
    {
        column.add_time_series("rr", "Surface rain rate", "kg m-2 s-1");
    }

    // Create cross sections
    // 1. Variables that this class can calculate/provide:
    std::vector<std::string> allowed_crossvars = {"rr_bot"};
    // 2. Cross-reference with the variables requested in the .ini file:
    crosslist = cross.get_enabled_variables(allowed_crossvars);
}

#ifndef USECUDA
template<typename TF>
void Microphys_2mom_warm<TF>::exec(Thermo<TF>& thermo, const double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    // Remove spurious negative values from qr and nr fields
    remove_negative_values(fields.sp.at("qr")->fld.data(), gd.istart, gd.jstart, gd.kstart,
                           gd.iend, gd.jend, gd.kend, gd.icells, gd.ijcells);
    remove_negative_values(fields.sp.at("nr")->fld.data(), gd.istart, gd.jstart, gd.kstart,
                           gd.iend, gd.jend, gd.kend, gd.icells, gd.ijcells);

    // Get cloud liquid water specific humidity from thermodynamics
    auto ql = fields.get_tmp();
    thermo.get_thermo_field(*ql, "qlqi", false, false);

    // Get pressure and exner function from thermodynamics
    std::vector<TF> p     = thermo.get_basestate_vector("p");
    std::vector<TF> exner = thermo.get_basestate_vector("exner");

    // Microphysics is handled in XZ slices, to
    // (1) limit the required number of tmp fields
    // (2) re-use some expensive calculations used in multiple microphysics routines.
    const int ikcells    = gd.icells * gd.kcells;                           // Size of XZ slice
    const int n_slices   = 12;                                              // Number of XZ slices required
    const int n_tmp_flds = std::ceil(static_cast<TF>(n_slices)/gd.jcells);  // Number of required tmp fields

    // Load the required number of tmp fields:
    std::vector<std::shared_ptr<Field3d<TF>>> tmp_fields;
    for (int n=0; n<n_tmp_flds; ++n)
        tmp_fields.push_back(fields.get_tmp());

    // Get pointers to slices in tmp fields:
    int slice_counter = 0;

    TF* w_qr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
    TF* w_nr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

    TF* c_qr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
    TF* c_nr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

    TF* slope_qr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
    TF* slope_nr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

    TF* flux_qr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
    TF* flux_nr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

    TF* rain_mass = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
    TF* rain_diam = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

    TF* lambda_r = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
    TF* mu_r     = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

    // ---------------------------------
    // Calculate microphysics tendencies
    // ---------------------------------

    // Autoconversion; formation of rain drop by coagulating cloud droplets
    mp3d::autoconversion(fields.st.at("qr")->fld.data(), fields.st.at("nr")->fld.data(), fields.st.at("qt")->fld.data(), fields.st.at("thl")->fld.data(),
                         fields.sp.at("qr")->fld.data(), ql->fld.data(), fields.rhoref.data(), exner.data(), Nc0,
                         gd.istart, gd.jstart, gd.kstart,
                         gd.iend,   gd.jend,   gd.kend,
                         gd.icells, gd.ijcells);

    // Accretion; growth of raindrops collecting cloud droplets
    mp3d::accretion(fields.st.at("qr")->fld.data(), fields.st.at("qt")->fld.data(), fields.st.at("thl")->fld.data(),
                    fields.sp.at("qr")->fld.data(), ql->fld.data(), fields.rhoref.data(), exner.data(),
                    gd.istart, gd.jstart, gd.kstart,
                    gd.iend,   gd.jend,   gd.kend,
                    gd.icells, gd.ijcells);

    // Rest of the microphysics is handled per XZ slice
    for (int j=gd.jstart; j<gd.jend; ++j)
    {
        // Prepare the XZ slices which are used in all routines
        mp2d::prepare_microphysics_slice(rain_mass, rain_diam, mu_r, lambda_r,
                                         fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(), fields.rhoref.data(),
                                         gd.istart, gd.iend, gd.kstart, gd.kend, gd.icells, gd.ijcells, j);

        // Evaporation; evaporation of rain drops in unsaturated environment
        mp2d::evaporation(fields.st.at("qr")->fld.data(), fields.st.at("nr")->fld.data(),  fields.st.at("qt")->fld.data(), fields.st.at("thl")->fld.data(),
                          fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(),  ql->fld.data(),
                          fields.sp.at("qt")->fld.data(), fields.sp.at("thl")->fld.data(), fields.rhoref.data(), exner.data(), p.data(),
                          rain_mass, rain_diam,
                          gd.istart, gd.jstart, gd.kstart,
                          gd.iend,   gd.jend,   gd.kend,
                          gd.icells, gd.ijcells, j);

        // Self collection and breakup; growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
        mp2d::selfcollection_breakup(fields.st.at("nr")->fld.data(), fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(), fields.rhoref.data(),
                                     rain_mass, rain_diam, lambda_r,
                                     gd.istart, gd.jstart, gd.kstart,
                                     gd.iend,   gd.jend,   gd.kend,
                                     gd.icells, gd.ijcells, j);

        // Sedimentation; sub-grid sedimentation of rain
        mp2d::sedimentation_ss08(fields.st.at("qr")->fld.data(), fields.st.at("nr")->fld.data(), rr_bot.data(),
                                 w_qr, w_nr, c_qr, c_nr, slope_qr, slope_nr, flux_qr, flux_nr, mu_r, lambda_r,
                                 fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(),
                                 fields.rhoref.data(), fields.rhorefh.data(), gd.dzi.data(), gd.dz.data(), dt,
                                 gd.istart, gd.jstart, gd.kstart,
                                 gd.iend,   gd.jend,   gd.kend,
                                 gd.icells, gd.kcells, gd.ijcells, j);
    }

    // Release all local tmp fields in use
    for (auto& it: tmp_fields)
        fields.release_tmp(it);

    fields.release_tmp(ql);

    stats.calc_tend(*fields.st.at("thl"), tend_name);
    stats.calc_tend(*fields.st.at("qt"),  tend_name);
    stats.calc_tend(*fields.st.at("qr"),  tend_name);
    stats.calc_tend(*fields.st.at("nr"),  tend_name);
}
#endif

template<typename TF>
void Microphys_2mom_warm<TF>::exec_stats(Stats<TF>& stats, Thermo<TF>& thermo, const double dt)
{
    auto& gd = grid.get_grid_data();

    // Time series
    const TF no_offset = 0.;
    const TF no_threshold = 0.;
    const TF threshold_qr = 1.e-6;

    stats.calc_stats_2d("rr", rr_bot, no_offset);
    stats.calc_stats("qr", *fields.sp.at("qr"), no_offset, threshold_qr);

    if (swmicrobudget)
    {
        // Vertical profiles. The statistics of qr & nr are handled by fields.cxx
        // Get cloud liquid water specific humidity from thermodynamics
        auto ql = fields.get_tmp();
        ql->loc = gd.sloc;
        thermo.get_thermo_field(*ql, "ql", false, false);

        // Get pressure and exner function from thermodynamics
        std::vector<TF> p     = thermo.get_basestate_vector("p");
        std::vector<TF> exner = thermo.get_basestate_vector("exner");

        // Microphysics is (partially) handled in XZ slices, to
        // (1) limit the required number of tmp fields
        // (2) re-use some expensive calculations used in multiple microphysics routines.
        const int ikcells    = gd.icells * gd.kcells;                           // Size of XZ slice
        const int n_slices   = 12;                                              // Number of XZ slices required
        const int n_tmp_flds = std::ceil(static_cast<TF>(n_slices)/gd.jcells);  // Number of required tmp fields

        // Load the required number of tmp fields:
        std::vector<std::shared_ptr<Field3d<TF>>> tmp_fields;
        for (int n=0; n<n_tmp_flds; ++n)
            tmp_fields.push_back(fields.get_tmp());

        // Get pointers to slices in tmp fields:
        int slice_counter = 0;

        TF* w_qr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
        TF* w_nr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

        TF* c_qr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
        TF* c_nr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

        TF* slope_qr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
        TF* slope_nr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

        TF* flux_qr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
        TF* flux_nr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

        TF* rain_mass = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
        TF* rain_diam = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

        TF* lambda_r = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
        TF* mu_r     = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

        // Get 4 tmp fields for all tendencies (qrt, nrt, thlt, qtt) :-(
        auto qrt  = fields.get_tmp();
        auto nrt  = fields.get_tmp();
        auto thlt = fields.get_tmp();
        auto qtt  = fields.get_tmp();
        qrt->loc  = gd.sloc;
        nrt->loc  = gd.sloc;
        thlt->loc = gd.sloc;
        qtt->loc  = gd.sloc;

        // Calculate tendencies
        // Autoconversion; formation of rain drop by coagulating cloud droplets
        // -------------------------
        zero_field(qrt->fld.data(),  gd.ncells);
        zero_field(nrt->fld.data(),  gd.ncells);
        zero_field(thlt->fld.data(), gd.ncells);
        zero_field(qtt->fld.data(),  gd.ncells);

        mp3d::autoconversion(qrt->fld.data(), nrt->fld.data(), qtt->fld.data(), thlt->fld.data(),
                             fields.sp.at("qr")->fld.data(), ql->fld.data(), fields.rhoref.data(), exner.data(), Nc0,
                             gd.istart, gd.jstart, gd.kstart,
                             gd.iend,   gd.jend,   gd.kend,
                             gd.icells, gd.ijcells);

        stats.calc_stats("auto_qrt" , *qrt , no_offset, no_threshold);
        stats.calc_stats("auto_nrt" , *nrt , no_offset, no_threshold);
        stats.calc_stats("auto_thlt", *thlt, no_offset, no_threshold);
        stats.calc_stats("auto_qtt" , *qtt , no_offset, no_threshold);

        // Accretion; growth of raindrops collecting cloud droplets
        // -------------------------
        zero_field(qrt->fld.data(),  gd.ncells);
        zero_field(thlt->fld.data(), gd.ncells);
        zero_field(qtt->fld.data(),  gd.ncells);

        mp3d::accretion(qrt->fld.data(), qtt->fld.data(), thlt->fld.data(),
                        fields.sp.at("qr")->fld.data(), ql->fld.data(), fields.rhoref.data(), exner.data(),
                        gd.istart, gd.jstart, gd.kstart,
                        gd.iend,   gd.jend,   gd.kend,
                        gd.icells, gd.ijcells);

        stats.calc_stats("accr_qrt" , *qrt , no_offset, no_threshold);
        stats.calc_stats("accr_thlt", *thlt, no_offset, no_threshold);
        stats.calc_stats("accr_qtt" , *qtt , no_offset, no_threshold);

        // Rest of the microphysics is handled per XZ slice
        // Evaporation; evaporation of rain drops in unsaturated environment
        // -------------------------
        zero_field(qrt->fld.data(),  gd.ncells);
        zero_field(nrt->fld.data(),  gd.ncells);
        zero_field(thlt->fld.data(), gd.ncells);
        zero_field(qtt->fld.data(),  gd.ncells);

        for (int j=gd.jstart; j<gd.jend; ++j)
        {
            mp2d::prepare_microphysics_slice(rain_mass, rain_diam, mu_r, lambda_r,
                                             fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(), fields.rhoref.data(),
                                             gd.istart, gd.iend, gd.kstart, gd.kend, gd.icells, gd.ijcells, j);

            mp2d::evaporation(qrt->fld.data(), nrt->fld.data(),  qtt->fld.data(), thlt->fld.data(),
                              fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(),  ql->fld.data(),
                              fields.sp.at("qt")->fld.data(), fields.sp.at("thl")->fld.data(), fields.rhoref.data(), exner.data(), p.data(),
                              rain_mass, rain_diam,
                              gd.istart, gd.jstart, gd.kstart,
                              gd.iend,   gd.jend,   gd.kend,
                              gd.icells, gd.ijcells, j);
        }

        stats.calc_stats("evap_qrt" , *qrt , no_offset, no_threshold);
        stats.calc_stats("evap_nrt" , *nrt , no_offset, no_threshold);
        stats.calc_stats("evap_thlt", *thlt, no_offset, no_threshold);
        stats.calc_stats("evap_qtt" , *qtt , no_offset, no_threshold);

        // Self collection and breakup; growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
        // -------------------------
        zero_field(nrt->fld.data(),  gd.ncells);

        for (int j=gd.jstart; j<gd.jend; ++j)
        {
            mp2d::prepare_microphysics_slice(rain_mass, rain_diam, mu_r, lambda_r,
                                             fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(), fields.rhoref.data(),
                                             gd.istart, gd.iend, gd.kstart, gd.kend, gd.icells, gd.ijcells, j);

            mp2d::selfcollection_breakup(nrt->fld.data(), fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(), fields.rhoref.data(),
                                         rain_mass, rain_diam, lambda_r,
                                         gd.istart, gd.jstart, gd.kstart,
                                         gd.iend,   gd.jend,   gd.kend,
                                         gd.icells, gd.ijcells, j);
        }

        stats.calc_stats("scbr_nrt" , *nrt , no_offset, no_threshold);

        // Sedimentation; sub-grid sedimentation of rain
        // -------------------------
        zero_field(qrt->fld.data(),  gd.ncells);
        zero_field(nrt->fld.data(),  gd.ncells);

        for (int j=gd.jstart; j<gd.jend; ++j)
        {
            mp2d::prepare_microphysics_slice(rain_mass, rain_diam, mu_r, lambda_r,
                                             fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(), fields.rhoref.data(),
                                             gd.istart, gd.iend, gd.kstart, gd.kend, gd.icells, gd.ijcells, j);

            mp2d::sedimentation_ss08(qrt->fld.data(), nrt->fld.data(), rr_bot.data(),
                                     w_qr, w_nr, c_qr, c_nr, slope_qr, slope_nr, flux_qr, flux_nr, mu_r, lambda_r,
                                     fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(),
                                     fields.rhoref.data(), fields.rhorefh.data(), gd.dzi.data(), gd.dz.data(), dt,
                                     gd.istart, gd.jstart, gd.kstart,
                                     gd.iend,   gd.jend,   gd.kend,
                                     gd.icells, gd.kcells, gd.ijcells, j);
        }

        stats.calc_stats("sed_qrt" , *qrt , no_offset, no_threshold);
        stats.calc_stats("sed_nrt" , *nrt , no_offset, no_threshold);

        // Release all local tmp fields in use
        for (auto& it: tmp_fields)
            fields.release_tmp(it);

        fields.release_tmp(ql  );
        fields.release_tmp(qrt );
        fields.release_tmp(nrt );
        fields.release_tmp(thlt);
        fields.release_tmp(qtt );
    }
}

#ifndef USECUDA
template<typename TF>
void Microphys_2mom_warm<TF>::exec_column(Column<TF>& column)
{
    const TF no_offset = 0.;
    column.calc_time_series("rr", rr_bot.data(), no_offset);
}
#endif

template<typename TF>
void Microphys_2mom_warm<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    TF no_offset = 0.;
    if (cross.get_switch())
    {
        for (auto& it : crosslist)
        {
            if (it == "rr_bot")
                cross.cross_plane(rr_bot.data(), no_offset, "rr_bot", iotime);
        }
    }
}

#ifndef USECUDA
template<typename TF>
unsigned long Microphys_2mom_warm<TF>::get_time_limit(unsigned long idt, const double dt)
{
    auto& gd = grid.get_grid_data();

    // Calculate the maximum sedimentation CFL number
    auto w_qr = fields.get_tmp();
    TF cfl = mp3d::calc_max_sedimentation_cfl(w_qr->fld.data(), fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(),
                                              fields.rhoref.data(), gd.dzi.data(), dt,
                                              gd.istart, gd.jstart, gd.kstart,
                                              gd.iend,   gd.jend,   gd.kend,
                                              gd.icells, gd.ijcells);
    fields.release_tmp(w_qr);

    // Get maximum CFL across all MPI tasks
    master.max(&cfl, 1);

    return idt * cflmax / cfl;
}
#endif

template<typename TF>
bool Microphys_2mom_warm<TF>::has_mask(std::string name)
{
    if (std::find(available_masks.begin(), available_masks.end(), name) != available_masks.end())
        return true;
    else
        return false;
}

template<typename TF>
void Microphys_2mom_warm<TF>::get_mask(Stats<TF>& stats, std::string mask_name)
{
    auto& gd = grid.get_grid_data();

    if (mask_name == "qr")
    {
        TF threshold = 1e-6;

        // Interpolate qr to half level:
        auto qrh = fields.get_tmp();
        grid.interpolate_2nd(qrh->fld.data(), fields.sp.at("qr")->fld.data(), gd.sloc.data(), gd.wloc.data());

        // Calculate masks
        stats.set_mask_thres(mask_name, *fields.sp.at("qr"), *qrh, threshold, Stats_mask_type::Plus);

        fields.release_tmp(qrh);
    }
    else
    {
        std::string message = "Double moment warm microphysics can not provide mask: \"" + mask_name +"\"";
        throw std::runtime_error(message);
    }
}

template<typename TF>
void Microphys_2mom_warm<TF>::get_surface_rain_rate(std::vector<TF>& field)
{
    // Make a hard copy of the surface precipitation field
    field = rr_bot;
}


#ifdef FLOAT_SINGLE
template class Microphys_2mom_warm<float>;
#else
template class Microphys_2mom_warm<double>;
#endif
