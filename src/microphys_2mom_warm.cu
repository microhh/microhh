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

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"

#include "boundary_cyclic.h"
#include "thermo_moist_functions.h"
#include "constants.h"
#include "tools.h"
#include "field3d_operators.h"
#include "stats.h"
#include "column.h"
#include "microphys.h"
#include "microphys_2mom_warm.h"
#include "microphys_sedi_kernels.h"

using namespace Constants;
using namespace Thermo_moist_functions;
using namespace Micro_2mom_warm_constants;
using namespace Micro_2mom_warm_functions;

namespace micro
{
    template<typename TF> __global__
    void remove_negative_values_g(TF* __restrict__ field,
                                  int istart, int jstart, int kstart,
                                  int iend,   int jend,   int kend,
                                  int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            field[ijk] = fmax(field[ijk], TF(0));
        }
    }

    // Autoconversion: formation of rain drop by coagulating cloud droplets
    template<typename TF> __global__
    void autoconversion_g(TF* const __restrict__ qrt, TF* const __restrict__ nrt,
                          TF* const __restrict__ qtt, TF* const __restrict__ thlt,
                          const TF* const __restrict__ qr,  const TF* const __restrict__ ql,
                          const TF* const __restrict__ rho, const TF* const __restrict__ exner, const TF nc,
                          const int istart, const int jstart, const int kstart,
                          const int iend,   const int jend,   const int kend,
                          const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        // BvS with the cpu and cuda implementation, these should really go to the header
        const TF x_star = 2.6e-10;       // SB06, list of symbols, same as UCLA-LES
        const TF k_cc   = 9.44e9;        // UCLA-LES (Long, 1974), 4.44e9 in SB06, p48
        const TF nu_c   = 1;             // SB06, Table 1., same as UCLA-LES
        const TF kccxs  = k_cc / (TF(20.) * x_star) * (nu_c+2)*(nu_c+4) / pow(nu_c+1, TF(2));

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            if (ql[ijk] > ql_min<TF>)
            {
                const TF xc      = rho[k] * ql[ijk] / nc;    // Mean mass of cloud drops [kg]
                const TF tau     = TF(1.) - ql[ijk] / (ql[ijk] + qr[ijk] + dsmall);    // SB06, Eq 5
                const TF phi_au  = TF(600.) * pow(tau, TF(0.68)) * pow(TF(1.) - pow(tau, TF(0.68)), TF(3));    // UCLA-LES
                //const TF phi_au  = 400. * pow(tau, 0.7) * pow(1. - pow(tau, 0.7), 3);    // SB06, Eq 6
                const TF au_tend = rho_0<TF> * kccxs * pow(ql[ijk], TF(2)) * pow(xc, TF(2)) *
                                       (TF(1.) + phi_au / pow(TF(1.)-tau, TF(2))); // SB06, eq 4

                qrt[ijk]  += au_tend;
                nrt[ijk]  += au_tend * rho[k] / x_star;
                qtt[ijk]  -= au_tend;
                thlt[ijk] += Lv<TF> / (cp<TF> * exner[k]) * au_tend;
            }
        }
    }

    // Accreation: growth of raindrops collecting cloud droplets
    template<typename TF> __global__
    void accretion_g(TF* const __restrict__ qrt, TF* const __restrict__ qtt, TF* const __restrict__ thlt,
                     const TF* const __restrict__ qr,  const TF* const __restrict__ ql,
                     const TF* const __restrict__ rho, const TF* const __restrict__ exner,
                     const int istart, const int jstart, const int kstart,
                     const int iend,   const int jend,   const int kend,
                     const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const TF k_cr  = 5.25; // SB06, p49

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            if (ql[ijk] > ql_min<TF> && qr[ijk] > qr_min<TF>)
            {
                const TF tau     = TF(1.) - ql[ijk] / (ql[ijk] + qr[ijk]); // SB06, Eq 5
                const TF phi_ac  = pow(tau / (tau + TF(5e-5)), TF(4)); // SB06, Eq 8
                const TF ac_tend = k_cr * ql[ijk] *  qr[ijk] * phi_ac * pow(rho_0<TF> / rho[k], TF(0.5)); // SB06, Eq 7

                qrt[ijk]  += ac_tend;
                qtt[ijk]  -= ac_tend;
                thlt[ijk] += Lv<TF> / (cp<TF> * exner[k]) * ac_tend;
            }
        }
    }

    // Evaporation: evaporation of rain drops in unsaturated environment
    template<typename TF> __global__
    void evaporation_g(TF* const restrict qrt, TF* const restrict nrt,
                       TF* const restrict qtt, TF* const restrict thlt,
                       const TF* const restrict qr, const TF* const restrict nr,
                       const TF* const restrict ql, const TF* const restrict qt, const TF* const restrict thl,
                       const TF* const restrict rho, const TF* const restrict exner, const TF* const restrict p,
                       const int istart, const int jstart, const int kstart,
                       const int iend,   const int jend,   const int kend,
                       const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const TF lambda_evap = 1.; // 1.0 in UCLA, 0.7 in DALES

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            if (qr[ijk] > qr_min<TF>)
            {
                const TF mr = calc_rain_mass(qr[ijk], nr[ijk], rho[k]);
                const TF dr = calc_rain_diameter(mr);

                const TF T   = thl[ijk] * exner[k] + (Lv<TF> * ql[ijk]) / (cp<TF> * exner[k]); // Absolute temperature [K]
                const TF Glv = pow(Rv<TF> * T / (esat_liq(T) * D_v<TF>) +
                                   (Lv<TF> / (K_t<TF> * T)) * (Lv<TF> / (Rv<TF> * T) - TF(1)), TF(-1)); // Cond/evap rate (kg m-1 s-1)?

                const TF S   = (qt[ijk] - ql[ijk]) / qsat_liq(p[k], T) - TF(1); // Saturation
                const TF F   = TF(1.); // Evaporation excludes ventilation term from SB06 (like UCLA, unimportant term? TODO: test)

                const TF ev_tend = TF(2.) * pi<TF> * dr * Glv * S * F * nr[ijk] / rho[k];

                qrt[ijk]  += ev_tend;
                nrt[ijk]  += lambda_evap * ev_tend * rho[k] / mr;
                qtt[ijk]  -= ev_tend;
                thlt[ijk] += Lv<TF> / (cp<TF> * exner[k]) * ev_tend;
            }
        }
    }

    // Selfcollection & breakup: growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
    template<typename TF> __global__
    void selfcollection_breakup_g(TF* const __restrict__ nrt, const TF* const __restrict__ qr, const TF* const __restrict__ nr, const TF* const __restrict__ rho,
                                  const int istart, const int jstart, const int kstart,
                                  const int iend,   const int jend,   const int kend,
                                  const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const TF k_rr     = 7.12;   // SB06, p49
        const TF kappa_rr = 60.7;   // SB06, p49
        const TF D_eq     = 0.9e-3; // SB06, list of symbols
        const TF k_br1    = 1.0e3;  // SB06, p50, for 0.35e-3 <= Dr <= D_eq
        const TF k_br2    = 2.3e3;  // SB06, p50, for Dr > D_eq

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            if (qr[ijk] > qr_min<TF>)
            {
                // Calculate mean rain drop mass and diameter
                const TF mr      = calc_rain_mass(qr[ijk], nr[ijk], rho[k]);
                const TF dr      = calc_rain_diameter(mr);
                const TF mu_r    = calc_mu_r(dr);
                const TF lambdar = calc_lambda_r(mu_r, dr);

                // Selfcollection
                const TF sc_tend = -k_rr * nr[ijk] * qr[ijk]*rho[k] * pow(TF(1.) + kappa_rr /
                                   lambdar * pow(pirhow<TF>, TF(1.)/TF(3.)), TF(-9)) * pow(rho_0<TF> / rho[k], TF(0.5));
                nrt[ijk] += sc_tend;

                // Breakup
                const TF dDr = dr - D_eq;
                if(dr > 0.35e-3)
                {
                    TF phi_br;
                    if(dr <= D_eq)
                        phi_br = k_br1 * dDr;
                    else
                        phi_br = TF(2.) * exp(k_br2 * dDr) - TF(1.);

                    const TF br_tend = -(phi_br + TF(1.)) * sc_tend;
                    nrt[ijk] += br_tend;
                }
            }
        }
    }
}

namespace sedimentation_2mom_warm
{
    // Sedimentation from Stevens and Seifert (2008)
    template<typename TF> __global__
    void calc_velocity_g(TF* const __restrict__ w_qr, TF* const __restrict__ w_nr,
                         const TF* const __restrict__ qr, const TF* const __restrict__ nr,
                         const TF* const __restrict__ rho,
                         const int istart, const int jstart, const int kstart,
                         const int iend,   const int jend,   const int kend,
                         const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const TF w_max = 9.65; // 9.65=UCLA, 20=SS08, appendix A
        const TF a_R   = 9.65;   // SB06, p51
        const TF c_R   = 600;    // SB06, p51
        const TF Dv    = 25.0e-6;
        const TF b_R   = a_R * exp(c_R*Dv); // UCLA-LES

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*icells + k*ijcells;

            if (qr[ijk] > qr_min<TF>)
            {
                // SS08:
                const TF rho_n   = pow(TF(1.2) / rho[k], TF(0.5));
                const TF mr      = calc_rain_mass(qr[ijk], nr[ijk], rho[k]);
                const TF dr      = calc_rain_diameter(mr);
                const TF mu_r    = calc_mu_r(dr);
                const TF lambdar = calc_lambda_r(mu_r, dr);

                w_qr[ijk] = fmin(w_max, fmax(TF(0.1), rho_n * a_R - b_R * TF(pow(TF(1.) + c_R/lambdar, TF(-1.)*(mu_r+TF(4.))))));
                w_nr[ijk] = fmin(w_max, fmax(TF(0.1), rho_n * a_R - b_R * TF(pow(TF(1.) + c_R/lambdar, TF(-1.)*(mu_r+TF(1.))))));
            }
            else
            {
                w_qr[ijk] = TF(0.);
                w_nr[ijk] = TF(0.);
            }
        }
    }

    // Sedimentation from Stevens and Seifert (2008)
    // Overloaded version for only w_qr (could also be done with templates, but then you
    // have to pass something for w_nr...
    template<typename TF> __global__
    void calc_velocity_g(TF* const __restrict__ w_qr,
                         const TF* const __restrict__ qr, const TF* const __restrict__ nr,
                         const TF* const __restrict__ rho,
                         const int istart, const int jstart, const int kstart,
                         const int iend,   const int jend,   const int kend,
                         const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const TF w_max = 9.65; // 9.65=UCLA, 20=SS08, appendix A
        const TF a_R   = 9.65;   // SB06, p51
        const TF c_R   = 600;    // SB06, p51
        const TF Dv    = 25.0e-6;
        const TF b_R   = a_R * exp(c_R*Dv); // UCLA-LES

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*icells + k*ijcells;

            if (qr[ijk] > qr_min<TF>)
            {
                // SS08:
                const TF rho_n   = pow(TF(1.2) / rho[k], TF(0.5));
                const TF mr      = calc_rain_mass(qr[ijk], nr[ijk], rho[k]);
                const TF dr      = calc_rain_diameter(mr);
                const TF mu_r    = calc_mu_r(dr);
                const TF lambdar = calc_lambda_r(mu_r, dr);

                w_qr[ijk] = fmin(w_max, fmax(TF(0.1), rho_n * a_R - b_R * TF(pow(TF(1.) + c_R/lambdar, TF(-1.)*(mu_r+TF(4.))))));
            }
            else
            {
                w_qr[ijk] = TF(0.);
            }
        }
    }

    template<typename TF> __global__
    void copy_rain_rate(
        TF* const __restrict__ rr_out,
        const TF* const __restrict__ rr_in,
        const int istart, const int iend,
        const int jstart, const int jend,
        const int icells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij = i + j*icells;
            rr_out[ij] = rr_in[ij];
        }
    }
}

#ifdef USECUDA
template<typename TF>
void Microphys_2mom_warm<TF>::exec(Thermo<TF>& thermo, const double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax+1);
    dim3 blockGPU(blocki, blockj, 1);

    dim3 grid2dGPU (gridi, gridj);
    dim3 block2dGPU(blocki, blockj);

    // Remove negative values from the qr and nr fields
    micro::remove_negative_values_g<<<gridGPU, blockGPU>>>(fields.ap.at("qr")->fld_g,
        gd.istart, gd.jstart, gd.kstart, gd.iend, gd.jend, gd.kend, gd.icells, gd.ijcells);
    cuda_check_error();

    micro::remove_negative_values_g<<<gridGPU, blockGPU>>>(fields.ap.at("nr")->fld_g,
        gd.istart, gd.jstart, gd.kstart, gd.iend, gd.jend, gd.kend, gd.icells, gd.ijcells);
    cuda_check_error();

    // Get cloud liquid water from thermodynamics
    auto ql = fields.get_tmp_g();
    thermo.get_thermo_field_g(*ql, "ql_qi", false);

    // Get GPU pointers basestate pressure and exner from thermo
    TF* p = thermo.get_basestate_fld_g("pref");
    TF* exner = thermo.get_basestate_fld_g("exner");

    // ---------------------------------
    // Calculate microphysics tendencies
    // ---------------------------------

    // Autoconversion; formation of rain drop by coagulating cloud droplets
    micro::autoconversion_g<<<gridGPU, blockGPU>>>(
        fields.st.at("qr")->fld_g, fields.st.at("nr")->fld_g,
        fields.st.at("qt")->fld_g, fields.st.at("thl")->fld_g,
        fields.sp.at("qr")->fld_g, ql->fld_g, fields.rhoref_g, exner, Nc0<TF>,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // Accretion; growth of raindrops collecting cloud droplets
    micro::accretion_g<<<gridGPU, blockGPU>>>(
        fields.st.at("qr")->fld_g, fields.st.at("qt")->fld_g, fields.st.at("thl")->fld_g,
        fields.sp.at("qr")->fld_g, ql->fld_g, fields.rhoref_g, exner,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // Evaporation; evaporation of rain drops in unsaturated environment
    micro::evaporation_g<<<gridGPU, blockGPU>>>(
        fields.st.at("qr")->fld_g, fields.st.at("nr")->fld_g,  fields.st.at("qt")->fld_g, fields.st.at("thl")->fld_g,
        fields.sp.at("qr")->fld_g, fields.sp.at("nr")->fld_g,  ql->fld_g,
        fields.sp.at("qt")->fld_g, fields.sp.at("thl")->fld_g, fields.rhoref_g, exner, p,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    fields.release_tmp_g(ql);

    // Self collection and breakup; growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
    micro::selfcollection_breakup_g<<<gridGPU, blockGPU>>>(
        fields.st.at("nr")->fld_g, fields.sp.at("qr")->fld_g, fields.sp.at("nr")->fld_g, fields.rhoref_g,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // Sedimentation; sub-grid sedimentation of rain
    // ---------------------------------------------
    // 1. Calculate sedimentation velocity of qr and nr
    auto w_qr = fields.get_tmp_g();
    auto w_nr = fields.get_tmp_g();

    sedimentation_2mom_warm::calc_velocity_g<<<gridGPU, blockGPU>>>(
        w_qr->fld_g, w_nr->fld_g,
        fields.sp.at("qr")->fld_g, fields.sp.at("nr")->fld_g, fields.rhoref_g,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // 2. Set the top and bottom velocity ghost cells
    Micro_sedimentation_kernels::set_bc_g<<<grid2dGPU, block2dGPU>>>(
        w_qr->fld_g,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    Micro_sedimentation_kernels::set_bc_g<<<grid2dGPU, block2dGPU>>>(
        w_nr->fld_g,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // qr --------------------
    // 3. Calculate CFL numbers for qr and nr
    auto cfl_qr = fields.get_tmp_g();

    Micro_sedimentation_kernels::calc_cfl_g<<<gridGPU, blockGPU>>>(
        cfl_qr->fld_g, w_qr->fld_g, gd.dzi_g, TF(dt),
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    fields.release_tmp_g(w_qr);

    // 4. Calculate slopes
    auto slope_qr = fields.get_tmp_g();

    Micro_sedimentation_kernels::calc_slope_g<<<gridGPU, blockGPU>>>(
        slope_qr->fld_g, fields.sp.at("qr")->fld_g,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // 5. Calculate sedimentation fluxes
    auto flux_qr = fields.get_tmp_g();

    Micro_sedimentation_kernels::calc_flux_g<<<gridGPU, blockGPU>>>(
        flux_qr->fld_g, fields.sp.at("qr")->fld_g,
        slope_qr->fld_g, cfl_qr->fld_g,
        gd.dz_g, gd.dzi_g, fields.rhoref_g,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // 6. Limit fluxes such that qr stays >= 0.
    Micro_sedimentation_kernels::limit_flux_g<<<grid2dGPU, block2dGPU>>>(
        flux_qr->fld_g, fields.sp.at("qr")->fld_g,
        gd.dz_g, fields.rhoref_g, TF(dt),
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // 7. Calculate flux divergence
    Micro_sedimentation_kernels::calc_flux_div_g<<<gridGPU, blockGPU>>>(
        fields.st.at("qr")->fld_g, flux_qr->fld_g,
        gd.dzi_g, fields.rhoref_g,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // 8. Store the surface precipitation rate
    Micro_sedimentation_kernels::copy_precip_rate_bot_g<<<grid2dGPU, block2dGPU>>>(
        rr_bot_g, flux_qr->fld_g,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart,
        gd.icells, gd.ijcells);
    cuda_check_error();

    fields.release_tmp_g(cfl_qr);
    fields.release_tmp_g(slope_qr);
    fields.release_tmp_g(flux_qr);

    // nr --------------------
    // 3. Calculate CFL numbers for qr and nr
    auto cfl_nr = fields.get_tmp_g();

    Micro_sedimentation_kernels::calc_cfl_g<<<gridGPU, blockGPU>>>(
        cfl_nr->fld_g, w_nr->fld_g, gd.dzi_g, TF(dt),
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    fields.release_tmp_g(w_nr);

    // 4. Calculate slopes
    auto slope_nr = fields.get_tmp_g();

    Micro_sedimentation_kernels::calc_slope_g<<<gridGPU, blockGPU>>>(
        slope_nr->fld_g, fields.sp.at("nr")->fld_g,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // 5. Calculate sedimentation fluxes
    auto flux_nr = fields.get_tmp_g();

    Micro_sedimentation_kernels::calc_flux_g<<<gridGPU, blockGPU>>>(
        flux_nr->fld_g, fields.sp.at("nr")->fld_g,
        slope_nr->fld_g, cfl_nr->fld_g,
        gd.dz_g, gd.dzi_g, fields.rhoref_g,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // 6. Limit fluxes such that nr stays >= 0.
    Micro_sedimentation_kernels::limit_flux_g<<<grid2dGPU, block2dGPU>>>(
        flux_nr->fld_g, fields.sp.at("nr")->fld_g,
        gd.dz_g, fields.rhoref_g, TF(dt),
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // 7. Calculate flux divergence
    Micro_sedimentation_kernels::calc_flux_div_g<<<gridGPU, blockGPU>>>(
        fields.st.at("nr")->fld_g, flux_nr->fld_g,
        gd.dzi_g, fields.rhoref_g,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    fields.release_tmp_g(cfl_nr);
    fields.release_tmp_g(slope_nr);
    fields.release_tmp_g(flux_nr);

    cudaDeviceSynchronize();

    stats.calc_tend(*fields.st.at("thl"), tend_name);
    stats.calc_tend(*fields.st.at("qt"),  tend_name);
    stats.calc_tend(*fields.st.at("qr"),  tend_name);
    stats.calc_tend(*fields.st.at("nr"),  tend_name);
}

template<typename TF>
unsigned long Microphys_2mom_warm<TF>::get_time_limit(unsigned long idt, const double dt)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    dim3 grid2dGPU (gridi, gridj);
    dim3 block2dGPU(blocki, blockj);

    // Calcuate sedimentation velocity qr, which is always higher than w_nr
    auto w_qr = fields.get_tmp_g();

    sedimentation_2mom_warm::calc_velocity_g<<<gridGPU, blockGPU>>>(
        w_qr->fld_g, fields.sp.at("qr")->fld_g, fields.sp.at("nr")->fld_g, fields.rhoref_g,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // Set the top and bottom velocity ghost cells
    Micro_sedimentation_kernels::set_bc_g<<<grid2dGPU, block2dGPU>>>(
        w_qr->fld_g,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // Calculate CFL number
    auto cfl_qr = fields.get_tmp_g();

    Micro_sedimentation_kernels::calc_cfl_g<<<gridGPU, blockGPU>>>(
        cfl_qr->fld_g, w_qr->fld_g, gd.dzi_g, TF(dt),
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    fields.release_tmp_g(w_qr);

    // Find maximum CFL across grid
    TF cfl = field3d_operators.calc_max_g(cfl_qr->fld_g);

    fields.release_tmp_g(cfl_qr);

    // Limit CFL at some small non-zero number
    cfl = fmax(cfl, cfl_min<TF>);

    const unsigned long idt_lim = idt * cflmax / cfl;
    return idt_lim;
}

template<typename TF>
void Microphys_2mom_warm<TF>::get_surface_rain_rate_g(
    TF* const __restrict__ rr_out)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj);
    dim3 blockGPU(blocki, blockj);

    sedimentation_2mom_warm::copy_rain_rate<<<gridGPU, blockGPU>>>(
        rr_out, rr_bot_g,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.icells);
    cuda_check_error();
}

template<typename TF>
void Microphys_2mom_warm<TF>::exec_column(Column<TF>& column)
{
    const TF no_offset = 0.;
    column.calc_time_series("rr", rr_bot_g, no_offset);
}

template<typename TF>
void Microphys_2mom_warm<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();
    const int memsize2d = gd.ijcells*sizeof(TF);

    cuda_safe_call(cudaMalloc(&rr_bot_g,  memsize2d));
}

template<typename TF>
void Microphys_2mom_warm<TF>::backward_device()
{
    auto& gd = grid.get_grid_data();
    const int dimemsize = gd.icells * sizeof(TF);

    cuda_safe_call(cudaMemcpy2D(rr_bot.data(), dimemsize, rr_bot_g, dimemsize, dimemsize, gd.jcells, cudaMemcpyDeviceToHost));
}

template<typename TF>
void Microphys_2mom_warm<TF>::clear_device()
{
    cuda_safe_call(cudaFree(rr_bot_g));
}
#endif


#ifdef FLOAT_SINGLE
template class Microphys_2mom_warm<float>;
#else
template class Microphys_2mom_warm<double>;
#endif
