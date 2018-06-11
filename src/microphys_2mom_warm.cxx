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

#include <iostream>
#include <cstdio>
#include <cmath>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "thermo_moist_functions.h"

#include "constants.h"
#include "microphys.h"
#include "microphys_2mom_warm.h"

using namespace Constants;
using namespace Thermo_moist_functions;

namespace
{
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
    TF* get_tmp_slice(std::vector<std::shared_ptr<Field3d<TF>>> &tmp_fields, int &slice_counter,
                      const int jcells, const int ikcells)
    {
        const int tmp_index   = slice_counter / jcells;     // Which tmp field in tmp_fields vector?
        const int fld_index   = slice_counter % jcells;     // Which slice in tmp field?
        const int slice_start = fld_index * ikcells;        // Start index of slice

        slice_counter++;

        return &(tmp_fields[tmp_index]->fld[slice_start]);
    }

    // Rational tanh approximation
    template<typename TF>
    inline TF tanh2(const TF x)
    {
        return x * (TF(27) + x * x) / (TF(27) + TF(9) * x * x);
    }

    // Given rain water content (qr), number density (nr) and density (rho)
    // calculate mean mass of rain drop
    template<typename TF>
    inline TF calc_rain_mass(const TF qr, const TF nr, const TF rho, const TF mr_min, const TF mr_max)
    {
        //TF mr = rho * qr / (nr + dsmall);
        TF mr = rho * qr / std::max(nr, TF(1.));
        return std::min(std::max(mr, mr_min), mr_max);
    }

    // Given mean mass rain drop, calculate mean diameter
    template<typename TF>
    inline TF calc_rain_diameter(const TF mr, const TF pirhow)
    {
        return pow(mr/pirhow, TF(1.)/TF(3.));
    }

    // Shape parameter mu_r
    template<typename TF>
    inline TF calc_mu_r(const TF dr)
    {
        //return 1./3.; // SB06
        //return 10. * (1. + tanh(1200 * (dr - 0.0015))); // SS08 (Milbrandt&Yau, 2005) -> similar as UCLA
        return 10. * (TF(1.) + tanh2(TF(1200.) * (dr - TF(0.0015)))); // SS08 (Milbrandt&Yau, 2005) -> similar as UCLA
    }

    // Slope parameter lambda_r
    template<typename TF>
    inline TF calc_lambda_r(const TF mur, const TF dr)
    {
        return pow((mur+3)*(mur+2)*(mur+1), TF(1.)/TF(3.)) / dr;
    }

    template<typename TF>
    inline TF minmod(const TF a, const TF b)
    {
        return copysign(TF(1.), a) * std::max(TF(0.), std::min(std::abs(a), TF(copysign(TF(1.), a))*b));
    }

    template<typename TF>
    inline TF min3(const TF a, const TF b, const TF c)
    {
        return std::min(a, std::min(b, c));
    }

    template<typename TF>
    inline TF max3(const TF a, const TF b, const TF c)
    {
        return std::max(a, std::max(b, c));
    }
}

// Microphysics calculated over entire 3D field
namespace mp3d
{
    // Autoconversion: formation of rain drop by coagulating cloud droplets
    template<typename TF>
    void autoconversion(TF* const restrict qrt, TF* const restrict nrt,
                        TF* const restrict qtt, TF* const restrict thlt,
                        const TF* const restrict qr,  const TF* const restrict ql,
                        const TF* const restrict rho, const TF* const restrict exner,
                        const int istart, const int jstart, const int kstart,
                        const int iend,   const int jend,   const int kend,
                        const int jj, const int kk,
                        Micro_2mom_warm_constants<TF> constants)
    {
        const TF x_star = 2.6e-10;       // SB06, list of symbols, same as UCLA-LES
        const TF k_cc   = 9.44e9;        // UCLA-LES (Long, 1974), 4.44e9 in SB06, p48
        const TF nu_c   = 1;             // SB06, Table 1., same as UCLA-LES
        const TF kccxs  = k_cc / (TF(20.) * x_star) * (nu_c+2)*(nu_c+4) / pow(nu_c+1, 2);

        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    if(ql[ijk] > constants.ql_min)
                    {
                        const TF xc      = rho[k] * ql[ijk] / constants.Nc0;    // Mean mass of cloud drops [kg]
                        const TF tau     = TF(1.) - ql[ijk] / (ql[ijk] + qr[ijk] + dsmall);    // SB06, Eq 5
                        const TF phi_au  = TF(600.) * pow(tau, TF(0.68)) * pow(TF(1.) - pow(tau, TF(0.68)), 3);    // UCLA-LES
                        //const TF phi_au  = 400. * pow(tau, 0.7) * pow(1. - pow(tau, 0.7), 3);    // SB06, Eq 6
                        const TF au_tend = rho[k] * kccxs * pow(ql[ijk], 2) * pow(xc, 2) *
                                               (TF(1.) + phi_au / pow(TF(1.)-tau, 2)); // SB06, eq 4

                        qrt[ijk]  += au_tend;
                        nrt[ijk]  += au_tend * rho[k] / x_star;
                        qtt[ijk]  -= au_tend;
                        thlt[ijk] += Lv / (cp * exner[k]) * au_tend;
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
                   const int jj, const int kk,
                   Micro_2mom_warm_constants<TF> constants)
    {
        const TF k_cr  = 5.25; // SB06, p49

        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    if(ql[ijk] > constants.ql_min && qr[ijk] > constants.qr_min)
                    {
                        const TF tau     = TF(1.) - ql[ijk] / (ql[ijk] + qr[ijk]); // SB06, Eq 5
                        const TF phi_ac  = pow(tau / (tau + TF(5e-5)), 4); // SB06, Eq 8
                        const TF ac_tend = k_cr * ql[ijk] *  qr[ijk] * phi_ac * pow(constants.rho_0 / rho[k], TF(0.5)); // SB06, Eq 7

                        qrt[ijk]  += ac_tend;
                        qtt[ijk]  -= ac_tend;
                        thlt[ijk] += Lv / (cp * exner[k]) * ac_tend;
                    }
                }
    }
}

// Microphysics calculated over 2D (xz) slices
namespace mp2d
{
    // Calculate microphysics properties which are used in multiple routines
    template<typename TF>
    void prepare_microphysics_slice(TF* const restrict rain_mass, TF* const restrict rain_diameter,
                                    TF* const restrict mu_r, TF* const restrict lambda_r,
                                    const TF* const restrict qr, const TF* const restrict nr,
                                    const TF* const restrict rho,
                                    const int istart, const int iend,
                                    const int kstart, const int kend,
                                    const int icells, const int ijcells, const int j,
                                    Micro_2mom_warm_constants<TF> constants)
    {

        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*icells + k*ijcells;
                const int ik  = i + k*icells;

                if(qr[ijk] > constants.qr_min)
                {
                    rain_mass[ik]     = calc_rain_mass(qr[ijk], nr[ijk], rho[k], constants.mr_min, constants.mr_max);
                    rain_diameter[ik] = calc_rain_diameter(rain_mass[ik], constants.pirhow);
                    mu_r[ik]          = calc_mu_r(rain_diameter[ik]);
                    lambda_r[ik]      = calc_lambda_r(mu_r[ik], rain_diameter[ik]);
                }
                else
                {
                    rain_mass[ik]     = 0;
                    rain_diameter[ik] = 0;
                    mu_r[ik]          = 0;
                    lambda_r[ik]      = 0;
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
                     const int jj, const int kk, const int j,
                     Micro_2mom_warm_constants<TF> constants)
    {
        const TF lambda_evap = 1.; // 1.0 in UCLA, 0.7 in DALES

        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ik  = i + k*jj;
                const int ijk = i + j*jj + k*kk;

                if(qr[ijk] > constants.qr_min)
                {
                    const TF mr  = rain_mass[ik];
                    const TF dr  = rain_diameter[ik];

                    const TF T   = thl[ijk] * exner[k] + (Lv * ql[ijk]) / (cp * exner[k]); // Absolute temperature [K]
                    const TF Glv = pow(Rv * T / (esat(T) * constants.D_v) +
                                       (Lv / (constants.K_t * T)) * (Lv / (Rv * T) - 1), -1); // Cond/evap rate (kg m-1 s-1)?

                    const TF S   = (qt[ijk] - ql[ijk]) / qsat(p[k], T) - 1; // Saturation
                    const TF F   = 1.; // Evaporation excludes ventilation term from SB06 (like UCLA, unimportant term? TODO: test)

                    const TF ev_tend = TF(2.) * constants.pi * dr * Glv * S * F * nr[ijk] / rho[k];

                    qrt[ijk]  += ev_tend;
                    nrt[ijk]  += lambda_evap * ev_tend * rho[k] / mr;
                    qtt[ijk]  -= ev_tend;
                    thlt[ijk] += Lv / (cp * exner[k]) * ev_tend;
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
                                const int jj, const int kk, const int j,
                                Micro_2mom_warm_constants<TF> constants)
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

                if(qr[ijk] > constants.qr_min)
                {
                    // Calculate mean rain drop mass and diameter
                    const TF dr      = rain_diameter[ik];
                    const TF lambdar = lambda_r[ik];

                    // Selfcollection
                    const TF sc_tend = -k_rr * nr[ijk] * qr[ijk]*rho[k] * pow(TF(1.) + kappa_rr /
                                       lambdar * pow(constants.pirhow, TF(1.)/TF(3.)), -9) * pow(constants.rho_0 / rho[k], TF(0.5));
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

    // Sedimentation from Stevens and Seifert (2008)
    template<typename TF>
    void sedimentation_ss08(TF* const restrict qrt, TF* const restrict nrt,
                            TF* const restrict w_qr, TF* const restrict w_nr,
                            TF* const restrict c_qr, TF* const restrict c_nr,
                            TF* const restrict slope_qr, TF* const restrict slope_nr,
                            TF* const restrict flux_qr, TF* const restrict flux_nr,
                            const TF* const restrict mu_r, const TF* const restrict lambda_r,
                            const TF* const restrict qr, const TF* const restrict nr,
                            const TF* const restrict rho, const TF* const restrict dzi,
                            const TF* const restrict dz, const double dt,
                            const int istart, const int jstart, const int kstart,
                            const int iend,   const int jend,   const int kend,
                            const int icells, const int kcells, const int ijcells, const int j,
                            Micro_2mom_warm_constants<TF> constants)
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
            const TF rho_n = pow(TF(1.2) / rho[k], TF(0.5));
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*icells + k*ijcells;
                const int ik  = i + k*icells;

                if(qr[ijk] > constants.qr_min)
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
            w_qr[ik2] = 0;
            w_nr[ik2] = 0;
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
            flux_qr[ik] = 0;
            flux_nr[ik] = 0;
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
                ftot  = 0;  // cumulative 'flux' (kg m-2)
                dzz   = 0;  // distance from zh[k]
                cc    = std::min(TF(1.), c_qr[ik]);
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
                ftot  = 0;  // cumulative 'flux'
                dzz   = 0;  // distance from zh[k]
                cc    = std::min(TF(1.), c_nr[ik]);
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
    }

}


template<typename TF>
Microphys_2mom_warm<TF>::Microphys_2mom_warm(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Microphys<TF>(masterin, gridin, fieldsin, inputin)
{
    swmicrophys = Microphys_type::Warm_2mom;

    #ifdef USECUDA
    throw std::runtime_error("swmicro = \"2mom_warm\" not (yet) implemented in CUDA\n");
    #endif

    // Read microphysics switches and settings
    swmicrobudget = inputin.get_item<bool>("micro", "swmicrobudget", "", false);
    cflmax        = inputin.get_item<TF>("micro", "cflmax", "", 2.);

    // Initialize the qr (rain water specific humidity) and nr (droplot number concentration) fields
    fields.init_prognostic_field("qr", "Rain water mixing ratio", "kg kg-1");
    fields.init_prognostic_field("nr", "Number density rain", "m-3");
}

template<typename TF>
Microphys_2mom_warm<TF>::~Microphys_2mom_warm()
{
}

template<typename TF>
void Microphys_2mom_warm<TF>::exec(Thermo<TF>& thermo, const double dt)
{
    auto& gd = grid.get_grid_data();

    // Remove spurious negative values from qr and nr fields
    remove_negative_values(fields.sp.at("qr")->fld.data(), gd.istart, gd.jstart, gd.kstart,
                           gd.iend, gd.jend, gd.kend, gd.icells, gd.ijcells);
    remove_negative_values(fields.sp.at("nr")->fld.data(), gd.istart, gd.jstart, gd.kstart,
                           gd.iend, gd.jend, gd.kend, gd.icells, gd.ijcells);

    // Get cloud liquid water specific humidity from thermodynamics
    auto ql = fields.get_tmp();
    thermo.get_thermo_field(*ql, "ql", false, false);

    // Get pressure and exner function from thermodynamics
    std::vector<TF> p     = thermo.get_p_vector();
    std::vector<TF> exner = thermo.get_exner_vector();

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
                         fields.sp.at("qr")->fld.data(), ql->fld.data(), fields.rhoref.data(), exner.data(),
                         gd.istart, gd.jstart, gd.kstart,
                         gd.iend,   gd.jend,   gd.kend,
                         gd.icells, gd.ijcells,
                         micro_constants);

    // Accretion; growth of raindrops collecting cloud droplets
    mp3d::accretion(fields.st.at("qr")->fld.data(), fields.st.at("qt")->fld.data(), fields.st.at("thl")->fld.data(),
                    fields.sp.at("qr")->fld.data(), ql->fld.data(), fields.rhoref.data(), exner.data(),
                    gd.istart, gd.jstart, gd.kstart,
                    gd.iend,   gd.jend,   gd.kend,
                    gd.icells, gd.ijcells,
                    micro_constants);

    // Rest of the microphysics is handled per XZ slice
    for (int j=gd.jstart; j<gd.jend; ++j)
    {
        // Prepare the XZ slices which are used in all routines
        mp2d::prepare_microphysics_slice(rain_mass, rain_diam, mu_r, lambda_r,
                                         fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(), fields.rhoref.data(),
                                         gd.istart, gd.iend, gd.kstart, gd.kend, gd.icells, gd.ijcells, j, micro_constants);

        // Evaporation; evaporation of rain drops in unsaturated environment
        mp2d::evaporation(fields.st.at("qr")->fld.data(), fields.st.at("nr")->fld.data(),  fields.st.at("qt")->fld.data(), fields.st.at("thl")->fld.data(),
                          fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(),  ql->fld.data(),
                          fields.sp.at("qt")->fld.data(), fields.sp.at("thl")->fld.data(), fields.rhoref.data(), exner.data(), p.data(),
                          rain_mass, rain_diam,
                          gd.istart, gd.jstart, gd.kstart,
                          gd.iend,   gd.jend,   gd.kend,
                          gd.icells, gd.ijcells, j,
                          micro_constants);

        // Self collection and breakup; growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
        mp2d::selfcollection_breakup(fields.st.at("nr")->fld.data(), fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(), fields.rhoref.data(),
                                     rain_mass, rain_diam, lambda_r,
                                     gd.istart, gd.jstart, gd.kstart,
                                     gd.iend,   gd.jend,   gd.kend,
                                     gd.icells, gd.ijcells, j,
                                     micro_constants);

        // Sedimentation; sub-grid sedimentation of rain
        mp2d::sedimentation_ss08(fields.st.at("qr")->fld.data(), fields.st.at("nr")->fld.data(),
                                 w_qr, w_nr, c_qr, c_nr, slope_qr, slope_nr, flux_qr, flux_nr, mu_r, lambda_r,
                                 fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(),
                                 fields.rhoref.data(), gd.dzi.data(), gd.dz.data(), dt,
                                 gd.istart, gd.jstart, gd.kstart,
                                 gd.iend,   gd.jend,   gd.kend,
                                 gd.icells, gd.kcells, gd.ijcells, j,
                                 micro_constants);

        // Sedimentation; sub-grid sedimentation of rain
        mp2d::sedimentation_ss08(fields.st.at("qr")->fld.data(), fields.st.at("nr")->fld.data(),
                                 w_qr, w_nr, c_qr, c_nr, slope_qr, slope_nr, flux_qr, flux_nr, mu_r, lambda_r,
                                 fields.sp.at("qr")->fld.data(), fields.sp.at("nr")->fld.data(),
                                 fields.rhoref.data(), gd.dzi.data(), gd.dz.data(), dt,
                                 gd.istart, gd.jstart, gd.kstart,
                                 gd.iend,   gd.jend,   gd.kend,
                                 gd.icells, gd.kcells, gd.ijcells, j,
                                 micro_constants);
    }

    // Release all local tmp fields in use
    for (auto& it: tmp_fields)
        fields.release_tmp(it);

    fields.release_tmp(ql);
}

template<typename TF>
void Microphys_2mom_warm<TF>::exec_stats(Stats<TF>& stats, std::string mask_name,
                                         Field3d<TF>& mask_field, Field3d<TF>& mask_fieldh, const double dt)
{
}

template<typename TF>
void Microphys_2mom_warm<TF>::exec_dump(Dump<TF>& dump, unsigned long iotime)
{
}

template<typename TF>
void Microphys_2mom_warm<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
}

template<typename TF>
unsigned long Microphys_2mom_warm<TF>::get_time_limit(unsigned long idt, const double dt)
{
    return Constants::ulhuge;
}

template<typename TF>
bool Microphys_2mom_warm<TF>::has_mask(std::string name)
{
    return false;
}

template<typename TF>
void Microphys_2mom_warm<TF>::get_mask(Field3d<TF>& mfield, Field3d<TF>& mfieldh, Stats<TF>& stats, std::string mask_name)
{
}

template class Microphys_2mom_warm<double>;
template class Microphys_2mom_warm<float>;
