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
#include "column.h"
#include "thermo.h"
#include "thermo_moist_functions.h"
#include "constants.h"

#include "microphys.h"
#include "microphys_sb06.h"

namespace
{
    using namespace Constants;
    namespace fm = Fast_math;
    namespace tmf = Thermo_moist_functions;

    // For simplicity, define constants here for now. These should probably move to the header.
    template<typename TF> constexpr TF pi       = 3.14159265359;
    template<typename TF> constexpr TF K_t      = 2.5e-2;              // Conductivity of heat [J/(sKm)]
    template<typename TF> constexpr TF D_v      = 3.e-5;               // Diffusivity of water vapor [m2/s]
    template<typename TF> constexpr TF rho_w    = 1.e3;                // Density water
    template<typename TF> constexpr TF rho_0    = 1.225;               // SB06, p48
    template<typename TF> constexpr TF pirhow   = pi<TF>*rho_w<TF>/6.;
    template<typename TF> constexpr TF mc_min   = 4.2e-15;             // Min mean mass of cloud droplet
    template<typename TF> constexpr TF mc_max   = 2.6e-10;             // Max mean mass of cloud droplet
    template<typename TF> constexpr TF mr_min   = mc_max<TF>;          // Min mean mass of precipitation drop
    template<typename TF> constexpr TF mr_max   = 3e-6;                // Max mean mass of precipitation drop // as in UCLA-LES
    template<typename TF> constexpr TF ql_min   = 1.e-6;               // Min cloud liquid water for which calculations are performed
    template<typename TF> constexpr TF qr_min   = 1.e-15;              // Min rain liquid water for which calculations are performed
    template<typename TF> constexpr TF cfl_min  = 1.e-5;               // Small non-zero limit at the CFL number
    template<typename TF> constexpr TF rho_vel  = 0.4;                 // Exponent for density correction (value from ICON)
    template<typename TF> constexpr TF q_crit   = 1.e-9;               // Min rain liquid water for which calculations are performed (value from ICON)

    template<typename TF>
    inline TF tanh2(const TF x)
    {
        // Rational `tanh` approximation
        return x * (TF(27) + x * x) / (TF(27) + TF(9) * x * x);
    }

    template<typename TF>
    inline TF calc_rain_mass(const TF qr, const TF nr)
    {
        // Calculate mean mass of rain droplets (kg)
        TF mr = qr / std::max(nr, TF(1.));
        return std::min(std::max(mr, mr_min<TF>), mr_max<TF>);
    }

    template<typename TF>
    inline TF calc_rain_diameter(const TF mr)
    {
        // Given mean mass rain drop, calculate mean diameter (m)
        return pow(mr/pirhow<TF>, TF(1.)/TF(3.));
    }

    template<typename TF>
    inline TF calc_mu_r(const TF dr)
    {
        // Calculate shape parameter mu_r
        //return 1./3.; // SB06
        //return 10. * (1. + tanh(1200 * (dr - 0.0015))); // SS08 (Milbrandt&Yau, 2005) -> similar as UCLA
        return TF(10) * (TF(1) + tanh2(TF(1200) * (dr - TF(0.0015)))); // SS08 (Milbrandt&Yau, 2005) -> similar as UCLA
    }

    template<typename TF> CUDA_MACRO
    inline TF calc_lambda_r(const TF mur, const TF dr)
    {
        // Calculate Slope parameter lambda_r
        return pow((mur+3)*(mur+2)*(mur+1), TF(1.)/TF(3.)) / dr;
    }


    template<typename TF>
    void prepare_microphysics_slice(
            TF* const restrict rain_mass,
            TF* const restrict rain_diameter,
            TF* const restrict mu_r,
            TF* const restrict lambda_r,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jstride + k*kstride;
                const int ij  = i + j*jstride;

                if (qr[ijk] > qr_min<TF>)
                {
                    rain_mass[ij]     = calc_rain_mass(qr[ijk], nr[ijk]);
                    rain_diameter[ij] = calc_rain_diameter(rain_mass[ij]);
                    mu_r[ij]          = calc_mu_r(rain_diameter[ij]);
                    lambda_r[ij]      = calc_lambda_r(mu_r[ij], rain_diameter[ij]);
                }
                else
                {
                    rain_mass[ij]     = TF(0);
                    rain_diameter[ij] = TF(0);
                    mu_r[ij]          = TF(0);
                    lambda_r[ij]      = TF(0);
                }
            }
    }


    template<typename TF>
    void convert_units(
            TF* const restrict qt,
            TF* const restrict qr,
            TF* const restrict qs,
            TF* const restrict qg,
            TF* const restrict qh,
            TF* const restrict qi,
            TF* const restrict ql,
            const TF* const restrict rho,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride,
            bool to_kgm3)
    {
        TF fac;

        for (int k=kstart; k<kend; k++)
        {
            if (to_kgm3)
                fac = rho[k];
            else
                fac = TF(1)/rho[k];

            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jstride + k*kstride;

                    qt[ijk] *= fac;
                    ql[ijk] *= fac;
                    qi[ijk] *= fac;
                    qr[ijk] *= fac;
                    qs[ijk] *= fac;
                    qg[ijk] *= fac;
                    qh[ijk] *= fac;
                }
        }
    }


    template<typename TF>
    void autoconversion(
            TF* const restrict qrt,
            TF* const restrict nrt,
            TF* const restrict qtt,
            TF* const restrict thlt,
            const TF* const restrict qr,
            const TF* const restrict ql,
            const TF* const restrict rho,
            const TF* const restrict exner,
            const TF Nc0,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        /* Formation of rain by coagulating cloud droplets */

        const TF x_star = 2.6e-10;       // SB06, list of symbols, same as UCLA-LES (kg)
        const TF k_cc = 9.44e9;          // UCLA-LES (Long, 1974), 4.44e9 in SB06, p48, same as ICON (m3 kg-2 s-1)
        const TF nu_c = 1;               // SB06, Table 1., same as UCLA-LES (-)
        const TF kccxs = k_cc / (TF(20.) * x_star) * (nu_c+2)*(nu_c+4) / fm::pow2(nu_c+1);

        for (int k=kstart; k<kend; k++)
        {
            const TF rho_i = TF(1)/rho[k];

            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jstride + k*kstride;

                    if (ql[ijk] > ql_min<TF>)
                    {
                        // Mean mass of cloud drops (kg):
                        const TF xc = ql[ijk] / Nc0;
                        // Dimensionless internal time scale (SB06, Eq 5):
                        const TF tau = TF(1.) - ql[ijk] / (ql[ijk] + qr[ijk] + dsmall);
                        // Collection kernel (SB06, Eq 6). Constants are from UCLA-LES and differ from SB06:
                        const TF phi_au = TF(600.) * pow(tau, TF(0.68)) * fm::pow3(TF(1.) - pow(tau, TF(0.68)));
                        // Autoconversion tendency (SB06, Eq 4, kg m-3 s-1):
                        const TF au_tend = rho_0<TF>/rho[k] * kccxs * fm::pow2(ql[ijk]) * fm::pow2(xc) *
                                               (TF(1.) + phi_au / fm::pow2(TF(1.)-tau)); // SB06, eq 4

                        qrt[ijk]  += au_tend;
                        nrt[ijk]  += au_tend / x_star;
                        qtt[ijk]  -= au_tend;
                        thlt[ijk] += rho_i * Lv<TF> / (cp<TF> * exner[k]) * au_tend;
                    }
                }
        }
    }


    template<typename TF>
    void accretion(
            TF* const restrict qrt,
            TF* const restrict qtt,
            TF* const restrict thlt,
            const TF* const restrict qr,
            const TF* const restrict ql,
            const TF* const restrict rho,
            const TF* const restrict exner,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        /* Accreation: growth of raindrops collecting cloud droplets */

        const TF k_cr = 5.25; // SB06, p49 (m3 kg-1 s-1)

        for (int k=kstart; k<kend; k++)
        {
            const TF rho_i = TF(1)/rho[k];

            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jstride + k*kstride;

                    if (ql[ijk] > ql_min<TF> && qr[ijk] > qr_min<TF>)
                    {
                        // Dimensionless internal time scale (SB06, Eq 5):
                        const TF tau = TF(1.) - ql[ijk] / (ql[ijk] + qr[ijk]);
                        // Collection kernel (SB06, Eq 8):
                        const TF phi_ac = fm::pow4(tau / (tau + TF(5e-5)));
                        // Accreation tendency (SB06, Eq 7, kg m-3 s-1):
                        const TF ac_tend = k_cr * ql[ijk] *  qr[ijk] * phi_ac * sqrt(rho_0<TF> / rho[k]);

                        qrt[ijk]  += ac_tend;
                        qtt[ijk]  -= ac_tend;
                        thlt[ijk] += rho_i * Lv<TF> / (cp<TF> * exner[k]) * ac_tend;
                    }
                }
        }
    }


    template<typename TF>
    void evaporation(
            TF* const restrict qrt,
            TF* const restrict nrt,
            TF* const restrict qtt,
            TF* const restrict thlt,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF* const restrict ql,
            const TF* const restrict qi,
            const TF* const restrict qt,
            const TF* const restrict T,
            const TF* const restrict rho,
            const TF* const restrict exner,
            const TF* const restrict p,
            const TF* const restrict rain_mass,
            const TF* const restrict rain_diameter,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        // Evaporation: evaporation of rain drops in unsaturated environment
        const TF lambda_evap = TF(1.); // 1.0 in UCLA, 0.7 in DALES
        const TF rho_i = TF(1.) / rho[k];

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij  = i + j*jstride;
                const int ijk = i + j*jstride + k*kstride;

                if (qr[ijk] > qr_min<TF>)
                {
                    const TF mr  = rain_mass[ij];
                    const TF dr  = rain_diameter[ij];

                    // ...Condensation/evaporation rate...?
                    const TF Glv = TF(1.) / (Rv<TF> * T[ijk] / (tmf::esat_liq(T[ijk]) * D_v<TF>) +
                                       (Lv<TF> / (K_t<TF> * T[ijk])) * (Lv<TF> / (Rv<TF> * T[ijk]) - TF(1.)));

                    // Supersaturation over water (-).
                    // NOTE: copy-pasted from UCLA, can't find the definition in SB06.
                    const TF qv = qt[ijk] - ql[ijk] - qi[ijk];
                    const TF S   = qv / tmf::qsat_liq(p[k], T[ijk]) - TF(1.);

                    // Ventilation factor. UCLA-LES=1, calculated in SB06 = TODO..
                    const TF F   = TF(1.);

                    // Evaporation tendency (kg m-3 s-1).
                    const TF ev_tend = TF(2.) * pi<TF> * dr * Glv * S * F * nr[ijk];

                    qrt[ijk]  += ev_tend;
                    nrt[ijk]  += lambda_evap * ev_tend / mr;
                    qtt[ijk]  -= ev_tend;
                    thlt[ijk] += rho_i * Lv<TF> / (cp<TF> * exner[k]) * ev_tend;
                }
            }
    }


    template<typename TF>
    void selfcollection_breakup(
            TF* const restrict nrt,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF* const restrict rho,
            const TF* const restrict rain_mass,
            const TF* const restrict rain_diameter,
            const TF* const restrict lambda_r,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        // Selfcollection & breakup: growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
        const TF k_rr     = 7.12;   // SB06, p49
        const TF kappa_rr = 60.7;   // SB06, p49
        const TF D_eq     = 0.9e-3; // SB06, list of symbols
        const TF k_br1    = 1.0e3;  // SB06, p50, for 0.35e-3 <= Dr <= D_eq
        const TF k_br2    = 2.3e3;  // SB06, p50, for Dr > D_eq

        for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ij  = i + j*jstride;
                    const int ijk = i + j*jstride + k*kstride;

                    if (qr[ijk] > qr_min<TF>)
                    {
                        // Short-cuts...
                        const TF dr = rain_diameter[ij];

                        // Selfcollection tendency:
                        // NOTE: this is quite different in ICON, UCLA-LES had 4 different versions, ...
                        const TF sc_tend = -k_rr * nr[ijk] * qr[ijk] /
                            pow(TF(1.) + kappa_rr /
                            lambda_r[ij] * pow(pirhow<TF>, TF(1.)/TF(3.)), -9) * sqrt(rho_0<TF> / rho[k]);

                        nrt[ijk] += sc_tend;

                        // Breakup by collisions:
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

    template<typename TF>
    inline TF particle_meanmass(
            Particle<TF>& particle,
            const TF q, const TF n)
    {
        // Mean mass of particle, with limiters (SB06, Eq 94)
        return std::min( std::max( q/(n+TF(Constants::dsmall)), particle.x_min ), particle.x_max );
    }

    template<typename TF>
    inline TF particle_diameter(
            Particle<TF>& particle,
            const TF mean_mass)
    {
        // Mass-diameter relation (SB06, Eq 32)
        return particle.a_geo * std::exp(particle.b_geo * std::log(mean_mass));
    }

    template<typename TF>
    inline TF rain_mue_dm_relation(
            Particle_rain_coeffs<TF>& coeffs,
            const TF d_m)
    {
        // mue-Dm relation of raindrops.
        TF mue;
        const TF delta = coeffs.cmu2 * (d_m - coeffs.cmu3);

        if (d_m <= coeffs.cmu3)
           mue = coeffs.cmu0 * std::tanh(fm::pow2(TF(4) * delta)) + coeffs.cmu4;
        else
           mue = coeffs.cmu1 * std::tanh(fm::pow2(delta)) + coeffs.cmu4;

        return mue;
    }

    template<typename TF>
    void sedi_vel_rain(
            TF* const restrict vn,
            TF* const restrict vq,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF* const restrict ql,
            Particle<TF>& rain,
            Particle_rain_coeffs<TF>& coeffs,
            const TF rho_corr,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k,
            bool qc_present)
    {
        TF mue;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j * jstride;
                const int ijk = i + j * jstride + k * kstride;

                if (qr[ijk] > q_crit<TF>)
                {
                    const TF x = particle_meanmass(rain, qr[ijk], nr[ijk]);
                    const TF d_m = particle_diameter(rain, x);

                    if (qc_present)
                    {
                        if (ql[ijk] >= q_crit<TF>)
                            mue = (rain.nu + TF(1.0)) / rain.b_geo - TF(1.0);
                        else
                            mue = rain_mue_dm_relation(coeffs, d_m);
                    } else
                        mue = rain_mue_dm_relation(coeffs, d_m);

                    const TF d_p =
                            d_m * std::exp(TF(-1. / 3.) * std::log((mue + TF(3)) * (mue + TF(2)) * (mue + TF(1))));
                    vn[ij] = coeffs.alfa - coeffs.beta * std::exp(-(mue + TF(1)) * std::log(TF(1) + coeffs.gama * d_p));
                    vq[ij] = coeffs.alfa - coeffs.beta * std::exp(-(mue + TF(4)) * std::log(TF(1) + coeffs.gama * d_p));

                    vn[ij] *= rho_corr;
                    vq[ij] *= rho_corr;
                } else
                {
                    vn[ij] = TF(0);
                    vq[ij] = TF(0);
                }
            }
    }


    template<typename TF>
    void copy_slice(
            TF* const restrict fld_2d,
            const TF* const restrict fld_3d,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        for (int j = jstart; j < jend; j++)
                #pragma ivdep
                for (int i = istart; i < iend; i++)
                {
                    const int ij = i + j * jstride;
                    const int ijk = i + j * jstride + k * kstride;

                    fld_2d[ij] = fld_3d[ijk];
                }
    }


    template<typename TF>
    void implicit_core(
            TF* const restrict q_val,
            TF* const restrict q_sum,
            TF* const restrict q_impl,
            TF* const restrict vsed_new,
            TF* const restrict vsed_now,
            TF* const restrict flux_new,
            TF* const restrict flux_now,
            const TF rdzdt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        for (int j = jstart; j < jend; j++)
            #pragma ivdep
            for (int i = istart; i < iend; i++)
            {
                const int ij = i + j * jstride;

                // `new` on r.h.s. is new value from level above
                vsed_new[ij] = TF(0.5) * (vsed_now[ij] + vsed_new[ij]);

                // `flux_new` are the updated flux values from the level above
                // `flux_now` are here the old (current time step) flux values from the level above
                const TF flux_sum = flux_new[ij] + flux_now[ij];

                // `flux_now` are here overwritten with the current level
                flux_now[ij] = std::min(vsed_now[ij] * q_val[ij], flux_sum);   // loop dependency
                flux_now[ij] = std::max(flux_now[ij], TF(0));                  // Maybe not necessary

                // Time integrated value without implicit weight
                q_sum[ij] = q_val[ij] + rdzdt * (flux_sum - flux_now[ij]);

                // Implicit weight
                q_impl[ij] = TF(1) / (TF(1) + vsed_new[ij] * rdzdt);

                // prepare for source term calculation
                const TF q_star = q_impl[ij] * q_sum[ij];
                q_val[ij]  = q_star;       // source/sinks work on star-values
                q_sum[ij]  = q_sum[ij] - q_star;
            }
    }


    template<typename TF>
    void implicit_time(
            TF* const restrict q_val,
            TF* const restrict q_sum,
            TF* const restrict q_impl,
            TF* const restrict vsed_new,
            TF* const restrict vsed_now,
            TF* const restrict flux_new,
            const TF rdzdt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        for (int j = jstart; j < jend; j++)
                #pragma ivdep
                for (int i = istart; i < iend; i++)
                {
                    const int ij = i + j * jstride;

                    // Time integration
                    q_val[ij] = std::max(TF(0), q_impl[ij] * (q_sum[ij] + q_val[ij]));

                    // Prepare for next level
                    flux_new[ij] = q_val[ij] * vsed_new[ij];
                    vsed_new[ij] = vsed_now[ij];
                }
    }


}

template<typename TF>
Microphys_sb06<TF>::Microphys_sb06(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
        Microphys<TF>(masterin, gridin, fieldsin, inputin)
{
    auto& gd = grid.get_grid_data();
    swmicrophys = Microphys_type::SB06;

    // Read microphysics switches and settings
    cflmax = inputin.get_item<TF>("micro", "cflmax", "", 1.2);
    Nc0 = inputin.get_item<TF>("micro", "Nc0", "");

    // Initialize the qr (rain water specific humidity) and nr (droplot number concentration) fields
    const std::string group_name = "thermo";
    fields.init_prognostic_field("qi", "Ice specific humidity", "kg kg-1", group_name, gd.sloc);
    fields.init_prognostic_field("qr", "Rain specific humidity", "kg kg-1", group_name, gd.sloc);
    fields.init_prognostic_field("qs", "Snow specific humidity", "kg kg-1", group_name, gd.sloc);
    fields.init_prognostic_field("qg", "Graupel specific humidity", "kg kg-1", group_name, gd.sloc);
    fields.init_prognostic_field("qh", "Hail specific humidity", "kg kg-1", group_name, gd.sloc);

    fields.init_prognostic_field("ni", "Number density ice", "m-3", group_name, gd.sloc);
    fields.init_prognostic_field("nr", "Number density rain", "m-3", group_name, gd.sloc);
    fields.init_prognostic_field("ns", "Number density snow", "m-3", group_name, gd.sloc);
    fields.init_prognostic_field("ng", "Number density graupel", "m-3", group_name, gd.sloc);
    fields.init_prognostic_field("nh", "Number density hail", "m-3", group_name, gd.sloc);

    // Load the viscosity for both fields.
    fields.sp.at("qi")->visc = inputin.get_item<TF>("fields", "svisc", "qi");
    fields.sp.at("qr")->visc = inputin.get_item<TF>("fields", "svisc", "qr");
    fields.sp.at("qg")->visc = inputin.get_item<TF>("fields", "svisc", "qg");
    fields.sp.at("qs")->visc = inputin.get_item<TF>("fields", "svisc", "qs");
    fields.sp.at("qh")->visc = inputin.get_item<TF>("fields", "svisc", "qh");

    fields.sp.at("ni")->visc = inputin.get_item<TF>("fields", "svisc", "ni");
    fields.sp.at("nr")->visc = inputin.get_item<TF>("fields", "svisc", "nr");
    fields.sp.at("ng")->visc = inputin.get_item<TF>("fields", "svisc", "ng");
    fields.sp.at("ns")->visc = inputin.get_item<TF>("fields", "svisc", "ns");
    fields.sp.at("nh")->visc = inputin.get_item<TF>("fields", "svisc", "nh");

    // Set the rain and rain_coeff types to the default provided values
    rain = rainSBB;
    rain_coeffs = rainSBBcoeffs;
}

template<typename TF>
Microphys_sb06<TF>::~Microphys_sb06()
{
}

template<typename TF>
void Microphys_sb06<TF>::init()
{
    auto& gd = grid.get_grid_data();

    rr_bot.resize(gd.ijcells);
    rs_bot.resize(gd.ijcells);
    rg_bot.resize(gd.ijcells);
}

template<typename TF>
void Microphys_sb06<TF>::create(
        Input& inputin, Netcdf_handle& input_nc,
        Stats<TF>& stats, Cross<TF>& cross, Dump<TF>& dump, Column<TF>& column)
{
    const std::string group_name = "micro";

    // Add variables to the statistics
    if (stats.get_switch())
    {
        // Time series
        stats.add_time_series("rr", "Mean surface rain rate", "kg m-2 s-1", group_name);
        stats.add_time_series("rs", "Mean surface snow rate", "kg m-2 s-1", group_name);
        stats.add_time_series("rg", "Mean surface graupel rate", "kg m-2 s-1", group_name);

        stats.add_tendency(*fields.st.at("thl"), "z", tend_name, tend_longname);
        stats.add_tendency(*fields.st.at("qt") , "z", tend_name, tend_longname);
        stats.add_tendency(*fields.st.at("qr") , "z", tend_name, tend_longname);
        stats.add_tendency(*fields.st.at("qs") , "z", tend_name, tend_longname);
        stats.add_tendency(*fields.st.at("qg") , "z", tend_name, tend_longname);

        // Profiles
        stats.add_prof("vqr", "Fall velocity rain mass density", "m s-1", "z" , group_name);
        stats.add_prof("vnr", "Fall velocity rain number density", "m s-1", "z" , group_name);
    }

    if (column.get_switch())
    {
        column.add_time_series("rr", "Surface rain rate", "kg m-2 s-1");
        column.add_time_series("rs", "Surface snow rate", "kg m-2 s-1");
        column.add_time_series("rg", "Surface graupel rate", "kg m-2 s-1");
    }

    // Create cross sections
    // 1. Variables that this class can calculate/provide:
    const std::vector<std::string> allowed_crossvars = {"rr_bot", "rs_bot", "rg_bot"};

    // 2. Cross-reference with the variables requested in the .ini file:
    crosslist = cross.get_enabled_variables(allowed_crossvars);
}

#ifndef USECUDA
template<typename TF>
void Microphys_sb06<TF>::exec(Thermo<TF>& thermo, const double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    // Get thermodynamic variables
    bool cyclic = false;
    bool is_stat = false;

    auto ql = fields.get_tmp();
    thermo.get_thermo_field(*ql, "ql", cyclic, is_stat);

    auto T = fields.get_tmp();
    thermo.get_thermo_field(*T, "T", cyclic, is_stat);

    const std::vector<TF>& p = thermo.get_basestate_vector("p");
    const std::vector<TF>& exner = thermo.get_basestate_vector("exner");
    const std::vector<TF>& rho = thermo.get_basestate_vector("rho");

    // 2D slices for quantities shared between different kernels.
    auto rain_mass = fields.get_tmp_xy();
    auto rain_diameter = fields.get_tmp_xy();
    auto mu_r = fields.get_tmp_xy();
    auto lambda_r = fields.get_tmp_xy();

    // 2D slices for sedimentation.
    auto get_zero_tmp_xy = [&]()
    {
        auto fld_xy = fields.get_tmp_xy();
        std::fill((*fld_xy).begin(), (*fld_xy).end(), TF(0));
        return fld_xy;
    };

    auto vr_sedq_now = get_zero_tmp_xy();
    auto vr_sedq_new = get_zero_tmp_xy();
    auto qr_flux_now = get_zero_tmp_xy();
    auto qr_flux_new = get_zero_tmp_xy();
    auto qr_sum = get_zero_tmp_xy();
    auto qr_impl = get_zero_tmp_xy();
    auto qr_slice = get_zero_tmp_xy();

    auto vr_sedn_now = get_zero_tmp_xy();
    auto vr_sedn_new = get_zero_tmp_xy();
    auto nr_flux_now = get_zero_tmp_xy();
    auto nr_flux_new = get_zero_tmp_xy();
    auto nr_sum = get_zero_tmp_xy();
    auto nr_impl = get_zero_tmp_xy();
    auto nr_slice = get_zero_tmp_xy();

    // Convert all units from `kg kg-1 to `kg m-3`.
    bool to_kgm3 = true;
    convert_units(
        fields.ap.at("qt")->fld.data(),
        fields.ap.at("qr")->fld.data(),
        fields.ap.at("qs")->fld.data(),
        fields.ap.at("qg")->fld.data(),
        fields.ap.at("qh")->fld.data(),
        fields.ap.at("qi")->fld.data(),
        ql->fld.data(),
        rho.data(),
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells,
        to_kgm3);
   
    for (int k=gd.kend-1; k>=gd.kstart; --k)
    {
        // Copy 3D fields to 2D slices.
        copy_slice(
                (*qr_slice).data(),
                fields.sp.at("qr")->fld.data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells, gd.ijcells, k);

        copy_slice(
                (*nr_slice).data(),
                fields.sp.at("nr")->fld.data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells, gd.ijcells, k);

        // Sedimentation
        // Density correction fall speeds
        const TF hlp = std::log(std::max(rho[k], TF(1e-6)) / rho_0<TF>);
        const TF rho_corr = std::exp(-rho_vel<TF>*hlp);

        bool ql_present = true;
        sedi_vel_rain(
            (*vr_sedn_now).data(),
            (*vr_sedq_now).data(),
            fields.sp.at("qr")->fld.data(),
            fields.sp.at("nr")->fld.data(),
            ql->fld.data(),
            rain, rain_coeffs,
            rho_corr,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells, gd.ijcells,
            k, ql_present);

        const TF rdzdt = TF(0.5) * gd.dzi[k] * dt;

        implicit_core(
            (*qr_slice).data(),
            (*qr_sum).data(),
            (*qr_impl).data(),
            (*vr_sedq_new).data(),
            (*vr_sedq_now).data(),
            (*qr_flux_new).data(),
            (*qr_flux_now).data(),
            rdzdt,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

        implicit_core(
            (*nr_slice).data(),
            (*nr_sum).data(),
            (*nr_impl).data(),
            (*vr_sedn_new).data(),
            (*vr_sedn_now).data(),
            (*nr_flux_new).data(),
            (*nr_flux_now).data(),
            rdzdt,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

        // Calculate microphysics processes.
        // Autoconversion; formation of rain drop by coagulating cloud droplets.
        autoconversion(
            fields.st.at("qr")->fld.data(),
            fields.st.at("nr")->fld.data(),
            fields.st.at("qt")->fld.data(),
            fields.st.at("thl")->fld.data(),
            fields.sp.at("qr")->fld.data(),
            ql->fld.data(),
            rho.data(),
            exner.data(), Nc0,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            k, k+1,
            gd.icells, gd.ijcells);

        // Accretion; growth of rain droplets by collecting cloud droplets
        accretion(
            fields.st.at("qr")->fld.data(),
            fields.st.at("qt")->fld.data(),
            fields.st.at("thl")->fld.data(),
            fields.sp.at("qr")->fld.data(),
            ql->fld.data(),
            rho.data(),
            exner.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            k, k+1,
            gd.icells, gd.ijcells);

        // Calculate quantities used by multiple kernels.
        prepare_microphysics_slice(
            (*rain_mass).data(),
            (*rain_diameter).data(),
            (*mu_r).data(),
            (*lambda_r).data(),
            fields.sp.at("qr")->fld.data(),
            fields.sp.at("nr")->fld.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells, gd.ijcells, k);

        // Evaporation of rain droplets.
        evaporation(
            fields.st.at("qr")->fld.data(),
            fields.st.at("nr")->fld.data(),
            fields.st.at("qt")->fld.data(),
            fields.st.at("thl")->fld.data(),
            fields.sp.at("qr")->fld.data(),
            fields.sp.at("nr")->fld.data(),
            ql->fld.data(),
            fields.sp.at("qi")->fld.data(),
            fields.sp.at("qt")->fld.data(),
            T->fld.data(),
            rho.data(),
            exner.data(),
            p.data(),
            (*rain_mass).data(),
            (*rain_diameter).data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells, gd.ijcells, k);

        selfcollection_breakup(
            fields.st.at("nr")->fld.data(),
            fields.sp.at("qr")->fld.data(),
            fields.sp.at("nr")->fld.data(),
            rho.data(),
            (*rain_mass).data(),
            (*rain_diameter).data(),
            (*lambda_r).data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells, gd.ijcells, k);

        // Sedimentation
        implicit_time(
                (*qr_slice).data(),
                (*qr_sum).data(),
                (*qr_impl).data(),
                (*vr_sedq_new).data(),
                (*vr_sedq_now).data(),
                (*qr_flux_new).data(),
                rdzdt,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells);

        implicit_time(
                (*nr_slice).data(),
                (*nr_sum).data(),
                (*nr_impl).data(),
                (*vr_sedn_new).data(),
                (*vr_sedn_now).data(),
                (*nr_flux_new).data(),
                rdzdt,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells);

        // TO-DO: calculate tendency from difference between e.g. `qr_slice` and `qr`.
    }

    // Calculate tendencies.
    stats.calc_tend(*fields.st.at("thl"), tend_name);
    stats.calc_tend(*fields.st.at("qt" ), tend_name);
    stats.calc_tend(*fields.st.at("qi" ), tend_name);
    stats.calc_tend(*fields.st.at("qr" ), tend_name);
    stats.calc_tend(*fields.st.at("qs" ), tend_name);
    stats.calc_tend(*fields.st.at("qg" ), tend_name);

    // Convert all units from `kg m-3 to `kg kg-1`.
    to_kgm3 = false;
    convert_units(
        fields.ap.at("qt")->fld.data(),
        fields.ap.at("qr")->fld.data(),
        fields.ap.at("qs")->fld.data(),
        fields.ap.at("qg")->fld.data(),
        fields.ap.at("qh")->fld.data(),
        fields.ap.at("qi")->fld.data(),
        ql->fld.data(),
        rho.data(),
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells,
        to_kgm3);

    // Release temporary fields.
    fields.release_tmp(ql);
    fields.release_tmp(T);

    fields.release_tmp_xy(rain_mass);
    fields.release_tmp_xy(rain_diameter);
    fields.release_tmp_xy(mu_r);
    fields.release_tmp_xy(lambda_r);

    fields.release_tmp_xy(vr_sedq_now);
    fields.release_tmp_xy(vr_sedq_new);
    fields.release_tmp_xy(qr_flux_now);
    fields.release_tmp_xy(qr_flux_new);
    fields.release_tmp_xy(qr_sum);
    fields.release_tmp_xy(qr_impl);
    fields.release_tmp_xy(qr_slice);

    fields.release_tmp_xy(vr_sedn_now);
    fields.release_tmp_xy(vr_sedn_new);
    fields.release_tmp_xy(nr_flux_now);
    fields.release_tmp_xy(nr_flux_new);
    fields.release_tmp_xy(nr_sum);
    fields.release_tmp_xy(nr_impl);
    fields.release_tmp_xy(nr_slice);
}
#endif

template<typename TF>
void Microphys_sb06<TF>::exec_stats(Stats<TF>& stats, Thermo<TF>& thermo, const double dt)
{
    auto& gd = grid.get_grid_data();

    const TF no_offset = 0.;
    const TF no_threshold = 0.;
    const bool is_stat = true;
    const bool cyclic = false;

    // Time series
    stats.calc_stats_2d("rr", rr_bot, no_offset);
    stats.calc_stats_2d("rs", rs_bot, no_offset);
    stats.calc_stats_2d("rg", rg_bot, no_offset);

    // Profiles
    auto vq = fields.get_tmp();
    auto vn = fields.get_tmp();

    const std::vector<TF>& rho = thermo.get_basestate_vector("rho");
    auto ql = fields.get_tmp();
    thermo.get_thermo_field(*ql, "ql", cyclic, is_stat);

    for (int k=gd.kend-1; k>=gd.kstart; --k)
    {
        // Sedimentation
        // Density correction fall speeds
        const TF hlp = std::log(std::max(rho[k], TF(1e-6)) / rho_0<TF>);
        const TF rho_corr = std::exp(-rho_vel<TF> * hlp);

        bool ql_present = true;
        sedi_vel_rain(
                &vn->fld.data()[k * gd.ijcells],
                &vq->fld.data()[k * gd.ijcells],
                fields.sp.at("qr")->fld.data(),
                fields.sp.at("nr")->fld.data(),
                ql->fld.data(),
                rain, rain_coeffs,
                rho_corr,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells, gd.ijcells,
                k, ql_present);
    }

    stats.calc_stats("vqr", *vq, no_offset, no_threshold);
    stats.calc_stats("vnr", *vn, no_offset, no_threshold);

    fields.release_tmp(vq);
    fields.release_tmp(vn);
    fields.release_tmp(ql);
}

#ifndef USECUDA
template<typename TF>
void Microphys_sb06<TF>::exec_column(Column<TF>& column)
{
    const TF no_offset = 0.;
    column.calc_time_series("rr", rr_bot.data(), no_offset);
    column.calc_time_series("rs", rs_bot.data(), no_offset);
    column.calc_time_series("rg", rg_bot.data(), no_offset);
}
#endif

template<typename TF>
void Microphys_sb06<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    if (cross.get_switch())
    {
        for (auto& it : crosslist)
        {
            if (it == "rr_bot")
                cross.cross_plane(rr_bot.data(), "rr_bot", iotime);
            if (it == "rs_bot")
                cross.cross_plane(rs_bot.data(), "rs_bot", iotime);
            if (it == "rg_bot")
                cross.cross_plane(rg_bot.data(), "rg_bot", iotime);
        }
    }
}

#ifndef USECUDA
template<typename TF>
unsigned long Microphys_sb06<TF>::get_time_limit(unsigned long idt, const double dt)
{
    auto& gd = grid.get_grid_data();
    auto tmp = fields.get_tmp();

    double cfl = 0.;

    // TO-DO

    // Get maximum CFL across all MPI tasks
    master.max(&cfl, 1);
    fields.release_tmp(tmp);

    // Prevent zero division.
    cfl = std::max(cfl, 1.e-5);

    return idt * this->cflmax / cfl;
}
#endif

template<typename TF>
bool Microphys_sb06<TF>::has_mask(std::string name)
{
    return false;
}

template<typename TF>
void Microphys_sb06<TF>::get_mask(Stats<TF>& stats, std::string mask_name)
{
    std::string message = "SB06 microphysics scheme can not provide mask: \"" + mask_name +"\"";
    throw std::runtime_error(message);
}

template<typename TF>
void Microphys_sb06<TF>::get_surface_rain_rate(std::vector<TF>& field)
{
    // Make a hard copy of the surface rain precipitation field
    field = rr_bot;

    // Add snow and graupel surface precipitation
    std::transform(field.begin(), field.end(), rs_bot.begin(), field.begin(), std::plus<TF>());
    std::transform(field.begin(), field.end(), rg_bot.begin(), field.begin(), std::plus<TF>());
}


template class Microphys_sb06<double>;
template class Microphys_sb06<float>;
