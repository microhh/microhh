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
#include <stdlib.h>

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
#include "timeloop.h"

#include "microphys.h"

#include "microphys_sb06.h"         // Class definition
#include "microphys_sb_common.h"    // Warm/cold common kernels
#include "microphys_sb_cold.h"      // Kernels from ICON
#include "microphys_sb_init.h"      // Initialisation code.

namespace
{
    using namespace Constants;
    namespace fm = Fast_math;
    namespace tmf = Thermo_moist_functions;
}

namespace warm
{
    /*
       These are the "old" kernels used in the old `2mom_warm` scheme,
       ported (mainly) from UCLA-LES and DALES.
    */

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
            const TF* const restrict rho,
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

                if (qr[ij] > qr_min<TF> * rho[k])  // TODO: remove *rho...
                {
                    rain_mass[ij]     = calc_rain_mass(qr[ij], nr[ij]);
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
    void autoconversion(
            TF* const restrict qrt,
            TF* const restrict nrt,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF* const restrict ql,
            const TF* const restrict rho,
            const TF* const restrict exner,
            const TF Nc0, const double dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        /* Formation of rain by coagulating cloud droplets */

        const TF x_star = 2.6e-10;       // SB06, list of symbols, same as UCLA-LES (kg)
        const TF k_cc = 9.44e9;          // UCLA-LES (Long, 1974), 4.44e9 in SB06, p48, same as ICON (m3 kg-2 s-1)
        const TF nu_c = 1;               // SB06, Table 1., same as UCLA-LES (-)
        const TF kccxs = k_cc / (TF(20.) * x_star) * (nu_c+2)*(nu_c+4) / fm::pow2(nu_c+1);

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;
                const int ijk = i + j*jstride + k*kstride;

                if (ql[ijk] > ql_min<TF> * rho[k])  // TODO: remove *rho[k]...
                {
                    // Mean mass of cloud drops (kg):
                    const TF xc = ql[ijk] / Nc0;
                    // Dimensionless internal time scale (SB06, Eq 5):
                    const TF tau = TF(1.) - ql[ijk] / (ql[ijk] + qr[ij] + dsmall);
                    // Collection kernel (SB06, Eq 6). Constants are from UCLA-LES and differ from SB06:
                    const TF phi_au = TF(600.) * pow(tau, TF(0.68)) * fm::pow3(TF(1.) - pow(tau, TF(0.68)));
                    // Autoconversion tendency (SB06, Eq 4, kg m-3 s-1):
                    const TF au_tend = rho_0<TF>/rho[k] * kccxs * fm::pow2(ql[ijk]) * fm::pow2(xc) *
                                       (TF(1.) + phi_au / fm::pow2(TF(1.)-tau)); // SB06, eq 4

                    // Update 2D slices:
                    qrt[ij] += au_tend;
                    nrt[ij] += au_tend / x_star;
                }
            }
    }


    template<typename TF>
    void accretion(
            TF* const restrict qrt,
            const TF* const restrict qr,
            const TF* const restrict ql,
            const TF* const restrict rho,
            const TF* const restrict exner,
            const double dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        /* Accretion: growth of raindrops collecting cloud droplets */

        const TF k_cr = 5.25; // SB06, p49 (m3 kg-1 s-1)

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;
                const int ijk = i + j*jstride + k*kstride;

                if (ql[ijk] > ql_min<TF> * rho[k] && qr[ij] > qr_min<TF> * rho[k])  // TODO: remove *rho[k]...
                {
                    // Dimensionless internal time scale (SB06, Eq 5):
                    const TF tau = TF(1.) - ql[ijk] / (ql[ijk] + qr[ij]);
                    // Collection kernel (SB06, Eq 8):
                    const TF phi_ac = fm::pow4(tau / (tau + TF(5e-5)));
                    // Accreation tendency (SB06, Eq 7, kg m-3 s-1):
                    const TF ac_tend = k_cr * ql[ijk] *  qr[ij] * phi_ac * sqrt(rho_0<TF> / rho[k]);

                    // Update 2D slice:
                    qrt[ij] += ac_tend;
                }
            }
    }


    template<typename TF>
    void evaporation(
            TF* const restrict qrt,
            TF* const restrict nrt,
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
            const double dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        // Evaporation: evaporation of rain drops in unsaturated environment
        const TF lambda_evap = TF(1.); // 1.0 in UCLA, 0.7 in DALES

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij  = i + j*jstride;
                const int ijk = i + j*jstride + k*kstride;

                if (qr[ij] > qr_min<TF> * rho[k]) // TODO: remove *rho..
                {
                    // Supersaturation over water (-, KP97 Eq 4-37).
                    // TMP BvS, for comparison with old scheme:
                    //const TF qv = qt[ijk] - ql[ijk] - qi[ijk];
                    const TF qv = qt[ijk] - ql[ijk];

                    const TF qs = tmf::qsat_liq(p[k], T[ijk]) * rho[k];
                    const TF S  = qv / qs - TF(1.);

                    if (S < TF(0))
                    {
                        const TF mr  = rain_mass[ij];
                        const TF dr  = rain_diameter[ij];

                        // ...Condensation/evaporation rate...?
                        const TF es = tmf::esat_liq(T[ijk]);
                        const TF Glv = TF(1.) / (Rv<TF> * T[ijk] / (es * D_v<TF>) +
                                                 (Lv<TF> / (K_t<TF> * T[ijk])) * (Lv<TF> / (Rv<TF> * T[ijk]) - TF(1.)));

                        // Ventilation factor. UCLA-LES=1, calculated in SB06 = TODO..
                        const TF F   = TF(1.);

                        // Evaporation tendency (kg m-3 s-1).
                        const TF ev_tend = TF(2.) * pi<TF> * dr * Glv * S * F * nr[ij];

                        // Update 2D slices:
                        qrt[ij] += ev_tend;
                        nrt[ij] += lambda_evap * ev_tend / mr;
                    }
                }
            }
    }


    template<typename TF>
    void selfcollection_breakup(
            TF* const restrict nrt,
            const TF* const restrict nr,
            const TF* const restrict qr,
            const TF* const restrict rho,
            const TF* const restrict rain_mass,
            const TF* const restrict rain_diameter,
            const TF* const restrict lambda_r,
            const double dt,
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

                if (qr[ij] > qr_min<TF> * rho[k]) // TODO: remove *rho..
                {
                    // Short-cuts...
                    const TF dr = rain_diameter[ij];

                    // Selfcollection tendency:
                    // NOTE: this is quite different in ICON, UCLA-LES had 4 different versions, ...
                    const TF sc_tend = -k_rr * nr[ij] * qr[ij] /
                                       pow(TF(1.) + kappa_rr /
                                                    lambda_r[ij] * pow(pirhow<TF>, TF(1.)/TF(3.)), -9) * sqrt(rho_0<TF> / rho[k]);

                    // Update 2D slice:
                    nrt[ij] += sc_tend;

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

                        // Update 2D slice:
                        nrt[ij] += br_tend;
                    }
                }
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
    sw_warm = inputin.get_item<bool>("micro", "swwarm", "", false);
    sw_microbudget = inputin.get_item<bool>("micro", "swmicrobudget", "", false);
    sw_debug = inputin.get_item<bool>("micro", "swdebug", "", false);
    sw_integrate = inputin.get_item<bool>("micro", "swintegrate", "", false);

    cfl_max = inputin.get_item<TF>("micro", "cflmax", "", 1.2);
    Nc0 = inputin.get_item<TF>("micro", "Nc0", "");
    Ni0 = inputin.get_item<TF>("micro", "Ni0", "");

    auto add_type = [&](
            const std::string& symbol,
            const std::string& short_name,
            const std::string& long_name,
            const std::string& units,
            const bool is_mass)
    {
        Hydro_type<TF> tmp = {short_name, long_name, units, is_mass};
        hydro_types.emplace(symbol, tmp);
    };

    // Switch between mass and density for precipitation types.
    const bool is_mass = true;

    if (sw_warm)
    {
        add_type("qr", "rain", "rain specific humidity", "kg kg-1", is_mass);
        add_type("nr", "rain", "number density rain", "kg-1", !is_mass);
    }
    else
    {
        add_type("qi", "ice", "ice specific humidity", "kg kg-1", is_mass);
        add_type("ni", "ice", "number density ice", "kg-1", !is_mass);

        add_type("qr", "rain", "rain specific humidity", "kg kg-1", is_mass);
        add_type("nr", "rain", "number density rain", "kg-1", !is_mass);

        add_type("qs", "snow", "snow specific humidity", "kg kg-1", is_mass);
        add_type("ns", "snow", "number density snow", "kg-1", !is_mass);

        add_type("qg", "graupel", "graupel specific humidity", "kg kg-1", is_mass);
        add_type("ng", "graupel", "number density graupel", "kg-1", !is_mass);

        add_type("qh", "hail", "hail specific humidity", "kg kg-1", is_mass);
        add_type("nh", "hail", "number density hail", "kg-1", !is_mass);
    }

    const std::string group_name = "thermo";
    for (auto& it : hydro_types)
    {
        fields.init_prognostic_field(it.first, it.second.long_name, it.second.units, group_name, gd.sloc);
        fields.sp.at(it.first)->visc = inputin.get_item<TF>("fields", "svisc", it.first);
    }

    // Setup/calculate cloud/rain/particle/... coefficients.
    init_2mom_scheme_once();

    // Init lookup table for conversion graupel to hail (wet growth):
    // NOTE: ICON supports input through ASCII and NetCDF files, and also
    //       some hard-coded internal lookup table.. ICON by default uses
    //       the NetCDF table; that's the only option that is supported here.
    const std::string file_name = "dmin_wetgrowth_lookup_61.nc";

    Sb_init::init_dmin_wg_gr_ltab_equi(
            ltabdminwgg,
            file_name,
            graupel,
            master);
}

template<typename TF>
Microphys_sb06<TF>::~Microphys_sb06()
{
}

template<typename TF>
void Microphys_sb06<TF>::init_2mom_scheme()
{
    cloud = cloud_nue1mue1;
    rain = rainSBB;
    ice = ice_cosmo5;
    snow = snowSBB;

    graupel = graupelhail_cosmo5;

    //SELECT TYPE (graupel)
    //TYPE IS (particle_frozen)
    //  call particle_frozen_assign(graupel,graupelhail_cosmo5)
    //TYPE IS (particle_lwf)
    //  call particle_lwf_assign(graupel,graupel_vivek)
    //END SELECT

    hail = hail_cosmo5;

    //SELECT TYPE (hail)
    //TYPE IS (particle_frozen)
    //  call particle_frozen_assign(hail,hail_cosmo5)
    //TYPE IS (particle_lwf)
    //  call particle_lwf_assign(hail,hail_vivek)
    //END SELECT


}

template<typename TF>
void Microphys_sb06<TF>::init_2mom_scheme_once()
{
    if (sw_debug)
    {
        master.print_message("Init SB06 2 moment scheme\n");
        master.print_message("---------------------------------------\n");
    }

    // bulk ventilation coefficient, Eq. (88) of SB2006
    auto vent_coeff_a = [](
            const Particle<TF> &parti, const int n)
    {
        const TF vent_coeff_a = parti.a_ven
                                * std::tgamma((parti.nu + n + parti.b_geo) / parti.mu)
                                / std::tgamma((parti.nu + 1.0) / parti.mu)
                                * std::pow(std::tgamma((parti.nu + 1.0) / parti.mu)
                                           / std::tgamma((parti.nu + 2.0) / parti.mu), (parti.b_geo + n - 1.0));

        return vent_coeff_a;
    };

    // bulk ventilation coefficient, Eq. (89) of SB2006
    auto vent_coeff_b = [](
            const Particle<TF> &parti, const int n)
    {
        constexpr TF m_f = 0.500; // see PK, S.541. Do not change.

        const TF vent_coeff_b = parti.b_ven
                                * std::tgamma((parti.nu + n + (m_f + 1.0) * parti.b_geo + m_f * parti.b_vel) / parti.mu)
                                / std::tgamma((parti.nu + 1.0) / parti.mu)
                                * std::pow(std::tgamma((parti.nu + 1.0) / parti.mu)
                                           / std::tgamma((parti.nu + 2.0) / parti.mu),
                                           ((m_f + 1.0) * parti.b_geo + m_f * parti.b_vel + n - 1.0));

        return vent_coeff_b;
    };

    auto moment_gamma = [](
            const Particle<TF> &p, const int n)
    {
        const TF moment_gamma = std::tgamma((n + p.nu + 1.0) / p.mu) / std::tgamma((p.nu + 1.0) / p.mu)
                                * std::pow(std::tgamma((p.nu + 1.0) / p.mu) / std::tgamma((p.nu + 2.0) / p.mu), n);
        return moment_gamma;
    };

    // Setup the subset of the coefficients that does not follow from the copy.
    auto setup_particle_coeffs = [&](
            const Particle<TF> &ptype,
            Particle_coeffs<TF> &pcoeffs)
    {
        constexpr TF N_sc = 0.710;  // Schmidt-Zahl (PK, S.541)
        constexpr TF n_f = 0.333;   // Exponent von N_sc im Vent-koeff. (PK, S.541)
        constexpr TF nu_l = 1.5e-5; // Kinematic viscosity of air (added by CvH).

        pcoeffs.c_i = 1. / ptype.cap;
        pcoeffs.a_f = vent_coeff_a(ptype, 1);
        pcoeffs.b_f = vent_coeff_b(ptype, 1) * std::pow(N_sc, n_f) / std::sqrt(nu_l);
        pcoeffs.c_z = moment_gamma(ptype, 2);

        if (sw_debug)
        {
            char l1 = ptype.name[0];
            master.print_message("Setup_particle_coeffs: %s\n",ptype.name.c_str());
            master.print_message(" | a_geo = %f\n", ptype.a_geo);
            master.print_message(" | b_geo = %f\n", ptype.b_geo);
            master.print_message(" | a_vel = %f\n", ptype.a_vel);
            master.print_message(" | b_vel = %f\n", ptype.b_vel);
            master.print_message(" | c_%c = %f\n", l1, pcoeffs.c_i);
            master.print_message(" | a_f = %f\n", pcoeffs.a_f);
            master.print_message(" | b_f = %f\n", pcoeffs.b_f);
        }
    };

    // initialize coefficients for bulk sedimentation velocity.
    auto init_2mom_sedi_vel = [&](
            const Particle_frozen<TF>& particle,
            Particle_sphere<TF>& coeffs)
    {
        coeffs.coeff_alfa_n = particle.a_vel * std::tgamma((particle.nu + particle.b_vel + TF(1)) / particle.mu) /
                std::tgamma((particle.nu + TF(1)) / particle.mu);
        coeffs.coeff_alfa_q = particle.a_vel * std::tgamma((particle.nu + particle.b_vel + TF(2)) / particle.mu) /
                std::tgamma((particle.nu + TF(2)) / particle.mu);
        coeffs.coeff_lambda = std::tgamma((particle.nu + TF(1)) / particle.mu) /
                std::tgamma((particle.nu + TF(2)) / particle.mu);

        if (sw_debug)
        {
            master.print_message(" | name = %s\n", particle.name.c_str());
            master.print_message(" | c_lam = %f\n", coeffs.coeff_lambda);
            master.print_message(" | alf_n = %f\n", coeffs.coeff_alfa_n);
            master.print_message(" | alf_q = %f\n", coeffs.coeff_alfa_q);
        }
    };

    auto coll_delta = [](
            const Particle<TF>& p1, const int n)
    {
        return std::tgamma((TF(2) * p1.b_geo + p1.nu + TF(1) + n) / p1.mu)
             / std::tgamma((p1.nu + TF(1)) / p1.mu)
             * pow(std::tgamma((p1.nu + TF(1)) / p1.mu), (TF(2) * p1.b_geo + n))
             / pow(std::tgamma((p1.nu + TF(2)) / p1.mu), (TF(2) * p1.b_geo + n));
    };

    // wrapper for coll_delta (unnecessary and unused argument p2, but do not remove this)
    auto coll_delta_11 = [&](
            const Particle<TF>& p1,
            const Particle<TF>& p2,
            const int n)
    {
        return coll_delta(p1, n);
    };

    // wrapper for coll_delta (unnecessary and unused argument p1, but do not remove this)
    auto coll_delta_22 = [&](
            const Particle<TF>& p1,
            const Particle<TF>& p2,
            const int n)
    {
        return coll_delta(p2, n);
    };

    // coefficient for general collision integral, Eq. (91) of SB2006
    auto coll_delta_12 = [&](
            const Particle<TF>& p1,
            const Particle<TF>& p2,
            const int n)
    {
        return TF(2) * std::tgamma((p1.b_geo + p1.nu + TF(1)) / p1.mu)
                     / std::tgamma((p1.nu + TF(1)) / p1.mu)
                     * pow(std::tgamma((p1.nu + TF(1)) / p1.mu), (p1.b_geo))
                     / pow(std::tgamma((p1.nu + TF(2)) / p1.mu), (p1.b_geo))
                     * std::tgamma((p2.b_geo+p2.nu+TF(1) + n)/p2.mu)
                     / std::tgamma((p2.nu + TF(1)) / p2.mu)
                     * pow(std::tgamma((p2.nu + TF(1)) / p2.mu), (p2.b_geo + n))
                     / pow(std::tgamma((p2.nu + TF(2)) / p2.mu), (p2.b_geo + n));
    };

    // coefficient for general collision integral, Eq. (92) of SB2006
    auto coll_theta = [](
            const Particle<TF>& p1, const int n)
    {
        return std::tgamma((TF(2) * p1.b_vel + TF(2) * p1.b_geo + p1.nu + TF(1) + n) / p1.mu)
             / std::tgamma((TF(2) * p1.b_geo + p1.nu + TF(1) + n) / p1.mu)
             * pow(std::tgamma((p1.nu + TF(1)) / p1.mu), (TF(2) * p1.b_vel))
             / pow(std::tgamma((p1.nu + TF(2)) / p1.mu), (TF(2) * p1.b_vel));
    };

    // wrapper for coll_theta (unnecessary and unused argument p2, but do not remove this)
    auto coll_theta_11 = [&](
            const Particle<TF>& p1,
            const Particle<TF>& p2,
            const int n)
    {
        return coll_theta(p1, n);
    };

    // wrapper for coll_theta (unnecessary and unused argument p1, but do not remove this)
    auto coll_theta_22 = [&](
            const Particle<TF>& p1,
            const Particle<TF>& p2,
            const int n)
    {
        return coll_theta(p2, n);
    };

    // coefficient for general collision integral, Eq. (93) of SB2006
    auto coll_theta_12 = [&](
            const Particle<TF>& p1,
            const Particle<TF>& p2,
            const int n)
    {
        return TF(2) * std::tgamma((p1.b_vel + TF(2) * p1.b_geo + p1.nu + TF(1)) / p1.mu)
                     / std::tgamma((TF(2) * p1.b_geo + p1.nu + TF(1)) / p1.mu)
                     * pow(std::tgamma((p1.nu + TF(1)) / p1.mu), p1.b_vel)
                     / pow(std::tgamma((p1.nu + TF(2)) / p1.mu), p1.b_vel)
                     * std::tgamma((p2.b_vel + TF(2) * p2.b_geo + p2.nu + TF(1) + n) / p2.mu)
                     / std::tgamma((TF(2) * p2.b_geo + p2.nu + TF(1) + n) / p2.mu)
                     * pow(std::tgamma((p2.nu + TF(1)) / p2.mu), p2.b_vel)
                     / pow(std::tgamma((p2.nu + TF(2)) / p2.mu), p2.b_vel);
    };

    auto setup_particle_collection_type1 = [&](
            const Particle<TF>& ptype,
            const Particle<TF>& qtype,
            Collection_coeffs<TF>& coll_coeffs,
            const std::string pname,
            const std::string qname,
            const bool print_types = true)
    {
        coll_coeffs.delta_n_aa = coll_delta_11(ptype, qtype, 0);
        coll_coeffs.delta_n_ab = coll_delta_12(ptype, qtype, 0);
        coll_coeffs.delta_n_bb = coll_delta_22(ptype, qtype, 0);
        coll_coeffs.delta_q_aa = coll_delta_11(ptype, qtype, 0);
        coll_coeffs.delta_q_ab = coll_delta_12(ptype, qtype, 1);
        coll_coeffs.delta_q_bb = coll_delta_22(ptype, qtype, 1);

        coll_coeffs.theta_n_aa = coll_theta_11(ptype, qtype, 0);
        coll_coeffs.theta_n_ab = coll_theta_12(ptype, qtype, 0);
        coll_coeffs.theta_n_bb = coll_theta_22(ptype, qtype, 0);
        coll_coeffs.theta_q_aa = coll_theta_11(ptype, qtype, 0);
        coll_coeffs.theta_q_ab = coll_theta_12(ptype, qtype, 1);
        coll_coeffs.theta_q_bb = coll_theta_22(ptype, qtype, 1);

        if (sw_debug)
        {
            char l1 = pname[0];
            char l2 = qname[0];

            master.print_message("%s-%s riming:\n", pname.c_str(), qname.c_str());

            if (print_types)
            {
                master.print_message(" | a_%s = %f\n", pname.c_str(),ptype.a_geo);
                master.print_message(" | b_%s = %f\n", pname.c_str(),ptype.b_geo);
                master.print_message(" | alf_%s = %f\n", pname.c_str(),ptype.a_vel);
                master.print_message(" | bet_%s = %f\n", pname.c_str(),ptype.b_vel);
                master.print_message(" | a_%s = %f\n", qname.c_str(),qtype.a_geo);
                master.print_message(" | b_%s = %f\n", qname.c_str(),qtype.b_geo);
                master.print_message(" | alf_%s  = %f\n", qname.c_str(),qtype.a_vel);
                master.print_message(" | bet_%s  = %f\n", qname.c_str(),qtype.b_vel);
            }

            master.print_message(" | delta_n_%c%c = %f\n", l1, l1, coll_coeffs.delta_n_aa);
            master.print_message(" | delta_n_%c%c = %f\n", l1, l2, coll_coeffs.delta_n_ab);
            master.print_message(" | delta_n_%c%c = %f\n", l2, l2, coll_coeffs.delta_n_bb);
            master.print_message(" | theta_n_%c%c = %f\n", l1, l1, coll_coeffs.theta_n_aa);
            master.print_message(" | theta_n_%c%c = %f\n", l1, l2, coll_coeffs.theta_n_ab);
            master.print_message(" | theta_n_%c%c = %f\n", l2, l2, coll_coeffs.theta_n_bb);
            master.print_message(" | delta_q_%c%c = %f\n", l1, l1, coll_coeffs.delta_q_aa);
            master.print_message(" | delta_q_%c%c = %f\n", l1, l2, coll_coeffs.delta_q_ab);
            master.print_message(" | delta_q_%c%c = %f\n", l2, l2, coll_coeffs.delta_q_bb);
            master.print_message(" | theta_q_%c%c = %f\n", l1, l1, coll_coeffs.theta_q_aa);
            master.print_message(" | theta_q_%c%c = %f\n", l1, l2, coll_coeffs.theta_q_ab);
            master.print_message(" | theta_q_%c%c = %f\n", l2, l2, coll_coeffs.theta_q_bb);
        }
    };

    auto setup_particle_collection_type2 = [&](
            const Particle<TF>& ptype,
            const Particle<TF>& qtype,
            Rain_riming_coeffs<TF>& coll_coeffs,
            const std::string pname,
            const std::string qname)
    {
        coll_coeffs.delta_n_aa = coll_delta_11(ptype, qtype, 0);
        coll_coeffs.delta_n_ab = coll_delta_12(ptype, qtype, 0);
        coll_coeffs.delta_n_bb = coll_delta_22(ptype, qtype, 0);
        coll_coeffs.delta_q_aa = coll_delta_11(ptype, qtype, 1); // mass weighted
        coll_coeffs.delta_q_ab = coll_delta_12(ptype, qtype, 1);
        coll_coeffs.delta_q_ba = coll_delta_12(qtype, ptype, 1);
        coll_coeffs.delta_q_bb = coll_delta_22(ptype, qtype, 1);

        coll_coeffs.theta_n_aa = coll_theta_11(ptype, qtype, 0);
        coll_coeffs.theta_n_ab = coll_theta_12(ptype, qtype, 0);
        coll_coeffs.theta_n_bb = coll_theta_22(ptype, qtype, 0);
        coll_coeffs.theta_q_aa = coll_theta_11(ptype, qtype, 1); // mass weighted
        coll_coeffs.theta_q_ab = coll_theta_12(ptype, qtype, 1);
        coll_coeffs.theta_q_ba = coll_theta_12(qtype, ptype, 1);
        coll_coeffs.theta_q_bb = coll_theta_22(ptype, qtype, 1);

        if (sw_debug)
        {
            char l1 = pname[0];
            char l2 = qname[0];

            master.print_message("%s-%s riming (TEST-abc!):\n", pname.c_str(), qname.c_str());
            master.print_message(" | a_%s = %f\n", pname.c_str(), ptype.a_geo);
            master.print_message(" | b_%s = %f\n", pname.c_str(), ptype.b_geo);
            master.print_message(" | alf_%s = %f\n", pname.c_str(), ptype.a_vel);
            master.print_message(" | bet_%s = %f\n", pname.c_str(), ptype.b_vel);
            master.print_message(" | a_%s = %f\n", qname.c_str(), qtype.a_geo);
            master.print_message(" | b_%s = %f\n", qname.c_str(), qtype.b_geo);
            master.print_message(" | alf_%s = %f\n", qname.c_str(), qtype.a_vel);
            master.print_message(" | bet_%s = %f\n", qname.c_str(), qtype.b_vel);
            master.print_message(" | delta_n_%c%c = %f\n", l1, l1, coll_coeffs.delta_n_aa);
            master.print_message(" | delta_n_%c%c = %f\n", l1, l2, coll_coeffs.delta_n_ab);
            master.print_message(" | delta_n_%c%c = %f\n", l2, l2, coll_coeffs.delta_n_bb);
            master.print_message(" | theta_n_%c%c = %f\n", l1, l1, coll_coeffs.theta_n_aa);
            master.print_message(" | theta_n_%c%c = %f\n", l1, l2, coll_coeffs.theta_n_ab);
            master.print_message(" | theta_n_%c%c = %f\n", l2, l2, coll_coeffs.theta_n_bb);
            master.print_message(" | delta_q_%c%c = %f\n", l1, l1, coll_coeffs.delta_q_aa);
            master.print_message(" | delta_q_%c%c = %f\n", l1, l2, coll_coeffs.delta_q_ab);
            master.print_message(" | delta_q_%c%c = %f\n", l2, l1, coll_coeffs.delta_q_ba);
            master.print_message(" | delta_q_%c%c = %f\n", l2, l2, coll_coeffs.delta_q_bb);
            master.print_message(" | theta_q_%c%c = %f\n", l1, l1, coll_coeffs.theta_q_aa);
            master.print_message(" | theta_q_%c%c = %f\n", l1, l2, coll_coeffs.theta_q_ab);
            master.print_message(" | theta_q_%c%c = %f\n", l2, l1, coll_coeffs.theta_q_ba);
            master.print_message(" | theta_q_%c%c = %f\n", l2, l2, coll_coeffs.theta_q_bb);
        }
    };

    auto setup_graupel_selfcollection = [&](
            const Particle<TF>& graupel,
            Particle_graupel_coeffs<TF>& graupel_coeffs)
    {
        const TF delta_n_11 = coll_delta_11(graupel, graupel, 0);
        const TF delta_n_12 = coll_delta_12(graupel, graupel, 0);
        const TF theta_n_11 = coll_theta_11(graupel, graupel, 0);
        const TF theta_n_12 = coll_theta_12(graupel, graupel, 0);

        const TF delta_n = (TF(2)*delta_n_11 + delta_n_12);
        const TF theta_n = std::pow((TF(2)*theta_n_11 - theta_n_12), TF(0.5));

        graupel_coeffs.sc_coll_n = Sb_cold::pi8<TF> * delta_n * theta_n;

        if (sw_debug)
        {
            master.print_message("setup_graupel_selfcollection:\n");
            master.print_message(" | delta_n_11 = %f\n", delta_n_11);
            master.print_message(" | delta_n_12 = %f\n", delta_n_12);
            master.print_message(" | delta_n = %f\n", delta_n);
            master.print_message(" | theta_n_11 = %f\n", theta_n_11);
            master.print_message(" | theta_n_12 = %f\n", theta_n_12);
            master.print_message(" | theta_n = %f\n", theta_n);
            master.print_message(" | coll_n = %f\n", graupel_coeffs.sc_coll_n);
        }
    };

    auto setup_snow_selfcollection = [&](
            const Particle<TF>& snow,
            Particle_snow_coeffs<TF>& snow_coeffs)
    {

        const TF delta_n_11 = coll_delta_11(snow, snow, 0);
        const TF delta_n_12 = coll_delta_12(snow, snow, 0);
        const TF theta_n_11 = coll_theta_11(snow, snow, 0);
        const TF theta_n_12 = coll_theta_12(snow, snow, 0);

        snow_coeffs.sc_delta_n = (TF(2)*delta_n_11 + delta_n_12);
        snow_coeffs.sc_theta_n = (TF(2)*theta_n_11 - theta_n_12);

        if (sw_debug)
        {
            master.print_message("setup_snow_selfcollection:\n");
            // Only printed in ICON in full debug mode:
            //master.print_message(" | a_snow = %f\n", snow.a_geo);
            //master.print_message(" | b_snow = %f\n", snow.b_geo);
            //master.print_message(" | alf_snow = %f\n", snow.a_vel);
            //master.print_message(" | bet_snow = %f\n", snow.b_vel);
            //master.print_message(" | delta_n_11 = %f\n", delta_n_11);
            //master.print_message(" | delta_n_12 = %f\n", delta_n_12);
            //master.print_message(" | theta_n_11 = %f\n", theta_n_11);
            //master.print_message(" | theta_n_12 = %f\n", theta_n_12);
            master.print_message(" | delta_n = %f\n", snow_coeffs.sc_delta_n);
            master.print_message(" | theta_n = %f\n", snow_coeffs.sc_theta_n);
        }
    };

    auto setup_ice_selfcollection = [&](
            const Particle<TF>& ice,
            Particle_ice_coeffs<TF>& ice_coeffs)
    {
        const TF delta_n_11 = coll_delta_11(ice,ice,0);
        const TF delta_n_12 = coll_delta_12(ice,ice,0);
        const TF delta_n_22 = coll_delta_22(ice,ice,0);
        const TF delta_q_11 = coll_delta_11(ice,ice,0);
        const TF delta_q_12 = coll_delta_12(ice,ice,1);
        const TF delta_q_22 = coll_delta_22(ice,ice,1);

        const TF theta_n_11 = coll_theta_11(ice,ice,0);
        const TF theta_n_12 = coll_theta_12(ice,ice,0);
        const TF theta_n_22 = coll_theta_22(ice,ice,0);
        const TF theta_q_11 = coll_theta_11(ice,ice,0);
        const TF theta_q_12 = coll_theta_12(ice,ice,1);
        const TF theta_q_22 = coll_theta_22(ice,ice,1);

        ice_coeffs.sc_delta_n = delta_n_11 + delta_n_12 + delta_n_22;
        ice_coeffs.sc_delta_q = delta_q_11 + delta_q_12 + delta_q_22;
        ice_coeffs.sc_theta_n = theta_n_11 - theta_n_12 + theta_n_22;
        ice_coeffs.sc_theta_q = theta_q_11 - theta_q_12 + theta_q_22;

        if(sw_debug)
        {
            master.print_message("setup_ice_selfcollection:\n");
            // Only printed in ICON in full debug mode:
            //master.print_message(" | a_ice      = %f\n",ice.a_geo);
            //master.print_message(" | b_ice      = %f\n",ice.b_geo);
            //master.print_message(" | alf_ice    = %f\n",ice.a_vel);
            //master.print_message(" | bet_ice    = %f\n",ice.b_vel);
            //master.print_message(" | delta_n_11 = %f\n",delta_n_11);
            //master.print_message(" | delta_n_12 = %f\n",delta_n_12);
            //master.print_message(" | delta_n_22 = %f\n",delta_n_22);
            //master.print_message(" | theta_n_11 = %f\n",theta_n_11);
            //master.print_message(" | theta_n_12 = %f\n",theta_n_12);
            //master.print_message(" | theta_n_22 = %f\n",theta_n_22);
            //master.print_message(" | delta_q_11 = %f\n",delta_q_11);
            //master.print_message(" | delta_q_12 = %f\n",delta_q_12);
            //master.print_message(" | delta_q_22 = %f\n",delta_q_22);
            //master.print_message(" | theta_q_11 = %f\n",theta_q_11);
            //master.print_message(" | theta_q_12 = %f\n",theta_q_12);
            //master.print_message(" | theta_q_22 = %f\n",theta_q_22);
            master.print_message(" | delta_n = %f\n",ice_coeffs.sc_delta_n);
            master.print_message(" | theta_n = %f\n",ice_coeffs.sc_theta_n);
            master.print_message(" | delta_q = %f\n",ice_coeffs.sc_delta_q);
            master.print_message(" | theta_q = %f\n",ice_coeffs.sc_theta_q);
        }
    };

    init_2mom_scheme();

    //ice_typ   = cloud_type/1000           ! (0) no ice, (1) no hail (2) with hail
    //nuc_i_typ = MOD(cloud_type/100,10)    ! choice of ice nucleation, see ice_nucleation_homhet()
    //nuc_c_typ = MOD(cloud_type/10,10)     ! choice of CCN assumptions, see cloud_nucleation()
    //auto_typ  = MOD(cloud_type,10)        ! choice of warm rain scheme, see clouds_twomoment()

    // Set the rain_coeff types to the default provided values
    rain_coeffs = rainSBBcoeffs;

    // Moved function to after assignment to prevent overwriting.
    setup_particle_coeffs(rain, rain_coeffs);
    // CvH: ICON overrides using the following code, but they are as far as I see the same values.
    // rain_coeffs%cmu0 = cfg_params%rain_cmu0
    // rain_coeffs%cmu1 = cfg_params%rain_cmu1
    // rain_coeffs%cmu3 = cfg_params%rain_cmu3kkk

    if (sw_debug)
    {
        master.print_message("Precipitation types:\n", this->cloud_type);
        master.print_message(" | cloud_type = %d\n", this->cloud_type);
        master.print_message(" | cloud = %s\n", cloud.name.c_str());
        master.print_message(" | rain = %s\n", rain.name.c_str());
        master.print_message(" | ice = %s\n", ice.name.c_str());
        master.print_message(" | snow = %s\n", snow.name.c_str());
        master.print_message(" | graupel = %s\n", graupel.name.c_str());
        master.print_message(" | hail = %s\n", hail.name.c_str());
    }

    // initialize bulk sedimentation velocities
    // calculates coeff_alfa_n, coeff_alfa_q, and coeff_lambda
    if (sw_debug)
        master.print_message("Sedimentation velocity coeffs:\n");
    init_2mom_sedi_vel(ice, ice_coeffs);
    init_2mom_sedi_vel(snow, snow_coeffs);
    init_2mom_sedi_vel(graupel, graupel_coeffs);
    init_2mom_sedi_vel(hail, hail_coeffs);

    // Look-up table and parameters for rain_freeze_gamlook
    const int nlookup   = 2000;      // Internal number of bins (low res part)
    const int nlookuphr_dummy = 10;  // Dummy if HR part is not needed.

    this->rain_nm1 = (rain.nu+TF(1)) / rain.mu;
    this->rain_nm2 = (rain.nu+TF(2)) / rain.mu;
    this->rain_nm3 = (rain.nu+TF(3)) / rain.mu;

    Sb_init::incgfct_lower_lookupcreate(rain_nm1, rain_ltable1, nlookup, nlookuphr_dummy);
    Sb_init::incgfct_lower_lookupcreate(rain_nm2, rain_ltable2, nlookup, nlookuphr_dummy);
    Sb_init::incgfct_lower_lookupcreate(rain_nm3, rain_ltable3, nlookup, nlookuphr_dummy);

    this->rain_g1 = rain_ltable1.igf[rain_ltable1.n-1]; // ordinary gamma function of nm1 is the last value in table 1
    this->rain_g2 = rain_ltable2.igf[rain_ltable2.n-1]; // ordinary gamma function of nm2 is the last value in table 2

    // Table and parameters for graupel_hail_conv_wet_gamlook
    this->graupel_nm1 = (graupel.nu+TF(1)) / graupel.mu;
    this->graupel_nm2 = (graupel.nu+TF(2)) / graupel.mu;

    Sb_init::incgfct_lower_lookupcreate(graupel_nm1, graupel_ltable1, nlookup, nlookuphr_dummy);
    Sb_init::incgfct_lower_lookupcreate(graupel_nm2, graupel_ltable2, nlookup, nlookuphr_dummy);

    this->graupel_g1 = graupel_ltable1.igf[graupel_ltable1.n-1]; // ordinary gamma function of nm1 is the last value in table 1
    this->graupel_g2 = graupel_ltable2.igf[graupel_ltable2.n-1]; // ordinary gamma function of nm2 is the last value in table 2

    // Other options for mue-D-relation of raindrops (for sensitivity studies)
    if (this->mu_Dm_rain_typ == 0)
    {
        // Constant mue value.
        rain_coeffs.cmu0 = TF(0.0);
        rain_coeffs.cmu1 = TF(0.0);
        rain_coeffs.cmu2 = TF(1.0);
        rain_coeffs.cmu3 = TF(1.0);
        rain_coeffs.cmu4 = (rain.nu + TF(1)) / rain.b_geo - TF(1); // This is the (constant) mue value
        rain_coeffs.cmu5 = 1;
        rain_gfak = TF(-1);  // In this case gamma = 1 in rain_evaporation
    }
    else if (this->mu_Dm_rain_typ == 1)
    {
        // Axel's mu-Dm-relation for raindrops based on 1d-bin model
        // This is the default and the cmus are set in the particle constructor.
        rain_gfak = TF(1);
    }
    else if (this->mu_Dm_rain_typ == 2)
    {
        // Modification of mu-Dm-relation for experiments with increased evaporation
        rain_coeffs.cmu0 = TF(11.0);     // instead of 6.0
        rain_coeffs.cmu1 = TF(30.0);
        rain_coeffs.cmu2 = TF(1.00e+3);
        rain_coeffs.cmu3 = TF(1.10e-3);
        rain_coeffs.cmu4 = TF(4.0);      // instead of 1.0
        rain_coeffs.cmu5 = TF(2);
        rain_gfak = TF(0.5);             // instead of 1.0
    }
    else if (this->mu_Dm_rain_typ == 3)
    {
        // Jason Milbrandts mu-Dm-relation
        rain_coeffs.cmu0 = TF(19.0);
        rain_coeffs.cmu1 = TF(19.0);     // Jason Milbrandt's mu-Dm-relation for rain_coeffs
        rain_coeffs.cmu2 = TF(0.60e+3);  // (Milbrandt&Yau 2005, JAS, Table 1)
        rain_coeffs.cmu3 = TF(1.80e-3);
        rain_coeffs.cmu4 = TF(17.0);
        rain_coeffs.cmu5 = TF(1);
        rain_gfak = TF(-1);              // In this case gamma = 1 in rain_evaporation
    }

    if (sw_debug)
    {
        master.print_message("Rain coeffs and sedimentation velocities:\n");
        master.print_message(" | name = %s\n", rain.name.c_str());
        master.print_message(" | alfa = %f\n", rain_coeffs.alfa);
        master.print_message(" | beta = %f\n", rain_coeffs.beta);
        master.print_message(" | gama = %f\n", rain_coeffs.gama);
        master.print_message(" | cmu0 = %f\n", rain_coeffs.cmu0);
        master.print_message(" | cmu1 = %f\n", rain_coeffs.cmu1);
        master.print_message(" | cmu2 = %f\n", rain_coeffs.cmu2);
        master.print_message(" | cmu3 = %f\n", rain_coeffs.cmu3);
        master.print_message(" | cmu4 = %f\n", rain_coeffs.cmu4);
        master.print_message(" | cmu5 = %f\n", rain_coeffs.cmu5);

        auto get_sedi_vel = [&](const TF x_r, const bool use_ql)
        {
            const TF nr_val = 1e-3 / x_r;
            std::vector<TF> qr = {1e-3};
            std::vector<TF> qc = {1e-3};
            std::vector<TF> nr = {nr_val};
            std::vector<TF> vq(1);
            std::vector<TF> vn(1);
            std::vector<TF> rho = {1};
            const TF rho_corr = 1;

            Sb_cold::sedi_vel_rain<TF>(
                    vq.data(),
                    vn.data(),
                    qr.data(),
                    nr.data(),
                    qc.data(),
                    rho.data(),
                    rain, rain_coeffs,
                    rho_corr,
                    0, 1,
                    0, 1,
                    0, 0, 0, use_ql);

            return std::pair<TF, TF>{vq[0], vn[0]};
        };

        std::pair<TF,TF> out_cloud_min = get_sedi_vel(rain.x_min, false);
        std::pair<TF,TF> out_cloud_max = get_sedi_vel(rain.x_max, false);

        master.print_message("Sedimentation out-of-cloud:\n");
        master.print_message(" | vn_rain_min = %f\n", out_cloud_min.second);
        master.print_message(" | vn_rain_max = %f\n", out_cloud_max.second);
        master.print_message(" | vq_rain_min = %f\n", out_cloud_min.first);
        master.print_message(" | vq_rain_max = %f\n", out_cloud_max.first);

        std::pair<TF,TF> in_cloud_min = get_sedi_vel(rain.x_min, true);
        std::pair<TF,TF> in_cloud_max = get_sedi_vel(rain.x_max, true);

        master.print_message("Sedimentation in-cloud:\n");
        master.print_message(" | vn_rain_min = %f\n", in_cloud_min.second);
        master.print_message(" | vn_rain_max = %f\n", in_cloud_max.second);
        master.print_message(" | vq_rain_min = %f\n", in_cloud_min.first);
        master.print_message(" | vq_rain_max = %f\n", in_cloud_max.first);
    }

    // Setup riming coefficients.
    setup_particle_collection_type1(snow, cloud, scr_coeffs, "snow", "cloud");
    setup_particle_collection_type2(snow, rain, srr_coeffs, "snow", "rain");
    setup_particle_collection_type2(ice, rain, irr_coeffs, "ice", "rain");
    setup_particle_collection_type1(ice, cloud, icr_coeffs, "ice", "cloud");
    setup_particle_collection_type1(hail, rain, hrr_coeffs, "hail", "rain");
    setup_particle_collection_type1(graupel, rain, grr_coeffs, "graupel", "rain", false);
    setup_particle_collection_type1(hail, cloud, hcr_coeffs, "hail", "cloud", false);
    setup_particle_collection_type1(graupel, cloud, gcr_coeffs, "graupel", "cloud", false);
    setup_particle_collection_type1(snow, ice, sic_coeffs, "snow", "ice", false);
    setup_particle_collection_type1(hail, ice, hic_coeffs, "hail", "ice", false);
    setup_particle_collection_type1(graupel, ice, gic_coeffs, "graupel", "ice");
    setup_particle_collection_type1(hail, snow, hsc_coeffs, "hail", "snow", false);
    setup_particle_collection_type1(graupel, snow, gsc_coeffs, "graupel", "snow", false);

    // Setup particle coefficients
    setup_particle_coeffs(ice, ice_coeffs);
    setup_particle_coeffs(graupel, graupel_coeffs);
    setup_particle_coeffs(hail, hail_coeffs);
    setup_particle_coeffs(snow, snow_coeffs);

    // Setup selfcollection of ice particles, coeffs are stored in their derived types
    setup_graupel_selfcollection(graupel, graupel_coeffs);
    setup_snow_selfcollection(snow, snow_coeffs);
    setup_ice_selfcollection(ice, ice_coeffs);

    // Setup run-time coeffs for cloud, e.g., used in cloud_freeze and autoconversionSB
    setup_particle_coeffs(cloud, cloud_coeffs);
    Sb_cold::setup_cloud_autoconversion(cloud, cloud_coeffs);

    if (sw_debug)
    {
        master.print_message("rain_coeffs:\n");
        master.print_message(" | c_z = %f\n", rain_coeffs.c_z);
        master.print_message("cloud_coeffs:\n");
        master.print_message(" | c_z = %f\n", cloud_coeffs.c_z);
    }

    //! Init SK Activation table
    //IF (nuc_c_typ > 5) THEN
    //  call ccn_activation_sk_4d()
    //  IF (isprint) THEN
    //    CALL message(routine,"Equidistant lookup table for Segal-Khain created")
    //  ENDIF
    //END IF


    if (sw_debug)
        master.print_message("---------------------------------------\n");
}

template<typename TF>
void Microphys_sb06<TF>::init()
{
    auto& gd = grid.get_grid_data();

    for (auto& it : hydro_types)
        if (it.second.is_mass)
            it.second.precip_rate.resize(gd.ijcells);
}

template<typename TF>
void Microphys_sb06<TF>::create(
        Input& inputin, Netcdf_handle& input_nc, Timeloop<TF>& timeloop,
        Stats<TF>& stats, Cross<TF>& cross, Dump<TF>& dump, Column<TF>& column, const std::string& sim_name)
{
    const std::string group_name = "thermo";

    // Add variables to the statistics
    if (stats.get_switch())
    {
        for (auto& it : hydro_types)
        {
            // Time series
            if (it.second.is_mass)
                stats.add_time_series(it.second.name + "_rate", "Mean surface " + it.second.name + "rate", "kg m-2 s-1", group_name);
                //stats.add_time_series("r" + it.first.substr(1,1), "Mean surface " + it.second.name + "rate", "kg m-2 s-1", group_name);

            // Profiles
            stats.add_prof("v" + it.first, "Fall velocity " + it.second.name + "mass density", "m s-1", "z" , group_name);

            // Tendencies
            stats.add_tendency(*fields.st.at(it.first), "z", tend_name, tend_longname);

            // Microphysics budget
            if (sw_microbudget)
            {
                const std::string group_name_budget = "micro_budget";

                /*
                 * Warm processes
                */
                stats.add_prof("auto_qr" , "Autoconversion tendency qr",  "kg kg-1 s-1", "z", group_name_budget);
                stats.add_prof("auto_nr" , "Autoconversion tendency nr",  "kg-1 s-1", "z", group_name_budget);

                stats.add_prof("accr_qr" , "Accretion tendency qr",  "kg kg-1 s-1", "z", group_name_budget);

                stats.add_prof("scbr_nr" , "Selfcollection/breakup tendency nr",  "kg-1 s-1", "z", group_name_budget);

                stats.add_prof("evap_qr" , "Evaporation tendency qr",  "kg kg-1 s-1", "z", group_name_budget);
                stats.add_prof("evap_nr" , "Evaporation tendency nr",  "kg-1 s-1", "z", group_name_budget);

                /*
                 * Ice cold processes.
                */
            }
        }

        // Thermo tendencies
        stats.add_tendency(*fields.st.at("thl"), "z", tend_name, tend_longname);
        stats.add_tendency(*fields.st.at("qt") , "z", tend_name, tend_longname);
    }

    if (column.get_switch())
    {
        for (auto& it : hydro_types)
            if (it.second.is_mass)
                column.add_time_series("r" + it.first.substr(1,1), "Surface " + it.second.name + " rate", "kg m-2 s-1");
    }

    // Create cross-sections
    // Variables that this class can calculate/provide:
    std::vector<std::string> allowed_crossvars;
    for (auto& it : hydro_types)
    {
        allowed_crossvars.push_back("r" + it.first.substr(1,1) + "_bot");
    }

    // Cross-reference with the variables requested in the .ini file:
    crosslist = cross.get_enabled_variables(allowed_crossvars);

    // Add timers for individual kernels.
    timer.create(timeloop.get_iotime(), sim_name);
    timer.add_timing("exec_total");
    timer.add_timing("set_default_n");
    timer.add_timing("limit_sizes");
    timer.add_timing("qr_sedi_vel");
    timer.add_timing("qs_sedi_vel");
    timer.add_timing("qi_sedi_vel");
    timer.add_timing("qg_sedi_vel");
    timer.add_timing("qh_sedi_vel");
    timer.add_timing("implicit_core");
    timer.add_timing("vapor_dep");
    timer.add_timing("qi_selfc");
    timer.add_timing("qs_selfc");
    timer.add_timing("qg_selfc");
    timer.add_timing("qiqs_coll");
    timer.add_timing("qiqg_coll");
    timer.add_timing("qsqg_coll");
    timer.add_timing("qiqh_coll");
    timer.add_timing("qsqh_coll");
    timer.add_timing("qgqh_conv");
    timer.add_timing("qi_riming");
    timer.add_timing("qs_riming");
    timer.add_timing("qhqc_riming");
    timer.add_timing("qhqr_riming");
    timer.add_timing("qgqc_riming");
    timer.add_timing("qgqr_riming");
    timer.add_timing("qr_freeze");
    timer.add_timing("qi_melt");
    timer.add_timing("qs_melt");
    timer.add_timing("qg_melt");
    timer.add_timing("qh_melt");
    timer.add_timing("qi_evap");
    timer.add_timing("qg_evap");
    timer.add_timing("qh_evap");
    timer.add_timing("qr_auto");
    timer.add_timing("qr_evap");
    timer.add_timing("qr_accr");
    timer.add_timing("qr_selfc");
}

#ifndef USECUDA
template<typename TF>
void Microphys_sb06<TF>::exec(Thermo<TF>& thermo, Timeloop<TF>& timeloop, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();
    const double dt = timeloop.get_sub_time_step();

    timer.start("exec_total");

    // Get thermodynamic variables
    bool cyclic = false;
    bool is_stat = false;

    auto ql = fields.get_tmp();
    auto T = fields.get_tmp();

    thermo.get_thermo_field(*ql, "ql", cyclic, is_stat);
    thermo.get_thermo_field(*T, "T", cyclic, is_stat);

    // Hack 1; get diagnostic qi from saturation adjustment,
    // as long as we don't have prognostic ice.
    thermo.get_thermo_field(*fields.ap.at("qi"), "qi", cyclic, is_stat);

    // Hack 2; set ice number concentration to fixed value from .ini file.
    std::fill(
        fields.ap.at("ni")->fld.begin(),
        fields.ap.at("ni")->fld.end(), Ni0);

    // Hack 3; calculate q_vapor as qt-ql-qi for now...
    auto qv = fields.get_tmp_xy();

    const std::vector<TF>& p = thermo.get_basestate_vector("p");
    const std::vector<TF>& exner = thermo.get_basestate_vector("exner");

    // TMP/HACK BvS
    //const std::vector<TF>& rho = thermo.get_basestate_vector("rho");
    const std::vector<TF>& rho = fields.rhoref;

    // 2D slices for quantities shared between different kernels.
    auto rain_mass = fields.get_tmp_xy();
    auto rain_diameter = fields.get_tmp_xy();
    auto mu_r = fields.get_tmp_xy();
    auto lambda_r = fields.get_tmp_xy();

    // Setup 2D slices for implicit solver
    const int n_slices = hydro_types.size() * 9;
    if (n_slices > gd.kcells)
        throw std::runtime_error("TODO.... :-)");

    auto tmp_slices = fields.get_tmp();
    std::fill(tmp_slices->fld.begin(), tmp_slices->fld.end(), TF(0));
    int n = 0;
    for (auto& it : hydro_types)
    {
        it.second.v_sed_now = &tmp_slices->fld.data()[n*gd.ijcells]; n+=1;
        it.second.v_sed_new = &tmp_slices->fld.data()[n*gd.ijcells]; n+=1;
        it.second.flux_now = &tmp_slices->fld.data()[n*gd.ijcells]; n+=1;
        it.second.flux_new = &tmp_slices->fld.data()[n*gd.ijcells]; n+=1;
        it.second.sum = &tmp_slices->fld.data()[n*gd.ijcells]; n+=1;
        it.second.impl = &tmp_slices->fld.data()[n*gd.ijcells]; n+=1;
        it.second.slice = &tmp_slices->fld.data()[n*gd.ijcells]; n+=1;
        it.second.conversion_tend = &tmp_slices->fld.data()[n*gd.ijcells]; n+=1;
        it.second.tmp1 = &tmp_slices->fld.data()[n*gd.ijcells]; n+=1;
    }

    // Diagnostic change in qt by warm (cloud) and cold (ice) processes.
    // Used for bookkeeping difference between e.g. evaporation and sublimation.
    auto qtt_liq = fields.get_tmp_xy();
    auto qtt_ice = fields.get_tmp_xy();

    // More tmp slices :-)
    auto tmpxy1 = fields.get_tmp_xy();
    auto tmpxy2 = fields.get_tmp_xy();

    // Deposition rate ice/snow; shared between kernels.
    auto dep_rate_ice  = fields.get_tmp_xy();
    auto dep_rate_snow = fields.get_tmp_xy();

    // Dummy fields for qcloud tendencies...
    auto qct_dummy = fields.get_tmp_xy();
    auto nct_dummy = fields.get_tmp_xy();
    auto nc_dummy = fields.get_tmp_xy();

    std::fill((*nc_dummy).begin(), (*nc_dummy).end(), this->Nc0);

    // Tmp slices for rime rates
    auto rime_rate_qc = fields.get_tmp_xy();
    auto rime_rate_nc = fields.get_tmp_xy();
    auto rime_rate_qi = fields.get_tmp_xy();
    auto rime_rate_qs = fields.get_tmp_xy();
    auto rime_rate_qr = fields.get_tmp_xy();
    auto rime_rate_nr = fields.get_tmp_xy();

    // Help functions to zero XY fields.
    auto zero_tmp_xy = [&](std::shared_ptr<std::vector<TF>>& fld_xy)
    {
        std::fill((*fld_xy).begin(), (*fld_xy).end(), TF(0));
    };

    auto convert_units_short = [&](TF* data_ptr, const bool is_to_kgm3)
    {
        Sb_common::convert_unit(
                data_ptr,
                rho.data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells,
                is_to_kgm3);
    };

    auto check = [&](const std::string& name, const int k)
    {
        /*
           After each process, check if sum of all tendencies is zero.
           If not, total water mass is not conserved...
        */

        if (!sw_debug)
            return;

        TF dtq_sum = TF(0);

        for (int j=gd.jstart; j<gd.jend; ++j)
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ij = i + j*gd.icells;

                dtq_sum += (*qtt_liq)[ij] + (*qtt_ice)[ij];

                for (auto& it : hydro_types)
                    if (it.first != "qi" && it.second.is_mass)
                        dtq_sum += it.second.conversion_tend[ij];
            }

        if (std::abs(dtq_sum) > 1e-16)
        {
            std::cout << "ERROR, SB06 water not conserved after " << name << ", sum dqx/dt = " << dtq_sum << std::endl;
            throw 1;
        }
    };

    // Convert all units from `kg kg-1` to `kg m-3` (mass) and `kg-1` to `m-3` (density).
    const bool to_kgm3 = true;
    convert_units_short(ql->fld.data(), to_kgm3);
    convert_units_short(fields.ap.at("qi")->fld.data(), to_kgm3);
    convert_units_short(fields.ap.at("qt")->fld.data(), to_kgm3);

    for (auto& it : hydro_types)
        convert_units_short(fields.ap.at(it.first)->fld.data(), to_kgm3);

    // Set to default values where qnx=0 and qx0
    timer.start("set_default_n");
    Sb_cold::set_default_n(
            fields.ap.at("qi")->fld.data(),
            fields.ap.at("ni")->fld.data(),
            fields.ap.at("qr")->fld.data(),
            fields.ap.at("nr")->fld.data(),
            fields.ap.at("qs")->fld.data(),
            fields.ap.at("ns")->fld.data(),
            fields.ap.at("qg")->fld.data(),
            fields.ap.at("ng")->fld.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);
    timer.stop("set_default_n");

    // NOTE BvS: in ICON, the size limits are set at the end of the chain of micro routines.
    // We have to do it at the start, since we don't integrate the fields in this exec() function.
    auto limit_sizes_wrapper = [&](
            TF* const restrict nx, const TF* const restrict qx, Particle<TF>& particle)
    {
        Sb_cold::limit_sizes(
                nx, qx, particle,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    };

    // size limits for all hydrometeors
    //IF (nuc_c_typ > 0) THEN
    //   DO k=kstart,kend
    //    DO i=istart,iend
    //      cloud%n(i,k) = MIN(cloud%n(i,k), cloud%q(i,k)/cloud%x_min)
    //      cloud%n(i,k) = MAX(cloud%n(i,k), cloud%q(i,k)/cloud%x_max)
    //      ! Hard upper limit for cloud number conc.
    //      cloud%n(i,k) = MIN(cloud%n(i,k), 5000d6)
    //    END DO
    //   END DO
    //END IF

    timer.start("limit_sizes");
    limit_sizes_wrapper(fields.ap.at("ni")->fld.data(), fields.ap.at("qi")->fld.data(), ice);
    limit_sizes_wrapper(fields.ap.at("nr")->fld.data(), fields.ap.at("qr")->fld.data(), rain);
    limit_sizes_wrapper(fields.ap.at("ns")->fld.data(), fields.ap.at("qs")->fld.data(), snow);
    limit_sizes_wrapper(fields.ap.at("ng")->fld.data(), fields.ap.at("qg")->fld.data(), graupel);
    limit_sizes_wrapper(fields.ap.at("nh")->fld.data(), fields.ap.at("qh")->fld.data(), hail);
    timer.stop("limit_sizes");

    for (int k=gd.kend-1; k>=gd.kstart; --k)
    {
        // Diagnose qv into 2D slice.
        Sb_cold::diagnose_qv(
            (*qv).data(),
            fields.ap.at("qt")->fld.data(),
            ql->fld.data(),
            fields.ap.at("qi")->fld.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells, gd.ijcells,
            k);

        for (auto& it : hydro_types)
        {
            // Copy 3D fields to 2D slices, and do partial
            // integration of dynamics tendencies.
            Sb_common::copy_slice_and_integrate(
                    it.second.slice,
                    fields.sp.at(it.first)->fld.data(),
                    fields.st.at(it.first)->fld.data(),
                    rho.data(),
                    TF(dt),
                    sw_integrate,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells, k);

            // Zero slice which gathers the conversion tendencies
            for (int n=0; n<gd.ijcells; ++n)
                it.second.conversion_tend[n] = TF(0);
        }

        // Zero diagnostic qt tendencies.
        zero_tmp_xy(qtt_liq);
        zero_tmp_xy(qtt_ice);

        // Density correction fall speeds
        // In ICON, `rhocorr` is written into the cloud/rain/etc particle types as `rho_v`.
        const TF hlp = std::log(std::max(rho[k], TF(1e-6)) / Sb_cold::rho_0<TF>);
        const TF rho_corr = std::exp(-Sb_cold::rho_vel<TF>*hlp);

        // Sedimentation velocity rain species.
        timer.start("qr_sedi_vel");
        Sb_cold::sedi_vel_rain<TF>(
                hydro_types.at("qr").v_sed_now,
                hydro_types.at("nr").v_sed_now,
                hydro_types.at("qr").slice,
                hydro_types.at("nr").slice,
                ql->fld.data(),
                rho.data(),
                rain, rain_coeffs,
                rho_corr,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells, gd.ijcells,
                k, use_ql_sedi_rain);
        timer.stop("qr_sedi_vel");

        timer.start("qi_sedi_vel");
        Sb_cold::sedi_vel_sphere(
                hydro_types.at("qi").v_sed_now,
                hydro_types.at("ni").v_sed_now,
                hydro_types.at("qi").slice,
                hydro_types.at("ni").slice,
                ice, ice_coeffs,
                rho_corr,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells);
        timer.stop("qi_sedi_vel");

        timer.start("qs_sedi_vel");
        Sb_cold::sedi_vel_sphere(
                hydro_types.at("qs").v_sed_now,
                hydro_types.at("ns").v_sed_now,
                hydro_types.at("qs").slice,
                hydro_types.at("ns").slice,
                snow, snow_coeffs,
                rho_corr,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells);
        timer.stop("qs_sedi_vel");

        //if (lprogmelt) then
        //  call sedi_vel_lwf(graupel_lwf,graupel_coeffs,  &
        //       & qg(:,k),qgl(:,k),xg_now,rhocorr(:,k),vg_sedn_now,vg_sedq_now,vg_sedl_now,its,ite)
        //  call sedi_vel_lwf(hail_lwf,hail_coeffs,        &
        //       & qh(:,k),qhl(:,k),xh_now,rhocorr(:,k),vh_sedn_now,vh_sedq_now,vh_sedl_now,its,ite)
        //else
        timer.start("qg_sedi_vel");
        Sb_cold::sedi_vel_sphere(
                hydro_types.at("qg").v_sed_now,
                hydro_types.at("ng").v_sed_now,
                hydro_types.at("qg").slice,
                hydro_types.at("ng").slice,
                graupel, graupel_coeffs,
                rho_corr,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells);
        timer.stop("qg_sedi_vel");

        timer.start("qh_sedi_vel");
        Sb_cold::sedi_vel_sphere(
                hydro_types.at("qh").v_sed_now,
                hydro_types.at("nh").v_sed_now,
                hydro_types.at("qh").slice,
                hydro_types.at("nh").slice,
                hail, hail_coeffs,
                rho_corr,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells);
        timer.stop("qh_sedi_vel");
        //end if

        const TF rdzdt = TF(0.5) * gd.dzi[k] * dt;

        timer.start("implicit_core");
        for (auto& it : hydro_types)
            Sb_common::implicit_core(
                    it.second.slice,
                    it.second.sum,
                    it.second.impl,
                    it.second.v_sed_new,
                    it.second.v_sed_now,
                    it.second.flux_new,
                    it.second.flux_now,
                    rdzdt,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
        timer.stop("implicit_core");

        if (sw_warm)
        {
            /*
               Calculate microphysics processes.
               These are the "old" `2mom_warm` kernels ported from UCLA-LES/DALES.
            */

            throw std::runtime_error("\"sw_warm\" is no longer (or at least for now...) supported in SB06\n");

            //// Autoconversion; formation of rain drop by coagulating cloud droplets.
            //warm::autoconversion(
            //        hydro_types.at("qr").conversion_tend,
            //        hydro_types.at("nr").conversion_tend,
            //        hydro_types.at("qr").slice,
            //        hydro_types.at("nr").slice,
            //        ql->fld.data(),
            //        rho.data(),
            //        exner.data(),
            //        Nc0, dt,
            //        gd.istart, gd.iend,
            //        gd.jstart, gd.jend,
            //        gd.icells, gd.ijcells, k);

            //// Accretion; growth of rain droplets by collecting cloud droplets
            //warm::accretion(
            //        hydro_types.at("qr").conversion_tend,
            //        hydro_types.at("qr").slice,
            //        ql->fld.data(),
            //        rho.data(),
            //        exner.data(),
            //        dt,
            //        gd.istart, gd.iend,
            //        gd.jstart, gd.jend,
            //        gd.icells, gd.ijcells, k);

            //// Calculate quantities used by multiple kernels.
            //warm::prepare_microphysics_slice(
            //        (*rain_mass).data(),
            //        (*rain_diameter).data(),
            //        (*mu_r).data(),
            //        (*lambda_r).data(),
            //        hydro_types.at("qr").slice,
            //        hydro_types.at("nr").slice,
            //        rho.data(),
            //        gd.istart, gd.iend,
            //        gd.jstart, gd.jend,
            //        gd.icells, gd.ijcells, k);

            //// Evaporation of rain droplets.
            //warm::evaporation(
            //        hydro_types.at("qr").conversion_tend,
            //        hydro_types.at("nr").conversion_tend,
            //        hydro_types.at("qr").slice,
            //        hydro_types.at("nr").slice,
            //        ql->fld.data(),
            //        // CvH, this is temporarily switched back to the sat-adjust qi until ice is enabled.
            //        // fields.sp.at("qi")->fld.data(),
            //        qi->fld.data(),
            //        fields.sp.at("qt")->fld.data(),
            //        T->fld.data(),
            //        rho.data(),
            //        exner.data(),
            //        p.data(),
            //        (*rain_mass).data(),
            //        (*rain_diameter).data(),
            //        dt,
            //        gd.istart, gd.iend,
            //        gd.jstart, gd.jend,
            //        gd.icells, gd.ijcells, k);

            //// Selfcollection & breakup: growth of raindrops by mutual (rain-rain)
            //// coagulation, and breakup by collisions.
            //warm::selfcollection_breakup(
            //        hydro_types.at("nr").conversion_tend,
            //        hydro_types.at("nr").slice,
            //        hydro_types.at("qr").slice,
            //        rho.data(),
            //        (*rain_mass).data(),
            //        (*rain_diameter).data(),
            //        (*lambda_r).data(),
            //        dt,
            //        gd.istart, gd.iend,
            //        gd.jstart, gd.jend,
            //        gd.icells, gd.ijcells, k);
        }
        else
        {
            /*
               Calculate microphysics processes.
               These are the new kernels ported from ICON.
            */

            check("start", k);

            zero_tmp_xy(dep_rate_ice);
            zero_tmp_xy(dep_rate_snow);

                    //IF (isdebug) CALL message(TRIM(routine),'cloud_nucleation')

                    //IF (nuc_c_typ .EQ. 0) THEN
                    //   IF (isdebug) CALL message(TRIM(routine),'  ... force constant cloud droplet number')
                    //   cloud%n(:,:) = qnc_const
                    //ELSEIF (nuc_c_typ < 6) THEN
                    //   IF (isdebug) CALL message(TRIM(routine),'  ... Hande et al CCN activation')
                    //   IF (PRESENT(n_cn)) THEN
                    //      CALL finish(TRIM(routine),&
                    //           & 'Error in two_moment_mcrph: Hande et al activation not supported for progn. aerosol')
                    //   ELSE
                    //      CALL ccn_activation_hdcp2(ik_slice,atmo,cloud)
                    //   END IF
                    //ELSE
                    //   IF (isdebug) CALL message(TRIM(routine), &
                    //        & '  ... CCN activation using look-up tables according to Segal& Khain')
                    //   IF (PRESENT(n_cn)) THEN
                    //      CALL ccn_activation_sk_4d(ik_slice,ccn_coeffs,atmo,cloud,n_cn)
                    //   ELSE
                    //      CALL ccn_activation_sk_4d(ik_slice,ccn_coeffs,atmo,cloud)
                    //   END IF
                    //END IF

                    //IF (ischeck) CALL check(ik_slice,'start',cloud,rain,ice,snow,graupel,hail)

                    //! homogeneous and heterogeneous ice nucleation
                    //CALL ice_nucleation_homhet(ik_slice, use_prog_in, atmo, cloud, ice, n_inact, n_inpot)
                    //IF (ischeck) CALL check(ik_slice,'ice nucleation',cloud,rain,ice,snow,graupel,hail)

                    //! homogeneous freezing of cloud droplets
                    //CALL cloud_freeze(ik_slice, dt, cloud_coeffs, qnc_const, atmo, cloud, ice)
                    //IF (ischeck) CALL check(ik_slice,'cloud_freeze', cloud, rain, ice, snow, graupel,hail)

            // Depositional growth of all ice particles.
            // Store deposition rate of ice and snow for conversion calculation in ice_riming and snow_riming.
            timer.start("vapor_dep");
            Sb_cold::vapor_dep_relaxation(
                    // 2D Output tendencies:
                    hydro_types.at("qi").conversion_tend,
                    hydro_types.at("ni").conversion_tend,
                    hydro_types.at("qs").conversion_tend,
                    hydro_types.at("ns").conversion_tend,
                    hydro_types.at("qg").conversion_tend,
                    hydro_types.at("ng").conversion_tend,
                    hydro_types.at("qh").conversion_tend,
                    hydro_types.at("nh").conversion_tend,
                    (*qtt_ice).data(),
                    (*dep_rate_ice).data(),
                    (*dep_rate_snow).data(),
                    // 2D tmp fields:
                    (*tmpxy1).data(),
                    (*tmpxy2).data(),
                    hydro_types.at("qi").tmp1,
                    hydro_types.at("qs").tmp1,
                    hydro_types.at("qg").tmp1,
                    hydro_types.at("qh").tmp1,
                    // 2D input:
                    hydro_types.at("qi").slice,
                    hydro_types.at("ni").slice,
                    hydro_types.at("qs").slice,
                    hydro_types.at("ns").slice,
                    hydro_types.at("qg").slice,
                    hydro_types.at("ng").slice,
                    hydro_types.at("qh").slice,
                    hydro_types.at("nh").slice,
                    &T->fld.data()[k*gd.ijcells],
                    (*qv).data(),
                    p[k], TF(dt),
                    ice,
                    snow,
                    graupel,
                    hail,
                    ice_coeffs,
                    snow_coeffs,
                    graupel_coeffs,
                    hail_coeffs,
                    rho_corr,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("vapor_dep");
            check("vapor_dep_relaxation", k);

            // Ice-ice collisions -> forms snow.
            timer.start("qi_selfc");
            Sb_cold::ice_selfcollection(
                    hydro_types.at("qi").conversion_tend,
                    hydro_types.at("ni").conversion_tend,
                    hydro_types.at("qs").conversion_tend,
                    hydro_types.at("ns").conversion_tend,
                    (*qtt_ice).data(),
                    hydro_types.at("qi").slice,
                    hydro_types.at("ni").slice,
                    &T->fld.data()[k*gd.ijcells],
                    ice,
                    snow,
                    ice_coeffs,
                    t_cfg_2mom,
                    rho_corr,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qi_selfc");
            check("ice_selfcollection", k);

            // Selfcollection of snow
            timer.start("qs_selfc");
            Sb_cold::snow_selfcollection(
                    hydro_types.at("ns").conversion_tend,
                    hydro_types.at("qs").slice,
                    hydro_types.at("ns").slice,
                    &T->fld.data()[k*gd.ijcells],
                    snow,
                    snow_coeffs,
                    rho_corr,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qs_selfc");
            check("snow_selfcollection", k);

            // Selfcollection of graupel.
            timer.start("qg_selfc");
            Sb_cold::graupel_selfcollection(
                    hydro_types.at("ng").conversion_tend,
                    hydro_types.at("qg").slice,
                    hydro_types.at("ng").slice,
                    &T->fld.data()[k*gd.ijcells],
                    graupel,
                    graupel_coeffs,
                    rho_corr,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qg_selfc");
            check("graupel_selfcollection", k);

            // Collection of ice by snow.
            const bool save_ice_tendency = true;

            timer.start("qiqs_coll");
            Sb_cold::particle_particle_collection(
                    hydro_types.at("qs").conversion_tend,
                    hydro_types.at("qi").conversion_tend,
                    hydro_types.at("ni").conversion_tend,
                    (*qtt_ice).data(),
                    hydro_types.at("qi").slice,
                    hydro_types.at("qs").slice,
                    hydro_types.at("ni").slice,
                    hydro_types.at("ns").slice,
                    &T->fld.data()[k*gd.ijcells],
                    ice, snow,
                    sic_coeffs,
                    rho_corr,
                    save_ice_tendency,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qiqs_coll");
            check("particle_particle_collection snow-ice", k);

            // Collection of ice by graupel.
            timer.start("qiqg_coll");
            Sb_cold::particle_particle_collection(
                    hydro_types.at("qg").conversion_tend,
                    hydro_types.at("qi").conversion_tend,
                    hydro_types.at("ni").conversion_tend,
                    (*qtt_ice).data(),
                    hydro_types.at("qi").slice,
                    hydro_types.at("qg").slice,
                    hydro_types.at("ni").slice,
                    hydro_types.at("ng").slice,
                    &T->fld.data()[k*gd.ijcells],
                    ice, graupel,
                    gic_coeffs,
                    rho_corr,
                    save_ice_tendency,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qiqg_coll");
            check("particle_particle_collection graupel-ice", k);

            // Collection of snow by graupel.
            timer.start("qsqg_coll");
            Sb_cold::particle_particle_collection(
                    hydro_types.at("qg").conversion_tend,
                    hydro_types.at("qs").conversion_tend,
                    hydro_types.at("ns").conversion_tend,
                    (*qtt_ice).data(),
                    hydro_types.at("qs").slice,
                    hydro_types.at("qg").slice,
                    hydro_types.at("ns").slice,
                    hydro_types.at("ng").slice,
                    &T->fld.data()[k*gd.ijcells],
                    snow, graupel,
                    gsc_coeffs,
                    rho_corr,
                    !save_ice_tendency,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qsqg_coll");
            check("particle_particle_collection graupel-snow", k);

            // Collection of ice by hail.
            timer.start("qiqh_coll");
            Sb_cold::particle_particle_collection(
                    hydro_types.at("qh").conversion_tend,
                    hydro_types.at("qi").conversion_tend,
                    hydro_types.at("ni").conversion_tend,
                    (*qtt_ice).data(),
                    hydro_types.at("qi").slice,
                    hydro_types.at("qh").slice,
                    hydro_types.at("ni").slice,
                    hydro_types.at("nh").slice,
                    &T->fld.data()[k*gd.ijcells],
                    ice, hail,
                    hic_coeffs,
                    rho_corr,
                    save_ice_tendency,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qiqh_coll");
            check("particle_particle_collection hail-ice", k);

            // Collection of snow by hail.
            timer.start("qsqh_coll");
            Sb_cold::particle_particle_collection(
                    hydro_types.at("qh").conversion_tend,
                    hydro_types.at("qs").conversion_tend,
                    hydro_types.at("ns").conversion_tend,
                    (*qtt_ice).data(),
                    hydro_types.at("qs").slice,
                    hydro_types.at("qh").slice,
                    hydro_types.at("ns").slice,
                    hydro_types.at("nh").slice,
                    &T->fld.data()[k*gd.ijcells],
                    snow, hail,
                    hsc_coeffs,
                    rho_corr,
                    !save_ice_tendency,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qsqh_coll");
            check("particle_particle_collection hail-snow", k);

            // Conversion of graupel to hail in wet growth regime
            timer.start("qgqh_conv");
            Sb_cold::graupel_hail_conv_wet_gamlook(
                    hydro_types.at("qg").conversion_tend,
                    hydro_types.at("ng").conversion_tend,
                    hydro_types.at("qh").conversion_tend,
                    hydro_types.at("nh").conversion_tend,
                    &ql->fld.data()[k*gd.ijcells],
                    hydro_types.at("qr").slice,
                    hydro_types.at("qg").slice,
                    hydro_types.at("ng").slice,
                    hydro_types.at("qi").slice,
                    hydro_types.at("qs").slice,
                    &T->fld.data()[k*gd.ijcells],
                    graupel_ltable1,
                    graupel_ltable2,
                    graupel,
                    ltabdminwgg,
                    graupel_nm1,
                    graupel_nm2,
                    graupel_g1,
                    graupel_g2,
                    p[k], TF(dt),
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qgqh_conv");
            check("graupel_hail_conv_wet_gamlook", k);

            // Riming of ice with cloud droplets and rain drops, and conversion to graupel
            timer.start("qi_riming");
            Sb_cold::ice_riming(
                    (*qct_dummy).data(),
                    (*nct_dummy).data(),
                    hydro_types.at("qi").conversion_tend,
                    hydro_types.at("ni").conversion_tend,
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("nr").conversion_tend,
                    hydro_types.at("qg").conversion_tend,
                    hydro_types.at("ng").conversion_tend,
                    (*dep_rate_ice).data(),
                    (*rime_rate_qc).data(),
                    (*rime_rate_nc).data(),
                    (*rime_rate_qi).data(),
                    (*rime_rate_qr).data(),
                    (*rime_rate_nr).data(),
                    (*qtt_liq).data(),
                    (*qtt_ice).data(),
                    hydro_types.at("qi").slice,
                    hydro_types.at("ni").slice,
                    &ql->fld.data()[k*gd.ijcells],
                    (*nc_dummy).data(),
                    hydro_types.at("qr").slice,
                    hydro_types.at("nr").slice,
                    &T->fld.data()[k*gd.ijcells],
                    ice, cloud, rain, graupel,
                    icr_coeffs, irr_coeffs,
                    t_cfg_2mom,
                    rho_corr,
                    this->ice_multiplication,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qi_riming");
            check("ice_riming", k);

            // Riming of snow with cloud droplets and rain drops, and conversion to graupel
            timer.start("qs_riming");
            Sb_cold::snow_riming(
                    (*qct_dummy).data(),
                    (*nct_dummy).data(),
                    hydro_types.at("qs").conversion_tend,
                    hydro_types.at("ns").conversion_tend,
                    hydro_types.at("qi").conversion_tend,
                    hydro_types.at("ni").conversion_tend,
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("nr").conversion_tend,
                    hydro_types.at("qg").conversion_tend,
                    hydro_types.at("ng").conversion_tend,
                    (*dep_rate_snow).data(),
                    (*rime_rate_qc).data(),
                    (*rime_rate_nc).data(),
                    (*rime_rate_qs).data(),
                    (*rime_rate_qr).data(),
                    (*rime_rate_nr).data(),
                    (*qtt_liq).data(),
                    (*qtt_ice).data(),
                    hydro_types.at("qs").slice,
                    hydro_types.at("ns").slice,
                    &ql->fld.data()[k*gd.ijcells],
                    (*nc_dummy).data(),
                    hydro_types.at("qr").slice,
                    hydro_types.at("nr").slice,
                    &T->fld.data()[k*gd.ijcells],
                    snow, ice, cloud, rain, graupel,
                    scr_coeffs, srr_coeffs,
                    t_cfg_2mom,
                    rho_corr,
                    this->ice_multiplication,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qs_riming");
            check("snow_riming", k);

            // Hail-cloud riming
            timer.start("qhqc_riming");
            Sb_cold::particle_cloud_riming(
                    hydro_types.at("qh").conversion_tend,
                    hydro_types.at("nh").conversion_tend,
                    (*qct_dummy).data(),
                    (*nct_dummy).data(),
                    hydro_types.at("qi").conversion_tend,
                    hydro_types.at("ni").conversion_tend,
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("nr").conversion_tend,
                    (*qtt_ice).data(),
                    &ql->fld.data()[k*gd.ijcells],
                    (*nc_dummy).data(),
                    hydro_types.at("qh").slice,
                    hydro_types.at("nh").slice,
                    &T->fld.data()[k*gd.ijcells],
                    ice, hail, cloud, rain,
                    hcr_coeffs,
                    rho_corr,
                    this->ice_multiplication,
                    this->enhanced_melting,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qhqc_riming");
            check("particle_cloud_riming hail-cloud", k);

            // Hail-rain riming
            timer.start("qhqr_riming");
            Sb_cold::particle_rain_riming(
                    hydro_types.at("qh").conversion_tend,
                    hydro_types.at("nh").conversion_tend,
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("nr").conversion_tend,
                    hydro_types.at("qi").conversion_tend,
                    hydro_types.at("ni").conversion_tend,
                    (*qtt_ice).data(),
                    hydro_types.at("qr").slice,
                    hydro_types.at("nr").slice,
                    hydro_types.at("qh").slice,
                    hydro_types.at("nh").slice,
                    &T->fld.data()[k*gd.ijcells],
                    rain, ice, hail,
                    hrr_coeffs,
                    rho_corr,
                    this->ice_multiplication,
                    this->enhanced_melting,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qhqr_riming");
            check("particle_rain_riming hail-rain", k);

            // Graupel-cloud riming
            timer.start("qgqc_riming");
            Sb_cold::particle_cloud_riming(
                    hydro_types.at("qg").conversion_tend,
                    hydro_types.at("ng").conversion_tend,
                    (*qct_dummy).data(),
                    (*nct_dummy).data(),
                    hydro_types.at("qi").conversion_tend,
                    hydro_types.at("ni").conversion_tend,
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("nr").conversion_tend,
                    (*qtt_ice).data(),
                    &ql->fld.data()[k*gd.ijcells],
                    (*nc_dummy).data(),
                    hydro_types.at("qg").slice,
                    hydro_types.at("ng").slice,
                    &T->fld.data()[k*gd.ijcells],
                    ice, graupel, cloud, rain,
                    gcr_coeffs,
                    rho_corr,
                    this->ice_multiplication,
                    this->enhanced_melting,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qgqc_riming");
            check("particle_cloud_riming graupel-cloud", k);

            // Graupel-rain riming
            timer.start("qgqr_riming");
            Sb_cold::particle_rain_riming(
                    hydro_types.at("qg").conversion_tend,
                    hydro_types.at("ng").conversion_tend,
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("nr").conversion_tend,
                    hydro_types.at("qi").conversion_tend,
                    hydro_types.at("ni").conversion_tend,
                    (*qtt_ice).data(),
                    hydro_types.at("qr").slice,
                    hydro_types.at("nr").slice,
                    hydro_types.at("qg").slice,
                    hydro_types.at("ng").slice,
                    &T->fld.data()[k*gd.ijcells],
                    rain, ice, graupel,
                    grr_coeffs,
                    rho_corr,
                    this->ice_multiplication,
                    this->enhanced_melting,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qgqr_riming");
            check("particle_rain_riming graupel-rain", k);

            // Freezing of rain and conversion to ice/graupel/hail
            timer.start("qr_freeze");
            Sb_cold::rain_freeze_gamlook(
                    hydro_types.at("qi").conversion_tend,
                    hydro_types.at("ni").conversion_tend,
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("nr").conversion_tend,
                    hydro_types.at("qg").conversion_tend,
                    hydro_types.at("ng").conversion_tend,
                    hydro_types.at("qh").conversion_tend,
                    hydro_types.at("nh").conversion_tend,
                    (*qtt_ice).data(),
                    hydro_types.at("qr").slice,
                    hydro_types.at("nr").slice,
                    &T->fld.data()[k*gd.ijcells],
                    rain_ltable1,
                    rain_ltable2,
                    rain_ltable3,
                    rain_coeffs,
                    rain,
                    t_cfg_2mom,
                    this->rain_nm1,
                    this->rain_nm2,
                    this->rain_nm3,
                    this->rain_g1,
                    this->rain_g2,
                    TF(dt),
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qr_freeze");
            check("rain_freeze_gamlook", k);

            // Melting of ice
            timer.start("qi_melt");
            Sb_cold::ice_melting(
                    hydro_types.at("qi").conversion_tend,
                    hydro_types.at("ni").conversion_tend,
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("nr").conversion_tend,
                    (*qct_dummy).data(),
                    (*nct_dummy).data(),
                    (*qtt_ice).data(),
                    (*qtt_liq).data(),
                    hydro_types.at("qi").slice,
                    hydro_types.at("ni").slice,
                    &T->fld.data()[k*gd.ijcells],
                    ice, cloud,
                    TF(dt),
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qi_melt");
            check("ice_melting", k);

            timer.start("qs_melt");
            Sb_cold::snow_melting(
                    hydro_types.at("qs").conversion_tend,
                    hydro_types.at("ns").conversion_tend,
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("nr").conversion_tend,
                    hydro_types.at("qs").slice,
                    hydro_types.at("ns").slice,
                    &T->fld.data()[k*gd.ijcells],
                    snow_coeffs,
                    snow,
                    rho_corr,
                    TF(dt),
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qs_melt");
            check("snow_melting", k);

            // Melting of graupel and hail can be simple or LWF-based
            //SELECT TYPE (graupel)
            //TYPE IS (particle_frozen)

            timer.start("qg_melt");
            Sb_cold::graupel_melting(
                    hydro_types.at("qg").conversion_tend,
                    hydro_types.at("ng").conversion_tend,
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("nr").conversion_tend,
                    hydro_types.at("qg").slice,
                    hydro_types.at("ng").slice,
                    &T->fld.data()[k*gd.ijcells],
                    graupel_coeffs,
                    graupel,
                    rho_corr,
                    TF(dt),
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qg_melt");
            check("graupel_melting", k);

            //TYPE IS (particle_lwf)
            //  CALL prepare_melting_lwf(ik_slice, atmo, gmelting)
            //  CALL particle_melting_lwf(ik_slice, dt, graupel, rain, gmelting)
            //END SELECT

            //SELECT TYPE (hail)
            //TYPE IS (particle_frozen)

            timer.start("qh_melt");
            Sb_cold::hail_melting_simple(
                    hydro_types.at("qh").conversion_tend,
                    hydro_types.at("nh").conversion_tend,
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("nr").conversion_tend,
                    hydro_types.at("qh").slice,
                    hydro_types.at("nh").slice,
                    &T->fld.data()[k*gd.ijcells],
                    hail_coeffs,
                    hail,
                    t_cfg_2mom,
                    rho_corr,
                    TF(dt),
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qh_melt");
            check("hail_melting", k);

            //TYPE IS (particle_lwf)
            //  CALL particle_melting_lwf(ik_slice, dt, hail, rain, gmelting)
            //END SELECT

            // Evaporation from melting ice particles
            timer.start("qi_evap");
            Sb_cold::evaporation(
                    hydro_types.at("qs").conversion_tend,
                    hydro_types.at("ns").conversion_tend,
                    (*qtt_liq).data(),
                    hydro_types.at("qs").slice,
                    hydro_types.at("ns").slice,
                    (*qv).data(),
                    &T->fld.data()[k*gd.ijcells],
                    snow,
                    snow_coeffs,
                    rho_corr,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qi_evap");
            check("evaporation of snow", k);

            timer.start("qg_evap");
            Sb_cold::evaporation(
                    hydro_types.at("qg").conversion_tend,
                    hydro_types.at("ng").conversion_tend,
                    (*qtt_liq).data(),
                    hydro_types.at("qg").slice,
                    hydro_types.at("ng").slice,
                    (*qv).data(),
                    &T->fld.data()[k*gd.ijcells],
                    graupel,
                    graupel_coeffs,
                    rho_corr,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qg_evap");
            check("evaporation of graupel", k);

            timer.start("qh_evap");
            Sb_cold::evaporation(
                    hydro_types.at("qh").conversion_tend,
                    hydro_types.at("nh").conversion_tend,
                    (*qtt_liq).data(),
                    hydro_types.at("qh").slice,
                    hydro_types.at("nh").slice,
                    (*qv).data(),
                    &T->fld.data()[k*gd.ijcells],
                    hail,
                    hail_coeffs,
                    rho_corr,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);
            timer.stop("qh_evap");
            check("evaporation of hail", k);

            //! warm rain processes
            //! (using something other than SB is somewhat inconsistent and not recommended)
            //IF (auto_typ == 1) THEN
            //   CALL autoconversionKB(ik_slice, dt, cloud, rain)   ! Beheng (1994)
            //   CALL accretionKB(ik_slice, dt, cloud, rain)
            //   CALL rain_selfcollectionSB(ik_slice, dt, rain)
            //ELSE IF (auto_typ == 2) THEN
            //   ! Khairoutdinov and Kogan (2000)
            //   ! (KK2000 originally assume a 25 micron size threshold)
            //   CALL autoconversionKK(ik_slice, dt, cloud, rain)
            //   CALL accretionKK(ik_slice, dt, cloud, rain)
            //   CALL rain_selfcollectionSB(ik_slice, dt, rain)
            //ELSE IF (auto_typ == 3) THEN

            // Autoconversion; formation of rain drop by coagulating cloud droplets.
            timer.start("qr_auto");
            Sb_cold::autoconversionSB(
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("nr").conversion_tend,
                    (*qtt_liq).data(),
                    hydro_types.at("qr").slice,
                    hydro_types.at("nr").slice,
                    &ql->fld.data()[k*gd.ijcells],
                    cloud_coeffs,
                    cloud, rain,
                    rho_corr,
                    Nc0,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells,
                    k);
            timer.stop("qr_auto");
            check("autoconversionSB", k);

            timer.start("qr_accr");
            Sb_cold::accretionSB(
                    hydro_types.at("qr").conversion_tend,
                    (*qtt_liq).data(),
                    hydro_types.at("qr").slice,
                    &ql->fld.data()[k*gd.ijcells],
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells,
                    k);
            timer.stop("qr_accr");
            check("accretionSB", k);

            timer.start("qr_selfc");
            Sb_cold::rain_selfcollectionSB(
                    hydro_types.at("nr").conversion_tend,
                    hydro_types.at("qr").slice,
                    hydro_types.at("nr").slice,
                    rain,
                    rho_corr,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells,
                    k);
            timer.stop("qr_selfc");
            check("rain_selfcollectionSB", k);

            //ENDIF

            // Evaporation of rain following Seifert (2008)
            timer.start("qr_evap");
            Sb_cold::rain_evaporation(
                    hydro_types.at("qr").conversion_tend,
                    hydro_types.at("nr").conversion_tend,
                    (*qtt_liq).data(),
                    hydro_types.at("qr").slice,
                    hydro_types.at("nr").slice,
                    (*qv).data(),
                    &ql->fld.data()[k*gd.ijcells],
                    &T->fld.data()[k*gd.ijcells],
                    p.data(),
                    rain_coeffs,
                    cloud,
                    rain,
                    t_cfg_2mom,
                    rain_gfak,
                    rho_corr,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells,
                    k);
            timer.stop("qr_evap");
            check("rain_evaporation", k);
        }

        for (auto& it : hydro_types)
        {
            // Integrate conversion tendencies into qr/Nr slices before implicit step.
            Sb_common::integrate_process(
                    it.second.slice,
                    it.second.conversion_tend,
                    dt,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);

            // Implicit sedimentation step
            Sb_common::implicit_time(
                    it.second.slice,
                    it.second.sum,
                    it.second.impl,
                    it.second.v_sed_new,
                    it.second.v_sed_now,
                    it.second.flux_new,
                    rdzdt,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells);

            // Evaluate total tendencies back from implicit solver.
            Sb_common::diagnose_tendency(
                    fields.st.at(it.first)->fld.data(),
                    fields.sp.at(it.first)->fld.data(),
                    it.second.slice,
                    rho.data(),
                    dt,
                    sw_integrate,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells,
                    k);
        }

        // Calculate thermodynamic tendencies `thl` and `qt`,
        // from microphysics tendencies excluding sedimentation.
        if (sw_warm)
            Sb_common::calc_thermo_tendencies_cloud(
                    fields.st.at("thl")->fld.data(),
                    fields.st.at("qt")->fld.data(),
                    hydro_types.at("qr").conversion_tend,
                    rho.data(),
                    exner.data(),
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells,
                    k);
        else
            Sb_common::calc_thermo_tendencies_cloud_ice(
                    fields.st.at("thl")->fld.data(),
                    fields.st.at("qt")->fld.data(),
                    (*qtt_liq).data(),
                    (*qtt_ice).data(),
                    rho.data(),
                    exner.data(),
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.icells, gd.ijcells,
                    k);
    }

    for (auto& it : hydro_types)
    {
        // Store surface precipitation rates.
        if (it.second.is_mass)
            for (int n=0; n<gd.ijcells; ++n)
                it.second.precip_rate[n] = it.second.flux_new[n];

        // Convert all units back from `kg m-3` to `kg kg-1` (mass) and `m-3` to `kg-1` (density)
        convert_units_short(fields.ap.at(it.first)->fld.data(), !to_kgm3);
    }

    // Convert specific humidity from `kg m-3` to `kg kg-1`
    convert_units_short(fields.ap.at("qt")->fld.data(), !to_kgm3);

    // Calculate tendencies.
    stats.calc_tend(*fields.st.at("thl"), tend_name);
    stats.calc_tend(*fields.st.at("qt" ), tend_name);
    for (auto& it : hydro_types)
        stats.calc_tend(*fields.st.at(it.first), tend_name);

    // Save timings.
    timer.stop("exec_total");

    // Only save on non-substeps. The timings are gathered over all substeps..
    if (!timeloop.in_substep())
        timer.save(timeloop.get_time());

    // Release temporary fields.
    fields.release_tmp(ql);
    fields.release_tmp(T);
    fields.release_tmp(tmp_slices);

    fields.release_tmp_xy(rain_mass);
    fields.release_tmp_xy(rain_diameter);
    fields.release_tmp_xy(mu_r);
    fields.release_tmp_xy(lambda_r);

    fields.release_tmp_xy(qtt_liq);
    fields.release_tmp_xy(qtt_ice);
    fields.release_tmp_xy(qv);

    fields.release_tmp_xy(tmpxy1);
    fields.release_tmp_xy(tmpxy2);

    fields.release_tmp_xy(dep_rate_ice);
    fields.release_tmp_xy(dep_rate_snow);

    fields.release_tmp_xy(qct_dummy);
    fields.release_tmp_xy(nct_dummy);
    fields.release_tmp_xy(nc_dummy);

    fields.release_tmp_xy(rime_rate_qc);
    fields.release_tmp_xy(rime_rate_nc);
    fields.release_tmp_xy(rime_rate_qi);
    fields.release_tmp_xy(rime_rate_qs);
    fields.release_tmp_xy(rime_rate_qr);
    fields.release_tmp_xy(rime_rate_nr);
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
    for (auto& it : hydro_types)
        if (it.second.is_mass)
            stats.calc_stats_2d(it.second.name + "_rate", it.second.precip_rate, no_offset);

    // Profiles
    auto vq = fields.get_tmp();
    auto vn = fields.get_tmp();

    const std::vector<TF>& rho = thermo.get_basestate_vector("rho");

    for (int k=gd.kend-1; k>=gd.kstart; --k)
    {
        // Sedimentation rain
        // Density correction fall speeds
        const TF hlp = std::log(std::max(rho[k], TF(1e-6)) / Sb_cold::rho_0<TF>);
        const TF rho_corr = std::exp(-Sb_cold::rho_vel<TF> * hlp);

        Sb_cold::sedi_vel_rain<TF>(
                &vq->fld.data()[k * gd.ijcells],
                &vn->fld.data()[k * gd.ijcells],
                &fields.sp.at("qr")->fld.data()[k*gd.ijcells],
                &fields.sp.at("nr")->fld.data()[k*gd.ijcells],
                nullptr,
                rho.data(),
                rain, rain_coeffs,
                rho_corr,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells, gd.ijcells,
                k, use_ql_sedi_rain);
    }

    stats.calc_stats("vqr", *vq, no_offset, no_threshold);
    stats.calc_stats("vnr", *vn, no_offset, no_threshold);

    for (int k=gd.kend-1; k>=gd.kstart; --k)
    {
        // Sedimentation ice types
        // Density correction fall speeds
        const TF hlp = std::log(std::max(rho[k], TF(1e-6)) / Sb_cold::rho_0<TF>);
        const TF rho_corr = std::exp(-Sb_cold::rho_vel<TF> * hlp);

        Sb_cold::sedi_vel_sphere(
                &vq->fld.data()[k * gd.ijcells],
                &vn->fld.data()[k * gd.ijcells],
                &fields.sp.at("qr")->fld.data()[k*gd.ijcells],
                &fields.sp.at("nr")->fld.data()[k*gd.ijcells],
                snow, snow_coeffs,
                rho_corr,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells);
    }

    stats.calc_stats("vqs", *vq, no_offset, no_threshold);
    stats.calc_stats("vns", *vn, no_offset, no_threshold);


    fields.release_tmp(vq);
    fields.release_tmp(vn);

    //if (sw_microbudget)
    //{
    //    auto qrt = fields.get_tmp();
    //    auto nrt = fields.get_tmp();
    //    auto qtt = fields.get_tmp();

    //    auto qt_xy = fields.get_tmp_xy();
    //    auto qr_xy = fields.get_tmp_xy();
    //    auto nr_xy = fields.get_tmp_xy();

    //    auto ql = fields.get_tmp();
    //    auto qi = fields.get_tmp();
    //    auto T = fields.get_tmp();

    //    thermo.get_thermo_field(*ql, "ql", cyclic, is_stat);
    //    thermo.get_thermo_field(*qi, "qi", cyclic, is_stat);
    //    thermo.get_thermo_field(*T, "T", cyclic, is_stat);

    //    // Transform ql en qi from `kg kg-1` to `kg m-3`.
    //    for (int k=gd.kstart; k<gd.kend; ++k)
    //        for (int j=gd.jstart; j<gd.jend; ++j)
    //            for (int i=gd.istart; i<gd.iend; ++i)
    //            {
    //                const int ijk = i + j*gd.icells + k*gd.ijcells;
    //                ql->fld[ijk] *= rho[k];
    //                qi->fld[ijk] *= rho[k];
    //            }

    //    const std::vector<TF>& p = thermo.get_basestate_vector("p");
    //    const std::vector<TF>& exner = thermo.get_basestate_vector("exner");

    //    // TMP/HACK BvS
    //    //const std::vector<TF>& rho = thermo.get_basestate_vector("rho");
    //    const std::vector<TF>& rho = fields.rhoref;

    //    auto zero_fields = [&]()
    //    {
    //        std::fill(qrt->fld.begin(), qrt->fld.end(), TF(0));
    //        std::fill(nrt->fld.begin(), nrt->fld.end(), TF(0));
    //        std::fill(qtt->fld.begin(), qtt->fld.end(), TF(0));
    //    };

    //    auto set_moisture_slices = [&](const int k)
    //    {
    //         // Copy xy slices moisture, and transform from
    //         // `kg kg-1` to `kg m-3` and from `kg-1` to `m-3`.
    //        for (int j=gd.jstart; j<gd.jend; ++j)
    //            for (int i=gd.istart; i<gd.iend; ++i)
    //            {
    //                const int ij  = i + j * gd.icells;
    //                const int ijk = ij + k * gd.ijcells;

    //                (*qt_xy)[ij] = fields.sp.at("qt")->fld[ijk] * rho[k];
    //                (*qr_xy)[ij] = fields.sp.at("qr")->fld[ijk] * rho[k];
    //                (*nr_xy)[ij] = fields.sp.at("nr")->fld[ijk] * rho[k];
    //            }
    //    };

    //    auto to_kgkg = [&](std::shared_ptr<Field3d<TF>>& fld)
    //    {
    //        for (int k=gd.kstart; k<gd.kend; ++k)
    //            for (int j=gd.jstart; j<gd.jend; ++j)
    //                for (int i=gd.istart; i<gd.iend; ++i)
    //                {
    //                    const int ijk = i + j * gd.icells + k * gd.ijcells;
    //                    fld->fld[ijk] /= rho[k];
    //                }
    //    };

    //    std::vector<TF> rho_corr(gd.kcells);
    //    for (int k=gd.kstart; k<gd.kend; ++k)
    //    {
    //        const TF hlp = std::log(std::max(rho[k], TF(1e-6)) / Sb_cold::rho_0<TF>);
    //        rho_corr[k] = std::exp(-Sb_cold::rho_vel<TF>*hlp);
    //    }

    //    // Autoconversion.
    //    zero_fields();

    //    for (int k=gd.kstart; k<gd.kend; ++k)
    //    {
    //        set_moisture_slices(k);

    //        Sb_cold::autoconversionSB(
    //                &qrt->fld.data()[k*gd.ijcells],
    //                &nrt->fld.data()[k*gd.ijcells],
    //                &qtt->fld.data()[k*gd.ijcells],
    //                (*qr_xy).data(),
    //                (*nr_xy).data(),
    //                &ql->fld.data()[k*gd.ijcells],
    //                cloud_coeffs,
    //                cloud, rain,
    //                rho_corr[k],
    //                Nc0,
    //                gd.istart, gd.iend,
    //                gd.jstart, gd.jend,
    //                gd.icells, gd.ijcells,
    //                k);
    //    }

    //    to_kgkg(qrt);
    //    to_kgkg(nrt);

    //    stats.calc_stats("auto_qr", *qrt, no_offset, no_threshold);
    //    stats.calc_stats("auto_nr", *nrt, no_offset, no_threshold);

    //    // Accretion
    //    zero_fields();

    //    for (int k=gd.kstart; k<gd.kend; ++k)
    //    {
    //        set_moisture_slices(k);

    //        Sb_cold::accretionSB(
    //                &qrt->fld.data()[k*gd.ijcells],
    //                &qtt->fld.data()[k*gd.ijcells],
    //                (*qr_xy).data(),
    //                &ql->fld.data()[k*gd.ijcells],
    //                gd.istart, gd.iend,
    //                gd.jstart, gd.jend,
    //                gd.icells, gd.ijcells,
    //                k);
    //    }

    //    to_kgkg(qrt);
    //    stats.calc_stats("accr_qr", *qrt, no_offset, no_threshold);

    //    // Selfcollection and breakup
    //    zero_fields();

    //    for (int k=gd.kstart; k<gd.kend; ++k)
    //    {
    //        set_moisture_slices(k);

    //        Sb_cold::rain_selfcollectionSB(
    //                &nrt->fld.data()[k*gd.ijcells],
    //                (*qr_xy).data(),
    //                (*nr_xy).data(),
    //                rain,
    //                rho_corr[k],
    //                gd.istart, gd.iend,
    //                gd.jstart, gd.jend,
    //                gd.icells, gd.ijcells,
    //                k);
    //    }

    //    to_kgkg(nrt);
    //    stats.calc_stats("scbr_nr", *nrt, no_offset, no_threshold);

    //    // Evaporation
    //    zero_fields();

    //    for (int k=gd.kstart; k<gd.kend; ++k)
    //    {
    //        set_moisture_slices(k);

    //        Sb_cold::rain_evaporation(
    //                &qrt->fld.data()[k*gd.ijcells],
    //                &nrt->fld.data()[k*gd.ijcells],
    //                &qtt->fld.data()[k*gd.ijcells],
    //                (*qr_xy).data(),
    //                (*nr_xy).data(),
    //                (*qt_xy).data(),
    //                &ql->fld.data()[k*gd.ijcells],
    //                &qi->fld.data()[k*gd.ijcells],
    //                &T->fld.data()[k*gd.ijcells],
    //                p.data(),
    //                rain_coeffs,
    //                cloud,
    //                rain,
    //                t_cfg_2mom,
    //                rain_gfak,
    //                rho_corr[k],
    //                gd.istart, gd.iend,
    //                gd.jstart, gd.jend,
    //                gd.icells, gd.ijcells,
    //                k);
    //    }

    //    to_kgkg(qrt);
    //    to_kgkg(nrt);

    //    stats.calc_stats("evap_qr", *qrt, no_offset, no_threshold);
    //    stats.calc_stats("evap_nr", *nrt, no_offset, no_threshold);




    //    fields.release_tmp(qrt);
    //    fields.release_tmp(nrt);
    //    fields.release_tmp(qtt);

    //    fields.release_tmp_xy(qt_xy);
    //    fields.release_tmp_xy(qr_xy);
    //    fields.release_tmp_xy(nr_xy);

    //    fields.release_tmp(ql);
    //    fields.release_tmp(qi);
    //    fields.release_tmp(T);
    //}
}

#ifndef USECUDA
template<typename TF>
void Microphys_sb06<TF>::exec_column(Column<TF>& column)
{
    const TF no_offset = 0.;

    for (auto& it : hydro_types)
        if (it.second.is_mass)
            column.calc_time_series("r" + it.first.substr(1,1), it.second.precip_rate.data(), no_offset);
}
#endif

template<typename TF>
void Microphys_sb06<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    if (cross.get_switch())
    {
        for (auto& it : crosslist)
        {
            // Yikes (BvS), this is easy to break.... Any better ideas? Should we stick to
            // using the full name, like e.g. `rain_rate` or `qr_flux` or ...?
            const std::string letter = it.substr(1,1);
            cross.cross_plane(hydro_types.at("q" + letter).precip_rate.data(), it, iotime);
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

    return idt * this->cfl_max / cfl;
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
    auto& gd = grid.get_grid_data();

    for (int n=0; n<gd.ijcells; ++n)
        field[n] = TF(0);

    for (auto& it : hydro_types)
    {
        if (it.second.is_mass)
        {
            for (int n=0; n<gd.ijcells; ++n)
                field[n] += it.second.precip_rate[n];
        }
    }
}

template class Microphys_sb06<double>;
template class Microphys_sb06<float>;
