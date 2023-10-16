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

#include <cmath>

#include "constants.h"
#include "fast_math.h"
#include "thermo_moist_functions.h"
#include "timer.h"

namespace Sb_cold
{
    /*
       Kernels ported from ICON.
       All kernels use the same name as the ICON subroutines
    */

    using namespace Constants;
    namespace fm = Fast_math;
    namespace tmf = Thermo_moist_functions;

    template<typename TF> constexpr TF pi = 3.14159265359;
    template<typename TF> constexpr TF pi4 = pi<TF>/TF(4);
    template<typename TF> constexpr TF pi6 = pi<TF>/TF(6);
    template<typename TF> constexpr TF pi8 = pi<TF>/TF(8);
    template<typename TF> constexpr TF rho_w = 1.e3;                   // Density water
    template<typename TF> constexpr TF rho_i = 916.7;                  // Density ice (ICON)
    template<typename TF> constexpr TF rho_0 = 1.225;                  // SB06, p48
    template<typename TF> constexpr TF pirhow = pi<TF>*rho_w<TF>/6.;
    template<typename TF> constexpr TF rho_vel = 0.4;                  // Exponent for density correction
    template<typename TF> constexpr TF kc_autocon = 9.44e+9;           // Long-Kernel
    template<typename TF> constexpr TF K_T = 2.40e-2;                  // Thermal/heat conductivity of dry air (J/m/s/K)
    template<typename TF> constexpr TF N_sc = 0.710;                   // Schmidt-Zahl (PK, S.541)
    template<typename TF> constexpr TF n_f  = 0.333;                   // Exponent von N_sc im Vent-koeff. (PK, S.541)
    template<typename TF> constexpr TF nu_l = 1.50E-5;                 // kinematic viscosity of dry air (m^2/s)
    template<typename TF> constexpr TF D_v = 2.22e-5;                  // diff coeff of H2O vapor in dry air at tmelt (m^2/s)
    template<typename TF> constexpr TF rcpl = 3.1733;                  // cp_d / cp_l - 1
    template<typename TF> constexpr TF clw = (rcpl<TF> + TF(1)) * Constants::cp<TF>; // cp_d / cp_l - 1
    template<typename TF> constexpr TF e_3 = 6.10780000e2;             // Saettigungsdamppfdruck bei T = T_3
    template<typename TF> constexpr TF N_avo = 6.02214179e23;          // Avogadro constant [1/mo]
    template<typename TF> constexpr TF k_b = 1.3806504e-23;            // Boltzmann constant [J/K]

    // Various parameters for collision and conversion rates
    template<typename TF> constexpr TF ecoll_min = 0.01;               // min. eff. for graupel_cloud, ice_cloud and snow_cloud
    template<typename TF> constexpr TF ecoll_gg = 0.10;                // collision efficiency for graupel selfcollection
    template<typename TF> constexpr TF ecoll_gg_wet = 0.40;            // ...... in case of wet graupel

    // Hallet-Mossop ice multiplication
    template<typename TF> constexpr TF C_mult     = 3.5e8;             // Koeff. fuer Splintering
    template<typename TF> constexpr TF T_mult_min = 265.0;             // Minimale Temp. Splintering
    template<typename TF> constexpr TF T_mult_max = 270.0;             // Maximale Temp. Splintering
    template<typename TF> constexpr TF T_mult_opt = 268.0;             // Optimale Temp. Splintering

    // Phillips et al. ice nucleation scheme, see ice_nucleation_homhet() for more details
    //template<typename TF> constexpr TF na_dust    = 160.e4;    // initial number density of dust [1/m], Phillips08 value 162e3 (never used, reset later)
    //template<typename TF> constexpr TF na_soot    =  25.e6;    // initial number density of soot [1/m], Phillips08 value 15e6 (never used, reset later)
    //template<typename TF> constexpr TF na_orga    =  30.e6;    // initial number density of organics [1/m3], Phillips08 value 177e6 (never used, reset later)
    template<typename TF> constexpr TF ni_het_max = 100.0e3;   // max number of IN between 1-10 per liter, i.e. 1d3-10d3
    template<typename TF> constexpr TF ni_hom_max = 5000.0e3;  // number of liquid aerosols between 100-5000 per liter

    // Look-up table for Phillips et al. nucleation
    const int ttmax  = 30;   // sets limit for temperature in look-up table
    const int ssmax  = 60;   // sets limit for ice supersaturation in look-up table
    const int ttstep = 2;    // increment for temperature in look-up table
    const int ssstep = 1;    // increment for ice supersaturation in look-up table

    // Even more parameters for collision and conversion rates
    template<typename TF> constexpr TF q_crit_ii = 1.000e-6;           // q-threshold for ice_selfcollection
    template<typename TF> constexpr TF D_crit_ii = 5.0e-6;             // D-threshold for ice_selfcollection
    template<typename TF> constexpr TF q_crit_r  = 1.000e-5;           // q-threshold for ice_rain_riming and snow_rain_riming
    template<typename TF> constexpr TF D_crit_r  = 100.0e-6;           // D-threshold for ice_rain_riming and snow_rain_riming
    template<typename TF> constexpr TF q_crit_fr = 1.000e-6;           // q-threshold for rain_freeze
    template<typename TF> constexpr TF q_crit_c  = 1.000e-6;           // q-threshold for cloud water
    template<typename TF> constexpr TF q_crit    = 1.000e-9;           // q-threshold elsewhere 1e-7 kg/m3 = 1e-4 g/m3 = 0.1 mg/m3
    template<typename TF> constexpr TF D_conv_sg = 200.0e-6;           // D-threshold for conversion of snow to graupel
    template<typename TF> constexpr TF D_conv_ig = 200.0e-6;           // D-threshold for conversion of ice to graupel
    template<typename TF> constexpr TF x_conv    = 0.100e-9;           // minimum mass of conversion due to riming
    template<typename TF> constexpr TF D_crit_c  = 10.00e-6;           // D-threshold for cloud drop collection efficiency
    template<typename TF> constexpr TF D_coll_c  = 40.00e-6;           // upper bound for diameter in collision efficiency

    template<typename TF> constexpr TF T_nuc  = 268.15;                // lower temperature threshold for ice nucleation, -5 C
    template<typename TF> constexpr TF T_freeze  = 273.15;             // lower temperature threshold for raindrop freezing
    template<typename TF> constexpr TF T_f  = 233.0;                   // below this temperature there is no liquid water


    template<typename TF>
    inline TF particle_meanmass(
            Particle<TF>& particle,
            const TF q, const TF n)
    {
        // Mean mass of particle, with limiters (SB06, Eq 94)
        return std::min( std::max( q/(n+TF(Constants::dsmall)), particle.x_min ), particle.x_max );
    }

    template<typename TF>
    inline TF diffusivity(const TF T, const TF p)
    {
        // Molecular diffusivity of water vapor, from ICON.
        return TF(8.7602e-5) * exp(TF(1.81)*log(T)) / p;
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
    inline TF particle_velocity(
            Particle<TF>& particle,
            const TF mean_mass)
    {
        // Mass-diameter relation (SB06, Eq 32)
        return particle.a_vel * exp(particle.b_vel * log(mean_mass)); // v = a_vel * x**b_vel
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
    void setup_cloud_autoconversion(
            Particle<TF>& cloud,
            Particle_cloud_coeffs<TF>& cloud_coeffs)
    {
        const TF nu = cloud.nu;
        const TF mu = cloud.mu;

        if (mu == 1.0)  // NOTE BvS: is this safe?
        {
            //.. see SB2001
            cloud_coeffs.k_au =
                    kc_autocon<TF> / cloud.x_max * (TF(1) / TF(20)) * (nu + TF(2)) * (nu + TF(4)) / fm::pow2(nu + TF(1));
            cloud_coeffs.k_sc = kc_autocon<TF> * (nu + TF(2)) / (nu + TF(1));
        }
        else
        {
            throw std::runtime_error("Non SB01 autoconversion/selfcollection constants not (yet) setup...");
            //!.. see Eq. (3.44) of Seifert (2002)
            //cloud_coeffs%k_au = kc_autocon / cloud%x_max * (1.0_wp / 20.0_wp)  &
            //     & * ( 2.0_wp * gfct((nu+4.0)/mu)**1                           &
            //     &            * gfct((nu+2.0)/mu)**1 * gfct((nu+1.0)/mu)**2    &
            //     &   - 1.0_wp * gfct((nu+3.0)/mu)**2 * gfct((nu+1.0)/mu)**2 )  &
            //     &   / gfct((nu+2.0)/mu)**4
            //cloud_coeffs%k_sc = kc_autocon * cloud_coeffs%c_z
        }
    }

    template<typename TF, bool sw_prognostic_ice>
    void diagnose_qv(
            TF* const restrict qv,
            const TF* const restrict qt,
            const TF* const restrict ql,
            const TF* const restrict qi,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;
                const int ijk = i + j*jstride + k*kstride;

                if (sw_prognostic_ice)
                    qv[ij] = qt[ijk] - ql[ijk];
                else
                    qv[ij] = qt[ijk] - ql[ijk] - qi[ijk];
            }
    }

    template<typename TF>
    inline TF set_qnc(const TF qc)
    {
        const TF Dmean = TF(10e-6); // Diameter of mean particle mass

        //!set_qnc = qc * 6.0_wp / (pi * rhoh2o * Dmean**3.0_wp)
        return qc * TF(6) / (Constants::pi<TF> * Constants::rho_w<TF> * std::exp(std::log(Dmean) * TF(3)));
    }

    template<typename TF>
    inline TF set_qni(const TF qi)
    {
        //const TF Dmean = TF(100e-6);
        //!set_qni  = qi / 1e-10   !  qiin / ( ( Dmean / ageo) ** (1.0_wp / bgeo) )
        //!set_qni =  5.0E+0_wp * EXP(0.304_wp *  (T_3 - T))   ! FR: Cooper (1986) used by Greg Thompson(2008)
        return qi / TF(1e-10);  //  qiin / ( exp(log(( Dmean / ageo)) * (1.0_wp / bgeo)) )
    }

    template<typename TF>
    inline TF set_qnr(const TF qr)
    {
        const TF N0r = TF(8000e3); // intercept of MP distribution

        //!set_qnr = N0r * ( qr * 6.0_wp / (pi * rhoh2o * N0r * gfct(4.0_wp)))**(0.25_wp)

        if (qr >= TF(1e-20))
            return N0r * std::exp( std::log( qr * TF(6) / (Constants::pi<TF> * Constants::rho_w<TF> * N0r * std::tgamma(TF(4)))) * TF(0.25) );
        else
            return TF(0);
    }

    template<typename TF>
    inline TF set_qns(const TF qs)
    {
        const TF N0s = TF(800.0e3);
        const TF ams = TF(0.038);  // needs to be connected to snow-type
        const TF bms = TF(2.0);

        //!set_qns = N0s * ( qs / ( ams * N0s * gfct(bms+1.0_wp)))**( 1.0_wp/(1.0_wp+bms) )

        if (qs >= TF(1e-20))
            return N0s * std::exp( std::log( qs / ( ams * N0s * std::tgamma(bms+TF(1)))) * ( TF(1)/(TF(1)+bms) ) );
        else
            return TF(0);
    }

    template<typename TF>
    inline TF set_qng(const TF qg)
    {
        const TF N0g = TF(4000.0e3);
        const TF amg = TF(169.6);     // needs to be connected to graupel-type
        const TF bmg = TF(3.1);

        //!set_qng = N0g * ( qg / ( amg * N0g * gfct(bmg+1.0_wp)))**( 1.0_wp/(1.0_wp+bmg) )

        if (qg >= TF(1e-20))
            return N0g * std::exp( std::log ( qg / ( amg * N0g * std::tgamma(bmg + TF(1)))) * ( TF(1)/(TF(1)+bmg) ) );
        else
            return TF(0);
    }

    template<typename TF>
    inline TF set_qnh(const TF qh)
    {
        const TF Dmean = TF(5e-3);    // Diameter of mean particle mass
        const TF rhob_hail = TF(750); // assumed bulk density of hail

        //!set_qnh = qh * 6.0_wp / (pi * rhob_hail * Dmean**3.0_wp)
        return qh * TF(6) / (Constants::pi<TF> * rhob_hail * std::exp(std::log(Dmean)*TF(3)) );
    }


    template<typename TF>
    void set_default_n_warm(
            TF* const restrict qr,
            TF* const restrict nr,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        // Set to a default number concentration in places with qnx = 0 and qx !=0

        const TF eps = TF(1e-3);

        // TODO for prognostic qc/nc:
        //if (qc && nc)
        //{
        //    for (int j=jstart; j<jend; ++j)
        //        for (int i=istart; i<iend; ++i)
        //        {
        //            const int ij = i + j*jstride;
        //            if (qc[ij] > TF(0) && nc[ij] < eps)
        //                nc[ij] = set_qnc(qc[ij]);
        //        }
        //}

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jstride + k*kstride;

                    if (qr[ijk] > TF(0) && nr[ijk] < eps)
                        nr[ijk] = set_qnr(qr[ijk]);
                }
    }

    template<typename TF>
    void set_default_n_cold(
            TF* const restrict qi,
            TF* const restrict ni,
            TF* const restrict qs,
            TF* const restrict ns,
            TF* const restrict qg,
            TF* const restrict ng,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        // Set to a default number concentration in places with qnx = 0 and qx !=0

        const TF eps = TF(1e-3);

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jstride + k*kstride;

                    if (qi[ijk] > TF(0) && ni[ijk] < eps)
                        ni[ijk] = set_qni(qi[ijk]);

                    if (qs[ijk] > TF(0) && ns[ijk] < eps)
                        ns[ijk] = set_qns(qs[ijk]);

                    if (qg[ijk] > TF(0) && ng[ijk] < eps)
                        ng[ijk] = set_qng(qg[ijk]);

                    // BvS: What about qh/nh? There is a `set_qnh()` function,
                    //      but is never called in ICON.
                }
    }

    template<typename TF>
    void limit_sizes(
            TF* const restrict nx,
            const TF* const restrict qx,
            Particle<TF>& particle,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jstride + k*kstride;

                    nx[ijk] = std::min(nx[ijk], qx[ijk] / particle.x_min);
                    nx[ijk] = std::max(nx[ijk], qx[ijk] / particle.x_max);
                }
    }

    template<typename TF>
    void sedi_vel_rain(
            TF* const restrict vq,
            TF* const restrict vn,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF* const restrict ql,
            const TF* const restrict rho,
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

                if (qr[ij] > q_crit<TF>)
                {
                    const TF x = particle_meanmass(rain, qr[ij], nr[ij]);
                    const TF d_m = particle_diameter(rain, x);

                    if (qc_present)
                    {
                        if (ql[ijk] >= q_crit<TF>)
                            mue = (rain.nu + TF(1.0)) / rain.b_geo - TF(1.0);
                        else
                            mue = rain_mue_dm_relation(coeffs, d_m);
                    }
                    else
                        mue = rain_mue_dm_relation(coeffs, d_m);

                    const TF d_p =
                            d_m * std::exp(TF(-1. / 3.) * std::log((mue + TF(3)) * (mue + TF(2)) * (mue + TF(1))));
                    vn[ij] = coeffs.alfa - coeffs.beta * std::exp(-(mue + TF(1)) * std::log(TF(1) + coeffs.gama * d_p));
                    vq[ij] = coeffs.alfa - coeffs.beta * std::exp(-(mue + TF(4)) * std::log(TF(1) + coeffs.gama * d_p));

                    vn[ij] *= rho_corr;
                    vq[ij] *= rho_corr;
                }
                else
                {
                    vn[ij] = TF(0);
                    vq[ij] = TF(0);
                }
            }
    }


    template<typename TF>
    void sedi_vel_sphere(
            TF* const restrict vqx,
            TF* const restrict vnx,
            const TF* const restrict qx,
            const TF* const restrict nx,
            Particle<TF>& particle,
            Particle_sphere<TF>& particle_coeffs,
            const TF rho_corr,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {

        for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ij = i + j*jstride;

                    if (qx[ij] > q_crit<TF>)
                    {
                        const TF x = particle_meanmass(particle, qx[ij], nx[ij]);
                        const TF lam = std::exp(particle.b_vel * std::log(particle_coeffs.coeff_lambda * x));

                        TF v_n = particle_coeffs.coeff_alfa_n * lam;
                        TF v_q = particle_coeffs.coeff_alfa_q * lam;

                        v_n = std::max(v_n, particle.vsedi_min);
                        v_q = std::max(v_q, particle.vsedi_min);
                        v_n = std::min(v_n, particle.vsedi_max);
                        v_q = std::min(v_q, particle.vsedi_max);

                        vnx[ij] = v_n * rho_corr;
                        vqx[ij] = v_q * rho_corr;
                    }
                    else
                    {
                        vqx[ij] = TF(0);
                        vnx[ij] = TF(0);
                    }
                }
    }


    template<typename TF>
    void autoconversionSB(
            TF* const restrict qct,
            TF* const restrict qrt,
            TF* const restrict nrt,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF* const restrict qc,
            Particle_cloud_coeffs<TF>& cloud_coeffs,
            Particle<TF>& cloud,
            Particle<TF>& rain,
            const TF cloud_rho_v,  // cloud%rho_v(i,k) in ICON, constant per layer in uHH.
            const TF Nc0,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        /*
           Autoconversion of Seifert and Beheng (2001, Atmos. Res.)
           Formation of rain by coagulating cloud droplets
        */
        const TF k_1  = 6.00e+2;   //..Parameter for Phi
        const TF k_2  = 0.68e+0;   //..Parameter fof Phi
        const TF eps  = 1.00e-25;  // NOTE BvS: dangerous for float version MicroHH.
        const TF x_s_i = TF(1) / cloud.x_max;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j * jstride;

                if (qc[ij] > q_crit<TF>)
                {
                    const TF n_c = Nc0; // cloud%n(i,k) in ICON
                    const TF x_c = particle_meanmass(cloud, qc[ij], n_c);

                    TF au = cloud_coeffs.k_au * fm::pow2(qc[ij]) * fm::pow2(x_c) * cloud_rho_v;    // NOTE `*dt` in ICON..
                    const TF tau = std::min(std::max(TF(1) - qc[ij] / (qc[ij] + qr[ij] + eps), eps), TF(0.9));
                    const TF phi = k_1 * std::pow(tau, k_2) * fm::pow3(TF(1) - std::pow(tau, k_2));
                    au = au * (TF(1) + phi / fm::pow2(TF(1) - tau));

                    nrt[ij] += au * x_s_i;
                    qrt[ij] += au;
                    qct[ij] -= au;

                    //au  = MAX(MIN(q_c,au),0.0_wp)
                    //sc  = cloud_coeffs%k_sc * q_c**2 * dt * cloud%rho_v(i,k)
                    //rain%n(i,k)  = rain%n(i,k)  + au * x_s_i
                    //rain%q(i,k)  = rain%q(i,k)  + au
                    //cloud%n(i,k) = cloud%n(i,k) - MIN(n_c,sc)
                    //cloud%q(i,k) = cloud%q(i,k) - au
                }
            }
    }


    template<typename TF>
    void accretionSB(
            TF* const restrict qct,
            TF* const restrict qrt,
            const TF* const restrict qr,
            const TF* const restrict qc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        /*
           Accretion: growth of raindrops collecting cloud droplets
           Accretion of Seifert and Beheng (2001, Atmos. Res.)
        */
        const TF k_r = 5.78;       // kernel
        const TF k_1 = 5.00e-04;   // Phi function
        const TF eps = 1.00e-25;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j * jstride;

                if (qc[ij] > TF(0) && qr[ij] > TF(0))
                {
                    // ..accretion rate of SB2001
                    const TF tau = std::min(std::max(TF(1) - qc[ij] / (qc[ij] + qr[ij] + eps), eps), TF(1));
                    const TF phi = fm::pow4(tau/(tau+k_1));
                    const TF ac  = k_r *  qc[ij] * qr[ij] * phi;  // NOTE: `*dt` in ICON..

                    qrt[ij] += ac;
                    qct[ij] -= ac;

                    //ac = MIN(q_c,ac)
                    //x_c = particle_meanmass(cloud, q_c,n_c)
                    //rain%q(i,k)  = rain%q(i,k)  + ac
                    //cloud%q(i,k) = cloud%q(i,k) - ac
                    //cloud%n(i,k) = cloud%n(i,k) - MIN(n_c,ac/x_c)
                }
            }
    }


    template<typename TF>
    void rain_selfcollectionSB(
            TF* const restrict nrt,
            const TF* const restrict qr,
            const TF* const restrict nr,
            Particle<TF>& rain,
            const TF rain_rho_v,  // rain%rho_v(i,k) in ICON, constant per layer in uHH.
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        /*
           Selfcollection & breakup: growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
           Selfcollection of Seifert and Beheng (2001, Atmos. Res.)
        */
        // Parameters based on Seifert (2008, JAS)
        const TF D_br = 1.10e-3;
        const TF k_rr = 4.33e+0;
        const TF k_br = 1.00e+3;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j * jstride;

                if (qr[ij] > TF(0))
                {
                    const TF xr = particle_meanmass(rain, qr[ij], nr[ij]);
                    const TF Dr = particle_diameter(rain, xr);

                    // Selfcollection as in SB2001
                    const TF sc = k_rr * nr[ij] * qr[ij] * rain_rho_v;  // `*dt` in ICON

                    // Breakup as in Seifert (2008, JAS), Eq. (A13)
                    TF br = TF(0);
                    if (Dr > TF(0.30e-3))
                        br = (k_br * (Dr - D_br) + TF(1)) * sc;

                    nrt[ij] += sc - br;
                }
            }
    }

    template<typename TF>
    void rain_evaporation(
            TF* const restrict qvt,
            TF* const restrict qrt,
            TF* const restrict nrt,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF* const restrict qv,
            const TF* const restrict ql,
            const TF* const restrict T,
            const TF* const restrict p,
            Particle_rain_coeffs<TF>& rain_coeffs,
            Particle<TF>& cloud,
            Particle<TF>& rain,
            const T_cfg_2mom<TF>& cfg_params,
            const TF rain_gfak,
            const TF rain_rho_v,  // rain%rho_v(i,k) in ICON, constant per layer in uHH.
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        /*
           Evaporation of rain based on Seifert (2008, J. Atmos. Sci.)
        */
        const TF D_br = TF(1.1e-3);
        const bool reduce_evaporation = true;

        const TF aa = rain_coeffs.alfa;
        const TF bb = rain_coeffs.beta;
        const TF cc = rain_coeffs.gama;

        const TF Lv2 = fm::pow2(Constants::Lv<TF>);

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j * jstride;

                const TF e_d = qv[ij] * Constants::Rv<TF> * T[ij];
                const TF e_sw = tmf::esat_liq(T[ij]); // in ICON, this calls `sat_pres_water`
                const TF s_sw = e_d / e_sw - TF(1);

                if (s_sw < TF(0) && qr[ij] > TF(0) && ql[ij] < q_crit<TF>)
                {
                    const TF D_vtp = diffusivity(T[ij], p[k]);

                    // Note that 2*pi is correct, because c_r = 1/2 is assumed
                    const TF g_d =
                        TF(2) * pi<TF> / ( Lv2 / (K_T<TF> * Constants::Rv<TF> * fm::pow2(T[ij]))
                            + Constants::Rv<TF> * T[ij] / (D_vtp * e_sw) );

                    const TF x_r = particle_meanmass(rain, qr[ij], nr[ij]);
                    const TF D_m = particle_diameter(rain, x_r);

                    // Eq. (20) of Seifert (2008)
                    TF mue;
                    if (D_m <= rain_coeffs.cmu3)
                       mue = rain_coeffs.cmu0 * tanh(
                        pow(TF(4) * rain_coeffs.cmu2 * (D_m - rain_coeffs.cmu3), rain_coeffs.cmu5)) + rain_coeffs.cmu4;
                    else
                        mue = rain_coeffs.cmu1 * tanh(
                        pow(TF(1) * rain_coeffs.cmu2 * (D_m - rain_coeffs.cmu3), rain_coeffs.cmu5)) + rain_coeffs.cmu4;

                    // Eq. (A8)
                    const TF lam = exp(TF(1)/TF(3) * log(pi6<TF> * rho_w<TF>*(mue+3.0)*(mue+2.0)*(mue+1.0)/x_r));

                    // Chebyshev approximation of Gamma(mue+5/2)/Gamma(mue+2)
                    const TF gfak =
                                   TF(0.1357940435E+01)
                        + mue * ( +TF(0.3033273220E+00)
                        + mue * ( -TF(0.1299313363E-01)
                        + mue * ( +TF(0.4002257774E-03)
                        - mue *    TF(0.4856703981E-05) ) ) );

                    const TF mm = mue + TF(5)/TF(2);

                    // Eq. (A7) rewritten with (A5) and (A9)
                    const TF f_v  =
                        rain.a_ven + rain.b_ven * pow(N_sc<TF>, n_f<TF>) * gfak *
                        sqrt(aa/nu_l<TF> * rain_rho_v / lam) *
                        (TF(1) - TF(1)/TF(2)   * (bb/aa)         * exp(mm*log(lam/(TF(1)*cc+lam)))
                               - TF(1)/TF(8)   * fm::pow2(bb/aa) * exp(mm*log(lam/(TF(2)*cc+lam)))
                               - TF(1)/TF(16)  * fm::pow3(bb/aa) * exp(mm*log(lam/(TF(3)*cc+lam)))
                               - TF(5)/TF(127) * fm::pow4(bb/aa) * exp(mm*log(lam/(TF(4)*cc+lam))) );

                    TF gamma_eva;
                    if (rain_gfak > TF(0))
                       gamma_eva = rain_gfak * (D_br/D_m) * exp(-TF(0.2)*mue); // Eq. (23)
                    else
                       gamma_eva = TF(1);

                    // Eq (A5) with (A9) and Gamma(mue+2) from (A7)
                    TF eva_q = g_d * nr[ij] * (mue + TF(1.0)) / lam * f_v * s_sw; // * dt in ICON

                    // UB: empirical correction factor to reduce evaporation of drizzle-like rain (D_m < D_br):
                    if (reduce_evaporation)
                    {
                        const TF eva_q_fak =
                            std::min(
                                std::max(
                                    (cfg_params.eva_q_fak_low +
                                    (cfg_params.eva_q_fak_high - cfg_params.eva_q_fak_low) /
                                    (cfg_params.eva_q_fak_Dbr_maxfak * D_br -
                                     cfg_params.eva_q_fak_Dbr_minfak * D_br) *
                                    (D_m - cfg_params.eva_q_fak_Dbr_minfak * D_br)),
                                cfg_params.eva_q_fak_low),
                            cfg_params.eva_q_fak_high);
                        eva_q *= eva_q_fak;
                    }

                    // Note to self: `eva_q` is a negative number at this point.
                    // Sign diff compared to ICON is caused by their `eva_q = MAX(-eva_q,0.0_wp)`.
                    eva_q = -eva_q;

                    qrt[ij] -= eva_q;
                    nrt[ij] -= gamma_eva * eva_q / x_r;
                    qvt[ij] += eva_q;

                    //const TF eva_q = MAX(-eva_q,0.0_wp)
                    //const TF eva_n = MAX(-eva_n,0.0_wp)

                    //const TF eva_q = MIN(eva_q,q_r)
                    //const TF eva_n = MIN(eva_n,n_r)

                    //const TF atmo%qv(i,k)     = atmo%qv(i,k)     + eva_q
                    //const TF rain%q(i,k) = rain%q(i,k) - eva_q
                    //const TF rain%n(i,k) = rain%n(i,k) - eva_n
                }
            }
    }


    template<typename TF>
    void evaporation(
            TF* const restrict qvt,
            TF* const restrict qxt,
            TF* const restrict nxt,
            const TF* const restrict qx,
            const TF* const restrict nx,
            const TF* const restrict qv,
            const TF* const restrict T,
            Particle<TF>& particle,
            Particle_coeffs<TF>& coeffs,
            const TF rho_v,  // rain%rho_v(i,k) in ICON, constant per layer in uHH.
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        /*
            Evaporation of melting snow/graupel/hail, see SB2006                                      *
        */

        bool reduce_melting = true;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j * jstride;

                if (qx[ij] > TF(0) && T[ij] > Constants::T0<TF>)
                {
                    const TF e_d = qv[ij] * Constants::Rv<TF> * T[ij];
                    const TF e_sw = tmf::esat_liq(T[ij]); // in ICON, this calls `sat_pres_water`
                    const TF s_sw = e_d / e_sw - TF(1);

                    //.. Eq. (37) of SB2006, note that 4*pi is correct because c is used below
                    const TF g_d = TF(4)*Constants::pi<TF> / ( fm::pow2(Constants::Lv<TF>) /
                            (K_T<TF> * Constants::Rv<TF> * fm::pow2(Constants::T0<TF>)) +
                            Constants::Rv<TF> * Constants::T0<TF> / (D_v<TF> * e_sw) );

                    const TF x = particle_meanmass(particle, qx[ij], nx[ij]);
                    const TF D = particle_diameter(particle, x);
                    const TF v = particle_velocity(particle, x) * rho_v;

                    const TF f_v  = coeffs.a_f + coeffs.b_f * sqrt(v*D);
                    TF eva_q = g_d * nx[ij] * coeffs.c_i * D * f_v * s_sw; // * dt in ICON
                    eva_q = std::max(-eva_q, TF(0));

                    //.. Complete evaporation of some of the melting frozen particles: parameterized in a way
                    //   to conserve the mean mass, similar to the case "gamma_eva = 1.0" in rain_evaporation() above:
                    if (reduce_melting)
                        nxt[ij] -= eva_q / x;

                    qxt[ij] -= eva_q;
                    qvt[ij] += eva_q;
                }
            }
    }


    template<typename TF>
    void particle_particle_collection(
            TF* const restrict qpt,
            TF* const restrict qit,
            TF* const restrict nit,
            const TF* const restrict qi,
            const TF* const restrict qp,
            const TF* const restrict ni,
            const TF* const restrict np,
            const TF* const restrict Ta,
            Particle_frozen<TF>& itype,
            Particle_frozen<TF>& ptype,
            Collection_coeffs<TF>& coeffs,
            const TF rho_v,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        /*  Most simple particle-particle collection for ice particles, e.g.,
            - graupel+ice -> graupel
            - graupel+snow -> graupel
            - hail+ice -> hail
            - hail+snow -> hail
            - snow+ice -> snow

            mass from `i/itype` is collected by `p/ptype`, such that
            both the mass and concentration of `i/itype` decreases,
            and only the mass of `p/ptype` increases.
        */

        for (int j = jstart; j < jend; j++)
            #pragma ivdep
            for (int i = istart; i < iend; i++)
            {
                const int ij = i + j * jstride;

                if (qi[ij] > Sb_cold::q_crit<TF> && qp[ij] > Sb_cold::q_crit<TF>)
                {
                    //.. Sticking efficiency of Lin (1983)
                    const TF e_coll = std::min(std::exp(TF(0.09) * (Ta[ij] - Constants::T0<TF>)), TF(1));

                    const TF xp = particle_meanmass(ptype, qp[ij], np[ij]);
                    const TF dp = particle_diameter(ptype, xp);
                    const TF vp = particle_velocity(ptype, xp) * rho_v;

                    const TF xi = particle_meanmass(itype, qi[ij], ni[ij]);
                    const TF di = particle_diameter(itype, xi);
                    const TF vi = particle_velocity(itype, xi) * rho_v;

                    // Both these terms have a `* dt` in ICON; left out since we need the tendency.
                    const TF coll_n = pi4<TF> * np[ij] * ni[ij] * e_coll
                                      * (coeffs.delta_n_aa * fm::pow2(dp)
                                       + coeffs.delta_n_ab * dp * di
                                       + coeffs.delta_n_bb * fm::pow2(di))
                                  * sqrt(coeffs.theta_n_aa * fm::pow2(vp)
                                       - coeffs.theta_n_ab * vp * vi
                                       + coeffs.theta_n_bb * fm::pow2(vi)
                                       + fm::pow2(itype.s_vel));

                    const TF coll_q = pi4<TF> * np[ij] * qi[ij] * e_coll
                                      * (coeffs.delta_q_aa * fm::pow2(dp)
                                       + coeffs.delta_q_ab * dp * di
                                       + coeffs.delta_q_bb * fm::pow2(di))
                                  * sqrt(coeffs.theta_q_aa * fm::pow2(vp)
                                       - coeffs.theta_q_ab * vp * vi
                                       + coeffs.theta_q_bb * fm::pow2(vi)
                                       + fm::pow2(itype.s_vel));

                    // Gain by the destination.
                    qpt[ij] += coll_q;

                    // Loss by the source.
                    qit[ij] -= coll_q;
                    nit[ij] -= coll_n;

                    //coll_n = MIN(n_i, coll_n)
                    //coll_q = MIN(q_i, coll_q)

                    //ptype%q(i,k) = ptype%q(i,k) + coll_q
                    //itype%q(i,k)  = itype%q(i,k)  - coll_q
                    //itype%n(i,k)  = itype%n(i,k)  - coll_n
                }
            }
    }

    template<typename TF>
    void vapor_deposition_generic(
            TF* const restrict dep_q,
            const TF* const restrict prtcl_q,
            const TF* const restrict prtcl_n,
            const TF* const restrict g_i,
            const TF* const restrict s_si,
            Particle<TF>& prtcl,
            Particle_coeffs<TF>& coeffs,
            const TF rho_v,
            const TF dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {

        for (int j = jstart; j < jend; j++)
            #pragma ivdep
            for (int i = istart; i < iend; i++)
            {
                const int ij = i + j*jstride;

                if (prtcl_q[ij] == TF(0))
                    dep_q[ij] = TF(0);
                else
                {
                    const TF x = particle_meanmass(prtcl, prtcl_q[ij], prtcl_n[ij]);
                    const TF d = particle_diameter(prtcl, x);
                    const TF v = particle_velocity(prtcl, x) * rho_v;

                    TF f_v  = coeffs.a_f + coeffs.b_f * sqrt(d*v);
                    f_v = std::max(f_v, coeffs.a_f/prtcl.a_ven);

                    dep_q[ij] = g_i[ij] * prtcl_n[ij] * coeffs.c_i * d * f_v * s_si[ij] * dt;
                }
            }
    }

    template<typename TF>
    void vapor_dep_relaxation(
            // 2D Output tendencies:
            TF* const restrict qvt,
            TF* const restrict qit,
            TF* const restrict nit,
            TF* const restrict qst,
            TF* const restrict nst,
            TF* const restrict qgt,
            TF* const restrict ngt,
            TF* const restrict qht,
            TF* const restrict nht,
            // Integrated deposition (so *dt):
            TF* const restrict dep_rate_ice,
            TF* const restrict dep_rate_snow,
            // 2D tmp fields:
            TF* const restrict s_si,
            TF* const restrict g_i,
            TF* const restrict dep_ice,
            TF* const restrict dep_snow,
            TF* const restrict dep_graupel,
            TF* const restrict dep_hail,
            // 2D input:
            const TF* const restrict qi,
            const TF* const restrict ni,
            const TF* const restrict qs,
            const TF* const restrict ns,
            const TF* const restrict qg,
            const TF* const restrict ng,
            const TF* const restrict qh,
            const TF* const restrict nh,
            const TF* const restrict T,
            const TF* const restrict qv,
            const TF p,
            const TF dt,
            Particle<TF>& ice,
            Particle<TF>& snow,
            Particle<TF>& graupel,
            Particle<TF>& hail,
            Particle_sphere<TF>& ice_coeffs,
            Particle_sphere<TF>& snow_coeffs,
            Particle_sphere<TF>& graupel_coeffs,
            Particle_sphere<TF>& hail_coeffs,
            const TF rho_v,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        const TF eps  = Constants::dtiny;   // 1e-20 in ICON

        // UB: if this new parameterization of n-reduction during sublimation
        //     really makes sense, move to a global constant or into the particle types.
        const bool reduce_sublimation = true;
        const TF dep_n_fac = TF(0.5);

        for (int j = jstart; j < jend; j++)
            #pragma ivdep
            for (int i = istart; i < iend; i++)
            {
                const int ij = i + j*jstride;

                if (T[ij] > Constants::T0<TF>)
                {
                    const TF e_d  = qv[ij] * Constants::Rv<TF> * T[ij];
                    const TF e_si = tmf::esat_ice(T[ij]);    // e_es(T_a) in ICON = sat_pres_ice
                    s_si[ij] = e_d / e_si - TF(1);           // Supersaturation over ice
                    const TF D_vtp = diffusivity(T[ij], p);  // D_v = 8.7602e-5 * T_a**(1.81) / p_a
                    g_i[ij] = TF(4) * Constants::pi<TF> / ( fm::pow2(Constants::Ls<TF>) /
                            (K_T<TF> * Constants::Rv<TF> * fm::pow2(T[ij])) +
                            Constants::Rv<TF> * T[ij] / (D_vtp * e_si) );
                }
                else
                {
                    g_i[ij]  = TF(0);
                    s_si[ij] = TF(0);
                }

                // Zero tmp slices
                dep_ice[ij] = TF(0);
                dep_snow[ij] = TF(0);
                dep_graupel[ij] = TF(0);
                dep_hail[ij] = TF(0);
            }

        Sb_cold::vapor_deposition_generic(
                dep_ice, qi, ni,
                g_i, s_si,
                ice, ice_coeffs,
                rho_v, dt,
                istart, iend,
                jstart, jend,
                jstride);

        Sb_cold::vapor_deposition_generic(
                dep_snow, qs, ns,
                g_i, s_si,
                snow, snow_coeffs,
                rho_v, dt,
                istart, iend,
                jstart, jend,
                jstride);

        Sb_cold::vapor_deposition_generic(
                dep_graupel, qg, ng,
                g_i, s_si,
                graupel, graupel_coeffs,
                rho_v, dt,
                istart, iend,
                jstart, jend,
                jstride);

        Sb_cold::vapor_deposition_generic(
                dep_hail, qh, nh,
                g_i, s_si,
                hail, hail_coeffs,
                rho_v, dt,
                istart, iend,
                jstart, jend,
                jstride);

        const TF zdt = TF(1) / dt;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                // Deposition only below T_3, evaporation of melting
                // particles at warmer T is treated elsewhere
                if(T[ij] > Constants::T0<TF>)
                {
                    // Depositional growth with relaxation time-scale approach based on:
                    // "A New Double-Moment Microphysics Parameterization for Application in Cloud and
                    // Climate Models. Part 1: Description" by H. Morrison, J.A.Curry, V.I. Khvorostyanov

                    const TF qvsidiff = qv[ij] - tmf::esat_ice(T[ij]) / (Constants::Rv<TF> * T[ij]);

                    if (std::abs(qvsidiff) > eps)
                    {
                        // Deposition rates are already multiplied with dt_local, therefore divide them here
                        const TF tau_i_i  = zdt / qvsidiff * dep_ice[ij];
                        const TF tau_s_i  = zdt / qvsidiff * dep_snow[ij];
                        const TF tau_g_i  = zdt / qvsidiff * dep_graupel[ij];
                        const TF tau_h_i  = zdt / qvsidiff * dep_hail[ij];

                        const TF Xi_i = ( tau_i_i + tau_s_i + tau_g_i + tau_h_i );

                        TF Xfac;
                        if (Xi_i < eps)
                            Xfac = TF(0);
                        else
                            Xfac = qvsidiff / Xi_i * (TF(1) - std::exp(-dt * Xi_i));

                        dep_ice[ij] = Xfac * tau_i_i;
                        dep_snow[ij] = Xfac * tau_s_i;
                        dep_graupel[ij] = Xfac * tau_g_i;
                        dep_hail[ij] = Xfac * tau_h_i;

                        // This limiter should not be necessary
                        if (qvsidiff < TF(0))
                        {
                            dep_ice[ij] = std::max(dep_ice[ij], -qi[ij]);
                            dep_snow[ij] = std::max(dep_snow[ij], -qs[ij]);
                            dep_graupel[ij] = std::max(dep_graupel[ij], -qg[ij]);
                            dep_hail[ij] = std::max(dep_hail[ij], -qh[ij]);
                        }

                        const TF dep_sum = dep_ice[ij] + dep_graupel[ij] + dep_snow[ij] + dep_hail[ij];

                        // NOTE BvS: this was commented out in ICON.
                        //! this limiter should not be necessary
                        //!IF (qvsidiff > 0.0_wp .and. dep_sum > qvsidiff) then
                        //!   weight = qvsidiff / dep_sum
                        //!   dep_sum          = weight * dep_sum
                        //!   dep_ice(i,k)     = weight * dep_ice(i,k)
                        //!   dep_snow(i,k)    = weight * dep_snow(i,k)
                        //!   dep_graupel(i,k) = weight * dep_graupel(i,k)
                        //!   dep_hail(i,k)    = weight * dep_hail(i,k)
                        //!END IF

                        qit[ij] += dep_ice[ij] * zdt;
                        qst[ij] += dep_snow[ij] * zdt;
                        qgt[ij] += dep_graupel[ij] * zdt;
                        qht[ij] += dep_hail[ij] * zdt;
                        qvt[ij] -= dep_sum * zdt;

                        // If deposition rate is negative, parameterize the complete evaporation of some of the
                        // particles in a way that mean size is conserved times a tuning factor < 1:
                        if (reduce_sublimation)
                        {
                            const TF x_i = particle_meanmass(ice, qi[ij], ni[ij]);
                            const TF x_s = particle_meanmass(snow, qs[ij], ns[ij]);
                            const TF x_g = particle_meanmass(graupel, qg[ij], ng[ij]);
                            const TF x_h = particle_meanmass(hail, qh[ij], nh[ij]);

                            const TF dep_ice_n = std::min(dep_ice[ij], TF(0)) / x_i;
                            const TF dep_snow_n = std::min(dep_snow[ij], TF(0)) / x_s;
                            const TF dep_graupel_n = std::min(dep_graupel[ij], TF(0)) / x_g;
                            const TF dep_hail_n = std::min(dep_hail[ij], TF(0)) / x_h;

                            nit[ij] += dep_n_fac * dep_ice_n * zdt;
                            nst[ij] += dep_n_fac * dep_snow_n * zdt;
                            ngt[ij] += dep_n_fac * dep_graupel_n * zdt;
                            nht[ij] += dep_n_fac * dep_hail_n * zdt;
                        }

                        dep_rate_ice[ij] += dep_ice[ij] * zdt;
                        dep_rate_snow[ij] += dep_snow[ij] * zdt;
                    }
                }
            }
    }

    template<typename TF>
    void ice_selfcollection(
            TF* const restrict qit,
            TF* const restrict nit,
            TF* const restrict qst,
            TF* const restrict nst,
            const TF* const restrict qi,
            const TF* const restrict ni,
            const TF* const restrict T,
            Particle_frozen<TF>& ice,
            Particle_frozen<TF>& snow,
            Particle_ice_coeffs<TF>& ice_coeffs,
            const T_cfg_2mom<TF>& cfg_params,
            const TF rho_v,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        // selfcollection of ice crystals, see SB2006 or Seifert (2002)                 *

        const TF x_conv_ii = pow(cfg_params.D_conv_ii/snow.a_geo, TF(1)/snow.b_geo);

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                const TF x_i = particle_meanmass(ice, qi[ij], ni[ij]);
                const TF D_i = particle_diameter(ice, x_i);

                if ( ni[ij] > TF(0) && qi[ij] > q_crit_ii<TF> && D_i > D_crit_ii<TF> )
                {
                    // Temperaturabhaengige Efficiency nach Cotton et al. (1986)
                    // (siehe auch Straka, 1989; S. 53)
                    const TF e_coll = std::min(
                            std::pow(TF(10), (TF(0.035) * (T[ij] - Constants::T0<TF>) - TF(0.7))), TF(0.2));

                    const TF v_i = ice.a_vel * pow(x_i, ice.b_vel) * rho_v;

                    const TF self_n = pi4<TF> * e_coll * ice_coeffs.sc_delta_n * ni[ij] * ni[ij] * D_i * D_i *
                         sqrt( ice_coeffs.sc_theta_n * v_i * v_i + TF(2) * pow(ice.s_vel, TF(2)) ); // * dt in ICON

                    const TF self_q = pi4<TF> * e_coll * ice_coeffs.sc_delta_q * ni[ij] * qi[ij] * D_i * D_i *
                         sqrt( ice_coeffs.sc_theta_q * v_i * v_i + TF(2) * pow(ice.s_vel, TF(2)) ); // * dt in ICON

                    qit[ij] -= self_q;
                    qst[ij] += self_q;

                    nit[ij] -= self_n;
                    nst[ij] += self_n / TF(2);      // BvS; why /2?

                    //self_q = MIN(self_q,q_i)
                    //self_n = MIN(MIN(self_n,self_q/x_conv_ii),n_i)

                    //ice%q(i,k)  = ice%q(i,k)  - self_q
                    //snow%q(i,k) = snow%q(i,k) + self_q

                    //ice%n(i,k)  = ice%n(i,k)  - self_n
                    //snow%n(i,k) = snow%n(i,k) + self_n / 2.0
                }
            }

    }


    template<typename TF>
    void snow_selfcollection(
            TF* const restrict nst,
            const TF* const restrict qs,
            const TF* const restrict ns,
            const TF* const restrict T,
            Particle_frozen<TF>& snow,
            Particle_snow_coeffs<TF>& snow_coeffs,
            const TF rho_v,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                if (qs[ij] > q_crit<TF>)
                {
                    //.. Temperaturabhaengige sticking efficiency nach Lin (1983)
                    const TF e_coll = std::max(TF(0.1),
                        std::min(std::exp(TF(0.09)*(T[ij]-Constants::T0<TF>)), TF(1.0)));

                    const TF x_s = particle_meanmass(snow, qs[ij], ns[ij]);
                    const TF D_s = particle_diameter(snow, x_s);
                    const TF v_s = particle_velocity(snow, x_s) * rho_v;

                    const TF self_n =
                            pi8<TF> * e_coll * ns[ij] * ns[ij] *
                            snow_coeffs.sc_delta_n * D_s * D_s *
                            sqrt(snow_coeffs.sc_theta_n * v_s * v_s + TF(2) *
                            fm::pow2(snow.s_vel) ); // * dt in ICON

                    nst[ij] -= self_n;
                }
            }
    }


    template<typename TF>
    void graupel_selfcollection(
            TF* const restrict ngt,
            const TF* const restrict qg,
            const TF* const restrict ng,
            const TF* const restrict Ta,
            Particle<TF>& graupel,
            Particle_graupel_coeffs<TF>& graupel_coeffs,
            const TF rho_v,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                if (qg[ij] > q_crit<TF>)
                {
                    //const TF n_g = graupel.n[ij];
                    const TF x_g = particle_meanmass(graupel, qg[ij], ng[ij]);
                    const TF D_g = particle_diameter(graupel, x_g);
                    const TF v_g = particle_velocity(graupel, x_g) * rho_v;

                    // Times dt in ICON
                    TF self_n = graupel_coeffs.sc_coll_n * fm::pow2(ng[ij]) * fm::pow2(D_g) * v_g;

                    // Sticking efficiency does only distinguish dry and wet based on T_3;
                    self_n = Ta[ij] > Constants::T0<TF> ? self_n * ecoll_gg_wet<TF> : self_n * ecoll_gg<TF>;

                    ngt[ij] -= self_n;
                }
            }
    }


    template<typename TF>
    void riming_cloud_core(
            TF* const restrict rime_rate_qb,
            TF* const restrict rime_rate_nb,
            const TF* const restrict qp,
            const TF* const restrict np,
            const TF* const restrict qc,
            const TF* const restrict nc,
            const TF rho_v,
            Particle_frozen<TF>& ptype,
            Particle<TF>& cloud,
            Collection_coeffs<TF>& coeffs,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        const TF const0 = TF(1)/(D_coll_c<TF> - D_crit_c<TF>);
        const TF const1 = const0 * ptype.ecoll_c;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                const TF x_p = particle_meanmass(ptype, qp[ij], np[ij]);
                const TF D_p = particle_diameter(ptype, x_p);
                const TF x_c = particle_meanmass(cloud, qc[ij], nc[ij]);
                const TF D_c = particle_diameter(cloud, x_c);

                if (qc[ij] > q_crit_c<TF> &&
                    qp[ij] > ptype.q_crit_c &&
                    D_p    > ptype.D_crit_c &&
                    D_c    > D_crit_c<TF> )
                {
                    const TF v_c = particle_velocity(cloud, x_c) * rho_v;
                    const TF v_p = particle_velocity(ptype, x_p) * rho_v;

                    const TF e_coll = std::min(
                            ptype.ecoll_c, std::max(const1*(D_c - D_crit_c<TF>), ecoll_min<TF>));

                    // Both terms are time integrated (* dt) in ICON...
                    const TF rime_n = pi4<TF> * e_coll * np[ij] * nc[ij] *
                                     ( coeffs.delta_n_aa * fm::pow2(D_p) +
                                       coeffs.delta_n_ab * D_p * D_c +
                                       coeffs.delta_n_bb * fm::pow2(D_c)) *
                                 sqrt( coeffs.theta_n_aa * fm::pow2(v_p) -
                                       coeffs.theta_n_ab * v_p * v_c +
                                       coeffs.theta_n_bb * fm::pow2(v_c) +
                                       fm::pow2(ptype.s_vel));

                    const TF rime_q = pi4<TF> * e_coll * np[ij] * qc[ij] *
                                     ( coeffs.delta_q_aa * fm::pow2(D_p) +
                                       coeffs.delta_q_ab * D_p * D_c +
                                       coeffs.delta_q_bb * fm::pow2(D_c)) *
                                 sqrt( coeffs.theta_q_aa * fm::pow2(v_p) -
                                       coeffs.theta_q_ab * v_p * v_c +
                                       coeffs.theta_q_bb * fm::pow2(v_c) +
                                       fm::pow2(ptype.s_vel));

                    rime_rate_qb[ij] = rime_q;
                    rime_rate_nb[ij] = rime_n;
                }
                else
                {
                    rime_rate_qb[ij] = TF(0);
                    rime_rate_nb[ij] = TF(0);
                }
            }
    }


    template<typename TF>
    void riming_rain_core(
            TF* const restrict rime_rate_qa,
            TF* const restrict rime_rate_qb,
            TF* const restrict rime_rate_nb,
            const TF* const restrict qa,
            const TF* const restrict na,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF rho_v,
            Particle_frozen<TF>& ptype,
            Particle<TF>& rain,
            Rain_riming_coeffs<TF>& coeffs,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                const TF x_a = particle_meanmass(ptype, qa[ij], na[ij]);
                const TF D_a = particle_diameter(ptype, x_a);

                if (qr[ij] > q_crit<TF> && qa[ij] > q_crit_r<TF> && D_a > D_crit_r<TF>)
                {
                    const TF x_r = particle_meanmass(rain, qr[ij], nr[ij]);
                    const TF D_r = particle_diameter(rain, x_r);
                    const TF v_r = particle_velocity(rain, x_r) * rho_v;
                    const TF v_a = particle_velocity(ptype, x_a) * rho_v;

                    // All three terms are time integrated (* dt) in ICON...
                    const TF rime_n = pi4<TF> * na[ij] * nr[ij] *
                                 (coeffs.delta_n_aa * D_a * D_a +
                                  coeffs.delta_n_ab * D_a * D_r +
                                  coeffs.delta_n_bb * D_r * D_r) *
                             sqrt(coeffs.theta_n_aa * v_a * v_a -
                                  coeffs.theta_n_ab * v_a * v_r +
                                  coeffs.theta_n_bb * v_r * v_r +
                                  fm::pow2(ptype.s_vel));

                    const TF rime_qr = pi4<TF> * na[ij] * qr[ij] *
                                 (coeffs.delta_n_aa * D_a * D_a +
                                  coeffs.delta_q_ab * D_a * D_r +
                                  coeffs.delta_q_bb * D_r * D_r) *
                             sqrt(coeffs.theta_n_aa * v_a * v_a -
                                  coeffs.theta_q_ab * v_a * v_r +
                                  coeffs.theta_q_bb * v_r * v_r +
                                  fm::pow2(ptype.s_vel));

                    const TF rime_qi = pi4<TF> * nr[ij] * qa[ij] *
                                 (coeffs.delta_q_aa * D_a * D_a +
                                  coeffs.delta_q_ba * D_a * D_r +
                                  coeffs.delta_n_bb * D_r * D_r) *
                             sqrt(coeffs.theta_q_aa * v_a * v_a -
                                  coeffs.theta_q_ba * v_a * v_r +
                                  coeffs.theta_n_bb * v_r * v_r +
                                  fm::pow2(ptype.s_vel));

                    rime_rate_nb[ij] = rime_n;
                    rime_rate_qa[ij] = rime_qi;
                    rime_rate_qb[ij] = rime_qr;
                }
                else
                {
                    rime_rate_nb[ij] = TF(0);
                    rime_rate_qa[ij] = TF(0);
                    rime_rate_qb[ij] = TF(0);
                }
            }
    }

    template<typename TF>
    void ice_riming(
            TF* const restrict qct,
            TF* const restrict nct,
            TF* const restrict qit,
            TF* const restrict nit,
            TF* const restrict qrt,
            TF* const restrict nrt,
            TF* const restrict qgt,
            TF* const restrict ngt,
            TF* const restrict dep_rate_ice,
            TF* const restrict rime_rate_qc,
            TF* const restrict rime_rate_nc,
            TF* const restrict rime_rate_qi,
            TF* const restrict rime_rate_qr,
            TF* const restrict rime_rate_nr,
            const TF* const restrict qi,
            const TF* const restrict ni,
            const TF* const restrict qc,
            const TF* const restrict nc,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF* const restrict Ta,
            Particle_frozen<TF>& ice,
            Particle<TF>& cloud,
            Particle<TF>& rain,
            Particle_frozen<TF>& graupel,
            Collection_coeffs<TF>& icr_coeffs,
            Rain_riming_coeffs<TF>& irr_coeffs,
            const T_cfg_2mom<TF>& cfg_params,
            const TF rho_v,
            const bool ice_multiplication,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        /*
           Riming of ice with cloud droplet and rain drops. First the process rates
           are calculated in
               snow_cloud_riming ()
               snow_rain_riming ()
           using those rates and the previously calculated and stored deposition
           rate the conversion of snow to graupel and rain is done.
        */

        const TF const3 = TF(1) / (T_mult_opt<TF> - T_mult_min<TF>);
        const TF const4 = TF(1) / (T_mult_opt<TF> - T_mult_max<TF>);

        Sb_cold::riming_cloud_core(
                rime_rate_qc,
                rime_rate_nc,
                qi, ni,
                qc, nc,
                rho_v,
                ice, cloud,
                icr_coeffs,
                istart, iend,
                jstart, jend,
                jstride);

        Sb_cold::riming_rain_core(
                rime_rate_qi,   // qa
                rime_rate_qr,   // qb
                rime_rate_nr,   // nb
                qi, ni,
                qr, nr,
                rho_v,
                ice, rain,
                irr_coeffs,
                istart, iend,
                jstart, jend,
                jstride);

        // This changes the results: !!!    const5 = rho_w/rho_ice * cfg_params%alpha_spacefilling
        const TF const5 = cfg_params.alpha_spacefilling * rho_w<TF>/rho_i<TF>;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j * jstride;

                if (dep_rate_ice[ij] > TF(0) && dep_rate_ice[ij] >= rime_rate_qc[ij]+rime_rate_qr[ij])
                {
                    // 1) Depositional growth is stronger than riming growth, therefore ice stays ice

                    // .. Ice_cloud_riming
                    if (rime_rate_qc[ij] > TF(0))
                    {
                        qit[ij] += rime_rate_qc[ij];
                        qct[ij] -= rime_rate_qc[ij];
                        nct[ij] -= rime_rate_nc[ij];

                        if (Ta[ij] < Constants::T0<TF> && ice_multiplication)
                        {
                            TF mult_1 = (Ta[ij] - T_mult_min<TF>) * const3;
                            TF mult_2 = (Ta[ij] - T_mult_max<TF>) * const4;
                            mult_1 = std::max(TF(0), std::min(mult_1, TF(1)));
                            mult_2 = std::max(TF(0), std::min(mult_2, TF(1)));
                            const TF mult_n = C_mult<TF> * mult_1 * mult_2 * rime_rate_qc[ij];

                            nit[ij] += mult_n;
                        }
                    }

                    // .. Ice_rain_riming
                    if (rime_rate_qr[ij] > TF(0))
                    {
                        qit[ij] += rime_rate_qr[ij];
                        qrt[ij] -= rime_rate_qr[ij];
                        nrt[ij] -= rime_rate_nr[ij];

                        // .. Ice multiplication
                        if (Ta[ij] < Constants::T0<TF> && ice_multiplication)
                        {
                            TF mult_1 = (Ta[ij] - T_mult_min<TF>) * const3;
                            TF mult_2 = (Ta[ij] - T_mult_max<TF>) * const4;
                            mult_1 = std::max(TF(0), std::min(mult_1, TF(1)));
                            mult_2 = std::max(TF(0), std::min(mult_2, TF(1)));
                            const TF mult_n = C_mult<TF> * mult_1 * mult_2 * rime_rate_qr[ij];

                            nit[ij] += mult_n;
                        }
                    }
                }
                else
                {
                    // 2) Depositional growth negative or smaller than riming growth, therefore ice is
                    //    allowed to convert to graupel and / or hail.

                    // Ice_cloud_riming
                    if (rime_rate_qc[ij] > TF(0))
                    {
                        const TF x_i = particle_meanmass(ice, qi[ij], ni[ij]);
                        const TF D_i = particle_diameter(ice, x_i);

                        qit[ij] += rime_rate_qc[ij];
                        qct[ij] -= rime_rate_qc[ij];
                        nct[ij] -= rime_rate_nc[ij];

                        // Ice multiplication;
                        const TF mult_q = TF(0);
                        if (Ta[ij] < Constants::T0<TF> && ice_multiplication)
                        {
                            TF mult_1 = (Ta[ij] - T_mult_min<TF>) * const3;
                            TF mult_2 = (Ta[ij] - T_mult_max<TF>) * const4;
                            mult_1 = std::max(TF(0), std::min(mult_1, TF(1)));
                            mult_2 = std::max(TF(0), std::min(mult_2, TF(1)));
                            const TF mult_n = C_mult<TF> * mult_1 * mult_2 * rime_rate_qc[ij];

                            nit[ij] += mult_n;
                        }

                        // Conversion ice -> graupel (depends on alpha_spacefilling);
                        if (D_i > D_conv_ig<TF> && Ta[ij] < cfg_params.Tmax_gr_rime)
                        {
                            const TF conv_q = (rime_rate_qc[ij] - mult_q) /
                                (const5 * (pi6<TF> * rho_i<TF> * fm::pow3(D_i) / x_i - TF(1)) );
                            const TF x_i = particle_meanmass(ice, qi[ij], ni[ij]);
                            const TF conv_n = conv_q / std::max(x_i, x_conv<TF>);

                            qit[ij] -= conv_q;
                            qgt[ij] += conv_q;

                            nit[ij] -= conv_n;
                            ngt[ij] += conv_n;
                        }
                    }

                    // Ice_rain_riming;
                    if (rime_rate_qi[ij] > TF(0))
                    {
                        nit[ij] -= rime_rate_nr[ij];
                        nrt[ij] -= rime_rate_nr[ij];

                        // BvS: I'm not sure about this part. If:
                        //      rime_rate_qi != -rime_rate_qr,
                        //      moisture is not conserved...?
                        qit[ij] -= rime_rate_qi[ij];
                        qrt[ij] -= rime_rate_qr[ij];

                        // Ice multiplication;
                        TF mult_q = TF(0);
                        TF mult_n = TF(0);

                        if (Ta[ij] < Constants::T0<TF> && ice_multiplication)
                        {
                            TF mult_1 = (Ta[ij] - T_mult_min<TF>) * const3;
                            TF mult_2 = (Ta[ij] - T_mult_max<TF>) * const4;
                            mult_1 = std::max(TF(0), std::min(mult_1, TF(1)));
                            mult_2 = std::max(TF(0), std::min(mult_2, TF(1)));
                            mult_n = C_mult<TF> * mult_1 * mult_2 * rime_rate_qr[ij];
                            mult_q = mult_n * ice.x_min;
                        }

                        if (Ta[ij] >= Constants::T0<TF>)
                        {
                            // Shedding of rain at warm temperatures;
                            // i.e. undo time integration, but with modified rain.n;
                            const TF x_r = particle_meanmass(rain, qr[ij], nr[ij]);

                            nit[ij] += rime_rate_nr[ij];
                            nrt[ij] += rime_rate_qr[ij] / x_r;

                            // BvS: I'm not sure about this part. If:
                            //      rime_rate_qi != -rime_rate_qr,
                            //      moisture is not conserved...?
                            qit[ij] += rime_rate_qi[ij];
                            qrt[ij] += rime_rate_qr[ij];
                        }
                        else
                        {
                            // New ice particles from multiplication;
                            nit[ij] += mult_n;
                            qit[ij] += mult_q;

                            // Riming to graupel;
                            if (Ta[ij] < cfg_params.Tmax_gr_rime)
                            {
                                ngt[ij] += rime_rate_nr[ij];
                                qgt[ij] += rime_rate_qi[ij] + rime_rate_qr[ij] - mult_q;
                            }
                            else
                            {
                                // Ice + frozen liquid stays ice:;
                                const TF conv_q = rime_rate_qi[ij] + rime_rate_qr[ij] - mult_q;

                                nit[ij] += rime_rate_nr[ij];
                                qit[ij] += conv_q;
                            }
                        }
                    }
                } // outer if
            } // i
    } // function


    template<typename TF>
    void snow_riming(
            TF* const restrict qct,
            TF* const restrict nct,
            TF* const restrict qst,
            TF* const restrict nst,
            TF* const restrict qit,
            TF* const restrict nit,
            TF* const restrict qrt,
            TF* const restrict nrt,
            TF* const restrict qgt,
            TF* const restrict ngt,
            TF* const restrict dep_rate_snow,
            TF* const restrict rime_rate_qc,
            TF* const restrict rime_rate_nc,
            TF* const restrict rime_rate_qs,
            TF* const restrict rime_rate_qr,
            TF* const restrict rime_rate_nr,
            const TF* const restrict qs,
            const TF* const restrict ns,
            const TF* const restrict qc,
            const TF* const restrict nc,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF* const restrict Ta,
            Particle_frozen<TF>& snow,
            Particle_frozen<TF>& ice,
            Particle<TF>& cloud,
            Particle<TF>& rain,
            Particle_frozen<TF>& graupel,
            Collection_coeffs<TF>& scr_coeffs,
            Rain_riming_coeffs<TF>& srr_coeffs,
            const T_cfg_2mom<TF>& cfg_params,
            const TF rho_v,
            const bool ice_multiplication,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        /*
            Riming of snow with cloud droplet and rain drops. First the process rates
            are calculated in
                snow_cloud_riming ()
                snow_rain_riming ()
            using those rates and the previously calculated and stored deposition
            rate the conversion of snow to graupel and rain is done.
        */
        const TF const3 = TF(1)/(T_mult_opt<TF> - T_mult_min<TF>);
        const TF const4 = TF(1)/(T_mult_opt<TF> - T_mult_max<TF>);

        Sb_cold::riming_cloud_core(
                rime_rate_qc,
                rime_rate_nc,
                qs, ns,
                qc, nc,
                rho_v,
                snow, cloud,
                scr_coeffs,
                istart, iend,
                jstart, jend,
                jstride);

        Sb_cold::riming_rain_core(
                rime_rate_qs,   // qa
                rime_rate_qr,   // qb
                rime_rate_nr,   // nb
                qs, ns,
                qr, nr,
                rho_v,
                snow, rain,
                srr_coeffs,
                istart, iend,
                jstart, jend,
                jstride);

        // This changes the results: !!!    const5 = rho_w/rho_ice * cfg_params%alpha_spacefilling
        const TF const5 = cfg_params.alpha_spacefilling * rho_w<TF>/rho_i<TF>;

        for (int j = jstart; j < jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                if (dep_rate_snow[ij] > TF(0) && dep_rate_snow[ij] >= rime_rate_qc[ij] + rime_rate_qr[ij])
                {
                    // 1) Depositional growth is stronger than riming growth, therefore snow stays snow:;

                    // Snow_cloud_riming;
                    if (rime_rate_qc[ij] > TF(0))
                    {
                        qst[ij] += rime_rate_qc[ij];
                        qct[ij] -= rime_rate_qc[ij];
                        nct[ij] -= rime_rate_nc[ij];

                        // Ice multiplication;
                        if (Ta[ij] < Constants::T0<TF> && ice_multiplication)
                        {
                            TF mult_1 = (Ta[ij] - T_mult_min<TF>) * const3;
                            TF mult_2 = (Ta[ij] - T_mult_max<TF>) * const4;

                            mult_1 = std::max(TF(0), std::min(mult_1, TF(1)));
                            mult_2 = std::max(TF(0), std::min(mult_2, TF(1)));

                            const TF mult_n = C_mult<TF> * mult_1 * mult_2 * rime_rate_qc[ij];
                            const TF mult_q = mult_n * ice.x_min;

                            nit[ij] += mult_n;
                            qit[ij] += mult_q;
                            qst[ij] -= mult_q;
                        }
                    }

                    //.. Snow_rain_riming;
                    if (rime_rate_qr[ij] > TF(0))
                    {
                        qst[ij] += rime_rate_qr[ij];
                        qrt[ij] -= rime_rate_qr[ij];
                        nrt[ij] -= rime_rate_nr[ij];

                        // Ice multiplication;
                        if (Ta[ij] < Constants::T0<TF> && ice_multiplication)
                        {
                            std::cout << "2.1" << std::endl;
                            TF mult_1 = (Ta[ij] - T_mult_min<TF>) * const3;
                            TF mult_2 = (Ta[ij] - T_mult_max<TF>) * const4;

                            mult_1 = std::max(TF(0), std::min(mult_1, TF(1)));
                            mult_2 = std::max(TF(0), std::min(mult_2, TF(1)));

                            const TF mult_n = C_mult<TF> * mult_1 * mult_2 * rime_rate_qr[ij];
                            const TF mult_q = mult_n * ice.x_min;

                            nit[ij] += mult_n;
                            qit[ij] += mult_q;
                            qst[ij] -= mult_q;
                        }
                    }
                }
                else
                {
                    // 2) Depositional growth is negative or smaller than riming growth, therefore snow is;
                    //    allowed to convert to graupel and / or hail:;

                    // snow_cloud_riming;
                    if (rime_rate_qc[ij] > TF(0))
                    {
                        const TF x_s = particle_meanmass(snow, qs[ij], ns[ij]);
                        const TF D_s = particle_diameter(snow, x_s);

                        qst[ij] += rime_rate_qc[ij];
                        qct[ij] -= rime_rate_qc[ij];
                        nct[ij] -= rime_rate_nc[ij];

                        // ice multiplication;
                        TF mult_q = TF(0);
                        if (Ta[ij] < Constants::T0<TF> && ice_multiplication)
                        {
                            TF mult_1 = (Ta[ij] - T_mult_min<TF>) * const3;
                            TF mult_2 = (Ta[ij] - T_mult_max<TF>) * const4;

                            mult_1 = std::max(TF(0), std::min(mult_1, TF(1)));
                            mult_2 = std::max(TF(0), std::min(mult_2, TF(1)));

                            const TF mult_n = C_mult<TF> * mult_1 * mult_2 * rime_rate_qc[ij];
                            mult_q = mult_n * ice.x_min;

                            nit[ij] += mult_n;
                            qit[ij] += mult_q;
                            qst[ij] -= mult_q;
                        }

                        // Conversion of snow to graupel, depends on alpha_spacefilling;
                        if (D_s > D_conv_sg<TF> && Ta[ij] < cfg_params.Tmax_gr_rime)
                        {
                            const TF conv_q = (rime_rate_qc[ij] - mult_q) /
                                ( const5 * (pi6<TF> * rho_i<TF> * fm::pow3(D_s)/x_s - TF(1)) );
                            const TF x_s = particle_meanmass(snow, qs[ij], ns[ij]);
                            const TF conv_n = conv_q / std::max(x_s, x_conv<TF>);

                            qst[ij] -= conv_q;
                            qgt[ij] += conv_q;
                            nst[ij] -= conv_n;
                            ngt[ij] += conv_n;
                        }
                    }

                    // snow_rain_riming;
                    if (rime_rate_qs[ij] > TF(0))
                    {
                        nst[ij] -= rime_rate_nr[ij];
                        nrt[ij] -= rime_rate_nr[ij];

                        // BvS: I'm not sure about this part. If:
                        //      rime_rate_qs != -rime_rate_qr,
                        //      moisture is not conserved...?
                        qst[ij] -= rime_rate_qs[ij];
                        qrt[ij] -= rime_rate_qr[ij];

                        // Ice multiplication;
                        TF mult_q = TF(0);
                        TF mult_n = TF(0);

                        if (Ta[ij] < Constants::T0<TF> && ice_multiplication)
                        {
                            TF mult_1 = (Ta[ij] - T_mult_min<TF>) * const3;
                            TF mult_2 = (Ta[ij] - T_mult_max<TF>) * const4;

                            mult_1 = std::max(TF(0), std::min(mult_1, TF(1)));
                            mult_2 = std::max(TF(0), std::min(mult_2, TF(1)));

                            mult_n = C_mult<TF> * mult_1 * mult_2 * rime_rate_qr[ij];
                            mult_q = mult_n * ice.x_min;
                        }

                        if (Ta[ij] >= Constants::T0<TF>)
                        {
                            // Shedding of rain at warm temperatures;
                            // i.e. undo time integration, but with modified rain.n;
                            const TF x_r = particle_meanmass(rain, qr[ij], nr[ij]);

                            nst[ij] += rime_rate_nr[ij];
                            nrt[ij] += rime_rate_qr[ij] / x_r;

                            // BvS: I'm not sure about this part. If:
                            //      rime_rate_qs != -rime_rate_qr,
                            //      moisture is not conserved...?
                            qst[ij] += rime_rate_qs[ij];
                            qrt[ij] += rime_rate_qr[ij];
                        }
                        else
                        {
                            // new ice particles from multiplication;
                            nit[ij] += mult_n;
                            qit[ij] += mult_q;

                            // riming to graupel;
                            if (Ta[ij] < cfg_params.Tmax_gr_rime)
                            {
                                ngt[ij] += rime_rate_nr[ij];
                                qgt[ij] += rime_rate_qr[ij] + rime_rate_qs[ij] - mult_q;
                            }
                            else
                            {
                                // Snow + frozen liquid stays snow:;
                                nst[ij] += rime_rate_nr[ij];
                                qst[ij] += rime_rate_qr[ij] + rime_rate_qs[ij] - mult_q;
                            }
                        }
                    }
                } // Outer if
            } // i
    } // function


    template<typename TF>
    void particle_cloud_riming(
            TF* const restrict qct,
            TF* const restrict nct,
            TF* const restrict qpt,
            TF* const restrict npt,
            TF* const restrict qit,
            TF* const restrict nit,
            TF* const restrict qrt,
            TF* const restrict nrt,
            const TF* const restrict qc,
            const TF* const restrict nc,
            const TF* const restrict qp,
            const TF* const restrict np,
            const TF* const restrict Ta,
            Particle_frozen<TF>& ice,
            Particle_frozen<TF>& ptype,
            Particle<TF>& cloud,
            Particle<TF>& rain,
            Collection_coeffs<TF>& coeffs,
            const TF rho_v,
            const bool ice_multiplication,
            const bool enhanced_melting,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        /*
            Riming of graupel or hail with cloud droplets                                *
        */
        const TF const0 = TF(1)/(D_coll_c<TF> - D_crit_c<TF>);
        const TF const2 = TF(1)/(T_mult_opt<TF> - T_mult_min<TF>);
        const TF const3 = TF(1)/(T_mult_opt<TF> - T_mult_max<TF>);
        const TF const4 = clw<TF> / Constants::Lf<TF>;
        const TF const1 = const0 * ptype.ecoll_c;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                const TF x_p = particle_meanmass(ptype, qp[ij], np[ij]);
                const TF D_p = particle_diameter(ptype, x_p);
                const TF x_c = particle_meanmass(cloud, qc[ij], nc[ij]);
                const TF D_c = particle_diameter(cloud, x_c);

                if (qc[ij] > q_crit_c<TF> && qp[ij] > ptype.q_crit_c && D_p > ptype.D_crit_c && D_c > D_crit_c<TF>)
                {
                    const TF v_p = particle_velocity(ptype, x_p) * rho_v;
                    const TF v_c = particle_velocity(cloud, x_c) * rho_v;

                    const TF e_coll_n = std::min(ptype.ecoll_c, std::max(const1*(D_c - D_crit_c<TF>), ecoll_min<TF>));
                    const TF e_coll_q = e_coll_n;

                    // Both terms are multiplied by dt in ICON
                    const TF rime_n = pi4<TF> * e_coll_n * np[ij] * nc[ij]
                        *     (coeffs.delta_n_aa * D_p*D_p + coeffs.delta_n_ab * D_p*D_c + coeffs.delta_n_bb * D_c*D_c)
                        * sqrt(coeffs.theta_n_aa * v_p*v_p - coeffs.theta_n_ab * v_p*v_c + coeffs.theta_n_bb * v_c*v_c);

                    const TF rime_q = pi4<TF> * e_coll_q * np[ij] * qc[ij]
                        *     (coeffs.delta_q_aa * D_p*D_p + coeffs.delta_q_ab * D_p*D_c + coeffs.delta_q_bb * D_c*D_c)
                        * sqrt(coeffs.theta_q_aa * v_p*v_p - coeffs.theta_q_ab * v_p*v_c + coeffs.theta_q_bb * v_c*v_c);

                    qpt[ij] += rime_q;
                    qct[ij] -= rime_q;
                    nct[ij] -= rime_n;

                    // Ice multiplication based on Hallet and Mossop;
                    if (Ta[ij] < Constants::T0<TF> && ice_multiplication)
                    {
                        TF mult_1 = const2 * (Ta[ij] - T_mult_min<TF>);
                        TF mult_2 = const3 * (Ta[ij] - T_mult_max<TF>);

                        mult_1 = std::max(TF(0), std::min(mult_1, TF(1)));
                        mult_2 = std::max(TF(0), std::min(mult_2, TF(1)));

                        const TF mult_n = C_mult<TF> * mult_1 * mult_2 * rime_q;
                        const TF mult_q = mult_n * ice.x_min;

                        nit[ij] += mult_n;
                        qit[ij] += mult_q;
                        qpt[ij] -= mult_q;
                    }

                    // Enhancement of melting;
                    if (Ta[ij] > Constants::T0<TF> && enhanced_melting)
                    {
                        const TF melt_q = const4 * (Ta[ij] - Constants::T0<TF>) * rime_q;
                        const TF melt_n = melt_q / x_p;

                        qpt[ij] -= melt_q;
                        qrt[ij] += melt_q;
                        npt[ij] -= melt_n;
                        nrt[ij] += melt_n;
                    }
                }
            } // i
    } // function


    template<typename TF>
    void particle_rain_riming(
            TF* const restrict qpt,
            TF* const restrict npt,
            TF* const restrict qrt,
            TF* const restrict nrt,
            TF* const restrict qit,
            TF* const restrict nit,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF* const restrict qp,
            const TF* const restrict np,
            const TF* const restrict Ta,
            Particle<TF>& rain,
            Particle<TF>& ice,
            Particle<TF>& ptype,
            Collection_coeffs<TF>& coeffs,
            const TF rho_v,
            const bool ice_multiplication,
            const bool enhanced_melting,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        /*
            Riming of graupel or hail with rain droplets
        */

        const TF const2 = TF(1) / (T_mult_opt<TF> - T_mult_min<TF>);
        const TF const3 = TF(1) / (T_mult_opt<TF> - T_mult_max<TF>);
        const TF const4 = clw<TF> / Constants::Lf<TF>;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                if (qr[ij] > q_crit<TF> && qp[ij] > q_crit<TF>)
                {
                    const TF x_p = particle_meanmass(ptype, qp[ij], np[ij]);
                    const TF D_p = particle_diameter(ptype, x_p);
                    const TF v_p = particle_velocity(ptype, x_p) * rho_v;

                    const TF x_r = particle_meanmass(rain, qr[ij], nr[ij]);
                    const TF D_r = particle_diameter(rain, x_r);
                    const TF v_r = particle_velocity(rain, x_r) * rho_v;

                    // Both terms are `* dt` in ICON.
                    const TF rime_n = pi4<TF> * np[ij] * nr[ij]
                             * (coeffs.delta_n_aa * fm::pow2(D_p)
                                + coeffs.delta_n_ab * D_p*D_r + coeffs.delta_n_bb * fm::pow2(D_r))
                             * sqrt(coeffs.theta_n_aa * fm::pow2(v_p)
                                - coeffs.theta_n_ab * v_p*v_r + coeffs.theta_n_bb * fm::pow2(v_r));

                    const TF rime_q = pi4<TF> * np[ij] * qr[ij]
                             * (coeffs.delta_n_aa * fm::pow2(D_p)
                                + coeffs.delta_q_ab * D_p*D_r + coeffs.delta_q_bb * fm::pow2(D_r))
                             * sqrt(coeffs.theta_n_aa * fm::pow2(v_p)
                                - coeffs.theta_q_ab * v_p*v_r + coeffs.theta_q_bb * fm::pow2(v_r));

                    qpt[ij] += rime_q;
                    qrt[ij] -= rime_q;
                    nrt[ij] -= rime_n;

                    // Ice multiplication based on Hallet and Mossop;
                    if (Ta[ij] < Constants::T0<TF> && ice_multiplication)
                    {
                        TF mult_1 = (Ta[ij] - T_mult_min<TF>) * const2;
                        TF mult_2 = (Ta[ij] - T_mult_max<TF>) * const3;

                        mult_1 = std::max(TF(0), std::min(mult_1, TF(1)));
                        mult_2 = std::max(TF(0), std::min(mult_2, TF(1)));

                        const TF mult_n = C_mult<TF> * mult_1 * mult_2 * rime_q;
                        TF mult_q = mult_n * ice.x_min;
                        mult_q = std::min(rime_q, mult_q);

                        nit[ij] += mult_n;
                        qit[ij] += mult_q;
                        qpt[ij] -= mult_q;
                    }

                    // Enhancement of melting of ptype;
                    if (Ta[ij] > Constants::T0<TF> && enhanced_melting)
                    {
                        const TF melt_q = const4 * (Ta[ij] - Constants::T0<TF>) * rime_q;
                        const TF melt_n = melt_q / x_p;

                        qpt[ij] -= melt_q;
                        qrt[ij] += melt_q;

                        npt[ij] -= melt_n;
                        nrt[ij] += melt_n;
                    }
                }
            } // i
    } // function


    template<typename TF>
    inline TF incgfct_lower_lookup(
            const TF x,
            Gamlookuptable<TF>& ltable)
    {
        /*
            Retrieve values from a lookup table of the lower incomplete gamma function,
            as function of x at a constant a, for which the lookup table has been
            created.

            The last value in the table has to correspond to x = infinity, so that
            during the reconstruction of incgfct-values from the table,
            the x-value can safely be truncated at the maximum table x-value:

            ltable%igf( ltable%x(ltable%n),...) = gfct(a)

            Profiling with ifort on a Linux-PC shows, that table lookup for the
            incompl. gamma-Funktion is faster by a factor of about 15 compared
            to the original function without optimization (-O0). Using optimization
            could change this ratio (we encoutered up to 300 depending on function inlining).

            Concerning the accuracy, comparisons show that the results of table lookup
            are accurate to within better than 0.1 % or even much less, except for
            very small values of X, for which the absolute values are however very
            close to 0. For X -> infinity (X > 99.5 % - value), accuracy may be
            somewhat reduced up to about 0.5 % ,
            because the table is truncated at the 99.5 % value (second-last value)
            and the last value is set to the ordinary gamma function.

            This function only uses the low resolution part of the table!
        */

        // Trunkcate x to the range of the table:
        const TF xt = std::max( std::min(x, ltable.x[ltable.n-1]), TF(0));

        // Calculate indices of the neighbouring regular x-values in the table:
        const int i0 = int(xt * ltable.odx);
        const int iu = std::min(i0, ltable.n-2);
        const int io = iu + 1;

        // Interpolate linearly and subtract from the ordinary
        // gamma function to get the upper incomplete gamma function:
        const TF res = ltable.igf[iu] + (ltable.igf[io] - ltable.igf[iu]) * ltable.odx * (xt-ltable.x[iu]);

        return res;
    }


    template<typename TF>
    inline TF incgfct_upper_lookup(
            const TF x,
            Gamlookuptable<TF>& ltable)
    {
        // Trunkcate x to the range of the table:
        const TF xt = std::max( std::min(x, ltable.x[ltable.n-1]), TF(0));

        // Calculate indices of the neighbouring regular x-values in the table:
        const int i0 = int(xt * ltable.odx);
        const int iu = std::min(i0, ltable.n-2);
        const int io = iu + 1;

        // interpolate lower inc. gamma function linearily and subtract from
        // the ordinary gamma function to get the upper incomplete gamma function:
        const TF res = ltable.igf[ltable.n-1] - ltable.igf[iu] -
            (ltable.igf[io] - ltable.igf[iu]) * ltable.odx * (xt - ltable.x[iu]);

        // Aufgrund von Rundungsfehlern (Differenz von 2 fast gleichen Zahlen) kann es beim table lookup passieren,
        // dass incgfct_upper_lookup(x, ltable) kleiner 0 wird, wenn eigentlich nahezu 0.0 herauskommen muesste.
        // Dies kommt vor allem dann vor, wenn x sehr gross ist.
        // Deswegen Begrenzung:

        return std::max(res, TF(0));
    }


    template<typename TF>
    void rain_freeze_gamlook(
            TF* const restrict qit,
            TF* const restrict nit,
            TF* const restrict qrt,
            TF* const restrict nrt,
            TF* const restrict qgt,
            TF* const restrict ngt,
            TF* const restrict qht,
            TF* const restrict nht,
            const TF* const restrict qr,
            const TF* const restrict nr,
            const TF* const restrict Ta,
            Gamlookuptable<TF>& rain_ltable1,
            Gamlookuptable<TF>& rain_ltable2,
            Gamlookuptable<TF>& rain_ltable3,
            Particle_rain_coeffs<TF>& rain_coeffs,
            Particle<TF>& rain,
            const T_cfg_2mom<TF>& cfg_params,
            const TF rain_nm1,
            const TF rain_nm2,
            const TF rain_nm3,
            const TF rain_g1,
            const TF rain_g2,
            const TF dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        const TF a_HET = 6.5e-1;      // Data of Barklie and Gokhale (PK S.350)
        const TF b_HET = 2.0e+2;      //         Barklie and Gokhale (PK S.350)

        //const TF eps = 1e-15;         // for clipping
        //const bool lclipping = true;

        const TF xmax_ice = std::pow( std::pow(cfg_params.D_rainfrz_ig / rain.a_geo, TF(1) / rain.b_geo), rain.mu);
        const TF xmax_gr  = std::pow( std::pow(cfg_params.D_rainfrz_gh / rain.a_geo, TF(1) / rain.b_geo), rain.mu);

        const TF zdt = TF(1) / dt;

        TF fr_q;
        TF fr_n;
        TF fr_n_i;
        TF fr_q_i;
        TF fr_n_g;
        TF fr_q_g;
        TF fr_n_h;
        TF fr_q_h;
        TF fr_n_tmp;
        TF fr_q_tmp;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                // Make copy; n_r is potentially changed..
                TF q_r = qr[ij];
                TF n_r = nr[ij];

                if (Ta[ij] < T_freeze<TF>)
                {
                    if (q_r <= q_crit_fr<TF>)
                    {
                        if (Ta[ij] < T_f<TF>)
                        {
                            // Instantaneous freezing below T_f or -40 C
                            fr_q = q_r;
                            fr_n = n_r;
                            fr_n_i = n_r;
                            fr_q_i = q_r;
                            fr_n_g = TF(0);
                            fr_q_g = TF(0);
                            fr_n_h = TF(0);
                            fr_q_h = TF(0);
                            fr_n_tmp = TF(1);
                            fr_q_tmp = TF(1);
                        }
                        else
                        {
                            fr_q = TF(0);
                            fr_n = TF(0);
                            fr_n_i = TF(0);
                            fr_q_i = TF(0);
                            fr_n_g = TF(0);
                            fr_q_g = TF(0);
                            fr_n_h = TF(0);
                            fr_q_h = TF(0);
                            fr_n_tmp = TF(0);
                            fr_q_tmp = TF(0);
                        }
                    }
                    else
                    {
                        const TF x_r = particle_meanmass(rain, q_r, n_r);
                        n_r = q_r / x_r;

                        if (Ta[ij] < T_f<TF>)
                        {
                            // This branch could also be omitted. While it is quantitatively correct, it is not
                            // consistent with the limit case for complete freezing of the
                            // calculation in the T_a >= T_f branch below.
                            fr_q = q_r;                 //  Ausfrieren unterhalb T_f \approx -40 C;
                            fr_n = n_r;

                            // Depending on the size, the frozen raindrops are added to the cloud ice, or to the
                            // graupel or hail. For this purpose, a partial integration of the spectrum from 0;
                            // up to a first separation mass xmax_ice (--> ice), from there to xmax_gr (--> graupel);
                            // and from xmax_gr to infinity (--> hail).
                            const TF lam = std::exp( std::log( rain_g1/rain_g2 * x_r ) * (-rain.mu) );
                            const TF lam_rnm1 = std::exp(rain_nm1 * std::log(lam));  // lam**rain_nm1;
                            const TF lam_rnm2 = std::exp(rain_nm2 * std::log(lam));  // lam**rain_nm2;

                            const TF n_0 = rain.mu * n_r * lam_rnm1 / rain_g1;
                            fr_n_i = n_0 / (rain.mu * lam_rnm1) * incgfct_lower_lookup(lam * xmax_ice, rain_ltable1);
                            fr_q_i = n_0 / (rain.mu * lam_rnm2) * incgfct_lower_lookup(lam * xmax_ice, rain_ltable2);
                            fr_n_g = n_0 / (rain.mu * lam_rnm1) * incgfct_lower_lookup(lam * xmax_gr,  rain_ltable1);
                            fr_q_g = n_0 / (rain.mu * lam_rnm2) * incgfct_lower_lookup(lam * xmax_gr,  rain_ltable2);

                            fr_n_h = fr_n - fr_n_g;
                            fr_q_h = fr_q - fr_q_g;
                            fr_n_g = fr_n_g - fr_n_i;
                            fr_q_g = fr_q_g - fr_q_i;
                            fr_n_tmp = n_r / std::max(fr_n, n_r);
                            fr_q_tmp = q_r / std::max(fr_q, q_r);
                        }
                        else
                        {
                            //..Heterogeneous freezing;
                            const TF j_het = std::max(b_HET * (
                                std::exp( a_HET * (Constants::T0<TF> - Ta[ij]))- TF(1) ), TF(0)) / rho_w<TF> * dt;

                            //if (use_prog_in) j_het = std::min(j_het, n_inact[ij]/q_r);

                            // Depending on the size, the frozen raindrops are added to cloud ice; or to graupel or hail.
                            // This is achieved by partial integration of the spectrum from 0; up to a first separation
                            // mass xmax_ice (--> ice), from there up to xmax_gr (--> graupel);
                            // and from xmax_gr to infinity (--> hail).
                            if (j_het >= TF(1.0e-20))
                            {
                                fr_n  = j_het * q_r;
                                fr_q  = j_het * q_r * x_r * rain_coeffs.c_z;

                                // lam = ( rain_g1 / rain_g2 * x_r)**(-rain.mu);
                                const TF lam = std::exp( std::log( rain_g1/rain_g2*x_r ) * (-rain.mu) );
                                const TF lam_rnm1 = std::exp(rain_nm1*std::log(lam));  // lam**rain_nm1;
                                const TF lam_rnm2 = std::exp(rain_nm2*std::log(lam));  // lam**rain_nm2;
                                const TF lam_rnm3 = std::exp(rain_nm3*std::log(lam));  // lam**rain_nm3;

                                const TF n_0 = rain.mu * n_r * lam_rnm1 / rain_g1;
                                fr_n_i = j_het * n_0 / (rain.mu * lam_rnm2) *
                                        incgfct_lower_lookup(lam*xmax_ice,rain_ltable2);
                                fr_q_i = j_het * n_0 / (rain.mu * lam_rnm3) *
                                        incgfct_lower_lookup(lam*xmax_ice,rain_ltable3);
                                fr_n_g = j_het * n_0 / (rain.mu * lam_rnm2) *
                                        incgfct_lower_lookup(lam*xmax_gr, rain_ltable2);
                                fr_q_g = j_het * n_0 / (rain.mu * lam_rnm3) *
                                        incgfct_lower_lookup(lam*xmax_gr, rain_ltable3);

                                fr_n_h = fr_n - fr_n_g;
                                fr_q_h = fr_q - fr_q_g;
                                fr_n_g = fr_n_g - fr_n_i;
                                fr_q_g = fr_q_g - fr_q_i;
                                fr_n_tmp = n_r / std::max(fr_n, n_r);
                                fr_q_tmp = q_r / std::max(fr_q, q_r);
                            }
                            else
                            {
                                fr_n = TF(0);
                                fr_q = TF(0);
                                fr_n_i = TF(0);
                                fr_q_i = TF(0);
                                fr_n_g = TF(0);
                                fr_q_g = TF(0);
                                fr_n_h = TF(0);
                                fr_q_h = TF(0);
                                fr_n_tmp = TF(0);
                                fr_q_tmp = TF(0);
                            }
                        }

                        fr_n = fr_n * fr_n_tmp;
                        fr_q = fr_q * fr_q_tmp;
                        fr_n_i = fr_n_i * fr_n_tmp;
                        fr_n_g = fr_n_g * fr_n_tmp;
                        fr_n_h = fr_n_h * fr_n_tmp;
                        fr_q_i = fr_q_i * fr_q_tmp;
                        fr_q_g = fr_q_g * fr_q_tmp;
                        fr_q_h = fr_q_h * fr_q_tmp;
                    }

                    qrt[ij] -= fr_q * zdt;
                    nrt[ij] -= fr_n * zdt;

                    //if (use_prog_in) then;
                    //   n_inact[ij] = n_inact[ij] + fr_n;
                    //end if;

                    // mit Hagelklasse, gefrierender Regen wird Eis, Graupel oder Hagel;
                    //snow.q[ij] = snow.q[ij]  + fr_q_i;
                    //snow.n[ij] = snow.n[ij]  + fr_n_i   // put this into snow;

                    qit[ij] += fr_q_i * zdt; // ... or into ice? --> UB: original idea was to put it into ice;
                    nit[ij] += fr_n_i * zdt;

                    qgt[ij] += fr_q_g * zdt;
                    ngt[ij] += fr_n_g * zdt;

                    qht[ij] += fr_q_h * zdt;
                    nht[ij] += fr_n_h * zdt;

                    //! clipping of small negatives is necessary here
                    //if (lclipping) then
                    //    IF (rain%q(i,k) < 0.0 .and. abs(rain%q(i,k)) < eps) rain%q(i,k) = 0.0_wp
                    //IF (rain%n(i,k) < 0.0 .and. abs(rain%n(i,k)) < eps) rain%n(i,k) = 0.0_wp
                    //IF (graupel%q(i,k) < 0.0 .and. abs(graupel%q(i,k)) < eps) graupel%q(i,k) = 0.0_wp
                    //IF (graupel%n(i,k) < 0.0 .and. abs(graupel%q(i,k)) < eps) graupel%n(i,k) = 0.0_wp
                    //IF (hail%q(i,k) < 0.0 .and. abs(hail%q(i,k)) < eps) hail%q(i,k) = 0.0_wp
                    //IF (hail%n(i,k) < 0.0 .and. abs(hail%n(i,k)) < eps) hail%n(i,k) = 0.0_wp
                    //end if
                }
            } // i
    }


    template<typename TF>
    TF dmin_wg_gr_ltab_equi(
            const TF p_a,
            const TF T_a,
            const TF qw_a,
            const TF qi_a,
            Lookupt_4D<TF>& ltab)
    {
        // Quadrilinear interpolation of lookup table.
        // LUT is equidistant to simplify indexing.

        if (T_a >= ltab.x2[ltab.n2-1])
            return TF(0);
        else if (T_a < ltab.x2[0])
            return TF(999.99);
        else
        {
            // Limit on lookup table bounds.
            const TF p_lok  = std::min( std::max(p_a,  ltab.x1[0]), ltab.x1[ltab.n1-1]);
            const TF T_lok  = std::min( std::max(T_a,  ltab.x2[0]), ltab.x2[ltab.n2-1]);
            const TF qw_lok = std::min( std::max(qw_a, ltab.x3[0]), ltab.x3[ltab.n3-1]);
            const TF qi_lok = std::min( std::max(qi_a, ltab.x4[0]), ltab.x4[ltab.n4-1]);

            // Index left of interpolation value.
            const int i0 = std::min( int((p_lok  - ltab.x1[0]) * ltab.odx1), ltab.n1-2);
            const int j0 = std::min( int((T_lok  - ltab.x2[0]) * ltab.odx2), ltab.n2-2);
            const int k0 = std::min( int((qw_lok - ltab.x3[0]) * ltab.odx3), ltab.n3-2);
            const int l0 = std::min( int((qi_lok - ltab.x4[0]) * ltab.odx4), ltab.n4-2);

            // Index right of interpolation value.
            const int i1 = i0 + 1;
            const int j1 = j0 + 1;
            const int k1 = k0 + 1;
            const int l1 = l0 + 1;

            // Interpolation factors.
            const TF fi0 = TF(1) - ((p_lok  - ltab.x1[i0]) / (ltab.x1[i1] - ltab.x1[i0]));
            const TF fj0 = TF(1) - ((T_lok  - ltab.x2[j0]) / (ltab.x2[j1] - ltab.x2[j0]));
            const TF fk0 = TF(1) - ((qw_lok - ltab.x3[k0]) / (ltab.x3[k1] - ltab.x3[k0]));
            const TF fl0 = TF(1) - ((qi_lok - ltab.x4[l0]) / (ltab.x4[l1] - ltab.x4[l0]));

            const TF fi1 = TF(1) - fi0;
            const TF fj1 = TF(1) - fj0;
            const TF fk1 = TF(1) - fk0;
            const TF fl1 = TF(1) - fl0;

            const int ls = ltab.n4_stride;
            const int ks = ltab.n3_stride;
            const int js = ltab.n2_stride;

            auto index = [](
                    const int l,  const int k,  const int j,  const int i,
                    const int ls, const int ks, const int js)
            {
                return i + j*js + k*ks + l*ls;
            };

            // Interpolate.
            const TF val =
                    fl0 * fk0 * fj0 * fi0 * ltab.ltable[ index(l0, k0, j0, i0, ls, ks, js) ] +
                    fl0 * fk0 * fj0 * fi1 * ltab.ltable[ index(l0, k0, j0, i1, ls, ks, js) ] +
                    fl0 * fk0 * fj1 * fi0 * ltab.ltable[ index(l0, k0, j1, i0, ls, ks, js) ] +
                    fl0 * fk0 * fj1 * fi1 * ltab.ltable[ index(l0, k0, j1, i1, ls, ks, js) ] +
                    fl0 * fk1 * fj0 * fi0 * ltab.ltable[ index(l0, k1, j0, i0, ls, ks, js) ] +
                    fl0 * fk1 * fj0 * fi1 * ltab.ltable[ index(l0, k1, j0, i1, ls, ks, js) ] +
                    fl0 * fk1 * fj1 * fi0 * ltab.ltable[ index(l0, k1, j1, i0, ls, ks, js) ] +
                    fl0 * fk1 * fj1 * fi1 * ltab.ltable[ index(l0, k1, j1, i1, ls, ks, js) ] +
                    fl1 * fk0 * fj0 * fi0 * ltab.ltable[ index(l1, k0, j0, i0, ls, ks, js) ] +
                    fl1 * fk0 * fj0 * fi1 * ltab.ltable[ index(l1, k0, j0, i1, ls, ks, js) ] +
                    fl1 * fk0 * fj1 * fi0 * ltab.ltable[ index(l1, k0, j1, i0, ls, ks, js) ] +
                    fl1 * fk0 * fj1 * fi1 * ltab.ltable[ index(l1, k0, j1, i1, ls, ks, js) ] +
                    fl1 * fk1 * fj0 * fi0 * ltab.ltable[ index(l1, k1, j0, i0, ls, ks, js) ] +
                    fl1 * fk1 * fj0 * fi1 * ltab.ltable[ index(l1, k1, j0, i1, ls, ks, js) ] +
                    fl1 * fk1 * fj1 * fi0 * ltab.ltable[ index(l1, k1, j1, i0, ls, ks, js) ] +
                    fl1 * fk1 * fj1 * fi1 * ltab.ltable[ index(l1, k1, j1, i1, ls, ks, js) ];

            return val;
        }
    }


    template<typename TF>
    void graupel_hail_conv_wet_gamlook(
            TF* const restrict qgt,
            TF* const restrict ngt,
            TF* const restrict qht,
            TF* const restrict nht,
            const TF* const restrict qc,
            const TF* const restrict qr,
            const TF* const restrict qg,
            const TF* const restrict ng,
            const TF* const restrict qi,
            const TF* const restrict qs,
            const TF* const restrict Ta,
            Gamlookuptable<TF>& graupel_ltable1,
            Gamlookuptable<TF>& graupel_ltable2,
            Particle_frozen<TF>& graupel,
            Lookupt_4D<TF>& ltabdminwgg,
            const TF graupel_nm1,
            const TF graupel_nm2,
            const TF graupel_g1,
            const TF graupel_g2,
            const TF pa,
            const TF dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        /*
            Wet growth and conversion of graupel to hail
            (uses look-up table for incomplete gamma functions)
            by Uli Blahak
        */
        const TF zdt = TF(1) / dt;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                const TF x_g = particle_meanmass(graupel, qg[ij], ng[ij]);
                const TF D_g = particle_diameter(graupel, x_g);

                const TF q_g = qg[ij];
                const TF n_g = q_g / x_g; // for consistency for limiters, n_g is used explicitly below;

                // Supercooled liquid water in the cloud environment = sum of rain and cloud water;
                const TF qw_a = qr[ij] + qc[ij];

                if (Ta[ij] < Constants::T0<TF> && q_g > graupel.q_crit_c && qw_a > 1e-3)
                {
                    // Umgebungsgehalt Eispartikel (vernachl. werden Graupel und Hagel wg. geringer Kollisionseff.);
                    // koennte problematisch sein, weil in konvekt. Wolken viel mehr Graupel und Hagel enthalten ist;
                    const TF qi_a = qi[ij] + qs[ij];

                    //if (luse_dmin_wetgrowth_table)
                    const TF d_trenn = dmin_wg_gr_ltab_equi(pa, Ta[ij], qw_a, qi_a, ltabdminwgg);
                    //else
                    //    d_trenn = dmin_wetgrowth_fun(p_a,T_a,qw_a,qi_a);
                    //end if;

                    if (d_trenn > TF(0) && d_trenn < TF(10) * D_g)
                    {
                        TF xmin = std::exp( std::log( d_trenn/graupel.a_geo ) * (TF(1) / graupel.b_geo) );
                        const TF lam  = std::exp( std::log( graupel_g2 / (graupel_g1 * x_g)) * (graupel.mu) );
                        xmin = std::exp( std::log(xmin) * graupel.mu );
                        const TF n_0  = graupel.mu * n_g * std::exp( std::log(lam) * graupel_nm1 ) / graupel_g1;

                        const TF conv_n = n_0 / (graupel.mu * std::exp( std::log(lam)*graupel_nm1))
                                * incgfct_upper_lookup(lam*xmin, graupel_ltable1);
                        const TF conv_q = n_0 / (graupel.mu * std::exp( std::log(lam)*graupel_nm2))
                                * incgfct_upper_lookup(lam*xmin, graupel_ltable2);

                        // BvS: this implies that conv_X is already time integrated here???
                        //conv_n = std::min(conv_n, n_g);
                        //conv_q = std::min(conv_q, q_g);

                        qgt[ij] -= conv_q * zdt;
                        ngt[ij] -= conv_n * zdt;

                        qht[ij] += conv_q * zdt;
                        nht[ij] += conv_n * zdt;
                    }
                }
            }
    }


    template<typename TF>
    void ice_melting(
            TF* const restrict qct,
            TF* const restrict nct,
            TF* const restrict qit,
            TF* const restrict nit,
            TF* const restrict qrt,
            TF* const restrict nrt,
            const TF* const restrict qi,
            const TF* const restrict ni,
            const TF* const restrict Ta,
            Particle<TF>& ice,
            Particle<TF>& cloud,
            const TF dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        const TF zdt = TF(1)/dt;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                if (Ta[ij] > Constants::T0<TF> && qi[ij] > TF(0))
                {
                    const TF x_i = particle_meanmass(ice, qi[ij], ni[ij]);

                    // Complete melting within this time step;
                    const TF melt_q = qi[ij] * zdt;
                    const TF melt_n = ni[ij] * zdt;

                    qit[ij] -= melt_q;
                    nit[ij] -= melt_n;

                    // Ice either melts into cloud droplets or rain depending on x_i;
                    if (x_i > cloud.x_max)
                    {
                        qrt[ij] += melt_q;
                        nrt[ij] += melt_n;
                    }
                    else
                    {
                        qct[ij] += melt_q;
                        nct[ij] += melt_n;
                    }
                }
            } // i
    } // function


    template<typename TF>
    void snow_melting(
            TF* const restrict qst,
            TF* const restrict nst,
            TF* const restrict qrt,
            TF* const restrict nrt,
            const TF* const restrict qs,
            const TF* const restrict ns,
            const TF* const restrict Ta,
            Particle_sphere<TF>& snow_coeffs,
            Particle<TF>& snow,
            const TF rho_v,
            const TF dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        const TF zdt = TF(1)/dt;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                if (Ta[ij] > Constants::T0<TF> && qs[ij] > TF(0))
                {
                    const TF e_a = tmf::esat_liq(Ta[ij]);

                    const TF x_s = particle_meanmass(snow, qs[ij], ns[ij]);
                    const TF D_s = particle_diameter(snow, x_s);
                    const TF v_s = particle_velocity(snow, x_s) * rho_v;

                    const TF fv_q = snow_coeffs.a_f + snow_coeffs.b_f * sqrt(v_s * D_s);

                    // UB: Based on Rasmussen and Heymsfield (1987) the ratio fh_q / fv_q
                    // is approx. 1.05 for a wide temperature- and pressure range;
                    const TF fh_q = TF(1.05) * fv_q;

                    const TF melt = TF(2) * pi<TF>/Constants::Lf<TF> * D_s * ns[ij] * dt;

                    // NOTE: ICON uses R_d here, which is `gas constant of water vapor (dampf)`, so Rv in our case...
                    const TF melt_h = melt * K_T<TF> * (Ta[ij] - Constants::T0<TF>);
                    const TF melt_v = melt * D_v<TF> * Constants::Lv<TF> /
                            Constants::Rv<TF> * (e_a/Ta[ij] - e_3<TF>/Constants::T0<TF>);
                    TF melt_q = (melt_h * fh_q + melt_v * fv_q);

                    // UB: for melt_n we assume that x_s is constant during melting;
                    TF melt_n = std::min(std::max( (melt_q - qs[ij]) / x_s + ns[ij], TF(0)), ns[ij]);

                    melt_q = std::min(qs[ij], std::max(melt_q, TF(0)));
                    melt_n = std::min(ns[ij], std::max(melt_n, TF(0)));

                    // UB: snow melts instantaneously at 10 C;
                    if (Ta[ij] - Constants::T0<TF> > TF(10))
                    {
                        melt_q = qs[ij];
                        melt_n = ns[ij];
                    }

                    qst[ij] -= melt_q * zdt;
                    qrt[ij] += melt_q * zdt;

                    nst[ij] -= melt_n * zdt;
                    nrt[ij] += melt_n * zdt;

                    // Oiiiii
                    //snow.n[ij] = std::max(snow.n[ij], snow.q[ij]/snow.x_max);
                }
            } // i
    } // function


    template<typename TF>
    void graupel_melting(
            TF* const restrict qgt,
            TF* const restrict ngt,
            TF* const restrict qrt,
            TF* const restrict nrt,
            const TF* const restrict qg,
            const TF* const restrict ng,
            const TF* const restrict Ta,
            Particle_sphere<TF>& graupel_coeffs,
            Particle<TF>& graupel,
            const TF rho_v,
            const TF dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        const TF zdt = TF(1)/dt;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                if (Ta[ij] > Constants::T0<TF> && qg[ij] > TF(0))
                {
                    const TF e_a = tmf::esat_liq(Ta[ij]);

                    const TF x_g = particle_meanmass(graupel, qg[ij], ng[ij]);
                    const TF D_g = particle_diameter(graupel, x_g);
                    const TF v_g = particle_velocity(graupel, x_g) * rho_v;

                    const TF fv_q = graupel_coeffs.a_f + graupel_coeffs.b_f * sqrt(v_g * D_g);
                    const TF fh_q = TF(1.05) * fv_q;

                    const TF melt   = TF(2) * pi<TF> / Constants::Lf<TF> * D_g * ng[ij] * dt;

                    const TF melt_h = melt * K_T<TF> * (Ta[ij] - Constants::T0<TF>);
                    const TF melt_v = melt * D_v<TF> * Constants::Lv<TF> /
                                      Constants::Rv<TF> * (e_a/Ta[ij] - e_3<TF>/Constants::T0<TF>);
                    TF melt_q = (melt_h * fh_q + melt_v * fv_q);

                    // UB: assume that x_g is constant during melting;
                    TF melt_n = std::min(std::max( (melt_q - qg[ij]) / x_g + ng[ij], TF(0)), ng[ij]);

                    melt_q = std::min(qg[ij], std::max(melt_q, TF(0)));
                    melt_n = std::min(ng[ij], std::max(melt_n, TF(0)));

                    qgt[ij] -= melt_q * zdt;
                    qrt[ij] += melt_q * zdt;

                    ngt[ij] -= melt_n * zdt;
                    nrt[ij] += melt_n * zdt;
                }
            } // i
    } // function


    template<typename TF>
    void hail_melting_simple(
            TF* const restrict qht,
            TF* const restrict nht,
            TF* const restrict qrt,
            TF* const restrict nrt,
            const TF* const restrict qh,
            const TF* const restrict nh,
            const TF* const restrict Ta,
            Particle_sphere<TF>& hail_coeffs,
            Particle<TF>& hail,
            const T_cfg_2mom<TF>& cfg_params,
            const TF rho_v,
            const TF dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        const TF zdt = TF(1)/dt;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                if (Ta[ij] > Constants::T0<TF> && qh[ij] > TF(0))
                {
                    const TF e_a = tmf::esat_liq(Ta[ij]);

                    const TF x_h = particle_meanmass(hail, qh[ij], nh[ij]);
                    const TF D_h = particle_diameter(hail, x_h);
                    const TF v_h = particle_velocity(hail, x_h) * rho_v;

                    const TF fv_q = hail_coeffs.a_f + hail_coeffs.b_f * sqrt(v_h * D_h);
                    const TF fh_q = TF(1.05) * fv_q;    // UB: based on Rasmussen and Heymsfield

                    const TF melt = TF(2) * pi<TF> / Constants::Lf<TF> * D_h * nh[ij] * dt;

                    const TF melt_h = melt * K_T<TF> * (Ta[ij] - Constants::T0<TF>);
                    const TF melt_v = melt * D_v<TF> * Constants::Lv<TF> /
                                      Constants::Rv<TF> * (e_a/Ta[ij] - e_3<TF>/Constants::T0<TF>);
                    TF melt_q = cfg_params.melt_h_tune_fak * (melt_h * fh_q + melt_v * fv_q);

                    // UB: assume that x_h is constant during melting;
                    TF melt_n = std::min(std::max( (melt_q - qh[ij]) / x_h + nh[ij], TF(0)), nh[ij]);

                    melt_q = std::min(qh[ij], std::max(melt_q, TF(0)));
                    melt_n = std::min(nh[ij], std::max(melt_n, TF(0)));

                    qht[ij] -= melt_q * zdt;
                    qrt[ij] += melt_q * zdt;

                    nht[ij] -= melt_n * zdt;
                    nrt[ij] += melt_n * zdt;
                }
            } // i
    } // function


    template<typename TF>
    void ice_nucleation_het_philips(
            TF* const restrict qit,
            TF* const restrict nit,
            TF* const restrict qvt,
            TF* const restrict n_inact,
            const TF* const restrict qv,
            const TF* const restrict ql,
            const TF* const restrict Ta,
            const TF* const restrict afrac_dust,
            const TF* const restrict afrac_soot,
            const TF* const restrict afrac_orga,
            Particle<TF>& ice,
            bool use_prog_in,
            const TF na_dust,
            const TF na_soot,
            const TF na_orga,
            const TF dt,
            const int afrac_stride,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        const TF eps = TF(1.0e-20);     // DANGEROUS for SP.
        const TF zdt = TF(1) / dt;

        // (i,j) strides in lookup table.
        const int iis = afrac_stride;
        const int jjs = 1;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                const TF e_si = tmf::esat_ice(Ta[ij]);    // e_es(T_a) in ICON = sat_pres_ice
                const TF ssi  = qv[ij] * Constants::Rv<TF> * Ta[ij] / e_si;

                TF infrac1, infrac2, infrac3; // Interpolated values from LUT

                if (Ta[ij] < T_nuc<TF> && Ta[ij] > TF(180) && ssi > TF(1) && n_inact[ij] < ni_het_max<TF>)
                {
                    // Interpolation index x-direction
                    TF xt = (TF(274.) - Ta[ij]) / ttstep;
                    xt = std::min(xt, TF(ttmax-1));
                    const int ii = int(xt);

                    // Interpolation factors x-direction
                    const TF f0x = TF(ii) + TF(1) - xt;
                    const TF f1x = xt - TF(ii);

                    if (ql[ij] > eps)
                    {
                        // Immersion freezing at water saturation.
                        const int jj = 99;
                        const int ij  = jj + ii*iis;

                        infrac1 = f0x * afrac_dust[ij]
                                + f1x * afrac_dust[ij+iis];

                        infrac2 = f0x * afrac_soot[ij]
                                + f1x * afrac_soot[ij+iis];

                        infrac3 = f0x * afrac_orga[ij]
                                + f1x * afrac_orga[ij+iis];
                    }
                    else
                    {
                        // Deposition nucleation below water saturation;
                        // Calculate indices used for 2D look-up tables;
                        TF xs = 100. * (ssi - TF(1.0)) / ssstep;
                        xs = std::min(xs, TF(ssmax-1));
                        const int jj = std::max(1, int(xs));
                        const TF ssr = std::max(TF(1), std::floor(xs));

                        // Interpolation factor y-direction.
                        const TF f0y = ssr + TF(1) - xs;
                        const TF f1y = xs - ssr;

                        const int ij  = jj + ii*iis;

                        infrac1 = f0x * f0y * afrac_dust[ij]
                                + f1x * f0y * afrac_dust[ij+iis]
                                + f0x * f1y * afrac_dust[ij+jjs]
                                + f1x * f1y * afrac_dust[ij+iis+jjs];

                        infrac2 = f0x * f0y * afrac_soot[ij]
                                + f1x * f0y * afrac_soot[ij+iis]
                                + f0x * f1y * afrac_soot[ij+jjs]
                                + f1x * f1y * afrac_soot[ij+iis+jjs];

                        infrac3 = f0x * f0y * afrac_orga[ij]
                                + f1x * f0y * afrac_orga[ij+iis]
                                + f0x * f1y * afrac_orga[ij+jjs]
                                + f1x * f1y * afrac_orga[ij+iis+jjs];
                    }

                    // Sum up the three modes;
                    TF ndiag;
                    if (use_prog_in)
                    {
                        throw std::runtime_error("\"use_prog_in\" not (yet) supported..");
                        // `n_inpot` replaces `na_dust`; `na_soot` and `na_orga` are assumed to be constant;
                        //ndiag  = n_inpot[ij] * infrac1 + na_soot * infrac2 + na_orga * infrac3;
                        //ndiag_dust = n_inpot[ij] * infrac1;
                        //ndiag_all = ndiag;
                    }
                    else
                    {
                        // All aerosol species are diagnostic;
                        ndiag = na_dust * infrac1 + na_soot * infrac2 + na_orga * infrac3;
                    }

                    ndiag = std::min(ndiag, ni_het_max<TF>);
                    TF nuc_n = std::max(ndiag - n_inact[ij],TF(0));
                    TF nuc_q = std::min(nuc_n * ice.x_min, qv[ij]);
                    nuc_n = nuc_q / ice.x_min;

                    // From absolute change -> tendency.
                    // NOTE: there is no `*dt` in ICON, but `nuc_q` and `nuc_n` are assumed to
                    // be total increments of `ni` and `qi`....?
                    nuc_q *= zdt;
                    nuc_n *= zdt;

                    // Store tendencies.
                    qit[ij] += nuc_q;
                    nit[ij] += nuc_n;
                    qvt[ij] -= nuc_q;
                    n_inact[ij] += nuc_n;

                    //lwrite_n_inpot = use_prog_in && ndiag .GT. 1.0e-12_wp;
                    //ndiag_mask(i, k) = lwrite_n_inpot;
                    //if (lwrite_n_inpot)
                    //{
                    //nuc_n = nuc_n * ndiag_dust / ndiag_all;
                    //}
                    //nuc_n_a(i, k) = nuc_n;
                    //}
                    //else
                    //{
                    //nuc_n_a(i, k) = TF(0);
                    //ndiag_mask(i, k) = .FALSE.;
                    //}
                } // if
            } // loop
    } // function


    template<typename TF>
    void ice_nucleation_homhet(
            TF* const restrict qvt,
            TF* const restrict qit,
            TF* const restrict nit,
            TF* const restrict n_inact,
            //TF* const restrict n_inpot,
            const TF* const restrict qi,
            const TF* const restrict ni,
            const TF* const restrict qv,
            const TF* const restrict ql,
            const TF* const restrict Ta,
            const TF* const restrict w,
            const TF* const restrict afrac_dust,
            const TF* const restrict afrac_soot,
            const TF* const restrict afrac_orga,
            Particle<TF>& ice,
            bool use_prog_in,
            const int nuc_i_typ,
            const TF p,
            const TF dt,
            const int afrac_stride,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        /*
            Homogeneous and heterogeneous ice nucleation.

            Nucleation scheme is based on the papers:
            - "A parametrization of cirrus cloud formation: Homogenous
              freezing of supercooled aerosols" by B. Kaercher and
              U. Lohmann 2002 (KL02 hereafter)
            - "Physically based parameterization of cirrus cloud formation
              for use in global atmospheric models" by B. Kaercher, J. Hendricks
              and U. Lohmann 2006 (KHL06 hereafter)
            - Phillips et al. (2008) with extensions
        */

        const TF zdt = TF(1) / dt;

        // Switch for version of Phillips et al. scheme.
        // Only `2010` is supported (also in ICON)..
        const int iphillips = 2010;

        // Other switches, set below based on certain criteria.
        bool use_homnuc;

        // Variables set dynamically below
        TF na_dust, na_soot, na_orga;

        // Some more constants needed for homogeneous nucleation scheme.
        const TF r_0     = 0.25e-6;                     // aerosol particle radius prior to freezing
        const TF alpha_d = 0.5;                         // deposition coefficient (KL02; Spichtinger & Gierens 2009)
        const TF M_w     = 18.01528e-3;                 // molecular mass of water [kg/mol]
        const TF M_a     = 28.96e-3;                    // molecular mass of air [kg/mol]
        const TF ma_w    = M_w / N_avo<TF>;             // mass of water molecule [kg]
        const TF svol    = ma_w / Constants::rho_i<TF>; // specific volume of a water molecule in ice

        // Parameters for Hande et al. nucleation parameterization for HDCP2 simulations
        //TYPE(dep_imm_coeffs), PARAMETER :: hdcp2_nuc_coeffs(5) &
        //     = (/ dep_imm_coeffs(&  ! Spring of Table 1
        //     nin_imm = 1.5684e5, &
        //     alf_imm = 0.2466,   &
        //     bet_imm = 1.2293,   &
        //     nin_dep = 1.7836e5, &
        //     alf_dep = 0.0075,   &
        //     bet_dep = 2.0341),  &
        //     dep_imm_coeffs(     &  ! Summer
        //     nin_imm = 2.9694e4, &
        //     alf_imm = 0.2813,   &
        //     bet_imm = 1.1778,   &
        //     nin_dep = 2.6543e4, &
        //     alf_dep = 0.0020,   &
        //     bet_dep = 2.5128),  &
        //     dep_imm_coeffs(     &  ! Autumn
        //     nin_imm = 4.9920e4, &
        //     alf_imm = 0.2622,   &
        //     bet_imm = 1.2044,   &
        //     nin_dep = 7.7167e4, &
        //     alf_dep = 0.0406,   &
        //     bet_dep = 1.4705),  &
        //     dep_imm_coeffs(     &  ! Winter
        //     nin_imm = 1.0259e5, &
        //     alf_imm = 0.2073,   &
        //     bet_imm = 1.2873,   &
        //     nin_dep = 1.1663e4, &
        //     alf_dep = 0.0194,   &
        //     bet_dep = 1.6943),  &
        //     dep_imm_coeffs(     &  ! Spring with 95th percentile scaling factor
        //     nin_imm = 1.5684e5 * 17.82, &
        //     alf_imm = 0.2466,   &
        //     bet_imm = 1.2293,   &
        //     nin_dep = 1.7836e5 * 5.87, &
        //     alf_dep = 0.0075,   &
        //     bet_dep = 2.0341) /)

        if (nuc_i_typ == 0)
            use_homnuc = false;
        else if (nuc_i_typ < 10)
            use_homnuc = true;
        else
            throw std::runtime_error("Invalid \"nuc_i_typ\".");

        // switch for Hande et al. ice nucleation, if `true`, this turns off Phillips scheme.
        const bool use_hdcp2_het = (nuc_i_typ <= 5);
        if (use_hdcp2_het)
            throw std::runtime_error("Hande et al. ice nucleation not (yet) implemented.");

        // Heterogeneous nucleation using Hande et al. scheme;
        if (use_hdcp2_het)
        {
            //IF (nuc_typ < 1 .OR. nuc_typ > 5) THEN
            //    CALL finish(TRIM(routine), &
            //        & 'Error in two_moment_mcrph: Invalid value nuc_typ in case of &
            //        &use_hdcp2_het=.true.')
            //END IF
            //CALL ice_nucleation_het_hdcp2(ik_slice, atmo, ice, cloud, &
            //    use_prog_in, hdcp2_nuc_coeffs(nuc_typ), n_inact, ndiag_mask, nuc_n_a)
        }
        else
        {
            // Heterogeneous nucleation using Phillips et al. scheme;
            if (iphillips == 2010)
            {
                // Possible pre-defined choices.
                if (nuc_i_typ == 6)
                {
                    // With no organics and rather high soot, coming close to Meyers formula at -20 C
                    na_dust  = TF(160.e4);  // initial number density of dust [1/m3];
                    na_soot  = TF(30.e6);   // initial number density of soot [1/m3];
                    na_orga  = TF(0.e0);    // initial number density of organics [1/m3];
                }
                else if (nuc_i_typ == 7)
                {
                    // With some organics and rather high soot.
                    na_dust  = TF(160.e4);  // Coming close to Meyers formula at -20 C;
                    na_soot  = TF( 25.e6);
                    na_orga  = TF( 30.e6);
                }
                else if (nuc_i_typ == 8)
                {
                    // No organics, no soot, coming close to DeMott et al. 2010 at -20 C
                    na_dust  = TF(70.e4);  // i.e. roughly one order in magnitude lower than Meyers;
                    na_soot  = TF(0.e6);
                    na_orga  = TF(0.e6);
                }
                else
                    throw std::runtime_error("Invalid value nuc_typ in case of use_hdcp2_het=false");
            }

            ice_nucleation_het_philips(
                qit, nit, qvt,
                n_inact,
                qv, ql, Ta,
                afrac_dust, afrac_soot, afrac_orga,
                ice, use_prog_in,
                na_dust, na_soot, na_orga,
                dt, afrac_stride,
                istart, iend,
                jstart, jend,
                jstride);
        }

        //IF (use_prog_in) THEN
        //  DO k = kstart, kend
        //    DO i = istart, iend
        //      n_inpot(i,k) = MERGE(MAX(n_inpot(i,k) - nuc_n_a(i, k), 0.0_wp), &
        //           &               n_inpot(i, k), &
        //           &               ndiag_mask(i, k))
        //    END DO
        //  END DO
        //END IF

        if (use_homnuc)
        {
            // Homogeneous nucleation using KHL06 approach
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ij = i + j*jstride;

                    const TF e_si = tmf::esat_ice(Ta[ij]);    // e_es(T_a) in ICON = sat_pres_ice
                    const TF ssi  = qv[ij] * Constants::Rv<TF> * Ta[ij] / e_si;

                    // Critical supersaturation for homogeneous nucleation;
                    const TF scr = TF(2.349) - Ta[ij] * (TF(1) / TF(259));

                    if (ssi > scr && Ta[ij] < TF(235) && ni[ij] < ni_hom_max<TF>)
                    {
                        const TF x_i = particle_meanmass(ice, qi[ij], ni[ij]);

                        //r_i = (x_i/(4./3.*pi*rho_ice))**(1./3.);
                        const TF r_i =
                            std::exp( (TF(1)/TF(3)) * std::log(x_i / (TF(4)/TF(3) * pi<TF> * Constants::rho_i<TF>)) );

                        const TF v_th = sqrt(TF(8) * k_b<TF> * Ta[ij] / (pi<TF> * ma_w));
                        const TF flux = alpha_d * v_th / TF(4);
                        const TF n_sat = e_si / (k_b<TF> * Ta[ij]);

                        // Coeffs of supersaturation equation;
                        const TF acoeff1 =
                            (Constants::Ls<TF> * Constants::grav<TF>) /
                            (Constants::cp<TF> * Constants::Rv<TF> * fm::pow2(Ta[ij])) - Constants::grav<TF> /
                                    (Constants::Rd<TF> * Ta[ij]);
                        const TF acoeff2 = TF(1) / n_sat;
                        const TF acoeff3 =
                            (fm::pow2(Constants::Ls<TF>) * M_w * ma_w) / (Constants::cp<TF> * p * Ta[ij] * M_a);

                        // Coeffs of depositional growth equation;
                        const TF bcoeff1 = flux * svol * n_sat * (ssi - TF(1));
                        const TF bcoeff2 = flux / diffusivity(Ta[ij], p);

                        // Pre-existing ice crystals included as reduced updraft speed;
                        const TF ri_dot = bcoeff1 / (TF(1) + bcoeff2 * r_i);
                        const TF R_ik = (TF(4) * pi<TF>) / svol * ni[ij] * fm::pow2(r_i) * ri_dot;
                        TF w_pre = (acoeff2 + acoeff3 * ssi)/(acoeff1 * ssi) * R_ik;  // KHL06 Eq. 19;
                        w_pre = std::max(w_pre, TF(0));

                        // NOTE BvS: using `w[ij]` is a bit cheap (staggered grid); this assumes
                        // that nucleation only happens if there is a positive vertical velocity
                        // at the bottom of the current grid box...
                        // Would interpolation of `w` to full levels be better?
                        if (w[ij] > w_pre)
                        {
                            // Homogenous nucleation event
                            // Timescales of freezing event (see KL02, RM05, KHL06);
                            const TF cool = Constants::grav<TF> / Constants::cp<TF> * w[ij];
                            const TF ctau = Ta[ij] * (TF(0.004) * Ta[ij] - TF(2)) + TF(304.4);
                            const TF tau = TF(1) / (ctau * cool);   // freezing timescale, eq. (5);
                            const TF delta = (bcoeff2 * r_0);       // dimless aerosol radius, eq.(4);
                            const TF phi = acoeff1 * ssi / ( acoeff2 + acoeff3 * ssi) * (w[ij] - w_pre);

                            // Monodisperse approximation following KHL06;
                            const TF kappa = TF(2) * bcoeff1 * bcoeff2 * tau / fm::pow2(TF(1) + delta);  // kappa, Eq. 8 KHL06;
                            const TF sqrtkap = sqrt(kappa);
                            // Analytical approx. of erfc by RM05:
                            const TF ren = TF(3) * sqrtkap / (TF(2) + sqrt(TF(1) + TF(9) * kappa / pi<TF>));
                            const TF R_imfc  = TF(4) * pi<TF> * bcoeff1/fm::pow2(bcoeff2) / svol;
                            const TF R_im    = R_imfc / (TF(1) + delta) * ( fm::pow2(delta) - TF(1)
                                + (TF(1) + TF(0.5) * kappa * fm::pow2(TF(1) + delta)) * ren/sqrtkap);    // RIM Eq. 6 KHL06

                            // Number concentration and radius of ice particles;
                            const TF ni_hom  = phi / R_im;                                  // ni Eq.9 KHL06;
                            const TF ri_0 = TF(1) + TF(0.5) * sqrtkap * ren;                // for Eq. 3 KHL06;
                            const TF ri_hom  = (ri_0 * (TF(1) + delta) - TF(1)) / bcoeff2;  // Eq. 3 KHL06 * REN = Eq.23 KHL06;
                            TF mi_hom  = (TF(4)/TF(3) * pi<TF> * Constants::rho_i<TF>) * ni_hom * fm::pow3(ri_hom);
                            mi_hom  = std::max(mi_hom, ice.x_min);

                            TF nuc_n = std::max(std::min(ni_hom, ni_hom_max<TF>), TF(0));
                            //const TF nuc_q = std::min(nuc_n * mi_hom, atmo.qv[ij]);
                            TF nuc_q = nuc_n * mi_hom;

                            // From absolute change -> tendency.
                            // NOTE: there is no `*dt` in ICON, but `nuc_q` and `nuc_n` are assumed to
                            // be total increments of `ni` and `qi`....?
                            nuc_n *= zdt;
                            nuc_q *= zdt;

                            qit[ij] += nuc_q;
                            nit[ij] += nuc_n;
                            qvt[ij] -= nuc_q;
                        }
                    }
                } // loop
        } // if
    } // function

    template<typename TF>
    void cloud_freeze(
            TF* const restrict qct,
            TF* const restrict nct,
            TF* const restrict qit,
            TF* const restrict nit,
            const TF* const restrict qc,
            const TF* const restrict nc,
            const TF* const restrict ni,
            const TF* const restrict Ta,
            const int nuc_c_typ,
            const TF dt,
            Particle<TF>& cloud,
            Particle_cloud_coeffs<TF>& cloud_coeffs,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        const TF log_10 = std::log(10);
        const TF rhow_i = TF(1)/Constants::rho_w<TF>;
        const TF dti = TF(1)/dt;

        TF fr_q, fr_n, j_hom;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j*jstride;

                if (Ta[ij] < Constants::T0<TF>)
                {
                    const TF T_c = Ta[ij] - Constants::T0<TF>;

                    if (qc[ij] > TF(0) && T_c < TF(-30))
                    {
                        if (T_c < TF(-50))
                        {
                            fr_q = qc[ij];  //..instantaneous freezing
                            fr_n = nc[ij];  //..below -50 C
                        }
                        else
                        {
                            const TF x_c = particle_meanmass(cloud, qc[ij], nc[ij]);

                            //..Hom. freezing based on Jeffrey und Austin (1997), see also Cotton und Field (2001)
                            //  (note that log in Cotton and Field is log10, not ln)
                            if (T_c > TF(-30.0))
                            {
                                // j_hom = 1.0e6_wp/rho_w * 10**(-7.63-2.996*(T_c+30.0))           !..J in 1/(kg s)
                                j_hom = TF(1.0e6) * rhow_i * std::exp((TF(-7.63) - TF(2.996) * (T_c+TF(30)))*log_10);
                            }
                            else
                            {
                                //j_hom = 1.0e6_wp/rho_w &
                                //     &  * 10**(-243.4-14.75*T_c-0.307*T_c**2-0.00287*T_c**3-0.0000102*T_c**4)
                                j_hom = TF(1.0e6) * rhow_i * std::exp((TF(-243.4) - TF(14.75) * T_c - TF(0.307) *
                                    fm::pow2(T_c) - TF(0.00287) * fm::pow3(T_c) - TF(0.0000102) * fm::pow4(T_c))*log_10);
                            }

                            fr_n  = j_hom * qc[ij] * dt;
                            fr_q  = j_hom * qc[ij] * dt * x_c * cloud_coeffs.c_z;
                            fr_q  = std::min(fr_q, qc[ij]);
                            fr_n  = std::min(fr_n, nc[ij]);
                        }

                        qct[ij] -= fr_q * dti;
                        nct[ij] -= fr_n * dti;

                        fr_n = std::max(fr_n, fr_q/cloud.x_max);
                        // Special treatment for constant drop number; Force upper bound in cloud_freeze
                        if (nuc_c_typ == 0)
                            fr_n = std::max(std::min(fr_n, nc[ij]-ni[ij]), TF(0));

                        qit[ij] += fr_q * dti;
                        nit[ij] += fr_n * dti;
                    }
                }
            }
    }
} // namespace
