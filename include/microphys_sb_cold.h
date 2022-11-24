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

#include "constants.h"
#include "fast_math.h"

namespace Sb_cold
{
    /*
       Kernels ported from ICON.
       All kernels use the same name as the ICON subroutines
    */

    using namespace Constants;
    namespace fm = Fast_math;

    template<typename TF> constexpr TF pi = 3.14159265359;
    template<typename TF> constexpr TF rho_w = 1.e3;                   // Density water
    template<typename TF> constexpr TF rho_0 = 1.225;                  // SB06, p48
    template<typename TF> constexpr TF pirhow = pi<TF>*rho_w<TF>/6.;
    template<typename TF> constexpr TF rho_vel = 0.4;                  // Exponent for density correction (value from ICON)
    template<typename TF> constexpr TF kc_autocon = 9.44e+9;           // Long-Kernel

    // Limiters on ql/qr/etc.
    template<typename TF> constexpr TF q_crit = 1.e-9;                 // Min rain liquid water for which calculations are performed (value from ICON)


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
                //if (qr[ij] > qr_min<TF> * rho[k]) // TODO: remove *rho..
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


    template<typename TF>
    void autoconversionSB(
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
        // Autoconversion of Seifert and Beheng (2001, Atmos. Res.)
        // Formation of rain by coagulating cloud droplets

        const TF k_1  = 6.00e+2;   //..Parameter for Phi
        const TF k_2  = 0.68e+0;   //..Parameter fof Phi
        const TF eps  = 1.00e-25;  // NOTE BvS: dangerous for float version MicroHH.
        const TF x_s_i = TF(1) / cloud.x_max;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j * jstride;
                const int ijk = i + j * jstride + k * kstride;

                if (qc[ijk] > q_crit<TF>)
                {
                    const TF n_c = Nc0; // cloud%n(i,k) in ICON
                    const TF x_c = particle_meanmass(cloud, qc[ijk], n_c);

                    TF au = cloud_coeffs.k_au * fm::pow2(qc[ijk]) * fm::pow2(x_c) * cloud_rho_v;    // NOTE `*dt` in ICON..
                    const TF tau = std::min(std::max(TF(1) - qc[ijk] / (qc[ijk] + qr[ij] + eps), eps), TF(0.9));
                    const TF phi = k_1 * std::pow(tau, k_2) * fm::pow3(TF(1) - std::pow(tau, k_2));
                    au = au * (TF(1) + phi / fm::pow2(TF(1) - tau));

                    nrt[ij] += au * x_s_i;
                    qrt[ij] += au;

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
            TF* const restrict qrt,
            const TF* const restrict qr,
            const TF* const restrict qc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride, const int kstride,
            const int k)
    {
        // Accretion: growth of raindrops collecting cloud droplets
        // Accretion of Seifert and Beheng (2001, Atmos. Res.)

        const TF k_r = 5.78;       // kernel
        const TF k_1 = 5.00e-04;   // Phi function
        const TF eps = 1.00e-25;

        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij = i + j * jstride;
                const int ijk = i + j * jstride + k * kstride;

                if (qc[ijk] > TF(0) && qr[ij] > TF(0))
                {

                    // ..accretion rate of SB2001
                    const TF tau = std::min(std::max(TF(1) - qc[ijk] / (qc[ijk] + qr[ij] + eps), eps), TF(1));
                    const TF phi = fm::pow4(tau/(tau+k_1));
                    const TF ac  = k_r *  qc[ijk] * qr[ij] * phi;  // NOTE: `*dt` in ICON..

                    qrt[ij] += ac;

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
        // Selfcollection & breakup: growth of raindrops by mutual (rain-rain) coagulation, and breakup by collisions
        // Selfcollection of Seifert and Beheng (2001, Atmos. Res.)                     *

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
}
