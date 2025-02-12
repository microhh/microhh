/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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

#ifndef LAND_SURFACE_KERNELS_GPU_H
#define LAND_SURFACE_KERNELS_GPU_H

#include "tools.h"
#include "constants.h"
#include "boundary_surface_kernels_gpu.h"
#include "thermo_moist_functions.h"
#include "monin_obukhov.h"

namespace Land_surface_kernels_g
{
    namespace bsk = Boundary_surface_kernels_g;
    namespace most = Monin_obukhov;

    template<typename TF> __global__
    void calc_tile_fractions_g(
            TF* const __restrict__ tile_frac_veg,
            TF* const __restrict__ tile_frac_soil,
            TF* const __restrict__ tile_frac_wet,
            const TF* const __restrict__ wl,
            const TF* const __restrict__ c_veg,
            const TF* const __restrict__ lai,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*icells;
            const TF wlm = Constants::wlmax<TF> * (TF(1)-c_veg[ij] + c_veg[ij]*lai[ij]);

            tile_frac_wet[ij]  = fmin(TF(1), wl[ij]/wlm);
            tile_frac_veg[ij]  = (TF(1)-tile_frac_wet[ij]) * c_veg[ij];
            tile_frac_soil[ij] = (TF(1)-tile_frac_wet[ij]) * (TF(1)-c_veg[ij]);
        }
    }

    template<typename TF> __global__
    void calc_liquid_water_reservoir_g(
            TF* const __restrict__ wl_tend,
            TF* const __restrict__ interception,
            TF* const __restrict__ throughfall,
            const TF* const __restrict__ wl,
            const TF* const __restrict__ LE_veg,
            const TF* const __restrict__ LE_soil,
            const TF* const __restrict__ LE_wet,
            const TF* const __restrict__ tile_frac_veg,
            const TF* const __restrict__ tile_frac_soil,
            const TF* const __restrict__ tile_frac_wet,
            const TF* const __restrict__ rain_rate,
            const TF* const __restrict__ c_veg,
            const TF* const __restrict__ lai,
            const double subdt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        const TF intercept_eff = TF(0.5);
        const TF to_ms  = TF(1) / (Constants::rho_w<TF> * Constants::Lv<TF>);
        const TF subdti = TF(1) / subdt;

        if (i < iend && j < jend)
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
            const TF wl_tend_liq = -fmax(TF(0), tile_frac_wet[ij] * LE_wet[ij] * to_ms);

            // Tendency due to dewfall into vegetation/soil/liquid water tiles
            const TF wl_tend_dew = -( fmin(TF(0), tile_frac_wet[ij]  * LE_wet[ij]  * to_ms)
                                    + fmin(TF(0), tile_frac_veg[ij]  * LE_veg[ij]  * to_ms)
                                    + fmin(TF(0), tile_frac_soil[ij] * LE_soil[ij] * to_ms) );

            // Tendency due to interception of precipitation by vegetation
            // Rain rate is positive downwards, so minus is excluded.
            const TF wl_tend_precip = intercept_eff * c_veg[ij] * rr_ms;

            // Total and limited tendencies
            const TF wl_tend_sum = wl_tend_liq + wl_tend_dew + wl_tend_precip;
            const TF wl_tend_lim = fmin(wl_tend_max, fmax(wl_tend_min,  wl_tend_sum));

            // Diagnose throughfall and interception
            throughfall[ij] =
                -(TF(1)-c_veg[ij]) * rr_ms
                -(TF(1)-intercept_eff) * c_veg[ij] * rr_ms +
                fmin(TF(0), wl_tend_lim - wl_tend_sum);

            interception[ij] = fmax(TF(0), wl_tend_lim);

            wl_tend[ij] += wl_tend_lim;
        }
    }

    template<typename TF> __global__
    void calc_resistance_functions_g(
            TF* const __restrict__ f1,
            TF* const __restrict__ f2,
            TF* const __restrict__ f3,
            const TF* const __restrict__ sw_dn,
            const TF* const __restrict__ theta_mean_n,
            const TF* const __restrict__ vpd,
            const TF* const __restrict__ gD,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells, const int ijcells)
    {
        // Constants f1 calculation:
        const TF a_f1 = 0.81;
        const TF b_f1 = 0.004;
        const TF c_f1 = 0.05;

        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*icells;

            // f1: reduction vegetation resistance as f(sw_in):
            const TF sw_dn_lim = fmax(TF(0), sw_dn[ij]);
            f1[ij] = TF(1)/fmin( TF(1), (b_f1*sw_dn_lim + c_f1) / (a_f1 * (b_f1*sw_dn_lim + TF(1))) );

            // f2: reduction vegetation resistance as f(theta):
            f2[ij] = TF(1)/fmin( TF(1), fmax(TF(1e-9), theta_mean_n[ij]) );

            // f3: reduction vegetation resistance as f(VPD):
            f3[ij] = TF(1)/exp(-gD[ij] * vpd[ij]);
        }
    }

    template<typename TF> __global__
    void calc_soil_resistance_function_g(
            TF* const __restrict__ f2b,
            const TF* const __restrict__ theta,
            const TF* const __restrict__ c_veg,
            const TF* const __restrict__ theta_wp,
            const TF* const __restrict__ theta_fc,
            const TF* const __restrict__ theta_res,
            const int* const __restrict__ soil_index,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kend,
            const int icells, const int ijcells)
    {
        // Constants f1 calculation:
        const TF a_f1 = 0.81;
        const TF b_f1 = 0.004;
        const TF c_f1 = 0.05;

        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*icells;
            const int ijk = i + j*icells + (kend-1)*ijcells;    // Top soil layer
            const int si  = soil_index[ijk];

            // f2b: reduction soil resistance as f(theta)
            const TF theta_min = c_veg[ij] * theta_wp[si] + (TF(1)-c_veg[ij]) * theta_res[si];
            const TF theta_rel = (theta[ijk] - theta_min) / (theta_fc[si] - theta_min);
            f2b[ij] = TF(1)/fmin(TF(1), fmax(TF(1e-9), theta_rel));
        }
    }


    template<typename TF> __global__
    void calc_canopy_resistance_g(
            TF* const __restrict__ rs,
            const TF* const __restrict__ rs_min,
            const TF* const __restrict__ lai,
            const TF* const __restrict__ f1,
            const TF* const __restrict__ f2,
            const TF* const __restrict__ f3,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*icells;
            rs[ij] = rs_min[ij] / (lai[ij]+Constants::dsmall) * f1[ij] * f2[ij] * f3[ij];
        }
    }


    template<typename TF> __device__
    inline TF e1_approx_g(const TF x)
    {
        // E1() approximation, from: https://doi.org/10.1016/S0022-1694(99)00184-5
        const TF euler = TF(0.5772156649015329);
        const TF G = exp(-euler);
        const TF b = pow(TF(2) * ((TF(1) - G) / (G * (TF(2) - G))), TF(0.5));
        const TF h_inf = (TF(1) - G) * (fm::pow2(G) - TF(6) * G + TF(12)) / (TF(3) * G * fm::pow2(TF(2) - G) * b);
        const TF q = TF(20/47.) * pow(x, pow(TF(31/26.), TF(0.5)));
        const TF h = TF(1) / (TF(1) + x * pow(x, TF(0.5)));
        const TF E1 = exp(-x) / (G + (TF(1) - G) * exp(-x / (TF(1) - G))) * log(TF(1) + G / x - (TF(1) - G) / fm::pow2(h + b * x));
        return E1;
    }


    template<typename TF> __global__
    void calc_canopy_resistance_ags_g(
            TF* const restrict rs,
            TF* const restrict an_co2,
            const TF* const restrict lai,
            const TF* const restrict t_bot,
            const TF* const restrict ra,
            const TF* const restrict co2,
            const TF* const restrict thl,
            const TF* const restrict qt,
            const TF* const restrict sw_flux_dn,
            const TF* const restrict theta_rel_mean,
            //const TF* const restrict albedo,
            const TF* const restrict vpds,
            const TF* const restrict alpha0,
            const TF* const restrict t1gm,
            const TF* const restrict t2gm,
            const TF* const restrict t1am,
            const TF* const restrict gm298,
            const TF* const restrict gmin,
            const TF* const restrict ammax298,
            const TF* const restrict f0,
            const TF* const restrict co2_comp298,
            //const TF cos_sza,
            const TF rho,
            const TF rho_bot,
            const bool sw_splitleaf,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int jstride,
            const int kstride)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*jstride;
            const int ijk = ij + kstart*kstride;

            // Fixed constants.
            const TF q10gm = TF(2);      // Parameter to calculate the mesophyll conductance
            const TF q10am = TF(2);      // Parameter to calculate max primary productivity
            const TF q10lambda = TF(1);  // Parameter to calculate the CO2 compensation concentration. (2 in IFS, 1.5 in DALES)

            const TF t2am = TF(311);     // Reference temperatures calculation max primary productivity [K]

            const TF nuco2q = TF(1.6);   // Ratio molecular viscosity water to carbon dioxide [-]
            const TF ad = TF(0.07);      // Regression coefficient to calculate Cfrac [kPa-1]
            const TF kx = TF(0.7);       // Extinction coefficient PAR in canopy [m_ground m_leaf-1]

            // For sw_splitleaf.
            //nr_gauss = 3                                   # Amount of bins to use for Gaussian integrations
            //weight_g = np.array([0.2778, 0.4444, 0.2778])  # Weights of the Gaussian bins (must add up to 1)
            //angle_g  = np.array([0.1127,    0.5, 0.8873])  # Sines of the leaf angles compared to the sun in the first Gaussian integration
            //angle_weight_sum = np.sum(weight_g * angle_g)
            //LAI_g    = np.array([0.1127,    0.5, 0.8873])  # Ratio of integrated LAI at locations where shaded leaves are evaluated in the second Gaussian integration
            //sigma    = 0.2                                 # Scattering coefficient
            //kdfbl    = 0.8                                 # Diffuse radiation extinction coefficient for black leaves

            // Fixed conversion factors.
            const TF to_mgm3 = Constants::xmco2<TF> / Constants::xmair<TF> * TF(1e6) * rho;
            const TF to_mol  = Constants::xmair<TF> / Constants::xmco2<TF> * TF(1e-6) * rho_bot;

            // Calculate the CO2 compensation concentration (IFS eq. 8.92).
            // "The compensation point Γ is defined as the CO2 concentration at which the net CO2 assimilation of a fully lit leaf becomes zero".
            const TF co2_comp = rho_bot * co2_comp298[ij] * pow(q10lambda, TF(0.1) * (t_bot[ij] - TF(298)));

            // Calculate the mesophyll conductance (IFS eq. 8.93).
            // "The mesophyll conductance gm describes the transport of CO2 from the substomatal cavities to the mesophyll cells where the carbon is fixed".
            const TF gm = gm298[ij] * pow(q10gm, TF(0.1) * (t_bot[ij] - TF(298))) / ((TF(1) + exp(TF(0.3) * (t1gm[ij] - t_bot[ij]))) * (TF(1) + exp(TF(0.3) * (t_bot[ij] - t2gm[ij])))) / TF(1000);

            // Calculate CO2 concentration inside the leaf (ci).
            // NOTE: Differs from IFS
            const TF fmin0 = gmin[ij] / nuco2q - TF(1./9.) * gm;
            const TF fmin = pow(max(-fmin0 + (fm::pow2(fmin0) + TF(4) * gmin[ij] / nuco2q * gm), TF(Constants::dsmall)), TF(0.5)) / (TF(2) * gm);

            // Calculate atmospheric moisture deficit.
            // "Therefore Ci/Cs is specified as a function of atmospheric moisture deficit Ds at the leaf surface".
            // Based on other literature, I think this indeed has to be the "leaf to air VPD", so esat(Ts) - e(Ta).
            const TF ds = vpds[ij] / 1000.;

            // This seems to differ from IFS?
            const TF dmax = (f0[ij] - fmin) / ad;

            // Coupling factor (IFS eq. 8.101).
            const TF cfrac = max(TF(0.01), f0[ij] * (TF(1) - ds / dmax) + fmin * (ds / dmax));

            // Absolute CO2 concentration. Input = [mol(CO2) mol-1(air)], co2_abs = [mg(CO2) m-3(air)].
            const TF co2_abs  = co2[ijk] * to_mgm3;

            // CO2 concentration in leaf (IFS eq. ???).
            //if (lrelaxci) then
            //  if (ci_old_set) then
            //    ci_inf        = cfrac * (co2_abs - co2_comp) + co2_comp
            //    ci            = ci_old(i,j) + min(kci*rk3coef, 1.0) * (ci_inf - ci_old(i,j))
            //    if (rk3step  == 3) then
            //      ci_old(i,j) = ci
            //    endif
            //  else
            //    ci            = cfrac * (co2_abs - co2_comp) + co2_comp
            //    ci_old(i,j)   = ci
            //  endif
            //else
            const TF ci = cfrac * (co2_abs - co2_comp) + co2_comp;
            //endif

            // Max gross primary production in high light conditions (Ag) (IFS eq. 8.94).
            const TF ammax = ammax298[ij] * pow(q10am, TF(0.1) * (t_bot[ij] - TF(298))) / ((TF(1) + exp(TF(0.3) * (t1am[ij] - t_bot[ij]))) * (TF(1) + exp(TF(0.3) * (t_bot[ij] - t2am))));

            // Effect of soil moisture stress on gross assimilation rate.
            // NOTE: this seems to be different in IFS...
            const TF fstr = max(TF(1e-3), min(TF(1), theta_rel_mean[ij]));

            // Gross assimilation rate (Am, IFS eq. 8.97).
            const TF am = ammax * (TF(1) - exp(-(gm * (ci - co2_comp) / ammax)));

            // Autotrophic dark respiration (IFS eq. 8.99).
            const TF rdark = am / TF(9);

            // Photosynthetically active radiation [W m-2].
            const TF par = TF(0.5) * max(TF(0.1), sw_flux_dn[ij]);

            // Light use efficiency.
            const TF alphac = alpha0[ij] * (co2_abs - co2_comp) / (co2_abs + 2 * co2_comp);

            TF an, gc_inf;
            //if (sw_splitleaf)
            //{
            //    PARdir   = 0.5 * max(0.1, sw_in[j,i])
            //    PARdif   = 0.5 * max(0.1, sw_in[j,i])
            //    cos_sza  = max(1e-10, cos_sza)

            //    kdrbl = 0.5 / cos_sza   # Direct radiation extinction coefficient for black leaves
            //    kdf = kdfbl * np.sqrt(1.0 - sigma)
            //    kdr = kdrbl * np.sqrt(1.0 - sigma)
            //    ref = (1.0 - np.sqrt(1.0 - sigma)) / (1.0 + np.sqrt(1.0 - sigma)) # Reflection coefficient
            //    ref_dir = 2 * ref / (1.0 + 1.6 * cos_sza)

            //    Hleaf = np.zeros(nr_gauss+1)
            //    Fleaf = np.zeros(nr_gauss+1)
            //    Agl   = np.zeros(nr_gauss+1)

            //    Fnet  = np.zeros(nr_gauss)
            //    gnet  = np.zeros(nr_gauss)

            //    # Loop over different LAI locations
            //    for it in range(nr_gauss):

            //        iLAI = LAI[j,i] * LAI_g[it]   # Integrated LAI between here and canopy top; Gaussian distributed
            //        fSL = np.exp(-kdrbl * iLAI)  # Fraction of sun-lit leaves

            //        PARdfD = PARdif * (1.0-ref)     * np.exp(-kdf * iLAI)             # Total downward PAR due to diffuse radiation at canopy top
            //        PARdrD = PARdir * (1.0-ref_dir) * np.exp(-kdr * iLAI)             # Total downward PAR due to direct radiation at canopy top
            //        PARdfU = PARdif * (1.0-ref)     * np.exp(-kdf * LAI[j,i]) * \
            //            albedo[j,i] * (1.0-ref)     * np.exp(-kdf * (LAI[j,i]-iLAI))  # Total upward (reflected) PAR that originates as diffuse radiation
            //        PARdrU = PARdir * (1.0-ref_dir) * np.exp(-kdr * LAI[j,i]) * \
            //            albedo[j,i] * (1.0-ref)     * np.exp(-kdf * (LAI[j,i]-iLAI))  # Total upward (reflected) PAR that originates as direct radiation
            //        PARdfT = PARdfD + PARdfU                                          # Total PAR due to diffuse radiation at canopy top
            //        PARdrT = PARdrD + PARdrU                                          # Total PAR due to direct radiation at canopy top

            //        dirPAR = (1.0-sigma) * PARdir * fSL   # Purely direct PAR (can only be downward)
            //        difPAR = PARdfT + PARdrT - dirPAR     # Total diffuse radiation

            //        HdfT = kdf * PARdfD + kdf * PARdfU
            //        HdrT = kdr * PARdrD + kdf * PARdrU
            //        dirH = kdrbl * dirPAR
            //        Hshad = HdfT + HdrT - dirH

            //        Hsun = Hshad + angle_g * (1.0-sigma) * kdrbl * PARdir / angle_weight_sum

            //        Hleaf[0] = Hshad
            //        Hleaf[1:] = Hsun

            //        Agl = fstr * (Am + Rdark) * (1 - np.exp(-alphac*Hleaf/(Am + Rdark)))
            //        gleaf = gmin[ij]/nuco2q +  Agl/(co2_abs-ci)
            //        Fleaf = Agl - Rdark

            //        Fshad  = Fleaf[0]
            //        Fsun   = np.sum(weight_g * Fleaf[1:])
            //        gshad  = gleaf[0]
            //        gsun   = np.sum(weight_g * gleaf[1:])

            //        Fnet[it] = Fsun * fSL + Fshad * (1 - fSL)
            //        gnet[it] = gsun * fSL + gshad * (1 - fSL)

            //    An     = LAI[j,i] * np.sum(weight_g * Fnet)
            //    gc_inf = LAI[j,i] * np.sum(weight_g * gnet)
            //}
            //else
            //{
                // Calculate upscaling from leaf to canopy: net flow CO2 into the plant (An)
                const TF agsa1  = TF(1) / (TF(1) - f0[ij]);
                const TF dstar  = dmax / (agsa1 * (f0[ij] - fmin));
                const TF tempy  = alphac * kx * par / (am + rdark);

                an     = (am + rdark) * (TF(1) - TF(1) / (kx * lai[ij]) * (e1_approx_g(tempy * exp(-kx * lai[ij])) - e1_approx_g(tempy)));
                gc_inf = lai[ij] * (gmin[ij] / nuco2q + agsa1 * fstr * an / ((co2_abs - co2_comp) * (TF(1) + ds / dstar)));
            //}

            //if (lrelaxgc) then
            //  if (gc_old_set) then
            //    gcco2       = gc_old(i,j) + min(kgc*rk3coef, 1.0) * (gc_inf - gc_old(i,j))
            //    if (rk3step ==3) then
            //      gc_old(i,j) = gcco2
            //    endif
            //  else
            //    gcco2 = gc_inf
            //    gc_old(i,j) = gcco2
            //  endif
            //else
            //    gcco2 = gc_inf
            const TF gcco2 = gc_inf;

            // Output values from kernel:
            // Canopy resistances for moisture.
            rs[ij] = TF(1) / (TF(1.6) * gcco2);

            // Net flux of CO2 into the plant (An).
            // `co2_abs` has units [mg(CO2) m-3(air)], so flux has units [mg(CO2) m-2 s-1)].
            // Conversion with `to_mol` results in flux in [kmol(CO2) kmol-1(air) m s-1].
            an_co2[ij] = (ci - co2_abs) / (ra[ij] + (TF(1) / gcco2)) * to_mol;


            // Check!
            if (isnan(rs[ij]))
            {
                printf("Canopy resistance contains NaNs!\n");
                asm("trap;");
            }
        }
    }


    template<typename TF> __global__
    void calc_soil_respiration_jacobs_g(
            TF* const restrict resp_co2,
            const TF* const restrict t_soil_mean,
            const TF* const restrict r10,
            const TF* const restrict ea,
            const TF rho_bot,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        const TF to_mol  = Constants::xmair<TF> / Constants::xmco2<TF> * TF(1e-6) * rho_bot;

        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*jstride;
            resp_co2[ij] = r10[ij] * exp(ea[ij] / (TF(283.15) * TF(8.314)) * (TF(1) - TF(283.15) / t_soil_mean[ij])) * to_mol;
        }
    }


    template<typename TF> __global__
    void set_net_co2_flux_g(
            TF* const restrict co2_flux,
            const TF* const restrict co2_an,
            const TF* const restrict co2_resp,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jstride)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*jstride;
            co2_flux[ij] = co2_an[ij] + co2_resp[ij];
        }
    }


    template<typename TF> __global__
    void calc_soil_resistance_g(
            TF* const __restrict__ rs,
            const TF* const __restrict__ rs_min,
            const TF* const __restrict__ f2b,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*icells;
            rs[ij] = rs_min[ij] * f2b[ij];
        }
    }

    template<typename TF, bool sw_constant_z0> __global__
    void calc_stability_g(
            TF* const __restrict__ ustar,
            TF* const __restrict__ obuk,
            TF* const __restrict__ bfluxbot,
            TF* const __restrict__ ra,
            int* const __restrict__ nobuk,
            const TF* const __restrict__ dutot,
            const TF* const __restrict__ b,
            const TF* const __restrict__ bbot,
            const TF* const __restrict__ z0m,
            const TF* const __restrict__ z0h,
            const float* const __restrict__ zL_sl,
            const float* const __restrict__ f_sl,
            const TF db_ref,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int jcells,
            const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*icells;
            const int ijk = i + j*icells + kstart*ijcells;

            const TF db = b[ijk] - bbot[ij] + db_ref;

            if (sw_constant_z0)
                obuk[ij] = bsk::calc_obuk_noslip_dirichlet_lookup_g(
                        zL_sl, f_sl, nobuk[ij], dutot[ij], db, zsl);
            else
                obuk[ij] = bsk::calc_obuk_noslip_dirichlet_iterative_g(
                        obuk[ij], dutot[ij], db, zsl, z0m[ij], z0h[ij]);

            ustar[ij] = dutot[ij] * most::fm(zsl, z0m[ij], obuk[ij]);
            bfluxbot[ij] = -ustar[ij] * db * most::fh(zsl, z0h[ij], obuk[ij]);
            ra[ij]  = TF(1) / (ustar[ij] * most::fh(zsl, z0h[ij], obuk[ij]));
        }
    }

    template<typename TF> __global__
    void calc_fluxes_g(
            TF* const __restrict__ H,
            TF* const __restrict__ LE,
            TF* const __restrict__ G,
            TF* const __restrict__ S,
            TF* const __restrict__ thl_bot,
            TF* const __restrict__ qt_bot,
            const TF* const __restrict__ T,
            const TF* const __restrict__ qt,
            const TF* const __restrict__ T_soil,
            const TF* const __restrict__ qsat_bot,
            const TF* const __restrict__ dqsatdT_bot,
            const TF* const __restrict__ ra,
            const TF* const __restrict__ rs,
            const TF* const __restrict__ lambda_stable,
            const TF* const __restrict__ lambda_unstable,
            const TF* const __restrict__ cs_veg,
            const TF* const __restrict__ sw_dn,
            const TF* const __restrict__ sw_up,
            const TF* const __restrict__ lw_dn,
            const TF* const __restrict__ lw_up,
            const TF* const __restrict__ b,
            const TF* const __restrict__ b_bot,
            const TF rho_bot,
            const TF exner_bot,
            const TF db_ref,
            const TF emis_sfc,
            const TF dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend_soil,
            const int icells, const int ijcells,
            bool use_cs_veg)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
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
            const TF fH  = rho_bot * Constants::cp<TF> / ra[ij];
            const TF fLE = rho_bot * Constants::Lv<TF> / (ra[ij] + rs_lim);

            // Net radiation; negative sign = net input of energy at surface
            const TF Qnet = sw_dn[ij] - sw_up[ij] + lw_dn[ij] - lw_up[ij];

            // Solve for the new surface temperature
            const TF num =
                Qnet + lw_up[ij] + fH*T[ij] +
                fLE*(qt[ijk] + dqsatdT_bot[ij]*T_bot - qsat_bot[ij]) +
                lambda*T_soil[ijk_s] + TF(3)*emis_sfc*Constants::sigma_b<TF>*fm::pow4(T_bot) - (TF(1)-emis_sfc) * lw_dn[ij];
            const TF denom = fH + fLE*dqsatdT_bot[ij] + lambda + TF(4)*emis_sfc*Constants::sigma_b<TF>*fm::pow3(T_bot);
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
            qt_bot[ij]  = qt[ijk] + LE[ij] * ra[ij] / (rho_bot * Constants::Lv<TF>);
        }
    }

    template<typename TF> __global__
    void calc_bulk_obuk_g(
            TF* const __restrict__ obuk,
            const TF* const __restrict__ bfluxbot,
            const TF* const __restrict__ ustar,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*icells;
            obuk[ij] = -fm::pow3(ustar[ij]) / (Constants::kappa<TF> * bfluxbot[ij]);
            obuk[ij] = zsl/fmin(fmax(zsl/obuk[ij], Constants::zL_min<TF>), Constants::zL_max<TF>);
        }
    }

    template<typename TF> __global__
    void set_bcs_momentum_g(
            TF* const restrict ufluxbot,
            TF* const restrict vfluxbot,
            TF* const restrict ugradbot,
            TF* const restrict vgradbot,
            const TF* const restrict ustar,
            const TF* const restrict u,
            const TF* const restrict ubot,
            const TF* const restrict v,
            const TF* const restrict vbot,
            const TF* const restrict z0m,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int jcells,
            const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        const int ii = 1;
        const int jj = icells;

        const TF minval = 1.e-2;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*ijcells;

            const TF vonu2 = fmax(minval, TF(0.25)*(
                        fm::pow2(v[ijk-ii]-vbot[ij-ii]) + fm::pow2(v[ijk-ii+jj]-vbot[ij-ii+jj])
                      + fm::pow2(v[ijk   ]-vbot[ij   ]) + fm::pow2(v[ijk   +jj]-vbot[ij   +jj])) );
            const TF uonv2 = fmax(minval, TF(0.25)*(
                        fm::pow2(u[ijk-jj]-ubot[ij-jj]) + fm::pow2(u[ijk+ii-jj]-ubot[ij+ii-jj])
                      + fm::pow2(u[ijk   ]-ubot[ij   ]) + fm::pow2(u[ijk+ii   ]-ubot[ij+ii   ])) );

            const TF u2 = fmax(minval, fm::pow2(u[ijk]-ubot[ij]) );
            const TF v2 = fmax(minval, fm::pow2(v[ijk]-vbot[ij]) );

            const TF ustaronu4 = TF(0.5)*(fm::pow4(ustar[ij-ii]) + fm::pow4(ustar[ij]));
            const TF ustaronv4 = TF(0.5)*(fm::pow4(ustar[ij-jj]) + fm::pow4(ustar[ij]));

            ufluxbot[ij] = -copysign(TF(1), u[ijk]-ubot[ij]) * pow(ustaronu4 / (TF(1) + vonu2 / u2), TF(0.5));
            vfluxbot[ij] = -copysign(TF(1), v[ijk]-vbot[ij]) * pow(ustaronv4 / (TF(1) + uonv2 / v2), TF(0.5));

            ugradbot[ij] = (u[ijk]-ubot[ij])/zsl;
            vgradbot[ij] = (v[ijk]-vbot[ij])/zsl;
        }
    }

    template<typename TF> __global__
    void set_bcs_thl_qt_g(
            TF* const restrict thl_gradbot,
            TF* const restrict qt_gradbot,
            const TF* const restrict thl,
            const TF* const restrict qt,
            const TF* const restrict thl_bot,
            const TF* const restrict qt_bot,
            const TF zsl, const int kstart,
            const int icells, const int jcells,
            const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*icells;
            const int ijk = i + j*icells + kstart*ijcells;

            // Use the linearly interpolated grad, rather than the MO grad,
            // to prevent giving unresolvable gradients to advection schemes
            thl_gradbot[ij] = (thl[ijk]-thl_bot[ij])/zsl;
            qt_gradbot[ij]  = (qt [ijk]-qt_bot [ij])/zsl;
        }
    }

    template<typename TF> __global__
    void set_bcs_scalars_dirichlet_g(
            TF* const restrict varbot,
            TF* const restrict vargradbot,
            TF* const restrict varfluxbot,
            const TF* const restrict ustar,
            const TF* const restrict obuk,
            const TF* const restrict var,
            const TF* const restrict z0h,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend, const int kstart,
            const int icells, const int jcells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*icells;
            const int ijk = i + j*icells + kstart*ijcells;

            varfluxbot[ij] = -(var[ijk]-varbot[ij])*ustar[ij]*most::fh(zsl, z0h[ij], obuk[ij]);
            vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
        }
    }

template<typename TF> __global__
    void set_bcs_scalars_flux_g(
            TF* const restrict varbot,
            TF* const restrict vargradbot,
            TF* const restrict varfluxbot,
            const TF* const restrict ustar,
            const TF* const restrict obuk,
            const TF* const restrict var,
            const TF* const restrict z0h,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend, const int kstart,
            const int icells, const int jcells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*icells;
            const int ijk = i + j*icells + kstart*ijcells;

            varbot[ij] = varfluxbot[ij] / (ustar[ij]*most::fh(zsl, z0h[ij], obuk[ij])) + var[ijk];
            vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
        }
    }


//    template<typename TF>
//    void calc_bcs(
//            TF* const __restrict__ thl_bot,
//            TF* const __restrict__ qt_bot,
//            const TF* const __restrict__ thl,
//            const TF* const __restrict__ qt,
//            const TF* const __restrict__ H_veg,
//            const TF* const __restrict__ H_soil,
//            const TF* const __restrict__ H_wet,
//            const TF* const __restrict__ LE_veg,
//            const TF* const __restrict__ LE_soil,
//            const TF* const __restrict__ LE_wet,
//            const TF* const __restrict__ tile_frac_veg,
//            const TF* const __restrict__ tile_frac_soil,
//            const TF* const __restrict__ tile_frac_wet,
//            const TF* const __restrict__ ra,
//            const TF* const __restrict__ rhorefh,
//            const int istart, const int iend,
//            const int jstart, const int jend,
//            const int kstart,
//            const int icells, const int ijcells)
//    {
//        const TF rhocp_i = TF(1) / (rhorefh[kstart] * cp<TF>);
//        const TF rholv_i = TF(1) / (rhorefh[kstart] * Lv<TF>);
//
//        for (int j=jstart; j<jend; ++j)
//            #pragma ivdep
//            for (int i=istart; i<iend; ++i)
//            {
//                const int ij  = i + j*icells;
//                const int ijk = ij + kstart*ijcells;
//
//                // Tile averaged surface fluxes
//                const TF wthl = (
//                    tile_frac_veg [ij] * H_veg [ij] +
//                    tile_frac_soil[ij] * H_soil[ij] +
//                    tile_frac_wet [ij] * H_wet [ij] ) * rhocp_i;
//
//                const TF wqt = (
//                    tile_frac_veg [ij] * LE_veg [ij] +
//                    tile_frac_soil[ij] * LE_soil[ij] +
//                    tile_frac_wet [ij] * LE_wet [ij] ) * rholv_i;
//
//                // Calculate surface values
//                thl_bot[ij] = thl[ijk] + wthl * ra[ij];
//                qt_bot [ij] = qt[ijk]  + wqt  * ra[ij];
//            }
//    }

    template<typename TF> __global__
    void set_water_tiles(
            TF* const __restrict__ c_veg,
            TF* const __restrict__ c_soil,
            TF* const __restrict__ c_wet,
            TF* const __restrict__ H_veg,
            TF* const __restrict__ H_soil,
            TF* const __restrict__ H_wet,
            TF* const __restrict__ LE_veg,
            TF* const __restrict__ LE_soil,
            TF* const __restrict__ LE_wet,
            TF* const __restrict__ G_veg,
            TF* const __restrict__ G_soil,
            TF* const __restrict__ G_wet,
            TF* const __restrict__ rs_veg,
            TF* const __restrict__ rs_soil,
            TF* const __restrict__ rs_wet,
            TF* const __restrict__ thl_bot_wet,
            TF* const __restrict__ qt_bot_wet,
            const int* const __restrict__ water_mask,
            const TF* const __restrict__ t_bot_water,
            const TF* const __restrict__ thl,
            const TF* const __restrict__ qt,
            const TF* const __restrict__ thl_bot,
            const TF* const __restrict__ qt_bot,
            const TF* const __restrict__ ra,
            const TF rho_bot,
            const TF pref_bot,
            const TF exner_bot,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        const TF rhocp = rho_bot * Constants::cp<TF>;
        const TF rholv = rho_bot * Constants::Lv<TF>;

        if (i < iend && j < jend)
        {
            const int ij = i + j*icells;
            const int ijk = ij + kstart*ijcells;

            if (water_mask[ij])
            {
                thl_bot_wet[ij] = t_bot_water[ij] / exner_bot;
                qt_bot_wet[ij]  = Thermo_moist_functions::qsat(pref_bot, t_bot_water[ij]);

                c_veg[ij]  = TF(0);
                c_soil[ij] = TF(0);
                c_wet[ij]  = TF(1);

                H_veg[ij]  = TF(0);
                H_soil[ij] = TF(0);
                H_wet[ij]  = rhocp / ra[ij] * (thl_bot_wet[ij] - thl[ijk]);

                LE_veg[ij]  = TF(0);
                LE_soil[ij] = TF(0);
                LE_wet[ij]  = rholv / ra[ij] * (qt_bot_wet[ij] - qt[ijk]);

                G_veg[ij]  = TF(0);
                G_soil[ij] = TF(0);
                G_wet[ij]  = TF(0);

                rs_veg[ij]  = TF(0);
                rs_soil[ij] = TF(0);
                rs_wet[ij]  = TF(0);
            }
        }
    }

    template<typename TF> __global__
    void calc_tiled_mean_g(
            TF* const __restrict__ fld_mean,
            const TF* const __restrict__ fld_veg,
            const TF* const __restrict__ fld_soil,
            const TF* const __restrict__ fld_wet,
            const TF* const __restrict__ tile_frac_veg,
            const TF* const __restrict__ tile_frac_soil,
            const TF* const __restrict__ tile_frac_wet,
            const TF fac,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij = i + j*icells;

            fld_mean[ij] = (
                tile_frac_veg [ij] * fld_veg [ij] +
                tile_frac_soil[ij] * fld_soil[ij] +
                tile_frac_wet [ij] * fld_wet [ij] ) * fac;
        }
    }

    template<typename TF> __global__
    void scale_tile_with_fraction_g(
            TF* const __restrict__ fld_scaled,
            const TF* const __restrict__ fld,
            const TF* const __restrict__ tile_frac,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij = i + j*icells;
            fld_scaled[ij] = fld[ij] * tile_frac[ij];
        }
    }

    template<typename TF> __global__
    void print_ij(
            TF* const __restrict__ fld,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*icells;
            printf("%d %d %f\n", i, j, fld[ij]);
        }
    }

    template<typename TF>
    void init_tile(Surface_tile<TF>& tile, const int ijcells)
    {
        const int memsize_tf  = ijcells * sizeof(TF);
        const int memsize_int = ijcells * sizeof(int);

        cuda_safe_call(cudaMalloc(&tile.fraction_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.thl_bot_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.qt_bot_g, memsize_tf));

        cuda_safe_call(cudaMalloc(&tile.obuk_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.ustar_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.bfluxbot_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.ra_g, memsize_tf));

        cuda_safe_call(cudaMalloc(&tile.nobuk_g, memsize_int));

        cuda_safe_call(cudaMalloc(&tile.rs_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.H_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.LE_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.G_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.S_g, memsize_tf));
    }

    template<typename TF>
    void clear_tile(Surface_tile<TF>& tile)
    {
        cuda_safe_call(cudaFree(tile.fraction_g));
        cuda_safe_call(cudaFree(tile.thl_bot_g));
        cuda_safe_call(cudaFree(tile.qt_bot_g));

        cuda_safe_call(cudaFree(tile.obuk_g));
        cuda_safe_call(cudaFree(tile.ustar_g));
        cuda_safe_call(cudaFree(tile.bfluxbot_g));
        cuda_safe_call(cudaFree(tile.ra_g));

        cuda_safe_call(cudaFree(tile.nobuk_g));

        cuda_safe_call(cudaFree(tile.rs_g));
        cuda_safe_call(cudaFree(tile.H_g));
        cuda_safe_call(cudaFree(tile.LE_g));
        cuda_safe_call(cudaFree(tile.G_g));
        cuda_safe_call(cudaFree(tile.S_g));
    }

    template<typename TF>
    void forward_device_tile(Surface_tile<TF>& tile, const int ijcells)
    {
        const int memsize_tf  = ijcells * sizeof(TF);
        const int memsize_int = ijcells * sizeof(int);

        cuda_safe_call(cudaMemcpy(tile.fraction_g, tile.fraction.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.thl_bot_g, tile.thl_bot.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.qt_bot_g, tile.qt_bot.data(), memsize_tf, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(tile.obuk_g, tile.obuk.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.ustar_g, tile.ustar.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.bfluxbot_g, tile.bfluxbot.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.ra_g, tile.ra.data(), memsize_tf, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(tile.nobuk_g, tile.nobuk.data(), memsize_int, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(tile.rs_g, tile.rs.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.H_g, tile.H.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.LE_g, tile.LE.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.G_g, tile.G.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.S_g, tile.S.data(), memsize_tf, cudaMemcpyHostToDevice));
    }

    template<typename TF>
    void backward_device_tile(Surface_tile<TF>& tile, const int ijcells)
    {
        const int memsize_tf  = ijcells * sizeof(TF);
        const int memsize_int = ijcells * sizeof(int);

        cuda_safe_call(cudaMemcpy(tile.fraction.data(), tile.fraction_g, memsize_tf, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(tile.thl_bot.data(), tile.thl_bot_g, memsize_tf, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(tile.qt_bot.data(), tile.qt_bot_g, memsize_tf, cudaMemcpyDeviceToHost));

        cuda_safe_call(cudaMemcpy(tile.obuk.data(), tile.obuk_g, memsize_tf, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(tile.ustar.data(), tile.ustar_g, memsize_tf, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(tile.bfluxbot.data(), tile.bfluxbot_g, memsize_tf, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(tile.ra.data(), tile.ra_g, memsize_tf, cudaMemcpyDeviceToHost));

        cuda_safe_call(cudaMemcpy(tile.nobuk.data(), tile.nobuk_g, memsize_int, cudaMemcpyDeviceToHost));

        cuda_safe_call(cudaMemcpy(tile.rs.data(), tile.rs_g, memsize_tf, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(tile.H.data(), tile.H_g, memsize_tf, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(tile.LE.data(), tile.LE_g, memsize_tf, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(tile.G.data(), tile.G_g, memsize_tf, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(tile.S.data(), tile.S_g, memsize_tf, cudaMemcpyDeviceToHost));
    }
}
#endif
