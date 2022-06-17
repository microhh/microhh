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

#ifndef SOIL_KERNELS_H
#define SOIL_KERNELS_H

#include "constants.h"
#include "boundary_surface_lsm.h"

using namespace Constants;

namespace Soil_kernels
{
    template<typename TF>
    inline TF calc_diffusivity_vg(
            const TF vg_a, const TF vg_l, const TF vg_m, const TF gamma_sat,
            const TF theta_res, const TF theta_sat, const TF theta_norm)
    {
        const TF vg_mi = TF(1) / vg_m;

        return (TF(1) - vg_m) * gamma_sat / (vg_a * vg_m * (theta_sat - theta_res))
                    * pow(theta_norm, (vg_l - vg_mi))
                    * (pow((TF(1) - pow(theta_norm, vg_mi)), -vg_m)
                    + pow((TF(1) - pow(theta_norm, vg_mi)), vg_m) - TF(2));
    }

    template<typename TF>
    inline TF calc_conductivity_vg(
            const TF theta_norm, const TF vg_l, const TF vg_m, const TF gamma_sat)
    {
        return gamma_sat * pow(theta_norm, vg_l)
                    * pow((TF(1) - pow((TF(1) - pow(theta_norm, (1. / vg_m))), vg_m)), 2);
    }

    template<typename TF>
    void init_soil_homogeneous(
            TF* const restrict soil_fld,
            const TF* const restrict soil_prof,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int isize, const int ijsize)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i+j * isize + k*ijsize;
                    soil_fld[ijk] = soil_prof[k-kstart];
                }
    }

    template<typename TF>
    void calc_soil_properties(
            TF* const restrict kappa_theta_min, TF* const restrict kappa_theta_max,
            TF* const restrict gamma_theta_min, TF* const restrict gamma_theta_max,
            TF* const restrict vg_m,
            TF* const restrict gamma_T_dry, TF* const restrict rho_C,
            const TF* const restrict vg_a,
            const TF* const restrict vg_l,
            const TF* const restrict vg_n,
            const TF* const restrict gamma_theta_sat,
            const TF* const restrict theta_res, const TF* const restrict theta_sat,
            const TF* const restrict theta_fc,
            const int table_size)
    {
        for (int i=0; i<table_size; ++i)
        {
            // van Genuchten parameter `m`
            vg_m[i] = (TF(1) - (TF(1) / vg_n[i]));

            // Min/max values diffusivity soil moisture
            const TF theta_norm_min =
                (TF(1.001) * theta_res[i] - theta_res[i]) / (theta_sat[i] - theta_res[i]);
            const TF theta_norm_max =
                (TF(0.999) * theta_sat[i] - theta_res[i]) / (theta_sat[i] - theta_res[i]);

            kappa_theta_min[i] = calc_diffusivity_vg(
                    vg_a[i], vg_l[i], vg_m[i], gamma_theta_sat[i],
                    theta_res[i], theta_sat[i], theta_norm_min);
            kappa_theta_max[i] = calc_diffusivity_vg(
                    vg_a[i], vg_l[i], vg_m[i], gamma_theta_sat[i],
                    theta_res[i], theta_sat[i], theta_norm_max);

            // Min/max values conductivity soil moisture
            gamma_theta_min[i] = TF(0);
            gamma_theta_max[i] = gamma_theta_sat[i];

            // Conductivity temperature
            const TF rho_solid = TF(2700);  // Density of dry solid soil (kg m-3); PL98, eq. 6
            const TF rho_dry = (TF(1) - theta_sat[i]) * rho_solid;  // Density of soil (kg m-3)

            gamma_T_dry[i] = (TF(0.135) * rho_dry + TF(64.7)) / (rho_solid - TF(0.947) * rho_dry);
            rho_C[i] = (TF(1) - theta_sat[i]) * Constants::rho_C_matrix<TF>
                    + theta_fc[i] * Constants::rho_C_water<TF>;
        }
    }

    template<typename TF>
    void calc_root_column(
            TF* const restrict root_frac,
            const TF* const restrict zh,
            const TF a_root, const TF b_root,
            const int kstart, const int kend)
    {
        TF root_frac_sum = TF(0);

        for (int k=kstart+1; k<kend; ++k)
        {
            root_frac[k] = 0.5 * (exp(a_root * zh[k+1]) + \
                                  exp(b_root * zh[k+1]) - \
                                  exp(a_root * zh[k  ]) - \
                                  exp(b_root * zh[k  ]));

            root_frac_sum += root_frac[k];
        }

        // Make sure the root fraction sums to one.
        root_frac[kstart] = TF(1) - root_frac_sum;
    }

    template<typename TF>
    void calc_root_fraction(
            TF* const restrict root_frac,
            const TF* const restrict a_root,
            const TF* const restrict b_root,
            const TF* const restrict zh,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        for (int j=jstart; j<jend; ++j)
            for (int i=istart; i<iend; ++i)
            {
                TF root_frac_sum = TF(0);

                for (int k=kstart+1; k<kend; ++k)
                {
                    const int ij  = i + j*icells;
                    const int ijk = i + j*icells + k*ijcells;

                    root_frac[ijk] = 0.5 * (exp(a_root[ij] * zh[k+1]) + \
                                            exp(b_root[ij] * zh[k+1]) - \
                                            exp(a_root[ij] * zh[k  ]) - \
                                            exp(b_root[ij] * zh[k  ]));

                    root_frac_sum += root_frac[ijk];
                }

                const int ijk = i +j*icells + kstart*ijcells;
                // Make sure the root fraction sums to one.
                root_frac[ijk] = TF(1) - root_frac_sum;
            }
    }

    template<typename TF>
    void calc_root_weighted_mean_theta(
            TF* const restrict theta_mean,
            const TF* const restrict theta,
            const int* const restrict soil_index,
            const TF* const restrict root_fraction,
            const TF* const restrict theta_wp,
            const TF* const restrict theta_fc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij = i + j*icells;
                theta_mean[ij] = TF(0);
            }

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*icells;
                    const int ijk = i + j*icells + k*ijcells;
                    const int si  = soil_index[ijk];

                    const TF theta_lim = std::max(theta[ijk], theta_wp[si]);
                    theta_mean[ij] += root_fraction[ijk]
                            * (theta_lim - theta_wp[si]) / (theta_fc[si] - theta_wp[si]);
                }
    }

    template<typename TF>
    void calc_thermal_properties(
            TF* const restrict kappa,
            TF* const restrict gamma,
            const int* const restrict soil_index,
            const TF* const restrict theta,
            const TF* const restrict theta_sat,
            const TF* const restrict gamma_dry,
            const TF* const restrict rho_C,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    const int si = soil_index[ijk];

                    // Heat conductivity at saturation (from IFS code..)
                    const TF gamma_T_sat = pow(Constants::gamma_T_matrix<TF>, (TF(1) - theta_sat[si]))
                                            * pow(Constants::gamma_T_water<TF>, theta[ijk])
                                            * pow(TF(2.2), (theta_sat[si] - theta[ijk]));

                    // Kersten number for fine soils [IFS eq 8.64] (-)
                    const TF kersten = log10(std::max(TF(0.1), theta[ijk] / theta_sat[si])) + TF(1);

                    // Heat conductivity soil [IFS eq 8.62] (W m-1 K-1)
                    gamma[ijk] = kersten * (gamma_T_sat - gamma_dry[si]) + gamma_dry[si];

                    // Heat diffusivity (m2 s-1)
                    kappa[ijk] = gamma[ijk] / rho_C[si];
                }
    }

    template<typename TF>
    void calc_hydraulic_properties(
            TF* const restrict kappa,
            TF* const restrict gamma,
            const int* const restrict soil_index,
            const TF* const restrict theta,
            const TF* const restrict theta_sat,
            const TF* const restrict theta_res,
            const TF* const restrict vg_a,
            const TF* const restrict vg_l,
            const TF* const restrict vg_m,
            const TF* const restrict gamma_sat,
            const TF* const restrict gamma_min,
            const TF* const restrict gamma_max,
            const TF* const restrict kappa_min,
            const TF* const restrict kappa_max,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    const int si = soil_index[ijk];

                    // Limit soil moisture just above the residual soil moisture content
                    const TF theta_lim = std::max(theta[ijk], TF(1.001) * theta_res[si]);

                    // Dimensionless soil water content
                    const TF theta_norm = (theta_lim - theta_res[si]) / (theta_sat[si] - theta_res[si]);

                    // Calculate & limit the diffusivity
                    kappa[ijk] = calc_diffusivity_vg(
                            vg_a[si], vg_l[si], vg_m[si], gamma_sat[si],
                            theta_res[si], theta_sat[si], theta_norm);
                    kappa[ijk] = std::max(std::min(kappa_max[si], kappa[ijk]), kappa_min[si]);

                    // Calculate & limit the conductivity
                    gamma[ijk] = calc_conductivity_vg(
                            theta_norm, vg_l[si], vg_m[si], gamma_sat[si]);
                    gamma[ijk] = std::max(std::min(gamma_max[si], gamma[ijk]), gamma_min[si]);
                }
    }

    template<typename TF>
    void calc_root_water_extraction(
            TF* const restrict extraction,
            TF* const restrict tmp,
            const TF* const restrict theta,
            const TF* const restrict root_frac,
            const TF* const restrict LE_veg,
            const TF* const restrict dzi,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const TF fac = TF(1) / (rho_w<TF> * Lv<TF>);

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                tmp[ij] = TF(0);
            }

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*icells;
                    const int ijk = ij + k*ijcells;

                    tmp[ij] += root_frac[ijk] * theta[ijk];
                }

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*icells;
                    const int ijk = ij + k*ijcells;

                    const TF theta_frac = root_frac[ijk] * theta[ijk] / tmp[ij];
                    extraction[ijk] = -std::max(TF(0), LE_veg[ij]) * fac * dzi[k] * theta_frac;
                }
    }

    template<typename TF>
    void calc_infiltration(
            TF* const restrict infiltration,
            TF* const restrict runoff,
            const TF* const restrict throughfall,
            const TF* const restrict theta,
            const TF* const restrict theta_sat,
            const TF* const restrict kappa_max,
            const TF* const restrict gamma_max,
            const TF* const restrict dz,
            const int* const restrict soil_index,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kend,
            const int icells, const int ijcells)
    {
        const TF dz2i = TF(1)/(TF(0.5)*dz[kend-1]);

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                const int ijk = i + j*icells + (kend-1)*ijcells;
                const int si  = soil_index[ijk];

                const TF i_max = std::min(TF(0),
                        -(kappa_max[si] * (theta_sat[si] - theta[ijk]) * dz2i + gamma_max[si]));

                infiltration[ij] = std::min(TF(0), std::max(throughfall[ij], i_max));
                runoff[ij]       = std::min(TF(0), throughfall[ij] - infiltration[ij]);
            }
    }

    template<typename TF, Soil_interpolation_type interpolation_type>
    void interp_2_vertical(
            TF* const restrict fldh,
            const TF* const restrict fld,
            const TF* const restrict dz,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int kk = ijcells;

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;

                    if (interpolation_type == Soil_interpolation_type::Mean)
                        fldh[ijk] = TF(0.5) * (fld[ijk] + fld[ijk-kk]);
                    else if(interpolation_type == Soil_interpolation_type::Max)
                        fldh[ijk] = std::max(fld[ijk], fld[ijk-kk]);
                    else if(interpolation_type == Soil_interpolation_type::Harmonic_mean)
                        fldh[ijk] = (dz[k-1]+dz[k])*(fld[ijk-kk]*fld[ijk]) /
                                (fld[ijk-kk]*dz[k] + fld[ijk]*dz[k-1]);
                }
    }

    template<typename TF>
    void set_bcs_temperature(
            TF* const restrict flux_top,
            TF* const restrict flux_bot,
            const TF* const restrict G,
            const TF* const restrict rho_C,
            const int* const restrict soil_index,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kend,
            const int icells, const int ijcells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                const int ijk = ij + (kend-1)*ijcells;  // Top soil layer
                const int si  = soil_index[ijk];

                flux_top[ij] = -G[ij] / rho_C[si];
                flux_bot[ij] = TF(0);
            }
    }

    template<typename TF>
    void set_bcs_moisture(
            TF* const restrict flux_top,
            TF* const restrict flux_bot,
            TF* const restrict conductivity_h,
            const TF* const restrict LE_soil,
            const TF* const restrict tile_frac_soil,
            const TF* const restrict infiltration,
            const bool sw_free_drainage,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const TF fac = TF(1) / (rho_w<TF> * Lv<TF>);

        const int kk = ijcells;

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij = i + j*icells;

                flux_top[ij] = tile_frac_soil[ij] * LE_soil[ij] * fac + infiltration[ij];
                flux_bot[ij] = TF(0);
            }

        if (sw_free_drainage)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + kstart*ijcells;
                    conductivity_h[ijk] = conductivity_h[ijk+kk];
                }
        }
        else
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + kstart*ijcells;
                    conductivity_h[ijk] = TF(0);
                }
        }
    }

    template<typename TF, bool sw_source_term, bool sw_conductivity_term>
    void diff_explicit(
            TF* const restrict tend,
            const TF* const restrict fld,
            const TF* const restrict kappa_h,
            const TF* const restrict gamma_h,
            const TF* const restrict source,
            const TF* const restrict flux_top,
            const TF* const restrict flux_bot,
            const TF* const restrict dzi,
            const TF* const restrict dzhi,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int kk = ijcells;
        int k;

        // Bottom soil level
        k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                const int ijk = ij + k*ijcells;

                tend[ijk] += ((kappa_h[ijk+kk] * (fld[ijk+kk] - fld[ijk]) * dzhi[k+1]) + flux_bot[ij])*dzi[k];

                if (sw_conductivity_term)
                    tend[ijk] += (gamma_h[ijk+kk] - gamma_h[ijk]) * dzi[k];
                if (sw_source_term)
                    tend[ijk] += source[ijk];
            }

        // Top soil level
        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                const int ijk = ij + k*ijcells;

                tend[ijk] += (-flux_top[ij] - (kappa_h[ijk] * (fld[ijk] - fld[ijk-kk]) * dzhi[k]))*dzi[k];

                if (sw_conductivity_term)
                    tend[ijk] -= gamma_h[ijk] * dzi[k];
                if (sw_source_term)
                    tend[ijk] += source[ijk];
            }

        // Interior
        for (int k=kstart+1; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;

                    tend[ijk] += ((kappa_h[ijk+kk] * (fld[ijk+kk] - fld[ijk   ]) * dzhi[k+1])
                               - (kappa_h[ijk   ] * (fld[ijk   ] - fld[ijk-kk]) * dzhi[k  ])) * dzi[k];

                    if (sw_conductivity_term)
                        tend[ijk] += (gamma_h[ijk+kk] - gamma_h[ijk]) * dzi[k];
                    if (sw_source_term)
                        tend[ijk] += source[ijk];
                }
    }
}
#endif
