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

#ifndef SOIL_KERNELS_GPU_H
#define SOIL_KERNELS_GPU_H

#include "constants.h"
#include "boundary_surface_lsm.h"
#include "fast_math.h"


namespace Soil_kernels_g
{
    using namespace Constants;
    namespace fm = Fast_math;

    template<typename TF> __device__
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

    template<typename TF> __device__
    inline TF calc_conductivity_vg(
            const TF theta_norm, const TF vg_l, const TF vg_m, const TF gamma_sat)
    {
        return gamma_sat * pow(theta_norm, vg_l)
                    * fm::pow2((TF(1) - pow((TF(1) - pow(theta_norm, (1. / vg_m))), vg_m)));
    }

    template<typename TF> __global__
    void calc_root_weighted_mean_theta_g(
            TF* const __restrict__ theta_mean,
            const TF* const __restrict__ theta,
            const int* const __restrict__ soil_index,
            const TF* const __restrict__ root_fraction,
            const TF* const __restrict__ theta_wp,
            const TF* const __restrict__ theta_fc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*icells;
            theta_mean[ij] = TF(0);

            // Keep loop at GPU to simplify the reduction...
            for (int k=kstart; k<kend; ++k)
            {
                const int ijk = i + j*icells + k*ijcells;
                const int si  = soil_index[ijk];

                const TF theta_lim = fmax(theta[ijk], theta_wp[si]);
                theta_mean[ij] += root_fraction[ijk]
                        * (theta_lim - theta_wp[si]) / (theta_fc[si] - theta_wp[si]);
            }
        }
    }

    template<typename TF> __global__
    void calc_thermal_properties_g(
            TF* const __restrict__ kappa,
            TF* const __restrict__ gamma,
            const int* const __restrict__ soil_index,
            const TF* const __restrict__ theta,
            const TF* const __restrict__ theta_sat,
            const TF* const __restrict__ gamma_dry,
            const TF* const __restrict__ rho_C,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*icells + k*ijcells;
            const int si = soil_index[ijk];

            // Heat conductivity at saturation (from IFS code..)
            const TF gamma_T_sat = pow(Constants::gamma_T_matrix<TF>, (TF(1) - theta_sat[si]))
                                    * pow(Constants::gamma_T_water<TF>, theta[ijk])
                                    * pow(TF(2.2), (theta_sat[si] - theta[ijk]));

            // Kersten number for fine soils [IFS eq 8.64] (-)
            const TF kersten = log10(fmax(TF(0.1), theta[ijk] / theta_sat[si])) + TF(1);

            // Heat conductivity soil [IFS eq 8.62] (W m-1 K-1)
            gamma[ijk] = kersten * (gamma_T_sat - gamma_dry[si]) + gamma_dry[si];

            // Heat diffusivity (m2 s-1)
            kappa[ijk] = gamma[ijk] / rho_C[si];
        }
    }

    template<typename TF> __global__
    void calc_hydraulic_properties_g(
            TF* const __restrict__ kappa,
            TF* const __restrict__ gamma,
            const int* const __restrict__ soil_index,
            const TF* const __restrict__ theta,
            const TF* const __restrict__ theta_sat,
            const TF* const __restrict__ theta_res,
            const TF* const __restrict__ vg_a,
            const TF* const __restrict__ vg_l,
            const TF* const __restrict__ vg_m,
            const TF* const __restrict__ gamma_sat,
            const TF* const __restrict__ gamma_min,
            const TF* const __restrict__ gamma_max,
            const TF* const __restrict__ kappa_min,
            const TF* const __restrict__ kappa_max,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*icells + k*ijcells;
            const int si = soil_index[ijk];

            // Limit soil moisture just above the residual soil moisture content
            const TF theta_lim = fmax(theta[ijk], TF(1.001) * theta_res[si]);

            // Dimensionless soil water content
            const TF theta_norm = (theta_lim - theta_res[si]) / (theta_sat[si] - theta_res[si]);

            // Calculate & limit the diffusivity
            kappa[ijk] = calc_diffusivity_vg(
                    vg_a[si], vg_l[si], vg_m[si], gamma_sat[si],
                    theta_res[si], theta_sat[si], theta_norm);
            kappa[ijk] = fmax(fmin(kappa_max[si], kappa[ijk]), kappa_min[si]);

            // Calculate & limit the conductivity
            gamma[ijk] = calc_conductivity_vg(
                    theta_norm, vg_l[si], vg_m[si], gamma_sat[si]);
            gamma[ijk] = fmax(fmin(gamma_max[si], gamma[ijk]), gamma_min[si]);
        }
    }

    template<typename TF> __global__
    void calc_root_water_extraction_g(
            TF* const __restrict__ extraction,
            TF* const __restrict__ tmp,
            const TF* const __restrict__ theta,
            const TF* const __restrict__ root_frac,
            const TF* const __restrict__ LE_veg,
            const TF* const __restrict__ dzi,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        const TF fac = TF(1) / (rho_w<TF> * Lv<TF>);

        if (i < iend && j < jend)
        {
            const int ij  = i + j*icells;
            tmp[ij] = TF(0);

            // Keep loop at GPU to simplify reduction.
            for (int k=kstart; k<kend; ++k)
            {
                const int ijk = ij + k*ijcells;
                tmp[ij] += root_frac[ijk] * theta[ijk];
            }

            for (int k=kstart; k<kend; ++k)
            {
                const int ijk = ij + k*ijcells;
                const TF theta_frac = root_frac[ijk] * theta[ijk] / tmp[ij];
                extraction[ijk] = -fmax(TF(0), LE_veg[ij]) * fac * dzi[k] * theta_frac;
            }
        }
    }

    template<typename TF> __global__
    void calc_infiltration_g(
            TF* const __restrict__ infiltration,
            TF* const __restrict__ runoff,
            const TF* const __restrict__ throughfall,
            const TF* const __restrict__ theta,
            const TF* const __restrict__ theta_sat,
            const TF* const __restrict__ kappa_max,
            const TF* const __restrict__ gamma_max,
            const TF* const __restrict__ dz,
            const int* const __restrict__ soil_index,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        const TF dz2i = TF(1)/(TF(0.5)*dz[kend-1]);

        if (i < iend && j < jend)
        {
            const int ij = i + j*icells;
            const int ijk = i + j*icells + (kend-1)*ijcells;
            const int si = soil_index[ijk];

            const TF i_max =
                fmin(TF(0), -(kappa_max[si] * (theta_sat[si] - theta[ijk]) * dz2i + gamma_max[si]));

            infiltration[ij] = fmin(TF(0), fmax(throughfall[ij], i_max));
            runoff[ij]       = fmin(TF(0), throughfall[ij] - infiltration[ij]);
        }
    }

    template<typename TF, Soil_interpolation_type interpolation_type> __global__
    void interp_2_vertical_g(
            TF* const __restrict__ fldh,
            const TF* const __restrict__ fld,
            const TF* const __restrict__ dz,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart+1;

        const int kk = ijcells;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*icells + k*ijcells;

            if (interpolation_type == Soil_interpolation_type::Mean)
                fldh[ijk] = TF(0.5) * (fld[ijk] + fld[ijk-kk]);
            else if(interpolation_type == Soil_interpolation_type::Max)
                fldh[ijk] = fmax(fld[ijk], fld[ijk-kk]);
            else if(interpolation_type == Soil_interpolation_type::Harmonic_mean)
                fldh[ijk] = (dz[k-1]+dz[k])*(fld[ijk-kk]*fld[ijk]) /
                        (fld[ijk-kk]*dz[k] + fld[ijk]*dz[k-1]);
        }
    }

    template<typename TF> __global__
    void set_bcs_temperature_g(
            TF* const __restrict__ flux_top,
            TF* const __restrict__ flux_bot,
            const TF* const __restrict__ G,
            const TF* const __restrict__ rho_C,
            const int* const __restrict__ soil_index,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*icells;
            const int ijk = ij + (kend-1)*ijcells;  // Top soil layer
            const int si  = soil_index[ijk];

            flux_top[ij] = -G[ij] / rho_C[si];
            flux_bot[ij] = TF(0);
        }
    }

    template<typename TF> __global__
    void set_bcs_moisture_g(
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
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        const TF fac = TF(1) / (rho_w<TF> * Lv<TF>);
        const int kk = ijcells;

        if (i < iend && j < jend)
        {
            const int ij = i + j*icells;
            const int ijk = i + j*icells + kstart*ijcells;

            flux_top[ij] = tile_frac_soil[ij] * LE_soil[ij] * fac + infiltration[ij];
            flux_bot[ij] = TF(0);

            if (sw_free_drainage)
                conductivity_h[ijk] = conductivity_h[ijk+kk];
            else
                conductivity_h[ijk] = TF(0);
        }
    }

    template<typename TF, bool sw_source_term, bool sw_conductivity_term> __global__
    void diff_explicit_g(
            TF* const __restrict__ tend,
            const TF* const __restrict__ fld,
            const TF* const __restrict__ kappa_h,
            const TF* const __restrict__ gamma_h,
            const TF* const __restrict__ source,
            const TF* const __restrict__ flux_top,
            const TF* const __restrict__ flux_bot,
            const TF* const __restrict__ dzi,
            const TF* const __restrict__ dzhi,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int kk = ijcells;

        if (i < iend && j < jend && k < kend)
        {
            const int ij  = i + j*icells;
            const int ijk = ij + k*ijcells;

            if (k == kstart)
            {
                // Bottom soil level..
                tend[ijk] += ((kappa_h[ijk+kk] * (fld[ijk+kk] - fld[ijk]) * dzhi[k+1]) + flux_bot[ij])*dzi[k];

                if (sw_conductivity_term)
                    tend[ijk] += (gamma_h[ijk+kk] - gamma_h[ijk]) * dzi[k];
                if (sw_source_term)
                    tend[ijk] += source[ijk];
            }
            else if (k == kend-1)
            {
                // Top soil levels..
                tend[ijk] += (-flux_top[ij] - (kappa_h[ijk] * (fld[ijk] - fld[ijk-kk]) * dzhi[k]))*dzi[k];

                if (sw_conductivity_term)
                    tend[ijk] -= gamma_h[ijk] * dzi[k];
                if (sw_source_term)
                    tend[ijk] += source[ijk];
            }
            else
            {
                // Interior..
                tend[ijk] += ((kappa_h[ijk+kk] * (fld[ijk+kk] - fld[ijk   ]) * dzhi[k+1])
                           - (kappa_h[ijk   ] * (fld[ijk   ] - fld[ijk-kk]) * dzhi[k  ])) * dzi[k];

                if (sw_conductivity_term)
                    tend[ijk] += (gamma_h[ijk+kk] - gamma_h[ijk]) * dzi[k];
                if (sw_source_term)
                    tend[ijk] += source[ijk];
            }
        }
    }
}
#endif
