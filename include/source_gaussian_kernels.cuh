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

#ifndef SOURCE_GAUSSIAN_KERNELS_G_H
#define SOURCE_GAUSSIAN_KERNELS_G_H

#include "fast_math.h"
#include "constants.h"

namespace fm = Fast_math;

namespace Source_gaussian_kernels_g
{
    template<typename TF> __global__
    void calc_source_g(
            TF* const __restrict__ st,
            const TF* const restrict x,
            const TF* const restrict y,
            const TF* const restrict z,
            const TF x0, const TF sigma_x, const TF line_x,
            const TF y0, const TF sigma_y, const TF line_y,
            const TF z0, const TF sigma_z, const TF line_z,
            const TF strength, const TF norm,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k  = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*icells + k*ijcells;

            if (line_x != 0)
            {
                if (x[i] >= x0+line_x)
                    st[ijk] += strength/norm*exp(
                            - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                            - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                            - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                else if (x[i] <= x0)
                    st[ijk] += strength/norm*exp(
                            - fm::pow2(x[i]-x0)       /fm::pow2(sigma_x)
                            - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                            - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                else
                    st[ijk] += strength/norm*exp(
                            - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                            - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
            }
            else if (line_y != 0)
            {
                if (y[j] >= y0+line_y)
                    st[ijk] += strength/norm*exp(
                            - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                            - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                            - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                else if (y[j] <= y0)
                    st[ijk] += strength/norm*exp(
                            - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                            - fm::pow2(y[j]-y0)       /fm::pow2(sigma_y)
                            - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                else
                    st[ijk] += strength/norm*exp(
                            - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                            - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
            }
            else if (line_z != 0)
            {
                if (z[k] >= z0+line_z)
                    st[ijk] += strength/norm*exp(
                            - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                            - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                            - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                else if (z[k] <= z0)
                    st[ijk] += strength/norm*exp(
                            - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                            - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                            - fm::pow2(z[k]-z0)       /fm::pow2(sigma_z));
                else
                    st[ijk] += strength/norm*exp(
                            - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                            - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y));
            }
            else
                st[ijk] += strength/norm*exp(
                        - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                        - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                        - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
        }
    }


    template<typename TF> __global__
    void calc_norm_g(
            TF* sum,
            const TF* const restrict x, const TF x0, const TF sigma_x, const TF line_x,
            const TF* const restrict y, const TF y0, const TF sigma_y, const TF line_y,
            const TF* const restrict z, const TF z0, const TF sigma_z, const TF line_z,
            const TF* const restrict emission_profile,
            const int range_x0, const int range_x1,
            const int range_y0, const int range_y1,
            const int range_z0, const int range_z1,
            const TF* const restrict dz,
            const TF dx, const TF dy,
            const TF* const restrict rhoref, const bool sw_vmr, const bool use_profile,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        namespace fm = Fast_math;

        *sum = 0.;
        TF blob_norm = 0.;
        TF scaling;

        const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k  = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jstride + k*kstride;

            if (sw_vmr)
                // Emissions come in [kmol tracers s-1] and are added to grid boxes in [VMR s-1] unit.
                // rhoref [kg m-3] divided by xmair [kg kmol-1] transfers to units [kmol(tracer) / kmol(air) / s].
                scaling = rhoref[k]/Constants::xmair<TF>;
            else
                // Emissions come in [kg tracer s-1]. [kg tracer s-1 / (m3 * kg m-3)] results in
                // emissions in units [kg tracer / kg air / s].
                scaling = rhoref[k];

            if (i>=range_x0 && i<=range_x1 && j>=range_y0 && j<=range_y1 && k>=range_z0 && k<=range_z1)
            {
                if (line_x != 0)
                {
                    if (x[i] >= x0+line_x)
                        blob_norm = exp(
                                - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                    else if (x[i]<=x0)
                        blob_norm = exp(
                                - fm::pow2(x[i]-x0)       /fm::pow2(sigma_x)
                                - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                    else
                        blob_norm = exp(
                                - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                }
                else if (line_y != 0)
                {
                    if (y[j] >= y0+line_y)
                        blob_norm = exp(
                                - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                    else if (y[j]<=y0)
                        blob_norm = exp(
                                - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                - fm::pow2(y[j]-y0)       /fm::pow2(sigma_y)
                                - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                    else
                        blob_norm = exp(
                                - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));

                }
                else if (line_z != 0)
                {
                    if (z[k] >= z0+line_z)
                        blob_norm = exp(
                                - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                    else if (z[k]<=z0)
                        blob_norm = exp(
                                - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                - fm::pow2(z[k]-z0)       /fm::pow2(sigma_z));
                    else
                        blob_norm = exp(
                                - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y));
                }
                else
                {
                    if (use_profile)
                        blob_norm = exp(
                                - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y))
                                    * emission_profile[k];
                    else
                        blob_norm = exp(
                                - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                }
            }
            else
                blob_norm = TF(0);

            atomicAdd(sum, blob_norm * dx * dy * dz[k] * scaling);
        }
    }
}
#endif
