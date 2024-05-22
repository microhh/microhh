/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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
#include <array>

#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "source.h"
#include "source_kernels.h"
#include "defines.h"
#include "fast_math.h"
#include "tools.h"
#include "constants.h"

namespace
{
    namespace fm = Fast_math;

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

// Add the source to the fields. This function is called in the main time loop.
#ifdef USECUDA
template<typename TF>
void Source<TF>::exec(Timeloop<TF>& timeloop)
{
    auto& gd = grid.get_grid_data();

    if (!swsource)
        return;

    if (swtimedep_location)
    {
        const int blocki = gd.ithread_block;
        const int blockj = gd.jthread_block;
        const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
        const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

        dim3 gridGPU (gridi, gridj, gd.kcells);
        dim3 blockGPU(blocki, blockj, 1);

        TF* norm_g;
        cudaMalloc(&norm_g, sizeof(TF));

        // Update source locations, and calculate new norm's
        for (int n=0; n<sourcelist.size(); ++n)
        {
            std::string name_x = "source_x0_" + std::to_string(n);
            std::string name_y = "source_y0_" + std::to_string(n);
            std::string name_z = "source_z0_" + std::to_string(n);

            tdep_source_x0.at(name_x)->update_time_dependent(source_x0[n], timeloop);
            tdep_source_y0.at(name_y)->update_time_dependent(source_y0[n], timeloop);
            tdep_source_z0.at(name_z)->update_time_dependent(source_z0[n], timeloop);

            // Shape of the source in each direction
            shape[n].range_x = calc_shape(gd.x.data(), source_x0[n], sigma_x[n], line_x[n], gd.istart, gd.iend);
            shape[n].range_y = calc_shape(gd.y.data(), source_y0[n], sigma_y[n], line_y[n], gd.jstart, gd.jend);
            shape[n].range_z = calc_shape(gd.z.data(), source_z0[n], sigma_z[n], line_z[n], gd.kstart, gd.kend);

            if (sw_emission_profile)
                throw std::runtime_error("Emission profiles with time dependent location/strength are not (yet) supported!");
            else
            {
                calc_norm_g<TF><<<gridGPU, blockGPU>>>(
                        norm_g,
                        gd.x_g, source_x0[n], sigma_x[n], line_x[n],
                        gd.y_g, source_y0[n], sigma_y[n], line_y[n],
                        gd.z_g, source_z0[n], sigma_z[n], line_z[n],
                        nullptr,
                        shape[n].range_x[0], shape[n].range_x[1],
                        shape[n].range_y[0], shape[n].range_y[1],
                        shape[n].range_z[0], shape[n].range_z[1],
                        gd.dz_g, gd.dx, gd.dy,
                        fields.rhoref_g, sw_vmr[n], false,
                        gd.istart, gd.iend,
                        gd.jstart, gd.jend,
                        gd.kstart, gd.kend,
                        gd.icells, gd.ijcells);

                cudaMemcpy(&norm[n], norm_g, sizeof(long), cudaMemcpyDeviceToHost);
            }
        }

        // Take sum over MPI tasks, in CPU version done inside kernel.
        master.sum(norm.data(), sourcelist.size());
    }

    if (swtimedep_strength)
    {
        // Update source locations, and calculate new norm's
        for (int n=0; n<sourcelist.size(); ++n)
        {
            std::string name_strength = "source_strength_" + std::to_string(n);
            tdep_source_strength.at(name_strength)->update_time_dependent(strength[n], timeloop);
        }
    }

    for (int n=0; n<sourcelist.size(); ++n)
    {
        const int range[3] = {
                shape[n].range_x[1]-shape[n].range_x[0],
                shape[n].range_y[1]-shape[n].range_y[0],
                shape[n].range_z[1]-shape[n].range_z[0]};

        if (range[0] == 0 || range[1] == 0 || range[2] == 0)
            continue;

        const int blocki = 16;
        const int blockj = 16;

        const int gridi  = range[0]/blocki + (range[0]%blocki > 0);
        const int gridj  = range[1]/blockj + (range[1]%blockj > 0);

        dim3 gridGPU (gridi, gridj, range[2]);
        dim3 blockGPU(blocki, blockj, 1);

        calc_source_g<TF><<<gridGPU, blockGPU>>>(
                fields.st[sourcelist[n]]->fld_g,
                gd.x_g, gd.y_g, gd.z_g,
                source_x0[n], sigma_x[n], line_x[n],
                source_y0[n], sigma_y[n], line_y[n],
                source_z0[n], sigma_z[n], line_z[n],
                strength[n], norm[n],
                shape[n].range_x[0], shape[n].range_x[1],
                shape[n].range_y[0], shape[n].range_y[1],
                shape[n].range_z[0], shape[n].range_z[1],
                gd.icells, gd.ijcells);

        cuda_check_error();
    }
}
#endif

#ifdef FLOAT_SINGLE
template class Source<float>;
#else
template class Source<double>;
#endif
