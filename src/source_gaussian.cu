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

#include <iostream>
#include <cmath>
#include <vector>
#include <array>

#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "fast_math.h"
#include "tools.h"
#include "constants.h"

#include "source_gaussian.h"
#include "source_gaussian_kernels.h"
#include "source_gaussian_kernels.cuh"

namespace sgk = Source_gaussian_kernels;        // CPU kernels.
namespace sgk_g = Source_gaussian_kernels_g;    // GPU kernels.


#ifdef USECUDA
template<typename TF>
void Source_gaussian<TF>::exec(Thermo<TF>& thermo, Timeloop<TF>& timeloop)
{
    auto& gd = grid.get_grid_data();

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

        sgk_g::calc_source_g<TF><<<gridGPU, blockGPU>>>(
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


template<typename TF>
void Source_gaussian<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    auto& gd = grid.get_grid_data();

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
            shape[n].range_x = sgk::calc_shape(gd.x.data(), source_x0[n], sigma_x[n], line_x[n], gd.istart, gd.iend);
            shape[n].range_y = sgk::calc_shape(gd.y.data(), source_y0[n], sigma_y[n], line_y[n], gd.jstart, gd.jend);
            shape[n].range_z = sgk::calc_shape(gd.z.data(), source_z0[n], sigma_z[n], line_z[n], gd.kstart, gd.kend);

            if (sw_emission_profile)
                throw std::runtime_error("Emission profiles with time dependent location/strength are not (yet) supported!");
            else
            {
                sgk_g::calc_norm_g<TF><<<gridGPU, blockGPU>>>(
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
}


template<typename TF>
void Source_gaussian<TF>::prepare_device()
{
}
#endif


#ifdef FLOAT_SINGLE
template class Source_gaussian<float>;
#else
template class Source_gaussian<double>;
#endif
