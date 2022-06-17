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
#include <array>

#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "source.h"
#include "defines.h"
#include "fast_math.h"
#include "tools.h"

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
}

// Add the source to the fields. This function is called in the main time loop.
#ifdef USECUDA
template<typename TF>
void Source<TF>::exec(Timeloop<TF>& timeloop)
{
    auto& gd = grid.get_grid_data();

    if (!swsource)
        return;

    if (swtimedep_location || swtimedep_strength)
        throw std::runtime_error("Time varying and/or moving sources not (yet) implemented on GPU!");

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

template class Source<double>;
template class Source<float>;
