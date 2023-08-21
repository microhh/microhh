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

#include "grid.h"
#include "tools.h"
#include "boundary_cyclic.h"

namespace
{
    template<typename TF> __global__
    void boundary_cyclic_x_g(TF* const __restrict__ data,
                             const int icells, const int jcells, const int kcells,
                             const int istart, const int jstart,
                             const int iend,   const int jend,
                             const int igc,    const int jgc)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;
        const int k = blockIdx.z;

        const int jj = icells;
        const int kk = icells*jcells;

        // East-west
        if (k < kcells && j < jcells && i < igc)
        {
            const int ijk0 = i          + j*jj + k*kk;
            const int ijk1 = iend-igc+i + j*jj + k*kk;
            const int ijk2 = i+iend     + j*jj + k*kk;
            const int ijk3 = i+istart   + j*jj + k*kk;

            data[ijk0] = data[ijk1];
            data[ijk2] = data[ijk3];
        }
    }

    template<typename TF> __global__
    void boundary_cyclic_y_g(TF* const __restrict__ data,
                             const int icells, const int jcells, const int kcells,
                             const int istart, const int jstart,
                             const int iend,   const int jend,
                             const int igc,    const int jgc)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;
        const int k = blockIdx.z;

        const int jj = icells;
        const int kk = icells*jcells;

        // North-south
        if (jend-jstart == 1)
        {
            if (k < kcells && j < jgc && i < icells)
            {
                const int ijkref   = i + jstart*jj   + k*kk;
                const int ijknorth = i + j*jj        + k*kk;
                const int ijksouth = i + (jend+j)*jj + k*kk;
                data[ijknorth] = data[ijkref];
                data[ijksouth] = data[ijkref];
            }
        }
        else
        {
            if (k < kcells && j < jgc && i < icells)
            {
                const int ijk0 = i + j           *jj + k*kk;
                const int ijk1 = i + (jend-jgc+j)*jj + k*kk;
                const int ijk2 = i + (j+jend  )  *jj + k*kk;
                const int ijk3 = i + (j+jstart)  *jj + k*kk;

                data[ijk0] = data[ijk1];
                data[ijk2] = data[ijk3];
            }
        }
    }
}

template<typename TF>
void Boundary_cyclic<TF>::exec_g(TF* data)
{
    auto& gd = grid.get_grid_data();

    const int blocki_x = gd.igc;
    const int blockj_x = 256 / gd.igc + (256%gd.igc > 0);
    const int gridi_x  = 1;
    const int gridj_x  = gd.jcells/blockj_x + (gd.jcells%blockj_x > 0);

    const int blocki_y = 256 / gd.jgc + (256%gd.jgc > 0);
    const int blockj_y = gd.jgc;
    const int gridi_y  = gd.icells/blocki_y + (gd.icells%blocki_y > 0);
    const int gridj_y  = 1;

    dim3 gridGPUx (gridi_x,  gridj_x,  gd.kcells);
    dim3 blockGPUx(blocki_x, blockj_x, 1);

    dim3 gridGPUy (gridi_y,  gridj_y,  gd.kcells);
    dim3 blockGPUy(blocki_y, blockj_y, 1);

    boundary_cyclic_x_g<TF><<<gridGPUx,blockGPUx>>>(
        data, gd.icells, gd.jcells, gd.kcells,
        gd.istart, gd.jstart, gd.iend, gd.jend, gd.igc, gd.jgc);

    boundary_cyclic_y_g<TF><<<gridGPUy,blockGPUy>>>(
        data, gd.icells, gd.jcells, gd.kcells,
        gd.istart, gd.jstart, gd.iend, gd.jend, gd.igc, gd.jgc);

    cuda_check_error();
}

template<typename TF>
void Boundary_cyclic<TF>::exec_2d_g(TF* data)
{
    auto& gd = grid.get_grid_data();

    const int blocki_x = gd.igc;
    const int blockj_x = 256 / gd.igc + (256%gd.igc > 0);
    const int gridi_x  = 1;
    const int gridj_x  = gd.jcells/blockj_x + (gd.jcells%blockj_x > 0);

    const int blocki_y = 256 / gd.jgc + (256%gd.jgc > 0);
    const int blockj_y = gd.jgc;
    const int gridi_y  = gd.icells/blocki_y + (gd.icells%blocki_y > 0);
    const int gridj_y  = 1;

    dim3 gridGPUx (gridi_x,  gridj_x,  1);
    dim3 blockGPUx(blocki_x, blockj_x, 1);

    dim3 gridGPUy (gridi_y,  gridj_y,  1);
    dim3 blockGPUy(blocki_y, blockj_y, 1);

    boundary_cyclic_x_g<TF><<<gridGPUx,blockGPUx>>>(
        data, gd.icells, gd.jcells, gd.kcells,
        gd.istart, gd.jstart, gd.iend, gd.jend, gd.igc, gd.jgc);

    boundary_cyclic_y_g<TF><<<gridGPUy,blockGPUy>>>(
        data, gd.icells, gd.jcells, gd.kcells,
        gd.istart, gd.jstart, gd.iend, gd.jend, gd.igc, gd.jgc);

    cuda_check_error();
}

#ifdef FLOAT_SINGLE
template class Boundary_cyclic<float>;
#else
template class Boundary_cyclic<double>;
#endif
