/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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
#include "math.h"

namespace
{
    __global__ 
    void boundary_cyclic_x_g(double* const __restrict__ data,
                             const int icells, const int jcells, const int kcells,
                             const int icellsp,
                             const int istart, const int jstart,
                             const int iend,   const int jend, 
                             const int igc,    const int jgc)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;
        const int k = blockIdx.z;

        const int jj = icellsp;
        const int kk = icellsp*jcells;

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

    __global__ 
    void boundary_cyclic_y_g(double* const __restrict__ data,
                             const int icells, const int jcells, const int kcells,
                             const int icellsp,
                             const int istart, const int jstart,
                             const int iend,   const int jend, 
                             const int igc,    const int jgc)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;
        const int k = blockIdx.z;

        const int jj = icellsp;
        const int kk = icellsp*jcells;

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

void Grid::prepare_device()
{
    /* Align the interior of the grid (i.e. excluding ghost cells) with 
       the 128 byte memory blocks of the GPU's global memory */
    memoffset = 16 - igc;           // Padding at start of array 
    int padl  = 16-(int)imax%16;    // Elements left in last 128 byte block
    icellsp   = imax + padl + (padl < 2*igc) * 16;
    ijcellsp  = icellsp * jcells;  
    ncellsp   = ijcellsp * kcells + memoffset;

    // Calculate optimal size thread blocks based on grid
    ithread_block = min(256, 16 * ((itot / 16) + (itot % 16 > 0)));
    jthread_block = 256 / ithread_block;

    const int kmemsize = kcells*sizeof(double);

    cuda_safe_call(cudaMalloc((void**)&z_g    , kmemsize));
    cuda_safe_call(cudaMalloc((void**)&zh_g   , kmemsize));
    cuda_safe_call(cudaMalloc((void**)&dz_g   , kmemsize));
    cuda_safe_call(cudaMalloc((void**)&dzh_g  , kmemsize));
    cuda_safe_call(cudaMalloc((void**)&dzi_g  , kmemsize));
    cuda_safe_call(cudaMalloc((void**)&dzhi_g , kmemsize));
    cuda_safe_call(cudaMalloc((void**)&dzi4_g , kmemsize));
    cuda_safe_call(cudaMalloc((void**)&dzhi4_g, kmemsize));

    cuda_safe_call(cudaMemcpy(z_g    , z    , kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(zh_g   , zh   , kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(dz_g   , dz   , kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(dzh_g  , dzh  , kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(dzi_g  , dzi  , kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(dzhi_g , dzhi , kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(dzi4_g , dzi4 , kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(dzhi4_g, dzhi4, kmemsize, cudaMemcpyHostToDevice));
}

void Grid::clear_device()
{
    cuda_safe_call(cudaFree(z_g    ));
    cuda_safe_call(cudaFree(zh_g   ));
    cuda_safe_call(cudaFree(dz_g   ));
    cuda_safe_call(cudaFree(dzh_g  ));
    cuda_safe_call(cudaFree(dzi_g  ));
    cuda_safe_call(cudaFree(dzhi_g ));
    cuda_safe_call(cudaFree(dzi4_g ));
    cuda_safe_call(cudaFree(dzhi4_g));
}

void Grid::boundary_cyclic_g(double* data)
{
    const int blocki_x = igc;
    const int blockj_x = 256 / igc + (256%igc > 0);
    const int gridi_x  = 1;
    const int gridj_x  = jcells/blockj_x + (jcells%blockj_x > 0);

    const int blocki_y = 256 / jgc + (256%jgc > 0);
    const int blockj_y = jgc;
    const int gridi_y  = icells/blocki_y + (icells%blocki_y > 0);
    const int gridj_y  = 1;

    dim3 gridGPUx (gridi_x,  gridj_x,  kcells);
    dim3 blockGPUx(blocki_x, blockj_x, 1);

    dim3 gridGPUy (gridi_y,  gridj_y,  kcells);
    dim3 blockGPUy(blocki_y, blockj_y, 1);

    boundary_cyclic_x_g<<<gridGPUx,blockGPUx>>>(
        data, icells, jcells, kcells, icellsp,
        istart, jstart, iend, jend, igc, jgc);

    boundary_cyclic_y_g<<<gridGPUy,blockGPUy>>>(
        data, icells, jcells, kcells, icellsp,
        istart, jstart, iend, jend, igc, jgc);

    cuda_check_error();
}

void Grid::boundary_cyclic2d_g(double* data)
{
    const int blocki_x = igc;
    const int blockj_x = 256 / igc + (256%igc > 0);
    const int gridi_x  = 1;
    const int gridj_x  = jcells/blockj_x + (jcells%blockj_x > 0);

    const int blocki_y = 256 / jgc + (256%jgc > 0);
    const int blockj_y = jgc;
    const int gridi_y  = icells/blocki_y + (icells%blocki_y > 0);
    const int gridj_y  = 1;

    dim3 gridGPUx (gridi_x,  gridj_x,  1);
    dim3 blockGPUx(blocki_x, blockj_x, 1);

    dim3 gridGPUy (gridi_y,  gridj_y,  1);
    dim3 blockGPUy(blocki_y, blockj_y, 1);

    boundary_cyclic_x_g<<<gridGPUx,blockGPUx>>>(
        data, icells, jcells, kcells, icellsp,
        istart, jstart, iend, jend, igc, jgc);

    boundary_cyclic_y_g<<<gridGPUy,blockGPUy>>>(
        data, icells, jcells, kcells, icellsp,
        istart, jstart, iend, jend, igc, jgc);

    cuda_check_error();
}


double Grid::get_max_g(double* data, double* tmp)
{
    using namespace Tools_g;

    const double scalefac = 1.;
    double maxvalue;

    // Reduce 3D field excluding ghost cells and padding to jtot*ktot values
    reduce_interior(data, tmp, itot, istart, iend, jtot, jstart, jend, ktot, kstart, icellsp, ijcellsp, maxType);
    // Reduce jtot*ktot to ktot values
    reduce_all     (tmp, &tmp[jtot*ktot], jtot*ktot, ktot, jtot, maxType, scalefac);
    // Reduce ktot values to a single value
    reduce_all     (&tmp[jtot*ktot], tmp, ktot, 1, ktot, maxType, scalefac);
    // Copy back result from GPU
    cuda_safe_call(cudaMemcpy(&maxvalue, &tmp[0], sizeof(double), cudaMemcpyDeviceToHost));

    return maxvalue;
}

double Grid::get_sum_g(double* data, double* tmp)
{
    using namespace Tools_g;

    const double scalefac = 1.;
    double sumvalue;

    // Reduce 3D field excluding ghost cells and padding to jtot*ktot values
    reduce_interior(data, tmp, itot, istart, iend, jtot, jstart, jend, ktot, kstart, icellsp, ijcellsp, sumType);
    // Reduce jtot*ktot to ktot values
    reduce_all     (tmp, &tmp[jtot*ktot], jtot*ktot, ktot, jtot, sumType, scalefac);
    // Reduce ktot values to a single value
    reduce_all     (&tmp[jtot*ktot], tmp, ktot, 1, ktot, sumType, scalefac);
    // Copy back result from GPU
    cuda_safe_call(cudaMemcpy(&sumvalue, &tmp[0], sizeof(double), cudaMemcpyDeviceToHost));

    return sumvalue;
}

void Grid::calc_mean_g(double* prof, double* data, double* tmp)
{
    using namespace Tools_g;

    const double scalefac = 1./(itot*jtot);

    // Reduce 3D field excluding ghost cells and padding to jtot*kcells values
    reduce_interior(data, tmp, itot, istart, iend, jtot, jstart, jend, kcells, 0, icellsp, ijcellsp, sumType);
    // Reduce jtot*kcells to kcells values
    reduce_all     (tmp, prof, jtot*kcells, kcells, jtot, sumType, scalefac);
} 
