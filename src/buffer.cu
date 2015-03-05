/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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

#include <cstdio>
#include <cmath>
#include <stdlib.h>
#include "grid.h"
#include "fields.h"
#include "buffer.h"
#include "constants.h"
#include "tools.h"

namespace
{
    __global__ 
    void buffer_g(double* __restrict__ at,   double* __restrict__ a,
                  double* __restrict__ abuf, double* __restrict__ z,
                  double zstart, double zsizebufi, double sigma,  double beta,
                  int istart, int jstart, int bufferkstart,
                  int iend,   int jend,   int kend,
                  int jj, int kk)
    {
        __shared__ double sigmaz;

        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
        const int k = blockIdx.z + bufferkstart; 

        /* sigmaz only depends on height. Let one thread calculate it to shared memory,
           other threads re-use value */
        if (threadIdx.x == 0 && threadIdx.y == 0)
            sigmaz = sigma * pow((z[k]-zstart)*zsizebufi, beta);
        __syncthreads();

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            at[ijk] -= sigmaz*(a[ijk]-abuf[k]);
        }
    }
}

void Buffer::prepare_device()
{
    if (swbuffer == "1")
    {
        const int nmemsize = grid->kcells*sizeof(double);

        // Allocate the buffer arrays at GPU.
        for (FieldMap::const_iterator it=fields->mp.begin(); it!=fields->mp.end(); ++it)
            cuda_safe_call(cudaMalloc(&bufferprofs_g[it->first], nmemsize));
        for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
            cuda_safe_call(cudaMalloc(&bufferprofs_g[it->first], nmemsize));

        // Copy buffers to GPU.
        for (FieldMap::const_iterator it=fields->mp.begin(); it!=fields->mp.end(); ++it)
            cuda_safe_call(cudaMemcpy(bufferprofs_g[it->first], bufferprofs[it->first], nmemsize, cudaMemcpyHostToDevice));
        for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
            cuda_safe_call(cudaMemcpy(bufferprofs_g[it->first], bufferprofs[it->first], nmemsize, cudaMemcpyHostToDevice));
    }
}

void Buffer::clear_device()
{
    if(swbuffer == "1")
    {
        for (FieldMap::const_iterator it=fields->mp.begin(); it!=fields->mp.end(); ++it)
            cuda_safe_call(cudaFree(bufferprofs_g[it->first]));
        for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
            cuda_safe_call(cudaFree(bufferprofs_g[it->first]));
    }
}

#ifdef USECUDA
void Buffer::exec()
{
    if (swbuffer == "1")
    {
        const int blocki = grid->ithread_block;
        const int blockj = grid->jthread_block;
        const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
        const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);
        const int gridk  = grid->kmax - bufferkstart + 1;

        dim3 gridGPU (gridi, gridj, gridk);
        dim3 blockGPU(blocki, blockj, 1);

        const int offs = grid->memoffset;
        const double zsizebufi = 1./(grid->zsize-zstart);

        buffer_g<<<gridGPU, blockGPU>>>(
            &fields->mt["u"]->data_g[offs], &fields->mp["u"]->data_g[offs],
            bufferprofs_g["u"], grid->z_g, 
            zstart, zsizebufi, sigma, beta, 
            grid->istart,  grid->jstart, bufferkstart,
            grid->iend,    grid->jend,   grid->kend,
            grid->icellsp, grid->ijcellsp);
        cuda_check_error();

        buffer_g<<<gridGPU, blockGPU>>>(
            &fields->mt["v"]->data_g[offs], &fields->mp["v"]->data_g[offs],
            bufferprofs_g["v"], grid->z_g, 
            zstart, zsizebufi, sigma, beta, 
            grid->istart,  grid->jstart, bufferkstart,
            grid->iend,    grid->jend,   grid->kend,
            grid->icellsp, grid->ijcellsp);
        cuda_check_error();

        buffer_g<<<gridGPU, blockGPU>>>(
            &fields->mt["w"]->data_g[offs], &fields->mp["w"]->data_g[offs],
            bufferprofs_g["w"], grid->zh_g, 
            zstart, zsizebufi, sigma, beta, 
            grid->istart,  grid->jstart, bufferkstart,
            grid->iend,    grid->jend,   grid->kend,
            grid->icellsp, grid->ijcellsp);
        cuda_check_error();

        for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
            buffer_g<<<gridGPU, blockGPU>>>(
                &fields->st[it->first]->data_g[offs], &it->second->data_g[offs],
                bufferprofs_g[it->first], grid->z_g, 
                zstart, zsizebufi, sigma, beta, 
                grid->istart,  grid->jstart, bufferkstart,
                grid->iend,    grid->jend,   grid->kend,
                grid->icellsp, grid->ijcellsp);
        cuda_check_error();
    }
}
#endif
