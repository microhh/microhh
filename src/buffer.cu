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

#include <cstdio>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "grid.h"
#include "fields.h"
#include "buffer.h"
#include "constants.h"
#include "stats.h"
#include "tools.h"

namespace
{
    template<typename TF>__global__
    void buffer_g(TF* __restrict__ at,   TF* __restrict__ a,
                  TF* __restrict__ abuf, TF* __restrict__ z,
                  TF zstart, TF zsizebufi, TF sigma,  TF beta,
                  int istart, int jstart, int bufferkstart,
                  int iend,   int jend,   int kend,
                  int jj, int kk)
    {
        __shared__ TF sigmaz;

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

template<typename TF>
void Buffer<TF>::prepare_device()
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    if (swbuffer)
    {
        const int nmemsize = gd.kcells*sizeof(TF);

        // Allocate the buffer arrays at GPU.
        for (auto& it : fields.ap)
        {
            bufferprofs_g.emplace(it.first, nullptr);
            cuda_safe_call(cudaMalloc(&bufferprofs_g.at(it.first), nmemsize));
        }

        // Copy buffers to GPU.
        for (auto& it : fields.ap)
            cuda_safe_call(cudaMemcpy(bufferprofs_g.at(it.first), bufferprofs.at(it.first).data(), nmemsize, cudaMemcpyHostToDevice));
    }
}

template<typename TF>
void Buffer<TF>::clear_device()
{
    if (swbuffer)
    {
        for (auto& it : bufferprofs_g)
            cuda_safe_call(cudaFree(it.second));
    }
}

#ifdef USECUDA
template<typename TF>
void Buffer<TF>::exec(Stats<TF>& stats)
{
    if (swbuffer)
    {
        const Grid_data<TF>& gd = grid.get_grid_data();

        const int blocki = gd.ithread_block;
        const int blockj = gd.jthread_block;
        const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
        const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);
        const int gridk  = gd.kmax - (bufferkstart-gd.kgc);

        dim3 gridGPU (gridi, gridj, gridk);
        dim3 blockGPU(blocki, blockj, 1);

        const TF zsizebufi = 1./(gd.zsize-zstart);

        if (swupdate)
        {
            buffer_g<TF><<<gridGPU, blockGPU>>>(
                fields.mt.at("u")->fld_g, fields.mp.at("u")->fld_g,
                fields.mp.at("u")->fld_mean_g, gd.z_g,
                zstart, zsizebufi, sigma, beta,
                gd.istart, gd.jstart, bufferkstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();

            buffer_g<TF><<<gridGPU, blockGPU>>>(
                fields.mt.at("v")->fld_g, fields.mp.at("v")->fld_g,
                fields.mp.at("v")->fld_mean_g, gd.z_g,
                zstart, zsizebufi, sigma, beta,
                gd.istart, gd.jstart, bufferkstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();

            buffer_g<TF><<<gridGPU, blockGPU>>>(
                fields.mt.at("w")->fld_g, fields.mp.at("w")->fld_g,
                fields.mp.at("w")->fld_mean_g, gd.zh_g,
                zstart, zsizebufi, sigma, beta,
                gd.istart, gd.jstart, bufferkstarth,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();

            for (auto& it : fields.sp)
                buffer_g<TF><<<gridGPU, blockGPU>>>(
                    fields.st.at(it.first)->fld_g, fields.sp.at(it.first)->fld_g,
                    fields.sp.at(it.first)->fld_mean_g, gd.z_g,
                    zstart, zsizebufi, sigma, beta,
                    gd.istart, gd.jstart, bufferkstart,
                    gd.iend,   gd.jend,   gd.kend,
                    gd.icells, gd.ijcells);
            cuda_check_error();
        }
        else
        {
            buffer_g<TF><<<gridGPU, blockGPU>>>(
                fields.mt.at("u")->fld_g, fields.mp.at("u")->fld_g,
                bufferprofs_g.at("u"), gd.z_g,
                zstart, zsizebufi, sigma, beta,
                gd.istart, gd.jstart, bufferkstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();

            buffer_g<TF><<<gridGPU, blockGPU>>>(
                fields.mt.at("v")->fld_g, fields.mp.at("v")->fld_g,
                bufferprofs_g.at("v"), gd.z_g,
                zstart, zsizebufi, sigma, beta,
                gd.istart, gd.jstart, bufferkstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();

            buffer_g<TF><<<gridGPU, blockGPU>>>(
                fields.mt.at("w")->fld_g, fields.mp.at("w")->fld_g,
                bufferprofs_g.at("w"), gd.zh_g,
                zstart, zsizebufi, sigma, beta,
                gd.istart, gd.jstart, bufferkstarth,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();

            for (auto& it : fields.sp)
                buffer_g<TF><<<gridGPU, blockGPU>>>(
                    fields.st.at(it.first)->fld_g, fields.sp.at(it.first)->fld_g,
                    bufferprofs_g.at(it.first), gd.z_g,
                    zstart, zsizebufi, sigma, beta,
                    gd.istart, gd.jstart, bufferkstart,
                    gd.iend,   gd.jend,   gd.kend,
                    gd.icells, gd.ijcells);
            cuda_check_error();
        }

        cudaDeviceSynchronize();

        stats.calc_tend(*fields.mt.at("u"), tend_name);
        stats.calc_tend(*fields.mt.at("v"), tend_name);
        stats.calc_tend(*fields.mt.at("w"), tend_name);
        for (auto it : fields.st)
            stats.calc_tend(*it.second, tend_name);
    }
}
#endif

template class Buffer<double>;
template class Buffer<float>;
