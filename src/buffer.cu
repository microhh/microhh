/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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

__global__ void buffer_buffer(double * __restrict__ at,   double * __restrict__ a,
                              double * __restrict__ abuf, double * __restrict__ z,
                              double zstart, double zsizebufi, double sigma,  double beta,
                              int istart, int jstart, int bufferkstart,
                              int iend,   int jend,   int kend,
                              int jj, int kk)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
  int k = blockIdx.z + bufferkstart; 

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    double sigmaz = sigma * pow((z[k]-zstart)*zsizebufi, beta);

    at[ijk] -= sigmaz*(a[ijk]-abuf[k]);
  }
}

// TODO: (also for pressure), deallocate fields on GPU...
int cbuffer::prepareGPU()
{
  if(swbuffer == "1")
  {
    const int nmemsize = grid->kcells*sizeof(double);

    // Allocate the buffer arrays at GPU.
    for(fieldmap::const_iterator it=fields->mp.begin(); it!=fields->mp.end(); ++it)
      cudaMalloc(&bufferprofs_g[it->first], nmemsize);
    for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
      cudaMalloc(&bufferprofs_g[it->first], nmemsize);

    // Copy buffers to GPU.
    for(fieldmap::const_iterator it=fields->mp.begin(); it!=fields->mp.end(); ++it)
      cudaMemcpy(bufferprofs_g[it->first], bufferprofs[it->first], nmemsize, cudaMemcpyHostToDevice);
    for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
      cudaMemcpy(bufferprofs_g[it->first], bufferprofs[it->first], nmemsize, cudaMemcpyHostToDevice);
  }
  return 0;
}

/*
#ifdef USECUDA
int cbuffer::exec()
{
  if(swbuffer == "1")
  {
    const int blocki = 128;
    const int blockj = 2;
    const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
    const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);
    const int gridk  = grid->kmax - bufferkstart + 1;

    dim3 gridGPU (gridi, gridj, gridk);
    dim3 blockGPU(blocki, blockj, 1);

    const int offs = grid->memoffset;
    double zsizebufi = 1./(grid->zsize-zstart);

    buffer_buffer<<<gridGPU, blockGPU>>>(&fields->mt["u"]->data_g[offs], &fields->mp["u"]->data_g[offs],
                                         bufferprofs_g["u"], grid->z_g, 
                                         zstart, zsizebufi, sigma, beta, 
                                         grid->istart, grid->jstart, bufferkstart,
                                         grid->iend,   grid->jend, grid->kend,
                                         grid->icellsp, grid->ijcellsp);
    
    buffer_buffer<<<gridGPU, blockGPU>>>(&fields->mt["v"]->data_g[offs], &fields->mp["v"]->data_g[offs],
                                         bufferprofs_g["v"], grid->z_g, 
                                         zstart, zsizebufi, sigma, beta, 
                                         grid->istart, grid->jstart, bufferkstart,
                                         grid->iend,   grid->jend, grid->kend,
                                         grid->icellsp, grid->ijcellsp);

    buffer_buffer<<<gridGPU, blockGPU>>>(&fields->mt["w"]->data_g[offs], &fields->mp["w"]->data_g[offs],
                                         bufferprofs_g["w"], grid->zh_g, 
                                         zstart, zsizebufi, sigma, beta, 
                                         grid->istart, grid->jstart, bufferkstart,
                                         grid->iend,   grid->jend, grid->kend,
                                         grid->icellsp, grid->ijcellsp);

    for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
      buffer_buffer<<<gridGPU, blockGPU>>>(&fields->st[it->first]->data_g[offs], &it->second->data_g[offs],
                                           bufferprofs_g[it->first], grid->z_g, 
                                           zstart, zsizebufi, sigma, beta, 
                                           grid->istart, grid->jstart, bufferkstart,
                                           grid->iend,   grid->jend, grid->kend,
                                           grid->icellsp, grid->ijcellsp);
  }

  return 0;
}
#endif
*/

