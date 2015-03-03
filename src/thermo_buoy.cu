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
#include "grid.h"
#include "fields.h"
#include "thermo_buoy.h"
#include "fd.h"
#include "tools.h"

namespace ThermoBuoy_g
{
  __global__ void calcBuoyancyTend_2nd(double * __restrict__ wt, double * __restrict__ b, 
                                       int istart, int jstart, int kstart,
                                       int iend,   int jend,   int kend,
                                       int jj, int kk)
  {
    int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
    int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
    int k = blockIdx.z + kstart; 
  
    if(i < iend && j < jend && k < kend)
    {
      int ijk = i + j*jj + k*kk;
      wt[ijk] += fd::o2::interp2(b[ijk-kk], b[ijk]);
    }
  }
  
  __global__ void calcBuoyancyTend_4th(double * __restrict__ wt, double * __restrict__ b, 
                                       int istart, int jstart, int kstart,
                                       int iend,   int jend,   int kend,
                                       int jj, int kk)
  {
    const int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
    const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
    const int k = blockIdx.z + kstart;
  
    const int kk1 = 1*kk;
    const int kk2 = 2*kk;
  
    if(i < iend && j < jend && k < kend)
    {
      const int ijk = i + j*jj + k*kk;
      wt[ijk] += fd::o4::ci0*b[ijk-kk2] + fd::o4::ci1*b[ijk-kk1] + fd::o4::ci2*b[ijk] + fd::o4::ci3*b[ijk+kk1];
    }
  }
} // end namespace

#ifdef USECUDA
void ThermoBuoy::exec()
{
  const int blocki = grid->iThreadBlock;
  const int blockj = grid->jThreadBlock;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax-1);
  dim3 blockGPU(blocki, blockj, 1);
  
  const int offs = grid->memoffset;

  if(grid->swspatialorder== "2")
  {
    ThermoBuoy_g::calcBuoyancyTend_2nd<<<gridGPU, blockGPU>>>(&fields->wt->data_g[offs], &fields->sp["b"]->data_g[offs], 
                                                              grid->istart, grid->jstart, grid->kstart+1,
                                                              grid->iend,   grid->jend, grid->kend,
                                                              grid->icellsp, grid->ijcellsp);
    cudaCheckError();
  }
  else if(grid->swspatialorder== "4")
  {
    ThermoBuoy_g::calcBuoyancyTend_4th<<<gridGPU, blockGPU>>>(&fields->wt->data_g[offs], &fields->sp["b"]->data_g[offs], 
                                                              grid->istart, grid->jstart, grid->kstart+1,
                                                              grid->iend,   grid->jend, grid->kend,
                                                              grid->icellsp, grid->ijcellsp);
    cudaCheckError();
  }
}
#endif
