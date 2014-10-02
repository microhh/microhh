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

#include "grid.h"
#include "fields.h"
#include "diff_2.h"
#include "defines.h"
#include "constants.h"
#include "tools.h"

__global__ void diff_2_diffc(double * __restrict__ const at, const double * __restrict__ const a,
                             const double * __restrict__ const dzi, const double * __restrict__ const dzhi,
                             const double dxidxi, const double dyidyi, const double visc,
                             const int jj, const int kk,
                             const int istart, const int jstart, const int kstart,
                             const int iend, const int jend, const int kend)
{
  const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  const int k = blockIdx.z + kstart;

  if(i < iend && j < jend && k < kend)
  {
    const int ijk = i + j*jj + k*kk;
    const int ii = 1;

    at[ijk] += visc * (
          + (  (a[ijk+ii] - a[ijk   ]) 
             - (a[ijk   ] - a[ijk-ii]) ) * dxidxi 
          + (  (a[ijk+jj] - a[ijk   ]) 
             - (a[ijk   ] - a[ijk-jj]) ) * dyidyi
          + (  (a[ijk+kk] - a[ijk   ]) * dzhi[k+1]
             - (a[ijk   ] - a[ijk-kk]) * dzhi[k]   ) * dzi[k] );
  }
}

__global__ void diff_2_diffw(double * __restrict__ const at, const double * __restrict__ const a,
                             const double * __restrict__ const dzi, const double * __restrict__ const dzhi,
                             const double dxidxi, const double dyidyi, const double visc,
                             const int jj, const int kk,
                             const int istart, const int jstart, const int kstart,
                             const int iend, const int jend, const int kend)
{
  const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  const int k = blockIdx.z + kstart;

  if(i < iend && j < jend && k > kstart && k < kend)
  {
    const int ijk = i + j*jj + k*kk;
    const int ii = 1;

    at[ijk] += visc * (
          + (  (a[ijk+ii] - a[ijk   ])
             - (a[ijk   ] - a[ijk-ii]) ) * dxidxi
          + (  (a[ijk+jj] - a[ijk   ])
             - (a[ijk   ] - a[ijk-jj]) ) * dyidyi
          + (  (a[ijk+kk] - a[ijk   ]) * dzi[k]
             - (a[ijk   ] - a[ijk-kk]) * dzi[k-1] ) * dzhi[k] );
  }
}

#ifdef USECUDA
int cdiff_2::exec()
{
  const int blocki = cuda::blockSizeI;
  const int blockj = cuda::blockSizeJ;
  const int gridi = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxidxi = 1./(grid->dx*grid->dx);
  const double dyidyi = 1./(grid->dy*grid->dy);

  const int offs = grid->memoffset;

  diff_2_diffc<<<gridGPU, blockGPU>>>(&fields->ut->data_g[offs], &fields->u->data_g[offs],
                                      grid->dzi_g, grid->dzhi_g,
                                      dxidxi, dyidyi, fields->visc,
                                      grid->icellsp, grid->ijcellsp,
                                      grid->istart, grid->jstart, grid->kstart,
                                      grid->iend, grid->jend, grid->kend);
  cudaCheckError();

  diff_2_diffc<<<gridGPU, blockGPU>>>(&fields->vt->data_g[offs], &fields->v->data_g[offs],
                                      grid->dzi_g, grid->dzhi_g,
                                      dxidxi, dyidyi, fields->visc,
                                      grid->icellsp, grid->ijcellsp,
                                      grid->istart, grid->jstart, grid->kstart,
                                      grid->iend, grid->jend, grid->kend);
  cudaCheckError();

  diff_2_diffw<<<gridGPU, blockGPU>>>(&fields->wt->data_g[offs], &fields->w->data_g[offs],
                                      grid->dzi_g, grid->dzhi_g,
                                      dxidxi, dyidyi, fields->visc,
                                      grid->icellsp, grid->ijcellsp,
                                      grid->istart, grid->jstart, grid->kstart,
                                      grid->iend, grid->jend, grid->kend);
  cudaCheckError();


  for(fieldmap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    diff_2_diffc<<<gridGPU, blockGPU>>>(&it->second->data_g[offs], &fields->sp[it->first]->data_g[offs],
                                        grid->dzi_g, grid->dzhi_g,
                                        dxidxi, dyidyi, fields->sp[it->first]->visc,
                                        grid->icellsp, grid->ijcellsp,
                                        grid->istart, grid->jstart, grid->kstart,
                                        grid->iend, grid->jend, grid->kend);
  cudaCheckError();

  return 0;
}
#endif
