/*
 * MicroHH
 * Copyright (c) 2011-2013 Chiel van Heerwaarden
 * Copyright (c) 2011-2013 Thijs Heus
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
#include "master.h"
#include "diff_4.h"
#include "defines.h"
#include "model.h"

__global__ void diffc_kernel(double * __restrict__ const at, const double * __restrict__ const a,
                             const double * __restrict__ const dzi4, const double * __restrict__ const dzhi4,
                             const double dx, const double dy, const double visc,
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

    const int ii1 = 1;
    const int ii2 = 2;
    const int ii3 = 3;
    const int jj1 = 1*jj;
    const int jj2 = 2*jj;
    const int jj3 = 3*jj;
    const int kk1 = 1*kk;
    const int kk2 = 2*kk;
    const int kk3 = 3*kk;

    const double dxidxi = 1./(dx*dx);
    const double dyidyi = 1./(dy*dy);

    // bottom boundary
    if(k == kstart)
    {
      at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
      at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
      at[ijk] += visc * ( cg0*(bg0*a[ijk-kk2] + bg1*a[ijk-kk1] + bg2*a[ijk    ] + bg3*a[ijk+kk1]) * dzhi4[k-1]
                        + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzhi4[k  ]
                        + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[k+1]
                        + cg3*(cg0*a[ijk    ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzhi4[k+2] )
                        * dzi4[k];
    }
    // top boundary
    else if(k == kend-1)
    {
      at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
      at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
      at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*a[ijk    ]) * dzhi4[k-1]
                        + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzhi4[k  ]
                        + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[k+1]
                        + cg3*(tg0*a[ijk-kk1] + tg1*a[ijk    ] + tg2*a[ijk+kk1] + tg3*a[ijk+kk2]) * dzhi4[k+2] )
                        * dzi4[k];
    }
    // interior
    else
    {
      at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
      at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
      at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*a[ijk    ]) * dzhi4[k-1]
                        + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzhi4[k  ]
                        + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[k+1]
                        + cg3*(cg0*a[ijk    ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzhi4[k+2] )
                        * dzi4[k];
    }
  }
}

__global__ void diffw_kernel(double * __restrict__ const at, const double * __restrict__ const a,
                             const double * __restrict__ const dzi4, const double * __restrict__ const dzhi4,
                             const double dx, const double dy, const double visc,
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

    const int ii1 = 1;
    const int ii2 = 2;
    const int ii3 = 3;
    const int jj1 = 1*jj;
    const int jj2 = 2*jj;
    const int jj3 = 3*jj;
    const int kk1 = 1*kk;
    const int kk2 = 2*kk;
    const int kk3 = 3*kk;

    const double dxidxi = 1./(dx*dx);
    const double dyidyi = 1./(dy*dy);

    // bottom boundary
    if(k == kstart+1)
    {
      at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
      at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
      at[ijk] += visc * ( cg0*(bg0*a[ijk-kk2] + bg1*a[ijk-kk1] + bg2*a[ijk    ] + bg3*a[ijk+kk1]) * dzi4[k-1]
                        + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzi4[k  ]
                        + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzi4[k+1]
                        + cg3*(cg0*a[ijk    ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzi4[k+2] )
                        * dzhi4[k];
    }
    else if(k == kend-1)
    {
      // top boundary
      at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
      at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
      at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*a[ijk    ]) * dzi4[k-2]
                        + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzi4[k-1]
                        + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzi4[k  ]
                        + cg3*(tg0*a[ijk-kk1] + tg1*a[ijk    ] + tg2*a[ijk+kk1] + tg3*a[ijk+kk2]) * dzi4[k+1] )
                        * dzhi4[k];
    }
    else
    {
      // interior
      at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
      at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
      at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*a[ijk    ]) * dzi4[k-2]
                        + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzi4[k-1]
                        + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzi4[k  ]
                        + cg3*(cg0*a[ijk    ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzi4[k+1] )
                        * dzhi4[k];
    }
  }
}

#ifdef USECUDA
int cdiff_4::exec()
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  fields->forwardGPU();
  diffc_kernel<<<gridGPU, blockGPU>>>(fields->ut->data_g, fields->u->data_g,
                                      grid->dzi4_g, grid->dzhi4_g,
                                      grid->dx, grid->dy, fields->visc,
                                      grid->icells, grid->ijcells,
                                      grid->istart, grid->jstart, grid->kstart,
                                      grid->iend, grid->jend, grid->kend);

  diffc_kernel<<<gridGPU, blockGPU>>>(fields->vt->data_g, fields->v->data_g,
                                      grid->dzi4_g, grid->dzhi4_g,
                                      grid->dx, grid->dy, fields->visc,
                                      grid->icells, grid->ijcells,
                                      grid->istart, grid->jstart, grid->kstart,
                                      grid->iend, grid->jend, grid->kend);

  diffw_kernel<<<gridGPU, blockGPU>>>(fields->wt->data_g, fields->w->data_g,
                                      grid->dzi4_g, grid->dzhi_g,
                                      grid->dx, grid->dy, fields->visc,
                                      grid->icells, grid->ijcells,
                                      grid->istart, grid->jstart, grid->kstart,
                                      grid->iend, grid->jend, grid->kend);


  for(fieldmap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    diffc_kernel<<<gridGPU, blockGPU>>>(it->second->data_g, fields->sp[it->first]->data_g,
                                        grid->dzi4_g, grid->dzhi4_g,
                                        grid->dx, grid->dy, fields->sp[it->first]->visc,
                                        grid->icells, grid->ijcells,
                                        grid->istart, grid->jstart, grid->kstart,
                                        grid->iend, grid->jend, grid->kend);

  fields->backwardGPU();

  cudaError_t error = cudaGetLastError();
  if(error != cudaSuccess)
    printf("CUDA ERROR: %s\n", cudaGetErrorString(error));

  return 0;
}
#endif
