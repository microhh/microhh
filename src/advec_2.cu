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

#include "advec_2.h"
#include "grid.h"
#include "fields.h"
#include "tools.h"
#include "constants.h"
#include "tools.h"

__device__ double interp2(double a, double b)
{
  return 0.5*(a + b);
}

//__global__ void advec_2_advecu(double * __restrict__ ut, double * __restrict__ u, 
//                               double * __restrict__ v, double * __restrict__ w,
//                               double * __restrict__ dzi, double dxi, double dyi, 
//                               int jj, int kk, int istart, int jstart, int kstart,
//                               int iend,   int jend,   int kend)
//{
//  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
//  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
//  int k = blockIdx.z + kstart;
//  int ii = 1;
//
//  if(i < iend && j < jend && k < kend)
//  {
//    int ijk = i + j*jj + k*kk;
//    ut[ijk] += 
//          - (  interp2(u[ijk   ], u[ijk+ii]) * interp2(u[ijk   ], u[ijk+ii])
//             - interp2(u[ijk-ii], u[ijk   ]) * interp2(u[ijk-ii], u[ijk   ]) ) * dxi
//
//          - (  interp2(v[ijk-ii+jj], v[ijk+jj]) * interp2(u[ijk   ], u[ijk+jj])
//             - interp2(v[ijk-ii   ], v[ijk   ]) * interp2(u[ijk-jj], u[ijk   ]) ) * dyi 
//
//          - (  interp2(w[ijk-ii+kk], w[ijk+kk]) * interp2(u[ijk   ], u[ijk+kk])
//             - interp2(w[ijk-ii   ], w[ijk   ]) * interp2(u[ijk-kk], u[ijk   ]) ) * dzi[k];
//  }
//}
//
//__global__ void advec_2_advecv(double * __restrict__ vt, double * __restrict__ u, 
//                               double * __restrict__ v, double * __restrict__ w,
//                               double * __restrict__ dzi, double dxi, double dyi, 
//                               int jj, int kk, int istart, int jstart, int kstart,
//                               int iend,   int jend,   int kend)
//{
//  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
//  int j = blockIdx.y*blockDisdm.y + threadIdx.y + jstart;
//  int k = blockIdx.z + kstart;
//  int ii = 1;
//
//  if(i < iend && j < jend && k < kend)
//  {
//    int ijk = i + j*jj + k*kk;
//    vt[ijk] += 
//          - (  interp2(u[ijk+ii-jj], u[ijk+ii]) * interp2(v[ijk   ], v[ijk+ii])
//             - interp2(u[ijk   -jj], u[ijk   ]) * interp2(v[ijk-ii], v[ijk   ]) ) * dxi
//
//          - (  interp2(v[ijk   ], v[ijk+jj]) * interp2(v[ijk   ], v[ijk+jj])
//             - interp2(v[ijk-jj], v[ijk   ]) * interp2(v[ijk-jj], v[ijk   ]) ) * dyi
//
//          - (  interp2(w[ijk-jj+kk], w[ijk+kk]) * interp2(v[ijk   ], v[ijk+kk])
//             - interp2(w[ijk-jj   ], w[ijk   ]) * interp2(v[ijk-kk], v[ijk   ]) ) * dzi[k];
//  }
//}
//
//__global__ void advec_2_advecw(double * __restrict__ wt, double * __restrict__ u, 
//                               double * __restrict__ v, double * __restrict__ w,
//                               double * __restrict__ dzhi, double dxi, double dyi, 
//                               int jj, int kk, int istart, int jstart, int kstart,
//                               int iend,   int jend,   int kend)
//{
//  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
//  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
//  int k = blockIdx.z + kstart + 1;
//  int ii = 1;
//
//  if(i < iend && j < jend && k < kend)
//  {
//    int ijk = i + j*jj + k*kk;
//    wt[ijk] += 
//          - (  interp2(u[ijk+ii-kk], u[ijk+ii]) * interp2(w[ijk   ], w[ijk+ii])
//             - interp2(u[ijk   -kk], u[ijk   ]) * interp2(w[ijk-ii], w[ijk   ]) ) * dxi
//
//          - (  interp2(v[ijk+jj-kk], v[ijk+jj]) * interp2(w[ijk   ], w[ijk+jj])
//             - interp2(v[ijk   -kk], v[ijk   ]) * interp2(w[ijk-jj], w[ijk   ]) ) * dyi
//
//          - (  interp2(w[ijk   ], w[ijk+kk]) * interp2(w[ijk   ], w[ijk+kk])
//             - interp2(w[ijk-kk], w[ijk   ]) * interp2(w[ijk-kk], w[ijk   ]) ) * dzhi[k];
//  }
//}

__global__ void advec_2_advecuvw(double * __restrict__ ut, double * __restrict__ vt, double * __restrict__ wt, 
                                 double * __restrict__ u,  double * __restrict__ v,  double * __restrict__ w,
                                 double * __restrict__ rhoref, double * __restrict__ rhorefh,
                                 double * __restrict__ dzi, double * __restrict__ dzhi, double dxi, double dyi, 
                                 int jj, int kk, int istart, int jstart, int kstart,
                                 int iend,   int jend,   int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart;
  int ii = 1;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    ut[ijk] += 
          - (  interp2(u[ijk   ], u[ijk+ii]) * interp2(u[ijk   ], u[ijk+ii])
             - interp2(u[ijk-ii], u[ijk   ]) * interp2(u[ijk-ii], u[ijk   ]) ) * dxi

          - (  interp2(v[ijk-ii+jj], v[ijk+jj]) * interp2(u[ijk   ], u[ijk+jj])
             - interp2(v[ijk-ii   ], v[ijk   ]) * interp2(u[ijk-jj], u[ijk   ]) ) * dyi 

          - (  rhorefh[k+1] * interp2(w[ijk-ii+kk], w[ijk+kk]) * interp2(u[ijk   ], u[ijk+kk])
             - rhorefh[k  ] * interp2(w[ijk-ii   ], w[ijk   ]) * interp2(u[ijk-kk], u[ijk   ]) ) / rhoref[k] * dzi[k];

    vt[ijk] += 
          - (  interp2(u[ijk+ii-jj], u[ijk+ii]) * interp2(v[ijk   ], v[ijk+ii])
             - interp2(u[ijk   -jj], u[ijk   ]) * interp2(v[ijk-ii], v[ijk   ]) ) * dxi

          - (  interp2(v[ijk   ], v[ijk+jj]) * interp2(v[ijk   ], v[ijk+jj])
             - interp2(v[ijk-jj], v[ijk   ]) * interp2(v[ijk-jj], v[ijk   ]) ) * dyi

          - (  rhorefh[k+1] * interp2(w[ijk-jj+kk], w[ijk+kk]) * interp2(v[ijk   ], v[ijk+kk])
             - rhorefh[k  ] * interp2(w[ijk-jj   ], w[ijk   ]) * interp2(v[ijk-kk], v[ijk   ]) ) / rhoref[k] * dzi[k];

    if(k>kstart)
    {
      wt[ijk] += 
            - (  interp2(u[ijk+ii-kk], u[ijk+ii]) * interp2(w[ijk   ], w[ijk+ii])
               - interp2(u[ijk   -kk], u[ijk   ]) * interp2(w[ijk-ii], w[ijk   ]) ) * dxi

            - (  interp2(v[ijk+jj-kk], v[ijk+jj]) * interp2(w[ijk   ], w[ijk+jj])
               - interp2(v[ijk   -kk], v[ijk   ]) * interp2(w[ijk-jj], w[ijk   ]) ) * dyi

            - (  rhoref[k  ] * interp2(w[ijk   ], w[ijk+kk]) * interp2(w[ijk   ], w[ijk+kk])
               - rhoref[k-1] * interp2(w[ijk-kk], w[ijk   ]) * interp2(w[ijk-kk], w[ijk   ]) ) / rhorefh[k] * dzhi[k];
    }
  }
}

__global__ void advec_2_advecs(double * __restrict__ st, double * __restrict__ s, 
                               double * __restrict__ u, double * __restrict__ v, double * __restrict__ w,
                               double * __restrict__ rhoref, double * __restrict__ rhorefh,
                               double * __restrict__ dzi, double dxi, double dyi, 
                               int jj, int kk, int istart, int jstart, int kstart,
                               int iend,   int jend,   int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart;
  int ii = 1;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    st[ijk] += 
          - (  u[ijk+ii] * interp2(s[ijk   ], s[ijk+ii])
             - u[ijk   ] * interp2(s[ijk-ii], s[ijk   ]) ) * dxi

          - (  v[ijk+jj] * interp2(s[ijk   ], s[ijk+jj])
             - v[ijk   ] * interp2(s[ijk-jj], s[ijk   ]) ) * dyi 

          - (  rhorefh[k+1] * w[ijk+kk] * interp2(s[ijk   ], s[ijk+kk])
             - rhorefh[k  ] * w[ijk   ] * interp2(s[ijk-kk], s[ijk   ]) ) / rhoref[k] * dzi[k];
  }
}

__global__ void advec_2_calccfl(double * __restrict__ u, double * __restrict__ v, double * __restrict__ w, 
                                double * __restrict__ cfl,
                                double * __restrict__ dzi, double dxi, double dyi,
                                int jj, int kk,
                                int istart, int jstart, int kstart,
                                int iend, int jend, int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
  int k = blockIdx.z + kstart; 
  int ii = 1;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    cfl[ijk] = std::abs(interp2(u[ijk], u[ijk+ii]))*dxi + 
               std::abs(interp2(v[ijk], v[ijk+jj]))*dyi + 
               std::abs(interp2(w[ijk], w[ijk+kk]))*dzi[k];
  }
}

#ifdef USECUDA
void Advec_2::exec()
{
  const int blocki = cuda::blockSizeI;
  const int blockj = cuda::blockSizeJ;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxi = 1./grid->dx;
  const double dyi = 1./grid->dy;
  
  const int offs = grid->memoffset;

  //advec_2_advecu<<<gridGPU, blockGPU>>>(fields->ut->data_g, fields->u->data_g, fields->v->data_g, 
  //                                      fields->w->data_g, grid->dzi_g, dxi, dyi,
  //                                      grid->icells, grid->ijcells,
  //                                      grid->istart, grid->jstart, grid->kstart,
  //                                      grid->iend,   grid->jend, grid->kend);

  //advec_2_advecv<<<gridGPU, blockGPU>>>(fields->vt->data_g, fields->u->data_g, fields->v->data_g, 
  //                                      fields->w->data_g, grid->dzi_g, dxi, dyi,
  //                                      grid->icells, grid->ijcells,
  //                                      grid->istart, grid->jstart, grid->kstart,
  //                                      grid->iend,   grid->jend, grid->kend);

  //advec_2_advecw<<<gridGPU, blockGPU>>>(fields->wt->data_g, fields->u->data_g, fields->v->data_g, 
  //                                      fields->w->data_g, grid->dzhi_g, dxi, dyi,
  //                                      grid->icells, grid->ijcells,
  //                                      grid->istart, grid->jstart, grid->kstart,
  //                                      grid->iend,   grid->jend, grid->kend);

  advec_2_advecuvw<<<gridGPU, blockGPU>>>(&fields->ut->data_g[offs], &fields->vt->data_g[offs], &fields->wt->data_g[offs], 
                                          &fields->u->data_g[offs],  &fields->v->data_g[offs],  &fields->w->data_g[offs], 
                                          fields->rhoref_g, fields->rhorefh_g, 
                                          grid->dzi_g, grid->dzhi_g, 
                                          dxi, dyi, grid->icellsp, grid->ijcellsp,
                                          grid->istart, grid->jstart, grid->kstart,
                                          grid->iend,   grid->jend, grid->kend);
  cudaCheckError(); 
  
  for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    advec_2_advecs<<<gridGPU, blockGPU>>>(&it->second->data_g[offs], &fields->sp[it->first]->data_g[offs], 
                                          &fields->u->data_g[offs], &fields->v->data_g[offs], &fields->w->data_g[offs],
                                          fields->rhoref_g, fields->rhorefh_g, 
                                          grid->dzi_g, dxi, dyi,
                                          grid->icellsp, grid->ijcellsp,
                                          grid->istart, grid->jstart, grid->kstart,
                                          grid->iend,   grid->jend, grid->kend);
  cudaCheckError(); 
}
#endif

#ifdef USECUDA
double Advec_2::calccfl(double * u, double * v, double * w, double * dzi, double dt)
{
  const int blocki = cuda::blockSizeI;
  const int blockj = cuda::blockSizeJ;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);
  double cfl = 0;

  dim3 gridGPU (gridi, gridj, grid->kcells);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxi = 1./grid->dx;
  const double dyi = 1./grid->dy;

  const int offs = grid->memoffset;

  advec_2_calccfl<<<gridGPU, blockGPU>>>(&fields->u->data_g[offs], &fields->v->data_g[offs], &fields->w->data_g[offs], 
                                         &fields->atmp["tmp1"]->data_g[offs],
                                         grid->dzi_g, dxi, dyi,
                                         grid->icellsp, grid->ijcellsp,
                                         grid->istart,  grid->jstart, grid->kstart,
                                         grid->iend,    grid->jend,   grid->kend);
  cudaCheckError(); 

  cfl = grid->getmax_g(&fields->atmp["tmp1"]->data_g[offs], fields->atmp["tmp2"]->data_g); 
  grid->getmax(&cfl); 
  cfl = cfl*dt;

  return cfl;
}
#endif
