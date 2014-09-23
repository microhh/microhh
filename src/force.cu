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

#include "master.h"
#include "force.h"
#include "grid.h"
#include "fields.h"
#include "fd.h"

__global__ void force_flux_step1(double * const __restrict__ usum, double * const __restrict__ utsum,
                                 const double * const __restrict__ u, const double * const __restrict__ ut,
                                 const double * const __restrict__ dz,
                                 const int jj, const int kk, 
                                 const int istart, const int jstart, const int kstart,
                                 const int iend,   const int jend,   const int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    usum [ijk] = u [ijk]*dz[k];
    utsum[ijk] = ut[ijk]*dz[k];
  }
}

__global__ void force_flux_step2(double * const __restrict__ ut,
                                 const double fbody,
                                 const int jj, const int kk, 
                                 const int istart, const int jstart, const int kstart,
                                 const int iend,   const int jend,   const int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    ut[ijk] += fbody;
  }
}

__global__ void force_coriolis_2nd(double * const __restrict__ ut, double * const __restrict__ vt,
                                   double * const __restrict__ u,  double * const __restrict__ v, 
                                   double * const __restrict__ ug, double * const __restrict__ vg, 
                                   const double fc, const double ugrid, const double vgrid,
                                   const int jj, const int kk, 
                                   const int istart, const int jstart, const int kstart,
                                   const int iend,   const int jend,   const int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart;

  int ii = 1;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    ut[ijk] += fc * (0.25*(v[ijk-ii] + v[ijk] + v[ijk-ii+jj] + v[ijk+jj]) + vgrid - vg[k]);
    vt[ijk] -= fc * (0.25*(u[ijk-jj] + u[ijk] + u[ijk+ii-jj] + u[ijk+ii]) + ugrid - ug[k]);
  }
}

__global__ void force_coriolis_4th(double * const __restrict__ ut, double * const __restrict__ vt,
                                   double * const __restrict__ u,  double * const __restrict__ v, 
                                   double * const __restrict__ ug, double * const __restrict__ vg, 
                                   const double fc, const double ugrid, const double vgrid,
                                   const int jj, const int kk, 
                                   const int istart, const int jstart, const int kstart,
                                   const int iend,   const int jend,   const int kend)
{
  using namespace fd::o4;

  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart;

  int ii  = 1;
  int ii2 = 2;
  int jj2 = 2*jj;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    ut[ijk] += fc * ( ( ci0*(ci0*v[ijk-ii2-jj ] + ci1*v[ijk-ii-jj ] + ci2*v[ijk-jj    ] + ci3*v[ijk+ii-jj  ])
                      + ci1*(ci0*v[ijk-ii2    ] + ci1*v[ijk-ii    ] + ci2*v[ijk       ] + ci3*v[ijk+ii     ])
                      + ci2*(ci0*v[ijk-ii2+jj ] + ci1*v[ijk-ii+jj ] + ci2*v[ijk+jj    ] + ci3*v[ijk+ii+jj  ])
                      + ci3*(ci0*v[ijk-ii2+jj2] + ci1*v[ijk-ii+jj2] + ci2*v[ijk+jj2   ] + ci3*v[ijk+ii+jj2 ]) )
                      + vgrid - vg[k] );

    vt[ijk] -= fc * ( ( ci0*(ci0*u[ijk-ii-jj2 ] + ci1*u[ijk-jj2   ] + ci2*u[ijk+ii-jj2] + ci3*u[ijk+ii2-jj2])
                      + ci1*(ci0*u[ijk-ii-jj  ] + ci1*u[ijk-jj    ] + ci2*u[ijk+ii-jj ] + ci3*u[ijk+ii2-jj ])
                      + ci2*(ci0*u[ijk-ii     ] + ci1*u[ijk       ] + ci2*u[ijk+ii    ] + ci3*u[ijk+ii2    ])
                      + ci3*(ci0*u[ijk-ii+jj  ] + ci1*u[ijk+jj    ] + ci2*u[ijk+ii+jj ] + ci3*u[ijk+ii2+jj ]) )
                      + ugrid - ug[k]);
  }
}

__global__ void force_advecwls_2nd(double * const __restrict__ st, double * const __restrict__ s,
                                   const double * const __restrict__ wls, const double * const __restrict__ dzhi,
                                   const int istart, const int jstart, const int kstart,
                                   const int iend,   const int jend,   const int kend,
                                   const int jj, const int kk)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;

    if(wls[k] > 0.)
      st[ijk] -=  wls[k] * (s[k]-s[k-1])*dzhi[k];
    else
      st[ijk] -=  wls[k] * (s[k+1]-s[k])*dzhi[k+1];
  }
}



//int cforce::advecwls_2nd(double * const restrict st, const double * const restrict s,
//                         const double * const restrict wls, const double * const dzhi)
//{
//  int ijk,jj,kk;
//
//  jj = grid->icells;
//  kk = grid->ijcells;
//
//  // use an upwind differentiation
//  for(int k=grid->kstart; k<grid->kend; ++k)
//  {
//    if(wls[k] > 0.)
//    {
//      for(int j=grid->jstart; j<grid->jend; ++j)
//        for(int i=grid->istart; i<grid->iend; ++i)
//        {
//          ijk = i + j*jj + k*kk;
//          st[ijk] -=  wls[k] * (s[k]-s[k-1])*dzhi[k];
//        }
//    }
//    else
//    {
//      for(int j=grid->jstart; j<grid->jend; ++j)
//        for(int i=grid->istart; i<grid->iend; ++i)
//        {
//          ijk = i + j*jj + k*kk;
//          st[ijk] -=  wls[k] * (s[k+1]-s[k])*dzhi[k+1];
//        }
//    }
//  }


int cforce::prepareDevice()
{
  const int nmemsize = grid->kcells*sizeof(double);

  if(swlspres == "geo")
  {
    cudaMalloc(&ug_g, nmemsize);
    cudaMalloc(&vg_g, nmemsize);

    cudaMemcpy(ug_g, ug, nmemsize, cudaMemcpyHostToDevice);
    cudaMemcpy(vg_g, vg, nmemsize, cudaMemcpyHostToDevice);
  }

  //if(swls == "1")
  //{
  //  for(std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
  //    lsprofs[*it] = new double[grid->kcells];
  //}

  if(swwls == "1")
  {
    cudaMalloc(&wls_g, nmemsize);
    cudaMemcpy(wls_g, wls, nmemsize, cudaMemcpyHostToDevice);
  }

  return 0;
}

int cforce::clearDevice()
{
  if(swlspres == "geo")
  {
    cudaFree(ug_g);
    cudaFree(vg_g);
  }

  //if(swls == "1")
  //{
  //}

  if(swwls == "1")
    cudaFree(wls_g);

  return 0; 
}

#ifdef USECUDA
int cforce::exec(double dt)
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kcells);
  dim3 blockGPU(blocki, blockj, 1);

  const int offs = grid->memoffset;

  if(swlspres == "uflux")
  {
    force_flux_step1<<<gridGPU, blockGPU>>>(&fields->a["tmp1"]->data_g[offs], &fields->a["tmp2"]->data_g[offs],
                                            &fields->u->data_g[offs], &fields->ut->data_g[offs],
                                            grid->dz_g,
                                            grid->icellsp, grid->ijcellsp,
                                            grid->istart,  grid->jstart, grid->kstart,
                                            grid->iend,    grid->jend,   grid->kend);

    double uavg  = grid->getsum_g(&fields->a["tmp1"]->data_g[offs], fields->a["tmp3"]->data_g); 
    double utavg = grid->getsum_g(&fields->a["tmp2"]->data_g[offs], fields->a["tmp3"]->data_g); 

    uavg  = uavg  / (grid->itot*grid->jtot*grid->zsize);
    utavg = utavg / (grid->itot*grid->jtot*grid->zsize);

    double fbody = (uflux - uavg - grid->utrans) / dt - utavg;

    force_flux_step2<<<gridGPU, blockGPU>>>(&fields->ut->data_g[offs],
                                            fbody,
                                            grid->icellsp, grid->ijcellsp,
                                            grid->istart,  grid->jstart, grid->kstart,
                                            grid->iend,    grid->jend,   grid->kend);
  }
  else if(swlspres == "geo")
  {
    if(grid->swspatialorder == "2")
      force_coriolis_2nd<<<gridGPU, blockGPU>>>(&fields->ut->data_g[offs], &fields->vt->data_g[offs],
                                                &fields->u->data_g[offs],  &fields->v->data_g[offs],
                                                ug_g, vg_g, fc, grid->utrans, grid->vtrans, 
                                                grid->icellsp, grid->ijcellsp,
                                                grid->istart,  grid->jstart, grid->kstart,
                                                grid->iend,    grid->jend,   grid->kend);
    else if(grid->swspatialorder == "4")
      force_coriolis_4th<<<gridGPU, blockGPU>>>(&fields->ut->data_g[offs], &fields->vt->data_g[offs],
                                                &fields->u->data_g[offs],  &fields->v->data_g[offs],
                                                ug_g, vg_g, fc, grid->utrans, grid->vtrans, 
                                                grid->icellsp, grid->ijcellsp,
                                                grid->istart,  grid->jstart, grid->kstart,
                                                grid->iend,    grid->jend,   grid->kend);
  }

  //if(swls == "1")
  //{
  //  for(std::vector<std::string>::const_iterator it=lslist.begin(); it!=lslist.end(); ++it)
  //    lssource(fields->st[*it]->data, lsprofs[*it]);
  //}

  if(swwls == "1")
  {
    for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
      force_advecwls_2nd<<<gridGPU, blockGPU>>>(&it->second->data_g[offs], fields->s[it->first]->datamean_g, wls_g, grid->dzhi_g,
                                                grid->istart,  grid->jstart, grid->kstart,
                                                grid->iend,    grid->jend,   grid->kend,
                                                grid->icellsp, grid->ijcellsp);
  }

  return 0;
}
#endif

