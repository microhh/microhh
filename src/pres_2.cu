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
#include <algorithm>
#include <fftw3.h>
#include <cufft.h>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "pres.h"
#include "pres_2.h"
#include "defines.h"
#include "model.h"
#include "tools.h"
#include "constants.h"

namespace Pres_2_g
{
  __global__ void presin(double * __restrict__ p,
                         double * __restrict__ u ,  double * __restrict__ v ,     double * __restrict__ w ,
                         double * __restrict__ ut,  double * __restrict__ vt,     double * __restrict__ wt,
                         double * __restrict__ dzi, double * __restrict__ rhoref, double * __restrict__ rhorefh,
                         double dxi, double dyi, double dti,
                         const int jj, const int kk,
                         const int jjp, const int kkp,
                         const int imax, const int jmax, const int kmax,
                         const int igc, const int jgc, const int kgc)
  {
    const int ii = 1;
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    const int j = blockIdx.y*blockDim.y + threadIdx.y;
    const int k = blockIdx.z;

    if(i < imax && j < jmax && k < kmax)
    {
      const int ijkp = i + j*jjp + k*kkp;
      const int ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;

      p[ijkp] = rhoref [k+kgc]   * ( (ut[ijk+ii] + u[ijk+ii] * dti) - (ut[ijk] + u[ijk] * dti) ) * dxi
              + rhoref [k+kgc]   * ( (vt[ijk+jj] + v[ijk+jj] * dti) - (vt[ijk] + v[ijk] * dti) ) * dyi
            + ( rhorefh[k+kgc+1] * (  wt[ijk+kk] + w[ijk+kk] * dti)
              - rhorefh[k+kgc  ] * (  wt[ijk   ] + w[ijk   ] * dti) ) * dzi[k+kgc];
    }
  }

  __global__ void presout(double * __restrict__ ut, double * __restrict__ vt, double * __restrict__ wt,
                          double * __restrict__ p,
                          double * __restrict__ dzhi, const double dxi, const double dyi,
                          const int jj, const int kk,
                          const int istart, const int jstart, const int kstart,
                          const int iend, const int jend, const int kend)
  {
    const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
    const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
    const int k = blockIdx.z + kstart;

    const int ii = 1;

    if(i < iend && j < jend && k < kend)
    {
      int ijk = i + j*jj + k*kk;
      ut[ijk] -= (p[ijk] - p[ijk-ii]) * dxi;
      vt[ijk] -= (p[ijk] - p[ijk-jj]) * dyi;
      wt[ijk] -= (p[ijk] - p[ijk-kk]) * dzhi[k];
    }
  }

  __global__ void solveout(double * __restrict__ p, double * __restrict__ work3d,
                           const int jj, const int kk,
                           const int jjp, const int kkp,
                           const int istart, const int jstart, const int kstart,
                           const int imax, const int jmax, const int kmax)
  {
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    const int j = blockIdx.y*blockDim.y + threadIdx.y;
    const int k = blockIdx.z;

    if(i < imax && j < jmax && k < kmax)
    {
      const int ijk  = i + j*jj + k*kk;
      const int ijkp = i+istart + (j+jstart)*jjp + (k+kstart)*kkp;

      p[ijkp] = work3d[ijk];

      if(k == 0)
        p[ijkp-kkp] = p[ijkp];
    }
  }

  __global__ void solvein(double * __restrict__ p,
                          double * __restrict__ work3d, double * __restrict__ b,
                          double * __restrict__ a, double * __restrict__ c,
                          double * __restrict__ dz, double * __restrict__ rhoref,
                          double * __restrict__ bmati, double * __restrict__ bmatj,
                          const int jj, const int kk,
                          const int imax, const int jmax, const int kmax,
                          const int kstart)
  {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;
    int k = blockIdx.z;

    if(i < imax && j < jmax && k < kmax)
    {
      int ijk = i + j*jj + k*kk;

      // CvH this needs to be taken into account in case of an MPI run
      // iindex = mpi->mpicoordy * iblock + i;
      // jindex = mpi->mpicoordx * jblock + j;
      // b[ijk] = dz[k+kgc]*dz[k+kgc] * (bmati[iindex]+bmatj[jindex]) - (a[k]+c[k]);
      //  if(iindex == 0 && jindex == 0)

      b[ijk] = dz[k+kstart]*dz[k+kstart] * rhoref[k+kstart]*(bmati[i]+bmatj[j]) - (a[k]+c[k]);
      p[ijk] = dz[k+kstart]*dz[k+kstart] * p[ijk];

      if(k == 0)
      {
        // substitute BC's
        // ijk = i + j*jj;
        b[ijk] += a[0];
      }
      else if(k == kmax-1)
      {
        // for wave number 0, which contains average, set pressure at top to zero
        if(i == 0 && j == 0)
          b[ijk] -= c[k];
        // set dp/dz at top to zero
        else
          b[ijk] += c[k];
      }
    }
  }

  __global__ void tdma(double * __restrict__ a, double * __restrict__ b, double * __restrict__ c,
                       double * __restrict__ p, double * __restrict__ work3d,
                       const int jj, const int kk,
                       const int imax, const int jmax, const int kmax)
  {
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    const int j = blockIdx.y*blockDim.y + threadIdx.y;

    if(i < imax && j < jmax)
    {
      const int ij = i + j*jj;
      int k,ijk;

      double work2d = b[ij];
      p[ij] /= work2d;

      for(k=1; k<kmax; k++)
      {
        ijk = ij + k*kk;
        work3d[ijk] = c[k-1] / work2d;
        work2d = b[ijk] - a[k]*work3d[ijk];
        p[ijk] -= a[k]*p[ijk-kk];
        p[ijk] /= work2d;
      }

      for(k=kmax-2; k>=0; k--)
      {
        ijk = ij + k*kk;
        p[ijk] -= work3d[ijk+kk]*p[ijk+kk];
      }
    }
  }

  __global__ void calcdivergence(double * __restrict__ u, double * __restrict__ v, double * __restrict__ w,
                                 double * __restrict__ div, double * __restrict__ dzi,
                                 double * __restrict__ rhoref, double * __restrict__ rhorefh,
                                 double dxi, double dyi,
                                 int jj, int kk, int istart, int jstart, int kstart,
                                 int iend, int jend, int kend)
  {
    int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
    int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
    int k = blockIdx.z + kstart;
    int ii = 1;

    if(i < iend && j < jend && k < kend)
    {
      int ijk = i + j*jj + k*kk;
      div[ijk] = rhoref[k]*((u[ijk+ii]-u[ijk])*dxi + (v[ijk+jj]-v[ijk])*dyi)
               + (rhorefh[k+1]*w[ijk+kk]-rhorefh[k]*w[ijk])*dzi[k];
    }
  }
} // End namespace.

void Pres_2::prepare_device()
{
  const int kmemsize = grid->kmax*sizeof(double);
  const int imemsize = grid->itot*sizeof(double);
  const int jmemsize = grid->jtot*sizeof(double);

  const int ijmemsize = grid->imax*grid->jmax*sizeof(double);

  cudaSafeCall(cudaMalloc((void**)&bmati_g, imemsize  ));
  cudaSafeCall(cudaMalloc((void**)&bmatj_g, jmemsize  ));
  cudaSafeCall(cudaMalloc((void**)&a_g, kmemsize      ));
  cudaSafeCall(cudaMalloc((void**)&c_g, kmemsize      ));
  cudaSafeCall(cudaMalloc((void**)&work2d_g, ijmemsize));

  cudaSafeCall(cudaMemcpy(bmati_g, bmati, imemsize, cudaMemcpyHostToDevice   ));
  cudaSafeCall(cudaMemcpy(bmatj_g, bmatj, jmemsize, cudaMemcpyHostToDevice   ));
  cudaSafeCall(cudaMemcpy(a_g, a, kmemsize, cudaMemcpyHostToDevice           ));
  cudaSafeCall(cudaMemcpy(c_g, c, kmemsize, cudaMemcpyHostToDevice           ));
  cudaSafeCall(cudaMemcpy(work2d_g, work2d, ijmemsize, cudaMemcpyHostToDevice));

  makeCufftPlan();
}

void Pres_2::clear_device()
{
  cudaSafeCall(cudaFree(bmati_g ));
  cudaSafeCall(cudaFree(bmatj_g ));
  cudaSafeCall(cudaFree(a_g     ));
  cudaSafeCall(cudaFree(c_g     ));
  cudaSafeCall(cudaFree(work2d_g));
}

#ifdef USECUDA
void Pres_2::exec(double dt)
{
  const int blocki = grid->iThreadBlock;
  const int blockj = grid->jThreadBlock;
  const int gridi = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj = grid->jmax/blockj + (grid->jmax%blockj > 0);

  // 3D grid
  dim3 gridGPU (gridi,  gridj,  grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  // 2D grid
  dim3 grid2dGPU (gridi,  gridj);
  dim3 block2dGPU(blocki, blockj);

  const int offs = grid->memoffset;

  // calculate the cyclic BCs first
  grid->boundary_cyclic_g(&fields->ut->data_g[offs]);
  grid->boundary_cyclic_g(&fields->vt->data_g[offs]);
  grid->boundary_cyclic_g(&fields->wt->data_g[offs]);

  Pres_2_g::presin<<<gridGPU, blockGPU>>>(fields->sd["p"]->data_g,
                                         &fields->u->data_g[offs],  &fields->v->data_g[offs],  &fields->w->data_g[offs],
                                         &fields->ut->data_g[offs], &fields->vt->data_g[offs], &fields->wt->data_g[offs],
                                         grid->dzi_g, fields->rhoref_g, fields->rhorefh_g,
                                         1./grid->dx, 1./grid->dy, 1./dt,
                                         grid->icellsp, grid->ijcellsp, grid->imax, grid->imax*grid->jmax,
                                         grid->imax, grid->jmax, grid->kmax,
                                         grid->igc, grid->jgc, grid->kgc);
  cudaCheckError();

  fft_forward(fields->sd["p"]->data_g, fields->atmp["tmp1"]->data_g, fields->atmp["tmp2"]->data_g);

  Pres_2_g::solvein<<<gridGPU, blockGPU>>>(fields->sd["p"]->data_g,
                                          fields->atmp["tmp1"]->data_g, fields->atmp["tmp2"]->data_g,
                                          a_g, c_g,
                                          grid->dz_g, fields->rhoref_g, bmati_g, bmatj_g,
                                          grid->imax, grid->imax*grid->jmax,
                                          grid->imax, grid->jmax, grid->kmax,
                                          grid->kstart);
  cudaCheckError();

  Pres_2_g::tdma<<<grid2dGPU, block2dGPU>>>(a_g, fields->atmp["tmp2"]->data_g, c_g,
                                           fields->sd["p"]->data_g, fields->atmp["tmp1"]->data_g,
                                           grid->imax, grid->imax*grid->jmax,
                                           grid->imax, grid->jmax, grid->kmax);
  cudaCheckError();

  fft_backward(fields->sd["p"]->data_g, fields->atmp["tmp1"]->data_g, fields->atmp["tmp2"]->data_g);

  cudaSafeCall(cudaMemcpy(fields->atmp["tmp1"]->data_g, fields->sd["p"]->data_g, grid->ncellsp*sizeof(double), cudaMemcpyDeviceToDevice));

  Pres_2_g::solveout<<<gridGPU, blockGPU>>>(&fields->sd["p"]->data_g[offs], fields->atmp["tmp1"]->data_g,
                                           grid->imax, grid->imax*grid->jmax,
                                           grid->icellsp, grid->ijcellsp,
                                           grid->istart, grid->jstart, grid->kstart,
                                           grid->imax, grid->jmax, grid->kmax);
  cudaCheckError();

  grid->boundary_cyclic_g(&fields->sd["p"]->data_g[offs]);

  Pres_2_g::presout<<<gridGPU, blockGPU>>>(&fields->ut->data_g[offs], &fields->vt->data_g[offs], &fields->wt->data_g[offs],
                                          &fields->sd["p"]->data_g[offs],
                                          grid->dzhi_g, 1./grid->dx, 1./grid->dy,
                                          grid->icellsp, grid->ijcellsp,
                                          grid->istart, grid->jstart, grid->kstart,
                                          grid->iend, grid->jend, grid->kend);
  cudaCheckError();
}
#endif

#ifdef USECUDA
double Pres_2::checkDivergence()
{
  const int blocki = grid->iThreadBlock;
  const int blockj = grid->jThreadBlock;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  double divmax = 0;

  dim3 gridGPU (gridi, gridj, grid->kcells);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxi = 1./grid->dx;
  const double dyi = 1./grid->dy;

  const int offs = grid->memoffset;

  Pres_2_g::calcdivergence<<<gridGPU, blockGPU>>>(&fields->u->data_g[offs], &fields->v->data_g[offs], &fields->w->data_g[offs],
                                                 &fields->atmp["tmp1"]->data_g[offs], grid->dzi_g,
                                                 fields->rhoref_g, fields->rhorefh_g, dxi, dyi,
                                                 grid->icellsp, grid->ijcellsp,
                                                 grid->istart,  grid->jstart, grid->kstart,
                                                 grid->iend,    grid->jend,   grid->kend);
  cudaCheckError();

  divmax = grid->get_max_g(&fields->atmp["tmp1"]->data_g[offs], fields->atmp["tmp2"]->data_g);
  grid->get_max(&divmax);

  return divmax;
}
#endif
