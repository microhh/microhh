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

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fftw3.h>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "pres_2.h"
#include "defines.h"
#include "model.h"

__global__ void pres_2_presin(double * __restrict__ p,
                              double * __restrict__ u ,  double * __restrict__ v , double * __restrict__ w ,
                              double * __restrict__ ut,  double * __restrict__ vt, double * __restrict__ wt,
                              double * __restrict__ dzi, double dx, double dy, double dt,
                              const int jj, const int kk,
                              const int jjp, const int kkp,
                              const int imax, const int jmax, const int kmax,
                              const int igc, const int jgc, const int kgc)
{
  const int ii = 1;
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  const int j = blockIdx.y*blockDim.y + threadIdx.y;
  const int k = blockIdx.z;

  const double dxi = 1./dx;
  const double dyi = 1./dy;

  if(i < imax && j < jmax && k < kmax)
  {
    const int ijkp = i + j*jjp + k*kkp;
    const int ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;
    p[ijkp] = ( (ut[ijk+ii] + u[ijk+ii] / dt) - (ut[ijk] + u[ijk] / dt) ) * dxi
            + ( (vt[ijk+jj] + v[ijk+jj] / dt) - (vt[ijk] + v[ijk] / dt) ) * dyi
            + ( (wt[ijk+kk] + w[ijk+kk] / dt) - (wt[ijk] + w[ijk] / dt) ) * dzi[k+kgc];
  }
}

__global__ void pres_2_presout(double * __restrict__ ut, double * __restrict__ vt, double * __restrict__ wt,
                               double * __restrict__ p,
                               double * __restrict__ dzhi, const double dx, const double dy,
                               const int jj, const int kk,
                               const int istart, const int jstart, const int kstart,
                               const int iend, const int jend, const int kend)
{
  const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  const int k = blockIdx.z + kstart;

  const int ii = 1;

  const double dxi = 1./dx;
  const double dyi = 1./dy;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    ut[ijk] -= (p[ijk] - p[ijk-ii]) * dxi;
    vt[ijk] -= (p[ijk] - p[ijk-jj]) * dyi;
    wt[ijk] -= (p[ijk] - p[ijk-kk]) * dzhi[k];
  }
}

__global__ void pres_2_solveout(double * __restrict__ p, double * __restrict__ work3d,
                                const int jj, const int kk,
                                const int jjp, const int kkp,
                                const int istart, const int jstart, const int kstart,
                                const int imax, const int jmax, const int kmax)
{
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  const int j = blockIdx.y*blockDim.y + threadIdx.y;
  const int k = blockIdx.z;
  const int ijk  = i + j*jj + k*kk;
  const int ijkp = i+istart + (j+jstart)*jjp + (k+kstart)*kkp;

  if(i < imax && j < jmax && k < kmax)
  {
    p[ijkp] = work3d[ijk];

    if(k == 0)
      p[ijkp-kkp] = p[ijkp];
  }
}

/*
__global__ void pres_solve_in(double * __restrict__ a, double * __restrict__ b, double * __restrict__ c,
                              double * __restrict__ p, double * __restrict__ work3d,
                              double * __restrict__ dz, double * __restrict__ bmati, double * __restrict__ bmatj,
                              int kmax, int jj, int kk, int kgc)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  int k = blockIdx.z;
  int ijk = i + j*jj + k*kk;

  // CvH this needs to be taken into account in case of an MPI run
  // iindex = mpi->mpicoordy * iblock + i;
  // jindex = mpi->mpicoordx * jblock + j;
  // b[ijk] = dz[k+kgc]*dz[k+kgc] * (bmati[iindex]+bmatj[jindex]) - (a[k]+c[k]);
  //  if(iindex == 0 && jindex == 0)

  b[ijk] = dz[k+kgc]*dz[k+kgc] * (bmati[i]+bmatj[j]) - (a[k]+c[k]);
  p[ijk] = dz[k+kgc]*dz[k+kgc] * p[ijk];

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

__global__ void pres_tdma(double * __restrict__ a, double * __restrict__ b, double * __restrict__ c, 
                          double * __restrict__ p, double * __restrict__ work3d,
                          int kmax, int jj, int kk)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  int ij = i + j*jj;

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

__global__ void pres_solve_out(double * __restrict__ p, double * __restrict__ work3d,
                               int jj, int kk, int jjp, int kkp, int igc, int jgc, int kgc)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  int k = blockIdx.z;
  int ijk  = i + j*jj + k*kk;
  int ijkp = i+igc + (j+jgc)*jjp + (k+kgc)*kkp;

  p[ijkp] = work3d[ijk];

  if(k == 0)
    p[ijkp-kkp] = p[ijkp];

  //for(int j=grid->jstart; j<grid->jend; j++)
  //  for(int i=grid->istart; i<grid->iend; i++)
  //      ijk = i + j*jjp + grid->kstart*kkp;
}*/

#ifdef USECUDA
int cpres_2::exec(double dt)
{
  fields->forwardGPU();

  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  // calculate the cyclic BCs first
  grid->boundary_cyclic(fields->ut->data_g);
  grid->boundary_cyclic(fields->vt->data_g);
  grid->boundary_cyclic(fields->wt->data_g);

  pres_2_presin<<<gridGPU, blockGPU>>>(fields->sd["p"]->data_g,
                                       fields->u->data_g, fields->v->data_g, fields->w->data_g,
                                       fields->ut->data_g, fields->vt->data_g, fields->wt->data_g,
                                       grid->dzi_g, grid->dx, grid->dy, dt,
                                       grid->icells, grid->ijcells, grid->imax, grid->imax*grid->jmax, 
                                       grid->imax, grid->jmax, grid->kmax,
                                       grid->igc, grid->jgc, grid->kgc);
  fields->backwardGPU();

  // solve the system
  grid->fftforward(fields->sd["p"]->data, fields->sd["tmp1"]->data,
                   grid->fftini, grid->fftouti, grid->fftinj, grid->fftoutj);

  pres_solve(fields->sd["p"]->data, fields->sd["tmp1"]->data, fields->sd["tmp2"]->data, grid->dz,
             grid->fftini, grid->fftouti, grid->fftinj, grid->fftoutj);

  grid->fftbackward(fields->sd["p"]->data, fields->sd["tmp1"]->data,
                    grid->fftini, grid->fftouti, grid->fftinj, grid->fftoutj);

  fields->forwardGPU();
  pres_2_solveout<<<gridGPU, blockGPU>>>(fields->sd["p"]->data_g, fields->sd["tmp1"]->data_g,
                                         grid->imax, grid->imax*grid->jmax,
                                         grid->icells, grid->ijcells,
                                         grid->istart, grid->jstart, grid->kstart,
                                         grid->imax, grid->jmax, grid->kmax);

  grid->boundary_cyclic(fields->sd["p"]->data_g);

  pres_2_presout<<<gridGPU, blockGPU>>>(fields->ut->data_g, fields->vt->data_g, fields->wt->data_g,
                                        fields->sd["p"]->data_g,
                                        grid->dzhi_g, grid->dx, grid->dy,
                                        grid->icells, grid->ijcells,
                                        grid->istart, grid->jstart, grid->kstart,
                                        grid->iend, grid->jend, grid->kend);
  fields->backwardGPU();

  return 0;
}
#endif

#ifdef USECUDA
int cpres_2::pres_solve(double * restrict p, double * restrict work3d, double * restrict b, double * restrict dz,
                        double * restrict fftini, double * restrict fftouti, 
                        double * restrict fftinj, double * restrict fftoutj)

{
  int i,j,k,jj,kk,ijk;
  int imax,jmax,kmax;
  int itot,jtot;
  int iblock,jblock,kblock;
  int igc,jgc,kgc;
  int iindex,jindex;

  imax   = grid->imax;
  jmax   = grid->jmax;
  kmax   = grid->kmax;
  itot   = grid->itot;
  jtot   = grid->jtot;
  iblock = grid->iblock;
  jblock = grid->jblock;
  kblock = grid->kblock;
  igc    = grid->igc;
  jgc    = grid->jgc;
  kgc    = grid->kgc;

  jj = iblock;
  kk = iblock*jblock;

  // solve the tridiagonal system
  // create vectors that go into the tridiagonal matrix solver
  for(k=0; k<kmax; k++)
    for(j=0; j<jblock; j++)
#pragma ivdep
      for(i=0; i<iblock; i++)
      {
        // swap the mpicoords, because domain is turned 90 degrees to avoid two mpi transposes
        iindex = master->mpicoordy * iblock + i;
        jindex = master->mpicoordx * jblock + j;

        ijk  = i + j*jj + k*kk;
        b[ijk] = dz[k+kgc]*dz[k+kgc] * (bmati[iindex]+bmatj[jindex]) - (a[k]+c[k]);
        p[ijk] = dz[k+kgc]*dz[k+kgc] * p[ijk];
      }

  for(j=0; j<jblock; j++)
#pragma ivdep
    for(i=0; i<iblock; i++)
    {
      iindex = master->mpicoordy * iblock + i;
      jindex = master->mpicoordx * jblock + j;

      // substitute BC's
      ijk = i + j*jj;
      b[ijk] += a[0];

      // for wave number 0, which contains average, set pressure at top to zero
      ijk  = i + j*jj + (kmax-1)*kk;
      if(iindex == 0 && jindex == 0)
        b[ijk] -= c[kmax-1];
      // set dp/dz at top to zero
      else
        b[ijk] += c[kmax-1];
    }

  // call tdma solver
  tdma(a, b, c, p, work2d, work3d);

  return 0;
}
#endif
