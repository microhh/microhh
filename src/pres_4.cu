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

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "pres.h"
#include "pres_4.h"
#include "fd.h"
#include "tools.h"

using namespace fd::o4;

namespace Pres_4_g
{
  __global__ void gcwt(double * const __restrict__ wt,
                       const int jj, const int kk,
                       const int istart, const int jstart, const int kstart,
                       const int iend, const int jend, const int kend)
  {
    const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
    const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

    if(i < iend && j < jend)
    {
      int ijk = i + j*jj + kstart*kk;
      wt[ijk-kk] = -wt[ijk+kk];

      ijk = i + j*jj + kend*kk;
      wt[ijk+kk] = -wt[ijk-kk];
    }
  }

  __global__ void presin(double * const __restrict__ p,
                         const double * const __restrict__ u , const double * const __restrict__ v , const double * const __restrict__ w ,
                         const double * const __restrict__ ut, const double * const __restrict__ vt, const double * const __restrict__ wt,
                         const double * const __restrict__ dzi4,
                         const double dxi, const double dyi, const double dti,
                         const int jj, const int kk,
                         const int jjp, const int kkp,
                         const int imax, const int jmax, const int kmax,
                         const int igc, const int jgc, const int kgc)
  {
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    const int j = blockIdx.y*blockDim.y + threadIdx.y;
    const int k = blockIdx.z;

    const int ii1 = 1;
    const int ii2 = 2;
    const int jj1 = 1*jj;
    const int jj2 = 2*jj;
    const int kk1 = 1*kk;
    const int kk2 = 2*kk;

    if(i < imax && j < jmax && k < kmax)
    {
      const int ijkp = i + j*jjp + k*kkp;
      const int ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;

      p[ijkp] = (cg0*(ut[ijk-ii1] + u[ijk-ii1]*dti) + cg1*(ut[ijk] + u[ijk]*dti) + cg2*(ut[ijk+ii1] + u[ijk+ii1]*dti) + cg3*(ut[ijk+ii2] + u[ijk+ii2]*dti)) * cgi*dxi
              + (cg0*(vt[ijk-jj1] + v[ijk-jj1]*dti) + cg1*(vt[ijk] + v[ijk]*dti) + cg2*(vt[ijk+jj1] + v[ijk+jj1]*dti) + cg3*(vt[ijk+jj2] + v[ijk+jj2]*dti)) * cgi*dyi
              + (cg0*(wt[ijk-kk1] + w[ijk-kk1]*dti) + cg1*(wt[ijk] + w[ijk]*dti) + cg2*(wt[ijk+kk1] + w[ijk+kk1]*dti) + cg3*(wt[ijk+kk2] + w[ijk+kk2]*dti)) * dzi4[k+kgc];
    }
  }

  __global__ void solvein(const double * const __restrict__ p,
                          const double * const __restrict__ m1, const double * const __restrict__ m2, const double * const __restrict__ m3, const double * const __restrict__ m4,
                          const double * const __restrict__ m5, const double * const __restrict__ m6, const double * const __restrict__ m7,
                          double * const __restrict__ m1temp, double * const __restrict__ m2temp, double * __restrict__ const m3temp, double * const __restrict__ m4temp,
                          double * const __restrict__ m5temp, double * const __restrict__ m6temp, double * __restrict__ const m7temp, double * const __restrict__ ptemp,
                          const double * const __restrict__ bmati, const double * const __restrict__ bmatj,
                          const int mpicoordx, const int mpicoordy,
                          const int iblock, const int jblock,
                          const int kmax,
                          const int n, const int jslice)
  {
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    const int j = blockIdx.y*blockDim.y + threadIdx.y;

    const int jj = iblock;
    const int kk = iblock*jblock;

    const int kki1 = 1*iblock*jslice;
    const int kki2 = 2*iblock*jslice;
    const int kki3 = 3*iblock*jslice;

    int ik,ijk,iindex,jindex;

    if(i < iblock && j < jslice)
    {
      // Swap the mpicoords, because domain is turned 90 degrees to avoid two mpi transposes.
      iindex = mpicoordy*iblock + i;
      jindex = mpicoordx*jblock + n*jslice + j;

      // Set a zero gradient bc at the bottom.
      ik = i + j*jj;
      m1temp[ik] =  0.;
      m2temp[ik] =  0.;
      m3temp[ik] =  0.;
      m4temp[ik] =  1.;
      m5temp[ik] =  0.;
      m6temp[ik] =  0.;
      m7temp[ik] = -1.;
      ptemp [ik] =  0.;

      m1temp[ik+kki1] =  0.;
      m2temp[ik+kki1] =  0.;
      m3temp[ik+kki1] =  0.;
      m4temp[ik+kki1] =  1.;
      m5temp[ik+kki1] = -1.;
      m6temp[ik+kki1] =  0.;
      m7temp[ik+kki1] =  0.;
      ptemp [ik+kki1] =  0.;

      for(int k=0; k<kmax; ++k)
      {
        // Swap the mpicoords, because domain is turned 90 degrees to avoid two mpi transposes.
        ijk = i + (j + n*jslice)*jj + k*kk;
        ik  = i + j*jj + k*kki1;
        m1temp[ik+kki2] = m1[k];
        m2temp[ik+kki2] = m2[k];
        m3temp[ik+kki2] = m3[k];
        m4temp[ik+kki2] = m4[k] + bmati[iindex] + bmatj[jindex];
        m5temp[ik+kki2] = m5[k];
        m6temp[ik+kki2] = m6[k];
        m7temp[ik+kki2] = m7[k];
        ptemp [ik+kki2] = p[ijk];
      }

      // Set the top boundary.
      ik = i + j*jj + kmax*kki1;
      if(iindex == 0 && jindex == 0)
      {
        m1temp[ik+kki2] =    0.;
        m2temp[ik+kki2] = -1/3.;
        m3temp[ik+kki2] =    2.;
        m4temp[ik+kki2] =    1.;

        m1temp[ik+kki3] =   -2.;
        m2temp[ik+kki3] =    9.;
        m3temp[ik+kki3] =    0.;
        m4temp[ik+kki3] =    1.;
      }

      // Set dp/dz at top to zero.
      else
      {
        m1temp[ik+kki2] =  0.;
        m2temp[ik+kki2] =  0.;
        m3temp[ik+kki2] = -1.;
        m4temp[ik+kki2] =  1.;

        m1temp[ik+kki3] = -1.;
        m2temp[ik+kki3] =  0.;
        m3temp[ik+kki3] =  0.;
        m4temp[ik+kki3] =  1.;
      }

      // Set the top boundary.
      m5temp[ik+kki2] = 0.;
      m6temp[ik+kki2] = 0.;
      m7temp[ik+kki2] = 0.;
      ptemp [ik+kki2] = 0.;

      m5temp[ik+kki3] = 0.;
      m6temp[ik+kki3] = 0.;
      m7temp[ik+kki3] = 0.;
      ptemp [ik+kki3] = 0.;
    }
  }

  __global__ void solveputback(double * const __restrict__ p,
                               const double * const __restrict__ ptemp,
                               const int iblock, const int jblock,
                               const int kmax,
                               const int n, const int jslice)
  {
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    const int j = blockIdx.y*blockDim.y + threadIdx.y;

    const int jj = iblock;
    const int kk = iblock*jblock;

    const int kki1 = 1*iblock*jslice;
    const int kki2 = 2*iblock*jslice;

    if(i < iblock && j < jslice)
    {
      // Put back the solution.
      for(int k=0; k<kmax; ++k)
      {
        const int ik  = i + j*jj + k*kki1;
        const int ijk = i + (j + n*jslice)*jj + k*kk;
        p[ijk] = ptemp[ik+kki2];
      }
    }
  }

  __global__ void hdma(double * const __restrict__ m1, double * const __restrict__ m2, double * const __restrict__ m3, double * const __restrict__ m4,
                       double * const __restrict__ m5, double * const __restrict__ m6, double * const __restrict__ m7, double * const __restrict__ p,
                       const int iblock, const int kmax, const int jslice)
  {
    const int i = blockIdx.x*blockDim.x + threadIdx.x;
    const int j = blockIdx.y*blockDim.y + threadIdx.y;

    const int jj = iblock;

    const int kk1 = 1*iblock*jslice;
    const int kk2 = 2*iblock*jslice;
    const int kk3 = 3*iblock*jslice;

    int k,ik;

    if(i < iblock && j < jslice)
    {
      // Use LU factorization.
      k = 0;
      ik = i + j*jj;
      m1[ik] = 1.;
      m2[ik] = 1.;
      m3[ik] = 1.            / m4[ik];
      m4[ik] = 1.;
      m5[ik] = m5[ik]*m3[ik];
      m6[ik] = m6[ik]*m3[ik];
      m7[ik] = m7[ik]*m3[ik];

      k = 1;
      ik = i + j*jj + k*kk1;
      m1[ik] = 1.;
      m2[ik] = 1.;
      m3[ik] = m3[ik]                     / m4[ik-kk1];
      m4[ik] = m4[ik] - m3[ik]*m5[ik-kk1];
      m5[ik] = m5[ik] - m3[ik]*m6[ik-kk1];
      m6[ik] = m6[ik] - m3[ik]*m7[ik-kk1];

      k = 2;
      ik = i + j*jj + k*kk1;
      m1[ik] = 1.;
      m2[ik] =   m2[ik]                                           / m4[ik-kk2];
      m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] ) / m4[ik-kk1];
      m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2];
      m5[ik] =   m5[ik] - m3[ik]*m6[ik-kk1] - m2[ik]*m7[ik-kk2];
      m6[ik] =   m6[ik] - m3[ik]*m7[ik-kk1];

      for(k=3; k<kmax+2; ++k)
      {
        ik = i + j*jj + k*kk1;
        m1[ik] = ( m1[ik]                                                            ) / m4[ik-kk3];
        m2[ik] = ( m2[ik]                                         - m1[ik]*m5[ik-kk3]) / m4[ik-kk2];
        m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] - m1[ik]*m6[ik-kk3]) / m4[ik-kk1];
        m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2] - m1[ik]*m7[ik-kk3];
        m5[ik] =   m5[ik] - m3[ik]*m6[ik-kk1] - m2[ik]*m7[ik-kk2];
        m6[ik] =   m6[ik] - m3[ik]*m7[ik-kk1];
      }

      k = kmax+1;
      ik = i + j*jj + k*kk1;
      m7[ik] = 1.;

      k = kmax+2;
      ik = i + j*jj + k*kk1;
      m1[ik] = ( m1[ik]                                                            ) / m4[ik-kk3];
      m2[ik] = ( m2[ik]                                         - m1[ik]*m5[ik-kk3]) / m4[ik-kk2];
      m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] - m1[ik]*m6[ik-kk3]) / m4[ik-kk1];
      m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2] - m1[ik]*m7[ik-kk3];
      m5[ik] =   m5[ik] - m3[ik]*m6[ik-kk1] - m2[ik]*m7[ik-kk2];
      m6[ik] = 1.;
      m7[ik] = 1.;

      k = kmax+3;
      ik = i + j*jj + k*kk1;
      m1[ik] = ( m1[ik]                                                            ) / m4[ik-kk3];
      m2[ik] = ( m2[ik]                                         - m1[ik]*m5[ik-kk3]) / m4[ik-kk2];
      m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] - m1[ik]*m6[ik-kk3]) / m4[ik-kk1];
      m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2] - m1[ik]*m7[ik-kk3];
      m5[ik] = 1.;
      m6[ik] = 1.;
      m7[ik] = 1.;

      // Do the backward substitution.
      // First, solve Ly = p, forward.
      ik = i + j*jj;
      p[ik    ] =             p[ik    ]*m3[ik    ];
      p[ik+kk1] = p[ik+kk1] - p[ik    ]*m3[ik+kk1];
      p[ik+kk2] = p[ik+kk2] - p[ik+kk1]*m3[ik+kk2] - p[ik]*m2[ik+kk2];

      for(k=3; k<kmax+4; ++k)
      {
        ik = i + j*jj + k*kk1;
        p[ik] = p[ik] - p[ik-kk1]*m3[ik] - p[ik-kk2]*m2[ik] - p[ik-kk3]*m1[ik];
      }

      // Second, solve Ux=y, backward.
      k = kmax+3;
      ik = i + j*jj + k*kk1;
      p[ik    ] =   p[ik    ]                                             / m4[ik    ];
      p[ik-kk1] = ( p[ik-kk1] - p[ik    ]*m5[ik-kk1] )                    / m4[ik-kk1];
      p[ik-kk2] = ( p[ik-kk2] - p[ik-kk1]*m5[ik-kk2] - p[ik]*m6[ik-kk2] ) / m4[ik-kk2];

      for(k=kmax; k>=0; --k)
      {
        ik = i + j*jj + k*kk1;
        p[ik] = ( p[ik] - p[ik+kk1]*m5[ik] - p[ik+kk2]*m6[ik] - p[ik+kk3]*m7[ik] ) / m4[ik];
      }
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

    const int kkp1 = 1*kkp;
    const int kkp2 = 2*kkp;

    if(i < imax && j < jmax && k < kmax)
    {
      const int ijk  = i + j*jj + k*kk;
      const int ijkp = i+istart + (j+jstart)*jjp + (k+kstart)*kkp;

      p[ijkp] = work3d[ijk];

      // set the BC
      if(k == 0)
      {
        p[ijkp-kkp1] = p[ijkp     ];
        p[ijkp-kkp2] = p[ijkp+kkp1];
      }
      else if(k == kmax-1)
      {
        p[ijkp+kkp1] = p[ijkp     ];
        p[ijkp+kkp2] = p[ijkp-kkp1];
      }
    }
  }

  __global__ void presout(double * const __restrict__ ut, double * const __restrict__ vt, double * const __restrict__ wt,
                          const double * const __restrict__ p,
                          const double * const __restrict__ dzhi4,
                          const double dxi, const double dyi,
                          const int jj, const int kk,
                          const int istart, const int jstart, const int kstart,
                          const int iend, const int jend, const int kend)
  {
    const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
    const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
    const int k = blockIdx.z + kstart;

    const int ii1 = 1;
    const int ii2 = 2;
    const int jj1 = 1*jj;
    const int jj2 = 2*jj;
    const int kk1 = 1*kk;
    const int kk2 = 2*kk;

    if(i < iend && j < jend && k == kstart)
    {
      const int ijk = i + j*jj + k*kk;
      ut[ijk] -= (cg0*p[ijk-ii2] + cg1*p[ijk-ii1] + cg2*p[ijk] + cg3*p[ijk+ii1]) * cgi*dxi;
      vt[ijk] -= (cg0*p[ijk-jj2] + cg1*p[ijk-jj1] + cg2*p[ijk] + cg3*p[ijk+jj1]) * cgi*dyi;
    }
    else if(i < iend && j < jend && k < kend)
    {
      const int ijk = i + j*jj1 + k*kk1;
      ut[ijk] -= (cg0*p[ijk-ii2] + cg1*p[ijk-ii1] + cg2*p[ijk] + cg3*p[ijk+ii1]) * cgi*dxi;
      vt[ijk] -= (cg0*p[ijk-jj2] + cg1*p[ijk-jj1] + cg2*p[ijk] + cg3*p[ijk+jj1]) * cgi*dyi;
      wt[ijk] -= (cg0*p[ijk-kk2] + cg1*p[ijk-kk1] + cg2*p[ijk] + cg3*p[ijk+kk1]) * dzhi4[k];
    }
  }

  __global__ void calcdivergence(double * __restrict__ div,
                                 double * __restrict__ u, double * __restrict__ v, double * __restrict__ w,
                                 double * __restrict__ dzi4,
                                 double dxi, double dyi,
                                 int jj, int kk,
                                 int istart, int jstart, int kstart,
                                 int iend, int jend, int kend)
  {
    const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
    const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
    const int k = blockIdx.z + kstart;

    const int ii1 = 1;
    const int ii2 = 2;
    const int jj1 = 1*jj;
    const int jj2 = 2*jj;
    const int kk1 = 1*kk;
    const int kk2 = 2*kk;

    if(i < iend && j < jend && k < kend)
    {
      const int ijk = i + j*jj + k*kk;
      div[ijk] = (cg0*u[ijk-ii1] + cg1*u[ijk] + cg2*u[ijk+ii1] + cg3*u[ijk+ii2]) * cgi*dxi
               + (cg0*v[ijk-jj1] + cg1*v[ijk] + cg2*v[ijk+jj1] + cg3*v[ijk+jj2]) * cgi*dyi
               + (cg0*w[ijk-kk1] + cg1*w[ijk] + cg2*w[ijk+kk1] + cg3*w[ijk+kk2]) * dzi4[k];
    }
  }
} // End namespace.

#ifdef USECUDA
void Pres_4::exec(double dt)
{
  // 1. Create the input for the pressure solver.
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  dim3 grid2dGPU (gridi, gridj);
  dim3 block2dGPU(blocki, blockj);

  dim3 grid1dGPU (gridi);
  dim3 block1dGPU(blocki);

  const int offs = grid->memoffset;

  // calculate the cyclic BCs first
  grid->boundaryCyclic_g(&fields->ut->data_g[offs]);
  grid->boundaryCyclic_g(&fields->vt->data_g[offs]);
  grid->boundaryCyclic_g(&fields->wt->data_g[offs]);

  Pres_4_g::gcwt<<<grid2dGPU, block2dGPU>>>(&fields->wt->data_g[offs],
                                           grid->icellsp, grid->ijcellsp,
                                           grid->istart, grid->jstart, grid->kstart,
                                           grid->iend, grid->jend, grid->kend);
  cudaCheckError();

  Pres_4_g::presin<<<gridGPU, blockGPU>>>(fields->sd["p"]->data_g,
                                         &fields->u ->data_g[offs], &fields->v ->data_g[offs], &fields->w ->data_g[offs],
                                         &fields->ut->data_g[offs], &fields->vt->data_g[offs], &fields->wt->data_g[offs],
                                         grid->dzi4_g,
                                         1./grid->dx, 1./grid->dy, 1./dt,
                                         grid->icellsp, grid->ijcellsp,
                                         grid->imax, grid->imax*grid->jmax,
                                         grid->imax, grid->jmax, grid->kmax,
                                         grid->igc, grid->jgc, grid->kgc);
  cudaCheckError();

  fftForward(fields->sd["p"]->data_g, fields->atmp["tmp1"]->data_g, fields->atmp["tmp2"]->data_g);

  double *tmp1_g = fields->atmp["tmp1"]->data_g;
  double *tmp2_g = fields->atmp["tmp2"]->data_g;

  // Set jslice to a higher value
  const int jslice = std::max(grid->jblock/4, 1);

  const int blockis = 128;
  const int blockjs = 1;
  const int gridis  = grid->iblock/blockis + (grid->iblock%blockis > 0);
  const int gridjs  =       jslice/blockjs + (      jslice%blockjs > 0);

  dim3 grid2dsGPU (gridis , gridjs );
  dim3 block2dsGPU(blockis, blockjs);

  const int ns = grid->iblock*jslice*(grid->kmax+4);
  const int nj = grid->jblock/jslice;

  for(int n=0; n<nj; ++n)
  {
    // Prepare the fields that go into the matrix solver
    Pres_4_g::solvein<<<grid2dsGPU,block2dsGPU>>>(fields->sd["p"]->data_g,
                                                 m1_g, m2_g, m3_g, m4_g,
                                                 m5_g, m6_g, m7_g,
                                                 &tmp1_g[0*ns], &tmp1_g[1*ns], &tmp1_g[2*ns], &tmp1_g[3*ns],
                                                 &tmp2_g[0*ns], &tmp2_g[1*ns], &tmp2_g[2*ns], &tmp2_g[3*ns],
                                                 bmati_g, bmatj_g,
                                                 master->mpicoordx, master->mpicoordy,
                                                 grid->iblock, grid->jblock,
                                                 grid->kmax,
                                                 n, jslice);
    cudaCheckError();

    // Solve the sevenbanded matrix
    Pres_4_g::hdma<<<grid2dsGPU,block2dsGPU>>>(&tmp1_g[0*ns], &tmp1_g[1*ns], &tmp1_g[2*ns], &tmp1_g[3*ns],
                                              &tmp2_g[0*ns], &tmp2_g[1*ns], &tmp2_g[2*ns], &tmp2_g[3*ns],
                                              grid->iblock, grid->kmax, jslice);
    cudaCheckError();

    // Put the solution back into the pressure field
    Pres_4_g::solveputback<<<grid2dsGPU,block2dsGPU>>>(fields->sd["p"]->data_g,
                                                      &tmp2_g[3*ns],
                                                      grid->iblock, grid->jblock,
                                                      grid->kmax,
                                                      n, jslice);
    cudaCheckError();
  }

  fftBackward(fields->sd["p"]->data_g, fields->atmp["tmp1"]->data_g, fields->atmp["tmp2"]->data_g);

  cudaSafeCall(cudaMemcpy(fields->atmp["tmp1"]->data_g, fields->sd["p"]->data_g, grid->ncellsp*sizeof(double), cudaMemcpyDeviceToDevice));

  Pres_4_g::solveout<<<gridGPU, blockGPU>>>(&fields->sd["p"]->data_g[offs], fields->atmp["tmp1"]->data_g,
                                           grid->imax, grid->imax*grid->jmax,
                                           grid->icellsp, grid->ijcellsp,
                                           grid->istart, grid->jstart, grid->kstart,
                                           grid->imax, grid->jmax, grid->kmax);
  cudaCheckError();

  grid->boundaryCyclic_g(&fields->sd["p"]->data_g[offs]);

  // 3. Get the pressure tendencies from the pressure field.
  Pres_4_g::presout<<<gridGPU, blockGPU>>>(&fields->ut->data_g[offs], &fields->vt->data_g[offs], &fields->wt->data_g[offs],
                                          &fields->sd["p"]->data_g[offs],
                                          grid->dzhi4_g,
                                          1./grid->dx, 1./grid->dy,
                                          grid->icellsp, grid->ijcellsp,
                                          grid->istart, grid->jstart, grid->kstart,
                                          grid->iend, grid->jend, grid->kend);
  cudaCheckError();
}

double Pres_4::checkDivergence()
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const int offs = grid->memoffset;

  Pres_4_g::calcdivergence<<<gridGPU, blockGPU>>>(&fields->atmp["tmp1"]->data_g[offs],
                                                 &fields->u->data_g[offs], &fields->v->data_g[offs], &fields->w->data_g[offs],
                                                 grid->dzi4_g,
                                                 grid->dxi, grid->dyi,
                                                 grid->icellsp, grid->ijcellsp,
                                                 grid->istart, grid->jstart, grid->kstart,
                                                 grid->iend, grid->jend, grid->kend);
  cudaCheckError();

  double divmax = grid->getMax_g(&fields->atmp["tmp1"]->data_g[offs], fields->atmp["tmp2"]->data_g);
  grid->getMax(&divmax);

  return divmax;
}

void Pres_4::prepare_device()
{
  const int kmemsize = grid->kmax*sizeof(double);
  const int imemsize = grid->itot*sizeof(double);
  const int jmemsize = grid->jtot*sizeof(double);

  cudaSafeCall(cudaMalloc((void**)&bmati_g, imemsize));
  cudaSafeCall(cudaMalloc((void**)&bmatj_g, jmemsize));

  cudaSafeCall(cudaMalloc((void**)&m1_g, kmemsize));
  cudaSafeCall(cudaMalloc((void**)&m2_g, kmemsize));
  cudaSafeCall(cudaMalloc((void**)&m3_g, kmemsize));
  cudaSafeCall(cudaMalloc((void**)&m4_g, kmemsize));
  cudaSafeCall(cudaMalloc((void**)&m5_g, kmemsize));
  cudaSafeCall(cudaMalloc((void**)&m6_g, kmemsize));
  cudaSafeCall(cudaMalloc((void**)&m7_g, kmemsize));

  cudaSafeCall(cudaMemcpy(bmati_g, bmati, imemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(bmatj_g, bmatj, jmemsize, cudaMemcpyHostToDevice));

  cudaSafeCall(cudaMemcpy(m1_g, m1, kmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(m2_g, m2, kmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(m3_g, m3, kmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(m4_g, m4, kmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(m5_g, m5, kmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(m6_g, m6, kmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(m7_g, m7, kmemsize, cudaMemcpyHostToDevice));

  makeCufftPlan();
}

void Pres_4::clear_device()
{
  cudaSafeCall(cudaFree(bmati_g));
  cudaSafeCall(cudaFree(bmatj_g));

  cudaSafeCall(cudaFree(m1_g));
  cudaSafeCall(cudaFree(m2_g));
  cudaSafeCall(cudaFree(m3_g));
  cudaSafeCall(cudaFree(m4_g));
  cudaSafeCall(cudaFree(m5_g));
  cudaSafeCall(cudaFree(m6_g));
  cudaSafeCall(cudaFree(m7_g));
}
#endif
