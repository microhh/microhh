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
#include "pres_4.h"
#include "fd.h"

using namespace fd::o4;

__global__ void pres_4_gcwt(double * const __restrict__ wt,
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

__global__ void pres_4_presin(double * const __restrict__ p,
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

__global__ void pres_4_presout(double * const __restrict__ ut, double * const __restrict__ vt, double * const __restrict__ wt,
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

__global__ void pres_4_calcdivergence(double * __restrict__ div,
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

#ifdef USECUDA
void cpres_4::exec(double dt)
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

  const int offs = grid->memoffset;

  // calculate the cyclic BCs first
  grid->boundary_cyclic_g(&fields->ut->data_g[offs]);
  grid->boundary_cyclic_g(&fields->vt->data_g[offs]);
  grid->boundary_cyclic_g(&fields->wt->data_g[offs]);

  pres_4_gcwt<<<grid2dGPU, block2dGPU>>>(&fields->wt->data_g[offs],
                                         grid->icellsp, grid->ijcellsp,
                                         grid->istart, grid->jstart, grid->kstart,
                                         grid->iend, grid->jend, grid->kend);

  pres_4_presin<<<gridGPU, blockGPU>>>(fields->sd["p"]->data_g,
                                       &fields->u ->data_g[offs], &fields->v ->data_g[offs], &fields->w ->data_g[offs],
                                       &fields->ut->data_g[offs], &fields->vt->data_g[offs], &fields->wt->data_g[offs],
                                       grid->dzi4_g,
                                       1./grid->dx, 1./grid->dy, 1./dt,
                                       grid->icellsp, grid->ijcellsp,
                                       grid->imax, grid->imax*grid->jmax,
                                       grid->imax, grid->jmax, grid->kmax,
                                       grid->igc, grid->jgc, grid->kgc);

  // Forward FFT -> how to get rid of the loop at the host side....
  int nsize = sizeof(double)*grid->ncells;
  cudaMemcpy(fields->a["p"]->data, fields->a["p"]->data_g, nsize, cudaMemcpyDeviceToHost);

  // 2. Solve the Poisson equation using FFTs and a heptadiagonal solver
  // Take slices out of a temporary field to save memory. The temp arrays
  // are always big enough, this cannot fail.
  double *tmp2 = fields->sd["tmp2"]->data;
  double *tmp3 = fields->sd["tmp3"]->data;
  const int ns = grid->iblock*(grid->kmax+4);
  pres_solve(fields->sd["p"]->data, fields->sd["tmp1"]->data, grid->dz,
             m1, m2, m3, m4,
             m5, m6, m7,
             &tmp2[0*ns], &tmp2[1*ns], &tmp2[2*ns], &tmp2[3*ns], 
             &tmp3[0*ns], &tmp3[1*ns], &tmp3[2*ns], &tmp3[3*ns], 
             bmati, bmatj);

  // 3. Get the pressure tendencies from the pressure field.
  cudaMemcpy(fields->a["p"]->data_g, fields->a["p"]->data, nsize, cudaMemcpyHostToDevice);
  pres_4_presout<<<gridGPU, blockGPU>>>(&fields->ut->data_g[offs], &fields->vt->data_g[offs], &fields->wt->data_g[offs],
                                        &fields->sd["p"]->data_g[offs],
                                        grid->dzhi4_g,
                                        1./grid->dx, 1./grid->dy,
                                        grid->icellsp, grid->ijcellsp,
                                        grid->istart, grid->jstart, grid->kstart,
                                        grid->iend, grid->jend, grid->kend);
}

double cpres_4::check()
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const int offs = grid->memoffset;

  pres_4_calcdivergence<<<gridGPU, blockGPU>>>(&fields->a["tmp1"]->data_g[offs],
                                               &fields->u->data_g[offs], &fields->v->data_g[offs], &fields->w->data_g[offs],
                                               grid->dzi4_g,
                                               grid->dxi, grid->dyi,
                                               grid->icellsp, grid->ijcellsp,
                                               grid->istart,  grid->jstart, grid->kstart,
                                               grid->iend,    grid->jend,   grid->kend);

  double divmax = grid->getmax_g(&fields->a["tmp1"]->data_g[offs], fields->a["tmp2"]->data_g);
  grid->getmax(&divmax);

  return divmax;
}

int cpres_4::prepareGPU()
{
  const int kmemsize = grid->kmax*sizeof(double);
  const int imemsize = grid->itot*sizeof(double);
  const int jmemsize = grid->jtot*sizeof(double);

  cudaMalloc((void**)&bmati_g, imemsize);
  cudaMalloc((void**)&bmatj_g, jmemsize);
  cudaMalloc((void**)&m1_g, kmemsize);
  cudaMalloc((void**)&m2_g, kmemsize);
  cudaMalloc((void**)&m3_g, kmemsize);
  cudaMalloc((void**)&m4_g, kmemsize);
  cudaMalloc((void**)&m5_g, kmemsize);
  cudaMalloc((void**)&m6_g, kmemsize);
  cudaMalloc((void**)&m7_g, kmemsize);

  cudaMemcpy(bmati_g, bmati, imemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(bmatj_g, bmatj, jmemsize, cudaMemcpyHostToDevice);

  cudaMemcpy(m1_g, m1, kmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(m2_g, m2, kmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(m3_g, m3, kmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(m4_g, m4, kmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(m5_g, m5, kmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(m6_g, m6, kmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(m7_g, m7, kmemsize, cudaMemcpyHostToDevice);

  // cuFFT
  cudaMalloc((void **)&ffti_complex_g, sizeof(cufftDoubleComplex)*(grid->jtot * (grid->itot/2+1))); // sizeof(complex) = 16
  cudaMalloc((void **)&fftj_complex_g, sizeof(cufftDoubleComplex)*(grid->itot * (grid->jtot/2+1)));

  // Make cuFFT plan
  int rank      = 1;

  // Double input
  int i_ni[]    = {grid->itot}; 
  int i_nj[]    = {grid->jtot};  
  int i_istride = 1;
  int i_jstride = grid->itot;
  int i_idist   = grid->itot;
  int i_jdist   = 1;

  // Double-complex output
  int o_ni[]    = {grid->itot/2+1};
  int o_nj[]    = {grid->jtot/2+1};
  int o_istride = 1;
  int o_jstride = grid->itot;
  int o_idist   = grid->itot/2+1;
  int o_jdist   = 1;

  // Forward FFTs
  cufftPlanMany(&iplanf, rank, i_ni, i_ni, i_istride, i_idist, o_ni, o_istride, o_idist, CUFFT_D2Z, grid->jtot);
  cufftPlanMany(&jplanf, rank, i_nj, i_nj, i_jstride, i_jdist, o_nj, o_jstride, o_jdist, CUFFT_D2Z, grid->itot);

  // Backward FFTs
  // NOTE: input size is always the 'logical' size of the FFT, so itot or jtot, not itot/2+1 or jtot/2+1 
  cufftPlanMany(&iplanb, rank, i_ni, o_ni, o_istride, o_idist, i_ni, i_istride, i_idist, CUFFT_Z2D, grid->jtot);
  cufftPlanMany(&jplanb, rank, i_nj, o_nj, o_jstride, o_jdist, i_nj, i_jstride, i_jdist, CUFFT_Z2D, grid->itot);

  return 0;
}
#endif
