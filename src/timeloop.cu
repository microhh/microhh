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

#include "timeloop.h"
#include "grid.h"
#include "master.h"
#include "constants.h"

#define cA0 0.
#define cA1 -5./9.
#define cA2 -153./128.
#define cB0 1./3.
#define cB1 15./16.
#define cB2 8./15.

__global__ void rk3_kernel(double * __restrict__ a, double * __restrict__ at, double dt,
                           const int substep, const int jj, const int kk,
                           const int istart, const int jstart, const int kstart,
                           const int iend, const int jend, const int kend)
{
  const double cA[] = {0., -5./9., -153./128.};
  const double cB[] = {1./3., 15./16., 8./15.};

  const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  const int k = blockIdx.z + kstart;

  if(i < iend && j < jend && k < kend)
  {
    const int ijk = i + j*jj + k*kk;
    a[ijk] = a[ijk] + cB[substep]*dt*at[ijk];

    const int substepn = (substep+1) % 3;
    // substep 0 resets the tendencies, because cA[0] == 0
    at[ijk] = cA[substepn]*at[ijk];
  }
}

template<int substep>
__global__ void rk3_kernel2(double * __restrict__ a, double * __restrict__ at, double dt,
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

    switch(substep)
    {
      case 0:
        a[ijk]  = a[ijk] + cB0*dt*at[ijk];
        at[ijk] = cA1*at[ijk];
        break;
      case 1:
        a[ijk]  = a[ijk] + cB1*dt*at[ijk];
        at[ijk] = cA2*at[ijk];
        break;
      case 2:
        a[ijk]  = a[ijk] + cB2*dt*at[ijk];
        at[ijk] = 0.; 
        break;
    }
  }
}

__global__ void rk4_kernel(double * __restrict__ a, double * __restrict__ at, double dt,
                           const int substep, const int jj, const int kk,
                           const int istart, const int jstart, const int kstart,
                           const int iend, const int jend, const int kend)
{
  const double cA [] = {
      0.,
    - 567301805773./1357537059087.,
    -2404267990393./2016746695238.,
    -3550918686646./2091501179385.,
    -1275806237668./ 842570457699.};

  const double cB [] = {
    1432997174477./ 9575080441755.,
    5161836677717./13612068292357.,
    1720146321549./ 2090206949498.,
    3134564353537./ 4481467310338.,
    2277821191437./14882151754819.};
  
  const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  const int k = blockIdx.z + kstart;

  if(i < iend && j < jend && k < kend)
  {
    const int ijk = i + j*jj + k*kk;
    a[ijk] = a[ijk] + cB[substep]*dt*at[ijk];

    const int substepn = (substep+1) % 5;
    // substep 0 resets the tendencies, because cA[0] == 0
    at[ijk] = cA[substepn]*at[ijk];
  }
}

int ctimeloop::rk3_GPU(double *a, double *at, double dt)
{
  const int blocki = cuda::blockSizeI;
  const int blockj = cuda::blockSizeJ;
  const int gridi = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const int offs = grid->memoffset;

  if(substep==0) {
    rk3_kernel2<0><<<gridGPU, blockGPU>>>(&a[offs], &at[offs], dt,
                                      grid->icellsp, grid->ijcellsp,
                                      grid->istart,  grid->jstart, grid->kstart,
                                      grid->iend,    grid->jend,   grid->kend); }
  else if(substep==1) {
    rk3_kernel2<1><<<gridGPU, blockGPU>>>(&a[offs], &at[offs], dt,
                                      grid->icellsp, grid->ijcellsp,
                                      grid->istart,  grid->jstart, grid->kstart,
                                      grid->iend,    grid->jend,   grid->kend); }
  else if(substep==2) {
    rk3_kernel2<2><<<gridGPU, blockGPU>>>(&a[offs], &at[offs], dt,
                                      grid->icellsp, grid->ijcellsp,
                                      grid->istart,  grid->jstart, grid->kstart,
                                      grid->iend,    grid->jend,   grid->kend); }

  //rk3_kernel<<<gridGPU, blockGPU>>>(a, at, dt,
  //                                  substep, grid->icells, grid->ijcells,
  //                                  grid->istart, grid->jstart, grid->kstart,
  //                                  grid->iend, grid->jend, grid->kend);

  return 0;
}

int ctimeloop::rk4_GPU(double *a, double *at, double dt)
{
  const int blocki = cuda::blockSizeI;
  const int blockj = cuda::blockSizeJ;
  const int gridi = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  rk4_kernel<<<gridGPU, blockGPU>>>(a, at, dt,
                                    substep, grid->icells, grid->ijcells,
                                    grid->istart, grid->jstart, grid->kstart,
                                    grid->iend, grid->jend, grid->kend);

  return 0;
}
