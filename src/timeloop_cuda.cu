#include "timeloop.h"
#include "grid.h"
#include "master.h"

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
  const int blocki = 128;
  const int blockj = 2;
  const int gridi = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  rk3_kernel<<<gridGPU, blockGPU>>>(a, at, dt,
                                    substep, grid->icells, grid->ijcells,
                                    grid->istart, grid->jstart, grid->kstart,
                                    grid->iend, grid->jend, grid->kend);

  return 0;
}

int ctimeloop::rk4_GPU(double *a, double *at, double dt)
{
  const int blocki = 128;
  const int blockj = 2;
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
