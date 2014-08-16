#include "advec_2.h"
#include "grid.h"

__device__ double interp2(double a, double b)
{
  return 0.5*(a + b);
} 

__global__ void advecu_kernel(double * __restrict__ ut, double * __restrict__ u, 
                              double * __restrict__ v, double * __restrict__ w,
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
    ut[ijk] += 
          - (  interp2(u[ijk   ], u[ijk+ii]) * interp2(u[ijk   ], u[ijk+ii])
             - interp2(u[ijk-ii], u[ijk   ]) * interp2(u[ijk-ii], u[ijk   ]) ) * dxi

          - (  interp2(v[ijk-ii+jj], v[ijk+jj]) * interp2(u[ijk   ], u[ijk+jj])
             - interp2(v[ijk-ii   ], v[ijk   ]) * interp2(u[ijk-jj], u[ijk   ]) ) * dyi 

          - (  interp2(w[ijk-ii+kk], w[ijk+kk]) * interp2(u[ijk   ], u[ijk+kk])
             - interp2(w[ijk-ii   ], w[ijk   ]) * interp2(u[ijk-kk], u[ijk   ]) ) * dzi[k];
  }
}

int cadvec_2::advecu_GPU(double * ut, double * u, double * v, double * w, double * dzi)
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxi = 1./grid->dx;
  const double dyi = 1./grid->dy;

  advecu_kernel<<<gridGPU, blockGPU>>>(ut, u, v, w, dzi, dxi, dyi,
                                      grid->icells, grid->ijcells,
                                      grid->istart, grid->jstart, grid->kstart,
                                      grid->iend,   grid->jend, grid->kend);

  cudaError_t error = cudaGetLastError();
  if(error != cudaSuccess)
    printf("CUDA ERROR: %s\n", cudaGetErrorString(error));

  return 0;
}

