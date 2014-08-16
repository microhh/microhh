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

__global__ void advecv_kernel(double * __restrict__ vt, double * __restrict__ u, 
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
    vt[ijk] += 
          - (  interp2(u[ijk+ii-jj], u[ijk+ii]) * interp2(v[ijk   ], v[ijk+ii])
             - interp2(u[ijk   -jj], u[ijk   ]) * interp2(v[ijk-ii], v[ijk   ]) ) * dxi

          - (  interp2(v[ijk   ], v[ijk+jj]) * interp2(v[ijk   ], v[ijk+jj])
             - interp2(v[ijk-jj], v[ijk   ]) * interp2(v[ijk-jj], v[ijk   ]) ) * dyi

          - (  interp2(w[ijk-jj+kk], w[ijk+kk]) * interp2(v[ijk   ], v[ijk+kk])
             - interp2(w[ijk-jj   ], w[ijk   ]) * interp2(v[ijk-kk], v[ijk   ]) ) * dzi[k];
  }
}

__global__ void advecw_kernel(double * __restrict__ wt, double * __restrict__ u, 
                              double * __restrict__ v, double * __restrict__ w,
                              double * __restrict__ dzhi, double dxi, double dyi, 
                              int jj, int kk, int istart, int jstart, int kstart,
                              int iend,   int jend,   int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart + 1;
  int ii = 1;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    wt[ijk] += 
          - (  interp2(u[ijk+ii-kk], u[ijk+ii]) * interp2(w[ijk   ], w[ijk+ii])
             - interp2(u[ijk   -kk], u[ijk   ]) * interp2(w[ijk-ii], w[ijk   ]) ) * dxi

          - (  interp2(v[ijk+jj-kk], v[ijk+jj]) * interp2(w[ijk   ], w[ijk+jj])
             - interp2(v[ijk   -kk], v[ijk   ]) * interp2(w[ijk-jj], w[ijk   ]) ) * dyi

          - (  interp2(w[ijk   ], w[ijk+kk]) * interp2(w[ijk   ], w[ijk+kk])
             - interp2(w[ijk-kk], w[ijk   ]) * interp2(w[ijk-kk], w[ijk   ]) ) * dzhi[k];
  }
}

__global__ void advecs_kernel(double * __restrict__ st, double * __restrict__ s, 
                              double * __restrict__ u, double * __restrict__ v, double * __restrict__ w,
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

          - (  w[ijk+kk] * interp2(s[ijk   ], s[ijk+kk])
             - w[ijk   ] * interp2(s[ijk-kk], s[ijk   ]) ) * dzi[k];
  }
}

int cadvec_2::advecu_g(double * ut, double * u, double * v, double * w, double * dzi)
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

int cadvec_2::advecv_g(double * vt, double * u, double * v, double * w, double * dzi)
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxi = 1./grid->dx;
  const double dyi = 1./grid->dy;

  advecv_kernel<<<gridGPU, blockGPU>>>(vt, u, v, w, dzi, dxi, dyi,
                                      grid->icells, grid->ijcells,
                                      grid->istart, grid->jstart, grid->kstart,
                                      grid->iend,   grid->jend, grid->kend);

  cudaError_t error = cudaGetLastError();
  if(error != cudaSuccess)
    printf("CUDA ERROR: %s\n", cudaGetErrorString(error));

  return 0;
}

int cadvec_2::advecw_g(double * wt, double * u, double * v, double * w, double * dzhi)
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxi = 1./grid->dx;
  const double dyi = 1./grid->dy;

  advecw_kernel<<<gridGPU, blockGPU>>>(wt, u, v, w, dzhi, dxi, dyi,
                                      grid->icells, grid->ijcells,
                                      grid->istart, grid->jstart, grid->kstart,
                                      grid->iend,   grid->jend, grid->kend);

  cudaError_t error = cudaGetLastError();
  if(error != cudaSuccess)
    printf("CUDA ERROR: %s\n", cudaGetErrorString(error));

  return 0;
}

int cadvec_2::advecs_g(double * st, double * s, double * u, double * v, double * w, double * dzi)
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxi = 1./grid->dx;
  const double dyi = 1./grid->dy;

  advecs_kernel<<<gridGPU, blockGPU>>>(st, s, u, v, w, dzi, dxi, dyi,
                                      grid->icells, grid->ijcells,
                                      grid->istart, grid->jstart, grid->kstart,
                                      grid->iend,   grid->jend, grid->kend);

  cudaError_t error = cudaGetLastError();
  if(error != cudaSuccess)
    printf("CUDA ERROR: %s\n", cudaGetErrorString(error));

  return 0;
}
