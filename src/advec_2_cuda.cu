#include "advec_2.h"
#include "grid.h"
#include "fields.h"

#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <thrust/device_ptr.h>

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

__global__ void calccfl_kernel(double * __restrict__ u, double * __restrict__ v, double * __restrict__ w, 
                               double * __restrict__ tmp1, double * __restrict__ dzi, double dxi, double dyi, 
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
    tmp1[ijk] = std::abs(interp2(u[ijk], u[ijk+ii]))*dxi + std::abs(interp2(v[ijk], v[ijk+jj]))*dyi + std::abs(interp2(w[ijk], w[ijk+kk]))*dzi[k];

    //if(i==5 && j==5 && k==5)
    //  printf("cfl gpu=%f\n",tmp1[ijk]);

  }
}

#ifdef USECUDA
int cadvec_2::exec()
{
  fields->forwardGPU();

  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxi = 1./grid->dx;
  const double dyi = 1./grid->dy;

  advecu_kernel<<<gridGPU, blockGPU>>>(fields->ut->data_g, fields->u->data_g, fields->v->data_g, 
                                       fields->w->data_g, grid->dzi_g, dxi, dyi,
                                       grid->icells, grid->ijcells,
                                       grid->istart, grid->jstart, grid->kstart,
                                       grid->iend,   grid->jend, grid->kend);

  advecv_kernel<<<gridGPU, blockGPU>>>(fields->vt->data_g, fields->u->data_g, fields->v->data_g, 
                                       fields->w->data_g, grid->dzi_g, dxi, dyi,
                                       grid->icells, grid->ijcells,
                                       grid->istart, grid->jstart, grid->kstart,
                                       grid->iend,   grid->jend, grid->kend);

  advecw_kernel<<<gridGPU, blockGPU>>>(fields->wt->data_g, fields->u->data_g, fields->v->data_g, 
                                       fields->w->data_g, grid->dzhi_g, dxi, dyi,
                                       grid->icells, grid->ijcells,
                                       grid->istart, grid->jstart, grid->kstart,
                                       grid->iend,   grid->jend, grid->kend);

  for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    advecs_kernel<<<gridGPU, blockGPU>>>((*it->second).data_g, (*fields->s[it->first]).data_g, 
                                         fields->u->data_g, fields->v->data_g, fields->w->data_g, 
                                         grid->dzi_g, dxi, dyi,
                                         grid->icells, grid->ijcells,
                                         grid->istart, grid->jstart, grid->kstart,
                                         grid->iend,   grid->jend, grid->kend);

  cudaError_t error = cudaGetLastError();
  if(error != cudaSuccess)
    printf("CUDA ERROR: %s\n", cudaGetErrorString(error));

  fields->backwardGPU();
  return 0;
}
#endif

double cadvec_2::calccfl2(double * u, double * v, double * w, double * dzi, double dt)
{
  fields->forwardGPU();

  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);
  double cfl = 0;

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxi = 1./grid->dx;
  const double dyi = 1./grid->dy;

  calccfl_kernel<<<gridGPU, blockGPU>>>(fields->u->data_g, fields->v->data_g, fields->w->data_g, 
                                        (*fields->a["tmp1"]).data_g, grid->dzi_g, dxi, dyi,
                                        grid->icells, grid->ijcells,
                                        grid->istart, grid->jstart, grid->kstart,
                                        grid->iend,   grid->jend, grid->kend);

  //thrust::device_ptr<double> cfl_g = thrust::device_pointer_cast((*fields->a["tmp1"]).data_g);
  //cfl = thrust::reduce(cfl_g, cfl_g + grid->ncells, -1, thrust::maximum<double>()); 

  fields->backwardGPU();

  return cfl;
}


