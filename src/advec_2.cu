#include "advec_2.h"
#include "grid.h"
#include "fields.h"

template <unsigned int blockSize, bool nispow2>
__device__ void reduceBlocks(const double *data, double *idata, unsigned int n)
{
  extern __shared__ double sdata[];
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*(blockSize*2) + threadIdx.x;
  unsigned int gridSize = blockSize*2*gridDim.x;

  double max = -1e12;
  while (i < n)
  {
    max = fmax(max,data[i]);
    if (nispow2 || i + blockSize < n)
      max = fmax(max,data[i+blockSize]);
    i += gridSize;
  }
  sdata[tid] = max;
  __syncthreads();

  if (blockSize >= 512) { if (tid < 256) { sdata[tid] = fmax(sdata[tid],sdata[tid + 256]); } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { sdata[tid] = fmax(sdata[tid],sdata[tid + 128]); } __syncthreads(); }
  if (blockSize >= 128) { if (tid <  64) { sdata[tid] = fmax(sdata[tid],sdata[tid +  64]); } __syncthreads(); }

  if (tid < 32)
  {
    volatile double *smem = sdata;
    if (blockSize >=  64) { smem[tid] = fmax(smem[tid],smem[tid + 32]); }
    if (blockSize >=  32) { smem[tid] = fmax(smem[tid],smem[tid + 16]); }
    if (blockSize >=  16) { smem[tid] = fmax(smem[tid],smem[tid +  8]); }
    if (blockSize >=   8) { smem[tid] = fmax(smem[tid],smem[tid +  4]); }
    if (blockSize >=   4) { smem[tid] = fmax(smem[tid],smem[tid +  2]); }
    if (blockSize >=   2) { smem[tid] = fmax(smem[tid],smem[tid +  1]); }
  }

  if (tid == 0)
    idata[blockIdx.x] = sdata[0];
}

template <unsigned int blockSize, bool nispow2>
__global__ void reduceMultiPass(const double *data, double *idata, unsigned int n)
{
  reduceBlocks<blockSize, nispow2>(data, idata, n);
}

bool ispow2(unsigned int x)
{
  return ((x&(x-1))==0);
}

void max_multipass(int size, int threads, int blocks, double *data, double *idata)
{
  dim3 dimBlock(threads, 1, 1);
  dim3 dimGrid(blocks, 1, 1);
  unsigned int nsmem = max(64,threads);
  int smemSize = nsmem * sizeof(double);

  if (ispow2(size))
  {
    switch (threads)
    {
      case 512:
        reduceMultiPass<512, true><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case 256:
        reduceMultiPass<256, true><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case 128:
        reduceMultiPass<128, true><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case 64:
        reduceMultiPass< 64, true><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case 32:
        reduceMultiPass< 32, true><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case 16:
        reduceMultiPass< 16, true><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case  8:
        reduceMultiPass<  8, true><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case  4:
        reduceMultiPass<  4, true><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case  2:
        reduceMultiPass<  2, true><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case  1:
        reduceMultiPass<  1, true><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
    }
  }
  else
  {
    switch (threads)
    {
      case 512:
        reduceMultiPass<512, false><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case 256:
        reduceMultiPass<256, false><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case 128:
        reduceMultiPass<128, false><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case 64:
        reduceMultiPass< 64, false><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case 32:
        reduceMultiPass< 32, false><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case 16:
        reduceMultiPass< 16, false><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case  8:
        reduceMultiPass<  8, false><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case  4:
        reduceMultiPass<  4, false><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case  2:
        reduceMultiPass<  2, false><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
      case  1:
        reduceMultiPass<  1, false><<< dimGrid, dimBlock, smemSize >>>(data, idata, size); break;
    }
  }
}

__device__ double interp2(double a, double b)
{
  return 0.5*(a + b);
}

__global__ void advec_2_advecu(double * __restrict__ ut, double * __restrict__ u, 
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

__global__ void advec_2_advecv(double * __restrict__ vt, double * __restrict__ u, 
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

__global__ void advec_2_advecw(double * __restrict__ wt, double * __restrict__ u, 
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

__global__ void advec_2_advecs(double * __restrict__ st, double * __restrict__ s, 
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

__global__ void advec_2_calccfl(double * __restrict__ u, double * __restrict__ v, double * __restrict__ w, 
                                double * __restrict__ tmp1, double * __restrict__ dzi, double dxi, double dyi, 
                                int jj, int kk, int istart, int jstart, int kstart,
                                int iend, int jend, int kend,
                                int icells, int jcells, int kcells)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x; 
  int j = blockIdx.y*blockDim.y + threadIdx.y; 
  int k = blockIdx.z; 
  int ii = 1;
  int ijk = i + j*jj + k*kk;

  if(i >= istart && i < iend && j >= jstart && j < jend && k >= kstart && k < kend)
    tmp1[ijk] = std::abs(interp2(u[ijk], u[ijk+ii]))*dxi + 
                std::abs(interp2(v[ijk], v[ijk+jj]))*dyi + 
                std::abs(interp2(w[ijk], w[ijk+kk]))*dzi[k];
  else if(i < icells && j < jcells && k < kcells) 
    tmp1[ijk] = 0.;
}

#ifdef USECUDA
int cadvec_2::exec()
{
  //fields->forwardGPU();

  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxi = 1./grid->dx;
  const double dyi = 1./grid->dy;

  // Cuda timer:
  //cudaEvent_t start, stop;
  //cudaEventCreate(&start);
  //cudaEventCreate(&stop);
  //cudaEventRecord(start);

  advec_2_advecu<<<gridGPU, blockGPU>>>(fields->ut->data_g, fields->u->data_g, fields->v->data_g, 
                                        fields->w->data_g, grid->dzi_g, dxi, dyi,
                                        grid->icells, grid->ijcells,
                                        grid->istart, grid->jstart, grid->kstart,
                                        grid->iend,   grid->jend, grid->kend);

  advec_2_advecv<<<gridGPU, blockGPU>>>(fields->vt->data_g, fields->u->data_g, fields->v->data_g, 
                                        fields->w->data_g, grid->dzi_g, dxi, dyi,
                                        grid->icells, grid->ijcells,
                                        grid->istart, grid->jstart, grid->kstart,
                                        grid->iend,   grid->jend, grid->kend);

  advec_2_advecw<<<gridGPU, blockGPU>>>(fields->wt->data_g, fields->u->data_g, fields->v->data_g, 
                                        fields->w->data_g, grid->dzhi_g, dxi, dyi,
                                        grid->icells, grid->ijcells,
                                        grid->istart, grid->jstart, grid->kstart,
                                        grid->iend,   grid->jend, grid->kend);

  for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    advec_2_advecs<<<gridGPU, blockGPU>>>((*it->second).data_g, (*fields->s[it->first]).data_g, 
                                          fields->u->data_g, fields->v->data_g, fields->w->data_g, 
                                          grid->dzi_g, dxi, dyi,
                                          grid->icells, grid->ijcells,
                                          grid->istart, grid->jstart, grid->kstart,
                                          grid->iend,   grid->jend, grid->kend);

  
  //cudaEventRecord(stop);  
  //cudaEventSynchronize(stop);
  //float elapsed;
  //cudaEventElapsedTime(&elapsed, start, stop);
  //printf("advec_2 = %f ms\n",elapsed);

  cudaError_t error = cudaGetLastError();
  if(error != cudaSuccess)
    printf("CUDA ERROR: %s\n", cudaGetErrorString(error));

  //fields->backwardGPU();
  return 0;
}
#endif

#ifdef USECUDA
double cadvec_2::calccfl(double * u, double * v, double * w, double * dzi, double dt)
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->icells/blocki + (grid->icells%blocki > 0);
  const int gridj  = grid->jcells/blockj + (grid->jcells%blockj > 0);
  double cfl = 0;

  dim3 gridGPU (gridi, gridj, grid->kcells);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxi = 1./grid->dx;
  const double dyi = 1./grid->dy;

  //fields->forwardGPU();

  advec_2_calccfl<<<gridGPU, blockGPU>>>(fields->u->data_g, fields->v->data_g, fields->w->data_g, 
                                         (*fields->a["tmp1"]).data_g, grid->dzi_g, dxi, dyi,
                                         grid->icells, grid->ijcells,
                                         grid->istart, grid->jstart, grid->kstart,
                                         grid->iend,   grid->jend, grid->kend,
                                         grid->icells, grid->jcells, grid->kcells);

  // Reduce cfl in tmp field
  unsigned int blocksize = 256; // # threads / block
  unsigned int nblocks = (int)ceil((double)grid->ncells / (double)blocksize); // # of blocks

  double *idata; // Array holding max values from first reduction step
  double *dmax; // Array holding max value from second reduction step

  cudaMalloc((void**)&idata, (nblocks)*sizeof(double));
  cudaMalloc((void**)&dmax, sizeof(double));

  // Do the reduction in two steps
  max_multipass(grid->ncells, blocksize, nblocks, fields->a["tmp1"]->data_g,  idata);
  max_multipass(nblocks,      blocksize, 1,       idata,                      dmax);

  // Get max value from GPU
  cudaMemcpy(&cfl, dmax, sizeof(double), cudaMemcpyDeviceToHost);

  // For now allocate and free per call, should probably be done once in init function
  cudaFree(idata);
  cudaFree(dmax);

  grid->getmax(&cfl);

  cfl = cfl*dt;

  return cfl;
}
#endif


