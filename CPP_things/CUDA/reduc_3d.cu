// nvcc reduc_3d.cu -o reduc_3d -arch=sm_21 -I/opt/cuda/include/
// /usr/local/cuda/bin/nvcc reduc_3d.cu -o reduc_3d -arch=sm_21 

#include <stdio.h>
#include <stdlib.h> 
#include <iostream>

#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <thrust/device_ptr.h>

//#define ntest 50000
#define ntest 1

template <unsigned int blockSize>
__device__ void reduceBlock(double *sdata, const unsigned int tid)
{
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
}

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

  reduceBlock<blockSize>(sdata, tid);

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

int nextpow2(unsigned int x)
{
  return (int)pow(2,ceil(log(x)/log(2)));
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

int main(void) 
{
  double *ha;   // Host array 
  double *da;   // Device array

  // Grid setup
  const unsigned int itot = 32;
  const unsigned int jtot = 32;
  const unsigned int ktot = 32;
  const unsigned int ng   = 1;

  unsigned int icells = itot + 2*ng;
  unsigned int jcells = jtot + 2*ng;
  unsigned int kcells = ktot + 2*ng;
  unsigned int ncells = icells * jcells * kcells;

  // Allocate host array, fill with random stuff and calculate maximum
  ha = new double[ncells];
  double cpumax = -1e12;
  for (int i = 0; i<ncells; ++i) 
  {
    ha[i] = rand() % 100000;
    if(ha[i] > cpumax)
      cpumax = ha[i];
  }
  printf("max value cpu = %e\n",cpumax);

  // Put device array on GPU
  cudaMalloc((void **)&da, ncells*sizeof(double));
  cudaMemcpy(da, ha, ncells*sizeof(double), cudaMemcpyHostToDevice);

  // Maximum value using Thrust
  // --------------------------------------
  thrust::device_ptr<double> ta = thrust::device_pointer_cast(da);
  double hmax;
  for (int i = 0; i<ntest; ++i)
  {
    hmax = thrust::reduce(ta, ta+ncells, -1.0, thrust::maximum<double>()); 
  }
  printf("max value thrust = %e\n",hmax);

  // Maximum value "own" (NVIDIA) kernel;
  // --------------------------------------
  unsigned int blocksize = 512; //1*nextpow2(kcells);
  unsigned int nblocks   = (int)ceil((double)ncells / (double)blocksize);
  //printf("nthreads=%i, nblocks=%i\n",blocksize,nblocks);

  double *idata, *dmax;
  cudaMalloc((void**)&idata, (nblocks)*sizeof(double));
  cudaMalloc((void**)&dmax, sizeof(double));
  for (int i = 0; i<ntest; ++i)
  {
    max_multipass(ncells,   blocksize,  nblocks,  da,     idata);
    max_multipass(nblocks,  blocksize,  1,        idata,  dmax);
  }
  double hmax2;
  cudaMemcpy(&hmax2, dmax, sizeof(double), cudaMemcpyDeviceToHost);
  printf("max kernel = %e\n",hmax2);

  // Cleanup
  // --------------------------------------
  cudaFree(da);
  cudaFree(idata);
  cudaFree(dmax);
  delete[] ha;
 
  return 0;
}
