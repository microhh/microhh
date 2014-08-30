#include <stdio.h>
#include "tools.h"

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

// Get the maximum value from GPU field "field" of size "nvalues"
double maximum_gpu(double * field, int nvalues)
{
  unsigned int blocksize = 256; 
  unsigned int nblocks = (int)ceil((double)nvalues / (double)blocksize); 

  double *idata; // Array holding max values from first reduction step
  double *maxval_g; // Pointer holding max value from second reduction step
  double maxval;  

  // TODO: move to some init function...
  cudaMalloc((void**)&idata, (nblocks)*sizeof(double));
  cudaMalloc((void**)&maxval_g, sizeof(double));

  // Do the reduction in two steps
  max_multipass(nvalues, blocksize, nblocks, field, idata);
  max_multipass(nblocks, blocksize, 1,       idata, maxval_g);

  // Get max value from GPU
  cudaMemcpy(&maxval, maxval_g, sizeof(double), cudaMemcpyDeviceToHost);

  // Cleanup GPU 
  cudaFree(idata);
  cudaFree(maxval_g);

  return maxval;
}


// CUDA error checking. 
void CudaCheckError()
{
  cudaError err = cudaGetLastError();
  if(cudaSuccess != err)
    printf("CUDA error : %s\n",cudaGetErrorString(err));

  //err = cudaDeviceSynchronize();
  //if(cudaSuccess != err)
  //{
  //  fprintf( stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",
  //           file, line, cudaGetErrorString( err ) );
  //  exit( -1 );
  //}
 
  return;
}
