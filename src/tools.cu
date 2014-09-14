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

#include <stdio.h>
#include "float.h"
#include "tools.h"

#define MAXTHREADS 512 // Maximum number of threads used in reduce algoritms
#define SUM 0
#define MAX 1
#define MIN 2

int nextpow2(unsigned int x)
{
  return (int)pow(2,ceil(log(x)/log(2)));
}

template <int func>
__device__ double reduction(double v1, double v2)
{
  double rval;
  if (func == SUM)
    rval = v1+v2;
  else if (func == MAX)
    rval = fmax(v1,v2);
  else if (func == MIN)
    rval = fmin(v1,v2);
  return rval;
} 

// Reduce one block of data
template <int func, int blockSize> 
__device__ void reduceBlock(volatile double *as, const unsigned int tid)
{
  /* Loop is completely unrolled for performance */
  if (blockSize >= 512) { if (tid < 256) { as[tid] = reduction<func>(as[tid],as[tid + 256]); } __syncthreads(); }
  if (blockSize >= 256) { if (tid < 128) { as[tid] = reduction<func>(as[tid],as[tid + 128]); } __syncthreads(); }
  if (blockSize >= 128) { if (tid <  64) { as[tid] = reduction<func>(as[tid],as[tid +  64]); } __syncthreads(); }

  /* Once we get to the last 32 values (1 thread warp), the __syncthreads() is no longer necessary */
  if (tid < 32)
  {
    if (blockSize >=  64) { if (tid < 32) { as[tid] = reduction<func>(as[tid],as[tid + 32]); }}
    if (blockSize >=  32) { if (tid < 16) { as[tid] = reduction<func>(as[tid],as[tid + 16]); }}
    if (blockSize >=  16) { if (tid <  8) { as[tid] = reduction<func>(as[tid],as[tid +  8]); }}
    if (blockSize >=   8) { if (tid <  4) { as[tid] = reduction<func>(as[tid],as[tid +  4]); }}
    if (blockSize >=   4) { if (tid <  2) { as[tid] = reduction<func>(as[tid],as[tid +  2]); }}
    if (blockSize >=   2) { if (tid <  1) { as[tid] = reduction<func>(as[tid],as[tid +  1]); }}
  }
}

// Reduce field from 3D to 2D, excluding ghost cells and padding
template <int func, int blockSize> 
__global__ void deviceReduceInterior(const double *a, double *a2d, 
                                     unsigned int istart, unsigned int jstart, unsigned int kstart, 
                                     unsigned int iend,   unsigned int jend,   
                                     unsigned int icells, unsigned int ijcells)
{
  extern __shared__ double as[];

  unsigned int tid  = threadIdx.x;
  unsigned int i    = istart + threadIdx.x;
  unsigned int j    = jstart + blockIdx.y;
  unsigned int k    = kstart + blockIdx.z; 
  unsigned int jk   = blockIdx.y+blockIdx.z*(jend-jstart);   // Index in 2D "a2d"
  unsigned int ijk  = i + j*icells + k*ijcells;              // Index in 3D "a"
  unsigned int ijkm = ijkm = iend + j*icells + k*ijcells;    // Max index in X-direction

  double tmpval;
  if (func == MAX)
    tmpval = -DBL_MAX;
  else if (func == MIN)
    tmpval = DBL_MAX;
  else if (func == SUM)
    tmpval = 0;
  
  int ii = ijk;
  while (ii < ijkm)
  {
    tmpval = reduction<func>(tmpval,a[ii]);
    if(ii + blockDim.x < ijkm)
      tmpval = reduction<func>(tmpval,a[ii+blockDim.x]);
    ii += 2*blockDim.x;
  }
  as[tid] = tmpval;

  __syncthreads();

  reduceBlock<func, blockSize>(as, tid);

  if (tid == 0)
    a2d[jk] = as[0];
}

// Reduce array, not accounting from ghost cells or padding 
template <int func, int blockSize> 
__global__ void deviceReduceAll(const double *a, double *aout, unsigned int ncells, unsigned int nvaluesperblock)  
{
  extern __shared__ double as[];

  unsigned int tid  = threadIdx.x;
  unsigned int ii   = nvaluesperblock *  blockIdx.x + threadIdx.x;
  unsigned int iim  = nvaluesperblock * (blockIdx.x+1);

  double tmpval;
  if (func == MAX)
    tmpval = -DBL_MAX;
  else if (func == MIN)
    tmpval = DBL_MAX;
  else if (func == SUM)
    tmpval = 0;
  
  while (ii < iim)
  {
    tmpval = reduction<func>(tmpval,a[ii]);
    if(ii + blockDim.x < iim && ii + blockDim.x < ncells)
      tmpval = reduction<func>(tmpval,a[ii+blockDim.x]);
    ii += 2*blockDim.x;
  }
  as[tid] = tmpval;

  /* Make sure all threads are synchronised before reducing the shared array */
  __syncthreads();

  /* Reduce block in shared memory */
  reduceBlock<func, blockSize>(as, tid);

  /* First value in shared array now holds the reduced value. Write back to global memory */
  if (tid == 0)
    aout[blockIdx.x] = as[0];
}

void reduceInterior(double *a, double *a2d, 
                    int itot, int istart, int iend,
                    int jtot, int jstart, int jend,
                    int ktot, int kstart, int kend,
                    int icells, int ijcells)
{
  int nthreads = min(MAXTHREADS, nextpow2(itot/2));
  dim3 gridGPU (1, jtot, ktot);
  dim3 blockGPU(nthreads, 1, 1);

  // HACK BVS: for now mode hardcoded at MAX for testing......
  switch (nthreads)
  {
    case 512:
      deviceReduceInterior<MAX, 512><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, a2d, istart, jstart, kstart, iend, jend, icells, ijcells); break;
    case 256:
      deviceReduceInterior<MAX, 256><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, a2d, istart, jstart, kstart, iend, jend, icells, ijcells); break;
    case 128:
      deviceReduceInterior<MAX, 128><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, a2d, istart, jstart, kstart, iend, jend, icells, ijcells); break;
    case 64:
      deviceReduceInterior<MAX,  64><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, a2d, istart, jstart, kstart, iend, jend, icells, ijcells); break;
    case 32:
      deviceReduceInterior<MAX,  32><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, a2d, istart, jstart, kstart, iend, jend, icells, ijcells); break;
    case 16:
      deviceReduceInterior<MAX,  16><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, a2d, istart, jstart, kstart, iend, jend, icells, ijcells); break;
    case 8:
      deviceReduceInterior<MAX,   8><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, a2d, istart, jstart, kstart, iend, jend, icells, ijcells); break;
    case 4:
      deviceReduceInterior<MAX,   4><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, a2d, istart, jstart, kstart, iend, jend, icells, ijcells); break;
    case 2:
      deviceReduceInterior<MAX,   2><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, a2d, istart, jstart, kstart, iend, jend, icells, ijcells); break;
    case 1:
      deviceReduceInterior<MAX,   1><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, a2d, istart, jstart, kstart, iend, jend, icells, ijcells); break;
  }
}

void reduceAll(double *a, double *aout, int ncells, int nblocks, int nvaluesperblock)
{
  int nthreads = min(MAXTHREADS, nextpow2(nvaluesperblock/2));
  dim3 gridGPU (nblocks,  1, 1);
  dim3 blockGPU(nthreads, 1, 1);

  // HACK BVS: for now mode hardcoded at MAX for testing......
  switch (nthreads)
  {
    case 512:
      deviceReduceAll<MAX, 512><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, aout, ncells, nvaluesperblock); break;
    case 256:
      deviceReduceAll<MAX, 256><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, aout, ncells, nvaluesperblock); break;
    case 128:
      deviceReduceAll<MAX, 128><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, aout, ncells, nvaluesperblock); break;
    case 64:
      deviceReduceAll<MAX,  64><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, aout, ncells, nvaluesperblock); break;
    case 32:
      deviceReduceAll<MAX,  32><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, aout, ncells, nvaluesperblock); break;
    case 16:
      deviceReduceAll<MAX,  16><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, aout, ncells, nvaluesperblock); break;
    case 8:
      deviceReduceAll<MAX,   8><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, aout, ncells, nvaluesperblock); break;
    case 4:
      deviceReduceAll<MAX,   4><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, aout, ncells, nvaluesperblock); break;
    case 2:
      deviceReduceAll<MAX,   2><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, aout, ncells, nvaluesperblock); break;
    case 1:
      deviceReduceAll<MAX,   1><<<gridGPU, blockGPU, nthreads*sizeof(double)>>>(a, aout, ncells, nvaluesperblock); break;
  }
}

// CUDA error checking. 
void CudaCheckError()
{
  cudaError err = cudaGetLastError();
  if(cudaSuccess != err)
    printf("CUDA error : %s\n",cudaGetErrorString(err));

  err = cudaDeviceSynchronize();
  if(cudaSuccess != err)
    printf("CUDA error with sync : %s\n",cudaGetErrorString(err));
 
  return;
}
