#include "fields.h"
#include "grid.h"
#include "master.h"

__global__ void setid(double * __restrict__ a, int id, int jj, int kk)
{
  int ijk;
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  int k = blockIdx.z;

  ijk = i + j*jj + k*kk;
  a[ijk] = (double)id;
}

__global__ void printarray(double * __restrict__ a, int id, int jj, int kk)
{
  int ijk;
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;
  int k = blockIdx.z;

  ijk = i + j*jj + k*kk;
  if(j == 2 && k == 2) printf("CvH (%d) %d, %E\n", id, i, a[ijk]);
}

int cfields::prepareGPU()
{
  double *a_gpu;

  const int nmemsize = grid->ncells*sizeof(double);

  cudaMalloc((void**)&a_gpu, nmemsize);
  cudaMemcpy(a_gpu, sp["th"]->data, nmemsize, cudaMemcpyHostToDevice);

  // test CUDA-aware MPI
  dim3 gpugrid(grid->icells, grid->jcells, grid->kcells);
  dim3 gpublock(1,1,1);

  cudaError err;

  setid<<<gpugrid,gpublock>>>(a_gpu, master->mpiid, grid->icells, grid->ijcells);
  grid->boundary_cyclic(a_gpu);
  err = cudaGetLastError();
  if(cudaSuccess != err)
    master->printMessage("%s\n", cudaGetErrorString(err));

  printarray<<<gpugrid,gpublock>>>(a_gpu, master->mpiid, grid->icells, grid->ijcells);
  err = cudaGetLastError();
  cudaDeviceSynchronize();
  if(cudaSuccess != err)
    master->printMessage("%s\n", cudaGetErrorString(err));

  return 0;
}
