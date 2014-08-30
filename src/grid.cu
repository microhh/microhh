#include "grid.h"

__global__ void grid_cyclic(double * __restrict__ data, 
                            const int icells, const int jcells, const int kcells,
                            const int istart, const int jstart,
                            const int iend,   const int jend, 
                            const int igc,    const int jgc)
{
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  const int j = blockIdx.y*blockDim.y + threadIdx.y;
  const int k = blockIdx.z;

  const int jj = icells;
  const int kk = icells*jcells;

  int ijk0 = i          + j*jj + k*kk;
  int ijk1 = iend-igc+i + j*jj + k*kk;

  // East-west
  if(k < kcells && j < jcells && i < igc)
  {
    ijk0 = i          + j*jj + k*kk;
    ijk1 = iend-igc+i + j*jj + k*kk;
    data[ijk0] = data[ijk1];
  }

  if(k < kcells && j < jcells && i < igc)
  {
    ijk0 = i+iend   + j*jj + k*kk;
    ijk1 = i+istart + j*jj + k*kk;
    data[ijk0] = data[ijk1];
  }

  // North-east
  if(k < kcells && j < jgc && i < icells)
  {
    ijk0 = i + j           *jj + k*kk;
    ijk1 = i + (jend-jgc+j)*jj + k*kk;
    data[ijk0] = data[ijk1];
  }

  if(k < kcells && j < jgc && i < icells)
  {
    ijk0 = i + (j+jend  )*jj + k*kk;
    ijk1 = i + (j+jstart)*jj + k*kk;
    data[ijk0] = data[ijk1];
  }

}

int cgrid::prepareGPU()
{
  const int kmemsize = kcells*sizeof(double);

  cudaMalloc((void**)&dz_g   , kmemsize);
  cudaMalloc((void**)&dzh_g  , kmemsize);
  cudaMalloc((void**)&dzi_g  , kmemsize);
  cudaMalloc((void**)&dzhi_g , kmemsize);
  cudaMalloc((void**)&dzi4_g , kmemsize);
  cudaMalloc((void**)&dzhi4_g, kmemsize);

  cudaMemcpy(dz_g   , dz   , kmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dzh_g  , dzh  , kmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dzi_g  , dzi  , kmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dzhi_g , dzhi , kmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dzi4_g , dzi4 , kmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(dzhi4_g, dzhi4, kmemsize, cudaMemcpyHostToDevice);

  return 0;
}

int cgrid::boundary_cyclic_gpu(double * data)
{
  const int blocki  = 128;
  const int blockj  = 2;
  const int gridi   = icells/blocki + (icells%blocki > 0);
  const int gridj   = jcells/blockj + (jcells%blockj > 0);

  dim3 gridGPU (gridi, gridj, kcells);
  dim3 blockGPU(blocki, blockj, 1);

  grid_cyclic<<<gridGPU,blockGPU>>>(data, icells, jcells, kcells,
                                    istart, jstart,
                                    iend,   jend,
                                    igc,    jgc);
 
  return 0;
}
