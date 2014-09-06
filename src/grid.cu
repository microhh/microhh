#include "grid.h"

__global__ void grid_cyclic_x(double * __restrict__ data, 
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

  int ijk0, ijk1, ijk2, ijk3;

  // East-west
  if(k < kcells && j < jcells && i < igc)
  {
    ijk0 = i          + j*jj + k*kk;
    ijk1 = iend-igc+i + j*jj + k*kk;
    ijk2 = i+iend     + j*jj + k*kk;
    ijk3 = i+istart   + j*jj + k*kk;

    data[ijk0] = data[ijk1];
    data[ijk2] = data[ijk3];
  }
}

__global__ void grid_cyclic_y(double * __restrict__ data, 
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

  int ijk0, ijk1, ijk2, ijk3;

  // North-south
  if(k < kcells && j < jgc && i < icells)
  {
    ijk0 = i + j           *jj + k*kk;
    ijk1 = i + (jend-jgc+j)*jj + k*kk;
    ijk2 = i + (j+jend  )*jj + k*kk;
    ijk3 = i + (j+jstart)*jj + k*kk;
    data[ijk0] = data[ijk1];
    data[ijk2] = data[ijk3];
  }
}

int cgrid::prepareGPU()
{
  /* Align the interior of the grid (i.e. excluding ghost cells) with 
     the 128 byte memory blocks of the GPU's global memory */
  memoffset = 16 - igc;           // Padding at start of array 
  int padl  = 16-(int)itot%16;    // Elements left in last 128 byte block
  icellsp   = itot + padl + (padl < 2*igc) * 16;
  ijcellsp  = icellsp * jcells;  
  ncellsp   = ijcellsp * kcells;

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

int cgrid::clearGPU()
{
  cudaFree(dz_g   );
  cudaFree(dzh_g  );
  cudaFree(dzi_g  );
  cudaFree(dzhi_g );
  cudaFree(dzi4_g );
  cudaFree(dzhi4_g);

  return 0;
}

int cgrid::boundary_cyclic_gpu(double * data)
{
  const int blocki_x = igc;
  const int blockj_x = 256 / igc + (256%igc > 0);
  const int gridi_x  = 1;
  const int gridj_x  = jcells/blockj_x + (jcells%blockj_x > 0);

  const int blocki_y = 256 / igc + (256%jgc > 0);
  const int blockj_y = jgc;
  const int gridi_y  = icells/blocki_y + (icells%blocki_y > 0);
  const int gridj_y  = 1;

  dim3 gridGPUx (gridi_x, gridj_x, kcells);
  dim3 blockGPUx(blocki_x, blockj_x, 1);

  dim3 gridGPUy (gridi_y, gridj_y, kcells);
  dim3 blockGPUy(blocki_y, blockj_y, 1);

  grid_cyclic_x<<<gridGPUx,blockGPUx>>>(data, icells, jcells, kcells,
                                        istart, jstart,
                                        iend,   jend,
                                        igc,    jgc);

  grid_cyclic_y<<<gridGPUy,blockGPUy>>>(data, icells, jcells, kcells,
                                        istart, jstart,
                                        iend,   jend,
                                        igc,    jgc);

  return 0;
}
