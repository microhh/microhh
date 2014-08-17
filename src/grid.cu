#include "grid.h"

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
