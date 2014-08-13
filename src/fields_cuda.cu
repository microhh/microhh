#include "fields.h"
#include "grid.h"

int cfields::prepareGPU()
{
  double *a;

  const int nmemsize = grid->ncells*sizeof(double);

  cudaMalloc((void**)&a, nmemsize);
  cudaMemcpy(a, sp["th"]->data, nmemsize, cudaMemcpyHostToDevice);

  // test CUDA-aware MPI
  grid->boundary_cyclic(a);

  return 0;
}
