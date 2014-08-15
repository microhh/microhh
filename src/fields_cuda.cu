#include "fields.h"
#include "grid.h"
#include "master.h"

int cfields::prepareGPU()
{
  const int nmemsize = grid->ncells*sizeof(double);

  for(fieldmap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
    cudaMalloc(&it->second->data_g, nmemsize);

  for(fieldmap::const_iterator it=at.begin(); it!=at.end(); ++it)
    cudaMalloc(&it->second->data_g, nmemsize);

  // copy all the data to the GPU
  forwardGPU();

  return 0;
}

int cfields::forwardGPU()
{
  const int nmemsize = grid->ncells*sizeof(double);

  for(fieldmap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
    cudaMemcpy(it->second->data_g, it->second->data, nmemsize, cudaMemcpyHostToDevice);

  for(fieldmap::const_iterator it=at.begin(); it!=at.end(); ++it)
    cudaMemcpy(it->second->data_g, it->second->data, nmemsize, cudaMemcpyHostToDevice);

  return 0;
}

int cfields::backwardGPU()
{
  const int nmemsize = grid->ncells*sizeof(double);

  for(fieldmap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
    cudaMemcpy(it->second->data, it->second->data_g, nmemsize, cudaMemcpyDeviceToHost);

  for(fieldmap::const_iterator it=at.begin(); it!=at.end(); ++it)
    cudaMemcpy(it->second->data, it->second->data_g, nmemsize, cudaMemcpyDeviceToHost);

  return 0;
}
