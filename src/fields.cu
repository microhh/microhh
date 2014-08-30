#include "fields.h"
#include "grid.h"
#include "master.h"

int cfields::prepareGPU()
{
  const int nmemsize = grid->ncells*sizeof(double);
  const int nmemsize2d = grid->ijcells*sizeof(double);

  for(fieldmap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
  {
    cudaMalloc(&it->second->data_g,        nmemsize);
    cudaMalloc(&it->second->databot_g,     nmemsize2d);
    cudaMalloc(&it->second->datatop_g,     nmemsize2d);
    cudaMalloc(&it->second->datagradbot_g, nmemsize2d);
    cudaMalloc(&it->second->datagradtop_g, nmemsize2d);
    cudaMalloc(&it->second->datafluxbot_g, nmemsize2d);
    cudaMalloc(&it->second->datafluxtop_g, nmemsize2d);
  }

  for(fieldmap::const_iterator it=at.begin(); it!=at.end(); ++it)
  {
    cudaMalloc(&it->second->data_g, nmemsize);
  }

  cudaMalloc(&a["p"]->data_g, nmemsize);
  cudaMalloc(&a["tmp1"]->data_g, nmemsize);
  cudaMalloc(&a["tmp2"]->data_g, nmemsize);

  // copy all the data to the GPU
  forwardGPU();

  return 0;
}

int cfields::forwardGPU()
{
  const int nmemsize = grid->ncells*sizeof(double);
  const int nmemsize2d = grid->ijcells*sizeof(double);

  for(fieldmap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
  {
    cudaMemcpy(it->second->data_g,        it->second->data,        nmemsize,   cudaMemcpyHostToDevice);
    cudaMemcpy(it->second->databot_g,     it->second->databot,     nmemsize2d, cudaMemcpyHostToDevice);
    cudaMemcpy(it->second->datatop_g,     it->second->datatop,     nmemsize2d, cudaMemcpyHostToDevice);
    cudaMemcpy(it->second->datagradbot_g, it->second->datagradbot, nmemsize2d, cudaMemcpyHostToDevice);
    cudaMemcpy(it->second->datagradtop_g, it->second->datagradtop, nmemsize2d, cudaMemcpyHostToDevice);
    cudaMemcpy(it->second->datafluxbot_g, it->second->datafluxbot, nmemsize2d, cudaMemcpyHostToDevice);
    cudaMemcpy(it->second->datafluxtop_g, it->second->datafluxtop, nmemsize2d, cudaMemcpyHostToDevice);
  }

  for(fieldmap::const_iterator it=at.begin(); it!=at.end(); ++it)
    cudaMemcpy(it->second->data_g, it->second->data, nmemsize, cudaMemcpyHostToDevice);

  cudaMemcpy(a["p"]->data_g, a["p"]->data, nmemsize, cudaMemcpyHostToDevice);  
  cudaMemcpy(a["tmp1"]->data_g, a["tmp1"]->data, nmemsize, cudaMemcpyHostToDevice);  
  cudaMemcpy(a["tmp2"]->data_g, a["tmp2"]->data, nmemsize, cudaMemcpyHostToDevice);  

  printf("--> forwardGPU\n");

  return 0;
}

int cfields::backwardGPU()
{
  const int nmemsize = grid->ncells*sizeof(double);
  const int nmemsize2d = grid->ijcells*sizeof(double);

  for(fieldmap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
  {
    cudaMemcpy(it->second->data,        it->second->data_g,        nmemsize,   cudaMemcpyDeviceToHost);
    cudaMemcpy(it->second->databot,     it->second->databot_g,     nmemsize2d, cudaMemcpyDeviceToHost);
    cudaMemcpy(it->second->datatop,     it->second->datatop_g,     nmemsize2d, cudaMemcpyDeviceToHost);
    cudaMemcpy(it->second->datagradbot, it->second->datagradbot_g, nmemsize2d, cudaMemcpyDeviceToHost);
    cudaMemcpy(it->second->datagradtop, it->second->datagradtop_g, nmemsize2d, cudaMemcpyDeviceToHost);
    cudaMemcpy(it->second->datafluxbot, it->second->datafluxbot_g, nmemsize2d, cudaMemcpyDeviceToHost);
    cudaMemcpy(it->second->datafluxtop, it->second->datafluxtop_g, nmemsize2d, cudaMemcpyDeviceToHost);
  }

  for(fieldmap::const_iterator it=at.begin(); it!=at.end(); ++it)
    cudaMemcpy(it->second->data, it->second->data_g, nmemsize, cudaMemcpyDeviceToHost);

  cudaMemcpy(a["p"]->data, a["p"]->data_g, nmemsize, cudaMemcpyDeviceToHost);  
  cudaMemcpy(a["tmp1"]->data, a["tmp1"]->data_g, nmemsize, cudaMemcpyDeviceToHost);  
  cudaMemcpy(a["tmp2"]->data, a["tmp2"]->data_g, nmemsize, cudaMemcpyDeviceToHost);  

  printf("--> backwardGPU\n");

  return 0;
}

int cfields::clearGPU()
{
  for(fieldmap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
  {
    cudaFree(&it->second->data_g);
    cudaFree(&it->second->databot_g);
    cudaFree(&it->second->datatop_g);
    cudaFree(&it->second->datagradbot_g);
    cudaFree(&it->second->datagradtop_g);
    cudaFree(&it->second->datafluxbot_g);
    cudaFree(&it->second->datafluxtop_g);
  }

  for(fieldmap::const_iterator it=at.begin(); it!=at.end(); ++it)
  {
    cudaFree(&it->second->data_g);
  }

  cudaFree(&a["p"]->data_g);
  cudaFree(&a["tmp1"]->data_g);
  cudaFree(&a["tmp2"]->data_g);

  return 0;
}

