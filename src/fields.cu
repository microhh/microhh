#include "fields.h"
#include "grid.h"
#include "master.h"

int cfields::prepareGPU()
{
  const int nmemsize = grid->ncellsp*sizeof(double);
  const int nmemsize2d = grid->ijcellsp*sizeof(double);

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
  const int jcells    = grid->jcells;
  const int jkcells   = grid->jcells * grid->kcells;
  const int imemsizep = grid->icellsp * sizeof(double);
  const int imemsize  = grid->icells  * sizeof(double);

  for(fieldmap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
  {
    cudaMemcpy2D(&it->second->data_g[grid->memoffset],        imemsizep, it->second->data,        imemsize, imemsize, jkcells, cudaMemcpyHostToDevice);
    cudaMemcpy2D(&it->second->databot_g[grid->memoffset],     imemsizep, it->second->databot,     imemsize, imemsize, jcells,  cudaMemcpyHostToDevice);
    cudaMemcpy2D(&it->second->datatop_g[grid->memoffset],     imemsizep, it->second->datatop,     imemsize, imemsize, jcells,  cudaMemcpyHostToDevice);
    cudaMemcpy2D(&it->second->datagradbot_g[grid->memoffset], imemsizep, it->second->datagradbot, imemsize, imemsize, jcells,  cudaMemcpyHostToDevice);
    cudaMemcpy2D(&it->second->datagradtop_g[grid->memoffset], imemsizep, it->second->datagradtop, imemsize, imemsize, jcells,  cudaMemcpyHostToDevice);
    cudaMemcpy2D(&it->second->datafluxbot_g[grid->memoffset], imemsizep, it->second->datafluxbot, imemsize, imemsize, jcells,  cudaMemcpyHostToDevice);
    cudaMemcpy2D(&it->second->datafluxtop_g[grid->memoffset], imemsizep, it->second->datafluxtop, imemsize, imemsize, jcells,  cudaMemcpyHostToDevice);
  }

  for(fieldmap::const_iterator it=at.begin(); it!=at.end(); ++it)
    cudaMemcpy2D(&it->second->data_g[grid->memoffset],        imemsizep, it->second->data,        imemsize, imemsize, jkcells, cudaMemcpyHostToDevice);

  cudaMemcpy2D(&a["p"]->data_g[grid->memoffset],              imemsizep, a["p"]->data,            imemsize, imemsize, jkcells, cudaMemcpyHostToDevice);  
  cudaMemcpy2D(&a["tmp1"]->data_g[grid->memoffset],           imemsizep, a["tmp1"]->data,         imemsize, imemsize, jkcells, cudaMemcpyHostToDevice);  
  cudaMemcpy2D(&a["tmp2"]->data_g[grid->memoffset],           imemsizep, a["tmp2"]->data,         imemsize, imemsize, jkcells, cudaMemcpyHostToDevice);  

  //printf("--> forwardGPU\n");

  return 0;
}

int cfields::backwardGPU()
{
  const int jcells    = grid->jcells;
  const int jkcells   = grid->jcells * grid->kcells;
  const int imemsizep = grid->icellsp * sizeof(double);
  const int imemsize  = grid->icells  * sizeof(double);

  for(fieldmap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
  {
    cudaMemcpy2D(it->second->data,        imemsize, &it->second->data_g[grid->memoffset],        imemsizep, imemsize, jkcells, cudaMemcpyDeviceToHost);
    cudaMemcpy2D(it->second->databot,     imemsize, &it->second->databot_g[grid->memoffset],     imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost);
    cudaMemcpy2D(it->second->datatop,     imemsize, &it->second->datatop_g[grid->memoffset],     imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost);
    cudaMemcpy2D(it->second->datagradbot, imemsize, &it->second->datagradbot_g[grid->memoffset], imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost);
    cudaMemcpy2D(it->second->datagradtop, imemsize, &it->second->datagradtop_g[grid->memoffset], imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost);
    cudaMemcpy2D(it->second->datafluxbot, imemsize, &it->second->datafluxbot_g[grid->memoffset], imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost);
    cudaMemcpy2D(it->second->datafluxtop, imemsize, &it->second->datafluxtop_g[grid->memoffset], imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost);
  }

  for(fieldmap::const_iterator it=at.begin(); it!=at.end(); ++it)
    cudaMemcpy2D(it->second->data,        imemsize, &it->second->data_g[grid->memoffset],        imemsizep, imemsize, jkcells, cudaMemcpyDeviceToHost);

  cudaMemcpy2D(a["p"]->data,              imemsize, &a["p"]->data_g[grid->memoffset],            imemsizep, imemsize, jkcells, cudaMemcpyDeviceToHost);  
  cudaMemcpy2D(a["tmp1"]->data,           imemsize, &a["tmp1"]->data_g[grid->memoffset],         imemsizep, imemsize, jkcells, cudaMemcpyDeviceToHost);  
  cudaMemcpy2D(a["tmp2"]->data,           imemsize, &a["tmp2"]->data_g[grid->memoffset],         imemsizep, imemsize, jkcells, cudaMemcpyDeviceToHost);  

  //printf("--> backwardGPU\n");

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

