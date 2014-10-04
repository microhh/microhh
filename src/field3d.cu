#include "field3d.h"
#include "grid.h"
#include "master.h"
#include "tools.h"

#ifdef USECUDA
Field3d::~Field3d()
{
  cudaSafeCall(cudaFreeHost(data));
  cudaSafeCall(cudaFreeHost(databot));
  cudaSafeCall(cudaFreeHost(datatop));
  cudaSafeCall(cudaFreeHost(datagradbot));
  cudaSafeCall(cudaFreeHost(datagradtop));
  cudaSafeCall(cudaFreeHost(datafluxbot));
  cudaSafeCall(cudaFreeHost(datafluxtop));
  cudaSafeCall(cudaFreeHost(datamean));
}

int Field3d::init()
{
  // allocate the memory
  master->printMessage("Allocating %d bytes of memory for %s\n", grid->ncells*(int)sizeof(double), name.c_str());

  const int ijksize = grid->ncells *sizeof(double);
  const int ijsize  = grid->ijcells*sizeof(double);
  const int ksize   = grid->kcells *sizeof(double);

  // Allocate the 3d field.
  cudaSafeCall(cudaMallocHost(&data, ijksize));

  // Allocate the boundary cells.
  cudaSafeCall(cudaMallocHost(&databot, ijsize));
  cudaSafeCall(cudaMallocHost(&datatop, ijsize));
  cudaSafeCall(cudaMallocHost(&datagradbot, ijsize));
  cudaSafeCall(cudaMallocHost(&datagradtop, ijsize));
  cudaSafeCall(cudaMallocHost(&datafluxbot, ijsize));
  cudaSafeCall(cudaMallocHost(&datafluxtop, ijsize));
  cudaSafeCall(cudaMallocHost(&datamean, ksize));

  // Set all values to zero
  for(int n=0; n<grid->ncells; n++)
    data[n] = 0.;

  for(int n=0; n<grid->kcells; n++)
    datamean[n] = 0.;

  for(int n=0; n<grid->icells*grid->jcells; n++)
  {
    databot    [n] = 0.;
    datatop    [n] = 0.;
    datagradbot[n] = 0.;
    datagradtop[n] = 0.;
    datafluxbot[n] = 0.;
    datafluxtop[n] = 0.;
  }

  return 0;
}
#endif

void Field3d::initDevice()
{
  const int nmemsize   = grid->ncellsp*sizeof(double);
  const int nmemsize1d = grid->kcells *sizeof(double);
  const int nmemsize2d = (grid->ijcellsp+grid->memoffset)*sizeof(double);

  cudaSafeCall(cudaMalloc(&data_g,        nmemsize  ));
  cudaSafeCall(cudaMalloc(&databot_g,     nmemsize2d));
  cudaSafeCall(cudaMalloc(&datatop_g,     nmemsize2d));
  cudaSafeCall(cudaMalloc(&datagradbot_g, nmemsize2d));
  cudaSafeCall(cudaMalloc(&datagradtop_g, nmemsize2d));
  cudaSafeCall(cudaMalloc(&datafluxbot_g, nmemsize2d));
  cudaSafeCall(cudaMalloc(&datafluxtop_g, nmemsize2d));
  cudaSafeCall(cudaMalloc(&datamean_g,    nmemsize1d));
}
