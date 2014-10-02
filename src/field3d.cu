#include "field3d.h"
#include "grid.h"
#include "master.h"

#ifdef USECUDA
cfield3d::~cfield3d()
{
  cudaFreeHost(data);
  cudaFreeHost(databot);
  cudaFreeHost(datatop);
  cudaFreeHost(datagradbot);
  cudaFreeHost(datagradtop);
  cudaFreeHost(datafluxbot);
  cudaFreeHost(datafluxtop);
  cudaFreeHost(datamean);
}

int cfield3d::init()
{
  // allocate the memory
  master->printMessage("Allocating %d bytes of memory for %s\n", grid->ncells*(int)sizeof(double), name.c_str());

  const int ijksize = grid->ncells *sizeof(double);
  const int ijsize  = grid->ijcells*sizeof(double);
  const int ksize   = grid->kcells *sizeof(double);

  // Allocate the 3d field.
  cudaMallocHost(&data, ijksize);

  // Allocate the boundary cells.
  cudaMallocHost(&databot, ijsize);
  cudaMallocHost(&datatop, ijsize);
  cudaMallocHost(&datagradbot, ijsize);
  cudaMallocHost(&datagradtop, ijsize);
  cudaMallocHost(&datafluxbot, ijsize);
  cudaMallocHost(&datafluxtop, ijsize);
  cudaMallocHost(&datamean, ksize);

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
