#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "mpicheck.h"

cmpicheck::cmpicheck(cgrid *gridin, cmpi *mpiin)
{
  std::printf("Creating instance of object mpicheck\n");
  grid   = gridin;
  mpi    = mpiin;
}

cmpicheck::~cmpicheck()
{
  std::printf("Destroying instance of object mpicheck\n");
}

int cmpicheck::showLayout()
{
  
  std::printf("MPI id, mpicoordx, mpicoordy, neast, nwest, nnorth, nsouth, nprocs: %2d, %2d, %2d, %2d, %2d, %2d, %2d, %2d\n",
    mpi->mpiid, mpi->mpicoordx, mpi->mpicoordy, mpi->neast, mpi->nwest, mpi->nnorth, mpi->nsouth, mpi->nprocs);

  return 0;
}

int cmpicheck::create()
{
  s = new cfield3d(grid, "s" );

  s->init();

  for(int n=0; n<grid->ncells; n++)
    s->data[n] = (double)mpi->mpiid;
  
  return 0;
}

int cmpicheck::boundary()
{
  mpi->boundary_cyclic(s->data, grid);

  return 0;
}

int cmpicheck::showLine()
{
  int i,j,k;
  int ijk,ii,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(k=0; k<grid->kcells; k++)
    for(j=0; j<grid->jcells; j++)
      for(i=0; i<grid->icells; i++)
      {
        ijk = i + j*jj + k*kk;
        std::printf("MPI id %d, s(%d,%d,%d) = %4.0f\n", mpi->mpiid, i, j, k, s->data[ijk]);
      }

  return 0;
}
