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

int cmpicheck::checkLayout()
{
  
  std::printf("MPI id, mpicoordx, mpicoordy, neast, nwest, nnorth, nsouth, nprocs: %2d, %2d, %2d, %2d, %2d, %2d, %2d, %2d\n",
    mpi->mpiid, mpi->mpicoordx, mpi->mpicoordy, mpi->neast, mpi->nwest, mpi->nnorth, mpi->nsouth, mpi->nprocs);

  return 0;
}

int cmpicheck::create()
{
  s     = new cfield3d(grid, "s");
  temp1 = new cfield3d(grid, "temp1");
  temp2 = new cfield3d(grid, "temp2");

  s->init();
  temp1->init();
  temp2->init();

  for(int n=0; n<grid->ncells; n++)
    s->data[n] = (double)mpi->mpiid;
  
  return 0;
}

int cmpicheck::checkBoundary()
{
  mpi->boundary_cyclic(s->data);

  int i,j,k;
  int ijk,ii,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  k = grid->kstart;
  j = grid->jstart;
  for(i=0; i<grid->icells; i++)
  {
    ijk = i + j*jj + k*kk;
    std::printf("MPI i-line id %d, s(%d,%d,%d) = %4.0f\n", mpi->mpiid, i, j, k, s->data[ijk]);
  }

  k = grid->kstart;
  i = grid->istart;
  for(j=0; j<grid->jcells; j++)
  {
    ijk = i + j*jj + k*kk;
    std::printf("MPI j-line id %d, s(%d,%d,%d) = %4.0f\n", mpi->mpiid, i, j, k, s->data[ijk]);
  }

  return 0;
}

int cmpicheck::checkTranspose()
{
  int i,j,k,ijk,ijkw;
  int jj,kk,jjw,kkw;

  jj  = grid->icells;
  kk  = grid->icells*grid->jcells;
  jjw = grid->imax;
  kkw = grid->imax*grid->jmax;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj  + k*kk;
        ijkw = i + j*jjw + k*kkw;

        temp1->data[ijkw] = s->data[ijk];
      }

  mpi->transposezx(temp1->data, temp2->data);

  jj = grid->imax;
  kk = grid->imax*grid->jmax;

  k = 0;
  j = 0;

  for(i=0; i<grid->itot; i++)
  {
    ijk = i + j*jj + k*kk;
    std::printf("MPI transx id %d, s(%d,%d,%d) = %4.0f\n", mpi->mpiid, i, j, k, temp2->data[ijk]);
  }

  return 0;
}

