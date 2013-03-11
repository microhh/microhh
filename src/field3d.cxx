#include <cstdio>
#include <iostream>
#include "grid.h"
#include "field3d.h"
#include "defines.h"

cfield3d::cfield3d(cgrid *gridin, cmpi *mpiin, std::string namein)
{
  std::printf("Creating instance of object field3d\n");
  grid = gridin;
  name = namein;
  mpi  = mpiin;
}

cfield3d::~cfield3d()
{
  if(allocated)
  {
    delete[] data;
    delete[] databot;
    delete[] datatop;
    delete[] datagradbot;
    delete[] datagradtop;
  }

  // std::printf("Destroying instance of object field3d\n");
}

int cfield3d::init()
{
  // allocate the memory
  if(mpi->mpiid == 0) std::printf("Allocating %d bytes of memory for %s\n", grid->ncells*(int)sizeof(double), name.c_str());
  data    = new double[grid->ncells];

  // allocate the boundary cells
  databot = new double[grid->icells*grid->jcells];
  datatop = new double[grid->icells*grid->jcells];
  datagradbot = new double[grid->icells*grid->jcells];
  datagradtop = new double[grid->icells*grid->jcells];

  allocated = true;

  // set all values to zero
  for(int n=0; n<grid->ncells; n++)
    data[n] = 0.;

  for(int n=0; n<grid->icells*grid->jcells; n++)
  {
    databot    [n] = 0.;
    datatop    [n] = 0.;
    datagradbot[n] = 0.;
    datagradtop[n] = 0.;
  }

  return 0;
}

/*
int cfield3d::boundary_bottop(int sw)
{ 
  int ijk0,ijk1,jj,kk,kstart,kend;

  kstart = grid->kstart;
  kend   = grid->kend;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  if(sw == 0)
  {
    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
#pragma ivdep
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*jj + (kstart-k-1)*kk;
          ijk1 = i + j*jj + (kstart+k  )*kk;
          data[ijk0] = -1.*data[ijk1];
        }

    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
#pragma ivdep
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*jj + (kend+k  )*kk;
          ijk1 = i + j*jj + (kend-k-1)*kk;
          data[ijk0] = -1.*data[ijk1];
        }
  }
  else if(sw == 1)
  {
    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
#pragma ivdep
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*jj + (kstart-k-1)*kk;
          ijk1 = i + j*jj + (kstart+k  )*kk;
          data[ijk0] = data[ijk1];
        }

    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
#pragma ivdep
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*jj + (kend+k  )*kk;
          ijk1 = i + j*jj + (kend-k-1)*kk;
          data[ijk0] = data[ijk1];
        }
  }

  return 0;
}

int cfield3d::boundary_cyclic()
{ 
  int ijk0,ijk1,jj,kk,istart,iend,jstart,jend,igc,jgc;

  istart = grid->istart;
  iend   = grid->iend;
  igc    = grid->igc;
  jstart = grid->jstart;
  jend   = grid->jend;
  jgc    = grid->jgc;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  // east west boundaries
  for(int k=0; k<grid->kcells; k++)
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->igc; i++)
      {
        ijk0 = i          + j*jj + k*kk;
        ijk1 = iend-igc+i + j*jj + k*kk;
        data[ijk0] = data[ijk1];
      }

  for(int k=0; k<grid->kcells; k++)
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->igc; i++)
      {
        ijk0 = i+iend   + j*jj + k*kk;
        ijk1 = i+istart + j*jj + k*kk;
        data[ijk0] = data[ijk1];
      }

  // north south boundaries
  for(int k=0; k<grid->kcells; k++)
    for(int j=0; j<grid->jgc; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ijk0 = i + j           *jj + k*kk;
        ijk1 = i + (jend-jgc+j)*jj + k*kk;
        data[ijk0] = data[ijk1];
      }

  for(int k=0; k<grid->kcells; k++)
    for(int j=0; j<grid->jgc; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ijk0 = i + (j+jend  )*jj + k*kk;
        ijk1 = i + (j+jstart)*jj + k*kk;
        data[ijk0] = data[ijk1];
      }

  return 0;
}*/

int cfield3d::save(int n, double * restrict tmp1, double * restrict tmp2)
{
  char filename[256];

  /* 
  std::sprintf(filename, "%s.%07d.%07d", name.c_str(), n, mpiid);
  FILE *pFile;
  pFile = fopen(filename, "wb");

  if(pFile == NULL)
  {
    std::printf("ERROR \"%s\" cannot be written", filename);
    return 1;
  }
  else
    std::printf("Saving \"%s\"\n", filename);

  int ijk,istart,jj,kk;

  istart = grid->istart;
  jj     = grid->icells;
  kk     = grid->icells*grid->jcells;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      {
        ijk = istart + j*jj + k*kk;
        fwrite(&data[ijk], sizeof(double), grid->imax, pFile);
      }

  fclose(pFile);
  */

  std::sprintf(filename, "%s.%07d", name.c_str(), n);

  if(mpi->mpiid == 0) std::printf("Saving \"%s\"\n", filename);

  if(grid->savefield3d(data, tmp1, tmp2, filename))
  {
    if(mpi->mpiid == 0) std::printf("ERROR \"%s\" cannot be written\n", filename);
    return 1;
  }

  return 0;
}

int cfield3d::load(int n, double * restrict tmp1, double * restrict tmp2)
{
  char filename[256];

  /*
  FILE *pFile;
  std::sprintf(filename, "%s.%07d.%07d", name.c_str(), n, mpiid);
  pFile = fopen(filename, "rb");

  if(pFile == NULL)
  {
    std::printf("ERROR \"%s\" does not exist\n", filename);
    return 1;
  }
  else
    std::printf("Loading \"%s\"\n", filename);

  int ijk,istart,jj,kk;

  istart = grid->istart;
  jj     = grid->icells;
  kk     = grid->icells*grid->jcells;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      {
        ijk = istart + j*jj + k*kk;
        fread(&data[ijk], sizeof(double), grid->imax, pFile);
      }

  fclose(pFile);
  */

  std::sprintf(filename, "%s.%07d", name.c_str(), n);

  if(mpi->mpiid == 0) std::printf("Loading \"%s\"\n", filename);

  if(grid->loadfield3d(data, tmp1, tmp2, filename))
  {
    if(mpi->mpiid == 0) std::printf("ERROR \"%s\" does not exist\n", filename);
    return 1;
  }

  return 0;
}
