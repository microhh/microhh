#include <cstdio>
#include "grid.h"
#include "field3d.h"

cfield3d::cfield3d(cgrid *gridin, double *dataref, std::string namein)
{
  std::printf("Creating instance of object field3d\n");
  grid = gridin;
  data = dataref;
  name = namein;
}

cfield3d::~cfield3d()
{
  std::printf("Destroying instance of object field3d\n");
}

int cfield3d::boundary_bottop(int sw)
{ 
  int ijk0,ijk1,ii,jj,kk,kstart,kend;

  kstart = grid->kstart;
  kend   = grid->kend;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  if(sw == 0)
  {
    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*jj + (kstart-k-1)*kk;
          ijk1 = i + j*jj + (kstart+k  )*kk;
          data[ijk0] = -1.*data[ijk1];
        }

    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
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
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*jj + (kstart-k-1)*kk;
          ijk1 = i + j*jj + (kstart+k  )*kk;
          data[ijk0] = data[ijk1];
        }

    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*jj + (kend+k  )*kk;
          ijk1 = i + j*jj + (kend-k-1)*kk;
          data[ijk0] = data[ijk1];
        }
  }

  // check loop
  // int j = grid->jstart+17;
  // int i = grid->istart+11;

  // for(int k=0; k<grid->kcells; k++)
  // {
  //   ijk0 = i + j*icells + k*ijcells;
  //   std::printf("k : %d, %f\n", k, data[ijk0]);
  // }

  return 0;
}

int cfield3d::boundary_cyclic()
{ 
  int ijk0,ijk1,ii,jj,kk,istart,iend,jstart,jend,igc,jgc;

  istart = grid->istart;
  iend   = grid->iend;
  igc    = grid->igc;
  jstart = grid->jstart;
  jend   = grid->jend;
  jgc    = grid->jgc;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  // east west boundaries
  for(int k=0; k<grid->kcells; k++)
    for(int j=0; j<grid->jcells; j++)
      for(int i=0; i<grid->igc; i++)
      {
        ijk0 = i          + j*jj + k*kk;
        ijk1 = iend-igc+i + j*jj + k*kk;
        data[ijk0] = data[ijk1];
      }

  for(int k=0; k<grid->kcells; k++)
    for(int j=0; j<grid->jcells; j++)
      for(int i=0; i<grid->igc; i++)
      {
        ijk0 = i+iend   + j*jj + k*kk;
        ijk1 = i+istart + j*jj + k*kk;
        data[ijk0] = data[ijk1];
      }

  // north south boundaries
  for(int k=0; k<grid->kcells; k++)
    for(int j=0; j<grid->jgc; j++)
      for(int i=0; i<grid->icells; i++)
      {
        ijk0 = i + j           *jj + k*kk;
        ijk1 = i + (jend-jgc+j)*jj + k*kk;
        data[ijk0] = data[ijk1];
      }

  for(int k=0; k<grid->kcells; k++)
    for(int j=0; j<grid->jgc; j++)
      for(int i=0; i<grid->icells; i++)
      {
        ijk0 = i + (j+jend  )*jj + k*kk;
        ijk1 = i + (j+jstart)*jj + k*kk;
        data[ijk0] = data[ijk1];
      }

  return 0;
}

int cfield3d::save(int n)
{
  FILE *pFile;
  char filename[256];
  std::sprintf(filename, "%s.%06d", name.c_str(), n);
  pFile = fopen(filename, "wb");
  fwrite(data, sizeof(double), grid->ncells, pFile);
  fclose(pFile);

  return 0;
}

int cfield3d::load(int n)
{
  FILE *pFile;
  char filename[256];
  std::sprintf(filename, "%s.%06d", name.c_str(), n);
  pFile = fopen(filename, "rb");
  fread(data, sizeof(double), grid->ncells, pFile);
  fclose(pFile);

  return 0;
}
