#include <cstdio>
#include "grid.h"
#include "field3d.h"

cfield3d::cfield3d(cgrid *gridin, double *dataref, std::string namein)
{
  std::printf("Creating instance of object field3d\n");
  grid = gridin;
  data = dataref;
  name = new std::string(namein);
}

cfield3d::~cfield3d()
{
  delete name;
  std::printf("Destroying instance of object field3d\n");
}

int cfield3d::boundary_bottop(int sw)
{ 
  int ijk0,ijk1,icells,ijcells,ii,jj,kk,kstart,kend;

  icells  = grid->icells;
  ijcells = grid->icells*grid->jcells;

  kstart = grid->kstart;
  kend   = grid->kend;

  ii = 1;
  jj = 1*icells;
  kk = 1*ijcells;

  if(sw == 0)
  {
    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*icells + (kstart-k-1)*ijcells;
          ijk1 = i + j*icells + (kstart+k  )*ijcells;
          data[ijk0] = -1.*data[ijk1];
        }

    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*icells + (kend+k  )*ijcells;
          ijk1 = i + j*icells + (kend-k-1)*ijcells;
          data[ijk0] = -1.*data[ijk1];
        }
  }
  else if(sw == 1)
  {
    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*icells + (kstart-k-1)*ijcells;
          ijk1 = i + j*icells + (kstart+k  )*ijcells;
          data[ijk0] = data[ijk1];
        }

    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*icells + (kend+k  )*ijcells;
          ijk1 = i + j*icells + (kend-k-1)*ijcells;
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
  int ijk0,ijk1,icells,jcells,ijcells,ii,jj,kk,istart,iend,jstart,jend,igc,jgc;

  icells  = grid->icells;
  jcells  = grid->jcells;
  ijcells = grid->icells*grid->jcells;

  istart = grid->istart;
  iend   = grid->iend;
  igc    = grid->igc;
  jstart = grid->jstart;
  jend   = grid->jend;
  jgc    = grid->jgc;

  ii = 1;
  jj = 1*icells;
  kk = 1*ijcells;

  // east west boundaries
  for(int k=0; k<grid->kcells; k++)
    for(int j=0; j<grid->jcells; j++)
      for(int i=0; i<grid->igc; i++)
      {
        ijk0 = i          + j*icells + k*ijcells;
        ijk1 = iend-igc+i + j*icells + k*ijcells;
        data[ijk0] = data[ijk1];
      }

  for(int k=0; k<grid->kcells; k++)
    for(int j=0; j<grid->jcells; j++)
      for(int i=0; i<grid->igc; i++)
      {
        ijk0 = i+iend   + j*icells + k*ijcells;
        ijk1 = i+istart + j*icells + k*ijcells;
        data[ijk0] = data[ijk1];
      }

  // north south boundaries
  for(int k=0; k<grid->kcells; k++)
    for(int j=0; j<grid->jgc; j++)
      for(int i=0; i<grid->icells; i++)
      {
        ijk0 = i + j           *icells + k*ijcells;
        ijk1 = i + (jend-jgc+j)*icells + k*ijcells;
        data[ijk0] = data[ijk1];
      }

  for(int k=0; k<grid->kcells; k++)
    for(int j=0; j<grid->jgc; j++)
      for(int i=0; i<grid->icells; i++)
      {
        ijk0 = i + (j+jend  )*icells + k*ijcells;
        ijk1 = i + (j+jstart)*icells + k*ijcells;
        data[ijk0] = data[ijk1];
      }

  // check loop
  // int k = 10;
  // int i = istart;

  // for(int j=0; j<jcells; j++)
  // {
  //   ijk0 = i + j*icells + k*ijcells;
  //   std::printf("j : %d, %f\n", j, data[ijk0]);
  // }

  // int j = jstart;
  // for(int i=0; i<icells; i++)
  // {
  //   ijk0 = i + j*icells + k*ijcells;
  //   std::printf("i : %d, %f\n", i, data[ijk0]);
  // }

  return 0;
}

int cfield3d::dump()
{
  FILE *pFile;
  pFile = fopen(name->c_str(), "wb");
  fwrite(data, sizeof(double), grid->ncells, pFile);
  fclose(pFile);

  return 0;
}
