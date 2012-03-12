#include <cstdio>
#include "grid.h"
#include "field3d.h"

cfield3d::cfield3d(cgrid *gridin, double *dataref)
{
  std::printf("Creating instance of object field3d\n");
  grid = gridin;
  data = dataref;
}

cfield3d::~cfield3d()
{
  std::printf("Destroying instance of object field3d\n");
}

int cfield3d::index(int i, int j, int k)
{
  return i + j*grid->icells + k*grid->icells*grid->jcells;
}

int cfield3d::boundary_bottop(int sw)
{ 
  if(sw == 0)
  {
    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
        for(int i=0; i<grid->icells; i++)
          data[index(i,j,grid->kstart-k-1)] = -1.*data[index(i,j,grid->kstart+k)];

    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
        for(int i=0; i<grid->icells; i++)
          data[index(i,j,grid->kend+k)] = -1.*data[index(i,j,grid->kend-k-1)];
  }
  else if(sw == 1)
  {
    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
        for(int i=0; i<grid->icells; i++)
          data[index(i,j,grid->kstart-k-1)] = data[index(i,j,grid->kstart+k)];

    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
        for(int i=0; i<grid->icells; i++)
          data[index(i,j,grid->kend+k)] = data[index(i,j,grid->kend-k-1)];
  }

  return 0;
}

int cfield3d::boundary_cyclic()
{ 
  for(int k=0; k<grid->kgc; k++)
    for(int j=0; j<grid->jcells; j++)
      for(int i=0; i<grid->icells; i++)
        data[index(i,j,grid->kstart-k-1)] = -1.*data[index(i,j,grid->kstart+k)];

  for(int k=0; k<grid->kgc; k++)
    for(int j=0; j<grid->jcells; j++)
      for(int i=0; i<grid->icells; i++)
        data[index(i,j,grid->kend+k)] = -1.*data[index(i,j,grid->kend-k-1)];

  return 0;
}
