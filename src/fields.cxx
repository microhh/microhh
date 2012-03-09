#include <cstdlib>
#include <cstdio>
#include "grid.h"
#include "fields.h"

cfields::cfields(cgrid *gridin)
{
  grid = gridin;
  std::printf("Creating fields\n");
  // allocate memory for 3d arrays
  std::printf("Allocating %d bytes of memory for flow\n", grid->ncells*3*sizeof(double));
  flow  = new double[grid->ncells*3];
  std::printf("Allocating %d bytes of memory for flowt\n", grid->ncells*3*sizeof(double));
  flowt = new double[grid->ncells*3];

  // set pointers to correct location
  cfield u (grid, &flow[grid->ncells*0]);
  cfield v (grid, &flow[grid->ncells*1]);
  cfield w (grid, &flow[grid->ncells*2]);

  cfield ut(grid, &flowt[grid->ncells*0]);
  cfield vt(grid, &flowt[grid->ncells*1]);
  cfield wt(grid, &flowt[grid->ncells*2]);

  // set all values to 0
  for(int n=0; n<grid->ncells*3; n++)
    flow[n] = 0.;

  // set all tendencies to 0
  resettend();

  // set Moser180 as a default setup
  double dpdxls = -1.5e-6;
  double visc   =  1.0e-5;
  double rndamp =  1.e-5;
  int k;

  // put initial perturbation in u, v and w
  std::srand(0);
  for(int n=0; n<grid->ncells; n++)
    flow[n] = rndamp * (double)(std::rand() % 10000) / 10000.;

  for(int n=0; n<grid->ncells; n++)
  {
    k          = n / (grid->icells*grid->jcells);
    u.data[n] += 1./(2.*visc)*dpdxls*(grid->z[k]*grid->z[k] - grid->zsize*grid->z[k]);
  }

  u.boundary_bottop(0);
  v.boundary_bottop(0);
  //w.boundary_bottop(0);
  // end Moser180 setup

  for(int k=grid->kstart-grid->kgc; k<grid->kend+grid->kgc; k++)
    std::printf("%4d %9.6f %9.6f %9.6f %9.6f %9.6f \n", k-grid->kstart+1, grid->z[k], grid->zh[k], grid->dz[k], grid->dzh[k], u.data[k*grid->icells*grid->jcells]);
}

int cfields::resettend()
{
  for(int n=0; n<grid->ncells*3; n++)
    flowt[n] = 0.;
  return 0;
}

cfields::~cfields()
{
  std::printf("Deleting fields\n");
  delete[] flow;
  delete[] flowt;
}

cfield::cfield(cgrid *gridin, double *dataref)
{
  grid = gridin;
  data = dataref;
}

int cfield::index(int i, int j, int k)
{
  int n = i + j*grid->icells + k*grid->icells*grid->jcells;
  return n;
}


int cfield::boundary_bottop(int sw)
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
}

