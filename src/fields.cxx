#include <cstdlib>
#include <cstdio>
#include "grid.h"
#include "fields.h"

cfields::cfields(cgrid *gridin)
{
  std::printf("Creating instance of object fields\n");
  grid = gridin;
}

cfields::~cfields()
{
  delete u;
  delete v;
  delete w;
  delete ut;
  delete vt;
  delete wt;
  delete[] flow;
  delete[] flowt;
  std::printf("Destroying instance of object fields\n");
}

int cfields::initfields()
{
  std::printf("Initializing fields\n");
  // allocate memory for 3d arrays
  std::printf("Allocating %d bytes of memory for flow\n", grid->ncells*3*(int)sizeof(double));
  flow  = new double[grid->ncells*3];
  std::printf("Allocating %d bytes of memory for flowt\n", grid->ncells*3*(int)sizeof(double));
  flowt = new double[grid->ncells*3];

  // set pointers to correct location
  u  = new cfield3d(grid, &flow[grid->ncells*0]);
  v  = new cfield3d(grid, &flow[grid->ncells*1]);
  w  = new cfield3d(grid, &flow[grid->ncells*2]);

  ut = new cfield3d(grid, &flowt[grid->ncells*0]);
  vt = new cfield3d(grid, &flowt[grid->ncells*1]);
  wt = new cfield3d(grid, &flowt[grid->ncells*2]);

  // set all values to 0
  for(int n=0; n<grid->ncells*3; n++)
    flow[n] = 0.;

  // set all tendencies to 0
  resettend();

  return 0;
}

int cfields::createfields()
{
  std::printf("Creating fields\n");
  // set Moser180 as a default setup
  visc = 1.0e-5;

  double dpdxls = -1.5e-6;
  double rndamp =  1.e-5;
  int k;

  // put initial perturbation in u, v and w
  std::srand(0);
  for(int n=0; n<grid->ncells*3; n++)
    flow[n] = rndamp * (double)(std::rand() % 10000) / 10000.;

  for(int n=0; n<grid->ncells; n++)
  {
    k           = n / (grid->icells*grid->jcells);
    u->data[n] += 1./(2.*visc)*dpdxls*(grid->z[k]*grid->z[k] - grid->zsize*grid->z[k]);
  }

  // u.boundary_bottop(0);
  // v.boundary_bottop(0);
  //w.boundary_bottop(0);
  // end Moser180 setup

  for(int k=grid->kstart-grid->kgc; k<grid->kend+grid->kgc; k++)
    std::printf("%4d %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n", k-grid->kstart+1, grid->z[k], grid->zh[k], grid->dz[k], grid->dzh[k], u->data[k*grid->icells*grid->jcells], v->data[k*grid->icells*grid->jcells]);

  return 0;
}

int cfields::boundary_bottop()
{
  u->boundary_bottop(0);
  v->boundary_bottop(0);
  w->boundary_bottop(0);

  return 0;
}

int cfields::resettend()
{
  for(int n=0; n<grid->ncells*3; n++)
    flowt[n] = 0.;
  return 0;
}

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

