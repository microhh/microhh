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
  u  = &flow[grid->ncells*0];
  v  = &flow[grid->ncells*1];
  w  = &flow[grid->ncells*2];

  ut = &flowt[grid->ncells*0];
  vt = &flowt[grid->ncells*1];
  wt = &flowt[grid->ncells*2];

  // set all values to 0
  for(int n=0; n<grid->ncells*3; n++)
    flow[n] = 0.;

  // set all tendencies to 0
  resettend();

  // set Moser180 as a default setup
  double dpdxls = -1.5e-6;
  double visc   =  1.0e-5;
  int k;

  for(int n=0; n<grid->ncells; n++)
  {
    k    = n / (grid->icells*grid->jcells);
    u[n] = 1./(2.*visc)*dpdxls*(grid->z[k]*grid->z[k] - grid->zsize*grid->z[k]);
  }
  // end Moser180 setup

  for(int k=grid->kstart-grid->kgc; k<grid->kend+grid->kgc; k++)
    std::printf("%4d %9.6f %9.6f %9.6f %9.6f %9.6f \n", k-grid->kstart+1, grid->z[k], grid->zh[k], grid->dz[k], grid->dzh[k], u[k*grid->icells*grid->jcells]);
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

