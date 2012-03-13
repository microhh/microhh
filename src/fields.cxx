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
  double rndamp =  1.e-3;
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

int cfields::boundary()
{
  u->boundary_bottop(0);
  v->boundary_bottop(0);
  w->boundary_bottop(0);

  u->boundary_cyclic();
  v->boundary_cyclic();
  w->boundary_cyclic();

  return 0;
}

int cfields::resettend()
{
  for(int n=0; n<grid->ncells*3; n++)
    flowt[n] = 0.;
  return 0;
}

int cfields::check()
{
  double mom;

  mom = momentum(u->data, v->data, w->data, grid->dz);

  std::printf("Total momentum = %24.14f\n", mom);

  return 0;
}

double cfields::momentum(double * __restrict__ u, double * __restrict__ v, double * __restrict__ w, double * __restrict__ dz)
{
  int    ijk,icells,ijcells,ii,jj,kk;

  icells  = grid->icells;
  ijcells = grid->icells*grid->jcells;

  ii = 1;
  jj = 1*icells;
  kk = 1*ijcells;

  double momentum;
  momentum = 0;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*icells + k*ijcells;
        momentum += (interp2(u[ijk], u[ijk+ii]) + interp2(v[ijk], v[ijk+jj]) + interp2(w[ijk], w[ijk+kk]))*dz[k];
      }

  momentum /= (grid->imax*grid->jmax*grid->zsize);

  return momentum;
}

inline double cfields::interp2(const double a, const double b)
{
  return 0.5*(a + b);
}

