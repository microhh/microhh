#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "grid.h"
#include "fields.h"
#include "defines.h"

cfields::cfields(cgrid *gridin, cmpi *mpiin)
{
  std::printf("Creating instance of object fields\n");
  grid = gridin;
  mpi  = mpiin;

  allocated = false;
}

cfields::~cfields()
{
  if(allocated)
  {
    delete u;
    delete v;
    delete w;
    delete p;

    delete ut;
    delete vt;
    delete wt;

    delete s;
    delete st;
  }

  std::printf("Destroying instance of object fields\n");
}

int cfields::readinifile(cinput *inputin)
{
  // input parameters
  int n = 0;

  // obligatory parameters
  n += inputin->getItem(&visc , "fields", "visc" );
  n += inputin->getItem(&viscs, "fields", "viscs");

  // optional parameters
  n += inputin->getItem(&rndamp, "fields", "rndamp", 1.e-3);

  if(n > 0)
    return 1;

  return 0;
}

int cfields::initfields()
{
  std::printf("Initializing fields\n");

  // set pointers to correct location
  u  = new cfield3d(grid, mpi, "u" );
  v  = new cfield3d(grid, mpi, "v" );
  w  = new cfield3d(grid, mpi, "w" );
  p  = new cfield3d(grid, mpi, "p" );

  ut = new cfield3d(grid, mpi, "ut");
  vt = new cfield3d(grid, mpi, "vt");
  wt = new cfield3d(grid, mpi, "wt");

  s  = new cfield3d(grid, mpi, "s" );

  st = new cfield3d(grid, mpi, "st");

  u->init();
  v->init();
  w->init();
  p->init();

  ut->init();
  vt->init();
  wt->init();

  s->init();

  st->init();

  allocated = true;

  return 0;
}

int cfields::createfields()
{
  std::printf("Creating fields\n");
  
  /*// set Taylor-Green vortex as default setup
  const double pi = std::acos((double)-1.);

  visc = 1. / (8.*pi*pi*100.);

  int ijk,ii,jj,kk;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        u->data[ijk] =  std::sin(2.*pi*(i-grid->istart)/grid->itot) *std::cos(2.*pi*grid->z [k]);
        w->data[ijk] = -std::cos(2.*pi*(0.5*grid->dx+(i-grid->istart)*grid->dx))*std::sin(2.*pi*grid->zh[k]);
      }
  // end Taylor-Green vortex setup */

  // set Moser180 as a default setup
  // put initial perturbation in u, v and w, set mpiid as random seed to avoid having the same field at all procs
  std::srand(mpi->mpiid);

  for(int n=0; n<grid->ncells; n++)
    u->data[n] = rndamp * (double)(std::rand() % 10000 - 5000) / 10000.;
  for(int n=0; n<grid->ncells; n++)
    v->data[n] = rndamp * (double)(std::rand() % 10000 - 5000) / 10000.;
  for(int n=0; n<grid->ncells; n++)
    w->data[n] = rndamp * (double)(std::rand() % 10000 - 5000) / 10000.;


  // add a double vortex to the initial conditions
  const double pi = std::acos((double)-1.);

  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        v->data[ijk] =  0.002*std::sin(2.*pi*(grid->y[j])/grid->ysize)*std::cos(pi*grid->z[k]/grid->zsize);
        w->data[ijk] = -0.002*std::cos(2.*pi*(grid->y[j])/grid->ysize)*std::sin(pi*grid->z[k]/grid->zsize);
      }

  // set the mean profile
  double dpdxls = -8.e-7;
  int k;

  for(int n=0; n<grid->ncells; n++)
  {
    k           = n / (grid->icells*grid->jcells);
    u->data[n] += 1./(2.*visc)*dpdxls*(grid->z[k]*grid->z[k] - grid->zsize*grid->z[k]);
    s->data[n] += (double)(k+1-grid->kgc) * grid->zsize / (double)(grid->ktot+1);
  }
  // end Moser180 setup 

  // set w equal to zero at the boundaries
  int nbot = grid->kstart*grid->icells*grid->jcells;
  int ntop = grid->kend  *grid->icells*grid->jcells;
  for(int n=0; n<grid->icells*grid->jcells; n++)
  {
    w->data[nbot + n] = 0.;
    w->data[ntop + n] = 0.;
  }

  return 0;
}

int cfields::load(int n)
{
  // check them all before returning error
  int nerror = 0;
  nerror += u->load(n);
  nerror += v->load(n);
  nerror += w->load(n);
  nerror += s->load(n);

  if(nerror > 0)
    return 1;

  return 0;
}

int cfields::save(int n)
{
  u->save(n);
  v->save(n);
  w->save(n);
  p->save(n);
  s->save(n);

  return 0;
}

int cfields::boundary()
{
  u->boundary_bottop(0);
  v->boundary_bottop(0);
  // w->boundary_bottop(0);
  s->boundary_bottop(1);

  // u->boundary_cyclic();
  // v->boundary_cyclic();
  // w->boundary_cyclic();
  // s->boundary_cyclic();
  
  grid->boundary_cyclic(u->data);
  grid->boundary_cyclic(v->data);
  grid->boundary_cyclic(w->data);
  grid->boundary_cyclic(s->data);

  // for(int k=grid->kstart-grid->kgc; k<grid->kend+grid->kgc; k++)
  //   std::printf("%4d %9.6f %9.6f %9.6f %9.6f %9.6f\n", k, grid->z[k], grid->zh[k], u->data[k*grid->icells*grid->jcells], v->data[k*grid->icells*grid->jcells], w->data[k*grid->icells*grid->jcells]);

  return 0;
}

double cfields::check(int n)
{
  double checkval = 0.;

  if(n == 0)
    checkval = calcmom (u->data, v->data, w->data, grid->dz);
  else if(n == 1)
    checkval = calctke (u->data, v->data, w->data, grid->dz);
  else if(n == 2)
    checkval = calcmass(s->data, grid->dz);

  return checkval;
}

double cfields::calcmass(double * restrict s, double * restrict dz)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  double mass;
  mass = 0;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        mass += s[ijk]*dz[k];
      }

  grid->getsum(&mass);

  mass /= (grid->itot*grid->jtot*grid->zsize);

  return mass;
}


double cfields::calcmom(double * restrict u, double * restrict v, double * restrict w, double * restrict dz)
{
  int ijk,icells,ijcells,ii,jj,kk;

  icells  = grid->icells;
  ijcells = grid->icells*grid->jcells;

  ii = 1;
  jj = 1*icells;
  kk = 1*ijcells;

  double momentum;
  momentum = 0;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*icells + k*ijcells;
        momentum += (interp2(u[ijk], u[ijk+ii]) + interp2(v[ijk], v[ijk+jj]) + interp2(w[ijk], w[ijk+kk]))*dz[k];
      }

  grid->getsum(&momentum);

  momentum /= (grid->itot*grid->jtot*grid->zsize);

  return momentum;
}

double cfields::calctke(double * restrict u, double * restrict v, double * restrict w, double * restrict dz)
{
  int ijk,icells,ijcells,ii,jj,kk;

  icells  = grid->icells;
  ijcells = grid->icells*grid->jcells;

  ii = 1;
  jj = 1*icells;
  kk = 1*ijcells;

  double tke;
  tke = 0;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*icells + k*ijcells;
        tke += ( interp2(u[ijk]*u[ijk], u[ijk+ii]*u[ijk+ii]) 
               + interp2(v[ijk]*v[ijk], v[ijk+jj]*v[ijk+jj]) 
               + interp2(w[ijk]*w[ijk], w[ijk+kk]*w[ijk+kk]))*dz[k];
      }

  grid->getsum(&tke);

  tke /= (grid->itot*grid->jtot*grid->zsize);
  tke *= 0.5;

  return tke;
}

inline double cfields::interp2(const double a, const double b)
{
  return 0.5*(a + b);
}

