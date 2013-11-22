#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "thermo_dry.h"
#include "defines.h"

cthermo_dry::cthermo_dry(cgrid *gridin, cfields *fieldsin, cmpi *mpiin) : cthermo(gridin, fieldsin, mpiin)
{
}

cthermo_dry::~cthermo_dry()
{
}

int cthermo_dry::readinifile(cinput *inputin)
{
  int nerror = 0;
  nerror += inputin->getItem(&gravitybeta, "thermo", "gravitybeta", "");

  nerror += fields->initpfld("s");
  nerror += inputin->getItem(&fields->sp["s"]->visc, "fields", "svisc", "s");

  return nerror;
}

int cthermo_dry::exec()
{
  // extend later for gravity vector not normal to surface
  if(grid->swspatialorder== "2")
    buoyancy_2nd(fields->wt->data, fields->s["s"]->data);
  else if(grid->swspatialorder == "4")
    buoyancy_4th(fields->wt->data, fields->s["s"]->data);

  return 0;
}

int cthermo_dry::buoyancy_2nd(double * restrict wt, double * restrict s)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  // CvH check the usage of the gravity term here, in case of scaled DNS we use one. But thermal expansion coeff??
  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        wt[ijk] += gravitybeta * interp2(s[ijk-kk], s[ijk]);
      }

  return 0;
}

int cthermo_dry::buoyancy_4th(double * restrict wt, double * restrict s)
{
  int ijk,jj;
  int kk1,kk2;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk1;
        wt[ijk] += gravitybeta * interp4(s[ijk-kk2], s[ijk-kk1], s[ijk], s[ijk+kk1]);
      }

  return 0;
}

inline double cthermo_dry::interp2(const double a, const double b)
{
  return 0.5*(a + b);
}

inline double cthermo_dry::interp4(const double a, const double b, const double c, const double d)
{
  return (-a + 9.*b + 9.*c - d) / 16.;
}

