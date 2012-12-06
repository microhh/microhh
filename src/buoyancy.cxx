#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "buoyancy.h"
#include "defines.h"

cbuoyancy::cbuoyancy(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  // std::printf("Creating instance of object buoyancy\n");
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cbuoyancy::~cbuoyancy()
{
  // std::printf("Destroying instance of object buoyancy\n");
}

int cbuoyancy::readinifile(cinput *inputin)
{
  int n = 0;

  n += inputin->getItem(&ibuoyancy  , "physics", "ibuoyancy");
  n += inputin->getItem(&gravity    , "physics", "gravity", 9.81);
  n += inputin->getItem(&buoyancyref, "physics", "buoyancyref", 0.);

  if(n > 0)
    return 1;

  return 0;
}

int cbuoyancy::exec()
{
  if(ibuoyancy == 0)
    return 0;

  // extend later for gravity vector not normal to surface
  if(ibuoyancy == 2)
    buoyancy_2nd((*fields->wt).data, (*fields->s).data);
  else if(ibuoyancy == 4)
    buoyancy_4th((*fields->wt).data, (*fields->s).data);

  return 0;
}

int cbuoyancy::buoyancy_2nd(double * restrict wt, double * restrict s)
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
        wt[ijk] += gravity * interp2(s[ijk-kk], s[ijk]);
      }

  return 0;
}

int cbuoyancy::buoyancy_4th(double * restrict wt, double * restrict s)
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
        wt[ijk] += gravity * interp4(s[ijk-kk2], s[ijk-kk1], s[ijk], s[ijk+kk1]);
      }

  return 0;
}

inline double cbuoyancy::interp2(const double a, const double b)
{
  return 0.5*(a + b);
}

inline double cbuoyancy::interp4(const double a, const double b, const double c, const double d)
{
  return (-a + 9.*b + 9.*c - d) / 16.;
}

inline double cbuoyancy::interp4biasbot(const double a, const double b, const double c, const double d)
{
  return ((5./16.)*a + (15./16.)*b - (5./16.)*c + (1./16)*d);
}

inline double cbuoyancy::interp4biastop(const double a, const double b, const double c, const double d)
{
  return ((5./16.)*d + (15./16.)*c - (5./16.)*b + (1./16)*a);
}

