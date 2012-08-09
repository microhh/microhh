#include <cstdio>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "boundary.h"
#include "defines.h"

cboundary::cboundary(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  std::printf("Creating instance of object boundary\n");
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cboundary::~cboundary()
{
  std::printf("Destroying instance of object boundary\n");
}

int cboundary::readinifile(cinput *inputin)
{
  int n = 0;

  // obligatory parameters
  n += inputin->getItem(&iboundary, "physics", "iboundary");

  n += inputin->getItem(&bcbotmom, "fields", "bcbotmom");
  n += inputin->getItem(&bctopmom, "fields", "bctopmom");

  n += inputin->getItem(&bcbotscal, "fields", "bcbotscal");
  n += inputin->getItem(&bctopscal, "fields", "bctopscal");

  n += inputin->getItem(&sbot, "fields", "sbot");
  n += inputin->getItem(&stop, "fields", "stop");

    // if one argument fails, then crash
  if(n > 0)
    return 1;

  return 0;
}

int cboundary::exec()
{
  if(iboundary == 2)
  {
    // bottom boundary conditions
    setgcbot_2nd((*fields->u).data, bcbotmom, 0.);
    setgcbot_2nd((*fields->v).data, bcbotmom, 0.);
    // setgcbot((*fields->w).data);
    setgcbot_2nd((*fields->s).data, bcbotscal, sbot);

    // top boundary conditions
    setgctop_2nd((*fields->u).data, bctopmom, 0.);
    setgctop_2nd((*fields->v).data, bctopmom, 0.);
    // setgcbot((*fields->w).data);
    setgctop_2nd((*fields->s).data, bctopscal, stop);
  }
  else if(iboundary == 4)
  {
    // bottom boundary conditions
    setgcbot_4th((*fields->u).data, grid->z, bcbotmom, 0.);
    setgcbot_4th((*fields->v).data, grid->z, bcbotmom, 0.);
    // setgcbot((*fields->w).data);
    setgcbot_4th((*fields->s).data, grid->z, bcbotscal, sbot);

    // top boundary conditions
    setgctop_4th((*fields->u).data, grid->z, bctopmom, 0.);
    setgctop_4th((*fields->v).data, grid->z, bctopmom, 0.);
    // setgcbot((*fields->w).data);
    setgctop_4th((*fields->s).data, grid->z, bctopscal, stop);
  }
 
  // cyclic boundary conditions
  grid->boundary_cyclic((*fields->u).data);
  grid->boundary_cyclic((*fields->v).data);
  grid->boundary_cyclic((*fields->w).data);
  grid->boundary_cyclic((*fields->s).data);
  mpi->waitall();

  return 0;
}

int cboundary::setgcbot_2nd(double * restrict a, int sw, double abot)
{ 
  int ijk,jj,kk,kstart;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  kstart = grid->kstart;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ijk = i + j*jj + kstart*kk;
        a[ijk-kk] = 2.*abot - a[ijk];
      }
  }
  else if(sw == 1)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ijk = i + j*jj + kstart*kk;
        a[ijk-kk] = a[ijk];
      }
  }

  return 0;
}

int cboundary::setgctop_2nd(double * restrict a, int sw, double atop)
{ 
  int ijk,jj,kk,kend;

  kend = grid->kend;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        // add the bcvalues later
        ijk = i + j*jj + (kend-1)*kk;
        a[ijk+kk] = 2.*atop - a[ijk];
      }
  }
  else if(sw == 1)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        // add the bcvalues later
        ijk = i + j*jj + (kend-1)*kk;
        a[ijk+kk] = a[ijk];
      }
  }

  return 0;
}

int cboundary::setgcbot_4th(double * restrict a, double * restrict z, int sw, double abot)
{ 
  int ijk,jj,kk1,kk2,kk3,kstart;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

  kstart = grid->kstart;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ijk = i + j*jj + kstart*kk1;
        a[ijk-kk1] = (16./5.)*abot - 3.*a[ijk] + a[ijk+kk1] - (1./5.)*a[ijk+kk2];
      }
  }
  else if(sw == 1)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ijk = i + j*jj + kstart*kk1;
        a[ijk-kk1] = ( - grad4xbiasbot(z[kstart-1], z[kstart], z[kstart+1], z[kstart+2])*abot 
                     + 21.*a[ijk] + 3.*a[ijk+kk1] - a[ijk+kk2] ) / 23.;
      }
  }

  return 0;
}

int cboundary::setgctop_4th(double * restrict a, double * restrict z, int sw, double atop)
{ 
  int ijk,jj,kend,kk1,kk2,kk3;

  kend = grid->kend;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

  if(sw == 0)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ijk = i + j*jj + (kend-1)*kk1;
        a[ijk+kk1] = (16./5.)*atop - 3.*a[ijk] + a[ijk-kk1] - (1./5.)*a[ijk-kk2];
      }
  }
  else if(sw == 1)
  {
    for(int j=0; j<grid->jcells; j++)
#pragma ivdep
      for(int i=0; i<grid->icells; i++)
      {
        ijk = i + j*jj + (kend-1)*kk1;
        a[ijk+kk1] = ( grad4xbiastop(z[kend-3], z[kend-2], z[kend-1], z[kend]) * atop 
                     + 21.*a[ijk] + 3.*a[ijk-kk1] - a[ijk-kk2] ) / 23.;
      }
  }

  return 0;
}

inline double cboundary::grad4xbiasbot(const double a, const double b, const double c, const double d)
{
  return (-23.*a + 21.*b + 3.*c - d);
}

inline double cboundary::grad4xbiastop(const double a, const double b, const double c, const double d)
{
  return ( 23.*d - 21.*c - 3.*b + a);
}

