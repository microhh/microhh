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
  n += inputin->getItem(&bcbotmom, "fields", "bcbotmom");
  n += inputin->getItem(&bctopmom, "fields", "bctopmom");

  n += inputin->getItem(&bcbotscal, "fields", "bcbotscal");
  n += inputin->getItem(&bctopscal, "fields", "bctopscal");

    // if one argument fails, then crash
  if(n > 0)
    return 1;

  return 0;
}

int cboundary::exec()
{
  /*
  // bottom boundary conditions
  setgcbot_2nd((*fields->u).data, bcbotmom);
  setgcbot_2nd((*fields->v).data, bcbotmom);
  // setgcbot((*fields->w).data);
  setgcbot_2nd((*fields->s).data, bcbotscal);

  // top boundary conditions
  setgctop_2nd((*fields->u).data, bctopmom);
  setgctop_2nd((*fields->v).data, bctopmom);
  // setgcbot((*fields->w).data);
  setgctop_2nd((*fields->s).data, bctopscal);
  */
 
  // bottom boundary conditions
  setgcbot_4th((*fields->u).data, bcbotmom);
  setgcbot_4th((*fields->v).data, bcbotmom);
  // setgcbot((*fields->w).data);
  setgcbot_4th((*fields->s).data, bcbotscal);

  // top boundary conditions
  setgctop_4th((*fields->u).data, bctopmom);
  setgctop_4th((*fields->v).data, bctopmom);
  // setgcbot((*fields->w).data);
  setgctop_4th((*fields->s).data, bctopscal);
 
  // cyclic boundary conditions
  grid->boundary_cyclic((*fields->u).data);
  grid->boundary_cyclic((*fields->v).data);
  grid->boundary_cyclic((*fields->w).data);
  grid->boundary_cyclic((*fields->s).data);
  mpi->waitall();

  return 0;
}

int cboundary::setgcbot_2nd(double * restrict a, int sw)
{ 
  int ijk0,ijk1,jj,kk,kstart;

  kstart = grid->kstart;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  if(sw == 0)
  {
    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
#pragma ivdep
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*jj + (kstart-k-1)*kk;
          ijk1 = i + j*jj + (kstart+k  )*kk;
          a[ijk0] = -1.*a[ijk1];
        }
  }
  else if(sw == 1)
  {
    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
#pragma ivdep
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*jj + (kstart-k-1)*kk;
          ijk1 = i + j*jj + (kstart+k  )*kk;
          a[ijk0] = a[ijk1];
        }
  }

  return 0;
}

int cboundary::setgctop_2nd(double * restrict a, int sw)
{ 
  int ijk0,ijk1,jj,kk,kend;

  kend   = grid->kend;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  if(sw == 0)
  {
    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
#pragma ivdep
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*jj + (kend+k  )*kk;
          ijk1 = i + j*jj + (kend-k-1)*kk;
          a[ijk0] = -1.*a[ijk1];
        }
  }
  else if(sw == 1)
  {
    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
#pragma ivdep
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*jj + (kend+k  )*kk;
          ijk1 = i + j*jj + (kend-k-1)*kk;
          a[ijk0] = a[ijk1];
        }
  }

  return 0;
}

int cboundary::setgcbot_4th(double * restrict a, int sw)
{ 
  int ijk0,ijk1,jj,kk,kstart;

  kstart = grid->kstart;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  if(sw == 0)
  {
    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
#pragma ivdep
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*jj + (kstart-k-1)*kk;
          ijk1 = i + j*jj + (kstart+k  )*kk;
          a[ijk0] = -1.*a[ijk1];
        }
  }
  else if(sw == 1)
  {
    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
#pragma ivdep
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*jj + (kstart-k-1)*kk;
          ijk1 = i + j*jj + (kstart+k  )*kk;
          a[ijk0] = a[ijk1];
        }
  }

  return 0;
}

int cboundary::setgctop_4th(double * restrict a, int sw)
{ 
  int ijk0,ijk1,jj,kk,kend;

  kend   = grid->kend;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  if(sw == 0)
  {
    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
#pragma ivdep
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*jj + (kend+k  )*kk;
          ijk1 = i + j*jj + (kend-k-1)*kk;
          a[ijk0] = -1.*a[ijk1];
        }
  }
  else if(sw == 1)
  {
    for(int k=0; k<grid->kgc; k++)
      for(int j=0; j<grid->jcells; j++)
#pragma ivdep
        for(int i=0; i<grid->icells; i++)
        {
          ijk0 = i + j*jj + (kend+k  )*kk;
          ijk1 = i + j*jj + (kend-k-1)*kk;
          a[ijk0] = a[ijk1];
        }
  }

  return 0;
}

