#include <cstdio>
#include <cmath>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "buffer.h"
#include "defines.h"

cbuffer::cbuffer(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cbuffer::~cbuffer()
{
  if(allocated)
    delete[] bufferprofs;
}

int cbuffer::readinifile(cinput *inputin)
{
  int n = 0;

  // optional parameters
  n += inputin->getItem(&ibuffer,      "fields", "ibuffer"     , 0 );
  n += inputin->getItem(&bufferkstart, "fields", "bufferkstart", 0 );
  n += inputin->getItem(&buffersigma,  "fields", "buffersigma" , 2.);
  n += inputin->getItem(&bufferbeta,   "fields", "bufferbeta"  , 2.);

    // if one argument fails, then crash
  if(n > 0)
    return 1;

  return 0;
}

int cbuffer::init()
{
  // allocate the buffer array 
  bufferkcells = grid->kmax - bufferkstart;

  // CvH fix this later with flexible number of scalars
  bufferprofs = new double[4*bufferkcells];

  allocated = true;

  // add the ghost cells to the starting point
  bufferkstart += grid->kstart;

  return 0;
}

int cbuffer::setbuffers()
{
  if(ibuffer)
  {
    // set the buffers according to the initial profiles
    setbuffer((*fields->u).data, &bufferprofs[0*bufferkcells]);
    setbuffer((*fields->v).data, &bufferprofs[1*bufferkcells]);
    setbuffer((*fields->w).data, &bufferprofs[2*bufferkcells]);
    setbuffer((*fields->s).data, &bufferprofs[3*bufferkcells]);
  }

  return 0;
}

int cbuffer::exec()
{
  if(ibuffer == 1)
  {
    // calculate the buffer tendencies
    buffer((*fields->ut).data, (*fields->u).data, &bufferprofs[0*bufferkcells], grid->z );
    buffer((*fields->vt).data, (*fields->v).data, &bufferprofs[1*bufferkcells], grid->z );
    buffer((*fields->wt).data, (*fields->w).data, &bufferprofs[2*bufferkcells], grid->zh);
    buffer((*fields->st).data, (*fields->s).data, &bufferprofs[3*bufferkcells], grid->z );
  }

  return 0;
}

int cbuffer::buffer(double * restrict at, double * restrict a, double * restrict abuf, double * restrict z)
{ 
  int ijk,jj,kk;
  int kloopstart;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  kloopstart = bufferkstart+1;

  double sigma;
  double zsizebuf;

  zsizebuf = grid->zsize - z[bufferkstart];

  for(int k=kloopstart; k<grid->kend; k++)
  {
    sigma = buffersigma*std::pow((z[k]-z[bufferkstart])/zsizebuf, bufferbeta);
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        at[ijk] -= sigma*(a[ijk]-abuf[k-kloopstart]);
      }
  }

  return 0;
}

int cbuffer::setbuffer(double * restrict a, double * restrict abuf)
{
  int ijk,jj,kk;
  int kloopstart;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  kloopstart = bufferkstart+1;

  for(int k=kloopstart; k<grid->kend; k++)
  {
    abuf[k-kloopstart] = 0.;
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        abuf[k-kloopstart] += a[ijk];
      }

    abuf[k-kloopstart] /= grid->imax*grid->jmax;
  }

  grid->getprof(abuf, bufferkcells);

  return 0;
}

