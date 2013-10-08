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

  allocated = false;
}

cbuffer::~cbuffer()
{
  if(allocated)
    for (std::map<std::string,double*>::iterator it = bufferprofs.begin(); it!=bufferprofs.end(); it++)
      delete[] it->second;
}

int cbuffer::readinifile(cinput *inputin)
{
  int n = 0;

  // optional parameters
  n += inputin->getItem(&swbuffer, "buffer", "swbuffer", "", "0");

  if(swbuffer == "1")
  {
    n += inputin->getItem(&bufferz    , "buffer", "bufferz"     , "");
    n += inputin->getItem(&buffersigma, "buffer", "buffersigma" , "", 2.);
    n += inputin->getItem(&bufferbeta , "buffer", "bufferbeta"  , "", 2.);
  }

  // if one argument fails, then crash
  if(n > 0)
    return 1;

  return 0;
}

int cbuffer::init()
{
  if(swbuffer == "1")
  {
    // allocate the buffer arrays
    for(fieldmap::const_iterator it=fields->mp.begin(); it!=fields->mp.end(); ++it)
      bufferprofs[it->first] = new double[grid->kcells];

    for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
      bufferprofs[it->first] = new double[grid->kcells];

    allocated = true;
  }

  return 0;
}

int cbuffer::create(cinput *inputin)
{
  int nerror = 0;

  if(swbuffer == "1")
  {
    // set the buffers according to the initial profiles of the variables
    nerror += inputin->getProf(&bufferprofs["u"][grid->kstart], "u", grid->kmax);
    nerror += inputin->getProf(&bufferprofs["v"][grid->kstart], "v", grid->kmax);
 
    for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
      nerror += inputin->getProf(&bufferprofs[it->first][grid->kstart], it->first, grid->kmax);

    // find the starting points
    bufferkstart  = grid->kstart;
    bufferkstarth = grid->kstart;

    for(int k=grid->kstart; k<grid->kend; ++k)
    {
      // check if the cell center is in the buffer zone
      if(grid->z[k] < bufferz)
        ++bufferkstart;
      // check if the cell face is in the buffer zone
      if(grid->zh[k] < bufferz)
        ++bufferkstarth;

    }

    // check whether the lowest of the two levels is contained in the buffer layer
    if(bufferkstarth == grid->kend)
    {
      ++nerror;
      if(mpi->mpiid == 0) std::printf("ERROR buffer is too close to the model top\n");
    }
  }
  return nerror;
}

int cbuffer::exec()
{
  if(swbuffer == "1")
  {
    // calculate the buffer tendencies
    buffer((*fields->mt["u"]).data, (*fields->mp["u"]).data, bufferprofs["u"], grid->z );
    buffer((*fields->mt["v"]).data, (*fields->mp["v"]).data, bufferprofs["v"], grid->z );
    buffer((*fields->mt["w"]).data, (*fields->mp["w"]).data, bufferprofs["w"], grid->zh);
 
    for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
      buffer(fields->st[it->first]->data, it->second->data, bufferprofs[it->first], grid->z);
  }

  return 0;
}

int cbuffer::buffer(double * const restrict at, const double * const restrict a, 
                    const double * const restrict abuf, const double * const restrict z)
{ 
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  double sigma;
  double zsizebuf;

  zsizebuf = grid->zsize - bufferz;

  for(int k=bufferkstart; k<grid->kend; k++)
  {
    sigma = buffersigma*std::pow((z[k]-bufferz)/zsizebuf, bufferbeta);
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        at[ijk] -= sigma*(a[ijk]-abuf[k]);
      }
  }

  return 0;
}

