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

    bufferkstart  = grid->kstart;
    bufferkstarth = grid->kstart;

    // find the starting points
    for(int k=grid->kstart; k<grid->kend; ++k)
    {
      // check if the cell center is in the buffer zone
      if(grid->z[k] < bufferz)
        ++bufferkstart;
      // check if the cell face is in the buffer zone
      if(grid->zh[k] < bufferz)
        ++bufferkstarth;
    }
    std::printf("CvH: %d, %d\n", bufferkstart, bufferkstarth);
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
  }

  return 0;
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
        at[ijk] -= sigma*(a[ijk]-abuf[k]);
      }
  }

  return 0;
}

/*
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

int cbuffer::save()
{
  if(swbuffer != "1")
    return 0;

  char filename[256];
  std::sprintf(filename, "%s.%07d", "buffer", 0);

  if(mpi->mpiid == 0)
  {
    std::printf("Saving \"%s\"\n", filename);
    FILE *pFile;
    pFile = fopen(filename, "wb");

    if(pFile == NULL)
    {
      std::printf("ERROR \"%s\" cannot be written", filename);
      return 1;
    }

    for (std::map<std::string,double*>::iterator itBuffer = bufferprofs.begin(); itBuffer!=bufferprofs.end(); itBuffer++)
      fwrite(itBuffer->second, sizeof(double), bufferkcells, pFile);
    
    fclose(pFile);
  }

  return 0;
}

int cbuffer::load()
{
  int nerror = 0;

  if(swbuffer != "1")
    return 0;

  char filename[256];
  std::sprintf(filename, "%s.%07d", "buffer", 0);

  if(mpi->mpiid == 0)
  {
    std::printf("Loading \"%s\"\n", filename);

    FILE *pFile;
    pFile = fopen(filename, "rb");

    if(pFile == NULL)
    {
      std::printf("ERROR \"%s\" does not exist\n", filename);
      ++nerror;
    }
    else
    {
      for (std::map<std::string,double*>::iterator itBuffer = bufferprofs.begin(); itBuffer!=bufferprofs.end(); itBuffer++)
        fread(itBuffer->second, sizeof(double), bufferkcells, pFile);
    
      fclose(pFile);
    }
  }

  mpi->broadcast(&nerror, 1);
  if(nerror)
    return 1;

  // send the buffers to all processes
  for (std::map<std::string,double*>::iterator itBuffer = bufferprofs.begin(); itBuffer!=bufferprofs.end(); itBuffer++)
    mpi->broadcast(itBuffer->second, bufferkcells);

  return 0;
}
*/
