#include <cstdio>
#include <netcdfcpp.h>
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "defines.h"

cstats::cstats(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cstats::~cstats()
{
  delete dataFile;
}

int cstats::readinifile(cinput *inputin)
{
  int n = 0;

  n += inputin->getItem(&istats, "postproc", "stats");

  if(n > 0)
    return 1;

  return 0;
}

int cstats::init()
{
  // create a NetCDF file for the statistics
  dataFile = new NcFile("test.nc", NcFile::Replace);
  if(!dataFile->is_valid())
  {
    std::printf("ERROR: can't write statistics file");
    return 1;
  }

  // create dimensions
  NcDim* zDim  = dataFile->add_dim("z" , grid->kmax);
  NcDim* zhDim = dataFile->add_dim("zh", grid->kmax+1);
  NcDim* tDim  = dataFile->add_dim("t");

  // create variables
  NcVar *u = dataFile->add_var("u", ncDouble, tDim, zDim );
  NcVar *v = dataFile->add_var("v", ncDouble, tDim, zDim );
  NcVar *w = dataFile->add_var("w", ncDouble, tDim, zhDim);
  NcVar *s = dataFile->add_var("s", ncDouble, tDim, zhDim);

  NcVar *u2 = dataFile->add_var("u2", ncDouble, tDim, zDim );
  NcVar *v2 = dataFile->add_var("v2", ncDouble, tDim, zDim );
  NcVar *w2 = dataFile->add_var("w2", ncDouble, tDim, zhDim);
  NcVar *s2 = dataFile->add_var("s2", ncDouble, tDim, zDim );

  dataFile->sync();

  return 0;
}

int cstats::exec(int iteration, double time)
{
  if(mpi->mpiid == 0) std::printf("Saving stats for iteration %d\n", iteration);

  // PROFILES
  // calculate means
  // calcmean(*(fields->u).data);

  // save means
  // saveprof(
  
  // calculate variance and higher order moments
  // save variances

  // calculate flux
  // save flux
  
  // SCALARS
  return 0;
}

