#include <cstdio>
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
}

int cstats::readinifile(cinput *inputin)
{
  int n = 0;

  n += inputin->getItem(&istats, "postproc", "stats");

  if(n > 0)
    return 1;

  return 0;
}

int cstats::exec(int iteration)
{
  if(mpi->mpiid == 0) std::printf("Saving stats for iteration %d\n", iteration);
  return 0;
}

