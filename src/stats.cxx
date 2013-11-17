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
  return 0;
}

int cstats::init(double ifactor)
{
  return 0;
}

int cstats::create(int n)
{
  return 0;
}

unsigned long cstats::gettimelim(unsigned long itime)
{
  return ulhuge;
}

int cstats::exec(int iteration, double time, unsigned long itime)
{
  return 0;
}

