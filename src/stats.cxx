#include <cstdio>
#include <cmath>
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

  stats_dns = new cstats_dns(grid, fields, mpi);
  stats_les = new cstats_les(grid, fields, mpi);
}

cstats::~cstats()
{
  delete stats_dns;
  delete stats_les;
}

int cstats::readinifile(cinput *inputin)
{
  int n = 0;

  // optional, by default switch stats off
  n += inputin->getItem(&swstats  , "stats", "swstats"  , "", "0");
  n += inputin->getItem(&statstime, "stats", "statstime", "", 60.);

  if(n > 0)
    return 1;

  return 0;
}

int cstats::init(double ifactor)
{
  istatstime = (unsigned long)(ifactor * statstime);

  if(swstats == "0")
    return 0;

  else if(swstats == "4")
  {
    if(stats_dns->init())
      return 1;
  }

  else if(swstats == "22")
  {
    if(stats_les->init())
      return 1;
  }

  return 0;
}

int cstats::create(int n)
{
  // check if switched on
  if(swstats == "0")
    return 0;

  if(swstats == "4")
  {
    if(stats_dns->create(n))
      return 1;
  }
  else if(swstats == "22")
  {
    if(stats_les->create(n))
      return 1;
  }

  return 0;
}

unsigned long cstats::gettimelim(unsigned long itime)
{
  unsigned long idtlim = istatstime -  itime % istatstime;

  return idtlim;
}

int cstats::exec(int iteration, double time, int itime)
{
  if(swstats == "0")
    return 0;

  // check if time for execution
  if(itime % istatstime != 0)
    return 0;

  if(swstats == "4")
  {
    if(stats_dns->exec(iteration, time))
      return 1;
  }

  else if(swstats == "22")
  {
    if(stats_les->exec(iteration, time))
      return 1;
  }

  return 0;
}

