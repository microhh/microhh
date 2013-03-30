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
  n += inputin->getItem(&swstats, "stats", "swstats", "", "0");

  if(n > 0)
    return 1;

  return 0;
}

int cstats::init()
{
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
  if(swstats == "0")
    return 0;

  else if(swstats == "4")
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

int cstats::exec(int iteration, double time)
{
  if(swstats == "0")
    return 0;

  else if(swstats == "4")
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

