#include <cstdio>
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "advec.h"
#include "defines.h"

#include <map>

cadvec::cadvec(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cadvec::~cadvec()
{
}

int cadvec::readinifile(cinput *inputin)
{
  return 0;
}

double cadvec::getcfl(double dt)
{
  double cfl;

  cfl = dsmall;
  
  return cfl;
}

int cadvec::exec()
{
  return 0;
}
