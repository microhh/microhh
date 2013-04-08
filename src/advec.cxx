#include <cstdio>
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "advec.h"
#include "defines.h"

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
  int n = 0;
  n += inputin->getItem(&cflmax      , "advec", "cflmax"      , "", 1.   );

  return (n>0);
}

unsigned long cadvec::gettimelim(unsigned long idt, double dt)
{
  unsigned long idtlim = (unsigned long) dbig;

  return idtlim;
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

