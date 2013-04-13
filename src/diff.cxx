#include <cstdio>
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "diff.h"
#include "defines.h"

cdiff::cdiff(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cdiff::~cdiff()
{
}


int cdiff::readinifile(cinput *inputin)
{
  int n = 0;

  n += inputin->getItem(&dnmax      , "diff", "dnmax"      , "", 0.5   );

  return (n>0);
}

unsigned long cdiff::gettimelim(unsigned long idtlim, double dt)
{
  idtlim = (unsigned long) dbig;

  return idtlim;
}

int cdiff::setvalues()
{
  return 0;
}

double cdiff::getdn(double dt)
{
  double dn;

  dn = dsmall;

  return dn;
}

int cdiff::exec()
{
  return 0;
}

