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
  return 0;
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

