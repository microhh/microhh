#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fftw3.h>
#include "grid.h"
#include "fields.h"
#include "pres.h"
#include "defines.h"

cpres::cpres(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cpres::~cpres()
{
}

int cpres::readinifile(cinput *inputin)
{
  return 0;
}

int cpres::init()
{
  return 0;
}

int cpres::setvalues()
{
  return 0;
}

int cpres::exec(double dt)
{
  return 0;
}

double cpres::check()
{
  double divmax = 0.;

  return divmax;
}
