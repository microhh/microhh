#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "defines.h"

cthermo::cthermo(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cthermo::~cthermo()
{
}

int cthermo::readinifile(cinput *inputin)
{
  return 0;
}

int cthermo::create()
{
  return 0;
}


int cthermo::exec()
{
  return 0;
}
