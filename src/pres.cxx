#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "pres.h"

cpres::cpres(cgrid *gridin, cfields *fieldsin)
{
  std::printf("Creating instance of object pres\n");
  grid   = gridin;
  fields = fieldsin;
}

cpres::~cpres()
{
  std::printf("Destroying instance of object pres\n");
}

int cpres::exec()
{
  return 0;
}

