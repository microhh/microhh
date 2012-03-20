#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "force.h"

cforce::cforce(cgrid *gridin, cfields *fieldsin)
{
  std::printf("Creating instance of object force\n");
  grid   = gridin;
  fields = fieldsin;
}

cforce::~cforce()
{
  std::printf("Destroying instance of object force\n");
}

int cforce::exec()
{
  return 0;
}
