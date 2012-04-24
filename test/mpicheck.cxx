#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "mpicheck.h"

cmpicheck::cmpicheck(cgrid *gridin, cfields *fieldsin)
{
  std::printf("Creating instance of object mpicheck\n");
  grid   = gridin;
  fields = fieldsin;
}

cmpicheck::~cmpicheck()
{
  std::printf("Destroying instance of object mpicheck\n");
}

int cmpicheck::create()
{
  return 0;
}

