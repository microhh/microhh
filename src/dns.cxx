#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "dns.h"

cdns::cdns(cgrid *gridin, cfields *fieldsin)
{
  grid   = gridin;
  fields = fieldsin;
  std::printf("Creating instance of object dns\n");
}

cdns::~cdns()
{
  std::printf("Destroying instance of object dns\n");
}
