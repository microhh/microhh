#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "timeint.h"

ctimeint::ctimeint(cgrid *gridin, cfields *fieldsin)
{
  std::printf("Creating instance of object timeint\n");
  grid   = gridin;
  fields = fieldsin;

  substep = 0;
}

ctimeint::~ctimeint()
{
  std::printf("Destroying instance of object timeint\n");
}

int ctimeint::exec(double dt)
{
  return rk3(fields->flow, fields->flowt, dt);
}

int ctimeint::rk3(double * __restrict__ a, double * __restrict__ at, double dt)
{
  const double cA [] = {0., -5./9., -153./128.};
  const double cB [] = {1./3., 15./16., 8./15.};
  const double cdt[] = {0., 1./3., 3./4.};
  int n;

  for(n=0; n<grid->ncells; n++)
    a[n] = a[n] + cB[substep]*dt*at[n];

  substep = (substep+1) % 3;

  if(substep < 3)
    for(n=0; n<grid->ncells; n++)
      at[n] = cA[substep]*at[n];

  return substep;
}
