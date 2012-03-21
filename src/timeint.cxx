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
  // return rk3(fields->flow, fields->flowt, dt);
  return rk4(fields->flow, fields->flowt, dt);
}

double ctimeint::subdt(double dt)
{
  // const double cB [] = {1./3., 15./16., 8./15.};
  const double cB [] = {
    1432997174477./ 9575080441755.,
    5161836677717./13612068292357.,
    1720146321549./ 2090206949498.,
    3134564353537./ 4481467310338.,
    2277821191437./14882151754819.};
  return cB[substep]*dt;
}

int ctimeint::rk3(double * __restrict__ a, double * __restrict__ at, double dt)
{
  const double cA [] = {0., -5./9., -153./128.};
  const double cB [] = {1./3., 15./16., 8./15.};
  // const double cdt[] = {0., 1./3., 3./4.};
  
  int n;

  for(n=0; n<grid->ncells*3; n++)
    a[n] = a[n] + cB[substep]*dt*at[n];

  substep = (substep+1) % 3;

  // substep 0 resets the tendencies, because cA[0] == 0
  for(n=0; n<grid->ncells*3; n++)
    at[n] = cA[substep]*at[n];

  return substep;
}

int ctimeint::rk4(double * __restrict__ a, double * __restrict__ at, double dt)
{
  const double cA [] = {
      0.,
    - 567301805773./1357537059087.,
    -2404267990393./2016746695238.,
    -3550918686646./2091501179385.,
    -1275806237668./ 842570457699.};

  const double cB [] = {
    1432997174477./ 9575080441755.,
    5161836677717./13612068292357.,
    1720146321549./ 2090206949498.,
    3134564353537./ 4481467310338.,
    2277821191437./14882151754819.};

  int n;

  for(n=0; n<grid->ncells*3; n++)
    a[n] = a[n] + cB[substep]*dt*at[n];

  substep = (substep+1) % 5;

  // substep 0 resets the tendencies, because cA[0] == 0
  for(n=0; n<grid->ncells*3; n++)
    at[n] = cA[substep]*at[n];

  return substep;
}

bool ctimeint::insubstep()
{
  if(substep > 0)
    return true;
  else
    return false;
}

