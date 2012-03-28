#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "timeint.h"

#define restrict __restrict__

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
  // return rk3(fields->flow, fields->flowt, fields->scal, fields->scalt, dt);
  return rk4(fields->flow, fields->flowt, fields->scal, fields->scalt, dt);
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

int ctimeint::rk3(double * restrict flow, double * restrict flowt, double * restrict scal, double * restrict scalt, double dt)
{
  const double cA [] = {0., -5./9., -153./128.};
  const double cB [] = {1./3., 15./16., 8./15.};
  // const double cdt[] = {0., 1./3., 3./4.};
  
  int i,j,k,n;
  int ijkn,ii,jj,kk,nn;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  nn = grid->ncells;

  for(n=0; n<3; n++)
    for(k=grid->kstart; k<grid->kend; k++)
      for(j=grid->jstart; j<grid->jend; j++)
        for(i=grid->istart; i<grid->iend; i++)
        {
          ijkn = i + j*jj + k*kk + n*nn;
          flow[ijkn] = flow[ijkn] + cB[substep]*dt*flowt[ijkn];
        }

  for(n=0; n<1; n++)
    for(k=grid->kstart; k<grid->kend; k++)
      for(j=grid->jstart; j<grid->jend; j++)
        for(i=grid->istart; i<grid->iend; i++)
        {
          ijkn = i + j*jj + k*kk + n*nn;
          scal[ijkn] = scal[ijkn] + cB[substep]*dt*scalt[ijkn];
        }

  substep = (substep+1) % 3;

  // substep 0 resets the tendencies, because cA[0] == 0
  for(n=0; n<3; n++)
    for(k=grid->kstart; k<grid->kend; k++)
      for(j=grid->jstart; j<grid->jend; j++)
        for(i=grid->istart; i<grid->iend; i++)
        {
          ijkn = i + j*jj + k*kk + n*nn;
          flowt[ijkn] = cA[substep]*flowt[ijkn];
        }

  for(n=0; n<1; n++)
    for(k=grid->kstart; k<grid->kend; k++)
      for(j=grid->jstart; j<grid->jend; j++)
        for(i=grid->istart; i<grid->iend; i++)
        {
          ijkn = i + j*jj + k*kk + n*nn;
          scalt[ijkn] = cA[substep]*scalt[ijkn];
        }

  return substep;
}

int ctimeint::rk4(double * restrict flow, double * restrict flowt, double * restrict scal, double * restrict scalt, double dt)
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
    flow[n] = flow[n] + cB[substep]*dt*flowt[n];

  for(n=0; n<grid->ncells; n++)
    scal[n] = scal[n] + cB[substep]*dt*scalt[n];

  substep = (substep+1) % 5;

  // substep 0 resets the tendencies, because cA[0] == 0
  for(n=0; n<grid->ncells*3; n++)
    flowt[n] = cA[substep]*flowt[n];

  for(n=0; n<grid->ncells; n++)
    scalt[n] = cA[substep]*scalt[n];

  return substep;
}

bool ctimeint::insubstep()
{
  if(substep > 0)
    return true;
  else
    return false;
}

