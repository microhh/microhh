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
  // rk3((*fields->u).data, (*fields->ut).data, dt);
  // rk3((*fields->v).data, (*fields->vt).data, dt);
  // rk3((*fields->w).data, (*fields->wt).data, dt);
  // rk3((*fields->s).data, (*fields->st).data, dt);
  // substep = (substep+1) % 3;

  rk4((*fields->u).data, (*fields->ut).data, dt);
  rk4((*fields->v).data, (*fields->vt).data, dt);
  rk4((*fields->w).data, (*fields->wt).data, dt);
  rk4((*fields->s).data, (*fields->st).data, dt);
  substep = (substep+1) % 5;

  return substep;
}

double ctimeint::getsubdt(double dt)
{
  double subdt;
  // subdt = rk3subdt(dt);
  subdt = rk4subdt(dt);
  return subdt;
}

double ctimeint::rk3subdt(double dt)
{
  const double cB [] = {1./3., 15./16., 8./15.};
  return cB[substep]*dt;
}

double ctimeint::rk4subdt(double dt)
{
  const double cB [] = {
    1432997174477./ 9575080441755.,
    5161836677717./13612068292357.,
    1720146321549./ 2090206949498.,
    3134564353537./ 4481467310338.,
    2277821191437./14882151754819.};

  return cB[substep]*dt;
}

int ctimeint::rk3(double * restrict a, double * restrict at, double dt)
{
  const double cA [] = {0., -5./9., -153./128.};
  const double cB [] = {1./3., 15./16., 8./15.};
  
  int i,j,k;
  int ijk,ii,jj,kk;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(k=grid->kstart; k<grid->kend; k++)
    for(j=grid->jstart; j<grid->jend; j++)
      for(i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        a[ijk] = a[ijk] + cB[substep]*dt*at[ijk];
      }

  int substepn = (substep+1) % 3;

  // substep 0 resets the tendencies, because cA[0] == 0
  for(k=grid->kstart; k<grid->kend; k++)
    for(j=grid->jstart; j<grid->jend; j++)
      for(i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        at[ijk] = cA[substep]*at[ijk];
      }

  return 0;
}

int ctimeint::rk4(double * restrict a, double * restrict at, double dt)
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

  int i,j,k;
  int ijk,ii,jj,kk;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(k=grid->kstart; k<grid->kend; k++)
    for(j=grid->jstart; j<grid->jend; j++)
      for(i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        a[ijk] = a[ijk] + cB[substep]*dt*at[ijk];
      }

  int substepn = (substep+1) % 5;

  // substep 0 resets the tendencies, because cA[0] == 0
  for(k=grid->kstart; k<grid->kend; k++)
    for(j=grid->jstart; j<grid->jend; j++)
      for(i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        at[ijk] = cA[substepn]*at[ijk];
      }

  return 0;
}

bool ctimeint::insubstep()
{
  if(substep > 0)
    return true;
  else
    return false;
}

