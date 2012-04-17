#include <cstdio>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "timeloop.h"

ctimeloop::ctimeloop(cgrid *gridin, cfields *fieldsin)
{
  std::printf("Creating instance of object timeloop\n");
  grid   = gridin;
  fields = fieldsin;

  substep = 0;
}

ctimeloop::~ctimeloop()
{
  std::printf("Destroying instance of object timeloop\n");
}

int ctimeloop::readinifile(cinput *inputin)
{
  // input parameters
  int n = 0;

  // obligatory parameters
  n += inputin->getItem(&adaptivestep, "time", "adaptivestep", true);
  n += inputin->getItem(&runtime     , "time", "runtime"     , true);
  n += inputin->getItem(&cflmax      , "time", "cflmax"      , true);

  if(n > 0)
    return 1;

  // optional parameters
  n = inputin->getItem(&iteration, "time", "iteration", false);
  if(n > 0)
    iteration = 0;

  // initializations
  loop      = true;
  time      = 0.;
  dt        = 0.1;
  cflmax    = 1.5;

  const int ifactor = 1000;

  itime    = (int)(ifactor * time);
  iruntime = (int)(ifactor * runtime);
  idt      = (int)(ifactor * dt);

  gettimeofday(&start, NULL);

  return 0;
}

int ctimeloop::timestep()
{
  time  += dt;
  itime += idt;

  iteration++;

  if(time >= runtime)
    loop = false;

  if(iteration % 100 == 0) 
  {
    gettimeofday(&end, NULL);

    double timeelapsed;

    timeelapsed = (double)(end.tv_sec-start.tv_sec) + (double)(end.tv_usec-start.tv_usec) * 1.e-6;

    std::printf("Iteration = %6d, runtime = %7.1f, cputime = %10.5f\n", iteration, time, timeelapsed);

    start = end;
  }

  if(iteration % 500 == 0) 
  {
    (*fields->u).save(iteration);
    (*fields->v).save(iteration);
    (*fields->w).save(iteration);
    (*fields->p).save(iteration);
    (*fields->s).save(iteration);
  }

  return 0;
}

int ctimeloop::settimestep(double cfl)
{
  if(adaptivestep)
    dt = dt * cflmax/cfl;

  std::printf("CFL = %f, dt = %f\n", cfl, dt);

  return 0;
}

int ctimeloop::exec()
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

double ctimeloop::getsubdt()
{
  double subdt;
  // subdt = rk3subdt(dt);
  subdt = rk4subdt(dt);
  return subdt;
}

double ctimeloop::rk3subdt(double dt)
{
  const double cB [] = {1./3., 15./16., 8./15.};
  return cB[substep]*dt;
}

double ctimeloop::rk4subdt(double dt)
{
  const double cB [] = {
    1432997174477./ 9575080441755.,
    5161836677717./13612068292357.,
    1720146321549./ 2090206949498.,
    3134564353537./ 4481467310338.,
    2277821191437./14882151754819.};

  return cB[substep]*dt;
}

int ctimeloop::rk3(double * restrict a, double * restrict at, double dt)
{
  const double cA [] = {0., -5./9., -153./128.};
  const double cB [] = {1./3., 15./16., 8./15.};
  
  int i,j,k;
  int ijk,jj,kk;

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
        at[ijk] = cA[substepn]*at[ijk];
      }

  return 0;
}

int ctimeloop::rk4(double * restrict a, double * restrict at, double dt)
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
  int ijk,jj,kk;

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

bool ctimeloop::insubstep()
{
  if(substep > 0)
    return true;
  else
    return false;
}
