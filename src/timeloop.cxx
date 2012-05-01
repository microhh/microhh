#include <cstdio>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "timeloop.h"
#include "defines.h"

ctimeloop::ctimeloop(cgrid *gridin, cfields *fieldsin)
{
  std::printf("Creating instance of object timeloop\n");
  grid   = gridin;
  fields = fieldsin;

  substep = 0;
  ifactor = 10000.;
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
  n += inputin->getItem(&adaptivestep, "time", "adaptivestep");
  n += inputin->getItem(&runtime     , "time", "runtime"     );

  // optional parameters
  n += inputin->getItem(&cflmax    , "time", "cflmax"     , 1. );
  n += inputin->getItem(&maxiter   , "time", "maxiter"    , 1e9);
  n += inputin->getItem(&iteration , "time", "iteration"  , 0  );
  n += inputin->getItem(&rkorder   , "time", "rkorder"    , 4  );
  n += inputin->getItem(&outputiter, "time", "outputiter" , 100);
  n += inputin->getItem(&saveiter  , "time", "saveiter"   , 500);

  // if one argument fails, then crash
  if(n > 0)
    return 1;

  // the maximum iteration is relative to the start iteration
  maxiter += iteration;

  // 3 and 4 are the only valid values for the rkorder
  if(!(rkorder == 3 || rkorder == 4))
  {
    std::printf("ERROR \"%d\" is an illegal value for rkorder\n", rkorder);
    return 1;
  }

  // initializations
  loop      = true;
  time      = 0.;
  dt        = 0.1;

  itime    = (long)(ifactor * time);
  iruntime = (long)(ifactor * runtime);
  idt      = (long)(ifactor * dt);

  gettimeofday(&start, NULL);

  return 0;
}

int ctimeloop::timestep()
{
  time  += dt;
  itime += idt;

  iteration++;

  if(itime >= iruntime || iteration >= maxiter)
    loop = false;

  return 0;
}

int ctimeloop::docheck()
{
  if(iteration % outputiter == 0)
    return 1;

  return 0;
}

int ctimeloop::dosave()
{
  if(iteration % saveiter == 0)
    return 1;

  return 0;
}

double ctimeloop::check()
{
  gettimeofday(&end, NULL);

  double timeelapsed = (double)(end.tv_sec-start.tv_sec) + (double)(end.tv_usec-start.tv_usec) * 1.e-6;
  start = end;

  return timeelapsed;
}

int ctimeloop::settimestep(double cfl)
{
  if(adaptivestep)
  {
    idt = (long)((double)idt * cflmax/cfl);
    dt  = (double)idt / ifactor;
  }

  return 0;
}

int ctimeloop::exec()
{
  if(rkorder == 3)
  {
    rk3((*fields->u).data, (*fields->ut).data, dt);
    rk3((*fields->v).data, (*fields->vt).data, dt);
    rk3((*fields->w).data, (*fields->wt).data, dt);
    rk3((*fields->s).data, (*fields->st).data, dt);
    substep = (substep+1) % 3;
  }

  if(rkorder == 4)
  {
    rk4((*fields->u).data, (*fields->ut).data, dt);
    rk4((*fields->v).data, (*fields->vt).data, dt);
    rk4((*fields->w).data, (*fields->wt).data, dt);
    rk4((*fields->s).data, (*fields->st).data, dt);
    substep = (substep+1) % 5;
  }

  return substep;
}

double ctimeloop::getsubdt()
{
  double subdt;
  if(rkorder == 3)
    subdt = rk3subdt(dt);
  else if(rkorder == 4)
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
#pragma ivdep
      for(i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        a[ijk] = a[ijk] + cB[substep]*dt*at[ijk];
      }

  int substepn = (substep+1) % 3;

  // substep 0 resets the tendencies, because cA[0] == 0
  for(k=grid->kstart; k<grid->kend; k++)
    for(j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
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
#pragma ivdep
      for(i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        a[ijk] = a[ijk] + cB[substep]*dt*at[ijk];
      }

  int substepn = (substep+1) % 5;

  // substep 0 resets the tendencies, because cA[0] == 0
  for(k=grid->kstart; k<grid->kend; k++)
    for(j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
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

int ctimeloop::save(int n, int mpiid)
{
  FILE *pFile;
  char filename[256];
  std::sprintf(filename, "time.%07d.%07d", n, mpiid);
  pFile = fopen(filename, "wb");

  if(pFile == NULL)
  {
    std::printf("ERROR \"%s\" cannot be written", filename);
    return 1;
  }
  else
    std::printf("Saving \"%s\"\n", filename);

  fwrite(&itime, sizeof(long), 1, pFile);

  fclose(pFile);

  return 0;
}

int ctimeloop::load(int n, int mpiid)
{
  FILE *pFile;
  char filename[256];
  std::sprintf(filename, "time.%07d.%07d", n, mpiid);
  pFile = fopen(filename, "rb");

  if(pFile == NULL)
  {
    std::printf("ERROR \"%s\" does not exist\n", filename);
    return 1;
  }
  else
    std::printf("Loading \"%s\"\n", filename);

  fread(&itime, sizeof(long), 1, pFile);

  fclose(pFile);

  // calculate the double precision time from the integer time
  time = (double)itime / ifactor;

  return 0;
}
