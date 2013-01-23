#include <cstdio>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "timeloop.h"
#include "defines.h"

ctimeloop::ctimeloop(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  // std::printf("Creating instance of object timeloop\n");
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;

  substep = 0;
  ifactor = 1e6;
}

ctimeloop::~ctimeloop()
{
  // std::printf("Destroying instance of object timeloop\n");
}

int ctimeloop::readinifile(cinput *inputin)
{
  // input parameters
  int n = 0;

  // obligatory parameters
  n += inputin->getItem(&maxiter     , "time", "maxiter"  );
  n += inputin->getItem(&startiter   , "time", "startiter");

  // optional parameters
  n += inputin->getItem(&adaptivestep, "time", "adaptivestep", true);
  n += inputin->getItem(&dt          , "time", "dt"          , 0.1 );
  n += inputin->getItem(&dtmax       , "time", "dtmax"       , dbig);
  n += inputin->getItem(&cflmax      , "time", "cflmax"      , 1.  );
  n += inputin->getItem(&dnmax       , "time", "dnmax"       , 0.5 );
  n += inputin->getItem(&rkorder     , "time", "rkorder"     , 4   );
  n += inputin->getItem(&outputiter  , "time", "outputiter"  , 100 );
  n += inputin->getItem(&saveiter    , "time", "saveiter"    , 500 );
  n += inputin->getItem(&statsiter   , "time", "statsiter"   , 100 );
  n += inputin->getItem(&postprociter, "time", "postprociter", 100 );

  // if one argument fails, then crash
  if(n > 0)
    return 1;

  // the maximum iteration is relative to the start iteration
  maxiter += startiter;

  // the current iteration is set to the start iteration
  iteration = startiter;

  // 3 and 4 are the only valid values for the rkorder
  if(!(rkorder == 3 || rkorder == 4))
  {
    if(mpi->mpiid == 0) std::printf("ERROR \"%d\" is an illegal value for rkorder\n", rkorder);
    return 1;
  }

  // initializations
  loop      = true;
  time      = 0.;

  itime    = (unsigned long)(ifactor * time);
  idt      = (unsigned long)(ifactor * dt);
  idtmax   = (unsigned long)(ifactor * dtmax);

  gettimeofday(&start, NULL);

  return 0;
}

int ctimeloop::timestep()
{
  time  += dt;
  itime += idt;

  iteration++;

  if(iteration >= maxiter)
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
  // do not save directly after the start of the simulation
  if(iteration % saveiter == 0 && iteration != startiter)
    return 1;

  return 0;
}

int ctimeloop::dostats()
{
  if(iteration % statsiter == 0)
  {
    // do not save directly after the start of the simulation, because it has been done
    // at the end of the previous run, except for iteration 0
    if(iteration == startiter && iteration != 0)
      return 0;
    return 1;
  }

  return 0;
}


double ctimeloop::check()
{
  gettimeofday(&end, NULL);

  double timeelapsed = (double)(end.tv_sec-start.tv_sec) + (double)(end.tv_usec-start.tv_usec) * 1.e-6;
  start = end;

  return timeelapsed;
}

int ctimeloop::settimestep(double cfl, double dn)
{
  if(adaptivestep)
  {
    idt = (long)(std::min((double)idt * cflmax/cfl, (double)idt * dnmax/dn));
    idt = std::min(idt, idtmax);
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
  double subdt = 0.;
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

int ctimeloop::save(int n)
{
  if(mpi->mpiid == 0)
  {
    char filename[256];
    std::sprintf(filename, "time.%07d", n);

    std::printf("Saving \"%s\"\n", filename);

    FILE *pFile;
    pFile = fopen(filename, "wb");

    if(pFile == NULL)
    {
      std::printf("ERROR \"%s\" cannot be written", filename);
      return 1;
    }

    fwrite(&itime, sizeof(long), 1, pFile);
    fwrite(&idt  , sizeof(long), 1, pFile);

    fclose(pFile);
  }

  return 0;
}

int ctimeloop::load(int n)
{
  if(mpi->mpiid == 0)
  {
    char filename[256];
    std::sprintf(filename, "time.%07d", n);

    std::printf("Loading \"%s\"\n", filename);

    FILE *pFile;
    pFile = fopen(filename, "rb");

    if(pFile == NULL)
    {
      if(mpi->mpiid == 0) std::printf("ERROR \"%s\" does not exist\n", filename);

      return 1;
    }

    fread(&itime, sizeof(long), 1, pFile);
    fread(&idt  , sizeof(long), 1, pFile);

    fclose(pFile);
  }

  mpi->broadcast(&itime, 1);
  mpi->broadcast(&idt  , 1);

  // calculate the double precision time from the integer time
  time = (double)itime / ifactor;
  dt   = (double)idt   / ifactor;

  return 0;
}

int ctimeloop::postprocstep()
{
  iteration += postprociter;

  if(iteration > maxiter)
    loop = false;

  return 0;
}
