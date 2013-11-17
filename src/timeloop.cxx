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
  n += inputin->getItem(&runtime, "time", "runtime", "");
  if(mpi->mode == "init")
  {
    starttime = 0.;
  }
  else
  {
    n += inputin->getItem(&starttime, "time", "starttime", "");
  }
  // optional parameters
  n += inputin->getItem(&adaptivestep, "time", "adaptivestep", "", true );
  n += inputin->getItem(&dtmax       , "time", "dtmax"       , "", dbig );
  n += inputin->getItem(&dt          , "time", "dt"          , "", dtmax);
  n += inputin->getItem(&rkorder     , "time", "rkorder"     , "", 4    );
  n += inputin->getItem(&outputiter  , "time", "outputiter"  , "", 100  );
  n += inputin->getItem(&savetime    , "time", "savetime"    , "", 3600 );
  n += inputin->getItem(&precision   , "time", "precision"   , "", 1.   );
  n += inputin->getItem(&postproctime, "time", "postproctime", "", 3600 );

  // if one argument fails, then crash
  if(n > 0)
    return 1;

  // 3 and 4 are the only valid values for the rkorder
  if(!(rkorder == 3 || rkorder == 4))
  {
    if(mpi->mpiid == 0) std::printf("ERROR \"%d\" is an illegal value for rkorder\n", rkorder);
    return 1;
  }

  // initializations
  loop      = true;
  time      = 0.;
  iteration = 0;

  // calculate all the integer times
  iruntime      = (unsigned long)(ifactor * runtime);
  istarttime    = (unsigned long)(ifactor * starttime);
  itime         = (unsigned long) 0;
  idt           = (unsigned long)(ifactor * dt);
  idtmax        = (unsigned long)(ifactor * dtmax);
  isavetime     = (unsigned long)(ifactor * savetime);
  ipostproctime = (unsigned long)(ifactor * postproctime);

  idtlim = idt;

  iotime = (int)(starttime / precision);

  gettimeofday(&start, NULL);

  return 0;
}

int ctimeloop::settimelim()
{
  idtlim = idtmax;
  idtlim = std::min(idtlim,isavetime -  itime % isavetime);

  return 0;
}

int ctimeloop::timestep()
{
  time   += dt;
  itime  += idt;
  iotime = (int)(1./precision*itime/ifactor);

  iteration++;

  if(itime >= iruntime)
  {
    loop = false;
  }
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
  if(itime % isavetime == 0 && iteration != 0)
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

int ctimeloop::settimestep()
{
  if(adaptivestep)
  {
    idt = idtlim;
    dt  = (double)idt / ifactor;
  }

  return 0;
}

int ctimeloop::exec()
{

  if(rkorder == 3)
  {
    for(fieldmap::const_iterator it = fields->at.begin(); it!=fields->at.end(); ++it)
      rk3(fields->ap[it->first]->data, it->second->data, dt);

    substep = (substep+1) % 3;
  }

  if(rkorder == 4)
  {
    for(fieldmap::const_iterator it = fields->at.begin(); it!=fields->at.end(); ++it)
      rk4(fields->ap[it->first]->data, it->second->data, dt);

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

int ctimeloop::save(int starttime)
{
  if(mpi->mpiid == 0)
  {
    char filename[256];
    std::sprintf(filename, "time.%07d", starttime);

    std::printf("Saving \"%s\" ... ", filename);

    FILE *pFile;
    pFile = fopen(filename, "wb");

    if(pFile == NULL)
    {
      std::printf("ERROR \"%s\" cannot be written", filename);
      return 1;
    }

    fwrite(&itime    , sizeof(long), 1, pFile);
    fwrite(&idt      , sizeof(long), 1, pFile);
    fwrite(&iteration, sizeof(long), 1, pFile);

    fclose(pFile);
    std::printf("OK\n");
  }

  return 0;
}

int ctimeloop::load(int starttime)
{
  int nerror = 0;

  if(mpi->mpiid == 0)
  {
    char filename[256];
    std::sprintf(filename, "time.%07d", starttime);

    std::printf("Loading \"%s\" ... ", filename);

    FILE *pFile;
    pFile = fopen(filename, "rb");

    if(pFile == NULL)
    {
      std::printf("ERROR \"%s\" does not exist\n", filename);
      ++nerror;
    }
    else
    {
      fread(&itime    , sizeof(long), 1, pFile);
      fread(&idt      , sizeof(long), 1, pFile);
      fread(&iteration, sizeof(long), 1, pFile);

      fclose(pFile);
    }
    std::printf("OK\n");
  }

  mpi->broadcast(&nerror, 1);
  if(nerror)
    return 1;

  mpi->broadcast(&itime    , 1);
  mpi->broadcast(&idt      , 1);
  mpi->broadcast(&iteration, 1);

  // calculate the double precision time from the integer time
  time = (double)itime / ifactor;
  dt   = (double)idt   / ifactor;

  return 0;
}

int ctimeloop::postprocstep()
{
  itime += ipostproctime;
  iotime = (int)(1./precision*itime/ifactor);

  if(itime > iruntime)
    loop = false;

  return 0;
}
