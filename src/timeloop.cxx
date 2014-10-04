/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdio>
#include <cmath>
#include "input.h"
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "timeloop.h"
#include "defines.h"
#include "constants.h"
#include "model.h"

Timeloop::Timeloop(Model *modelin, Input *inputin)
{
  model  = modelin;
  grid   = model->grid;
  fields = model->fields;
  master = model->master;

  substep = 0;
  ifactor = 1e6;

  // input parameters
  int n = 0;

  // obligatory parameters
  if(master->mode == "init")
    starttime = 0.;
  else
    n += inputin->getItem(&starttime, "time", "starttime", "");

  n += inputin->getItem(&endtime , "time", "endtime" , "");
  n += inputin->getItem(&savetime, "time", "savetime", "");

  // optional parameters
  n += inputin->getItem(&adaptivestep, "time", "adaptivestep", "", true            );
  n += inputin->getItem(&dtmax       , "time", "dtmax"       , "", constants::dbig );
  n += inputin->getItem(&dt          , "time", "dt"          , "", dtmax           );
  n += inputin->getItem(&rkorder     , "time", "rkorder"     , "", 3               );
  n += inputin->getItem(&outputiter  , "time", "outputiter"  , "", 20              );
  n += inputin->getItem(&iotimeprec  , "time", "iotimeprec"  , "", 0               );

  if(master->mode == "post")
    n += inputin->getItem(&postproctime, "time", "postproctime", "");

  // if one argument fails, then crash
  if(n > 0)
    throw 1;

  // 3 and 4 are the only valid values for the rkorder
  if(!(rkorder == 3 || rkorder == 4))
  {
    master->printError("\"%d\" is an illegal value for rkorder\n", rkorder);
    throw 1;
  }

  // initializations
  loop      = true;
  time      = 0.;
  iteration = 0;

  // set or calculate all the integer times
  itime         = (unsigned long) 0;

  // add 0.5 to prevent roundoff errors
  iendtime      = (unsigned long)(ifactor * endtime + 0.5);
  istarttime    = (unsigned long)(ifactor * starttime + 0.5);
  idt           = (unsigned long)(ifactor * dt + 0.5);
  idtmax        = (unsigned long)(ifactor * dtmax + 0.5);
  isavetime     = (unsigned long)(ifactor * savetime + 0.5);
  if(master->mode == "post")
    ipostproctime = (unsigned long)(ifactor * postproctime + 0.5);

  idtlim = idt;

  // take the proper precision for the output files into account
  iiotimeprec = (unsigned long)(ifactor * std::pow(10., iotimeprec) + 0.5);

  // check whether starttime and savetime are an exact multiple of iotimeprec
  if((istarttime % iiotimeprec) || (isavetime % iiotimeprec))
  {
    master->printError("starttime or savetime is not an exact multiple of iotimeprec\n");
    throw 1;
  }

  iotime = (int)(istarttime / iiotimeprec);

  gettimeofday(&start, NULL);
}

Timeloop::~Timeloop()
{
}

int Timeloop::settimelim()
{
  idtlim = idtmax;
  idtlim = std::min(idtlim, isavetime - itime % isavetime);

  return 0;
}

void Timeloop::stepTime()
{
  // Only step forward in time if we are not in a substep
  if(inSubStep())
    return;

  time  += dt;
  itime += idt;
  iotime = (int)(itime/iiotimeprec);

  ++iteration;

  if(itime >= iendtime)
    loop = false;
}

bool Timeloop::doCheck()
{
  if(iteration % outputiter == 0 && !inSubStep())
    return true;

  return false;
}

bool Timeloop::doSave()
{
  // do not save directly after the start of the simulation and not in a substep
  if(itime % isavetime == 0 && iteration != 0 && !inSubStep())
    return true;

  return false;
}

double Timeloop::check()
{
  gettimeofday(&end, NULL);

  double timeelapsed = (double)(end.tv_sec-start.tv_sec) + (double)(end.tv_usec-start.tv_usec) * 1.e-6;
  start = end;

  return timeelapsed;
}

void Timeloop::setTimeStep()
{
  // Only set the time step if we are not in a substep
  if(inSubStep())
    return;

  if(adaptivestep)
  {
    idt = idtlim;
    dt  = (double)idt / ifactor;
  }
}

#ifdef USECUDA
int Timeloop::exec()
{
  if(rkorder == 3)
  {
    for(fieldmap::const_iterator it = fields->at.begin(); it!=fields->at.end(); ++it)
      rk3_GPU(fields->ap[it->first]->data_g, it->second->data_g, dt);

    substep = (substep+1) % 3;
  }

  if(rkorder == 4)
  {
    for(fieldmap::const_iterator it = fields->at.begin(); it!=fields->at.end(); ++it)
      rk4_GPU(fields->ap[it->first]->data_g, it->second->data_g, dt);

    substep = (substep+1) % 5;
  }
  return substep;
}
#else
int Timeloop::exec()
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
#endif

double Timeloop::getSubTimeStep()
{
  double subdt = 0.;
  if(rkorder == 3)
    subdt = rk3subdt(dt);
  else if(rkorder == 4)
    subdt = rk4subdt(dt);

  return subdt;
}

double Timeloop::rk3subdt(double dt)
{
  const double cB [] = {1./3., 15./16., 8./15.};
  return cB[substep]*dt;
}

double Timeloop::rk4subdt(double dt)
{
  const double cB [] = {
    1432997174477./ 9575080441755.,
    5161836677717./13612068292357.,
    1720146321549./ 2090206949498.,
    3134564353537./ 4481467310338.,
    2277821191437./14882151754819.};

  return cB[substep]*dt;
}

int Timeloop::rk3(double * restrict a, double * restrict at, double dt)
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

int Timeloop::rk4(double * restrict a, double * restrict at, double dt)
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

bool Timeloop::inSubStep()
{
  if(substep > 0)
    return true;
  else
    return false;
}

bool Timeloop::isStatsStep()
{
  // In case we are not in a substep and not at the first iteration
  // after a restart, we can could do statistics.
  if(!inSubStep() && !((iteration > 0) && (itime == istarttime)))
    return true;
  else
    return false;
}

void Timeloop::save(int starttime)
{
  if(master->mpiid == 0)
  {
    char filename[256];
    std::sprintf(filename, "time.%07d", starttime);

    master->printMessage("Saving \"%s\" ... ", filename);

    FILE *pFile;
    pFile = fopen(filename, "wb");

    if(pFile == NULL)
    {
      master->printError("\"%s\" cannot be written", filename);
      throw 1;
    }

    fwrite(&itime    , sizeof(long), 1, pFile);
    fwrite(&idt      , sizeof(long), 1, pFile);
    fwrite(&iteration, sizeof(long), 1, pFile);

    fclose(pFile);
    master->printMessage("OK\n");
  }
}

void Timeloop::load(int starttime)
{
  int nerror = 0;

  if(master->mpiid == 0)
  {
    char filename[256];
    std::sprintf(filename, "time.%07d", starttime);

    master->printMessage("Loading \"%s\" ... ", filename);

    FILE *pFile;
    pFile = fopen(filename, "rb");

    if(pFile == NULL)
    {
      master->printError("\"%s\" does not exist\n", filename);
      ++nerror;
    }
    else
    {
      fread(&itime    , sizeof(long), 1, pFile);
      fread(&idt      , sizeof(long), 1, pFile);
      fread(&iteration, sizeof(long), 1, pFile);

      fclose(pFile);
    }
    master->printMessage("OK\n");
  }

  master->broadcast(&nerror, 1);
  if(nerror)
    throw 1;

  master->broadcast(&itime    , 1);
  master->broadcast(&idt      , 1);
  master->broadcast(&iteration, 1);

  // calculate the double precision time from the integer time
  time = (double)itime / ifactor;
  dt   = (double)idt   / ifactor;
}

int Timeloop::postprocstep()
{
  itime += ipostproctime;
  iotime = (int)(itime/iiotimeprec);

  if(itime > iendtime)
    loop = false;

  return 0;
}
