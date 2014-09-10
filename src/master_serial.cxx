/*
 * MicroHH
 * Copyright (c) 2011-2013 Chiel van Heerwaarden
 * Copyright (c) 2011-2013 Thijs Heus
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

#ifndef PARALLEL
#include <sys/time.h>
#include "grid.h"
#include "defines.h"
#include "master.h"

cmaster::cmaster()
{
  initialized = false;
  allocated   = false;
}

cmaster::~cmaster()
{
  printMessage("Finished run on %d processes\n", nprocs);
}

void cmaster::startup(int argc, char *argv[])
{
  initialized = true;

  // set the rank of the only process to 0
  mpiid = 0;
  // set the number of processes to 1
  nprocs = 1;

  printMessage("Starting run on %d processes\n", nprocs);

  // process the command line options
  if(argc <= 1)
  {
    printError("Specify init, run or post mode\n");
    throw 1;
  }
  else
  {
    // check the execution mode
    mode = argv[1];
    if(mode != "init" && mode != "run" && mode != "post")
    {
      printError("Specify init, run or post mode\n");
      throw 1;
    }
    // set the name of the simulation
    if(argc > 2)
      simname = argv[2];
    else
      simname = "microhh";
  }
}

void cmaster::init(cinput *inputin)
{
  int nerror = 0;
  nerror += inputin->getItem(&npx, "mpi", "npx", "", 1);
  nerror += inputin->getItem(&npy, "mpi", "npy", "", 1);
  if(nerror)
    throw 1;

  if(nprocs != npx*npy)
  {
    printError("npx*npy = %d*%d has to be equal to 1*1 in serial mode\n", npx, npy);
    throw 1;
  }

  // set the coordinates to 0
  mpicoordx = 0;
  mpicoordy = 0;

  allocated = true;
}

double cmaster::gettime()
{
  timeval timestruct;
  gettimeofday(&timestruct, NULL);
  double time;
  time = (double)timestruct.tv_sec + (double)timestruct.tv_usec*1.e-6;
  return time;
}

int cmaster::waitall()
{
  return 0;
}

// all broadcasts return directly, because there is nothing to broadcast
int cmaster::broadcast(char *data, int datasize)
{
  return 0;
}

// overloaded broadcast functions
int cmaster::broadcast(int *data, int datasize)
{
  return 0;
}

int cmaster::broadcast(unsigned long *data, int datasize)
{
  return 0;
}

int cmaster::broadcast(double *data, int datasize)
{
  return 0;
}

int cmaster::sum(int *var, int datasize)
{
  return 0;
}

int cmaster::sum(double *var, int datasize)
{
  return 0;
}

int cmaster::max(double *var, int datasize)
{
  return 0;
}

int cmaster::min(double *var, int datasize)
{
  return 0;
}
#endif
