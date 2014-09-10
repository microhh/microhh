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

#ifdef PARALLEL
#include <mpi.h>
#include <stdexcept>
#include "grid.h"
#include "defines.h"
#include "master.h"

cmaster::cmaster()
{
  initialized = false;
  allocated   = false;

  // set the mpiid, to ensure that errors can be written if MPI init fails
  mpiid = 0;
}

cmaster::~cmaster()
{
  if(allocated)
  {
    delete[] reqs;
    MPI_Comm_free(&commxy);
    MPI_Comm_free(&commx);
    MPI_Comm_free(&commy);
  }

  printMessage("Finished run on %d processes\n", nprocs);

  if(initialized)
    MPI_Finalize();
}

void cmaster::startup(int argc, char *argv[])
{
  int n;

  // initialize the MPI
  n = MPI_Init(NULL, NULL);
  if(checkerror(n))
    throw 1;

  initialized = true;

  // get the rank of the current process
  n = MPI_Comm_rank(MPI_COMM_WORLD, &mpiid);
  if(checkerror(n))
    throw 1;

  // get the total number of processors
  n = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if(checkerror(n))
    throw 1;

  // store a temporary copy of COMM_WORLD in commxy
  n = MPI_Comm_dup(MPI_COMM_WORLD, &commxy);
  if(checkerror(n))
    throw 1;

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
    printError("nprocs = %d does not equal npx*npy = %d*%d\n", nprocs, npx, npy);
    throw 1;
  }

  int n;
  int dims    [2] = {npy, npx};
  int periodic[2] = {true, true};

  // define the dimensions of the 2-D grid layout
  n = MPI_Dims_create(nprocs, 2, dims);
  if(checkerror(n))
    throw 1;

  // create a 2-D grid communicator that is optimized for grid to grid transfer
  // first, free our temporary copy of COMM_WORLD
  n = MPI_Comm_free(&commxy);
  if(checkerror(n))
    throw 1;

  // for now, do not reorder processes, blizzard gives large performance loss
  n = MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, false, &commxy);
  if(checkerror(n))
    throw 1;

  n = MPI_Comm_rank(commxy, &mpiid);
  if(checkerror(n))
    throw 1;

  // retrieve the x- and y-coordinates in the 2-D grid for each process
  int mpicoords[2];
  n = MPI_Cart_coords(commxy, mpiid, 2, mpicoords);
  if(checkerror(n))
    throw 1;

  mpicoordx = mpicoords[1];
  mpicoordy = mpicoords[0];

  int dimx[2] = {false, true };
  int dimy[2] = {true , false};

  n = MPI_Cart_sub(commxy, dimx, &commx);
  if(checkerror(n))
    throw 1;

  n = MPI_Cart_sub(commxy, dimy, &commy);
  if(checkerror(n))
    throw 1;

  // find out who are the neighbors of this process to facilitate the communication routines
  n = MPI_Cart_shift(commxy, 1, 1, &nwest , &neast );
  if(checkerror(n))
    throw 1;

  n = MPI_Cart_shift(commxy, 0, 1, &nsouth, &nnorth);
  if(checkerror(n))
    throw 1;

  // create the requests arrays for the nonblocking sends
  int npmax;
  npmax = std::max(npx, npy);

  // have at least as many communicators as prognostic variables
  npmax = std::max(npmax, 8*4);
  reqs  = new MPI_Request[npmax*2];
  reqsn = 0;

  allocated = true;
}

double cmaster::gettime()
{
  return MPI_Wtime();
}

int cmaster::checkerror(int n)
{
  char errbuffer[MPI_MAX_ERROR_STRING];
  int errlen;

  if(n != MPI_SUCCESS)
  {
    MPI_Error_string(n, errbuffer, &errlen);
    printError("MPI: %s\n", errbuffer);
    return 1;
  }

  return 0;
}

int cmaster::waitall()
{
  // wait for MPI processes and reset the number of pending requests
  MPI_Waitall(reqsn, reqs, MPI_STATUSES_IGNORE);
  reqsn = 0;

  return 0;
}

// do all broadcasts over the MPI_COMM_WORLD, to avoid complications in the input file reading
int cmaster::broadcast(char *data, int datasize)
{
  MPI_Bcast(data, datasize, MPI_CHAR, 0, commxy);
  return 0;
}

// overloaded broadcast functions
int cmaster::broadcast(int *data, int datasize)
{
  MPI_Bcast(data, datasize, MPI_INT, 0, commxy);
  return 0;
}

int cmaster::broadcast(unsigned long *data, int datasize)
{
  MPI_Bcast(data, datasize, MPI_UNSIGNED_LONG, 0, commxy);
  return 0;
}

int cmaster::broadcast(double *data, int datasize)
{
  MPI_Bcast(data, datasize, MPI_DOUBLE, 0, commxy);
  return 0;
}

int cmaster::sum(int *var, int datasize)
{
  MPI_Allreduce(MPI_IN_PLACE, var, datasize, MPI_INT, MPI_SUM, commxy);
  return 0;
}

int cmaster::sum(double *var, int datasize)
{
  MPI_Allreduce(MPI_IN_PLACE, var, datasize, MPI_DOUBLE, MPI_SUM, commxy);
  return 0;
}

int cmaster::max(double *var, int datasize)
{
  MPI_Allreduce(MPI_IN_PLACE, var, datasize, MPI_DOUBLE, MPI_MAX, commxy);
  return 0;
}

int cmaster::min(double *var, int datasize)
{
  MPI_Allreduce(MPI_IN_PLACE, var, datasize, MPI_DOUBLE, MPI_MIN, commxy);
  return 0;
}
#endif
