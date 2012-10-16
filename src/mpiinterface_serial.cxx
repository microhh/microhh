#ifndef PARALLEL
#include <cstdio>
#include <sys/time.h>
#include "grid.h"
#include "defines.h"
#include "mpiinterface.h"

cmpi::cmpi()
{
  initialized = false;
  allocated   = false;
}

cmpi::~cmpi()
{
  if(mpiid == 0) std::printf("Finished run on %d processes\n", nprocs);
}

int cmpi::readinifile(cinput *inputin)
{
  int n = 0;

  n += inputin->getItem(&npx, "grid", "npx", 1);
  n += inputin->getItem(&npy, "grid", "npy", 1);

  if(n > 0)
    return 1;
  
  return 0;
}

int cmpi::startup()
{
  initialized = true;

  // set the rank of the only process to 0
  mpiid = 0;
  // set the number of processes to 1
  nprocs = 1;

  if(mpiid == 0) std::printf("Starting run on %d processes\n", nprocs);

  return 0;
}

int cmpi::init()
{
  if(nprocs != npx*npy)
  {
    if(mpiid == 0) std::printf("ERROR nprocs = %d does not equal npx*npy = %d*%d\n", nprocs, npx, npy);
    return 1;
  }

  // set the coordinates to 0
  mpicoordx = 0;
  mpicoordy = 0;

  allocated = true;

  return 0;
}

double cmpi::gettime()
{
  timeval timestruct;
  gettimeofday(&timestruct, NULL);
  double time;
  time = (double)timestruct.tv_sec + (double)timestruct.tv_usec*1.e-6;
  return time;
}

int cmpi::waitall()
{
  return 0;
}

// all broadcasts return directly, because there is nothing to broadcast
int cmpi::broadcast(char *data, int datasize)
{
  return 0;
}

// overloaded broadcast functions
int cmpi::broadcast(int *data, int datasize)
{
  return 0;
}

int cmpi::broadcast(unsigned long *data, int datasize)
{
  return 0;
}

int cmpi::broadcast(double *data, int datasize)
{
  return 0;
}
#endif
