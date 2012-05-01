#include <cstdio>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "pres.h"
#include "timeloop.h"

int main()
{
  // create the class objects
  cinput    input;
  cgrid     grid;
  cmpi      mpi     (&grid);
  cfields   fields  (&grid, &mpi);
  cpres     pres    (&grid, &fields);
  ctimeloop timeloop(&grid, &fields);
  
  // read the input data and terminate on error
  if(input.readinifile())
    return 1;
  if(grid.readinifile(&input))
    return 1;
  if(mpi.readinifile(&input))
    return 1;
  if(fields.readinifile(&input))
    return 1;
  if(timeloop.readinifile(&input))
    return 1;
  
  // initialize the objects, allocate the required memory
  grid.init(mpi.npx, mpi.npy);
  if(mpi.init())
    return 1;
  fields.initfields();

  // fill the fields with data
  grid.create();
  fields.createfields();

  // store the data on disk
  if(grid.save(mpi.mpiid))
    return 1;
  if(pres.save(mpi.mpiid))
    return 1;
  if(fields.save(0, mpi.mpiid))
    return 1;
  if(timeloop.save(0, mpi.mpiid))
    return 1;

  return 0;
}

