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
  cmpi      mpi;
  cgrid     grid    (&mpi);
  cfields   fields  (&grid, &mpi);
  cpres     pres    (&grid, &fields, &mpi);
  ctimeloop timeloop(&grid, &fields, &mpi);
  
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
  if(mpi.init())
    return 1;
  if(grid.init())
    return 1;
  fields.initfields();

  // fill the fields with data
  grid.create();
  fields.createfields();

  // store the data on disk
  if(grid.save())
    return 1;
  if(pres.save())
    return 1;
  if(fields.save(0))
    return 1;
  if(timeloop.save(0))
    return 1;

  return 0;
}

