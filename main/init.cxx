#include <cstdio>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "buffer.h"
#include "timeloop.h"

int main(int argc, char *argv[])
{
  // set the name of the simulation
  std::string simname("microhh");
  if(argc > 1)
    simname = argv[1];

  // start the MPI
  cmpi      mpi;
  mpi.startup();

  // create the class objects
  cinput    input   (&mpi);
  cgrid     grid    (&mpi);
  cfields   fields  (&grid, &mpi);
  cbuffer   buffer  (&grid, &fields, &mpi);
  ctimeloop timeloop(&grid, &fields, &mpi);
  
  // read the input data
  if(input.readinifile(simname))
    return 1;
  // read the profiles
  if(input.readproffile(simname))
    return 1;

  // process the settings data
  if(grid.readinifile(&input))
    return 1;
  if(mpi.readinifile(&input))
    return 1;
  if(fields.readinifile(&input))
    return 1;
  if(buffer.readinifile(&input))
    return 1;
  if(timeloop.readinifile(&input))
    return 1;
  
  // initialize the objects, allocate the required memory
  if(mpi.init())
    return 1;
  if(grid.init())
    return 1;
  if(fields.init())
    return 1;
  if(buffer.init())
    return 1;

  // read the grid from the input
  if(grid.create(&input))
    return 1;

  // create the random field
  // create the initial field
  if(fields.create(&input))
    return 1;

  if(buffer.setbuffers())
    return 1;

  // free the memory of the input data
  input.clear();

  // store the data on disk
  if(grid.save())
    return 1;
  if(buffer.save())
    return 1;
  if(fields.save(timeloop.iteration))
    return 1;
  if(timeloop.save(timeloop.iteration))
    return 1;

  return 0;
}

