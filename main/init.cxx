#include <cstdio>
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "pres.h"
#include "timeloop.h"

int main(int argc, char *argv[])
{
  // set the name of the simulation
  std::string simname("microhh");
  if(argc > 1)
    simname = argv[1];

  // create the class objects
  cinput    input;
  cmpi      mpi;
  cgrid     grid    (&mpi);
  cfields   fields  (&grid, &mpi);
  cpres     pres    (&grid, &fields, &mpi);
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
  if(timeloop.readinifile(&input))
    return 1;
  
  // initialize the objects, allocate the required memory
  if(mpi.init())
    return 1;
  if(grid.init())
    return 1;
  fields.initfields();

  // INPUT FILE CREATION
  // read the grid from the input
  grid.create();

  // create the random field
  // create the initial field
  fields.createfields(&input);

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

