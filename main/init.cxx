#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "model.h"

int main(int argc, char *argv[])
{
  // set the name of the simulation
  std::string simname("microhh");
  if(argc > 1)
    simname = argv[1];

  // start up the message passing interface
  cmpi mpi;
  mpi.startup();

  // create the instances of the objects
  cinput  input (&mpi);
  cgrid   grid  (&mpi);
  cfields fields(&grid, &mpi);

  cmodel  model (&grid, &fields, &mpi, simname);

  // read the input data
  if(input.readinifile(simname))
    return 1;
  if(input.readproffile(simname))
    return 1;

  if(mpi.readinifile(&input))
    return 1;
  if(grid.readinifile(&input))
    return 1;
  if(fields.readinifile(&input))
    return 1;
  if(model.readinifile(&input))
    return 1;

  // init the mpi 
  if(mpi.init())
    return 1;
  if(grid.init())
    return 1;
  if(fields.init())
    return 1;
  if(model.init())
    return 1;

  // read the grid from the input
  if(grid.create(&input))
    return 1;
  if(fields.create(&input))
    return 1;

  // save the data
  if(grid.save())
    return 1;
  // CvH, the model saves the fields now, not good...
  if(model.save())
    return 1;

  // free the memory of the input
  input.clear();

  return 0;
}
