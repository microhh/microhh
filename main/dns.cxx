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

  // free the memory of the input
  input.clear();

  // fill the fields with data
  if(grid.load())
    return 1;
  // CvH, the model loads the fields now, not good...
  if(model.load())
    return 1;

  // run the model
  model.exec();

  return 0;
}
