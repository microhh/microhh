#include <cstdio>
#include "grid.h"
#include "mpiinterface.h"
#include "model.h"

int main(int argc, char *argv[])
{
  // start up the message passing interface
  cmpi mpi;
  if(mpi.startup(argc, argv))
    return 1;

  // create the instances of the objects
  cinput  input (&mpi);
  cgrid   grid  (&mpi);
  cmodel  model (&grid, &mpi);

  // read the input data
  if(input.readinifile())
    return 1;
  if(mpi.mode == "init")
    if(input.readproffile())
      return 1;

  if(mpi.readinifile(&input))
    return 1;
  if(grid.readinifile(&input))
    return 1;
  if(model.readinifile(&input))
    return 1;

  // init the mpi 
  if(mpi.init())
    return 1;
  if(grid.init())
    return 1;
  if(model.init())
    return 1;

  if(mpi.mode == "init")
  {
    // read the grid from the input
    if(grid.create(&input))
      return 1;

    // save the data
    if(grid.save())
      return 1;
    if(model.save(&input))
      return 1;
  }
  else
  {
  if(grid.load())
    return 1;
  if(model.load())
    return 1;
  }

  // free the memory of the input
  input.clear();

  if(mpi.mode != "init")
  {
    // run the model
    model.exec();
  }

  return 0;
}
