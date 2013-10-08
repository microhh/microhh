#include <cstdio>
#include "grid.h"
#include "mpiinterface.h"
#include "model.h"

int main(int argc, char *argv[])
{
  cmpi mpi;
  if(mpi.startup(argc, argv))
    return 1;

  // start up the message passing interface
  if(mpi.mpiid == 0) std::printf("Microhh git-hash: " GITHASH "\n");

  // create the instances of the objects
  cinput  input (&mpi);
  cmodel  model (&mpi, &input);

  // read the input data
  if(input.readinifile())
    return 1;
  if(input.readproffile())
    return 1;

  if(mpi.readinifile(&input))
    return 1;
  if(model.readinifile())
    return 1;

  // init the mpi 
  if(mpi.init())
    return 1;
  if(model.init())
    return 1;

  if(mpi.mode == "init")
  {
    if(model.create())
      return 1;

    // save the data
    if(model.save())
      return 1;
  }
  else
  {
    // load the data
    if(model.load())
      return 1;
  }

  // check unused input
  input.printUnused(); 
  // free the memory of the input
  input.clear();

  // run the model
  if(mpi.mode != "init")
    if(model.exec())
      return 1;

  return 0;
}

