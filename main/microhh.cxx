#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "model.h"

int main(int argc, char *argv[])
{
  // start up the message passing interface
  cmpi mpi;
  mpi.startup();

  // set the mode and the name of the simulation
  std::string mode;
  std::string simname("microhh");
  if(argc <= 1)
  {
    if(mpi.mpiid == 0) std::printf("ERROR: specify init, run or post mode\n");
    return 1;
  }
  else
  {
    // check the execution mode
    mode = argv[1];
    if(mode != "init" && mode != "run" && mode != "post")
    {
      if(mpi.mpiid == 0) std::printf("ERROR: specify init, run or post mode\n");
      return 1;
    }
    // set the name of the simulation
    if(argc > 2)
      simname = argv[2];
  }

  // create the instances of the objects
  cinput  input (&mpi);
  cgrid   grid  (&mpi);
  cfields fields(&grid, &mpi);

  cmodel  model (&grid, &fields, &mpi, simname);

  // read the input data
  if(input.readinifile(simname))
    return 1;
  if(mode == "init")
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

  if(mode == "init")
  {
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
  }
  else
  {
  // fill the fields with data
  if(grid.load())
    return 1;
  // CvH, the model loads the fields now, not good...
  if(model.load())
    return 1;
  }

  // free the memory of the input
  input.clear();

  if(mode != "init")
  {
    // run the model
    model.exec(mode);
  }

  return 0;
}
