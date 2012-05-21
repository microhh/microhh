#include <cstdio>
#include "grid.h"
#include "mpiinterface.h"
#include "mpicheck.h"

int main(int argc, char *argv[])
{
  // set the name of the simulation
  std::string simname("microhh");
  if(argc > 1)
    simname = argv[1];

  // create the instances of the objects
  cinput    input;
  cmpi      mpi;
  cgrid     grid(&mpi);
  cmpicheck mpicheck(&grid, &mpi);

  if(input.readinifile(simname))
    return 1;
  if(grid.readinifile(&input))
    return 1;
  if(mpi.readinifile(&input))
    return 1;

  // initialize the objects, allocate the required memory
  if(mpi.init())
    return 1;
  if(grid.init())
    return 1;

  // check the layout
  mpicheck.checkLayout();

  // fill the fields with data
  if(grid.create())
    return 1;

  // fill the fields with the test data
  mpicheck.create();

  // trigger the boundary conditions
  mpicheck.checkBoundary();

  mpicheck.checkTranspose();

  return 0;
}

