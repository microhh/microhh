#include <cstdio>
#include "grid.h"
#include "mpiinterface.h"
#include "mpicheck.h"

int main()
{
  // create the instances of the objects
  cmpi      mpi;
  cinput    input;
  cgrid     grid;
  cmpicheck mpicheck(&grid, &mpi);

  // read the input data
  if(input.readinifile())
    return 1;
  if(mpi.readinifile(&input))
    return 1;
  if(grid.readinifile(&input))
    return 1;

  // initialize the objects, allocate the required memory
  if(mpi.init())
    return 1;
  if(grid.init(mpi.npx, mpi.npy))
    return 1;

  mpi.initTypes(&grid);

  // CHECK the layout
  mpicheck.showLayout();

  // fill the fields with data
  if(grid.create())
    return 1;

  // fill the fields with the test data
  mpicheck.create();

  // trigger the boundary conditions
  mpicheck.boundary();

  mpicheck.showLine();
  return 0;
}

