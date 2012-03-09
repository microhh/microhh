#include "grid.h"
#include "fields.h"

int main()
{
  // INIT
  // initialize the MPI interface
  // create the grid
  cgrid grid;

  // create the fields
  cfields fields(&grid);
  // END INIT

  // DNS
  // create the model
  
  // start the time loop
  // 1. boundary conditions
  // 2. advection
  // 3. diffusion
  // 4. gravity
  // 5. large scale forcings
  // 6. pressure

  // END DNS
  
  return 0;
}

