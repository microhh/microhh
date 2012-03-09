#include "grid.h"
#include "fields.h"

int main()
{
  // INIT
  // initialize the MPI interface
  // create the grid
  grid dnsgrid;

  // create the fields
  fields dnsfields(&dnsgrid);
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

