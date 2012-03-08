#include <iostream>
#include "grid.h"
#include "fields.h"

int main()
{
  // initialize the MPI interface
  // create the grid
  grid dnsgrid;

  // create the fields
  fields dnsfields(&dnsgrid);

  // create the model
  
  // start the time loop
  // 1. boundary conditions
  // 2. advection
  // 3. diffusion
  // 4. gravity
  // 5. large scale forcings
  // 6. pressure
  std::cout << "Hello world!" << std::endl;
  return 0;
}

