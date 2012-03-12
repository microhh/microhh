#include "grid.h"
#include "fields.h"
#include "dns.h"
#include "advec.h"
#include "diff.h"

int main()
{
  // INIT
  // initialize the MPI interface
  // create the objects, read the inputdata
  cgrid grid;
  cfields fields(&grid);

  // initialize the objects, allocate the required memory
  grid.initgrid();
  fields.initfields();

  // create the objects, fill the fields with data
  grid.creategrid();
  fields.createfields();
  // END INIT

  // DNS
  // create the model and the operators
  cdns   dns  (&grid, &fields);
  cadvec advec(&grid, &fields);
  cdiff  diff (&grid, &fields);
  
  // start the time loop
  while(dns.loop)
  {
    // 1. boundary conditions
    fields.boundary_bottop();
    // 2. advection
    advec.exec(dns.dt);
    // 3. diffusion
    diff.exec();
    // 4. gravity
    // 5. large scale forcings
    // 6. pressure
    dns.timestep();
    fields.resettend();
    fields.check();
  }

  // END DNS
  
  return 0;
}

