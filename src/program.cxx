#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "dns.h"
#include "advec.h"
#include "diff.h"
#include "pres.h"
#include "timeint.h"

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
  cdns     dns    (&grid, &fields);
  cadvec   advec  (&grid, &fields);
  cdiff    diff   (&grid, &fields);
  cpres    pres   (&grid, &fields);
  ctimeint timeint(&grid, &fields);
  
  // start the time loop
  while(dns.loop)
  {
    // 0. determine the time step
    if(not timeint.insubstep())
      dns.settimestep(advec.getcfl(dns.dt));
    // 1. boundary conditions
    fields.boundary();
    // 2. advection
    advec.exec();
    // 3. diffusion
    diff.exec();
    // 4. gravity
    // 5. large scale forcings
    // 6. pressure
    pres.exec();
    // 7. perform the timestepping substep
    timeint.exec(dns.dt);

    if(not timeint.insubstep())
    {
      dns.timestep();
      fields.check();
    }
  }

  // END DNS
  
  return 0;
}

