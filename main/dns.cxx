#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "timeloop.h"
#include "advec.h"
#include "diff.h"
#include "force.h"
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
  // grid.creategrid();
  grid.load();
  // fields.createfields();
  fields.load(0);
  // END INIT

  // DNS
  // create the model and the operators
  ctimeloop timeloop(&grid, &fields);
  cadvec    advec   (&grid, &fields);
  cdiff     diff    (&grid, &fields);
  cpres     pres    (&grid, &fields);
  cforce    force   (&grid, &fields);
  ctimeint  timeint (&grid, &fields);

  // initialize the pressure solver
  pres.init();

  // start the time loop
  while(timeloop.loop)
  {
    // 0. determine the time step
    if(not timeint.insubstep())
      timeloop.settimestep(advec.getcfl(timeloop.dt));
    // 1. boundary conditions
    fields.boundary();

    if(not timeint.insubstep())
    {
      fields.check();
      pres.divergence();
    }
    // 2. advection
    advec.exec();
    // 3. diffusion
    diff.exec();
    // 4. gravity
    // 5. large scale forcings
    force.exec(timeint.subdt(timeloop.dt));
    // 6. pressure
    pres.exec(timeint.subdt(timeloop.dt));
    // 7. perform the timestepping substep
    timeint.exec(timeloop.dt);

    if(not timeint.insubstep())
      timeloop.timestep();
  }
  // END DNS
  
  return 0;
}

