#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "timeloop.h"
#include "advec.h"
#include "diff.h"
#include "force.h"
#include "pres.h"

int main()
{
  // read the input data
  cinput input;
  if(input.readinifile())
    return 1;
  
  // create the objects, read the inputdata
  cgrid grid;
  if(grid.readinifile(&input))
    return 1;

  cfields fields(&grid);

  // initialize the objects, allocate the required memory
  grid.initgrid();
  fields.initfields();

  // create the model and the operators
  ctimeloop timeloop(&grid, &fields);
  cadvec    advec   (&grid, &fields);
  cdiff     diff    (&grid, &fields);
  cpres     pres    (&grid, &fields);
  cforce    force   (&grid, &fields);

  // read the inputdata
  if(timeloop.readinifile(&input))
    return 1;

  // fill the fields with data
  if(grid.load())
    return 1;
  if(fields.load(timeloop.iteration))
    return 1;
  if(pres.load())
    return 1;

  // initialize the diffusion to get the time step requirement
  diff.init();

  // initialize the pressure solver
  pres.init();

  // start the time loop
  while(timeloop.loop)
  {
    // determine the time step
    if(!timeloop.insubstep())
    {
      double cfl = advec.getcfl(timeloop.dt);
      timeloop.settimestep(cfl);
    }

    // boundary conditions
    fields.boundary();

    if(!timeloop.insubstep())
    {
      fields.check();
      pres.divergence();
    }
    // advection
    advec.exec();
    // diffusion
    diff.exec();
    // large scale forcings
    force.exec(timeloop.getsubdt());
    // pressure
    pres.exec(timeloop.getsubdt());
    // perform the timestepping substep
    timeloop.exec();

    if(!timeloop.insubstep())
      timeloop.timestep();
  }
  
  return 0;
}

