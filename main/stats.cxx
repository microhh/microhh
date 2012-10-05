#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "boundary.h"
#include "advec.h"
#include "diff.h"
#include "force.h"
#include "buoyancy.h"
#include "pres.h"
#include "timeloop.h"
#include "stats.h"

int main(int argc, char *argv[])
{
  // set the name of the simulation
  std::string simname("microhh");
  if(argc > 1)
    simname = argv[1];

  // start up the message passing interface
  cmpi mpi;
  mpi.startup();

  // create the instances of the objects
  cinput  input (&mpi);
  cgrid   grid  (&mpi);
  cfields fields(&grid, &mpi);

  // read the input data
  if(input.readinifile(simname))
    return 1;
  if(mpi.readinifile(&input))
    return 1;
  if(grid.readinifile(&input))
    return 1;
  if(fields.readinifile(&input))
    return 1;

  // init the mpi 
  if(mpi.init())
    return 1;
  // initialize the objects, allocate the required memory
  if(grid.init())
    return 1;
  if(fields.init())
    return 1;

  // create the boundary conditions class
  cboundary boundary(&grid, &fields, &mpi);

  // create the instances of the model operations
  ctimeloop timeloop(&grid, &fields, &mpi);
  cadvec    advec   (&grid, &fields, &mpi);
  cdiff     diff    (&grid, &fields, &mpi);
  cpres     pres    (&grid, &fields, &mpi);
  cforce    force   (&grid, &fields, &mpi);
  cbuoyancy buoyancy(&grid, &fields, &mpi);

  // read the inputdata
  if(boundary.readinifile(&input))
    return 1;
  if(advec.readinifile(&input))
    return 1;
  if(diff.readinifile(&input))
    return 1;
  if(force.readinifile(&input))
    return 1;
  if(buoyancy.readinifile(&input))
    return 1;
  if(pres.readinifile(&input))
    return 1;
  if(timeloop.readinifile(&input))
    return 1;

  // load the postprocessing modules
  cstats stats(&grid, &fields, &mpi);

  // free the memory of the input
  input.clear();

  // fill the fields with data
  if(grid.load())
    return 1;

  // initialize the pressure solver
  // CvH check this later, init is already using information from grid, should not happen...
  if(pres.init())
    return 1;
  if(pres.load())
    return 1;

  // initialize the diffusion to get the time step requirement
  if(diff.setvalues())
    return 1;

  // initialize the statistics
  if(stats.init())
    return 1;

  // initialize the check variables
  double cfl,dn;
  double div,mom,tke,mass;

  // start the time loop
  while(true)
  {
    // load the data
    if(timeloop.load(timeloop.iteration))
      return 1;
    if(fields.load(timeloop.iteration))
      return 1;

    // set the boundary conditions
    boundary.exec();

    // determine the time step
    timeloop.settimestep(cfl, dn);
    cfl = advec.getcfl(timeloop.dt);
    dn  = diff.getdn(timeloop.dt);

    // get diagnostic variables
    div  = pres.check();
    mom  = fields.checkmom();
    tke  = fields.checktke();
    mass = fields.checkmass();

    // advection
    advec.exec();
    // diffusion
    diff.exec();
    // large scale forcings
    force.exec(timeloop.getsubdt());
    // buoyancy
    buoyancy.exec();
    // pressure
    pres.exec(timeloop.getsubdt());

    if(timeloop.dostats() && !timeloop.insubstep())
      stats.exec(timeloop.iteration, timeloop.time);

    timeloop.postprocstep();

    // exit the simulation when the runtime has been hit after the pressure calculation
    if(!timeloop.loop)
      break;
  }

  return 0;
}

