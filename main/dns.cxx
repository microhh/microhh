#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
//#include "mpiinterface2.h"
#include "boundary.h"
#include "advec.h"
#include "diff.h"
#include "force.h"
#include "buoyancy.h"
#include "pres.h"
#include "timeloop.h"

int main(int argc, char *argv[])
{
  // set the name of the simulation
  std::string simname("microhh");
  if(argc > 1)
    simname = argv[1];

  // create the instances of the objects
  cinput  input;
  cmpi    mpi;
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

  // free the memory of the input
  input.clear();

  // fill the fields with data
  if(grid.load())
    return 1;
  if(timeloop.load(timeloop.iteration))
    return 1;
  if(fields.load(timeloop.iteration))
    return 1;

  // initialize the diffusion to get the time step requirement
  if(diff.init())
    return 1;

  // initialize the pressure solver
  // CvH check this later, init is already using information from grid, should not happen...
  if(pres.init())
    return 1;
  if(pres.load())
    return 1;

  // initialize the check variables
  int    iter;
  double time, dt;
  double mom, tke, mass;
  double div;
  double cfl, dn;
  double cputime, start, end;

  // write output file header to the main processor and set the time
  FILE *dnsout = NULL;
  if(mpi.mpiid == 0)
  {
    std::string outputname = simname + ".out";
    dnsout = std::fopen(outputname.c_str(), "a");
    std::setvbuf(dnsout, NULL, _IOLBF, 1024);
    std::fprintf(dnsout, "%8s %11s %10s %11s %8s %8s %11s %16s %16s %16s\n",
      "ITER", "TIME", "CPUDT", "DT", "CFL", "DNUM", "DIV", "MOM", "TKE", "MASS");
  }

  // set the boundary conditions
  boundary.exec();

  // print the initial information
  if(timeloop.docheck() && !timeloop.insubstep())
  {
    iter    = timeloop.iteration;
    time    = timeloop.time;
    cputime = 0;
    dt      = timeloop.dt;
    cfl     = advec.getcfl(timeloop.dt);
    dn      = diff.getdn(timeloop.dt);
    div     = pres.check();
    mom     = fields.checkmom();
    tke     = fields.checktke();
    mass    = fields.checkmass();

    // write the output to file
    if(mpi.mpiid == 0)
      std::fprintf(dnsout, "%8d %11.3E %10.4f %11.3E %8.4f %8.4f %11.3E %16.8E %16.8E %16.8E\n",
        iter, time, cputime, dt, cfl, dn, div, mom, tke, mass);
    }

  // catch the start time for the first iteration
  start = mpi.gettime();

  // start the time loop
  while(true)
  {
    // determine the time step
    if(!timeloop.insubstep())
    {
      cfl = advec.getcfl(timeloop.dt);
      dn  = diff.getdn(timeloop.dt);
      timeloop.settimestep(cfl, dn);
    }

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
    if(timeloop.dosave() && !timeloop.insubstep())
      fields.p->save(timeloop.iteration);

    // exit the simulation when the runtime has been hit after the pressure calculation
    if(!timeloop.loop)
      break;

    // perform the timestepping substep
    timeloop.exec();

    // boundary conditions
    boundary.exec();

    if(!timeloop.insubstep())
      timeloop.timestep();

    if(timeloop.docheck() && !timeloop.insubstep())
    {
      iter    = timeloop.iteration;
      time    = timeloop.time;
      dt      = timeloop.dt;
      div     = pres.check();
      mom     = fields.checkmom();
      tke     = fields.checktke();
      mass    = fields.checkmass();

      end     = mpi.gettime();
      cputime = end - start;
      start   = end;

      // write the output to file
      if(mpi.mpiid == 0)
        std::fprintf(dnsout, "%8d %11.3E %10.4f %11.3E %8.4f %8.4f %11.3E %16.8E %16.8E %16.8E\n",
          iter, time, cputime, dt, cfl, dn, div, mom, tke, mass);
    }

    if(timeloop.dosave() && !timeloop.insubstep())
    {
      timeloop.save(timeloop.iteration);
      fields.save  (timeloop.iteration);
    }
  }

  // close the output file
  if(mpi.mpiid == 0)
    std::fclose(dnsout);
  
  return 0;
}

