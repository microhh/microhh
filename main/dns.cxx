#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "timeloop.h"
#include "advec.h"
#include "diff.h"
#include "force.h"
#include "pres.h"

int main()
{
  // create the instances of the objects
  cinput  input;
  cmpi    mpi;
  cgrid   grid  (&mpi);
  cfields fields(&grid, &mpi);

  // read the input data
  if(input.readinifile())
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
  fields.initfields();

  // create the instances of the model operations
  ctimeloop timeloop(&grid, &fields, &mpi);
  cadvec    advec   (&grid, &fields, &mpi);
  cdiff     diff    (&grid, &fields, &mpi);
  cpres     pres    (&grid, &fields, &mpi);
  cforce    force   (&grid, &fields, &mpi);

  // read the inputdata
  if(timeloop.readinifile(&input))
    return 1;

  // fill the fields with data
  if(grid.load())
    return 1;
  if(pres.load())
    return 1;
  if(timeloop.load(timeloop.iteration))
    return 1;
  if(fields.load(timeloop.iteration))
    return 1;

  // initialize the diffusion to get the time step requirement
  diff.init();

  // initialize the pressure solver
  pres.init();

  // initialize the check variables
  int    iter;
  double time, dt;
  double mom, tke, mass;
  double div;
  double cfl, dnum;
  double cputime, start, end;

  // write output file header to the main processor and set the time
  FILE *dnsout;
  if(mpi.mpiid == 0)
  {
    dnsout = std::fopen("dns.out", "w");
    std::setvbuf(dnsout, NULL, _IOLBF, 1024);
    std::fprintf(dnsout, "%8s %12s %10s %10s %8s %8s %13s %13s %13s %13s\n", 
      "ITER", "TIME", "CPUDT", "DT", "CFL", "DNUM", "DIV", "MOM", "TKE", "MASS");

    start = mpi.gettime();
  }

  // set the boundary conditions
  fields.boundary();

  // start the time loop
  while(timeloop.loop)
  {
    // determine the time step
    if(!timeloop.insubstep())
    {
      cfl  = advec.getcfl(timeloop.dt);
      dnum = diff.getdnum(timeloop.dt);
      timeloop.settimestep(cfl, dnum);
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

    // boundary conditions
    fields.boundary();

    if(!timeloop.insubstep())
      timeloop.timestep();

    if(timeloop.docheck() && !timeloop.insubstep())
    {
      iter    = timeloop.iteration;
      time    = timeloop.time;
      dt      = timeloop.dt;
      div     = pres.check();
      mom     = fields.check(0);
      tke     = fields.check(1);
      mass    = fields.check(2);

      end     = mpi.gettime();
      cputime = end - start;
      start   = end;

      // write the output to file
      if(mpi.mpiid == 0)
        std::fprintf(dnsout, "%8d %12.3f %10.4f %10.4f %8.4f %8.4f %13.5E %13.5E %13.5E %13.5E\n", 
          iter, time, cputime, dt, cfl, dnum, div, mom, tke, mass);
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

