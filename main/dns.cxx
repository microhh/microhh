#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "timeloop.h"
// #include "advec_g2i2.h"
#include "advec_g2i4.h"
#include "diff.h"
#include "force.h"
#include "pres.h"

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

  // create the instances of the model operations
  ctimeloop timeloop(&grid, &fields, &mpi);
  cadvec    advec   (&grid, &fields, &mpi);
  cdiff     diff    (&grid, &fields, &mpi);
  cpres     pres    (&grid, &fields, &mpi);
  cforce    force   (&grid, &fields, &mpi);

  // read the inputdata
  if(timeloop.readinifile(&input))
    return 1;
  if(force.readinifile(&input))
    return 1;

  // free the memory of the input
  input.clear();

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
  double cfl, dn;
  double cputime, start, end;

  // write output file header to the main processor and set the time
  FILE *dnsout = NULL;
  if(mpi.mpiid == 0)
  {
    std::string outputname = simname + ".out";
    dnsout = std::fopen(outputname.c_str(), "a");
    std::setvbuf(dnsout, NULL, _IOLBF, 1024);
    std::fprintf(dnsout, "%8s %12s %10s %10s %8s %8s %13s %13s %13s %13s\n", 
      "ITER", "TIME", "CPUDT", "DT", "CFL", "DNUM", "DIV", "MOM", "TKE", "MASS");
  }

  // set the boundary conditions
  fields.boundary();

  // catch the start time for the first iteration
  start = mpi.gettime();

  // start the time loop
  while(timeloop.loop)
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
      mom     = fields.checkmom();
      tke     = fields.checktke();
      mass    = fields.checkmass();

      end     = mpi.gettime();
      cputime = end - start;
      start   = end;

      // write the output to file
      if(mpi.mpiid == 0)
        std::fprintf(dnsout, "%8d %12.3f %10.4f %10.4f %8.4f %8.4f %13.5E %13.5E %13.5E %13.5E\n", 
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

