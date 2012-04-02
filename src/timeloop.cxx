#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "timeloop.h"

ctimeloop::ctimeloop(cgrid *gridin, cfields *fieldsin)
{
  std::printf("Creating instance of object dns\n");
  grid   = gridin;
  fields = fieldsin;

  // input parameters
  adaptivestep = true;
  runtime      = 10000.;

  // initializations
  loop      = true;
  time      = 0.;
  dt        = 0.1;
  iteration = 0;
  cflmax    = 1.5;

  const int ifactor = 1000;

  itime    = (int)(ifactor * time);
  iruntime = (int)(ifactor * runtime);
  idt      = (int)(ifactor * dt);

  gettimeofday(&start, NULL);
}

ctimeloop::~ctimeloop()
{
  std::printf("Destroying instance of object dns\n");
}

int ctimeloop::timestep()
{
  time  += dt;
  itime += idt;

  iteration++;

  if(time >= runtime)
    loop = false;

  if(iteration % 100 == 0) 
  {
    gettimeofday(&end, NULL);

    double timeelapsed;

    timeelapsed = (double)(end.tv_sec-start.tv_sec) + (double)(end.tv_usec-start.tv_usec) * 1.e-6;

    std::printf("Iteration = %6d, runtime = %7.1f, cputime = %10.5f\n", iteration, time, timeelapsed);

    start = end;
  }

  if(iteration % 500 == 0) 
  {
    (*fields->u).save(iteration);
    (*fields->v).save(iteration);
    (*fields->w).save(iteration);
    (*fields->p).save(iteration);
    (*fields->s).save(iteration);
  }

  return 0;
}

int ctimeloop::settimestep(double cfl)
{
  if(adaptivestep)
    dt = dt * cflmax/cfl;

  std::printf("CFL = %f, dt = %f\n", cfl, dt);

  return 0;
}

