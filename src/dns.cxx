#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "dns.h"

cdns::cdns(cgrid *gridin, cfields *fieldsin)
{
  std::printf("Creating instance of object dns\n");
  grid   = gridin;
  fields = fieldsin;

  loop = true;
  adaptivestep = true;

  time      = 0.;
  runtime   = 20000.;
  dt        = 0.1;
  iteration = 0;
  cflmax    = 0.5;

  const int ifactor = 1000;

  itime    = (int)(ifactor * time);
  iruntime = (int)(ifactor * runtime);
  idt      = (int)(ifactor * dt);

  clock_gettime(CLOCK_REALTIME, &start);
}

cdns::~cdns()
{
  std::printf("Destroying instance of object dns\n");
}

int cdns::timestep()
{
  time  += dt;
  itime += idt;

  iteration++;

  if(time >= runtime)
    loop = false;

  if(iteration % 100 == 0) 
  {
    clock_gettime(CLOCK_REALTIME, &end);

    double timeelapsed;

    timeelapsed = (double)(end.tv_sec-start.tv_sec) + (double)(end.tv_nsec-start.tv_nsec) * 1.e-9;

    std::printf("Iteration = %6d, runtime = %7.1f, cputime = %10.7f\n", iteration, time, timeelapsed);

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

int cdns::settimestep(double cfl)
{
  if(adaptivestep)
    dt = dt * cflmax/cfl;

  std::printf("CFL = %f, dt = %f\n", cfl, dt);

  return 0;
}

