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
  adaptivestep = false;

  time      = 0.;
  runtime   = 1.;
  dt        = 0.01;
  iteration = 0;
  cflmax    = 0.8;

  const int ifactor = 1000;

  itime    = (int)(ifactor * time);
  iruntime = (int)(ifactor * runtime);
  idt      = (int)(ifactor * dt);
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
    std::printf("Iteration = %6d, time = %7.1f\n", iteration, time);

  return 0;
}

int cdns::settimestep(double cfl)
{
  if(adaptivestep)
    dt = dt * cflmax/cfl;

  std::printf("CFL = %f, dt = %f\n", cfl, dt);

  return 0;
}

