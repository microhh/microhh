#include <string>
#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "model.h"
#include "defines.h"

// advection schemes
#include "advec_g2.h"
#include "advec_g2i4.h"
#include "advec_g42.h"
#include "advec_g4.h"
#include "advec_g4m.h"

// diffusion schemes
#include "diff_g2.h"
#include "diff_g42.h"
#include "diff_g4.h"
#include "diff_les_g2.h"

// pressure schemes
#include "pres_g2.h"
#include "pres_g42.h"
#include "pres_g4.h"

cmodel::cmodel(cgrid *gridin, cmpi *mpiin)
{
  grid    = gridin;
  mpi     = mpiin;

  // create the fields class
  fields   = new cfields  (grid, mpi);

  // create the boundary conditions class
  boundary = new cboundary(grid, fields, mpi);

  // create the instances of the model operations
  timeloop = new ctimeloop(grid, fields, mpi);
  force    = new cforce   (grid, fields, mpi);
  buoyancy = new cbuoyancy(grid, fields, mpi);
  buffer   = new cbuffer  (grid, fields, mpi);

  // load the postprocessing moduls
  stats    = new cstats   (grid, fields, mpi);
  cross    = new ccross   (grid, fields, mpi);

  // set null pointers for classes that will be initialized later
  advec = NULL;
  diff  = NULL;
  pres  = NULL;
}

cmodel::~cmodel()
{
  // delete the components in reversed order
  delete cross;
  delete stats;
  delete buffer;
  delete buoyancy;
  delete force;
  delete pres;
  delete diff;
  delete advec;
  delete timeloop;

  delete boundary;
  delete fields;
}

int cmodel::readinifile(cinput *inputin)
{
  // input parameters
  int n = 0;

  // fields
  if(fields->readinifile(inputin))
    return 1;

  // check the advection scheme
  n += inputin->getItem(&swadvec, "advec", "swadvec", "", grid->swspatialorder);
  if(swadvec == "0")
    advec = new cadvec     (grid, fields, mpi);
  else if(swadvec == "2")
    advec = new cadvec_g2  (grid, fields, mpi);
  else if(swadvec == "24")
    advec = new cadvec_g2i4(grid, fields, mpi);
  else if(swadvec == "42")
    advec = new cadvec_g42 (grid, fields, mpi);
  else if(swadvec == "4")
    advec = new cadvec_g4  (grid, fields, mpi);
  else if(swadvec == "44")
    advec = new cadvec_g4m (grid, fields, mpi);
  else
  {
    std::printf("ERROR \"%s\" is an illegal value for swadvec\n", swadvec.c_str());
    return 1;
  }
  if(advec->readinifile(inputin))
    return 1;

  // check the diffusion scheme
  n += inputin->getItem(&swdiff, "diff", "swdiff", "", grid->swspatialorder);
  if(swdiff == "0")
    diff = new cdiff    (grid, fields, mpi);
  else if(swdiff == "2")
    diff = new cdiff_g2 (grid, fields, mpi);
  else if(swdiff == "42")
    diff = new cdiff_g42(grid, fields, mpi);
  else if(swdiff == "4")
    diff = new cdiff_g4 (grid, fields, mpi);
  // CvH move to new model file later
  else if(swdiff == "22")
    diff = new cdiff_les_g2(grid, fields, mpi);
  else
  {
    std::printf("ERROR \"%s\" is an illegal value for swdiff\n", swdiff.c_str());
    return 1;
  }
  if(diff->readinifile(inputin))
    return 1;


  // check the pressure scheme
  n += inputin->getItem(&swpres, "pres", "swpres", "", grid->swspatialorder);
  if(swpres == "0")
    pres = new cpres    (grid, fields, mpi);
  else if(swpres == "2")
    pres = new cpres_g2 (grid, fields, mpi);
  else if(swpres == "42")
    pres = new cpres_g42(grid, fields, mpi);
  else if(swpres == "4")
    pres = new cpres_g4 (grid, fields, mpi);
  else
  {
    std::printf("ERROR \"%s\" is an illegal value for swpres\n", swpres.c_str());
    return 1;
  }
  if(pres->readinifile(inputin))
    return 1;

  // model operations
  if(force->readinifile(inputin))
    return 1;
  if(buoyancy->readinifile(inputin))
    return 1;
  if(timeloop->readinifile(inputin))
    return 1;

  // model classes that need to know prognostic fields
  if(boundary->readinifile(inputin))
    return 1;
  if(buffer->readinifile(inputin))
    return 1;

  // statistics
  if(stats->readinifile(inputin))
    return 1;
  if(cross->readinifile(inputin))
    return 1;

  // if one or more arguments fails, then crash
  if(n > 0)
    return 1;

  return 0;
}

int cmodel::init()
{
  if(fields->init())
    return 1;
  // TODO change this to surface init
  if(boundary->init())
    return 1;
  if(buffer->init())
    return 1;
  if(pres->init())
    return 1;
  if(stats->init())
    return 1;

  return 0;
}

int cmodel::load()
{
  if(timeloop->load(timeloop->iteration))
    return 1;
  if(fields->load(timeloop->iteration))
    return 1;
  if(buffer->load())
    return 1;
  if(stats->create(timeloop->iteration))
    return 1;

  // initialize the diffusion to get the time step requirement
  if(boundary->setvalues())
    return 1;
  if(diff->setvalues())
    return 1;
  if(pres->setvalues())
    return 1;

  return 0;
}

int cmodel::create(cinput *inputin)
{
  if(fields->create(inputin))
    return 1;
  if(buffer->setbuffers())
    return 1;

  return 0;
}

int cmodel::save()
{
  if(fields->save(timeloop->iteration))
    return 1;
  if(buffer->save())
    return 1;
  if(timeloop->save(timeloop->iteration))
    return 1;

  return 0;
}

int cmodel::exec()
{
  // initialize the check variables
  int    iter;
  double time, dt;
  double mom, tke, mass;
  double div;
  double cfl, dn;
  double cputime, start, end;

  // write output file header to the main processor and set the time
  FILE *dnsout = NULL;
  if(mpi->mpiid == 0)
  {
    std::string outputname = mpi->simname + ".out";
    dnsout = std::fopen(outputname.c_str(), "a");
    std::setvbuf(dnsout, NULL, _IOLBF, 1024);
    std::fprintf(dnsout, "%8s %11s %10s %11s %8s %8s %11s %16s %16s %16s\n",
      "ITER", "TIME", "CPUDT", "DT", "CFL", "DNUM", "DIV", "MOM", "TKE", "MASS");
  }

  // set the boundary conditions
  boundary->exec();

  // TODO very ugly!!!
  if(swdiff == "22")
    (*(cdiff_les_g2 *)diff).execvisc(boundary);

  // set the initial cfl and dn
  cfl = advec->getcfl(timeloop->dt);
  dn  = diff->getdn(timeloop->dt);

  // print the initial information
  if(timeloop->docheck() && !timeloop->insubstep())
  {
    iter    = timeloop->iteration;
    time    = timeloop->time;
    cputime = 0;
    dt      = timeloop->dt;
    div     = pres->check();
    mom     = fields->checkmom();
    tke     = fields->checktke();
    mass    = fields->checkmass();

    // write the output to file
    if(mpi->mpiid == 0)
      std::fprintf(dnsout, "%8d %11.3E %10.4f %11.3E %8.4f %8.4f %11.3E %16.8E %16.8E %16.8E\n",
        iter, time, cputime, dt, cfl, dn, div, mom, tke, mass);
    }

  // catch the start time for the first iteration
  start = mpi->gettime();

  // start the time loop
  while(true)
  {
    // determine the time step
    if(!timeloop->insubstep())
    {
      cfl = advec->getcfl(timeloop->dt);
      dn  = diff->getdn(timeloop->dt);
      timeloop->settimestep(cfl, dn);
    }

    // advection
    advec->exec();
    // diffusion
    diff->exec();
    // large scale forcings
    force->exec(timeloop->getsubdt());
    // buoyancy
    buoyancy->exec();
    // buffer
    buffer->exec();

    // pressure
    pres->exec(timeloop->getsubdt());
    if(timeloop->dosave() && !timeloop->insubstep())
      fields->s["p"]->save(timeloop->iteration, fields->s["tmp1"]->data, fields->s["tmp2"]->data);

    if(timeloop->dostats() && !timeloop->insubstep())
    {
      stats->exec(timeloop->iteration, timeloop->time);
      cross->exec(timeloop->iteration);
    }

    // exit the simulation when the runtime has been hit after the pressure calculation
    if(!timeloop->loop)
      break;

    // RUN MODE
    if(mpi->mode == "run")
    {
      // integrate in time
      timeloop->exec();

      // step the time step
      if(!timeloop->insubstep())
        timeloop->timestep();

      // save the fields
      if(timeloop->dosave() && !timeloop->insubstep())
      {
        timeloop->save(timeloop->iteration);
        fields->save  (timeloop->iteration);
      }
    }

    // POST PROCESS MODE
    else if(mpi->mode == "post")
    {
      // step to the next time step
      timeloop->postprocstep();

      // if simulation is done break
      if(!timeloop->loop)
        break;

      // load the data
      if(timeloop->load(timeloop->iteration))
        return 1;
      if(fields->load(timeloop->iteration))
        return 1;
    }

    // boundary conditions
    // TODO UGLY!!
    if(swdiff == "22")
      (*(cdiff_les_g2 *)diff).execvisc(boundary);

    boundary->exec();

    if(timeloop->docheck() && !timeloop->insubstep())
    {
      iter    = timeloop->iteration;
      time    = timeloop->time;
      dt      = timeloop->dt;
      div     = pres->check();
      mom     = fields->checkmom();
      tke     = fields->checktke();
      mass    = fields->checkmass();

      end     = mpi->gettime();
      cputime = end - start;
      start   = end;

      // write the output to file
      if(mpi->mpiid == 0)
        std::fprintf(dnsout, "%8d %11.3E %10.4f %11.3E %8.4f %8.4f %11.3E %16.8E %16.8E %16.8E\n",
          iter, time, cputime, dt, cfl, dn, div, mom, tke, mass);
    }
  }

  // close the output file
  if(mpi->mpiid == 0)
    std::fclose(dnsout);
  
  return 0;
}
