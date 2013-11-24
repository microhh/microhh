#include <string>
#include <cstdio>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "model.h"
#include "defines.h"

// boundary schemes
#include "boundary_surface.h"
#include "boundary_user.h"

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

// thermo schemes
#include "thermo.h"
#include "thermo_dry.h"
#include "thermo_moist.h"

// stats schemes
#include "stats_dns.h"
#include "stats_les.h"

cmodel::cmodel(cmpi * mpiin, cinput * inputin)
{
  mpi   = mpiin;
  input = inputin;

  // create the grid class
  grid = new cgrid(mpi);

  // create the fields class
  fields = new cfields(grid, mpi);

  // create the instances of the model operations
  timeloop = new ctimeloop(grid, fields, mpi);
  force    = new cforce   (grid, fields, mpi);
  buffer   = new cbuffer  (grid, fields, mpi);

  // set null pointers for classes that will be initialized later
  boundary = NULL;
  advec    = NULL;
  diff     = NULL;
  pres     = NULL;

  // load the postprocessing moduls
  stats = NULL;
  cross = new ccross(grid, fields, mpi);
}

cmodel::~cmodel()
{
  // delete the components in reversed order
  delete cross;
  delete stats;
  delete buffer;
  delete thermo;
  delete force;
  delete pres;
  delete diff;
  delete advec;
  delete timeloop;

  delete boundary;
  delete fields;
  delete grid;
}

int cmodel::readinifile()
{
  // input parameters
  int n = 0;

  // grid
  if(grid->readinifile(input))
    return 1;

  // fields
  if(fields->readinifile(input))
    return 1;

  // first, get the switches for the schemes
  n += input->getItem(&swadvec   , "advec"   , "swadvec"   , "", grid->swspatialorder);
  n += input->getItem(&swdiff    , "diff"    , "swdiff"    , "", grid->swspatialorder);
  n += input->getItem(&swpres    , "pres"    , "swpres"    , "", grid->swspatialorder);
  n += input->getItem(&swboundary, "boundary", "swboundary", "", "default");
  n += input->getItem(&swstats   , "stats"   , "swstats"   , "", "0");
  n += input->getItem(&swthermo  , "thermo"  , "swthermo"  , "", "off");

  // if one or more arguments fails, then crash
  if(n > 0)
    return 1;

  // check the advection scheme
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
  if(advec->readinifile(input))
    return 1;

  // check the diffusion scheme
  if(swdiff == "0")
    diff = new cdiff    (grid, fields, mpi);
  else if(swdiff == "2")
    diff = new cdiff_g2 (grid, fields, mpi);
  else if(swdiff == "42")
    diff = new cdiff_g42(grid, fields, mpi);
  else if(swdiff == "4")
    diff = new cdiff_g4 (grid, fields, mpi);
  // TODO move to new model file later
  else if(swdiff == "22")
  {
    diff = new cdiff_les_g2(grid, fields, mpi);
    // the subgrid model requires a surface model because of the MO matching at first level
    if(swboundary != "surface")
    {
      std::printf("ERROR swdiff == \"22\" requires swboundary == \"surface\"\n");
      return 1;
    }
  }
  else
  {
    std::printf("ERROR \"%s\" is an illegal value for swdiff\n", swdiff.c_str());
    return 1;
  }
  if(diff->readinifile(input))
    return 1;

  // check the pressure scheme
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
  if(pres->readinifile(input))
    return 1;

  // model operations
  if(force->readinifile(input))
    return 1;
  if(timeloop->readinifile(input))
    return 1;

  if(swthermo== "moist")
    thermo = new cthermo_moist(grid, fields, mpi);
  else if(swthermo == "dry")
    thermo = new cthermo_dry(grid, fields, mpi);
  else if(swthermo == "off")
    thermo = new cthermo(grid, fields, mpi);
  else
  {
    std::printf("ERROR \"%s\" is an illegal value for swthermo\n", swthermo.c_str());
    return 1;
  }
  if(thermo->readinifile(input))
    return 1;

  // read the boundary and buffer in the end because they need to know the requested fields
  if(swboundary == "surface")
    boundary = new cboundary_surface(grid, fields, mpi);
  else if(swboundary == "user")
    boundary = new cboundary_user(grid, fields, mpi);
  else if(swboundary == "default")
    boundary = new cboundary(grid, fields, mpi);
  else
  {
    std::printf("ERROR \"%s\" is an illegal value for swboundary\n", swboundary.c_str());
    return 1;
  }
  if(boundary->readinifile(input))
    return 1;

  if(buffer->readinifile(input))
    return 1;

  // statistics
  if(swstats == "0")
    stats = new cstats    (grid, fields, mpi);
  else if(swstats == "dns")
    stats = new cstats_dns(grid, fields, mpi);
  else if(swstats == "les")
    stats = new cstats_les(grid, fields, mpi);
  else
  {
    std::printf("ERROR \"%s\" is an illegal value for swstats\n", swstats.c_str());
    return 1;
  }

  if(stats->readinifile(input))
    return 1;
  if(cross->readinifile(input))
    return 1;

  // set dependencies at the same hierarchy level
  if(swboundary == "surface")
    static_cast<cboundary_surface *>(boundary)->setdepends(thermo);

  return 0;
}

int cmodel::init()
{
  if(grid->init())
    return 1;
  if(fields->init())
    return 1;
  if(boundary->init())
    return 1;
  if(buffer->init())
    return 1;
  if(force->init())
    return 1;
  if(pres->init())
    return 1;

  if(stats->init(timeloop->ifactor))
    return 1;
  if(cross->init(timeloop->ifactor))
    return 1;

  return 0;
}

int cmodel::load()
{
  if(grid->load())
    return 1;
  if(timeloop->load(timeloop->iotime))
    return 1;
  if(fields->load(timeloop->iotime))
    return 1;
  if(boundary->load(timeloop->iotime))
    return 1;
  if(buffer->create(input))
    return 1;
  if(force->create(input))
    return 1;
  if(thermo->create())
    return 1;

  if(swstats == "les")
  {
    cstats_les *stats_lesptr = static_cast<cstats_les *>(stats);
    if(stats_lesptr->create(timeloop->iotime, thermo))
      return 1;
  }
  else
  {
    if(stats->create(timeloop->iotime))
      return 1;
  }


  if(boundary->setvalues())
    return 1;
  if(diff->setvalues())
    return 1;
  if(pres->setvalues())
    return 1;

  return 0;
}

int cmodel::create()
{
  if(grid->create(input))
    return 1;
  if(fields->create(input))
    return 1;

  return 0;
}

int cmodel::save()
{
  if(grid->save())
    return 1;
  if(fields->save(timeloop->iotime))
    return 1;
  if(timeloop->save(timeloop->iotime))
    return 1;
  if(boundary->save(timeloop->iotime))
    return 1;

  return 0;
}

int cmodel::exec()
{

  // set the boundary conditions
  boundary->exec();
  diff->execvisc(boundary);

  if(settimestep())
    return 1;

  // print the initial information
  if(outputfile(!timeloop->loop))
    return 1;

  // start the time loop
  while(true)
  {
    // determine the time step
    if(!timeloop->insubstep())
    {
      if(settimestep())
        return 1;
    }
    // advection
    advec->exec();
    // diffusion
    diff->exec();
    // thermo
    thermo->exec();
    // buffer
    buffer->exec();
    // large scale forcings
    force->exec(timeloop->getsubdt());

    // pressure
    pres->exec(timeloop->getsubdt());

    // statistics when not in substep and not directly after restart
    if(!timeloop->insubstep() && !((timeloop->iteration > 0) && (timeloop->itime == timeloop->istarttime)))
    {
      if(swstats == "les")
      {
        cstats_les *stats_lesptr = static_cast<cstats_les *>(stats);
        stats_lesptr->exec(timeloop->iteration, timeloop->time, timeloop->itime, thermo);
      }
      else
        stats->exec(timeloop->iteration, timeloop->time, timeloop->itime);

      cross->exec(timeloop->time, timeloop->itime, timeloop->iotime);
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

      // save the data for a restart
      if(timeloop->dosave() && !timeloop->insubstep())
      {
        // save the time data
        timeloop->save(timeloop->iotime);
        // save the fields
        fields->save  (timeloop->iotime);
        // save the boundary data
        boundary->save(timeloop->iotime);
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
      if(timeloop->load(timeloop->iotime))
        return 1;
      if(fields->load(timeloop->iotime))
        return 1;
      if(boundary->load(timeloop->iotime))
        return 1;
    }

    // boundary conditions
    boundary->exec();
    diff->execvisc(boundary);

    if(outputfile(!timeloop->loop))
      return 1;

  }

  return 0;
}

int cmodel::outputfile(bool doclose)
{
  // initialize the check variables
  int    nerror=0;
  int    iter;
  double time, dt;
  double mom, tke, mass;
  double div;
  double cfl, dn;
  double cputime, end;
  static double start;
  static FILE *dnsout = NULL;

  // write output file header to the main processor and set the time
  if(mpi->mpiid == 0 && dnsout == NULL)
  {
    std::string outputname = mpi->simname + ".out";
    dnsout = std::fopen(outputname.c_str(), "a");
    std::setvbuf(dnsout, NULL, _IOLBF, 1024);
    std::fprintf(dnsout, "%8s %11s %10s %11s %8s %8s %11s %16s %16s %16s\n",
      "ITER", "TIME", "CPUDT", "DT", "CFL", "DNUM", "DIV", "MOM", "TKE", "MASS");
    start = mpi->gettime();
  }

  if(timeloop->docheck() && !timeloop->insubstep())
  {
    iter    = timeloop->iteration;
    time    = timeloop->time;
    dt      = timeloop->dt;
    div     = pres->check();
    mom     = fields->checkmom();
    tke     = fields->checktke();
    mass    = fields->checkmass();
    cfl     = advec->getcfl(timeloop->dt);
    dn      = diff->getdn(timeloop->dt);

    end     = mpi->gettime();
    cputime = end - start;
    start   = end;

    // write the output to file
    if(mpi->mpiid == 0)
    {
      std::fprintf(dnsout, "%8d %11.3E %10.4f %11.3E %8.4f %8.4f %11.3E %16.8E %16.8E %16.8E\n",
        iter, time, cputime, dt, cfl, dn, div, mom, tke, mass);
    }
  }

  if(doclose)
  {
    // close the output file
    if(mpi->mpiid == 0)
    std::fclose(dnsout);
  }

  return(nerror>0);
}

int cmodel::settimestep()
{
  if(timeloop->settimelim())
    return 1;

  timeloop->idtlim = std::min(timeloop->idtlim, advec->gettimelim(timeloop->idt, timeloop->dt));
  timeloop->idtlim = std::min(timeloop->idtlim, diff ->gettimelim(timeloop->idt, timeloop->dt));
  timeloop->idtlim = std::min(timeloop->idtlim, stats->gettimelim(timeloop->itime));
  timeloop->idtlim = std::min(timeloop->idtlim, cross->gettimelim(timeloop->itime));
  timeloop->settimestep();

  return 0;
}
