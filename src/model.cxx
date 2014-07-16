/*
 * MicroHH
 * Copyright (c) 2011-2013 Chiel van Heerwaarden
 * Copyright (c) 2011-2013 Thijs Heus
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <string>
#include <cstdio>
#include <algorithm>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "model.h"
#include "defines.h"
#include "timeloop.h"
#include "buffer.h"
#include "force.h"
#include "stats.h"
#include "cross.h"
#include "budget.h"

// boundary schemes
#include "boundary.h"
#include "boundary_surface.h"
#include "boundary_user.h"

// advection schemes
#include "advec.h"
#include "advec_2.h"
#include "advec_2int4.h"
#include "advec_4.h"
#include "advec_4m.h"

// diffusion schemes
#include "diff.h"
#include "diff_2.h"
#include "diff_4.h"
#include "diff_les2s.h"

// pressure schemes
#include "pres.h"
#include "pres_2.h"
#include "pres_4.h"

// thermo schemes
#include "thermo.h"
#include "thermo_buoy.h"
#include "thermo_buoy_slope.h"
#include "thermo_dry.h"
#include "thermo_moist.h"

cmodel::cmodel(cmaster *masterin, cinput *inputin)
{
  master = masterin;
  input  = inputin;

  // create the grid class
  grid = new cgrid(this);

  // create the fields class
  fields = new cfields(this);

  // create the instances of the model operations
  timeloop = new ctimeloop(this);
  force    = new cforce(this);
  buffer   = new cbuffer(this);

  // set null pointers for classes that will be initialized later
  boundary = NULL;
  advec    = NULL;
  diff     = NULL;
  pres     = NULL;
  thermo   = NULL;

  // load the postprocessing modules
  stats  = new cstats(this);
  cross  = new ccross(this);
  budget = new cbudget(this);
}

cmodel::~cmodel()
{
  // delete the components in reversed order
  delete budget;
  delete cross;
  delete stats;
  delete buffer;
  delete force;
  delete pres;
  delete diff;
  delete advec;
  delete timeloop;
  delete thermo;

  delete boundary;
  delete fields;
  delete grid;
}

int cmodel::readinifile()
{
  // input parameters
  int nerror = 0;

  // grid
  if(grid->readinifile(input))
    return 1;

  // fields
  if(fields->readinifile(input))
    return 1;

  // first, get the switches for the schemes
  nerror += input->getItem(&swadvec     , "advec"   , "swadvec"   , "", grid->swspatialorder);
  nerror += input->getItem(&swdiff      , "diff"    , "swdiff"    , "", grid->swspatialorder);
  nerror += input->getItem(&swpres      , "pres"    , "swpres"    , "", grid->swspatialorder);
  nerror += input->getItem(&swboundary  , "boundary", "swboundary", "", "default");
  nerror += input->getItem(&swthermo    , "thermo"  , "swthermo"  , "", "0");

  // Get base state option (boussinesq or anelastic)
  nerror += input->getItem(&swbasestate , "grid"  , "swbasestate" , "", "");  // BvS: where to put switch??

  if(!(swbasestate == "boussinesq" || swbasestate == "anelastic"))
  {
    if(master->mpiid == 0) std::printf("ERROR \"%s\" is an illegal value for swbasestate\n", swbasestate.c_str());
    return 1;
  }

  // 2nd order with LES diffusion is only option supporting anelastic (for now) 
  //if(swdiff != "les2s" and swbasestate != "boussinesq")
  //{
  //  std::printf("ERROR swdiff=%s is not allowed with swbasestate=%s \n", swdiff.c_str(),swbasestate.c_str());
  //  return 1;
  //}
 
  // get the list of masks
  nerror += input->getList(&masklist, "stats", "masklist", "");
  for(std::vector<std::string>::const_iterator it=masklist.begin(); it!=masklist.end(); ++it)
  {
    if(*it != "wplus" &&
       *it != "wmin"  &&
       *it != "ql"    &&
       *it != "qlcore")
      if(master->mpiid == 0) std::printf("WARNING %s is an undefined mask for conditional statistics\n", it->c_str());
    else
      stats->addmask(*it);
  }



  // if one or more arguments fails, then crash
  if(nerror > 0)
    return 1;

  // check the advection scheme
  if(swadvec == "0")
    advec = new cadvec(this);
  else if(swadvec == "2")
    advec = new cadvec_2(this);
  else if(swadvec == "2int4")
    advec = new cadvec_2int4(this);
  else if(swadvec == "4")
    advec = new cadvec_4(this);
  else if(swadvec == "4m")
    advec = new cadvec_4m(this);
  else
  {
    if(master->mpiid == 0) std::printf("ERROR \"%s\" is an illegal value for swadvec\n", swadvec.c_str());
    return 1;
  }
  if(advec->readinifile(input))
    return 1;

  // check the diffusion scheme
  if(swdiff == "0")
    diff = new cdiff(this);
  else if(swdiff == "2")
    diff = new cdiff_2(this);
  else if(swdiff == "4")
    diff = new cdiff_4(this);
  // TODO move to new model file later?
  else if(swdiff == "les2s")
  {
    diff = new cdiff_les2s(this);
    // the subgrid model requires a surface model because of the MO matching at first level
    if(swboundary != "surface")
    {
      if(master->mpiid == 0) std::printf("ERROR swdiff == \"les2s\" requires swboundary == \"surface\"\n");
      return 1;
    }
  }
  else
  {
    if(master->mpiid == 0) std::printf("ERROR \"%s\" is an illegal value for swdiff\n", swdiff.c_str());
    return 1;
  }
  if(diff->readinifile(input))
    return 1;

  // check the pressure scheme
  if(swpres == "0")
    pres = new cpres(this);
  else if(swpres == "2")
    pres = new cpres_2(this);
  else if(swpres == "4")
    pres = new cpres_4(this);
  else
  {
    if(master->mpiid == 0) std::printf("ERROR \"%s\" is an illegal value for swpres\n", swpres.c_str());
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
    thermo = new cthermo_moist(this);
  else if(swthermo == "buoy")
    thermo = new cthermo_buoy(this);
  else if(swthermo == "dry")
    thermo = new cthermo_dry(this);
  else if(swthermo == "buoy_slope")
    thermo = new cthermo_buoy_slope(this);
  else if(swthermo == "0")
    thermo = new cthermo(this);
  else
  {
    if(master->mpiid == 0) std::printf("ERROR \"%s\" is an illegal value for swthermo\n", swthermo.c_str());
    return 1;
  }
  if(thermo->readinifile(input))
    return 1;

  // read the boundary and buffer in the end because they need to know the requested fields
  if(swboundary == "surface")
    boundary = new cboundary_surface(this);
  else if(swboundary == "user")
    boundary = new cboundary_user(this);
  else if(swboundary == "default")
    boundary = new cboundary(this);
  else
  {
    if(master->mpiid == 0) std::printf("ERROR \"%s\" is an illegal value for swboundary\n", swboundary.c_str());
    return 1;
  }
  if(boundary->readinifile(input))
    return 1;
  if(buffer->readinifile(input))
    return 1;

  if(stats->readinifile(input))
    return 1;
  if(cross->readinifile(input))
    return 1;
  if(budget->readinifile(input))
    return 1;

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
  if(thermo->init())
    return 1;
  if(buffer->init())
    return 1;
  if(force->init())
    return 1;
  if(pres->init())
    return 1;
  if(thermo->init())
    return 1;

  if(stats->init(timeloop->ifactor))
    return 1;
  if(cross->init(timeloop->ifactor))
    return 1;
  if(budget->init())
    return 1;

  return 0;
}

int cmodel::load()
{
  // first load the grid and time to make their information available
  if(grid->load())
    return 1;
  if(timeloop->load(timeloop->iotime))
    return 1;
  // initialize the statistics file to open the possiblity to add profiles
  if(stats->create(timeloop->iotime))
    return 1;

  if(fields->load(timeloop->iotime))
    return 1;

  // \TODO call boundary load for the data and then timedep, not nice...
  if(boundary->load(timeloop->iotime))
    return 1;
  if(boundary->create(input))
    return 1;

  if(buffer->create(input))
    return 1;
  if(force->create(input))
    return 1;
  if(thermo->create(input))
    return 1;

  if(budget->create())
    return 1;

  // end with modules that require all fields to be present
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
  if(master->mpiid == 0) std::printf("Starting time integration\n");
  // update the time dependent values
  boundary->settimedep();
  force->settimedep();

  // set the boundary conditions
  boundary->exec();

  // get the field means, in case needed
  fields->exec();
  // get the viscosity to be used in diffusion
  diff->execvisc();

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
      if(stats->dostats())
      {
        // always process the default mask
        stats->getmask(fields->sd["tmp3"], fields->sd["tmp4"], &stats->masks["default"]);
        calcstats("default");

        // work through the potential masks for the statistics
        for(std::vector<std::string>::const_iterator it=masklist.begin(); it!=masklist.end(); ++it)
        {
          if(*it == "wplus" || *it == "wmin")
          {
            fields->getmask(fields->sd["tmp3"], fields->sd["tmp4"], &stats->masks[*it]);
            calcstats(*it);
          }
          else if(*it == "ql" || *it == "qlcore")
          {
            thermo->getmask(fields->sd["tmp3"], fields->sd["tmp4"], &stats->masks[*it]);
            calcstats(*it);
          }
        }

        // store the stats data
        stats->exec(timeloop->iteration, timeloop->time, timeloop->itime);
      }

      if(cross->docross())
      {
        if(fields->execcross())
          return 1;
        if(thermo->execcross())
          return 1;
        if(boundary->execcross())
          return 1;
      }
    }

    // exit the simulation when the runtime has been hit after the pressure calculation
    if(!timeloop->loop)
      break;

    // RUN MODE
    if(master->mode == "run")
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
        fields->save(timeloop->iotime);
        // save the boundary data
        boundary->save(timeloop->iotime);
      }
    }

    // POST PROCESS MODE
    else if(master->mode == "post")
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
    // update the time dependent values
    boundary->settimedep();
    force->settimedep();

    // set the boundary conditions
    boundary->exec();
    // get the field means, in case needed
    fields->exec();
    // get the viscosity to be used in diffusion
    diff->execvisc();

    if(outputfile(!timeloop->loop))
      return 1;

  }

  return 0;
}

int cmodel::calcstats(std::string maskname)
{
  fields->execstats(&stats->masks[maskname]);
  thermo->execstats(&stats->masks[maskname]);
  budget->execstats(&stats->masks[maskname]);
  boundary->execstats(&stats->masks[maskname]);

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
  if(master->mpiid == 0 && dnsout == NULL)
  {
    std::string outputname = master->simname + ".out";
    dnsout = std::fopen(outputname.c_str(), "a");
    std::setvbuf(dnsout, NULL, _IOLBF, 1024);
    std::fprintf(dnsout, "%8s %11s %10s %11s %8s %8s %11s %16s %16s %16s\n",
      "ITER", "TIME", "CPUDT", "DT", "CFL", "DNUM", "DIV", "MOM", "TKE", "MASS");
    start = master->gettime();
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

    end     = master->gettime();
    cputime = end - start;
    start   = end;

    // write the output to file
    if(master->mpiid == 0)
    {
      std::fprintf(dnsout, "%8d %11.3E %10.4f %11.3E %8.4f %8.4f %11.3E %16.8E %16.8E %16.8E\n",
        iter, time, cputime, dt, cfl, dn, div, mom, tke, mass);
    }
  }

  if(doclose)
  {
    // close the output file
    if(master->mpiid == 0)
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
