/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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
#include "advec.h"
#include "diff.h"
#include "pres.h"
#include "thermo.h"
#include "boundary.h"
#include "buffer.h"
#include "force.h"
#include "stats.h"
#include "cross.h"
#include "budget.h"

#ifdef USECUDA
#include <cuda_runtime_api.h> // Needed for cudaDeviceReset(), to check mem leaks 
#endif

Model::Model(Master *masterin, Input *inputin)
{
  master = masterin;
  input  = inputin;

  // initialize the pointers at zero
  grid     = 0;
  fields   = 0;
  diff     = 0;
  pres     = 0;
  thermo   = 0;
  timeloop = 0;
  force    = 0;
  buffer   = 0;

  stats  = 0;
  cross  = 0;
  budget = 0;

  try
  {
    // create the grid class
    grid = new Grid(this, input);

    // create the fields class
    fields = new Fields(this, input);

    // create the model components
    boundary = Boundary::factory(master, input, this);
    advec    = Advec   ::factory(master, input, this, grid->swspatialorder);
    diff     = Diff    ::factory(master, input, this, grid->swspatialorder);
    pres     = Pres    ::factory(master, input, this, grid->swspatialorder);
    thermo   = Thermo  ::factory(master, input, this);

    timeloop = new Timeloop(this, input);
    force    = new Force   (this, input);
    buffer   = new Buffer  (this, input);

    // load the postprocessing modules
    stats  = new Stats (this, input);
    cross  = new Cross (this, input);
    budget = new Budget(this, input);

    // get the list of masks
    // TODO This is really UGLY: make an interface that takes this out of the main loops
    int nerror = 0;
    nerror += input->getList(&masklist, "stats", "masklist", "");
    for(std::vector<std::string>::const_iterator it=masklist.begin(); it!=masklist.end(); ++it)
    {
      if(*it != "wplus" &&
         *it != "wmin"  &&
         *it != "ql"    &&
         *it != "qlcore")
      {
        master->printWarning("%s is an undefined mask for conditional statistics\n", it->c_str());
      }
      else
        stats->addMask(*it);
    }

    // Get base state option (boussinesq or anelastic)
    nerror += input->getItem(&swbasestate , "grid"  , "swbasestate" , "", "");  // BvS: where to put switch??

    if(!(swbasestate == "boussinesq" || swbasestate == "anelastic"))
    {
      master->printError("\"%s\" is an illegal value for swbasestate\n", swbasestate.c_str());
      throw 1;
    }
    // 2nd order with LES diffusion is only option supporting anelastic (for now) 
    //if(swdiff != "smag2" and swbasestate != "boussinesq")
    //{
    //  std::printf("ERROR swdiff=%s is not allowed with swbasestate=%s \n", swdiff.c_str(),swbasestate.c_str());
    //  return 1;
    //}

    // if one or more arguments fails, then crash
    if(nerror > 0)
      throw 1;
  }
  catch (int &e)
  {
    // In case of a failing constructor, delete the class objects and rethrow.
    deleteObjects();
    throw;
  }
}

void Model::deleteObjects()
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

Model::~Model()
{
  deleteObjects();
#ifdef USECUDA
  cudaDeviceReset();
#endif
}

void Model::init()
{
  grid  ->init();
  fields->init();

  boundary->init(input);
  buffer  ->init();
  force   ->init();
  pres    ->init();
  thermo  ->init();

  stats ->init(timeloop->ifactor);
  cross ->init(timeloop->ifactor);
  budget->init();
}

void Model::load()
{
  // first load the grid and time to make their information available
  grid    ->load();
  timeloop->load(timeloop->iotime);

  // initialize the statistics file to open the possiblity to add profiles
  stats->create(timeloop->iotime);
  cross->create();

  fields->load(timeloop->iotime);

  // \TODO call boundary load for the data and then timedep, not nice...
  boundary->load(timeloop->iotime);
  boundary->create(input);

  buffer->create(input);
  force ->create(input);
  thermo->create(input);

  budget->create();

  // end with modules that require all fields to be present
  boundary->setValues();
  diff    ->setValues();
  pres    ->setValues();
}

void Model::save()
{
  // Initialize the grid and the fields from the input data.
  grid  ->create(input);
  fields->create(input);

  // Save the initialized data to disk for the run mode.
  grid    ->save();
  fields  ->save(timeloop->iotime);
  timeloop->save(timeloop->iotime);
  boundary->save(timeloop->iotime);
}

void Model::exec()
{
#ifdef USECUDA
  master  ->printMessage("Preparing the GPU\n");
  grid    ->prepareDevice();
  fields  ->prepareDevice();
  pres    ->prepareDevice();
  buffer  ->prepareDevice();
  thermo  ->prepareDevice();
  boundary->prepareDevice();
  diff    ->prepareDevice();
  force   ->prepareDevice();
#endif

  master->printMessage("Starting time integration\n");

  // update the time dependent values
  boundary->setTimeDep();
  force->setTimeDep();

  // set the boundary conditions
  boundary->exec();

  // get the field means, in case needed
  fields->exec();
  // get the viscosity to be used in diffusion
  diff->execViscosity();

  setTimeStep();

  // print the initial information
  printOutputFile(!timeloop->loop);

  // start the time loop
  while(true)
  {
    // Determine the time step.
    setTimeStep();

    // Calculate the advection tendency.
    advec->exec();
    // Calculate the diffusion tendency.
    diff->exec();
    // Calculate the thermodynamics and the buoyancy tendency.
    thermo->exec();
    // Calculate the tendency due to damping in the buffer layer.
    buffer->exec();

    // Apply the large scale forcings. Keep this one always right before the pressure.
    force->exec(timeloop->getSubTimeStep());

    // pressure
    pres->exec(timeloop->getSubTimeStep());

    // Do only data analysis statistics when not in substep and not directly after restart.
    if(timeloop->isStatsStep())
    {
      if(stats->doStats())
      {
        #ifdef USECUDA
        fields  ->backwardDevice();
        boundary->backwardDevice();
        #endif

        // always process the default mask
        stats->getMask(fields->atmp["tmp3"], fields->atmp["tmp4"], &stats->masks["default"]);
        calcStats("default");

        // Work through the potential masks for the statistics.
        for(std::vector<std::string>::const_iterator it=masklist.begin(); it!=masklist.end(); ++it)
        {
          if(*it == "wplus" || *it == "wmin")
          {
            fields->getMask(fields->atmp["tmp3"], fields->atmp["tmp4"], &stats->masks[*it]);
            calcStats(*it);
          }
          else if(*it == "ql" || *it == "qlcore")
          {
            thermo->getMask(fields->atmp["tmp3"], fields->atmp["tmp4"], &stats->masks[*it]);
            calcStats(*it);
          }
        }

        // Store the stats data.
        stats->exec(timeloop->iteration, timeloop->time, timeloop->itime);
      }

      if(cross->doCross())
      {
        // Copy back the data from the GPU
        #ifdef USECUDA
        fields  ->backwardDevice();
        boundary->backwardDevice();
        #endif
      
        fields  ->execCross();
        thermo  ->execCross();
        boundary->execCross();
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
      timeloop->stepTime();

      // save the data for a restart
      if(timeloop->doSave())
      {
        #ifdef USECUDA
        fields  ->backwardDevice();
        boundary->backwardDevice();
        #endif

        // Save data to disk.
        timeloop->save(timeloop->iotime);
        fields  ->save(timeloop->iotime);
        boundary->save(timeloop->iotime);
      }
    }

    // POST PROCESS MODE
    else if(master->mode == "post")
    {
      // step to the next time step
      timeloop->stepPostProcTime();

      // if simulation is done break
      if(!timeloop->loop)
        break;

      // load the data
      timeloop->load(timeloop->iotime);
      fields  ->load(timeloop->iotime);
      boundary->load(timeloop->iotime);
    }
    // update the time dependent values
    boundary->setTimeDep();
    force   ->setTimeDep();

    // set the boundary conditions
    boundary->exec();
    // get the field means, in case needed
    fields->exec();
    // get the viscosity to be used in diffusion
    diff->execViscosity();

    printOutputFile(!timeloop->loop);
  } // end time loop

  #ifdef USECUDA
  fields  ->backwardDevice();
  boundary->backwardDevice();
  #endif
}

void Model::calcStats(std::string maskname)
{
  fields  ->execStats(&stats->masks[maskname]);
  thermo  ->execStats(&stats->masks[maskname]);
  budget  ->execStats(&stats->masks[maskname]);
  boundary->execStats(&stats->masks[maskname]);
}

void Model::printOutputFile(bool doclose)
{
  // initialize the check variables
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
    start = master->getTime();
  }

  if(timeloop->doCheck())
  {
    iter    = timeloop->iteration;
    time    = timeloop->time;
    dt      = timeloop->dt;
    div     = pres->checkDivergence();
    mom     = fields->checkMomentum();
    tke     = fields->checkTke();
    mass    = fields->checkMass();
    cfl     = advec->getcfl(timeloop->dt);
    dn      = diff->getdn(timeloop->dt);

    end     = master->getTime();
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
}

void Model::setTimeStep()
{
  // Only set the time step if the model is not in a substep
  if(timeloop->inSubStep())
    return;

  timeloop->setTimeLimit();

  timeloop->idtlim = std::min(timeloop->idtlim, advec->getTimeLimit(timeloop->idt, timeloop->dt));
  timeloop->idtlim = std::min(timeloop->idtlim, diff ->getTimeLimit(timeloop->idt, timeloop->dt));
  timeloop->idtlim = std::min(timeloop->idtlim, stats->getTimeLimit(timeloop->itime));
  timeloop->idtlim = std::min(timeloop->idtlim, cross->getTimeLimit(timeloop->itime));
  timeloop->setTimeStep();
}
