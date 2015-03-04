/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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
#include "dump.h"
#include "budget.h"

#ifdef USECUDA
#include <cuda_runtime_api.h>
#endif

// In the constructor all classes are initialized and their input is read.
Model::Model(Master *masterin, Input *inputin)
{
  master = masterin;
  input  = inputin;

  // Initialize the pointers as zero
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
  dump   = 0;
  budget = 0;

  try
  {
    // Create an instance of the Grid class.
    grid = new Grid(this, input);

    // Create an instance of the Fields class.
    fields = new Fields(this, input);

    // Create instances of the other model classes.
    boundary = Boundary::factory(master, input, this);
    advec    = Advec   ::factory(master, input, this, grid->swspatialorder);
    diff     = Diff    ::factory(master, input, this, grid->swspatialorder);
    pres     = Pres    ::factory(master, input, this, grid->swspatialorder);
    thermo   = Thermo  ::factory(master, input, this);

    timeloop = new Timeloop(this, input);
    force    = new Force   (this, input);
    buffer   = new Buffer  (this, input);

    // Create instances of the statistics classes.
    stats  = new Stats (this, input);
    cross  = new Cross (this, input);
    dump   = new Dump  (this, input);
    budget = new Budget(this, input);

    // Get the list of masks.
    // TODO Make an interface that takes this out of the main loop.
    int nerror = 0;
    nerror += input->get_list(&masklist, "stats", "masklist", "");
    for (std::vector<std::string>::const_iterator it=masklist.begin(); it!=masklist.end(); ++it)
    {
      if (*it != "wplus" &&
         *it != "wmin"  &&
         *it != "ql"    &&
         *it != "qlcore")
      {
        master->print_warning("%s is an undefined mask for conditional statistics\n", it->c_str());
      }
      else if ((*it == "ql" || *it == "qlcore") && thermo->getSwitch() != "moist")
        master->print_warning("%s mask only works for swthermo=moist \n", it->c_str());
      else
        stats->addMask(*it);
    }

    // if one or more arguments fails, then crash
    if (nerror > 0)
      throw 1;
  }
  catch (int &e)
  {
    // In case of a failing constructor, delete the class objects and rethrow.
    delete_objects();
    throw;
  }
}

// In this function all instances of objects are deleted and the memory is freed.
void Model::delete_objects()
{
  // Delete the components in reversed order.
  delete budget;
  delete dump;
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

// In the destructor the deletion of all class instances is triggered.
Model::~Model()
{
  delete_objects();
  #ifdef USECUDA
  cudaDeviceReset();
  #endif
}

// In the init stage all class individual settings are known and the dynamic arrays are allocated.
void Model::init()
{
  grid  ->init();
  fields->init();

  boundary->init(input);
  buffer  ->init();
  force   ->init();
  pres    ->init();
  thermo  ->init();

  stats ->init(timeloop->get_ifactor());
  cross ->init(timeloop->get_ifactor());
  dump  ->init(timeloop->get_ifactor());
  budget->init();
}

// In these functions data necessary to start the model is loaded from disk.
void Model::load()
{
  // First load the grid and time to make their information available.
  grid    ->load();
  timeloop->load(timeloop->get_iotime());

  // Initialize the statistics file to open the possiblity to add profiles.
  stats->create(timeloop->get_iotime());
  cross->create();
  dump ->create();

  fields->load(timeloop->get_iotime());
  fields->create_stats();

  // Initialize data or load data from disk.
  boundary->create(input);

  buffer->create(input);
  force ->create(input);
  thermo->create(input);

  budget->create();

  // End with those modules that require all fields to be loaded.
  boundary->set_values();
  diff    ->set_values();
  pres    ->set_values();
}

// In these functions data necessary to start the model is saved to disk.
void Model::save()
{
  // Initialize the grid and the fields from the input data.
  grid  ->create(input);
  fields->create(input);

  // Save the initialized data to disk for the run mode.
  grid    ->save();
  fields  ->save(timeloop->get_iotime());
  timeloop->save(timeloop->get_iotime());
}

void Model::exec()
{
  #ifdef USECUDA
  // Load all the necessary data to the GPU.
  master  ->print_message("Preparing the GPU\n");
  grid    ->prepare_device();
  fields  ->prepare_device();
  buffer  ->prepare_device();
  thermo  ->prepare_device();
  boundary->prepare_device();
  diff    ->prepare_device();
  force   ->prepare_device();
  // Prepare pressure last, for memory check
  pres    ->prepare_device(); 
  #endif

  master->print_message("Starting time integration\n");

  // Update the time dependent parameters.
  boundary->update_time_dependent();
  force   ->update_time_dependent();

  // Set the boundary conditions.
  boundary->exec();

  // Calculate the field means, in case needed.
  fields->exec();

  // Get the viscosity to be used in diffusion.
  diff->exec_viscosity();

  // Set the time step.
  set_time_step();

  // Print the initial status information.
  print_status();

  // start the time loop
  while (true)
  {
    // Determine the time step.
    set_time_step();

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

    // Solve the poisson equation for pressure.
    pres->exec(timeloop->getSubTimeStep());

    // Allow only for statistics when not in substep and not directly after restart.
    if (timeloop->isStatsStep())
    {
      #ifdef USECUDA
      // Copy fields from device to host
      if (stats->doStats() || cross->do_cross() || dump->do_dump())
      {
        fields  ->backward_device();
        boundary->backward_device();
      }
      #endif

      // Do the statistics.
      if (stats->doStats())
      {
        // Always process the default mask (the full field)
        stats->get_mask(fields->atmp["tmp3"], fields->atmp["tmp4"], &stats->masks["default"]);
        calc_stats("default");

        // Work through the potential masks for the statistics.
        for (std::vector<std::string>::const_iterator it=masklist.begin(); it!=masklist.end(); ++it)
        {
          if (*it == "wplus" || *it == "wmin")
          {
            fields->get_mask(fields->atmp["tmp3"], fields->atmp["tmp4"], &stats->masks[*it]);
            calc_stats(*it);
          }
          else if (*it == "ql" || *it == "qlcore")
          {
            thermo->get_mask(fields->atmp["tmp3"], fields->atmp["tmp4"], &stats->masks[*it]);
            calc_stats(*it);
          }
        }

        // Store the stats data.
        stats->exec(timeloop->get_iteration(), timeloop->get_time(), timeloop->get_itime());
      }

      // Save the selected cross sections to disk, cross sections are handled on CPU.
      if (cross->do_cross())
      {
        fields  ->exec_cross();
        thermo  ->exec_cross();
        boundary->exec_cross();
      }

      // Save the 3d dumps to disk
      if (dump->do_dump())
      {
        fields->exec_dump();
        thermo->exec_dump();
      }
    }

    // Exit the simulation when the runtime has been hit.
    if (timeloop->isFinished())
      break;

    // RUN MODE: In case of run mode do the time stepping.
    if (master->mode == "run")
    {
      // Integrate in time.
      timeloop->exec();

      // Increase the time with the time step.
      timeloop->stepTime();

      // Save the data for restarts.
      if (timeloop->doSave())
      {
        #ifdef USECUDA
        fields  ->backward_device();
        boundary->backward_device();
        #endif

        // Save data to disk.
        timeloop->save(timeloop->get_iotime());
        fields  ->save(timeloop->get_iotime());
      }
    }

    // POST PROCESS MODE: In case of post-process mode, load a new set of files.
    else if (master->mode == "post")
    {
      // Step to the next time step.
      timeloop->stepPostProcTime();

      // In case the simulation is done, step out of the loop.
      if (timeloop->isFinished())
        break;

      // Load the data from disk.
      timeloop->load(timeloop->get_iotime());
      fields  ->load(timeloop->get_iotime());
    }

    // Update the time dependent parameters.
    boundary->update_time_dependent();
    force   ->update_time_dependent();

    // Set the boundary conditions.
    boundary->exec();

    // Calculate the field means, in case needed.
    fields->exec();

    // Get the viscosity to be used in diffusion.
    diff->exec_viscosity();

    // Write status information to disk.
    print_status();

  } // End time loop.

  #ifdef USECUDA
  // At the end of the run, copy the data back from the GPU.
  fields  ->backward_device();
  boundary->backward_device();
  #endif
}

void Model::set_time_step()
{
  // Only set the time step if the model is not in a substep.
  if (timeloop->inSubStep())
    return;

  // Retrieve the maximum allowed time step per class.
  timeloop->setTimeStepLimit();
  timeloop->setTimeStepLimit(advec->get_time_limit(timeloop->get_idt(), timeloop->get_dt()));
  timeloop->setTimeStepLimit(diff ->get_time_limit(timeloop->get_idt(), timeloop->get_dt()));
  timeloop->setTimeStepLimit(stats->get_time_limit(timeloop->get_itime()));
  timeloop->setTimeStepLimit(cross->get_time_limit(timeloop->get_itime()));
  timeloop->setTimeStepLimit(dump ->get_time_limit(timeloop->get_itime()));

  // Set the time step.
  timeloop->setTimeStep();
}

// Calculate the statistics for all classes that have a statistics function.
void Model::calc_stats(std::string maskname)
{
  fields  ->exec_stats(&stats->masks[maskname]);
  thermo  ->exec_stats(&stats->masks[maskname]);
  budget  ->exec_stats(&stats->masks[maskname]);
  boundary->exec_stats(&stats->masks[maskname]);
}

// Print the status information to the .out file.
void Model::print_status()
{
  // Initialize the check variables
  int    iter;
  double time, dt;
  double mom, tke, mass;
  double div;
  double cfl, dn;
  double cputime, end;
  static double start;
  static FILE *dnsout = NULL;

  // Write output file header on the main process and set the time of writing.
  if (master->mpiid == 0 && dnsout == NULL)
  {
    std::string outputname = master->simname + ".out";
    dnsout = std::fopen(outputname.c_str(), "a");
    std::setvbuf(dnsout, NULL, _IOLBF, 1024);
    std::fprintf(dnsout, "%8s %11s %10s %11s %8s %8s %11s %16s %16s %16s\n",
      "ITER", "TIME", "CPUDT", "DT", "CFL", "DNUM", "DIV", "MOM", "TKE", "MASS");
    start = master->get_wall_clock_time();
  }

  // Retrieve all the status information.
  if (timeloop->doCheck())
  {
    // Get status variables.
    iter = timeloop->get_iteration();
    time = timeloop->get_time();
    dt   = timeloop->get_dt();
    div  = pres->checkDivergence();
    mom  = fields->check_momentum();
    tke  = fields->check_tke();
    mass = fields->check_mass();
    cfl  = advec->get_cfl(timeloop->get_dt());
    dn   = diff->get_dn(timeloop->get_dt());

    // Store time interval in betwteen two writes.
    end     = master->get_wall_clock_time();
    cputime = end - start;
    start   = end;

    // Write the status information to disk.
    if (master->mpiid == 0)
      std::fprintf(dnsout, "%8d %11.3E %10.4f %11.3E %8.4f %8.4f %11.3E %16.8E %16.8E %16.8E\n",
        iter, time, cputime, dt, cfl, dn, div, mom, tke, mass);
  }

  if (timeloop->isFinished())
  {
    // Close the output file when the run is done.
    if (master->mpiid == 0)
      std::fclose(dnsout);
  }
}
