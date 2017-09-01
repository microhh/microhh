/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "data_block.h"
#include "timeloop.h"
#include "boundary.h"
#include "model.h"

#ifdef USECUDA
#include <cuda_runtime_api.h>
#endif

namespace
{
    void process_command_line_options(std::string& simmode, std::string& simname,
                                      int argc, char *argv[],
                                      Master& master)
    {
        // Process the command line options.
        if (argc <= 1)
        {
            master.print_error("Specify init, run or post mode\n");
            throw 1;
        }
        else
        {
            // Check the execution mode.
            simmode = argv[1];
            if (simmode != "init" && simmode != "run" && simmode != "post")
            {
                master.print_error("Specify init, run or post mode\n");
                throw 1;
            }
            // Set the name of the simulation.
            if (argc > 2)
                simname = argv[2];
            else
                simname = "microhh";
        }

        master.print_message("Simulation name: %s\n", simname.c_str());
        master.print_message("Simulation mode: %s\n", simmode.c_str());
    }
}

// In the constructor all classes are initialized and their input is read.
template<typename TF>
Model<TF>::Model(Master *masterin, int argc, char *argv[])
{
    master = masterin;

    process_command_line_options(simmode, simname, argc, argv, *master);

    input = new Input(simname + ".ini");
    profs = new Data_block(simname + ".prof");

    // Initialize the pointers as nullptr.
    grid = nullptr;
    fields = nullptr;

    try
    {
        int nerror = 0;

        grid = new Grid<TF>(*master, *input);

        fields = new Fields<TF>(*master, *grid, *input);

        timeloop = new Timeloop<TF>(*master, *grid, *fields, *input, simmode);

        boundary = Boundary<TF>::factory(*master, *grid, *fields, *input);

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
template<typename TF>
void Model<TF>::delete_objects()
{
    // Delete the components in reversed order.
    delete boundary;
    delete fields;
    delete grid;

    delete profs;
    delete input;
}

// In the destructor the deletion of all class instances is triggered.
template<typename TF>
Model<TF>::~Model()
{
    delete_objects();
}

// In the init stage all class individual settings are known and the dynamic arrays are allocated.
template<typename TF>
void Model<TF>::init()
{
    master->init(*input);
    grid->init();
    fields->init();

    boundary->init(*input);
}

template<typename TF>
void Model<TF>::load_or_save()
{
    if (simmode == "init")
    {
        // Initialize the allocated fields and save the data.
        save();
    }
    else if (simmode == "run" || simmode == "post")
    {
        // Initialize the allocated fields using data from disk.
        load();
    }
}

// In these functions data necessary to start the model is loaded from disk.
template<typename TF>
void Model<TF>::load()
{
    // First load the grid and time to make their information available.
    grid->load();
    timeloop->load(timeloop->get_iotime());

    fields->load(timeloop->get_iotime());

    boundary->create(*input);
}

// In these functions data necessary to start the model is saved to disk.
template<typename TF>
void Model<TF>::save()
{
    // Initialize the grid and the fields from the input data.
    grid->create(*profs);
    fields->create(*input, *profs);

    // Save the initialized data to disk for the run mode.
    grid->save();
    fields->save(timeloop->get_iotime());
    timeloop->save(timeloop->get_iotime());
}

template<typename TF>
void Model<TF>::exec()
{
    if (simmode == "init")
        return;

    // #ifdef USECUDA
    // // Load all the necessary data to the GPU.
    // master  ->print_message("Preparing the GPU\n");
    // grid    ->prepare_device();
    // fields  ->prepare_device();
    // buffer  ->prepare_device();
    // thermo  ->prepare_device();
    // boundary->prepare_device();
    // diff    ->prepare_device();
    // force   ->prepare_device();
    // // Prepare pressure last, for memory check
    // pres    ->prepare_device(); 
    // #endif

    master->print_message("Starting time integration\n");

    // Update the time dependent parameters.
    // boundary->update_time_dependent();
    // force   ->update_time_dependent();

    // Set the boundary conditions.
    boundary->exec();

    // Calculate the field means, in case needed.
    // fields->exec();

    // Get the viscosity to be used in diffusion.
    // diff->exec_viscosity();

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
        boundary->set_ghost_cells_w(Boundary_w_type::Conservation_type);
        // advec->exec();
        boundary->set_ghost_cells_w(Boundary_w_type::Normal_type);

        // Calculate the diffusion tendency.
        // diff->exec();

        // Calculate the thermodynamics and the buoyancy tendency.
        // thermo->exec();

        // Calculate the tendency due to damping in the buffer layer.
        // buffer->exec();

        // Apply the large scale forcings. Keep this one always right before the pressure.
        // force->exec(timeloop->get_sub_time_step());

        // Solve the poisson equation for pressure.
        boundary->set_ghost_cells_w(Boundary_w_type::Conservation_type);
        // pres->exec(timeloop->get_sub_time_step());
        boundary->set_ghost_cells_w(Boundary_w_type::Normal_type);

        // Allow only for statistics when not in substep and not directly after restart.
        if (timeloop->is_stats_step())
        {
            // #ifdef USECUDA
            // // Copy fields from device to host
            // if (stats->doStats() || cross->do_cross() || dump->do_dump())
            // {
            //     fields  ->backward_device();
            //     boundary->backward_device();
            // }
            // #endif

            // // Do the statistics.
            // if (stats->doStats())
            // {
            //     // Always process the default mask (the full field)
            //     stats->get_mask(fields->atmp["tmp3"], fields->atmp["tmp4"], &stats->masks["default"]);
            //     calc_stats("default");

            //     // Work through the potential masks for the statistics.
            //     for (std::vector<std::string>::const_iterator it=masklist.begin(); it!=masklist.end(); ++it)
            //     {
            //         if (*it == "wplus" || *it == "wmin")
            //         {
            //             fields->get_mask(fields->atmp["tmp3"], fields->atmp["tmp4"], &stats->masks[*it]);
            //             calc_stats(*it);
            //         }
            //         else if (*it == "ql" || *it == "qlcore")
            //         {
            //             thermo->get_mask(fields->atmp["tmp3"], fields->atmp["tmp4"], &stats->masks[*it]);
            //             calc_stats(*it);
            //         }
            //         else if (*it == "patch_high" || *it == "patch_low")
            //         {
            //             boundary->get_mask(fields->atmp["tmp3"], fields->atmp["tmp4"], &stats->masks[*it]);
            //             calc_stats(*it);
            //         }
            //     }

            //     // Store the stats data.
            //     stats->exec(timeloop->get_iteration(), timeloop->get_time(), timeloop->get_itime());
            // }

            // // Save the selected cross sections to disk, cross sections are handled on CPU.
            // if (cross->do_cross())
            // {
            //     fields  ->exec_cross();
            //     thermo  ->exec_cross();
            //     boundary->exec_cross();
            // }

            // // Save the 3d dumps to disk
            // if (dump->do_dump())
            // {
            //     fields->exec_dump();
            //     thermo->exec_dump();
            // }
        }

        // Exit the simulation when the runtime has been hit.
        if (timeloop->is_finished())
            break;

        // RUN MODE: In case of run mode do the time stepping.
        if (simmode == "run")
        {
            // Integrate in time.
            timeloop->exec();

            // Increase the time with the time step.
            timeloop->step_time();

            // Save the data for restarts.
            if (timeloop->do_save())
            {
                // #ifdef USECUDA
                // fields  ->backward_device();
                // boundary->backward_device();
                // #endif

                // Save data to disk.
                timeloop->save(timeloop->get_iotime());
                fields  ->save(timeloop->get_iotime());
            }
        }

        // POST PROCESS MODE: In case of post-process mode, load a new set of files.
        else if (simmode == "post")
        {
            // Step to the next time step.
            timeloop->step_post_proc_time();

            // In case the simulation is done, step out of the loop.
            if (timeloop->is_finished())
                break;

            // Load the data from disk.
            timeloop->load(timeloop->get_iotime());
            fields  ->load(timeloop->get_iotime());
        }

        // Update the time dependent parameters.
        // boundary->update_time_dependent();
        // force   ->update_time_dependent();

        // Set the boundary conditions.
        boundary->exec();

        // Calculate the field means, in case needed.
        // fields->exec();

        // Get the viscosity to be used in diffusion.
        // diff->exec_viscosity();

        // Write status information to disk.
        print_status();

    } // End time loop.

    // #ifdef USECUDA
    // // At the end of the run, copy the data back from the GPU.
    // fields  ->backward_device();
    // boundary->backward_device();
    // #endif
}

template<typename TF>
void Model<TF>::set_time_step()
{
    // Only set the time step if the model is not in a substep.
    if (timeloop->in_substep())
        return;

    // Retrieve the maximum allowed time step per class.
    timeloop->set_time_step_limit();
    // timeloop->set_time_step_limit(advec ->get_time_limit(timeloop->get_idt(), timeloop->get_dt()));
    // timeloop->set_time_step_limit(diff  ->get_time_limit(timeloop->get_idt(), timeloop->get_dt()));
    // timeloop->set_time_step_limit(thermo->get_time_limit(timeloop->get_idt(), timeloop->get_dt()));
    // timeloop->set_time_step_limit(stats ->get_time_limit(timeloop->get_itime()));
    // timeloop->set_time_step_limit(cross ->get_time_limit(timeloop->get_itime()));
    // timeloop->set_time_step_limit(dump  ->get_time_limit(timeloop->get_itime()));

    // Set the time step.
    timeloop->set_time_step();
}

// Calculate the statistics for all classes that have a statistics function.
template<typename TF>
void Model<TF>::calc_stats(std::string maskname)
{
}

// Print the status information to the .out file.
template<typename TF>
void Model<TF>::print_status()
{
}

template class Model<double>;
template class Model<float>;
