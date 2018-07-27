/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
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
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "buffer.h"
#include "data_block.h"
#include "timeloop.h"
#include "fft.h"
#include "boundary.h"
#include "advec.h"
#include "diff.h"
#include "pres.h"
#include "force.h"
#include "thermo.h"
#include "radiation.h"
#include "microphys.h"
#include "decay.h"
#include "stats.h"
#include "column.h"
#include "cross.h"
#include "dump.h"
#include "model.h"

#ifdef USECUDA
#include <cuda_runtime_api.h>
#endif

namespace
{
    void process_command_line_options(Sim_mode& sim_mode, std::string& sim_name,
                                      int argc, char *argv[],
                                      Master& master)
    {
        // Process the command line options.
        if (argc <= 1)
        {
            master.print_error("Specify init, run or post mode\n");
            throw std::runtime_error("No run mode specified");
        }
        else
        {
            // Check the execution mode.
            std::string sim_mode_str = argv[1];
            if (sim_mode_str != "init" && sim_mode_str != "run" && sim_mode_str != "post")
            {
                master.print_error("Specify init, run or post mode\n");
                throw std::runtime_error("Illegal run mode specified");
            }
            else
            {
                if (sim_mode_str == "init")
                    sim_mode = Sim_mode::Init;
                else if (sim_mode_str == "run")
                    sim_mode = Sim_mode::Run;
                else
                    sim_mode = Sim_mode::Post;
            }

            // Set the name of the simulation.
            if (argc > 2)
                sim_name = argv[2];
            else
                sim_name = "microhh";

            master.print_message("Simulation name: %s\n", sim_name.c_str());
            master.print_message("Simulation mode: %s\n", sim_mode_str.c_str());
        }
    }
}

// In the constructor all classes are initialized and their input is read.
template<typename TF>
Model<TF>::Model(Master& masterin, int argc, char *argv[]) :
    master(masterin)
{
    process_command_line_options(sim_mode, sim_name, argc, argv, master);

    input = std::make_shared<Input>(master, sim_name + ".ini");
    profs = std::make_shared<Data_block>(master, sim_name + ".prof");

    try
    {
        grid      = std::make_shared<Grid<TF>>(master, *input);
        fields    = std::make_shared<Fields<TF>>(master, *grid, *input);
        timeloop  = std::make_shared<Timeloop<TF>>(master, *grid, *fields, *input, sim_mode);
        fft       = std::make_shared<FFT<TF>>(master, *grid);

        boundary  = Boundary<TF> ::factory(master, *grid, *fields, *input);
        advec     = Advec<TF>    ::factory(master, *grid, *fields, *input);
        diff      = Diff<TF>     ::factory(master, *grid, *fields, *input);
        pres      = Pres<TF>     ::factory(master, *grid, *fields, *fft, *input);
        thermo    = Thermo<TF>   ::factory(master, *grid, *fields, *input);
        microphys = Microphys<TF>::factory(master, *grid, *fields, *input);

        radiation = std::make_shared<Radiation<TF>>(master, *grid, *fields, *input);
        force     = std::make_shared<Force    <TF>>(master, *grid, *fields, *input);
        buffer    = std::make_shared<Buffer   <TF>>(master, *grid, *fields, *input);
        decay     = std::make_shared<Decay    <TF>>(master, *grid, *fields, *input);
        stats     = std::make_shared<Stats    <TF>>(master, *grid, *fields, *input);
        column    = std::make_shared<Column   <TF>>(master, *grid, *fields, *input);
        dump      = std::make_shared<Dump     <TF>>(master, *grid, *fields, *input);
        cross     = std::make_shared<Cross    <TF>>(master, *grid, *fields, *input);

        // Parse the statistics masks
        add_statistics_masks();
    }
    catch (std::exception& e)
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
}

// In the destructor the deletion of all class instances is triggered.
template<typename TF>
Model<TF>::~Model()
{
    delete_objects();
    #pragma omp taskwait
}

// In the init stage all class individual settings are known and the dynamic arrays are allocated.
template<typename TF>
void Model<TF>::init()
{
    master.init(*input);

    grid->init();
    fields->init(*dump, *cross);

    fft->init();

    boundary->init(*input, *thermo);
    buffer->init();
    diff->init();
    pres->init();
    force->init();
    thermo->init();
    microphys->init();
    radiation->init();
    decay->init(*input);

    stats->init(timeloop->get_ifactor());
    column->init(timeloop->get_ifactor());
    cross->init(timeloop->get_ifactor());
    dump->init(timeloop->get_ifactor());
}

template<typename TF>
void Model<TF>::load_or_save()
{
    if (sim_mode == Sim_mode::Init)
    {
        // Initialize the allocated fields and save the data.
        save();
    }
    else if (sim_mode == Sim_mode::Run || sim_mode == Sim_mode::Post)
    {
        // Initialize the allocated fields using data from disk.
        load();
    }

    // This marks the end of the entire initialization.
    // Print warnings for input variables that are unused.
    input->print_unused_items();

    // Free the memory taken by the input fields.
    input.reset();

}

// In these functions data necessary to start the model is loaded from disk.
template<typename TF>
void Model<TF>::load()
{
    // First load the grid and time to make their information available.
    grid->load();
    fft->load();
    timeloop->load(timeloop->get_iotime());

    // Initialize the statistics file to open the possiblity to add profiles in other routines
    stats->create(timeloop->get_iotime(), sim_name);
    column->create(*input, timeloop->get_iotime(), sim_name);
    dump->create();

    // Load the fields, and create the field statistics
    fields->load(timeloop->get_iotime());
    fields->create_stats(*stats);
    fields->create_column(*column);

    boundary->create(*input, *stats);
    buffer->create(*input, *profs);
    force->create(*input, *profs);
    thermo->create(*input, *profs, *stats, *column, *cross, *dump);
    microphys->create(*input, *profs, *stats, *cross, *dump);
    radiation->create(*thermo); // Radiation needs to be created after thermo as it needs base profiles.
    decay->create(*input);

    cross->create();    // Cross needs to be called at the end!

    boundary->set_values();
    pres->set_values();
    diff->create(*stats);
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
    fft->save();
    fields->save(timeloop->get_iotime());
    timeloop->save(timeloop->get_iotime());
}

template<typename TF>
void Model<TF>::exec()
{
    if (sim_mode == Sim_mode::Init)
        return;

    #ifdef USECUDA
    prepare_gpu();
    #endif

    master.print_message("Starting time integration\n");

    // Update the time dependent parameters.
    boundary->update_time_dependent(*timeloop);
    force   ->update_time_dependent(*timeloop);
    thermo  ->update_time_dependent(*timeloop);

    // Set the boundary conditions.
    boundary->exec(*thermo);

    // Calculate the field means, in case needed.
    fields->exec();

    // Get the viscosity to be used in diffusion.
    diff->exec_viscosity(*boundary, *thermo);

    // Set the time step.
    set_time_step();

    // Print the initial status information.
    print_status();
    #ifdef USECUDA
        #ifdef _OPENMP
            omp_set_nested(1);
            const int nthreads_out=2;
            master.print_message("Running with %i OpenMP threads\n", omp_get_max_threads());
        #endif
    #else
        #ifdef _OPENMP
            omp_set_num_threads(1);
            const int nthreads_out=1;
        #endif
    #endif
    #pragma omp parallel num_threads(nthreads_out)
    {
    #pragma omp master
    {
        // start the time loop
        while (true)
        {
            // Determine the time step.
            set_time_step();

            // Calculate the advection tendency.
            boundary->set_ghost_cells_w(Boundary_w_type::Conservation_type);
            advec->exec();
            boundary->set_ghost_cells_w(Boundary_w_type::Normal_type);

            // Calculate the diffusion tendency.
            diff->exec(*boundary);

            // Calculate the thermodynamics and the buoyancy tendency.
            thermo->exec(timeloop->get_sub_time_step());

            // Calculate the microphysics.
            microphys->exec(*thermo, timeloop->get_sub_time_step());

            // Calculate the radiation fluxes and the related heating rate.
            radiation->exec(*thermo);

            // Calculate the tendency due to damping in the buffer layer.
            buffer->exec();

            // Apply the scalar decay.
            decay->exec(timeloop->get_sub_time_step());

            // Apply the large scale forcings. Keep this one always right before the pressure.
            force->exec(timeloop->get_sub_time_step());

            // Solve the poisson equation for pressure.
            boundary->set_ghost_cells_w(Boundary_w_type::Conservation_type);
            pres->exec(timeloop->get_sub_time_step());
            boundary->set_ghost_cells_w(Boundary_w_type::Normal_type);

            // Allow only for statistics when not in substep and not directly after restart.
            if (timeloop->is_stats_step())
            {
                if (stats->do_statistics(timeloop->get_itime()) || cross->do_cross(timeloop->get_itime()) ||
                    dump->do_dump(timeloop->get_itime()) || column->do_column(timeloop->get_itime()))
                {
                    #pragma omp taskwait
                    #ifdef USECUDA
                        fields  ->backward_device();
                        boundary->backward_device();
                        thermo  ->backward_device();
                    #endif
                    #pragma omp task default(shared)
                    calculate_statistics(timeloop->get_iteration(), timeloop->get_time(), timeloop->get_itime(), timeloop->get_iotime(), timeloop->get_dt());

                }

            }

            // Exit the simulation when the runtime has been hit.
            if (timeloop->is_finished())
                break;

            // RUN MODE: In case of run mode do the time stepping.
            if (sim_mode == Sim_mode::Run)
            {
                // Integrate in time.
                timeloop->exec();

                // Increase the time with the time step.
                timeloop->step_time();

                // Save the data for restarts.
                if (timeloop->do_save())
                {
                    #pragma omp taskwait
                    #ifdef USECUDA
                    fields  ->backward_device();
                    boundary->backward_device();
                    thermo  ->backward_device();
                    #endif

                    // Save data to disk.
                    timeloop->save(timeloop->get_iotime());
                    fields  ->save(timeloop->get_iotime());
                }
            }

            // POST PROCESS MODE: In case of post-process mode, load a new set of files.
            else if (sim_mode == Sim_mode::Post)
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
            boundary->update_time_dependent(*timeloop);
            thermo  ->update_time_dependent(*timeloop);
            force   ->update_time_dependent(*timeloop);

            // Set the boundary conditions.
            boundary->exec(*thermo);

            // Calculate the field means, in case needed.
            // fields->exec();

            // Get the viscosity to be used in diffusion.
            diff->exec_viscosity(*boundary, *thermo);

            // Write status information to disk.
            print_status();

        } // End time loop.
    }}  //End OpenMP master and parallel
    #ifdef USECUDA
    // At the end of the run, copy the data back from the GPU.
    fields  ->backward_device();
    // boundary->backward_device();
    thermo  ->backward_device();

    clear_gpu();
    #endif
}

#ifdef USECUDA
template<typename TF>
void Model<TF>::prepare_gpu()
{
    // Load all the necessary data to the GPU.
    master.print_message("Preparing the GPU\n");
    grid    ->prepare_device();
    fields  ->prepare_device();
    buffer  ->prepare_device();
    thermo  ->prepare_device();
    boundary->prepare_device();
    diff    ->prepare_device(*boundary);
    force   ->prepare_device();
    // decay   ->prepare_device();
    // // Prepare pressure last, for memory check
    pres    ->prepare_device();
}

template<typename TF>
void Model<TF>::clear_gpu()
{
    master.print_message("Clearing the GPU\n");
    grid    ->clear_device();
    fields  ->clear_device();
    // buffer  ->clear_device();
    thermo  ->clear_device();
    // boundary->clear_device();
    // diff    ->clear_device();
    force   ->clear_device();
    // decay   ->clear_device();
    // // Clear pressure last, for memory check
    pres    ->clear_device();
}
#endif

// Calculate the statistics for all classes that have a statistics function.
template<typename TF>
void Model<TF>::calculate_statistics(int iteration, double time, unsigned long itime, int iotime, double dt)
{
    // Do the statistics.
    if(stats->do_statistics(itime))
    {
        const std::vector<std::string>& mask_list = stats->get_mask_list();

        stats->initialize_masks();
        for (auto& mask_name : mask_list)
        {
            // Get the mask from one of the mask providing classes
            if (fields->has_mask(mask_name))
                fields->get_mask(*stats, mask_name);
            else if (thermo->has_mask(mask_name))
                thermo->get_mask(*stats, mask_name);
            else if (microphys->has_mask(mask_name))
                microphys->get_mask(*stats, mask_name);
            else
            {
                std::string error_message = "Can not calculate mask for \"" + mask_name + "\"";
                throw std::runtime_error(error_message);
            }
        }
        stats->finalize_masks();
        // Calculate statistics
        fields   ->exec_stats(*stats, *diff);
        thermo   ->exec_stats(*stats, *diff);
        microphys->exec_stats(*stats, *diff, *thermo, dt);
        diff     ->exec_stats(*stats);
        //budget  ->exec_stats(&stats->masks[maskname]);
        boundary ->exec_stats(*stats);

        // Store the statistics data.
        stats->exec(iteration, time, itime);
}

//    // Save the selected cross sections to disk, cross sections are handled on CPU.
   if(cross->do_cross(itime))
    {
        fields   ->exec_cross(*cross, iotime);
        thermo   ->exec_cross(*cross, iotime);
        microphys->exec_cross(*cross, iotime);
//        boundary->exec_cross(iotime);
    }
   // Save the 3d dumps to disk
    if(dump->do_dump(itime))
    {
        fields   ->exec_dump(*dump, iotime);
        thermo   ->exec_dump(*dump, iotime);
        microphys->exec_dump(*dump, iotime);
    }

    if(column->do_column(itime))
    {
        fields->exec_column(*column);
        thermo->exec_column(*column);
        column->exec(iteration, time, itime);
    }
}

template<typename TF>
void Model<TF>::set_time_step()
{
    // Only set the time step if the model is not in a substep.
    if (timeloop->in_substep())
        return;

    // Retrieve the maximum allowed time step per class.
    timeloop->set_time_step_limit();
    timeloop->set_time_step_limit(advec    ->get_time_limit(timeloop->get_idt(), timeloop->get_dt()));
    timeloop->set_time_step_limit(diff     ->get_time_limit(timeloop->get_idt(), timeloop->get_dt()));
    timeloop->set_time_step_limit(thermo   ->get_time_limit(timeloop->get_idt(), timeloop->get_dt()));
    timeloop->set_time_step_limit(microphys->get_time_limit(timeloop->get_idt(), timeloop->get_dt()));
    timeloop->set_time_step_limit(stats    ->get_time_limit(timeloop->get_itime()));
    timeloop->set_time_step_limit(cross    ->get_time_limit(timeloop->get_itime()));
    timeloop->set_time_step_limit(dump     ->get_time_limit(timeloop->get_itime()));
    timeloop->set_time_step_limit(column   ->get_time_limit(timeloop->get_itime()));

    // Set the time step.
    timeloop->set_time_step();
}

// Add all masks
template<typename TF>
void Model<TF>::add_statistics_masks()
{
    const std::vector<std::string>& mask_list = stats->get_mask_list();

    // Check whether the mask can be retrieved from any of the mask-providing classes
    for (auto& mask_name : mask_list)
    {
        if (mask_name == "default")
            stats->add_mask(mask_name);
        else if (fields->has_mask(mask_name))
            stats->add_mask(mask_name);
        else if (thermo->has_mask(mask_name))
            stats->add_mask(mask_name);
        else if (microphys->has_mask(mask_name))
            stats->add_mask(mask_name);
        else
        {
            std::string error_message = "Can not calculate mask for \"" + mask_name + "\"";
            throw std::runtime_error(error_message);
        }
    }
}

// Print the status information to the .out file.
template<typename TF>
void Model<TF>::print_status()
{
    double cputime, end;
    static double start;
    static FILE *dnsout = NULL;
    static bool first = true;

    if (first)
    {
        start = master.get_wall_clock_time();

        if (master.get_mpiid() == 0)
        {
            std::string outputname = sim_name + ".out";
            dnsout = std::fopen(outputname.c_str(), "a");
            std::setvbuf(dnsout, NULL, _IOLBF, 1024);
            std::fprintf(dnsout, "%8s %13s %10s %11s %8s %8s %11s %16s %16s %16s\n",
                    "ITER", "TIME", "CPUDT", "DT", "CFL", "DNUM", "DIV", "MOM", "TKE", "MASS");
        }
        first = false;
    }

    if (timeloop->do_check())
    {
        const double time = timeloop->get_time();
        const int iter = timeloop->get_iteration();
        const double dt = timeloop->get_dt();
        end     = master.get_wall_clock_time();
        cputime = end - start;
        start   = end;

        boundary->set_ghost_cells_w(Boundary_w_type::Conservation_type);
        const TF div = pres->check_divergence();
        boundary->set_ghost_cells_w(Boundary_w_type::Normal_type);
        TF mom  = fields->check_momentum();
        TF tke  = fields->check_tke();
        TF mass = fields->check_mass();
        TF cfl  = advec->get_cfl(timeloop->get_dt());
        TF dn   = diff->get_dn(timeloop->get_dt());

        if (master.get_mpiid() == 0)
        {
            std::fprintf(dnsout, "%8d %13.6G %10.4f %11.3E %8.4f %8.4f %11.3E %16.8E %16.8E %16.8E\n",
                    iter, time, cputime, dt, cfl, dn, div, mom, tke, mass);
            std::fflush(dnsout);
        }

        if ( !(std::isfinite(cfl) && std::isfinite(dn) && std::isfinite(mom) && std::isfinite(tke) && std::isfinite(mass)) )
        {
            std::string error_message = "Simulation has non-finite numbers";
            throw std::runtime_error(error_message);
        }
    }
}

template class Model<double>;
template class Model<float>;
