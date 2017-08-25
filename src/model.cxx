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
#include "grid.h"
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
    }
}

// In the constructor all classes are initialized and their input is read.
template<typename TF>
Model<TF>::Model(Master *masterin, Input *inputin, int argc, char *argv[])
{
    master = masterin;
    input  = inputin;

    process_command_line_options(simmode, simname, argc, argv, *master);

    // Initialize the pointers as nullptr.
    grid = nullptr;

    try
    {
        int nerror = 0;

        // Create an instance of the Grid class.
        grid = new Grid<TF>(master, input);

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
    delete grid;
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
    grid->init();
}

// In these functions data necessary to start the model is loaded from disk.
template<typename TF>
void Model<TF>::load()
{
    // First load the grid and time to make their information available.
    grid->load();
}

// In these functions data necessary to start the model is saved to disk.
template<typename TF>
void Model<TF>::save()
{
    // Initialize the grid and the fields from the input data.
    grid->create(input);

    // Save the initialized data to disk for the run mode.
    grid->save();
}

template<typename TF>
void Model<TF>::exec()
{
    master->print_message("Starting time integration\n");
}

template<typename TF>
void Model<TF>::set_time_step()
{
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
// template class Model<float>;
