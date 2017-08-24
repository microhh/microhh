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

// In the constructor all classes are initialized and their input is read.
Model::Model(Master *masterin, Input *inputin)
{
    master = masterin;
    input  = inputin;

    // Initialize the pointers as nullptr.
    grid = nullptr;

    try
    {
        int nerror = 0;

        // Create an instance of the Grid class.
        grid = new Grid(master, input);

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
    delete grid;
}

// In the destructor the deletion of all class instances is triggered.
Model::~Model()
{
    delete_objects();
}

// In the init stage all class individual settings are known and the dynamic arrays are allocated.
void Model::init()
{
    grid->init();
}

// In these functions data necessary to start the model is loaded from disk.
void Model::load()
{
    // First load the grid and time to make their information available.
    grid->load();
}

// In these functions data necessary to start the model is saved to disk.
void Model::save()
{
    // Initialize the grid and the fields from the input data.
    grid->create(input);

    // Save the initialized data to disk for the run mode.
    grid->save();
}

void Model::exec()
{
    master->print_message("Starting time integration\n");
}

void Model::set_time_step()
{
}

// Calculate the statistics for all classes that have a statistics function.
void Model::calc_stats(std::string maskname)
{
}

// Print the status information to the .out file.
void Model::print_status()
{
}
