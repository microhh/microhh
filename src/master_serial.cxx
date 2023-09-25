/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

#ifndef USEMPI

#include <sys/time.h>
#include <stdexcept>

#include "grid.h"
#include "defines.h"
#include "master.h"

Master::Master()
{
    initialized = false;
    allocated   = false;
}

Master::~Master()
{
    print_message("Finished run on %d processes\n", md.nprocs);
}

void Master::start()
{
    initialized = true;

    // Set the rank of the only process to 0.
    md.mpiid = 0;
    // Set the number of processes to 1.
    md.nprocs = 1;

    // Set the wall clock time at start.
    wall_clock_start = get_wall_clock_time();

    print_message("Starting run on %d processes\n", md.nprocs);

}

void Master::init(Input& input)
{
    md.npx = input.get_item<int>("master", "npx", "", 1);
    md.npy = input.get_item<int>("master", "npy", "", 1);

    // Get the wall clock limit with a default value of 1E8 hours, which will be never hit
    double wall_clock_limit = input.get_item<double>("master", "wallclocklimit", "", 1E8);

    wall_clock_end = wall_clock_start + 3600.*wall_clock_limit;

    if (md.nprocs != md.npx*md.npy)
    {
        std::string msg = "npx*npy = " + std::to_string(md.npy) + "*" + std::to_string(md.npy) + " has to be equal to 1*1 in serial mode";
        throw std::runtime_error(msg);
    }

    // set the coordinates to 0
    md.mpicoordx = 0;
    md.mpicoordy = 0;

    allocated = true;
}

double Master::get_wall_clock_time()
{
    timeval timestruct;
    gettimeofday(&timestruct, NULL);
    return (double)timestruct.tv_sec + (double)timestruct.tv_usec*1.e-6;
}

// void Master::wait_all() {}

// All broadcasts return directly, because there is nothing to broadcast.
void Master::broadcast(char* data, int datasize, int mpiid_to_send) {}
void Master::broadcast(signed char* data, int datasize, int mpiid_to_send) {}
void Master::broadcast(int* data, int datasize, int mpiid_to_send) {}
void Master::broadcast(unsigned long* data, int datasize, int mpiid_to_send) {}
void Master::broadcast(double* data, int datasize, int mpiid_to_send) {}
void Master::broadcast(float* data, int datasize, int mpiid_to_send) {}
void Master::sum(int* var, int datasize) {}
void Master::sum(double* var, int datasize) {}
void Master::sum(float* var, int datasize) {}
void Master::max(double* var, int datasize) {}
void Master::max(float* var, int datasize) {}
void Master::min(double* var, int datasize) {}
void Master::min(float* var, int datasize) {}
#endif
