/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
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

#ifdef USEMPI

#include <mpi.h>
#include <stdexcept>

#include "grid.h"
#include "defines.h"
#include "master.h"

Master::Master()
{
    initialized = false;
    allocated   = false;

    // set the mpiid, to ensure that errors can be written if MPI init fails
    md.mpiid = 0;
}

Master::~Master()
{
    if (allocated)
    {
        delete[] reqs;
        MPI_Comm_free(&md.commxy);
        MPI_Comm_free(&md.commx);
        MPI_Comm_free(&md.commy);
    }

    print_message("Finished run on %d processes\n", md.nprocs);

    if (initialized)
        MPI_Finalize();
}

void Master::start()
{
    // initialize the MPI
    int n = MPI_Init(NULL, NULL);
    if (check_error(n))
        throw std::runtime_error("MPI init error");

    wall_clock_start = get_wall_clock_time();

    initialized = true;

    // get the rank of the current process
    n = MPI_Comm_rank(MPI_COMM_WORLD, &md.mpiid);
    if (check_error(n))
        throw std::runtime_error("MPI init error");

    // get the total number of processors
    n = MPI_Comm_size(MPI_COMM_WORLD, &md.nprocs);
    if (check_error(n))
        throw std::runtime_error("MPI init error");

    // store a temporary copy of COMM_WORLD in commxy
    n = MPI_Comm_dup(MPI_COMM_WORLD, &md.commxy);
    if (check_error(n))
        throw std::runtime_error("MPI init error");

    print_message("Starting run on %d processes\n", md.nprocs);
}

void Master::init(Input& input)
{
    md.npx = input.get_item<int>("master", "npx", "", 1);
    md.npy = input.get_item<int>("master", "npy", "", 1);

    // Get the wall clock limit with a default value of 1E8 hours, which will be never hit.
    double wall_clock_limit = input.get_item<double>("master", "wallclocklimit", "", 1E8);

    wall_clock_end = wall_clock_start + 3600.*wall_clock_limit;

    if (md.nprocs != md.npx*md.npy)
    {
        std::string msg = "nprocs = " + std::to_string(md.nprocs) + " does not equal npx*npy = " + std::to_string(md.npx) + "*" + std::to_string(md.npy);
        throw std::runtime_error(msg);
    }

    int n;
    int dims    [2] = {md.npy, md.npx};
    int periodic[2] = {true, true};

    // define the dimensions of the 2-D grid layout
    n = MPI_Dims_create(md.nprocs, 2, dims);
    if (check_error(n))
        throw std::runtime_error("MPI init error");

    // create a 2-D grid communicator that is optimized for grid to grid transfer
    // first, free our temporary copy of COMM_WORLD
    n = MPI_Comm_free(&md.commxy);
    if (check_error(n))
        throw std::runtime_error("MPI init error");

    // for now, do not reorder processes, blizzard gives large performance loss
    n = MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, false, &md.commxy);
    if (check_error(n))
        throw std::runtime_error("MPI init error");

    n = MPI_Comm_rank(md.commxy, &md.mpiid);
    if (check_error(n))
        throw std::runtime_error("MPI init error");

    // retrieve the x- and y-coordinates in the 2-D grid for each process
    int mpicoords[2];
    n = MPI_Cart_coords(md.commxy, md.mpiid, 2, mpicoords);
    if (check_error(n))
        throw std::runtime_error("MPI init error");

    md.mpicoordx = mpicoords[1];
    md.mpicoordy = mpicoords[0];

    int dimx[2] = {false, true };
    int dimy[2] = {true , false};

    n = MPI_Cart_sub(md.commxy, dimx, &md.commx);
    if (check_error(n))
        throw std::runtime_error("MPI init error");

    n = MPI_Cart_sub(md.commxy, dimy, &md.commy);
    if (check_error(n))
        throw std::runtime_error("MPI init error");

    // find out who are the neighbors of this process to facilitate the communication routines
    n = MPI_Cart_shift(md.commxy, 1, 1, &md.nwest , &md.neast );
    if (check_error(n))
        throw std::runtime_error("MPI init error");

    n = MPI_Cart_shift(md.commxy, 0, 1, &md.nsouth, &md.nnorth);
    if (check_error(n))
        throw std::runtime_error("MPI init error");

    // create the requests arrays for the nonblocking sends
    int npmax;
    npmax = std::max(md.npx, md.npy);

    // have at least as many communicators as prognostic variables
    npmax = std::max(npmax, 8*4);
    reqs  = new MPI_Request[npmax*2];
    reqsn = 0;

    allocated = true;
}

double Master::get_wall_clock_time()
{
    return MPI_Wtime();
}

int Master::check_error(int n)
{
    char errbuffer[MPI_MAX_ERROR_STRING];
    int errlen;

    if (n != MPI_SUCCESS)
    {
        MPI_Error_string(n, errbuffer, &errlen);
        print_message("MPI: %s\n", errbuffer);
        return 1;
    }

    return 0;
}

MPI_Request* Master::get_request_ptr()
{
    // Get a request out of the pointer and increment the counter.
    MPI_Request* req = &reqs[reqsn];
    ++reqsn;
    return req;
}

void Master::wait_all()
{
    // Wait for MPI processes and reset the number of pending requests.
    MPI_Waitall(reqsn, reqs, MPI_STATUSES_IGNORE);
    reqsn = 0;
}

// CvH obsolete: do all broadcasts over the MPI_COMM_WORLD, to avoid complications in the input file reading
void Master::broadcast(char *data, int datasize, int mpiid_to_send)
{
    MPI_Bcast(data, datasize, MPI_CHAR, mpiid_to_send, md.commxy);
}

void Master::broadcast(signed char *data, int datasize, int mpiid_to_send)
{
    MPI_Bcast(data, datasize, MPI_SIGNED_CHAR, mpiid_to_send, md.commxy);
}

void Master::broadcast(int* data, int datasize, int mpiid_to_send)
{
    MPI_Bcast(data, datasize, MPI_INT, mpiid_to_send, md.commxy);
}

void Master::broadcast(unsigned long* data, int datasize, int mpiid_to_send)
{
    MPI_Bcast(data, datasize, MPI_UNSIGNED_LONG, mpiid_to_send, md.commxy);
}

void Master::broadcast(double* data, int datasize, int mpiid_to_send)
{
    MPI_Bcast(data, datasize, MPI_DOUBLE, mpiid_to_send, md.commxy);
}

void Master::broadcast(float* data, int datasize, int mpiid_to_send)
{
    MPI_Bcast(data, datasize, MPI_FLOAT, mpiid_to_send, md.commxy);
}

void Master::sum(int* var, int datasize)
{
    MPI_Allreduce(MPI_IN_PLACE, var, datasize, MPI_INT, MPI_SUM, md.commxy);
}

void Master::sum(double* var, int datasize)
{
    MPI_Allreduce(MPI_IN_PLACE, var, datasize, MPI_DOUBLE, MPI_SUM, md.commxy);
}

void Master::sum(float* var, int datasize)
{
    MPI_Allreduce(MPI_IN_PLACE, var, datasize, MPI_FLOAT, MPI_SUM, md.commxy);
}

void Master::max(double* var, int datasize)
{
    MPI_Allreduce(MPI_IN_PLACE, var, datasize, MPI_DOUBLE, MPI_MAX, md.commxy);
}

void Master::max(float* var, int datasize)
{
    MPI_Allreduce(MPI_IN_PLACE, var, datasize, MPI_FLOAT, MPI_MAX, md.commxy);
}

void Master::min(double* var, int datasize)
{
    MPI_Allreduce(MPI_IN_PLACE, var, datasize, MPI_DOUBLE, MPI_MIN, md.commxy);
}

void Master::min(float* var, int datasize)
{
    MPI_Allreduce(MPI_IN_PLACE, var, datasize, MPI_FLOAT, MPI_MIN, md.commxy);
}
#endif
