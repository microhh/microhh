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

#ifndef MASTER
#define MASTER

#ifdef USEMPI
#include <mpi.h>
#endif
#include <string>
#include "input.h"

class Input;

class Master
{
    public:
        Master();
        ~Master();

        void start();
        void init(Input&);

        double get_wall_clock_time();
        bool at_wall_clock_limit();

        // Overload the broadcast function.
        void broadcast(char*, int);
        void broadcast(int*, int);
        void broadcast(bool*, int);
        void broadcast(double*, int);
        void broadcast(float*, int);
        void broadcast(unsigned long*, int);

        // Overload the sum function.
        void sum(int*, int);
        void sum(double*, int);
        void sum(float*, int);

        // Overload the max function.
        void max(double*, int);
        void max(float*, int);

        // Overload the min function.
        void min(double*, int);
        void min(float*, int);

        void print_message(const char *format, ...);
        void print_message(const std::ostringstream&);
        void print_message(const std::string&);

        void print_warning(const char *format, ...);
        void print_warning(const std::ostringstream&);
        void print_warning(const std::string&);

        void print_error  (const char *format, ...);

        int nprocs;
        int npx;
        int npy;
        int mpiid;
        int mpicoordx;
        int mpicoordy;

        #ifdef USEMPI
        int nnorth;
        int nsouth;
        int neast;
        int nwest;

        MPI_Comm commxy;
        MPI_Comm commx;
        MPI_Comm commy;

        MPI_Request* get_request_ptr();
        void wait_all();
        #endif

    private:
        bool initialized;
        bool allocated;

        double wall_clock_start;
        double wall_clock_end;

        #ifdef USEMPI
        MPI_Request* reqs;
        int reqsn;

        int check_error(int);
        #endif
};
#endif
