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

#ifndef MASTER
#define MASTER

#ifdef PARALLEL
#include <mpi.h>
#endif
#include <string>
#include "input.h"

class cinput;

class cmaster
{
  public:
    cmaster();
    ~cmaster();

    void startup(int, char**);
    void init(cinput *);

    double gettime();
    int waitall();

    // overload the broadcast function
    int broadcast(char *, int);
    int broadcast(int *, int);
    int broadcast(double *, int);
    int broadcast(unsigned long *, int);

    // overload the sum function
    int sum(int *, int);
    int sum(double *, int);

    // overload the max function
    int max(double *, int);

    // overload the min function
    int min(double *, int);

    void printMessage(const char *format, ...);
    void printWarning(const char *format, ...);
    void printError  (const char *format, ...);

    std::string mode;
    std::string simname;

    int nprocs;
    int npx;
    int npy;
    int mpiid;
    int mpicoordx;
    int mpicoordy;

#ifdef PARALLEL
    int nnorth;
    int nsouth;
    int neast;
    int nwest;

    MPI_Comm commxy;
    MPI_Comm commx;
    MPI_Comm commy;

    MPI_Request *reqs;
    int reqsn;
#endif

  private:
    bool initialized;
    bool allocated;

#ifdef PARALLEL
    int checkerror(int);
#endif
};
#endif
