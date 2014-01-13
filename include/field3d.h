/*
 * MicroHH
 * Copyright (c) 2011-2013 Chiel van Heerwaarden
 * Copyright (c) 2011-2013 Thijs Heus
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

#ifndef FIELD3D
#define FIELD3D

#include <string>

// forward declarations to reduce compilation time
class cmpi;
class cgrid;

class cfield3d
{
  public:
    // functions
    cfield3d(cgrid *, cmpi *, std::string);
    ~cfield3d();

    int init();
    // int boundary_bottop(int);
    // int boundary_cyclic();
    // int save(int, double *, double *);
    // int load(int, double *, double *);
    int checkfornan();

    // variables
    double *data;
    double *databot;
    double *datatop;
    double *datamean;
    double *datagradbot;
    double *datagradtop;
    double *datafluxbot;
    double *datafluxtop;
    std::string name;
    double visc;

  private:
    cgrid *grid;
    cmpi  *mpi;
    bool allocated;
};
#endif

