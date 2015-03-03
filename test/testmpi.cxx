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

#include <cstdio>
#include "grid.h"
#include "mpiinterface.h"
#include "mpicheck.h"

int main(int argc, char *argv[])
{
  // set the name of the simulation
  std::string simname("microhh");
  if(argc > 1)
    simname = argv[1];

  // start the MPI
  cmpi      mpi;
  mpi.startup();

  // create the instances of the objects
  cinput    input   (&mpi);
  cgrid     grid    (&mpi);
  cmpicheck mpicheck(&grid, &mpi);

  if(input.readinifile(simname))
    return 1;
  if(grid.readinifile(&input))
    return 1;
  if(mpi.readinifile(&input))
    return 1;

  // initialize the objects, allocate the required memory
  if(mpi.init())
    return 1;
  if(grid.init())
    return 1;

  // check the layout
  mpicheck.checkLayout();

  // fill the fields with data
  if(grid.create(&input))
    return 1;

  // fill the fields with the test data
  mpicheck.create();

  // trigger the boundary conditions
  mpicheck.checkBoundary();

  mpicheck.checkTranspose();

  return 0;
}

