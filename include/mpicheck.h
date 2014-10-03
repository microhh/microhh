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

#ifndef MPICHECK
#define MPICHECK

#include "grid.h"
#include "field3d.h"
#include "master.h"

class Mpicheck
{
  public:
    Mpicheck(Grid *, cmpi *);
    ~Mpicheck();
    
    int create();
    int checkLayout();
    int checkBoundary();
    int checkTranspose();

  private:
    Grid   *grid;
    cmpi    *mpi;

    Field3d *s;
    Field3d *temp1;
    Field3d *temp2;
};
#endif
