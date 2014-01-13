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

#ifndef MPICHECK
#define MPICHECK

#include "grid.h"
#include "field3d.h"
#include "master.h"

class cmpicheck
{
  public:
    cmpicheck(cgrid *, cmpi *);
    ~cmpicheck();
    
    int create();
    int checkLayout();
    int checkBoundary();
    int checkTranspose();

  private:
    cgrid   *grid;
    cmpi    *mpi;

    cfield3d *s;
    cfield3d *temp1;
    cfield3d *temp2;
};
#endif
