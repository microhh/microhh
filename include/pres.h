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

#ifndef PRES
#define PRES

#include <fftw3.h>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

// forward declaration
class cmodel;

class cpres
{
  public:
    cpres(cmodel *);
    virtual ~cpres();

    virtual int readinifile(cinput *);

    virtual int init();
    virtual int setvalues();
    virtual int exec(double);
    virtual double check();

  protected:
    cmodel  *model;
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;
};
#endif
