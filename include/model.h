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

#ifndef MODEL
#define MODEL

#include <string>

#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "boundary.h"
#include "advec.h"
#include "diff.h"
#include "force.h"
#include "thermo.h"
#include "thermo_moist.h"
#include "pres.h"
#include "buffer.h"
#include "timeloop.h"
#include "stats.h"
#include "cross.h"

class cmodel
{
  public:
    cmodel(cmpi *, cinput *);
    ~cmodel();
    int readinifile();
    int init();
    int create();
    int load();
    int save();
    int exec();

    // make the pointers public for use in other classes
    // maybe safer to create get functions
    cmpi    *mpi;
    cinput  *input;
    cgrid   *grid;
    cfields *fields;

    // model operators
    cboundary *boundary;
    ctimeloop *timeloop;
    cadvec    *advec;
    cdiff     *diff;
    cpres     *pres;  
    cforce    *force;   
    cthermo   *thermo;
    cbuffer   *buffer;

    // postprocessing modules
    cstats    *stats;
    ccross    *cross;

  private:
    // switches for included schemes
    std::string swadvec;
    std::string swdiff;
    std::string swpres;
    std::string swboundary;
    std::string swthermo;

    std::string swstats;

    int outputfile(bool);
    int settimestep();
};
#endif
