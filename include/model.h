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

// forward declaration to prevent very long building times
class cmaster;
class cinput;
class cgrid;
class cfields;
class cboundary;
class ctimeloop;
class cadvec;
class cdiff;
class cpres;
class cforce;
class cthermo;
class cbuffer;
class cstats;
class ccross;
class cbudget;

class cmodel
{
  public:
    cmodel(cmaster *, cinput *);
    ~cmodel();

    void init();
    void create();
    void load();
    void save();
    void exec();

    // make the pointers public for use in other classes
    // TODO maybe it is safer to create get functions
    cmaster *master;
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
    cstats  *stats;
    ccross  *cross;
    cbudget *budget;

  private:
    // list of masks for statistics
    std::vector<std::string> masklist;

    void printOutputFile(bool);
    void calcstats(std::string);
    void settimestep();
};
#endif
