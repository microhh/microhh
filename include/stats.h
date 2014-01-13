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

#ifndef STATS
#define STATS

#include <netcdfcpp.h>
#include <stats.h>

// forward declarations to reduce compilation time
class cmpi;
class cmodel;
class cgrid;
class cfields;

struct statsvar
{
  NcVar *ncvar;
  double *data;
};
typedef std::map<std::string, statsvar> profmap;

class cstats
{
  public:
    cstats(cmodel *);
    virtual ~cstats();

    virtual int readinifile(cinput *);
    virtual int init(double);
    virtual int create(int);
    virtual unsigned long gettimelim(unsigned long);
    virtual int exec(int, double, unsigned long);

  protected:
    cmodel  *model;
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    double sampletime;
    unsigned long isampletime;
};
#endif

