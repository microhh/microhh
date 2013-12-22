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

#ifndef CROSS
#define CROSS

#include <netcdfcpp.h>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"

// forward declaration
class cmodel;

class ccross
{
  public:
    ccross(cmodel *);
    ~ccross();

    int readinifile(cinput *);
    int init(int);
    unsigned long gettimelim(unsigned long);
    int exec(double, unsigned long, int);

  private:
    cmodel  *model;
    cgrid   *grid;
    cfields *fields;
    cmpi    *mpi;

    double crosstime;
    unsigned long icrosstime;

    std::string swcross;

    std::vector<int> jxz;
    std::vector<int> kxy;

    std::vector<std::string> simple;
    std::vector<std::string> lngrad;

    int crosssimple(double *, double *, std::string, std::vector<int>, std::vector<int>, int);
    int crosslngrad(double *, double *, double *, double *, std::string, std::vector<int>, std::vector<int>, int);
};
#endif

