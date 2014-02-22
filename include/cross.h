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

// forward declarations to reduce compilation time
class cmaster;
class cmodel;
class cgrid;
class cfields;

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
    cmaster *master;
    cmodel  *model;
    cgrid   *grid;
    cfields *fields;

    double sampletime;
    unsigned long isampletime;

    std::string swcross;

    std::vector<int> jxz;
    std::vector<int> kxy;

    std::vector<std::string> simple;
    std::vector<std::string> bot;
    std::vector<std::string> fluxbot;
    std::vector<std::string> lngrad;
    std::vector<std::string> path;

    int checkList(std::vector<std::string> *, fieldmap *, std::string crossname);

    int crosssimple (double *, double *, std::string, std::vector<int>, std::vector<int>, int);
    int crossbot    (double *, double *, double *, std::string, int);
    int crossfluxbot(double *, double *, std::string, int);
    int crosslngrad (double *, double *, double *, double *, std::string, std::vector<int>, std::vector<int>, int);
    int crosspath   (double *, double *, double *, std::string, int);
};
#endif

