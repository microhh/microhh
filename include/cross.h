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
    ccross(cmodel *, cinput *);
    ~ccross();

    int init(int);
    int create();
    unsigned long gettimelim(unsigned long);
    //int exec(double, unsigned long, int);

    std::string swcross;
    int docross();

    int crosssimple (double *, double *, std::string);
    int crosslngrad (double *, double *, double *, double *, std::string);
    int crossplane  (double *, double *, std::string);
    int crosspath   (double *, double *, double *, std::string);

  private:
    cmaster *master;
    cmodel  *model;
    cgrid   *grid;
    cfields *fields;

    double sampletime;
    unsigned long isampletime;

    std::vector<int> jxz;   ///< Index of nearest full y position of xz input
    std::vector<int> kxy;   ///< Index of nearest full height level of xy input
    std::vector<int> jxzh;  ///< Index of nearest half y position of xz input
    std::vector<int> kxyh;  ///< Index of nearest half height level of xy input
    std::vector<double> xz; ///< Y-position [m] xz cross from ini file
    std::vector<double> xy; ///< Z-position [m] xy cross from ini file

    std::vector<std::string> simple;
    std::vector<std::string> bot;
    std::vector<std::string> fluxbot;
    std::vector<std::string> lngrad;
    std::vector<std::string> path;

    int checkList(std::vector<std::string> *, fieldmap *, std::string crossname);
    int checkSave(int, char *);
};
#endif

