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

// forward declarations to reduce compilation time
class cmaster;
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
    ~cstats();

    int readinifile(cinput *);
    int init(double);
    int create(int);
    unsigned long gettimelim(unsigned long);
    int exec(int, double, unsigned long);
    int dostats();

    // interface functions
    profmap profs;
    int addprof(std::string, std::string, std::string, std::string);

    int calcmean    (double *, double *, double);
    int calcmoment  (double *, double *, double *, double, int);
    int calcdiff_2nd(double *, double *, double *, double *, double *, double *, double);
    int calcdiff_4th(double *, double *, double *, double);
    int calcgrad_2nd(double *, double *, double *);
    int calcgrad_4th(double *, double *, double *);
    int calcflux_2nd(double *, double *, double *, double *, int, int);
    int calcflux_4th(double *, double *, double *, double *, int, int);
    int addfluxes   (double *, double *, double *);
    int calccount   (double* data, double* prof, double threshold);

  private:
    bool allocated;
    bool initialized;

    NcFile *dataFile;
    NcDim  *z_dim, *zh_dim, *t_dim;
    NcVar  *z_var, *zh_var, *t_var, *iter_var;

    double *umodel, *vmodel;

    int nstats;

  protected:
    cmodel  *model;
    cgrid   *grid;
    cfields *fields;
    cmaster *master;

    double sampletime;
    unsigned long isampletime;
};
#endif
