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

#ifndef THERMO_MOIST
#define THERMO_MOIST

#include "grid.h"
#include "fields.h"
#include "master.h"
#include "thermo.h"
#include <cmath>

#define rd 287.04
#define rv 461.5
#define ep rd/rv
#define cp 1005
#define lv 2.5e6
#define rhow    1.e3
#define tmelt   273.15
#define p0 1.e5
#define grav 9.81

class cthermo_moist : public cthermo
{
  public:
    cthermo_moist(cgrid *, cfields *, cmaster *);
    ~cthermo_moist();
    int readinifile(cinput *);
    int init(cinput *);
    int create();
    int exec();
    int getql(cfield3d *, cfield3d *);

    // functions to retrieve buoyancy properties, to be called from other classes
    int getbuoyancysurf   (cfield3d *);
    int getbuoyancyfluxbot(cfield3d *);
    int getbuoyancy(cfield3d *, cfield3d *);
    // int getbuoyancyh();
    // int getsat();
    // int getsath();

  private:
    double ps;
    double thvs;
    double rhos;
    double *pmn;

    bool allocated;
    
    int calcbuoyancytend_2nd(double *, double *, double *, double *);
    int calcbuoyancytend_4th(double *, double *, double *, double *);

    int calcbuoyancy(double *, double *, double *, double *);

    int calcpres(double *, double *, double *);
    int calcqlfield(double *, double *, double *, double *);
    int calcbuoyancybot(double *, double *,
                        double *, double *,
                        double *, double *);
    int calcbuoyancyfluxbot(double *, double *, double *, double *, double *);

    inline double calcql(const double, const double, const double);
    inline double bu(const double, const double, const double, const double);
    inline double bunoql(const double, const double, const double);
    inline double bufluxnoql(const double, const double, const double, const double, const double);
    inline double exner(const double);
    inline double rslf(const double, const double);
    inline double esl(const double);

    inline double interp2(const double, const double);
    inline double interp4(const double, const double, const double, const double);
};
#endif
