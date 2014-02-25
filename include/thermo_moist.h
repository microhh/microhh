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

#include "thermo.h"

// forward declarations to speed up build time
class cmaster;
class cgrid;
class cfields;

#define rd 287.04
#define rv 461.5
#define ep rd/rv
#define cp 1005
#define lv 2.5e6
#define rhow    1.e3
#define tmelt   273.15
#define p0 1.e5
#define grav 9.81

#define ex1 2.85611940298507510698e-06
#define ex2 -1.02018879928714644313e-11
#define ex3 5.82999832046362073082e-17
#define ex4 -3.95621945728655163954e-22
#define ex5 2.93898686274077761686e-27
#define ex6 -2.30925409555411170635e-32
#define ex7 1.88513914720731231360e-37

#define at 17.27
#define bt 35.86
#define es0 610.78

#define c0 0.6105851e+03
#define c1 0.4440316e+02
#define c2 0.1430341e+01
#define c3 0.2641412e-01
#define c4 0.2995057e-03
#define c5 0.2031998e-05
#define c6 0.6936113e-08
#define c7 0.2564861e-11
#define c8 -.3704404e-13

class cthermo_moist : public cthermo
{
  public:
    cthermo_moist(cmodel *);
    ~cthermo_moist();
    int readinifile(cinput *);
    int init();
    int create(cinput *);
    int exec();
    int getql(cfield3d *, cfield3d *);

    // functions to retrieve buoyancy properties, to be called from other classes
    int checkthermofield(std::string name);
    int getthermofield(cfield3d *, cfield3d *, std::string name);
    int getbuoyancysurf(cfield3d *);
    int getbuoyancyfluxbot(cfield3d *);
    int getprogvars(std::vector<std::string> *); ///< Retrieve a list of prognostic variables.

  private:
    double ps;
    double thvref;
    double rhos;
    double *pmn;

    bool allocated;
    
    int calcbuoyancytend_2nd(double *, double *, double *, double *, double *, double *, double *, double *);
    int calcbuoyancytend_4th(double *, double *, double *, double *, double *, double *, double *, double *);

    int calcbuoyancy(double *, double *, double *, double *, double *, double *);
    int calcN2(double *, double *, double *, double *); ///< Calculation of the Brunt-Vaissala frequency.

    int calchydropres_2nd(double *, double *, double *, double *, double *);
    int calchydropres_4th(double *, double *, double *, double *, double *);

    int calcqlfield(double *, double *, double *, double *);
    int calcbuoyancybot(double *, double *,
                        double *, double *,
                        double *, double *,
                        double *, double *);
    int calcbuoyancyfluxbot(double *, double *, double *, double *, double *, double *);

    inline double calcql(const double, const double, const double ,const double);
    inline double bu(const double, const double, const double, const double, const double);
    inline double bunoql(const double, const double, const double);
    inline double bufluxnoql(const double, const double, const double, const double, const double);
    inline double exn(const double);
    inline double exn2(const double);
    inline double rslf(const double, const double);
    inline double esl(const double);

    inline double interp2(const double, const double);
    inline double interp4(const double, const double, const double, const double);

    // REFERENCE PROFILES
    double *thref;
    double *pref;
    double *exner;
    // double *rhoref;

    double *threfh;
    double *prefh;
    double *exnerh;
    // double *rhorefh;
};
#endif
