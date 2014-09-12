/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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
class cstats;
struct mask;

class cthermo_moist : public cthermo
{
  public:
    cthermo_moist(cmodel *, cinput *);
    ~cthermo_moist();

    void init();
    int create(cinput *);
    int exec();
    int getmask(cfield3d *, cfield3d *, mask *);
    int execstats(mask *);
    int execcross();

    // functions to retrieve buoyancy properties, to be called from other classes
    int checkthermofield(std::string name);
    int getthermofield(cfield3d *, cfield3d *, std::string name);
    int getbuoyancysurf(cfield3d *);
    int getbuoyancyfluxbot(cfield3d *);
    int getprogvars(std::vector<std::string> *); ///< Retrieve a list of prognostic variables.

  private:

    int swupdatebasestate;

    // cross sections
    std::vector<std::string> crosslist;        // List with all crosses from ini file
    std::vector<std::string> allowedcrossvars; // List with allowed cross variables

    cstats *stats;

    // masks
    int calcmaskql    (double *, double *, double *, int *, int *, int *, double *);
    int calcmaskqlcore(double *, double *, double *, int *, int *, int *, double *, double *, double *);

    int calcbuoyancytend_2nd(double *, double *, double *, double *, double *, double *, double *, double *);
    int calcbuoyancytend_4th(double *, double *, double *, double *, double *, double *, double *, double *);

    int calcbuoyancy(double *, double *, double *, double *, double *, double *);
    int calcN2(double *, double *, double *, double *); ///< Calculation of the Brunt-Vaissala frequency.
    int calcbasestate(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

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

    double pbot;
    double thvref0; ///< Reference virtual potential temperature in case of Boussinesq

    // REFERENCE PROFILES
    double *thl0;    // Initial thl profile 
    double *qt0;     // Initial qt profile
    double *thvref; 
    double *thvrefh;
    double *exnref;
    double *exnrefh;
    double *pref;
    double *prefh;
};
#endif
