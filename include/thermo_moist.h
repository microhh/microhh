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

class Master;
class Grid;
class Fields;
class Stats;
struct Mask;

class ThermoMoist : public Thermo
{
  public:
    ThermoMoist(Model *, Input *);
    virtual ~ThermoMoist();

    virtual void init();
    virtual void create(Input *);
    virtual void exec();
    virtual void getMask(Field3d *, Field3d *, Mask *);
    virtual void execStats(Mask *);
    virtual void execCross();
    virtual void execDump();

    // functions to retrieve buoyancy properties, to be called from other classes
    virtual bool checkThermoField(std::string name);
    virtual void getThermoField(Field3d *, Field3d *, std::string name);
    virtual void getBuoyancySurf(Field3d *);
    virtual void getBuoyancyFluxbot(Field3d *);
    virtual void getProgVars(std::vector<std::string> *); ///< Retrieve a list of prognostic variables.

    #ifdef USECUDA
    // GPU functions and variables
    void prepareDevice();
    void clearDevice();
    #endif

  private:

    int swupdatebasestate;

    // cross sections
    std::vector<std::string> crosslist;        ///< List with all crosses from ini file
    std::vector<std::string> allowedcrossvars; ///< List with allowed cross variables
    std::vector<std::string> dumplist;         ///< List with all 3d dumps from the ini file.

    Stats *stats;

    // masks
    int calcMask_ql    (double *, double *, double *, int *, int *, int *, double *);
    int calcMask_qlcore(double *, double *, double *, int *, int *, int *, double *, double *, double *);

    int calcBuoyancyTend_2nd(double *, double *, double *, double *, double *, double *, double *, double *);
    int calcBuoyancyTend_4th(double *, double *, double *, double *, double *, double *, double *, double *);

    int calcBuoyancy(double *, double *, double *, double *, double *, double *);
    int calcN2(double *, double *, double *, double *); ///< Calculation of the Brunt-Vaissala frequency.
    int calcBaseState(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);

    int calcLiquidWater(double *, double *, double *, double *);
    int calcBuoyancyBot(double *, double *,
                        double *, double *,
                        double *, double *,
                        double *, double *);
    int calcBuoyancyFluxBot(double *, double *, double *, double *, double *, double *);

    inline double satAdjust(const double, const double, const double ,const double);
    inline double buoyancy(const double, const double, const double, const double, const double, const double);
    inline double buoyancyNoql(const double, const double, const double);
    inline double buoyancyFluxNoql(const double, const double, const double, const double, const double);
    inline double exner(const double);
    inline double exn2(const double);
    inline double qsat(const double, const double);
    inline double esat(const double);

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

    // GPU functions and variables
    double *thvref_g; 
    double *thvrefh_g;
    double *exnref_g;
    double *exnrefh_g;
    double *pref_g;
    double *prefh_g;
};
#endif
