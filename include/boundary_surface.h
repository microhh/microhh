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

#ifndef BOUNDARY_SURFACE
#define BOUNDARY_SURFACE

#include "boundary.h"
#include "stats.h"

class Model;
class Input;
class Stats;
struct Mask;

class BoundarySurface : public Boundary
{
  public:
    BoundarySurface(Model *, Input *);
    ~BoundarySurface();

    void init(Input *);
    void create(Input *);
    void setValues();

    void execStats(Mask *); ///< Execute statistics of surface
    void execCross();       ///< Execute cross sections of surface

    void save(int);
    void load(int);

    // Make these variables public for out-of-class usage.
    double *obuk;
    double *ustar;

    double z0m;
    double z0h;

    #ifdef USECUDA
    // GPU functions and variables
    void prepareDevice();
    void clearDevice();
    void forwardDevice();  // TMP BVS
    void backwardDevice(); // TMP BVS 

    double *obuk_g;
    double *ustar_g;
    #endif

  private:
    // cross sections
    std::vector<std::string> crosslist;        // List with all crosses from ini file
    std::vector<std::string> allowedcrossvars; // List with allowed cross variables

    // surface scheme
    void updateBcs();
    void stability(double *, double *, double *,
                   double *, double *, double *,
                   double *, double *, double *,
                   double *, double *);
    void stabilityNeutral(double *, double *,
                          double *, double *,
                          double *, double *,
                          double *, double *);
    void surfm(double *, double *,
               double *, double *, double *, double *,
               double *, double *, double *, double *,
               double, int);
    void surfs(double *, double *, double *,
               double *, double *, double *,
               double, int);

    double calcObukNoslipFlux     (double, double, double, double, double);
    double calcObukNoslipDirichlet(double, double, double, double);

    double ustarin;

    Stats *stats;

    double* zL_sl;
    double* f_sl;

    // typedef std::map<std::string, int> bcbotmap;
    // int surfmbcbot;
    // bcbotmap surfsbcbot;

    int thermobc;
};
#endif
