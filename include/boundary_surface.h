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

// forward declaration
class Model;
class Stats;
struct mask;

class Boundary_surface : public Boundary
{
  public:
    Boundary_surface(Model *, Input *);
    ~Boundary_surface();

    void init(Input *);
    void create(Input *);
    void setValues();

    int exec();

    void execCross(); ///< Execute cross sections of surface
    int execStats(mask *); ///< Execute statistics of surface

    void save(int);
    void load(int);

    // make these variables public for out-of-class usage
    double *obuk;
    double *ustar;

    double z0m;
    double z0h;

#ifdef USECUDA
    // GPU functions and variables
    int prepareDevice();
    int clearDevice();
    int forwardDevice();  // TMP BVS
    int backwardDevice(); // TMP BVS 

    double *obuk_g;
    double *ustar_g;
#endif

  private:
    // cross sections
    std::vector<std::string> crosslist;        // List with all crosses from ini file
    std::vector<std::string> allowedcrossvars; // List with allowed cross variables

    // surface scheme
    int bcvalues();
    int stability(double *, double *, double *,
                  double *, double *, double *,
                  double *, double *, double *,
                  double *, double *);
    int stability_neutral(double *, double *,
                          double *, double *,
                          double *, double *,
                          double *, double *);
    int surfm(double *, double *,
              double *, double *, double *, double *,
              double *, double *, double *, double *,
              double, int);
    int surfs(double *, double *, double *,
              double *, double *, double *,
              double, int);
    double calcobuk_noslip_flux     (double, double, double, double);
    double calcobuk_noslip_dirichlet(double, double, double, double);
    inline double fm(double, double, double);
    inline double fh(double, double, double);
    inline double psim(double);
    inline double psih(double);
    inline double phim(double);
    inline double phih(double);
    double ustarin;

    Stats *stats;

    typedef std::map<std::string, int> bcbotmap;
    // int surfmbcbot;
    // bcbotmap surfsbcbot;

    int thermobc;
};
#endif
