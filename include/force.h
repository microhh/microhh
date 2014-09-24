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

#ifndef FORCE
#define FORCE

#include <vector>
#include <string>
#include <map>

// forward declaration
class cmodel;
class cgrid;
class cfields;
class cmaster;
class cinput;

/**
 * Class for the right-hand side terms that contain large-scale forcings
 * This class contains the large-scale pressure forcings, either in flux for or through a
 * geostrophic wind and a coriolis force. Furthermore, a large scale vertical velocity can
 * be imposed that advects the scalars through the domain. Profiles of sources/sinks can be
 * assigned to all scalars.
 */
class cforce
{
  public:
    cforce(cmodel *, cinput *); ///< Constructor of the force class.
    ~cforce();                  ///< Destructor of the force class.
    void init();                ///< Initialize the arrays that contain the profiles.
    void create(cinput *);      ///< Read the profiles of the forces from the input.
    int exec(double);           ///< Add the tendencies belonging to the large-scale processes.
    int settimedep();           ///< Set the time dependent parameters.

    std::vector<std::string> lslist;         ///< List of variables that have large-scale forcings.
    std::map<std::string, double *> lsprofs; ///< Map of profiles with forcings stored by its name.

    // GPU functions and variables
    int prepareDevice();
    int clearDevice();

    std::map<std::string, double *> lsprofs_g; ///< Map of profiles with forcings stored by its name.

  private:
    cmaster *master; ///< Pointer to master class.
    cmodel  *model;  ///< Pointer to model class.
    cgrid   *grid;   ///< Pointer to grid class.
    cfields *fields; ///< Pointer to fields class.

    std::string swlspres; ///< Switch for the large scale pressure force.
    std::string swls;     ///< Switch for large scale scalar tendencies.
    std::string swwls;    ///< Switch for large-scale vertical transport of scalars.

    double uflux; ///< Mean velocity used to enforce constant flux.
    double fc;    ///< Coriolis parameter.

    double *ug;  ///< Pointer to array u-component geostrophic wind.
    double *vg;  ///< Pointer to array v-component geostrophic wind.
    double *wls; ///< Pointer to array large-scale vertical velocity.

    // time dependent variables
    std::string swtimedep;
    std::vector<double> timedeptime;
    std::vector<std::string> timedeplist;
    std::map<std::string, double *> timedepdata;

    int settimedepprofiles(double, double, int, int); ///< Set the time dependent profiles.

    int flux(double * const, const double * const,
             const double * const, const double);  ///< Calculates the pressure force to enforce a constant mass-flux.

    int coriolis_2nd(double * const, double * const,
                     const double * const, const double * const,
                     const double * const, const double * const); ///< Calculates Coriolis force with 2nd-order accuracy.
    int coriolis_4th(double * const, double * const,
                     const double * const, const double * const,
                     const double * const, const double * const); ///< Calculates Coriolis force with 4th-order accuracy.

    int lssource(double * const, const double * const); ///< Applies the large scale scalar tendency.

    int advecwls_2nd(double * const, const double * const,
                     const double * const, const double * const); ///< Calculates the large-scale vertical transport.

    // GPU functions and variables
    double *ug_g;  ///< Pointer to GPU array u-component geostrophic wind.
    double *vg_g;  ///< Pointer to GPU array v-component geostrophic wind.
    double *wls_g; ///< Pointer to GPU array large-scale vertical velocity.
    std::map<std::string, double *> timedepdata_g;

};
#endif
