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

#ifndef FORCE
#define FORCE

// forward declaration
class cmodel;
class cgrid;
class cfields;
class cmaster;

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
    int init();                 ///< Initialize the arrays that contain the profiles.
    int create(cinput *);       ///< Read the profiles of the forces from the input.
    int exec(double);           ///< Add the tendencies belonging to the large-scale processes.
    int settimedep();           ///< Set the time dependent parameters.

    std::vector<std::string> lslist;         ///< List of variables that have large-scale forcings.
    std::map<std::string, double *> lsprofs; ///< Map of profiles with forcings stored by its name.

  private:
    cmaster *master; ///< Pointer to master class.
    cmodel  *model;  ///< Pointer to model class.
    cgrid   *grid;   ///< Pointer to grid class.
    cfields *fields; ///< Pointer to fields class.

    bool allocated; ///< Boolean flag to indicate allocation of arrays.

    std::string swlspres; ///< Switch for the large scale pressure force.
    std::string swls;     ///< Switch for large scale scalar tendencies.
    std::string swwls;    ///< Switch for large-scale vertical transport of scalars.

    double uflux; ///< Mean velocity used to enforce constant flux.
    double fc;    ///< Coriolis parameter.

    double *ug; ///< Pointer to array u-component geostrophic wind.
    double *vg; ///< Pointer to array v-component geostrophic wind.

    double *wls; ///< Pointer to array large-scale vertical velocity.

    // time dependent variables
    std::string swtimedep;
    std::vector<double> timedeptime;
    std::vector<std::string> timedeplist;
    std::map<std::string, double *> timedepdata;

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

    inline double interp2(const double, const double); ///< 2nd order interpolation function.
};
#endif

