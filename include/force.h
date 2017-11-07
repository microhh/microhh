/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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

class Model;
class Grid;
class Fields;
class Master;
class Input;

/**
 * Class for the right-hand side terms that contain large-scale forcings
 * This class contains the large-scale pressure forcings, either in flux for or through a
 * geostrophic wind and a coriolis force. Furthermore, a large scale vertical velocity can
 * be imposed that advects the scalars through the domain. Profiles of sources/sinks can be
 * assigned to all scalars.
 */
class Force
{
    public:
        Force(Model*, Input*); ///< Constructor of the force class.
        ~Force();              ///< Destructor of the force class.

        void init();           ///< Initialize the arrays that contain the profiles.
        void create(Input*);   ///< Read the profiles of the forces from the input.
        void exec(double);     ///< Add the tendencies belonging to the large-scale processes.

        void update_time_dependent(); ///< Update the time dependent parameters.

        std::vector<std::string> lslist;        ///< List of variables that have large-scale forcings.
        std::map<std::string, double*> lsprofs; ///< Map of profiles with forcings stored by its name.

        std::vector<std::string> nudgelist;        ///< List of variables that are nudged to a provided profile
        std::map<std::string, double*> nudgeprofs; ///< Map of nudge profiles stored by its name.

        // GPU functions and variables
        void prepare_device();
        void clear_device();

        std::map<std::string, double*> lsprofs_g; ///< Map of profiles with forcings stored by its name.

        // Accessor functions
        std::string get_switch_lspres()      { return swlspres; }
        double      get_coriolis_parameter() { return fc;       }

    private:
        Master* master; ///< Pointer to master class.
        Model*  model;  ///< Pointer to model class.
        Grid*   grid;   ///< Pointer to grid class.
        Fields* fields; ///< Pointer to fields class.

        std::string swlspres; ///< Switch for the large scale pressure force.
        std::string swls;     ///< Switch for large scale scalar tendencies.
        std::string swwls;    ///< Switch for large-scale vertical transport of scalars.
        std::string swnudge;  ///< Switch for nudging to reference profiles

        double uflux; ///< Mean velocity used to enforce constant flux.
        double fc;    ///< Coriolis parameter.

        double* ug;  ///< Pointer to array u-component geostrophic wind.
        double* vg;  ///< Pointer to array v-component geostrophic wind.
        double* wls; ///< Pointer to array large-scale vertical velocity.

        double nudge_tau;                  ///< Nudging time scale (s)
        std::vector<double> nudge_factor;  ///< Height varying nudging factor (-)

        // Time dependence geostrophic wind
        std::string swtimedep_geo;
        std::map<std::string, std::vector<double>> timedeptime_geo;
        std::map<std::string, double*> timedepdata_geo;

        // Time dependence large-scale forcings
        std::string swtimedep_ls;
        std::vector<std::string> timedeplist_ls;
        std::map<std::string, std::vector<double>> timedeptime_ls;
        std::map<std::string, double*> timedepdata_ls;

        // Time dependence nudging profiles
        std::string swtimedep_nudge;
        std::vector<std::string> timedeplist_nudge;
        std::map<std::string, std::vector<double>> timedeptime_nudge;
        std::map<std::string, double*> timedepdata_nudge;

        // Time dependence subsidence
        std::string swtimedep_wls;
        std::vector<double> timedeptime_wls;
        double* timedepdata_wls;

        int create_timedep(std::map<std::string, double*>&, std::map<std::string, std::vector<double>>&,
                            std::vector<std::string>&, std::vector<std::string>, std::string);

        void update_time_dependent_profs(std::map<std::string, double*>&, std::map<std::string, double*>,
                                         std::map<std::string, std::vector<double>> times, std::string);

        void update_time_dependent_prof(double* const, const double* const, std::vector<double>);

        void calc_flux(double* const, const double* const,
                       const double* const, const double);  ///< Calculates the pressure force to enforce a constant mass-flux.

        void calc_coriolis_2nd(double* const, double* const,
                               const double* const, const double* const,
                               const double* const, const double* const); ///< Calculates Coriolis force with 2nd-order accuracy.

        void calc_coriolis_4th(double* const, double* const,
                               const double* const, const double* const,
                               const double* const, const double* const); ///< Calculates Coriolis force with 4th-order accuracy.

        void calc_large_scale_source(double* const, const double* const); ///< Applies the large scale scalar tendency.

        void calc_nudging_tendency(double* const, const double* const,
                                   const double* const, const double* const); ///< Calculate nudging tendency.

        void advec_wls_2nd(double* const, const double* const,
                           const double* const, const double* const); ///< Calculates the large-scale vertical transport.

        // GPU functions and variables
        double* ug_g;  ///< Pointer to GPU array u-component geostrophic wind.
        double* vg_g;  ///< Pointer to GPU array v-component geostrophic wind.
        double* wls_g; ///< Pointer to GPU array large-scale vertical velocity.
        std::map<std::string, double*> timedepdata_ls_g;

};
#endif
