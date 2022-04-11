/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
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

#ifndef FORCE_H
#define FORCE_H

#include <vector>
#include <string>
#include <map>

class Master;
class Input;
class Netcdf_handle;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Field3d_operators;
template<typename> class Timedep;
template<typename> class Stats;
template<typename> class Thermo;

/**
 * Class for the right-hand side terms that contain large-scale forcings
 * This class contains the large-scale pressure forcings, either in flux for or through a
 * geostrophic wind and a coriolis force. Furthermore, a large scale vertical velocity can
 * be imposed that advects the scalars through the domain. Profiles of sources/sinks can be
 * assigned to all scalars.
 */

enum class Large_scale_pressure_type {Disabled, Fixed_flux, Geo_wind, Pressure_gradient};
enum class Large_scale_tendency_type {Disabled, Enabled};
enum class Large_scale_subsidence_type {Disabled, Mean_field, Local_field};
enum class Nudging_type {Disabled, Enabled};

template<typename TF>
class Force
{
    public:
        Force(Master&, Grid<TF>&, Fields<TF>&, Input&); ///< Constructor of the force class.
        ~Force();                                       ///< Destructor of the force class.

        void init();           ///< Initialize the arrays that contain the profiles.
        void create(Input&, Netcdf_handle&, Stats<TF>&);   ///< Read the profiles of the forces from the input.
        void exec(double, Thermo<TF>&, Stats<TF>&);     ///< Add the tendencies belonging to the large-scale processes.

        void update_time_dependent(Timeloop<TF>&); ///< Update the time dependent parameters.

        std::vector<std::string> lslist;        ///< List of variables that have large-scale forcings.
        std::map<std::string, std::vector<TF>> lsprofs; ///< Map of profiles with forcings stored by its name.

        std::vector<std::string> nudgelist;        ///< List of variables that are nudged to a provided profile
        std::vector<std::string> scalednudgelist;        ///< List of variables that are nudged to a provided profile
        std::map<std::string, std::vector<TF>> nudgeprofs; ///< Map of nudge profiles stored by its name.

        // GPU functions and variables
        void prepare_device();
        void clear_device();

        std::map<std::string, TF*> lsprofs_g;    ///< Map of profiles with forcings stored by its name.
        std::map<std::string, TF*> nudgeprofs_g; ///< Map of nudging profiles stored by its name.

        // Accessor functions
        Large_scale_pressure_type get_switch_lspres() { return swlspres; }
        TF get_coriolis_parameter() const { return fc; }


    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Field3d_operators<TF> field3d_operators;

        // Internal switches for various forcings
        Large_scale_pressure_type swlspres;
        Large_scale_tendency_type swls;
        Large_scale_subsidence_type swwls;
        Nudging_type swnudge;
        bool swwls_mom;

        TF uflux; ///< Mean velocity used to enforce constant flux.
        TF dpdx;  ///< Large-scale pressure gradient
        TF fc;    ///< Coriolis parameter.

        std::vector<TF> ug;  ///< Pointer to array u-component geostrophic wind.
        std::vector<TF> vg;  ///< Pointer to array v-component geostrophic wind.
        std::vector<TF> wls; ///< Pointer to array large-scale vertical velocity.

        std::vector<TF> nudge_factor;  ///< Height varying nudging factor (1/s)

        std::map<std::string, Timedep<TF>*> tdep_ls;
        std::map<std::string, Timedep<TF>*> tdep_geo;
        std::map<std::string, Timedep<TF>*> tdep_nudge;
        std::unique_ptr<Timedep<TF>> tdep_wls;

        // GPU functions and variables
        TF* ug_g;  ///< Pointer to GPU array u-component geostrophic wind.
        TF* vg_g;  ///< Pointer to GPU array v-component geostrophic wind.
        TF* wls_g; ///< Pointer to GPU array large-scale vertical velocity.
        TF* nudge_factor_g; ///< Pointer to GPU array nudge factor.

        const std::string tend_name_pres      = "lspres";
        const std::string tend_longname_pres  = "Large Scale Pressure";
        const std::string tend_name_cor       = "cor";
        const std::string tend_longname_cor   = "Coriolis";
        const std::string tend_name_ls        = "ls";
        const std::string tend_longname_ls    = "Large Scale";
        const std::string tend_name_nudge     = "nudge";
        const std::string tend_longname_nudge = "Nudging";
        const std::string tend_name_subs      = "subs";
        const std::string tend_longname_subs  = "Subsidence";
};
#endif
