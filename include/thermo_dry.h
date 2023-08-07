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

#ifndef THERMO_DRY_H
#define THERMO_DRY_H

#include <stdexcept>

#include "boundary_cyclic.h"
#include "timedep.h"
#include "thermo.h"

class Master;
class Input;
class Netcdf_handle;
template<typename> class Grid;
template<typename> class Stats;
template<typename> class Diff;
template<typename> class Column;
template<typename> class Dump;
template<typename> class Cross;
template<typename> class Field3d;
template<typename> class Timedep;
template<typename> class Timeloop;

/**
 * Class for the dry thermodynamics.
 * This class is responsible for the computation of the right hand side term related to
 * the acceleration by buoyancy. In the dry thermodynamics temperature and buoyancy are
 * equivalent and no complex buoyancy function is required.
 */
template<typename TF>
class Thermo_dry : public Thermo<TF>
{
    public:
        Thermo_dry(Master&, Grid<TF>&, Fields<TF>&, Input&, const Sim_mode); // Constructor of the dry thermodynamics class.
        virtual ~Thermo_dry(); // Destructor of the dry thermodynamics class.

        void init();
        void create(Input&, Netcdf_handle&, Stats<TF>&, Column<TF>&, Cross<TF>&, Dump<TF>&, Timeloop<TF>&);
        void exec(const double, Stats<TF>&); // Add the tendencies belonging to the buoyancy.
        unsigned long get_time_limit(unsigned long, double); // Compute the time limit (n/a for thermo_dry).

        void exec_stats(Stats<TF>&);
        void exec_cross(Cross<TF>&, unsigned long);
        void exec_dump(Dump<TF>&, unsigned long);
        void exec_column(Column<TF>&);

        bool check_field_exists(std::string name);
        void get_thermo_field(
                Field3d<TF>&, const std::string&, const bool, const bool);
        void get_buoyancy_surf(std::vector<TF>&, std::vector<TF>&, bool);
        void get_buoyancy_surf(std::vector<TF>&, std::vector<TF>&, std::vector<TF>&)
            { throw std::runtime_error("Function get_buoyancy_surf not implemented"); }
        void get_buoyancy_fluxbot(std::vector<TF>&, bool);
        void get_temperature_bot(Field3d<TF>&, bool);
        const std::vector<TF>& get_p_vector() const;
        const std::vector<TF>& get_ph_vector() const;
        const std::vector<TF>& get_exner_vector() const;
        const std::vector<TF>& get_basestate_vector(std::string) const;
        TF get_db_ref() const;

        int get_bl_depth();
        TF get_buoyancy_diffusivity();

        void get_prog_vars(std::vector<std::string>&); // Retrieve a list of prognostic variables.

        #ifdef USECUDA
        // GPU functions and variables
        void prepare_device();
        void clear_device();
        void forward_device();
        void backward_device();
        void get_thermo_field_g(
                Field3d<TF>&, const std::string&, const bool);
        void get_buoyancy_surf_g(Field3d<TF>&);
        void get_buoyancy_fluxbot_g(Field3d<TF>&);
        TF* get_basestate_fld_g(std::string)
            { throw std::runtime_error("Function get_basestate_fld_g not implemented"); };
        void get_radiation_fields_g(Field3d<TF>&, Field3d<TF>&, Field3d<TF>&, Field3d<TF>&, Field3d<TF>&) const
            { throw std::runtime_error("Function get_radiation_fields_g not implemented"); }
        void get_radiation_columns_g(Field3d<TF>&, const int*, const int*, const int) const
            { throw std::runtime_error("Function get_radiation_columns_g not implemented"); }
        void get_land_surface_fields_g(TF*, TF*, TF*, TF*, TF*)
            { throw std::runtime_error("Function \"get_land_surface_fields_g\" not implemented in thermo_dry"); };
        void get_buoyancy_surf_g(TF*, TF*, TF*)
            { throw std::runtime_error("Function \"get_buoyancy_surf_g\" not implemented in thermo_disabled"); };
        #endif

        // Functions which are not implemented and throw error
        void get_radiation_fields(Field3d<TF>&, Field3d<TF>&, Field3d<TF>&, Field3d<TF>&, Field3d<TF>&) const
            { throw std::runtime_error("Function get_radiation_fields not implemented"); }
        void get_radiation_columns(Field3d<TF>&, std::vector<int>&, std::vector<int>&) const
            { throw std::runtime_error("Function get_radiation_columns not implemented"); }
        void get_land_surface_fields(
                std::vector<TF>&, std::vector<TF>&, std::vector<TF>&,
                std::vector<TF>&, std::vector<TF>&)
            { throw std::runtime_error("Function get_land_surface_fields not implemented"); }

        // Empty functions that are allowed to pass.
        void create_basestate(Input&, Netcdf_handle&) {};
        void load(const int) {};
        void save(const int) {};

        void get_mask(Stats<TF>&, std::string) {};
        bool has_mask(std::string) {return false;};

        void update_time_dependent(Timeloop<TF>&);

    private:
        using Thermo<TF>::swthermo;
        using Thermo<TF>::master;
        using Thermo<TF>::grid;
        using Thermo<TF>::fields;

        Boundary_cyclic<TF> boundary_cyclic;

        // cross sections
        std::vector<std::string> crosslist;        // List with all crosses from ini file
        std::vector<std::string> allowedcrossvars; // List with allowed cross variables
        bool swcross_b;
        std::vector<std::string> dumplist;         // List with all 3d dumps from the ini file.

        void create_stats(Stats<TF>&);   // Initialization of the statistics.
        void create_column(Column<TF>&); // Initialization of the single column output.
        void create_dump(Dump<TF>&);     // Initialization of the single column output.
        void create_cross(Cross<TF>&);   // Initialization of the single column output.

        enum class Basestate_type {anelastic, boussinesq};
        struct background_state
        {
            Basestate_type swbasestate;

            TF pbot;   // Surface pressure.
            TF thref0; // Reference potential temperature in case of Boussinesq

            std::vector<TF> thref;
            std::vector<TF> threfh;
            std::vector<TF> pref;
            std::vector<TF> prefh;
            std::vector<TF> exnref;
            std::vector<TF> exnrefh;

            // GPU functions and variables
            TF*  thref_g;
            TF*  threfh_g;
            TF*  pref_g;
            TF*  prefh_g;
            TF*  exnref_g;
            TF*  exnrefh_g;
        };
        background_state bs;
        background_state bs_stats;

        bool swbaroclinic;
        TF dthetady_ls;

        std::unique_ptr<Timedep<TF>> tdep_pbot;
        const std::string tend_name = "buoy";
        const std::string tend_longname = "Buoyancy";
};
#endif
