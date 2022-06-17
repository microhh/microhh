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

#ifndef THERMO_H
#define THERMO_H

class Master;
class Input;
class Netcdf_handle;
template<typename> class Grid;
template<typename> class Stats;
template<typename> class Advec;
template<typename> class Diff;
template<typename> class Column;
template<typename> class Dump;
template<typename> class Cross;
template<typename> class Field3d;
template<typename> class Timeloop;

enum class Sim_mode;

/**
 * Base class for the thermo scheme. This class is abstract and only
 * derived classes can be instantiated. Derived classes are
 * implemented that handle different thermodynamics.
 */
template<typename TF>
class Thermo
{
    public:
        Thermo(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Thermo();
        static std::shared_ptr<Thermo> factory(Master&, Grid<TF>&, Fields<TF>&, Input&, const Sim_mode);
        const std::string& get_switch() const;

        // Below are the functions that the derived class has to implement.
        virtual void init() = 0;
        virtual void create(
                Input&, Netcdf_handle&, Stats<TF>&, Column<TF>&, Cross<TF>&, Dump<TF>&) = 0;
        virtual void create_basestate(Input&, Netcdf_handle&) = 0;
        virtual unsigned long get_time_limit(unsigned long, double) = 0;
        virtual void load(const int) = 0;
        virtual void save(const int) = 0;

        virtual void exec(const double, Stats<TF>&) = 0;
        virtual void exec_stats(Stats<TF>&) = 0; ///< Calculate the statistics
        virtual void exec_column(Column<TF>&) = 0; ///< Output the column
        virtual void exec_dump(Dump<TF>&, unsigned long) = 0;
        virtual void exec_cross(Cross<TF>&, unsigned long) = 0;

        virtual void get_mask(Stats<TF>&, std::string) = 0;
        virtual bool has_mask(std::string) = 0;

        // Interfacing functions to get buoyancy properties from other classes.
        virtual bool check_field_exists(std::string name) = 0;
        virtual void get_thermo_field(
                Field3d<TF>&, const std::string&, const bool, const bool) = 0;
        virtual void get_buoyancy_surf(std::vector<TF>&, std::vector<TF>&, bool) = 0;
        virtual void get_buoyancy_surf(std::vector<TF>&, std::vector<TF>&, std::vector<TF>&) = 0;
        virtual void get_buoyancy_fluxbot(std::vector<TF>&, bool) = 0;
        virtual void get_temperature_bot(Field3d<TF>&, bool) = 0;
        virtual void get_prog_vars(std::vector<std::string>&) = 0;

        virtual void get_radiation_fields(
                Field3d<TF>&, Field3d<TF>&, Field3d<TF>&, Field3d<TF>&, Field3d<TF>&) const = 0;
        virtual void get_radiation_columns(Field3d<TF>&, std::vector<int>&, std::vector<int>&) const = 0;
        virtual void get_land_surface_fields(
                std::vector<TF>&, std::vector<TF>&, std::vector<TF>&, std::vector<TF>&, std::vector<TF>&) = 0;

        virtual const std::vector<TF>& get_basestate_vector(std::string) const = 0;
        virtual TF get_db_ref() const = 0;

        virtual int get_bl_depth() = 0;
        virtual TF get_buoyancy_diffusivity() = 0;

        virtual void update_time_dependent(Timeloop<TF>&) = 0;

        #ifdef USECUDA
        // GPU functions and variables.
        virtual void prepare_device() = 0;
        virtual void clear_device() = 0;
        virtual void forward_device() = 0;
        virtual void backward_device() = 0;
        virtual void get_thermo_field_g(Field3d<TF>&, const std::string&, const bool) = 0;
        virtual void get_buoyancy_surf_g(Field3d<TF>&)  = 0;
        virtual void get_buoyancy_surf_g(TF*, TF*, TF*)  = 0;
        virtual void get_buoyancy_fluxbot_g(Field3d<TF>&) = 0;
        virtual void get_land_surface_fields_g(TF*, TF*, TF*, TF*, TF*) = 0;
        virtual TF* get_basestate_fld_g(std::string) = 0;

        virtual void get_radiation_fields_g(
                Field3d<TF>&, Field3d<TF>&, Field3d<TF>&, Field3d<TF>&, Field3d<TF>&) const = 0;
        virtual void get_radiation_columns_g(Field3d<TF>&, const int*, const int*, const int) const = 0;
        #endif

    protected:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        std::string swthermo;
};
#endif
