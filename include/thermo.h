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

#ifndef THERMO
#define THERMO

class Master;
class Input;
template<typename> class Grid;
template<typename> class Stats;
template<typename> class Diff;
template<typename> class Column;
template<typename> class Dump;
template<typename> class Cross;
template<typename> class Field3d;
class Data_block;
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
        static std::shared_ptr<Thermo> factory(Master&, Grid<TF>&, Fields<TF>&, Input&); ///< Factory function for thermo class generation.
        std::string get_switch();

        // Below are the functions that the derived class has to implement.

        virtual void init() = 0;
        virtual void create(Input&, Data_block&, Stats<TF>&, Column<TF>&, Cross<TF>&, Dump<TF>&) = 0;
        virtual unsigned long get_time_limit(unsigned long, double) = 0;

        virtual void exec() = 0;
        virtual void exec_stats(Stats<TF>&, std::string, Field3d<TF>&, Field3d<TF>&, const Diff<TF>&) = 0;   ///< Calculate the statistics
        virtual void exec_column(Column<TF>&) = 0;   ///< Output the column
        virtual void exec_dump(Dump<TF>&, unsigned long) = 0;
        virtual void exec_cross(Cross<TF>&, unsigned long) = 0;

        virtual void get_mask(Field3d<TF>&, Field3d<TF>&, Stats<TF>&, std::string) = 0;

        // Interfacing functions to get buoyancy properties from other classes.
        virtual bool check_field_exists(std::string name) = 0;
        virtual void get_thermo_field(Field3d<TF>&, std::string, bool) = 0;
        virtual void get_buoyancy_surf(Field3d<TF>&) = 0;
        virtual void get_buoyancy_fluxbot(Field3d<TF>&) = 0;
        virtual void get_prog_vars(std::vector<std::string>) = 0;

        virtual TF get_buoyancy_diffusivity() = 0;

        virtual void update_time_dependent() = 0;

        #ifdef USECUDA
        // GPU functions and variables.
        virtual void prepare_device() = 0;
        virtual void clear_device() = 0;
        virtual void forward_device() = 0;
        virtual void backward_device() = 0;

        #endif

    protected:
        Master& master; ///< Pointer to master class.
        Grid<TF>& grid; ///< Pointer to grid class.
        Fields<TF>& fields; ///< Pointer to fields class.

        std::string swthermo;
};
#endif
