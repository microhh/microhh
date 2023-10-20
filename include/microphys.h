/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

#ifndef MICROPHYS_H
#define MICROPHYS_H

class Master;
class Input;
class Netcdf_handle;

template<typename> class Grid;
template<typename> class Stats;
template<typename> class Diff;
template<typename> class Dump;
template<typename> class Cross;
template<typename> class Thermo;
template<typename> class Field3d;

enum class Microphys_type {Disabled, Warm_2mom, Nsw6};

/**
 * Base class for the microphysics scheme. This class is abstract and only
 * derived classes can be instantiated. Derived classes are
 * implemented that handle different microphysiscs schemes.
 */
template<typename TF>
class Microphys
{
    public:
        Microphys(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Microphys();

        static std::shared_ptr<Microphys> factory(Master&, Grid<TF>&, Fields<TF>&, Input&);
        Microphys_type get_switch();

        // Below are the functions that the derived class has to implement.
        virtual void init() = 0;
        virtual void create(Input&, Netcdf_handle&, Stats<TF>&, Cross<TF>&, Dump<TF>&, Column<TF>&) = 0;
        virtual unsigned long get_time_limit(unsigned long, double) = 0;

        virtual void exec(Thermo<TF>&, const double, Stats<TF>&) = 0;
        virtual void exec_stats(Stats<TF>&, Thermo<TF>&, const double) = 0; ///< Calculate the statistics
        virtual void exec_column(Column<TF>&) = 0;

        virtual void exec_dump(Dump<TF>&, unsigned long) = 0;
        virtual void exec_cross(Cross<TF>&, unsigned long) = 0;

        virtual void get_mask(Stats<TF>&, std::string) = 0;
        virtual bool has_mask(std::string) = 0;

        virtual void get_surface_rain_rate(std::vector<TF>&) = 0;

        virtual TF get_Nc0() = 0;
        virtual TF get_Ni0() = 0;

        // GPU functions and variables.
        #ifdef USECUDA
        virtual void get_surface_rain_rate_g(TF*) = 0;
        virtual void prepare_device() = 0;
        virtual void clear_device() = 0;
        virtual void backward_device() = 0;
        #endif

    protected:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Field3d_operators<TF> field3d_operators;

        Microphys_type swmicrophys;
};
#endif
