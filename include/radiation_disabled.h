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

#ifndef RADIATION_DISABLED_H
#define RADIATION_DISABLED_H

#include <cstdio>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "radiation.h"

class Master;
class Input;
template<typename> class Grid;
template<typename> class Stats;
template<typename> class Diff;
template<typename> class Column;
template<typename> class Dump;
template<typename> class Cross;
template<typename> class Field3d;
template<typename> class Thermo;
template<typename> class Timeloop;

template<typename TF>
class Radiation_disabled : public Radiation<TF>
{
    public:
        Radiation_disabled(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Radiation_disabled();

        bool check_field_exists(std::string name);
        void init(Timeloop<TF>&) {};
        void create(
                Input&, Netcdf_handle&, Thermo<TF>&,
                Stats<TF>&, Column<TF>&, Cross<TF>&, Dump<TF>&) {};
        void exec(Thermo<TF>&, double, Timeloop<TF>&, Stats<TF>&) {};

        unsigned long get_time_limit(unsigned long);

        // Empty functions that should throw
        void get_radiation_field(Field3d<TF>&, std::string, Thermo<TF>&, Timeloop<TF>&)
            { throw std::runtime_error("\"get_radiation_field()\" is not implemented in radiation_disabled"); }

        std::vector<TF>& get_surface_radiation(std::string)
            { throw std::runtime_error("\"get_surface_radiation()\" is not implemented in radiation_disabled"); }

        // void exec_stats(Stats<TF>&, Thermo<TF>&, Timeloop<TF>&) {};
        // void exec_cross(Cross<TF>&, unsigned long, Thermo<TF>&, Timeloop<TF>&) {};
        // void exec_dump(Dump<TF>&, unsigned long, Thermo<TF>&, Timeloop<TF>&) {};

        void exec_all_stats(
                Stats<TF>&, Cross<TF>&, Dump<TF>&, Column<TF>&,
                Thermo<TF>&, Timeloop<TF>&,
                const unsigned long, const int) {};
        void exec_individual_column_stats(
                Column<TF>&, Thermo<TF>&, Timeloop<TF>&, Stats<TF>&) {};
        void exec_column(Column<TF>&, Thermo<TF>&, Timeloop<TF>&) {};

    private:
        using Radiation<TF>::swradiation;
};
#endif
