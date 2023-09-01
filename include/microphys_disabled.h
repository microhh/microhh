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

#ifndef MICROPHYS_DISABLED_H
#define MICROPHYS_DISABLED_H

#include <stdexcept>
#include "microphys.h"

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
template<typename> class Microphys;

template<typename TF>
class Microphys_disabled : public Microphys<TF>
{
    public:
        Microphys_disabled(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Microphys_disabled();

        void init() {};
        void create(Input&, Netcdf_handle&, Stats<TF>&, Cross<TF>&, Dump<TF>&, Column<TF>&) {};
        void exec(Thermo<TF>&, const double, Stats<TF>&) {};
        void exec_stats(Stats<TF>&, Thermo<TF>&, const double) {};
        void exec_column(Column<TF>&) {};
        void exec_dump(Dump<TF>&, unsigned long) {};
        void exec_cross(Cross<TF>&, unsigned long) {};
        void get_mask(Stats<TF>&, std::string) {};
        bool has_mask(std::string) {return false;};

        void get_surface_rain_rate(std::vector<TF>&);

        unsigned long get_time_limit(unsigned long, double);

        TF get_Nc0() { throw std::runtime_error("Microphys_disabled cannot provide Nc0"); };
        TF get_Ni0() { throw std::runtime_error("Microphys_disabled cannot provide Ni0"); };

        #ifdef USECUDA
        void get_surface_rain_rate_g(TF*);
        void prepare_device() {};
        void clear_device() {};
        void backward_device() {};
        #endif

    private:
        using Microphys<TF>::swmicrophys;
        using Microphys<TF>::grid;
};
#endif
