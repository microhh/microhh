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

#ifndef THERMO_DISABLED
#define THERMO_DISABLED

#include "thermo.h"

class Master;
class Input;
class Data_block;
template<typename> class Grid;
template<typename> class Stats;
template<typename> class Diff;
template<typename> class Column;
template<typename> class Dump;
template<typename> class Cross;
template<typename> class Field3d;
template<typename> class Thermo;


template<typename TF>
class Thermo_disabled : public Thermo<TF>
{
    public:
        Thermo_disabled(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Thermo_disabled();

        // Interfacing functions to get buoyancy properties from other classes.
        bool check_field_exists(std::string name);

        // Empty functions that are allowed to pass.

        // Interfacing functions to get buoyancy properties from other classes.
        void init() {};
        void create(Input&, Data_block&, Stats<TF>&, Column<TF>&, Cross<TF>&, Dump<TF>&) {};
        void exec() {};
        void exec_stats(Stats<TF>&, std::string, Field3d<TF>&, Field3d<TF>&, const Diff<TF>&) {};
        void exec_column(Column<TF>&) {};
        virtual void exec_dump(Dump<TF>&, unsigned long) {};
        virtual void exec_cross(Cross<TF>&, unsigned long) {};
        void get_mask(Field3d<TF>&, Field3d<TF>&, Stats<TF>&, std::string) {}
        void get_prog_vars(std::vector<std::string>) {};
        void update_time_dependent() {};
        TF get_buoyancy_diffusivity();

        unsigned long get_time_limit(unsigned long, double);

#ifdef USECUDA
        void prepare_device() {};
        void clear_device() {};
        void forward_device() {};
        void backward_device() {};
#endif

        // Empty functions that shall throw.
        void get_thermo_field(Field3d<TF>&, std::string, bool) { throw 1; }
        void get_buoyancy_surf(Field3d<TF>&) { throw 1; }
        void get_buoyancy_fluxbot(Field3d<TF>&) { throw 1; }
    private:
        using Thermo<TF>::swthermo;
};
#endif
