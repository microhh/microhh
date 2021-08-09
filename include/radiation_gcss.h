/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
 * Copyright (c) 2018-2019 Elynn Wu
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

#ifndef RADIATION_GCSS_H
#define RADIATION_GCSS_H

#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>

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
template<typename> class Timedep;
template<typename> class Timeloop;

template<typename TF>
class Radiation_gcss : public Radiation<TF>
{
    public:
        Radiation_gcss(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Radiation_gcss();
        void init(Timeloop<TF>&);
        void create(
                Input&, Netcdf_handle&, Thermo<TF>&,
                Stats<TF>&, Column<TF>&, Cross<TF>&, Dump<TF>&);
        void exec(Thermo<TF>&, double, Timeloop<TF>&, Stats<TF>&);

        unsigned long get_time_limit(unsigned long);

        bool check_field_exists(std::string name);
        void get_radiation_field(Field3d<TF>&, std::string, Thermo<TF>&, Timeloop<TF>&);

        std::vector<TF>& get_surface_radiation(std::string)
            { throw std::runtime_error("\"get_surface_radiation()\" is not implemented in radiation_disabled"); }

        void exec_all_stats(
                Stats<TF>&, Cross<TF>&, Dump<TF>&, Column<TF>&,
                Thermo<TF>&, Timeloop<TF>&,
                const unsigned long, const int);

        void exec_column(Column<TF>&, Thermo<TF>&, Timeloop<TF>&);

    private:
        void create_stats(Stats<TF>&);   ///< Initialization of the statistics.
        void create_column(Column<TF>&); ///< Initialization of the single column output.
        void create_dump(Dump<TF>&);     ///< Initialization of the single column output.
        void create_cross(Cross<TF>&);   ///< Initialization of the single column output.

        using Radiation<TF>::swradiation;
        using Radiation<TF>::master;
        using Radiation<TF>::grid;
        using Radiation<TF>::fields;
        using Radiation<TF>::field3d_operators;

        std::vector<std::string> available_masks;  // Vector with the masks

        std::vector<std::string> crosslist;        ///< List with all crosses from ini file.
        bool swcross_rflx;
        std::vector<std::string> dumplist;         ///< List with all 3d dumps from the ini file.

        TF lat;
        TF lon;
        TF xka;
        TF fr0;
        TF fr1;
        TF div;

        const TF mu_min = 0.035;

        #ifdef USECUDA
        // GPU functions and variables
        void get_radiation_field_g(Field3d<TF>&, std::string, Thermo<TF>&, Timeloop<TF>&);
        void prepare_device() {};
        void clear_device() {};
        #endif

        const std::string tend_name = "rad";
        const std::string tend_longname = "Radiation";
};
#endif
