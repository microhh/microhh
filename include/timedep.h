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

#ifndef TIMEDEP_H
#define TIMEDEP_H

#include "master.h"
#include "grid.h"
#include "timeloop.h"

class Master;
template<typename> class Grid;

enum class Timedep_switch {Disabled, Enabled};

template<typename TF>
class Timedep
{
    public:
        Timedep(Master&, Grid<TF>&, const std::string, const bool);
        ~Timedep();

        void create_timedep(Netcdf_handle&, const std::string);
        void update_time_dependent(TF&, Timeloop<TF>&);

        void create_timedep_prof(Netcdf_handle&, const TF, const std::string);
        void create_timedep_background_prof(Netcdf_handle&, const TF, const std::string, const TF);
        void update_time_dependent_prof(std::vector<TF>&, Timeloop<TF>&);
        void update_time_dependent_background_prof(std::vector<TF>&, Timeloop<TF>&, const TF);

        #ifdef USECUDA
        TF* data_g;
        void update_time_dependent_prof_g(TF*, Timeloop<TF>&);
        void prepare_device();
        void clear_device();
        #endif

    private:
        Master& master;
        Grid<TF>& grid;
        const std::string varname;

        Timedep_switch sw;

        std::vector<double> time;
        std::vector<TF> data;
};
#endif
