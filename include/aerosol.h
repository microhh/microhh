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

#ifndef AEROSOL_H
#define AEROSOL_H

#include <vector>
#include <string>
#include <map>

class Master;
class Input;
class Netcdf_handle;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Timedep;
template<typename> class Timeloop;
template<typename> class Stats;
template<typename> class Thermo;

template<typename TF>
class Aerosol
{
    public:
        Aerosol(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Aerosol();

        void init();
        void create(Input&, Netcdf_handle&, Stats<TF>&);
        void exec(Thermo<TF>&);
        void exec_stats(Stats<TF>&);
        void update_time_dependent(Timeloop<TF>&);

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        // Case switches
        bool sw_aerosol;
        bool sw_timedep;

        // Arrays
        std::vector<TF> bla;
};
#endif
