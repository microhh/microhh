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

#ifndef CANOPY_H
#define CANOPY_H

#include <vector>
#include <string>

class Master;
class Input;
class Netcdf_handle;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Field3d_operators;
template<typename> class Timedep;
template<typename> class Stats;
template<typename> class Thermo;
template<typename> class Field3d_io;

template<typename TF>
class Canopy
{
    public:
        Canopy(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Canopy();

        void init();
        void create(Input&, Netcdf_handle&, Stats<TF>&);
        void exec();

        // GPU functions and variables
        #ifdef USECUDA
        void prepare_device();
        void clear_device();
        #endif

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Field3d_operators<TF> field3d_operators;
        Field3d_io<TF> field3d_io;

        // Internal switches.
        bool sw_canopy;
        bool sw_3d_pad;

        // Canopy settings.
        int ktot_canopy;      // Vertical extent canopy (-).
        int kend_canopy;      // Vertical extent canopy (-).
        TF cd;                // Drag coeffient canopy (-).

        // Vertical profiles
        std::vector<TF> pad;   // Plant area density (m2 m-2)
        std::vector<TF> padh;  // Plant area density (m2 m-2)

        #ifdef USECUDA
        TF* pad_g;
        TF* padh_g;
        #endif
};
#endif
