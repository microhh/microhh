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
#include <algorithm>
#include <memory>
#include "Gas_concs.h"

using Aerosol_concs = Gas_concs;

class Master;
class Input;
class Netcdf_handle;
class Netcdf_file;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Timedep;
template<typename> class Timeloop;
template<typename> class Stats;
template<typename> class Thermo;
template<typename> class Field3d;

template<typename TF>
class Aerosol
{
    public:
        Aerosol(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Aerosol();

        void init();
        void create(Input&, Netcdf_handle&, Stats<TF>&);
        void exec_stats(Stats<TF>&);
        void update_time_dependent(Timeloop<TF>&);

        #ifdef USECUDA
        // GPU functions and variables
        void prepare_device();
        void clear_device();
        #endif

        void get_radiation_fields(Aerosol_concs&);

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        // Case switches
        bool sw_aerosol;
        bool sw_timedep;

        // Arrays
        // to fill with input values
        std::vector<TF> aermr01;
        std::vector<TF> aermr02;
        std::vector<TF> aermr03;
        std::vector<TF> aermr04;
        std::vector<TF> aermr05;
        std::vector<TF> aermr06;
        std::vector<TF> aermr07;
        std::vector<TF> aermr08;
        std::vector<TF> aermr09;
        std::vector<TF> aermr10;
        std::vector<TF> aermr11;

        std::unique_ptr<Timedep<TF>> tdep_aermr01;
        std::unique_ptr<Timedep<TF>> tdep_aermr02;
        std::unique_ptr<Timedep<TF>> tdep_aermr03;
        std::unique_ptr<Timedep<TF>> tdep_aermr04;
        std::unique_ptr<Timedep<TF>> tdep_aermr05;
        std::unique_ptr<Timedep<TF>> tdep_aermr06;
        std::unique_ptr<Timedep<TF>> tdep_aermr07;
        std::unique_ptr<Timedep<TF>> tdep_aermr08;
        std::unique_ptr<Timedep<TF>> tdep_aermr09;
        std::unique_ptr<Timedep<TF>> tdep_aermr10;
        std::unique_ptr<Timedep<TF>> tdep_aermr11;

        #ifdef USECUDA
        // GPU functions and variables
        TF* aermr01_g;  ///< Pointer to GPU array.
        TF* aermr02_g;
        TF* aermr03_g;
        TF* aermr04_g;
        TF* aermr05_g;
        TF* aermr06_g;
        TF* aermr07_g;
        TF* aermr08_g;
        TF* aermr09_g;
        TF* aermr10_g;
        TF* aermr11_g;
        #endif

};
#endif
