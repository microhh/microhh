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

#ifndef SOIL_DISABLED
#define SOIL_DISABLED

#include "soil.h"

template<typename TF>
class Soil_disabled : public Soil<TF>
{
    public:
        Soil_disabled(Master&, Grid<TF>&, Fields<TF>&, Input&);
        virtual ~Soil_disabled();

        void init() {};
        void create_cold_start(Input&, Netcdf_handle&) {};
        void create_fields_grid_stats(Input&, Netcdf_handle&, Stats<TF>&, Cross<TF>&) {};
        void save_prognostic_fields(int) {};
        void load_prognostic_fields(int) {};
        void calc_tendencies() {};

    private:
        using Soil<TF>::sw_soil;
};
#endif
