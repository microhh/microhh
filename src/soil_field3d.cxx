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

#include <iostream>

#include "master.h"
#include "grid.h"
#include "soil_grid.h"
#include "soil_field3d.h"

template<typename TF>
Soil_field3d<TF>::Soil_field3d(
        Master& masterin, Grid<TF>& gridin, Soil_grid<TF>& soilgridin,
        std::string namein, std::string longnamein, std::string unitin) :
    master(masterin),
    grid(gridin),
    soil_grid(soilgridin)
{
    name = namein;
    longname = longnamein;
    unit = unitin;
}

template<typename TF>
Soil_field3d<TF>::~Soil_field3d()
{
}

template<typename TF>
int Soil_field3d<TF>::init()
{
    // Allocate all fields belonging to the soil field
    auto& agd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    fld     .resize(sgd.ncells);
    fld_bot .resize(agd.ijcells);
    fld_top .resize(agd.ijcells);
    flux_bot.resize(agd.ijcells);
    flux_top.resize(agd.ijcells);

    // Set all values to zero
    for (int n=0; n<sgd.ncells; ++n)
        fld [n] = 0.;

    for (int n=0; n<agd.ijcells; ++n)
    {
        fld_bot[n]  = 0.;
        fld_top[n]  = 0.;
        flux_bot[n] = 0.;
        flux_top[n] = 0.;
    }

    return 0;
}


#ifdef FLOAT_SINGLE
template class Soil_field3d<float>;
#else
template class Soil_field3d<double>;
#endif
