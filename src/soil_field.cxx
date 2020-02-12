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

#include "master.h"
#include "grid.h"
#include "soil_field.h"

template<typename TF>
Soil_field<TF>::Soil_field(Master& masterin, Grid<TF>& gridin) :
    master(masterin),
    grid(gridin)
{
}

template<typename TF>
Soil_field<TF>::~Soil_field()
{
}

template<typename TF>
void Soil_field<TF>::init(const int ktot)
{
    // Allocate all fields belonging to the soil field
    const Grid_data<TF>& gd = grid.get_grid_data();
    const int ncells_soil   = gd.ijcells *  ktot;

    fld .resize(ncells_soil);
    tend.resize(ncells_soil);

    // Set all values to zero
    for (int n=0; n<ncells_soil; ++n)
    {
        fld [n] = 0.;
        tend[n] = 0.;
    }
}

template class Soil_field<double>;
template class Soil_field<float>;
