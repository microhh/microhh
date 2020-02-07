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

#include <cstdio>
#include <iostream>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "constants.h"
#include "netcdf_interface.h"

#include "soil.h"
#include "soil_enabled.h"


template<typename TF>
Soil_enabled<TF>::Soil_enabled(Master& masterin, Grid<TF>& gridin, Input& inputin) : 
    Soil<TF>(masterin, gridin, inputin)
{
    sw_soil = Soil_type::Enabled;

    // Read the soil settings from the ini file
    ktot = inputin.get_item<int>("soil", "ktot", "");
    sw_interactive = inputin.get_item<bool>("soil", "sw_interactive", "", false);
    sw_homogeneous = inputin.get_item<bool>("soil", "sw_homogeneous", "", true);

    // Checks on input & limitations
    if (sw_interactive)
        throw std::runtime_error("Interactive soil not (yet) implemented");
    if (!sw_homogeneous)
        throw std::runtime_error("Heterogeneous soil input not (yet) implemented");
}

template<typename TF>
Soil_enabled<TF>::~Soil_enabled()
{
}

template<typename TF>
void Soil_enabled<TF>::init()
{
    auto& gd = grid.get_grid_data();

    // Resize the vectors which hold the vertical grid info
    // Full level
    z.resize(ktot);
    dz.resize(ktot);
    dzi.resize(ktot);

    // Half level
    zh.resize(ktot+1);
    dzh.resize(ktot+1);
    dzhi.resize(ktot+1);

    // Resize the data vectors
    const int ncells_soil = gd.ijcells*ktot;
    theta.resize(ncells_soil);
    temperature.resize(ncells_soil);
}

template<typename TF>
void Soil_enabled<TF>::create(Input& input, Netcdf_handle& input_nc)
{
    // Calculate soil grid
    // 1. Read full level grid height from input NetCDF file
    Netcdf_group& soil_group = input_nc.get_group("soil");
    input_nc.get_variable(z, "z_soil", {0}, {ktot});

    // Etc. -> need to discuss with Chiel..
}

template class Soil_enabled<double>;
template class Soil_enabled<float>;
