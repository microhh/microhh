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

#include "soil_grid.h"
#include "grid.h"
#include "master.h"
#include "input.h"
#include "defines.h"
#include "netcdf_interface.h"

/**
 * This function constructs the grid class.
 * @param modelin Pointer to the model class.
 * @param inputin Pointer to the input class.
 */
template<typename TF>
Soil_grid<TF>::Soil_grid(Master& masterin, Grid<TF>& gridin, Input& input) :
    master(masterin), grid(gridin)
{
    sw_soil_enabled = input.get_item<bool>("soil", "sw_soil", "", false);

    if (sw_soil_enabled)
    {
        gd.zsize = input.get_item<TF> ("soil", "zsize", "");
        gd.ktot  = input.get_item<int>("soil", "ktot",  "");
    }
}

template<typename TF>
Soil_grid<TF>::~Soil_grid()
{
}

/**
 * This function allocates the dynamic arrays in the field class
 * variables and calculates the derived grid indices and dimensions.
 */
template<typename TF>
void Soil_grid<TF>::init()
{
    if (!sw_soil_enabled)
        return;

    auto& agd = grid.get_grid_data();   // Atmospheric grid data

    gd.kmax    = gd.ktot;
    gd.kmaxh   = gd.ktot+1;
    gd.kgc     = 0;

    gd.kcells  = gd.kmax +2*gd.kgc; 
    gd.kcellsh = gd.kmaxh+2*gd.kgc; 

    gd.ncells  = gd.kcells  * agd.ijcells;
    gd.ncellsh = gd.kcellsh * agd.ijcells;

    gd.kstart  = gd.kgc;
    gd.kend    = gd.kmax+gd.kgc;

    // Allocate all arrays
    gd.z   .resize(gd.kcells);
    gd.dz  .resize(gd.kcells);
    gd.dzi .resize(gd.kcells);

    gd.zh  .resize(gd.kcellsh);
    gd.dzh .resize(gd.kcellsh);
    gd.dzhi.resize(gd.kcellsh);
}

/**
 * This function initializes the fields containing the grid dimensions based
 * on the profiles in the input file.
 * @param inputin Pointer to the input class.
 */
template<typename TF>
void Soil_grid<TF>::create(Netcdf_handle& input_nc)
{
    if (!sw_soil_enabled)
        return;

    // Get the grid coordinates from the input. This one is read from the global scope.
    input_nc.get_variable(gd.z, "z_soil", {0}, {gd.ktot});
    std::rotate(gd.z.rbegin(), gd.z.rbegin() + gd.kstart, gd.z.rend());

    // Calculate the grid
    for (int k=gd.kstart+1; k<gd.kend; ++k)
        gd.zh[k] = 0.5*(gd.z[k-1] + gd.z[k]);

    gd.zh[gd.kend  ] = 0.;
    gd.zh[gd.kstart] = gd.zsize;

    // Calculate grid spacing
    for (int k=gd.kstart; k<gd.kend; ++k)
        gd.dz[k] = gd.zh[k+1] - gd.zh[k];

    for (int k=gd.kstart+1; k<gd.kend; ++k)
        gd.dzh[k] = gd.z[k] - gd.z[k-1];

    gd.dzh[gd.kend  ] = 2*-gd.z[gd.kend-1];
    gd.dzh[gd.kstart] = 2*(gd.z[gd.kstart] - gd.zh[gd.kstart]);

    // Inverse grid spacings
    for (int k=gd.kstart; k<gd.kend; ++k)
        gd.dzi[k] = 1./gd.dz[k];

    for (int k=gd.kstart; k<gd.kend+1; ++k)
        gd.dzhi[k] = 1./gd.dzh[k];
}

template<typename TF>
const Soil_grid_data<TF>& Soil_grid<TF>::get_grid_data()
{
    return gd;
}

template class Soil_grid<double>;
template class Soil_grid<float>;
