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
#include <algorithm>
#include <cmath>
#include <math.h>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "cross.h"
#include "netcdf_interface.h"
#include "constants.h"

#include "land_surface.h"

using namespace Constants;

namespace
{
    template<typename TF>
    void init_tile(Surface_tile<TF>& tile, const int ijcells)
    {
        tile.fraction.resize(ijcells);

        tile.H.resize(ijcells);
        tile.LE.resize(ijcells);
        tile.G.resize(ijcells);

        tile.T_bot.resize(ijcells);
        tile.thl_bot.resize(ijcells);
        tile.qt_bot.resize(ijcells);

        tile.thl_fluxbot.resize(ijcells);
        tile.qt_fluxbot.resize(ijcells);
    }
}

template<typename TF>
Land_surface<TF>::Land_surface(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin) 
{
    sw_land_surface = inputin.get_item<bool>("land_surface", "sw_land_surface", "", false);

    if (sw_land_surface)
    {
        // Create the surface tiles
        tiles.emplace("low_veg",   Surface_tile<TF>{});
        tiles.emplace("bare_soil", Surface_tile<TF>{});
        tiles.emplace("wet_skin",  Surface_tile<TF>{});
    }
}

template<typename TF>
Land_surface<TF>::~Land_surface()
{
}

template<typename TF>
void Land_surface<TF>::init()
{
    if (!sw_land_surface)
        return;

    auto& gd = grid.get_grid_data();

    // Allocate the surface tiles
    for (auto& tile : tiles)
        init_tile(tile.second, gd.ijcells);
}

template<typename TF>
void Land_surface<TF>::exec()
{
    if (!sw_land_surface)
        return;
}

template<typename TF>
void Land_surface<TF>::exec_stats(Stats<TF>& stats)
{
    if (!sw_land_surface)
        return;
}

template<typename TF>
void Land_surface<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    if (!sw_land_surface)
        return;
}

template class Land_surface<double>;
template class Land_surface<float>;
