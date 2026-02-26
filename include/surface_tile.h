/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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

#ifndef SURFACE_TILE_H
#define SURFACE_TILE_H

#include <map>
#include <string>
#include <vector>

template<typename TF>
struct Surface_tile
{
    std::string long_name;    // Descriptive name of tile

    // Shared
    std::vector<TF> fraction; // Grid point fraction tile (-)
    std::vector<TF> thl_bot;  // Skin (liquid water) potential temperature (K)
    std::vector<TF> qt_bot;   // Skin specific humidity (kg kg-1)

    // Surface layer
    std::vector<TF> obuk;     // Obukhov length (m)
    std::vector<TF> ustar;    // Friction velocity (m s-1)
    std::vector<TF> bfluxbot; // Buoyancy flux at bottom (m s-1)
    std::vector<int> nobuk;   // Index in LUT
    std::vector<TF> ra;       // Aerodynamic resistance (s m-1)

    // Land surface
    std::vector<TF> rs;       // Surface resistance (canopy or soil, s m-1)
    std::vector<TF> H;        // Sensible heat flux (W m-2)
    std::vector<TF> LE;       // Latent heat flux (W m-2)
    std::vector<TF> G;        // Soil heat flux (W m-2)
    std::vector<TF> S;        // Storage flux (W m-2)

    #ifdef USECUDA
    // Shared
    TF* fraction_g; // Grid point fraction tile (-)
    TF* thl_bot_g;  // Skin (liquid water) potential temperature (K)
    TF* qt_bot_g;   // Skin specific humidity (kg kg-1)

    // Surface layer
    TF* obuk_g;     // Obukhov length (m)
    TF* ustar_g;    // Friction velocity (m s-1)
    TF* bfluxbot_g; // Buoyancy flux at bottom (m s-1)
    int* nobuk_g;   // Index in LUT
    TF* ra_g;       // Aerodynamic resistance (s m-1)

    // Land surface
    TF* rs_g;       // Surface resistance (canopy or soil, s m-1)
    TF* H_g;        // Sensible heat flux (W m-2)
    TF* LE_g;       // Latent heat flux (W m-2)
    TF* G_g;        // Soil heat flux (W m-2)
    TF* S_g;        // Storage flux (W m-2)
    #endif
};

template<typename TF>
using Tile_map = std::map<std::string, Surface_tile<TF>>;

#endif
