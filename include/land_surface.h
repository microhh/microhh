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

#ifndef LAND_SURFACE
#define LAND_SURFACE

class Netcdf_file;

template<typename TF>
struct Surface_tile
{
    // Grid point fraction tile
    std::vector<TF> fraction;

    std::vector<TF> H;
    std::vector<TF> LE;
    std::vector<TF> G;

    std::vector<TF> T_bot;
    std::vector<TF> qt_bot;
    std::vector<TF> thl_bot;

    std::vector<TF> thl_fluxbot;
    std::vector<TF> qt_fluxbot;
};

template<typename TF>
using Tile_map = std::map<std::string, Surface_tile<TF>>;

template<typename TF>
class Land_surface
{
    public:
        Land_surface(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Land_surface();

        void init();
        void exec();
        void exec_stats(Stats<TF>&);
        void exec_cross(Cross<TF>&, unsigned long);

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;

        bool sw_land_surface;

        Tile_map<TF> tiles;
};
#endif
