/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
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

#ifndef SOIL_FIELD3D
#define SOIL_FIELD3D

#include <string>
#include <vector>
#include <array>

class Master;
template<typename> class Grid;
template<typename> class Soil_grid;

template<typename TF>
class Soil_field3d
{
    public:
        // Functions
        Soil_field3d(Master&, Grid<TF>&, Soil_grid<TF>&, std::string, std::string, std::string);
        ~Soil_field3d();

        int init();

        // Variables at CPU.
        std::vector<TF> fld;
        std::vector<TF> fld_bot;
        std::vector<TF> fld_top;
        std::vector<TF> flux_bot;
        std::vector<TF> flux_top;

        #ifdef USECUDA
        void init_device();
        void clear_device();

        TF* fld_g;
        TF* fld_bot_g;
        TF* fld_top_g;
        TF* flux_bot_g;
        TF* flux_top_g;
        #endif

        std::string name;
        std::string longname;
        std::string unit;

    private:
        Master& master;
        Grid<TF>& grid;
        Soil_grid<TF>& soil_grid;
};
#endif
