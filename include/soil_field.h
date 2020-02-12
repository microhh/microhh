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

#ifndef SOIL_FIELD
#define SOIL_FIELD

#include <string>
#include <vector>
#include <array>

class Master;
template<typename> class Grid;

template<typename TF>
class Soil_field
{
    public:
        // Functions
        Soil_field(Master&, Grid<TF>&);
        ~Soil_field();

        void init(int);

        // Variables at CPU.
        std::vector<TF> fld;
        std::vector<TF> tend;

    private:
        Master& master;
        Grid<TF>& grid;
};
#endif
