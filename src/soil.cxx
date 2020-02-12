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
#include "fields.h"
#include "input.h"

// Soil schemes
#include "soil.h"
#include "soil_disabled.h"
#include "soil_enabled.h"

template<typename TF>
Soil<TF>::Soil(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& input) : 
    master(masterin), grid(gridin), fields(fieldsin)
{
}

template<typename TF>
Soil<TF>::~Soil()
{
}

template<typename TF>
Soil_type Soil<TF>::get_switch()
{
    return sw_soil;
}

template<typename TF>
std::shared_ptr<Soil<TF>> Soil<TF>::factory(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin)
{
    std::string sw_soil_str = inputin.get_item<std::string>("soil", "swsoil", "", "0");

    if (sw_soil_str == "0")
        return std::make_shared<Soil_disabled<TF>>(masterin, gridin, fieldsin, inputin);
    else if (sw_soil_str == "1")
        return std::make_shared<Soil_enabled<TF>>(masterin, gridin, fieldsin, inputin);
    else
    {
        std::string msg = sw_soil_str + " is an illegal value for swsoil";
        throw std::runtime_error(msg);
    }
}

template class Soil<double>;
template class Soil<float>;
