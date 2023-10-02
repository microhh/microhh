/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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
#include <cmath>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "input.h"
// thermo schemes
#include "thermo.h"
#include "thermo_buoy.h"
#include "thermo_dry.h"
#include "thermo_moist.h"
#include "thermo_disabled.h"

template<typename TF>
Thermo<TF>::Thermo(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& input) :
    master(masterin), grid(gridin), fields(fieldsin)
{
}

template<typename TF>
Thermo<TF>::~Thermo()
{
}

template<typename TF>
Thermo_type Thermo<TF>::get_switch()
{
    return swthermo;
}

template<typename TF>
std::shared_ptr<Thermo<TF>> Thermo<TF>::factory(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin, const Sim_mode sim_mode)
{
    std::string swthermo = inputin.get_item<std::string>("thermo", "swthermo", "", "0");

    if (swthermo == "0")
        return std::make_shared<Thermo_disabled<TF>>(masterin, gridin, fieldsin, inputin);
    else if (swthermo == "dry")
        return std::make_shared<Thermo_dry<TF>>(masterin, gridin, fieldsin, inputin, sim_mode);
    else if (swthermo == "moist")
        return std::make_shared<Thermo_moist<TF>>(masterin, gridin, fieldsin, inputin, sim_mode);
    else if (swthermo == "buoy")
        return std::make_shared<Thermo_buoy<TF>>(masterin, gridin, fieldsin, inputin);
    else
    {
        std::string msg = swthermo + " is an illegal value for swthermo";
        throw std::runtime_error(msg);
    }
}


#ifdef FLOAT_SINGLE
template class Thermo<float>;
#else
template class Thermo<double>;
#endif
