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
#include "fields.h"
#include "input.h"

// Microphysics schemes
#include "microphys.h"
#include "microphys_disabled.h"


template<typename TF>
Microphys<TF>::Microphys(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& input) :
    master(masterin), grid(gridin), fields(fieldsin)
{
}

template<typename TF>
Microphys<TF>::~Microphys()
{
}

template<typename TF>
std::string Microphys<TF>::get_switch()
{
    return swmicro;
}

template<typename TF>
std::shared_ptr<Microphys<TF>> Microphys<TF>::factory(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin)
{
    std::string swmicro = inputin.get_item<std::string>("micro", "swmicro", "", "0");

    if (swmicro == "0")
        return std::make_shared<Microphys_disabled<TF>>(masterin, gridin, fieldsin, inputin);
    //else if (swthermo == "2mom_warm")
    //    return std::make_shared<Micro_2mom_warm<TF>>(masterin, gridin, fieldsin, inputin);
    else
    {
        masterin.print_error("\"%s\" is an illegal value for \"swmicro\"\n", swmicro.c_str());
        throw std::runtime_error("Illegal options \"swmicro\"");
    }
}

template class Microphys<double>;
template class Microphys<float>;
