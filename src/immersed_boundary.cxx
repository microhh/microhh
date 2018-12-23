/*
 * MicroHH
 * Copyright (c) 2011-2019 Chiel van Heerwaarden
 * Copyright (c) 2011-2019 Thijs Heus
 * Copyright (c) 2014-2019 Bart van Stratum
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

#include "immersed_boundary.h"

template<typename TF>
Immersed_boundary<TF>::Immersed_boundary(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
}

template <typename TF>
Immersed_boundary<TF>::~Immersed_boundary()
{
}

template <typename TF>
void Immersed_boundary<TF>::init(Input& inputin)
{
    std::string sw_ib_str = inputin.get_item<std::string>("immersed_boundary", "sw_ib", "", "0");

    if (sw_ib_str == "0")
        sw_ib = IB_type::disabled;
    else if (sw_ib_str == "dem")
        sw_ib = IB_type::dem;
    else if (sw_ib_str == "user")
        sw_ib = IB_type::user;
    else
    {
        std::string error = "\"" + sw_ib_str + "\" is an illegal value for \"sw_ib\"";
        throw std::runtime_error(error);
    }
}

template class Immersed_boundary<double>;
template class Immersed_boundary<float>;
