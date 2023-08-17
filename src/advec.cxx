/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
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

#include <memory>
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "constants.h"
#include "master.h"
#include "field3d_operators.h"
#include "stats.h"

#include "advec.h"
#include "advec_disabled.h"
#include "advec_2.h"
#include "advec_2i4.h"
#include "advec_2i5.h"
#include "advec_2i62.h"
#include "advec_2i53.h"
#include "advec_4.h"
#include "advec_4m.h"

template<typename TF>
Advec<TF>::Advec(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& input) :
    master(masterin), grid(gridin), fields(fieldsin), field3d_operators(master, grid, fields),
    cflmin(1.E-5)
{
    cflmax = input.get_item<TF>("advec", "cflmax", "", 1.);
}

template<typename TF>
Advec<TF>::~Advec()
{
}

template<typename TF>
std::shared_ptr<Advec<TF>> Advec<TF>::factory(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin)
{
    std::string swspatialorder = (gridin.get_spatial_order() == Grid_order::Second) ? "2" : "4";
    std::string swadvec = inputin.get_item<std::string>("advec", "swadvec", "", swspatialorder);

    if (swadvec == "0")
        return std::make_shared<Advec_disabled<TF>>(masterin, gridin, fieldsin, inputin);
    else if (gridin.get_spatial_order() == Grid_order::Second)
    {
        if (swadvec == "2")
            return std::make_shared<Advec_2<TF>>(masterin, gridin, fieldsin, inputin);
        else if (swadvec == "2i4")
            return std::make_shared<Advec_2i4<TF>>(masterin, gridin, fieldsin, inputin);
        else if (swadvec == "2i5")
            return std::make_shared<Advec_2i5<TF>>(masterin, gridin, fieldsin, inputin);
        else if (swadvec == "2i53")
            return std::make_shared<Advec_2i53<TF>>(masterin, gridin, fieldsin, inputin);
        else if (swadvec == "2i62")
            return std::make_shared<Advec_2i62<TF>>(masterin, gridin, fieldsin, inputin);
        else
        {
            std::string msg = "swadvec = \"" + swadvec +  "\" is an illegal value with swspatialorder = \"2\"";
            throw std::runtime_error(msg);
        }
    }
    else
    {
        if (swadvec == "4")
            return std::make_shared<Advec_4<TF>>(masterin, gridin, fieldsin, inputin);
        else if (swadvec == "4m")
            return std::make_shared<Advec_4m<TF>>(masterin, gridin, fieldsin, inputin);
        else
        {
            std::string msg = "swadvec = \"" + swadvec +  "\" is an illegal value with swspatialorder = \"4\"";
            throw std::runtime_error(msg);
        }
    }
}

template class Advec<double>;
template class Advec<float>;
