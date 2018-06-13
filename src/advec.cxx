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

#include <memory>
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "constants.h"
#include "master.h"
#include "field3d_operators.h"

#include "advec.h"
#include "advec_disabled.h"
#include "advec_2.h"
#include "advec_2i3.h"
#include "advec_2i4.h"
// #include "advec_4.h"
// #include "advec_4m.h"

template<typename TF>
Advec<TF>::Advec(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& input) :
    master(masterin), grid(gridin), fields(fieldsin), field3d_operators(master, grid, fields),
    cflmin(1.E-5)
{
    cflmax = input.get_item<TF>("advec", "cflmax", "", 1.);
    swadvec = "0";
}

template<typename TF>
Advec<TF>::~Advec()
{
}

template<typename TF>
std::shared_ptr<Advec<TF>> Advec<TF>::factory(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin)
{
    std::string swspatialorder = (gridin.get_spatial_order() == Grid_order::Second) ? "2nd" : "4th";
    std::string swadvec = inputin.get_item<std::string>("advec", "swadvec", "", swspatialorder);

    if (swadvec == "0")
        return std::make_shared<Advec_disabled<TF>>(masterin, gridin, fieldsin, inputin);
    else if (swadvec == "2")
        return std::make_shared<Advec_2<TF>>(masterin, gridin, fieldsin, inputin);
    else if (swadvec == "2i3")
        return std::make_shared<Advec_2i3<TF>>(masterin, gridin, fieldsin, inputin);
    else if (swadvec == "2i4")
        return std::make_shared<Advec_2i4<TF>>(masterin, gridin, fieldsin, inputin);
    // else if (swadvec == "4")
    //     return new Advec_4(modelin, inputin);
    // else if (swadvec == "4m")
    //     return new Advec_4m(modelin, inputin);
    else
    {
        masterin.print_error("\"%s\" is an illegal value for swadvec\n", swadvec.c_str());
        throw std::runtime_error("Illegal options swadvec");
    }
}

template<typename TF>
std::string Advec<TF>::get_switch()
{
    return swadvec;
}

template class Advec<double>;
template class Advec<float>;
