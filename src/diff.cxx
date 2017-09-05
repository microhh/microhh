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
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "grid.h"
#include "fields.h"
#include "master.h"
#include "defines.h"
#include "constants.h"
#include "model.h"

// diffusion schemes
#include "diff.h"
#include "diff_disabled.h"
//#include "diff_2.h"
//#include "diff_4.h"
//#include "diff_smag2.h"

template<typename TF>
Diff<TF>::Diff(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& input) :
    master(masterin), grid(gridin), fields(fieldsin)
{
}

template <typename TF>
Diff<TF>::~Diff()
{
}

template<typename TF>
std::shared_ptr<Diff<TF>> Diff<TF>::factory(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin, const std::string swspatialorder)
{
    std::string swdiff     = inputin.get_item<std::string>("diff",     "swdiff",     "", swspatialorder);
    std::string swboundary = inputin.get_item<std::string>("boundary", "swboundary", "", "default");

    if (swdiff == "0")
        return std::make_shared<Diff_disabled<TF>>(masterin, gridin, fieldsin, inputin);
    else
    {
        masterin.print_error("\"%s\" is an illegal value for swdiff\n", swdiff.c_str());
        throw std::runtime_error("Illegal options swdiff");
    }
}

template class Diff<double>;
template class Diff<float>;




//Diff* Diff::factory(Master* masterin, Input* inputin, Model* modelin, const std::string swspatialorder)
//{
//    std::string swdiff;
//    std::string swboundary;
//
//    int nerror = 0;
//    nerror += inputin->get_item(&swdiff, "diff", "swdiff", "", swspatialorder);
//    // load the boundary switch as well in order to be able to check whether the surface model is used
//    nerror += inputin->get_item(&swboundary, "boundary", "swboundary", "", "default");
//    if (nerror)
//        return 0;
//
//    if (swdiff == "0")
//        return new Diff_disabled(modelin, inputin);
//    else if (swdiff == "2")
//        return new Diff_2(modelin, inputin);
//    else if (swdiff == "4")
//        return new Diff_4(modelin, inputin);
//    else if (swdiff == "smag2")
//    {
//        // the subgrid model requires a surface model because of the MO matching at first level
//        if ((swboundary == "surface") || (swboundary == "surface_bulk") || (swboundary == "surface_patch"))
//            return new Diff_smag_2(modelin, inputin);
//        else
//        {
//            masterin->print_error("swdiff=\"smag2\" requires a surface model (swboundary = \"surface\", \"surface_bulk\" or \"surface_patch\")\n");
//            throw 1;
//        }
//    }
//    else
//    {
//        masterin->print_error("\"%s\" is an illegal value for swdiff\n", swdiff.c_str());
//        throw 1;
//    }
//}
//
//std::string Diff::get_switch()
//{
//    return swdiff;
//}
