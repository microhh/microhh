/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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

#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "master.h"
#include "field3d_operators.h"

#include "source.h"
#include "source_disabled.h"
#include "source_gaussian.h"
#include "source_3d.h"

template<typename TF>
Source<TF>::Source(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& input) :
    master(masterin), grid(gridin), fields(fieldsin)
{
}

template<typename TF>
Source<TF>::~Source()
{
}

template<typename TF>
std::shared_ptr<Source<TF>> Source<TF>::factory(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin)
{
    std::string sw_source = inputin.get_item<std::string>("source", "swsource", "", "0");

    if (sw_source == "1")
    {
        masterin.print_warning("swsource=1 is deprecated. Defaulting to swsource=gaussian\n");
        sw_source = "gaussian";
    }

    if (sw_source == "0")
        return std::make_shared<Source_disabled<TF>>(masterin, gridin, fieldsin, inputin);
    else if (sw_source == "gaussian")
        return std::make_shared<Source_gaussian<TF>>(masterin, gridin, fieldsin, inputin);
    else if (sw_source == "3d")
        return std::make_shared<Source_3d<TF>>(masterin, gridin, fieldsin, inputin);
    else
        throw std::runtime_error("Illegal option for \"swsource\".");
}


#ifdef FLOAT_SINGLE
template class Source<float>;
#else
template class Source<double>;
#endif
