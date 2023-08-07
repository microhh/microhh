/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
 * Copyright (c) 2018-2019 Elynn Wu
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
#include <stdexcept>

#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"

#include "radiation.h"
#include "radiation_disabled.h"
#include "radiation_gcss.h"
#include "radiation_rrtmgp.h"
#include "radiation_prescribed.h"

#include "Optical_props.h"

template<typename TF>
Radiation<TF>::Radiation(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), field3d_operators(master, grid, fields)
{
}

template<typename TF>
Radiation<TF>::~Radiation()
{
}

template<typename TF>
std::string Radiation<TF>::get_switch()
{
    return swradiation;
}

template<typename TF>
std::shared_ptr<Radiation<TF>> Radiation<TF>::factory(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin)
{
    std::string swradiation = inputin.get_item<std::string>("radiation", "swradiation", "", "0");
    if (swradiation == "0")
        return std::make_shared<Radiation_disabled<TF>>(masterin, gridin, fieldsin, inputin);
    else if (swradiation == "rrtmgp")
        return std::make_shared<Radiation_rrtmgp<TF>>(masterin, gridin, fieldsin, inputin);
    else if (swradiation == "gcss") // gcss - for Sc clouds.
        return std::make_shared<Radiation_gcss<TF>>(masterin, gridin, fieldsin, inputin);
    else if (swradiation == "prescribed")
        return std::make_shared<Radiation_prescribed<TF>>(masterin, gridin, fieldsin, inputin);
    else
    {
        std::string error_message = swradiation + " is an illegal value for swradiation";
        throw std::runtime_error(error_message);
    }
}

#ifdef FLOAT_SINGLE
template class Radiation<float>;
#else
template class Radiation<double>;
#endif

