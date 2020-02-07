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

#include <cstdio>
#include <iostream>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "constants.h"

#include "soil.h"
#include "soil_disabled.h"


template<typename TF>
Soil_disabled<TF>::Soil_disabled(Master& masterin, Grid<TF>& gridin, Input& inputin): 
    Soil<TF>(masterin, gridin, inputin)
{
    sw_soil = Soil_type::Disabled;
}

template<typename TF>
Soil_disabled<TF>::~Soil_disabled()
{
}

template class Soil_disabled<double>;
template class Soil_disabled<float>;
