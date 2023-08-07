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

#include <stdexcept>

#include "radiation_disabled.h"
#include "constants.h"

template<typename TF>
Radiation_disabled<TF>::Radiation_disabled(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Radiation<TF>(masterin, gridin, fieldsin, inputin)
{
    swradiation = "0";
}

template<typename TF>
Radiation_disabled<TF>::~Radiation_disabled()
{
}

template<typename TF>
unsigned long Radiation_disabled<TF>::get_time_limit(unsigned long itime)
{
    return Constants::ulhuge;
}

template<typename TF>
bool Radiation_disabled<TF>::check_field_exists(const std::string& name)
{
    return false;  // always returns error
}

#ifdef FLOAT_SINGLE
template class Radiation_disabled<float>;
#else
template class Radiation_disabled<double>;
#endif
