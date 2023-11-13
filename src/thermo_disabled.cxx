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
#include "master.h"
#include "grid.h"
#include "fields.h"

#include "constants.h"
#include "thermo.h"
#include "thermo_disabled.h"

template<typename TF>
Thermo_disabled<TF>::Thermo_disabled(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) : Thermo<TF>(masterin, gridin, fieldsin, inputin)
{
    swthermo = Thermo_type::Disabled;
}

template<typename TF>
Thermo_disabled<TF>::~Thermo_disabled()
{
}

template<typename TF>
unsigned long Thermo_disabled<TF>::get_time_limit(unsigned long idt, const double dt)
{
    return Constants::ulhuge;
}

template<typename TF>
bool Thermo_disabled<TF>::check_field_exists(std::string name)
{
    return false;  // always returns error
}

template<typename TF>
TF Thermo_disabled<TF>::get_buoyancy_diffusivity()
{
    return 0.;
}


#ifdef FLOAT_SINGLE
template class Thermo_disabled<float>;
#else
template class Thermo_disabled<double>;
#endif
