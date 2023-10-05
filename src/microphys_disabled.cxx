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
#include "microphys.h"
#include "microphys_disabled.h"

template<typename TF>
Microphys_disabled<TF>::Microphys_disabled(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) : 
    Microphys<TF>(masterin, gridin, fieldsin, inputin)
{
    swmicrophys = Microphys_type::Disabled;

    Nc0 = inputin.get_item<TF>("micro", "Nc0", "", -1);
    Ni0 = inputin.get_item<TF>("micro", "Ni0", "", -1);
}

template<typename TF>
Microphys_disabled<TF>::~Microphys_disabled()
{
}

template<typename TF>
unsigned long Microphys_disabled<TF>::get_time_limit(unsigned long idt, const double dt)
{
    return Constants::ulhuge;
}

template<typename TF>
void Microphys_disabled<TF>::get_surface_rain_rate(std::vector<TF>& field)
{
    std::fill(field.begin(), field.end(), TF(0));
}

template<typename TF>
TF Microphys_disabled<TF>::get_Nc0()
{
    if (this->Nc0 > 0)
        return this->Nc0;
    else
        throw std::runtime_error("Requested uninitialised Nc0!");
}

template<typename TF>
TF Microphys_disabled<TF>::get_Ni0()
{
    if (this->Ni0 > 0)
        return this->Ni0;
    else
        throw std::runtime_error("Requested uninitialised Ni0!");
}

#ifdef FLOAT_SINGLE
template class Microphys_disabled<float>;
#else
template class Microphys_disabled<double>;
#endif
