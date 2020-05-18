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

#include "grid.h"
#include "fields.h"
#include "constants.h"
#include "master.h"
#include "stats.h"

#include "advec.h"
#include "advec_disabled.h"

template<typename TF>
Advec_disabled<TF>::Advec_disabled(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Advec<TF>(masterin, gridin, fieldsin, inputin)
{
}

template<typename TF>
Advec_disabled<TF>::~Advec_disabled() {}

template<typename TF>
unsigned long Advec_disabled<TF>::get_time_limit(unsigned long idt, const double dt)
{
    return Constants::ulhuge;
}

template<typename TF>
double Advec_disabled<TF>::get_cfl(const double dt)
{
    return cflmin;
}

template<typename TF>
void Advec_disabled<TF>::create(Stats<TF>& stats)
{
}

template<typename TF>
void Advec_disabled<TF>::exec(Stats<TF>&) {}

template<typename TF>
void Advec_disabled<TF>::get_advec_flux(Field3d<TF>& advec_flux, const Field3d<TF>& fld)
{
    std::fill(advec_flux.fld.begin(), advec_flux.fld.end(), TF(0.));
}

template class Advec_disabled<double>;
template class Advec_disabled<float>;
