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

#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "timeloop.h"
#include "constants.h"

#include "particle_lagrangian.h"

namespace
{
}


template<typename TF>
Particle_lagrangian<TF>::Particle_lagrangian(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    sw_particle = inputin.get_item<bool>("particle_lagrangian", "sw_particle", "", false);
}


template<typename TF>
Particle_lagrangian<TF>::~Particle_lagrangian()
{
}


template<typename TF>
void Particle_lagrangian<TF>::create(Timeloop<TF>& timeloop)
{
    if (!sw_particle)
        return;
}


template<typename TF>
unsigned long Particle_lagrangian<TF>::get_time_limit()
{
    if (!sw_particle)
        return Constants::ulhuge;
}


#ifndef USECUDA
template<typename TF>
void Particle_lagrangian<TF>::exec(Stats<TF>& stats)
{
    if (!sw_particle)
        return;
}
#endif


#ifdef FLOAT_SINGLE
template class Particle_lagrangian<float>;
#else
template class Particle_lagrangian<double>;
#endif