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

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"

#include "advec_2i53.h"

#ifdef USECUDA
template<typename TF>
double Advec_2i53<TF>::get_cfl(double dt)
{
    throw std::runtime_error("Advec_2i53 is not (yet) implemented on the GPU");
}

template<typename TF>
unsigned long Advec_2i53<TF>::get_time_limit(unsigned long idt, double dt)
{
    throw std::runtime_error("Advec_2i53 is not (yet) implemented on the GPU");
}

template<typename TF>
void Advec_2i53<TF>::exec(Stats<TF>& stats)
{
    throw std::runtime_error("Advec_2i53 is not (yet) implemented on the GPU");
}
#endif

template class Advec_2i53<double>;
template class Advec_2i53<float>;
