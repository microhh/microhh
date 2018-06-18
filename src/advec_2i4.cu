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

#include "advec_2i4.h"
#include "grid.h"
#include "fields.h"
#include "constants.h"
#include "tools.h"
#include "finite_difference.h"
#include "field3d_operators.h"

using namespace Finite_difference::O2;

#ifdef USECUDA
template<typename TF>
unsigned long Advec_2i4<TF>::get_time_limit(unsigned long idt, double dt)
{
    throw std::runtime_error("swadvec=2i4 not (yet) implemented on GPU");
}

template<typename TF>
double Advec_2i4<TF>::get_cfl(const double dt)
{
    throw std::runtime_error("swadvec=2i4 not (yet) implemented on GPU");
}

template<typename TF>
void Advec_2i4<TF>::exec()
{
    throw std::runtime_error("swadvec=2i4 not (yet) implemented on GPU");
}
#endif

template class Advec_2i4<double>;
template class Advec_2i4<float>;
