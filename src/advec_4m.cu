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

#include "advec_4m.h"
#include "grid.h"
#include "fields.h"
#include "tools.h"
#include "constants.h"
#include "finite_difference.h"
#include "field3d_operators.h"

#ifdef USECUDA
template<typename TF>
unsigned long Advec_4m<TF>::get_time_limit(unsigned long idt, double dt)
{
    throw std::runtime_error("Advec_4m not yet implemented on GPU\n");
}

template<typename TF>
double Advec_4m<TF>::get_cfl(const double dt)
{
    throw std::runtime_error("Advec_4m not yet implemented on GPU\n");
}

template<typename TF>
void Advec_4m<TF>::exec()
{
    throw std::runtime_error("Advec_4m not yet implemented on GPU\n");
}
#endif

template class Advec_4m<double>;
template class Advec_4m<float>;
