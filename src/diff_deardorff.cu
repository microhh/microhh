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

#include <algorithm>
#include <cmath>
#include <iostream>

#include "grid.h"
#include "fields.h"
#include "master.h"
#include "thermo.h"
#include "boundary.h"
#include "stats.h"

#include "diff_deardorff.h"
#include "diff_kernels.h"

namespace
{
}

#ifdef USECUDA
template<typename TF>
unsigned long Diff_deardorff<TF>::get_time_limit(const unsigned long idt, const double dt)
{
}

template<typename TF>
double Diff_deardorff<TF>::get_dn(const double dt)
{
}

template<typename TF>
void Diff_deardorff<TF>::exec(Stats<TF>& stats)
{
}

template<typename TF>
void Diff_deardorff<TF>::exec_viscosity(Stats<TF>& stats, Thermo<TF>& thermo)
{
}

template<typename TF>
void Diff_deardorff<TF>::prepare_device(Boundary<TF>& boundary)
{
}

template<typename TF>
void Diff_deardorff<TF>::clear_device()
{
}
#endif

template class Diff_deardorff<double>;
template class Diff_deardorff<float>;
