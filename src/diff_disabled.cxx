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

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "constants.h"

// Diffusion schemes
#include "diff.h"
#include "diff_disabled.h"
//#include "diff_2.h"
//#include "diff_4.h"
//#include "diff_smag2.h"


template<typename TF>
Diff_disabled<TF>::Diff_disabled(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Diff<TF>(masterin, gridin, fieldsin, inputin)
{
}

template <typename TF>
Diff_disabled<TF>::~Diff_disabled()
{
}

template <typename TF>
Diffusion_type Diff_disabled<TF>::get_switch() const
{
    return swdiff;
}

template<typename TF>
unsigned long Diff_disabled<TF>::get_time_limit(const unsigned long idtlim, const double dt)
{
    return Constants::ulhuge;
}

template<typename TF>
double Diff_disabled<TF>::get_dn(const double dt)
{
    return Constants::dsmall;
}

template class Diff_disabled<double>;
template class Diff_disabled<float>;
