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
#include "pres_disabled.h"

template<typename TF>
Pres_disabled<TF>::Pres_disabled(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, FFT<TF>& fftin, Input& inputin) :
    Pres<TF>(masterin, gridin, fieldsin, fftin, inputin) {}

template<typename TF>
Pres_disabled<TF>::~Pres_disabled() {}

template<typename TF>
TF Pres_disabled<TF>::check_divergence()
{
    TF divmax = 0.;
    return divmax;
}

template<typename TF>
void Pres_disabled<TF>::init() {}

template<typename TF>
void Pres_disabled<TF>::create(Stats<TF>& stats) {}

template<typename TF>
void Pres_disabled<TF>::set_values() {}

template<typename TF>
void Pres_disabled<TF>::exec(const double dt, Stats<TF>& stats) {}

// BvS: mixing CUDA and CPU code; put in pres_disabled.cu? Might be a bit of overkill?
#ifdef USECUDA
template<typename TF>
void Pres_disabled<TF>::prepare_device() {}

template<typename TF>
void Pres_disabled<TF>::clear_device() {}
#endif


#ifdef FLOAT_SINGLE
template class Pres_disabled<float>;
#else
template class Pres_disabled<double>;
#endif
