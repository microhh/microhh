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

#include "grid.h"
#include "fields.h"
#include "master.h"

#include "source.h"
#include "source_disabled.h"


template<typename TF>
Source_disabled<TF>::Source_disabled(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Source<TF>(masterin, gridin, fieldsin, inputin)
{
}

template<typename TF>
Source_disabled<TF>::~Source_disabled()
{
}

template<typename TF>
void Source_disabled<TF>::init()
{
}

template<typename TF>
void Source_disabled<TF>::create(Input& input, Netcdf_handle& nc_in)
{
}

template<typename TF>
void Source_disabled<TF>::exec()
{
}

template<typename TF>
void Source_disabled<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
}

#ifdef USECUDA
template<typename TF>
void Source_disabled<TF>::prepare_device()
{
}
#endif

#ifdef FLOAT_SINGLE
template class Source_disabled<float>;
#else
template class Source_disabled<double>;
#endif
