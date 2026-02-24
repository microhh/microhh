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


#include <stdexcept>
#include "canopy.h"
#include "field3d_operators.h"
#include "field3d_io.h"

#ifdef USECUDA
template <typename TF>
void Canopy<TF>::exec()
{
    if (!sw_canopy)
        return;

    throw std::runtime_error("Canopy is not (yet) implemented on the GPU.");
}

template <typename TF>
void Canopy<TF>::prepare_device()
{
    throw std::runtime_error("Canopy is not (yet) implemented on the GPU.");
}

template <typename TF>
void Canopy<TF>::clear_device()
{
    throw std::runtime_error("Canopy is not (yet) implemented on the GPU.");
}
#endif

#ifdef FLOAT_SINGLE
template class Canopy<float>;
#else
template class Canopy<double>;
#endif
