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

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fftw3.h>
#include <cufft.h>
#include <iostream>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "pres.h"
#include "pres_4.h"
#include "defines.h"
#include "tools.h"
#include "constants.h"
#include "field3d_operators.h"

#ifdef USECUDA
template<typename TF>
void Pres_4<TF>::prepare_device()
{
    throw std::runtime_error("Pres_4 not yet implemented on GPU\n");
}

template<typename TF>
void Pres_4<TF>::clear_device()
{
    throw std::runtime_error("Pres_4 not yet implemented on GPU\n");
}

template<typename TF>
void Pres_4<TF>::exec(double dt)
{
    throw std::runtime_error("Pres_4 not yet implemented on GPU\n");
}

template<typename TF>
TF Pres_4<TF>::check_divergence()
{
    throw std::runtime_error("Pres_4 not yet implemented on GPU\n");
}
#endif

template class Pres_4<double>;
template class Pres_4<float>;
