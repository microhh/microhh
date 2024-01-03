/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

#include "advec_2i6.h"
//#include "advec_2i6_kernels.cuh"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "tools.h"
#include "constants.h"
#include "finite_difference.h"
#include "field3d_operators.h"
#include "cuda_launcher.h"


#ifdef USECUDA
template<typename TF>
unsigned long Advec_2i6<TF>::get_time_limit(unsigned long idt, double dt)
{
    throw std::runtime_error("advec_2i6 is not (yet) implemented on the GPU.");
}


template<typename TF>
double Advec_2i6<TF>::get_cfl(const double dt)
{
    throw std::runtime_error("advec_2i6 is not (yet) implemented on the GPU.");
}


template<typename TF>
void Advec_2i6<TF>::exec(Stats<TF>& stats)
{
    throw std::runtime_error("advec_2i6 is not (yet) implemented on the GPU.");
}
#endif


#ifdef FLOAT_SINGLE
template class Advec_2i6<float>;
#else
template class Advec_2i6<double>;
#endif
