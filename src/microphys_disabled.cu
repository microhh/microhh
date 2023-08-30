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

#include "fields.h"
#include "microphys.h"
#include "microphys_disabled.h"
#include "tools.h"
#include "grid.h"

#ifdef USECUDA
template<typename TF>
void Microphys_disabled<TF>::get_surface_rain_rate_g(
    TF* const __restrict__ rr_g)
{
    auto& gd = grid.get_grid_data();
    cuda_safe_call(cudaMemset(rr_g, TF(0), gd.ijcells*sizeof(TF)));
}
#endif


#ifdef FLOAT_SINGLE
template class Microphys_disabled<float>;
#else
template class Microphys_disabled<double>;
#endif
