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

#include "soil_field3d.h"
#include "soil_grid.h"
#include "grid.h"
#include "tools.h"

#ifdef USECUDA
template<typename TF>
void Soil_field3d<TF>::init_device()
{
    auto& agd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    const int ijmemsize = agd.ijcells * sizeof(TF);
    const int nmemsize = sgd.ncells * sizeof(TF);

    cuda_safe_call(cudaMalloc(&fld_g, nmemsize));
    cuda_safe_call(cudaMalloc(&fld_bot_g, ijmemsize));
    cuda_safe_call(cudaMalloc(&fld_top_g, ijmemsize));
    cuda_safe_call(cudaMalloc(&flux_bot_g, ijmemsize));
    cuda_safe_call(cudaMalloc(&flux_top_g, ijmemsize));
}

template<typename TF>
void Soil_field3d<TF>::clear_device()
{
    cuda_safe_call(cudaFree(fld_g));
    cuda_safe_call(cudaFree(fld_bot_g));
    cuda_safe_call(cudaFree(fld_top_g));
    cuda_safe_call(cudaFree(flux_bot_g));
    cuda_safe_call(cudaFree(flux_top_g));
}
#endif


#ifdef FLOAT_SINGLE
template class Soil_field3d<float>;
#else
template class Soil_field3d<double>;
#endif
