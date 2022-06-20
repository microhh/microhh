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

#include <map>
#include <vector>
#include "field3d.h"
#include "grid.h"
#include "master.h"
#include "tools.h"
#include <iostream>

template<typename TF>
void Field3d<TF>::init_device()
{
    const Grid_data<TF>& gd = grid.get_grid_data();
    const int nmemsize   = gd.ncells  * sizeof(TF);
    const int nmemsize1d = gd.kcells  * sizeof(TF);
    const int nmemsize2d = gd.ijcells * sizeof(TF);

    cuda_safe_call(cudaMalloc(&fld_g,      nmemsize  ));
    cuda_safe_call(cudaMalloc(&fld_bot_g,  nmemsize2d));
    cuda_safe_call(cudaMalloc(&fld_top_g,  nmemsize2d));
    cuda_safe_call(cudaMalloc(&grad_bot_g, nmemsize2d));
    cuda_safe_call(cudaMalloc(&grad_top_g, nmemsize2d));
    cuda_safe_call(cudaMalloc(&flux_bot_g, nmemsize2d));
    cuda_safe_call(cudaMalloc(&flux_top_g, nmemsize2d));
    cuda_safe_call(cudaMalloc(&fld_mean_g, nmemsize1d));
}

template<typename TF>
void Field3d<TF>::clear_device()
{
    cuda_safe_call(cudaFree(fld_g));
    cuda_safe_call(cudaFree(fld_bot_g));
    cuda_safe_call(cudaFree(fld_top_g));
    cuda_safe_call(cudaFree(grad_bot_g));
    cuda_safe_call(cudaFree(grad_top_g));
    cuda_safe_call(cudaFree(flux_bot_g));
    cuda_safe_call(cudaFree(flux_top_g));
    cuda_safe_call(cudaFree(fld_mean_g));
}

template class Field3d<double>;
template class Field3d<float>;
