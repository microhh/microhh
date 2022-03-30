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

#include <iostream>

#include "soil_grid.h"
#include "grid.h"
#include "tools.h"

#ifdef USECUDA
template<typename TF>
void Soil_grid<TF>::prepare_device()
{
    const int kmemsize = gd.kcells*sizeof(TF);
    const int khmemsize = gd.kcellsh*sizeof(TF);

    cuda_safe_call(cudaMalloc((void**)&gd.z_g,   kmemsize));
    cuda_safe_call(cudaMalloc((void**)&gd.dz_g,  kmemsize));
    cuda_safe_call(cudaMalloc((void**)&gd.dzi_g, kmemsize));

    cuda_safe_call(cudaMalloc((void**)&gd.zh_g,   khmemsize));
    cuda_safe_call(cudaMalloc((void**)&gd.dzh_g,  khmemsize));
    cuda_safe_call(cudaMalloc((void**)&gd.dzhi_g, khmemsize));

    cuda_safe_call(cudaMemcpy(gd.z_g,   gd.z.data(),   kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.dz_g,  gd.dz.data(),  kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.dzi_g, gd.dzi.data(), kmemsize, cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy(gd.zh_g,   gd.zh.data(),   khmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.dzh_g,  gd.dzh.data(),  khmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.dzhi_g, gd.dzhi.data(), khmemsize, cudaMemcpyHostToDevice));
}

template<typename TF>
void Soil_grid<TF>::clear_device()
{
    cuda_safe_call(cudaFree(gd.z_g));
    cuda_safe_call(cudaFree(gd.dz_g));
    cuda_safe_call(cudaFree(gd.dzi_g));

    cuda_safe_call(cudaFree(gd.zh_g));
    cuda_safe_call(cudaFree(gd.dzh_g));
    cuda_safe_call(cudaFree(gd.dzhi_g));
}
#endif

template class Soil_grid<double>;
template class Soil_grid<float>;
