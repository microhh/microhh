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

#include "grid.h"
#include "tools.h"
#include "math.h"

template<typename TF>
void Grid<TF>::prepare_device()
{
    // Calculate optimal size thread blocks based on grid
    gd.ithread_block = min(256, 16 * ((gd.itot / 16) + (gd.itot % 16 > 0)));
    gd.jthread_block = 256 / gd.ithread_block;

    const int imemsize = gd.icells*sizeof(TF);
    const int jmemsize = gd.jcells*sizeof(TF);
    const int kmemsize = gd.kcells*sizeof(TF);

    cuda_safe_call(cudaMalloc((void**)&gd.x_g,     imemsize));
    cuda_safe_call(cudaMalloc((void**)&gd.y_g,     jmemsize));
    cuda_safe_call(cudaMalloc((void**)&gd.z_g,     kmemsize));
    cuda_safe_call(cudaMalloc((void**)&gd.zh_g,    kmemsize));
    cuda_safe_call(cudaMalloc((void**)&gd.dz_g,    kmemsize));
    cuda_safe_call(cudaMalloc((void**)&gd.dzh_g,   kmemsize));
    cuda_safe_call(cudaMalloc((void**)&gd.dzi_g,   kmemsize));
    cuda_safe_call(cudaMalloc((void**)&gd.dzhi_g,  kmemsize));
    cuda_safe_call(cudaMalloc((void**)&gd.dzi4_g,  kmemsize));
    cuda_safe_call(cudaMalloc((void**)&gd.dzhi4_g, kmemsize));

    cuda_safe_call(cudaMemcpy(gd.x_g,     gd.x.data(),     imemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.y_g,     gd.y.data(),     jmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.z_g,     gd.z.data(),     kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.zh_g,    gd.zh.data(),    kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.dz_g,    gd.dz.data(),    kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.dzh_g,   gd.dzh.data(),   kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.dzi_g,   gd.dzi.data(),   kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.dzhi_g,  gd.dzhi.data(),  kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.dzi4_g,  gd.dzi4.data(),  kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gd.dzhi4_g, gd.dzhi4.data(), kmemsize, cudaMemcpyHostToDevice));
}

template<typename TF>
void Grid<TF>::clear_device()
{
    cuda_safe_call(cudaFree(gd.y_g    ));
    cuda_safe_call(cudaFree(gd.z_g    ));
    cuda_safe_call(cudaFree(gd.zh_g   ));
    cuda_safe_call(cudaFree(gd.dz_g   ));
    cuda_safe_call(cudaFree(gd.dzh_g  ));
    cuda_safe_call(cudaFree(gd.dzi_g  ));
    cuda_safe_call(cudaFree(gd.dzhi_g ));
    cuda_safe_call(cudaFree(gd.dzi4_g ));
    cuda_safe_call(cudaFree(gd.dzhi4_g));
}

template class Grid<double>;
template class Grid<float>;
