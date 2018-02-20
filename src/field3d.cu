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

#include "field3d.h"
#include "grid.h"
#include "master.h"
#include "tools.h"

template<typename TF>
void Field3d<TF>::release_cuda_fields() 
{
    cuda_safe_call(cudaFreeHost(fld.data()));
    cuda_safe_call(cudaFreeHost(fld_bot.data()));
    cuda_safe_call(cudaFreeHost(fld_top.data()));
    cuda_safe_call(cudaFreeHost(grad_bot.data()));
    cuda_safe_call(cudaFreeHost(grad_top.data()));
    cuda_safe_call(cudaFreeHost(flux_bot.data()));
    cuda_safe_call(cudaFreeHost(flux_top.data()));
    cuda_safe_call(cudaFreeHost(fld_mean.data()));
}

template<typename TF>
void Field3d<TF>::init_cuda()
{
    const Grid_data<TF>& gd = grid.get_grid_data();
    
    const int ijksize = gd.ncells *sizeof(TF);
    const int ijsize  = gd.ijcells*sizeof(TF);
    const int ksize   = gd.kcells *sizeof(TF);

    // Allocate the 3d field.
    cuda_safe_call(cudaMallocHost(&(fld.data()), ijksize));

    // Allocate the boundary cells.
    cuda_safe_call(cudaMallocHost(&(fld_bot.data()), ijsize)); //TH: Should we use the C++ api here, and ensure that it is pinned mem?
    cuda_safe_call(cudaMallocHost(&(fld_top.data()), ijsize)); 
    cuda_safe_call(cudaMallocHost(&(grad_bot.data()), ijsize)); 
    cuda_safe_call(cudaMallocHost(&(grad_top.data()), ijsize)); 
    cuda_safe_call(cudaMallocHost(&(flux_bot.data()), ijsize)); 
    cuda_safe_call(cudaMallocHost(&(flux_top.data()), ijsize)); 
    cuda_safe_call(cudaMallocHost(&(fld_mean.data()), ksize)); 
}

template<typename TF>
void Field3d<TF>::init_device()
{
    const Grid_data<TF>& gd = grid.get_grid_data();
    const int nmemsize   = gd.ncellsp*sizeof(TF);
    const int nmemsize1d = gd.kcells *sizeof(TF);
    const int nmemsize2d = (gd.ijcellsp+gd.memoffset)*sizeof(TF);

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

#ifdef USECUDA
template<typename TF>
void Field3d<TF>::calc_mean_profile()
{
    using namespace Tools_g;
    const Grid_data<TF>& gd = grid.get_grid_data();
    const TF scalefac = 1./(gd.itot*gd.jtot);

    // Reduce 3D field excluding ghost cells and padding to jtot*kcells values
    reduce_interior(data, tmp, gd.itot, gd.istart, gd.iend, gd.jtot, gd.jstart, gd.jend, gd.kcells, 0, gd.icellsp, gd.ijcellsp, sumType);
    // Reduce jtot*kcells to kcells values
    reduce_all     (tmp, fld_mean, gd.jtot*gd.kcells, gd.kcells, gd.jtot, sumType, scalefac);
}

template<typename TF>
TF Field3d<TF>::calc_mean()
{

    using namespace Tools_g;
    const Grid_data<TF>& gd = grid.get_grid_data();
    const TF scalefac = 1./(gd.itot*gd.jtot*gd.ktot);
    TF sumvalue;

    TF* tmp = atmp.at(tmp1)->fld.data();

    // Reduce 3D field excluding ghost cells and padding to jtot*ktot values
    reduce_interior(fld.data(), atmp.at(tmp1)->fld.data(), gd.itot, gd.istart, gd.iend, gd.jtot, gd.jstart, gd.jend, gd.kcells, 0, gd.icellsp, gd.ijcellsp, sumType);
    // Reduce jtot*ktot to ktot values
    reduce_all     (tmp, tmp[jtot*ktot], gd.jtot*gd.ktot, gd.ktot, gd.jtot, sumType, 1);
    // Reduce ktot values to a single value
    reduce_all     (tmp[jtot*ktot], tmp, ktot, 1, ktot, sumType, scalefac);
    // Copy back result from GPU
    cuda_safe_call(cudaMemcpy(&sumvalue, &tmp[0], sizeof(TF), cudaMemcpyDeviceToHost));

    return sumvalue;
}

#endif

