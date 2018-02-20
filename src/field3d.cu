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
#include "data_block.h"

template<typename TF>
void Field3d<TF>::release_cuda_fields() 
{
    cuda_safe_call(cudaFreeHost(data));
    cuda_safe_call(cudaFreeHost(databot));
    cuda_safe_call(cudaFreeHost(datatop));
    cuda_safe_call(cudaFreeHost(datagradbot));
    cuda_safe_call(cudaFreeHost(datagradtop));
    cuda_safe_call(cudaFreeHost(datafluxbot));
    cuda_safe_call(cudaFreeHost(datafluxtop));
    cuda_safe_call(cudaFreeHost(datamean));
}

template<typename TF>
void Field3d<TF>::init_cuda()
{
    const Grid_data<TF>& gd = grid.get_grid_data();
    
    const int ijksize = gd.ncells *sizeof(TF);
    const int ijsize  = gd.ijcells*sizeof(TF);
    const int ksize   = gd.kcells *sizeof(TF);

    // Allocate the 3d field.
    cuda_safe_call(cudaMallocHost(&data, ijksize));

    // Allocate the boundary cells.
    cuda_safe_call(cudaMallocHost(&databot, ijsize));
    cuda_safe_call(cudaMallocHost(&datatop, ijsize));
    cuda_safe_call(cudaMallocHost(&datagradbot, ijsize));
    cuda_safe_call(cudaMallocHost(&datagradtop, ijsize));
    cuda_safe_call(cudaMallocHost(&datafluxbot, ijsize));
    cuda_safe_call(cudaMallocHost(&datafluxtop, ijsize));
    cuda_safe_call(cudaMallocHost(&datamean, ksize));
}

template<typename TF>
void Field3d<TF>::init_device()
{
    const int nmemsize   = grid->ncellsp*sizeof(TF);
    const int nmemsize1d = grid->kcells *sizeof(TF);
    const int nmemsize2d = (grid->ijcellsp+grid->memoffset)*sizeof(TF);

    cuda_safe_call(cudaMalloc(&fld_g,        nmemsize  ));
    cuda_safe_call(cudaMalloc(&dfld_bot_g,     nmemsize2d));
    cuda_safe_call(cudaMalloc(&fld_top_g,     nmemsize2d));
    cuda_safe_call(cudaMalloc(&dgrad_bot_g, nmemsize2d));
    cuda_safe_call(cudaMalloc(&grad_top_g, nmemsize2d));
    cuda_safe_call(cudaMalloc(&flux_bot_g, nmemsize2d));
    cuda_safe_call(cudaMalloc(&flux_top_g, nmemsize2d));
    cuda_safe_call(cudaMalloc(&datamean_g,    nmemsize1d));
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
    cuda_safe_call(cudaFree(datamean_g));
}

#ifdef USECUDA
template<typename TF>
void Field3d<TF>::calc_mean_profile()
{
    using namespace Tools_g;

    const TF scalefac = 1./(itot*jtot);

    // Reduce 3D field excluding ghost cells and padding to jtot*kcells values
    reduce_interior(data, tmp, itot, istart, iend, jtot, jstart, jend, kcells, 0, icellsp, ijcellsp, sumType);
    // Reduce jtot*kcells to kcells values
    reduce_all     (tmp, prof, jtot*kcells, kcells, jtot, sumType, scalefac);
}

template<typename TF>
TF Field3d<TF>::calc_mean()
{
}
#endif

template class Field3d<double>;
template class Field3d<float>;
