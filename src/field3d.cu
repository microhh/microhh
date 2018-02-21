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

#include <map>
#include <vector>
#include "field3d.h"
#include "grid.h"
#include "master.h"
#include "tools.h"

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
void Field3d<TF>::calc_mean_profile(TF * tmp)
{
    using namespace Tools_g;
    const Grid_data<TF>& gd = grid.get_grid_data();
    const TF scalefac = 1./(gd.itot*gd.jtot);

    // Reduce 3D field excluding ghost cells and padding to jtot*kcells values
    reduce_interior<TF>(fld.data(), tmp, gd.itot, gd.istart, gd.iend, gd.jtot, gd.jstart, gd.jend, gd.kcells, 0, gd.icellsp, gd.ijcellsp, sumType);
    // Reduce jtot*kcells to kcells values
    reduce_all<TF>     (tmp, fld_mean.data(), gd.jtot*gd.kcells, gd.kcells, gd.jtot, sumType, scalefac);
}

template<typename TF>
TF Field3d<TF>::calc_mean(TF * tmp)
{

    using namespace Tools_g;
    const Grid_data<TF>& gd = grid.get_grid_data();
    const TF scalefac = 1./(gd.itot*gd.jtot*gd.ktot);
    TF sumvalue;

    // Reduce 3D field excluding ghost cells and padding to jtot*ktot values
    reduce_interior<TF>(fld.data(), tmp, gd.itot, gd.istart, gd.iend, gd.jtot, gd.jstart, gd.jend, gd.kcells, 0, gd.icellsp, gd.ijcellsp, sumType);
    // Reduce jtot*ktot to ktot values
    reduce_all<TF>     (tmp, &tmp[gd.jtot*gd.ktot], gd.jtot*gd.ktot, gd.ktot, gd.jtot, sumType, scalefac);
    // Reduce ktot values to a single value
    reduce_all<TF>     (&tmp[gd.jtot*gd.ktot], tmp, gd.ktot, 1, gd.ktot, sumType, tmp[0]);
    // Copy back result from GPU
    cuda_safe_call(cudaMemcpy(&sumvalue, &tmp[0], sizeof(TF), cudaMemcpyDeviceToHost));

    return sumvalue;
}

#endif

template class Field3d<double>;
template class Field3d<float>;
