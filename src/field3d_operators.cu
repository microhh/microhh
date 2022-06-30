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

#include <cstdio>
#include <iostream>
#include <cmath>
#include "master.h"
#include "grid.h"
#include "field3d.h"
#include "fields.h"
#include "defines.h"
#include "tools.h"
#include "field3d_operators.h"

namespace
{
    template<typename TF> __global__
    void get_mean_profile(
            TF* const __restrict__ prof,
            const TF* const __restrict__ fld,
            const TF scalefac,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kcells,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z;

        if (i < iend && j < jend && k < kcells)
        {
            const int ijk = i + j*icells + k*ijcells;
            atomicAdd(&prof[k], fld[ijk]*scalefac);
        }
    }

    template<typename TF> __global__
    void get_mean_2d(
            TF* value, const TF* const __restrict__ fld,
            const TF scalefac,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij = i + j*icells;
            atomicAdd(value, fld[ij]*scalefac);
        }
    }
}

#ifdef USECUDA
template<typename TF>
void Field3d_operators<TF>::calc_mean_profile_g(TF* const restrict prof, const TF* const restrict fld)
{
    using namespace Tools_g;

    const auto& gd = grid.get_grid_data();
    const TF scalefac = 1./(gd.itot*gd.jtot);

    /*
    // Naive reduction with `atomicAdd`; slow, but this one gives perfect identical
    // results compared to the CPU reduction, so it can be helpful for debugging...
    const int blocki_1d = 128;
    const int gridi_1d  = gd.kcells/blocki_1d + (gd.kcells%blocki_1d > 0);
    dim3 gridGPU_1d (gridi_1d, 1, 1);
    dim3 blockGPU_1d(blocki_1d, 1, 1);

    set_to_val<<<gridGPU_1d, blockGPU_1d>>>(prof, gd.kcells, TF(0));

    const int blocki_3d = gd.ithread_block;
    const int blockj_3d = gd.jthread_block;
    const int gridi_3d  = gd.imax/blocki_3d + (gd.imax%blocki_3d > 0);
    const int gridj_3d  = gd.jmax/blockj_3d + (gd.jmax%blockj_3d > 0);

    dim3 gridGPU_3d(gridi_3d, gridj_3d, gd.kcells);
    dim3 blockGPU_3d(blocki_3d, blockj_3d, 1);

    get_mean_profile<<<gridGPU_3d, blockGPU_3d>>>(
        prof, fld, scalefac,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kcells,
        gd.icells, gd.ijcells);
    */

    // Optimized reduction method. This gives slightly different results compared to the CPU...
    auto tmp = fields.get_tmp_g();

    // Reduce 3D field excluding ghost cells and padding to jtot*kcells values
    reduce_interior<TF>(
        fld, tmp->fld_g, gd.itot, gd.istart, gd.iend, gd.jtot,
        gd.jstart, gd.jend, gd.kcells, 0, gd.icells, gd.ijcells, Sum_type);

    // Reduce jtot*kcells to kcells values
    reduce_all<TF>(
        tmp->fld_g, prof, gd.jtot*gd.kcells, gd.kcells, gd.jtot, Sum_type, scalefac);

    fields.release_tmp_g(tmp);
}

template<typename TF>
TF Field3d_operators<TF>::calc_mean_2d_g(const TF* const restrict fld)
{
    const auto& gd = grid.get_grid_data();
    const TF scalefac = 1./(gd.itot*gd.jtot);

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.itot/blocki + (gd.itot%blocki > 0);
    const int gridj  = gd.jtot/blockj + (gd.jtot%blockj > 0);

    dim3 gridGPU(gridi, gridj);
    dim3 blockGPU(blocki, blockj);

    TF mean_value = 0;
    TF* mean_value_g;
    cuda_safe_call(cudaMalloc(&mean_value_g, sizeof(TF)));
    cudaMemcpy(mean_value_g, &mean_value, sizeof(TF), cudaMemcpyHostToDevice);

    // Very naive reduction from itot*jtot to single value.
    get_mean_2d<<<gridGPU, blockGPU>>>(
        mean_value_g, fld, scalefac,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.icells);

    cudaDeviceSynchronize();

    cudaMemcpy(&mean_value, mean_value_g, sizeof(TF), cudaMemcpyDeviceToHost);
    cudaFree(mean_value_g);

    return mean_value;
}

template<typename TF>
TF Field3d_operators<TF>::calc_mean_g(const TF* const restrict fld)
{
    using namespace Tools_g;

    const Grid_data<TF>& gd = grid.get_grid_data();
    const TF scalefac = 1./(gd.itot*gd.jtot*gd.zsize);
    TF mean_value;

    auto tmp = fields.get_tmp_g();

    // Reduce 3D field excluding ghost cells and padding to jtot*ktot values
    reduce_interior<TF>(fld, tmp->fld_g, gd.itot, gd.istart, gd.iend, gd.jtot, gd.jstart, gd.jend, gd.ktot, gd.kstart, gd.icells, gd.ijcells, Sum_type);
    // Reduce jtot*ktot to ktot values
    for (int k=0; k<gd.ktot; ++k)
    {
        reduce_all<TF> (&tmp->fld_g[gd.jtot*k], &tmp->fld_g[gd.jtot*gd.ktot+k], gd.jtot, 1., gd.jtot, Sum_type, gd.dz[k+gd.kstart]);
    }
    // Reduce ktot values to a single value
    reduce_all<TF>     (&tmp->fld_g[gd.jtot*gd.ktot], tmp->fld_g, gd.ktot, 1, gd.ktot, Sum_type, scalefac);
    // Copy back result from GPU
    cuda_safe_call(cudaMemcpy(&mean_value, tmp->fld_g, sizeof(TF), cudaMemcpyDeviceToHost));

    fields.release_tmp_g(tmp);
    return mean_value;
}

template<typename TF>
TF Field3d_operators<TF>::calc_max_g(const TF* const restrict fld)
{
    using namespace Tools_g;

    const Grid_data<TF>& gd = grid.get_grid_data();
    const TF scalefac = 1.;
    TF max_value;

    auto tmp = fields.get_tmp_g();

    // Reduce 3D field excluding ghost cells and padding to jtot*ktot values
    reduce_interior<TF>(fld, tmp->fld_g, gd.itot, gd.istart, gd.iend, gd.jtot, gd.jstart, gd.jend, gd.ktot, gd.kstart, gd.icells, gd.ijcells, Max_type);
    // Reduce jtot*ktot to ktot values
    reduce_all<TF>     (tmp->fld_g, &tmp->fld_g[gd.jtot*gd.ktot], gd.jtot*gd.ktot, gd.ktot, gd.jtot, Max_type, scalefac);
    // Reduce ktot values to a single value
    reduce_all<TF>     (&tmp->fld_g[gd.jtot*gd.ktot], tmp->fld_g, gd.ktot, 1, gd.ktot, Max_type, scalefac);
    // Copy back result from GPU
    cuda_safe_call(cudaMemcpy(&max_value, tmp->fld_g, sizeof(TF), cudaMemcpyDeviceToHost));

    fields.release_tmp_g(tmp);

    return max_value;
}
#endif

template class Field3d_operators<double>;
template class Field3d_operators<float>;
