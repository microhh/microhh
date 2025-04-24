/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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
    template<typename TF> __global__
    void get_projected_sum(
            TF* const __restrict__ fldxy,
            const TF* const __restrict__ fld,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kcells,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij = i + j*icells;
            fldxy[ij] = 0;
            for (int k=kstart; k<kcells; ++k)
            {
                const int ijk = ij + k*ijcells;
                fldxy[ij] += fld[ijk];
            }
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
    // Reduce 3D field excluding ghost cells and padding to jtot*kcells values
    TF* tmp = grid.get_tmp_3d_g();

    reduce_interior<TF>(
        fld, tmp, gd.imax, gd.istart, gd.iend, gd.jmax,
        gd.jstart, gd.jend, gd.kcells, 0, gd.icells, gd.ijcells, Sum_type);

    // Reduce jtot*kcells to kcells values
    reduce_all<TF>(
        tmp, prof, gd.jmax*gd.kcells, gd.kcells, gd.jmax, Sum_type, scalefac);

    grid.release_tmp_3d_g(tmp);
}

template<typename TF>
void Field3d_operators<TF>::calc_sum_profile_g(TF* const restrict prof, const TF* const restrict fld)
{
    using namespace Tools_g;

    const auto& gd = grid.get_grid_data();
    // const TF scalefac = 1./(gd.itot*gd.jtot);

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
        prof, fld, 1,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kcells,
        gd.icells, gd.ijcells);
    */

    // Optimized reduction method. This gives slightly different results compared to the CPU...
    // Reduce 3D field excluding ghost cells and padding to jtot*kcells values
    TF* tmp = grid.get_tmp_3d_g();

    reduce_interior<TF>(
        fld, tmp, gd.imax, gd.istart, gd.iend, gd.jmax,
        gd.jstart, gd.jend, gd.kcells, 0, gd.icells, gd.ijcells, Sum_type);

    // Reduce jtot*kcells to kcells values
    reduce_all<TF>(
        tmp, prof, gd.jmax*gd.kcells, gd.kcells, gd.jmax, Sum_type, 1.);

    grid.release_tmp_3d_g(tmp);
}

template<typename TF>
TF Field3d_operators<TF>::calc_mean_2d_g(const TF* const restrict fld)
{
    const auto& gd = grid.get_grid_data();
    const TF scalefac = 1./(gd.itot*gd.jtot);

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU(gridi, gridj);
    dim3 blockGPU(blocki, blockj);

    TF mean_value = 0;
    TF* mean_value_g;
    cuda_safe_call(cudaMalloc(&mean_value_g, sizeof(TF)));
    cudaMemcpy(mean_value_g, &mean_value, sizeof(TF), cudaMemcpyHostToDevice);

    // Very naive reduction from imax*jmax to single value.
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
TF Field3d_operators<TF>::calc_sum_2d_g(const TF* const restrict fld)
{
    const auto& gd = grid.get_grid_data();
    // const TF scalefac = 1./(gd.itot*gd.jtot);

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU(gridi, gridj);
    dim3 blockGPU(blocki, blockj);

    TF sum_value = 0;
    TF* sum_value_g;
    cuda_safe_call(cudaMalloc(&sum_value_g, sizeof(TF)));
    cudaMemcpy(sum_value_g, &sum_value, sizeof(TF), cudaMemcpyHostToDevice);

    // Very naive reduction from imax*jmax to single value.
    get_mean_2d<<<gridGPU, blockGPU>>>(
        sum_value_g, fld, TF(1.),
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.icells);

    cudaDeviceSynchronize();

    cudaMemcpy(&sum_value, sum_value_g, sizeof(TF), cudaMemcpyDeviceToHost);
    cudaFree(sum_value_g);

    return sum_value;
}

template<typename TF>
TF Field3d_operators<TF>::calc_mean_g(const TF* const restrict fld)
{
    using namespace Tools_g;

    const Grid_data<TF>& gd = grid.get_grid_data();
    const TF scalefac = 1./(gd.itot*gd.jtot*gd.zsize);
    TF mean_value;

    TF* tmp = grid.get_tmp_3d_g();

    // Reduce 3D field excluding ghost cells and padding to jtot*ktot values
    reduce_interior<TF>(fld, tmp, gd.imax, gd.istart, gd.iend, gd.jmax, gd.jstart, gd.jend, gd.ktot, gd.kstart, gd.icells, gd.ijcells, Sum_type);
    // Reduce jtot*ktot to ktot values
    for (int k=0; k<gd.ktot; ++k)
    {
        reduce_all<TF> (&tmp[gd.jmax*k], &tmp[gd.jmax*gd.ktot+k], gd.jmax, 1., gd.jmax, Sum_type, gd.dz[k+gd.kstart]);
    }
    // Reduce ktot values to a single value
    reduce_all<TF>     (&tmp[gd.jmax*gd.ktot], tmp, gd.ktot, 1, gd.ktot, Sum_type, scalefac);
    // Copy back result from GPU
    cuda_safe_call(cudaMemcpy(&mean_value, tmp, sizeof(TF), cudaMemcpyDeviceToHost));

    grid.release_tmp_3d_g(tmp);

    return mean_value;
}


template<typename TF>
TF Field3d_operators<TF>::calc_sum_g(const TF* const restrict fld)
{
    using namespace Tools_g;

    const Grid_data<TF>& gd = grid.get_grid_data();
    // const TF scalefac = 1./(gd.itot*gd.jtot*gd.zsize);
    TF sum_value;

    TF* tmp = grid.get_tmp_3d_g();

    // Reduce 3D field excluding ghost cells and padding to jtot*ktot values
    reduce_interior<TF>(fld, tmp, gd.imax, gd.istart, gd.iend, gd.jmax, gd.jstart, gd.jend, gd.ktot, gd.kstart, gd.icells, gd.ijcells, Sum_type);
    // Reduce jtot*ktot to ktot values
    for (int k=0; k<gd.ktot; ++k)
    {
        reduce_all<TF> (&tmp[gd.jmax*k], &tmp[gd.jmax*gd.ktot+k], gd.jmax, 1., gd.jmax, Sum_type, 1);
    }
    // Reduce ktot values to a single value
    reduce_all<TF>     (&tmp[gd.jmax*gd.ktot], tmp, gd.ktot, 1, gd.ktot, Sum_type, 1);
    // Copy back result from GPU
    cuda_safe_call(cudaMemcpy(&sum_value, tmp, sizeof(TF), cudaMemcpyDeviceToHost));

    grid.release_tmp_3d_g(tmp);

    return sum_value;
}

template<typename TF>
TF Field3d_operators<TF>::calc_max_g(const TF* const restrict fld)
{
    using namespace Tools_g;

    const Grid_data<TF>& gd = grid.get_grid_data();
    const TF scalefac = 1.;
    TF max_value;

    TF* tmp = grid.get_tmp_3d_g();

    // Reduce 3D field excluding ghost cells and padding to jtot*ktot values
    reduce_interior<TF>(fld, tmp, gd.imax, gd.istart, gd.iend, gd.jmax, gd.jstart, gd.jend, gd.ktot, gd.kstart, gd.icells, gd.ijcells, Max_type);
    // Reduce jtot*ktot to ktot values
    reduce_all<TF>     (tmp, &tmp[gd.jmax*gd.ktot], gd.jmax*gd.ktot, gd.ktot, gd.jmax, Max_type, scalefac);
    // Reduce ktot values to a single value
    reduce_all<TF>     (&tmp[gd.jmax*gd.ktot], tmp, gd.ktot, 1, gd.ktot, Max_type, scalefac);
    // Copy back result from GPU
    cuda_safe_call(cudaMemcpy(&max_value, tmp, sizeof(TF), cudaMemcpyDeviceToHost));

    grid.release_tmp_3d_g(tmp);

    return max_value;
}

template<typename TF>
void Field3d_operators<TF>::calc_proj_sum_g(TF* const restrict fldxy, const TF* const restrict fld)
{
    using namespace Tools_g;

    const Grid_data<TF>& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU(gridi, gridj);
    dim3 blockGPU(blocki, blockj);

    get_projected_sum<<<gridGPU, blockGPU>>>(
        fldxy, fld,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kcells,
        gd.icells, gd.ijcells);
}

#endif


#ifdef FLOAT_SINGLE
template class Field3d_operators<float>;
#else
template class Field3d_operators<double>;
#endif
