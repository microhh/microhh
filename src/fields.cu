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
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>
#include "data_block.h"
#include "fields.h"
#include "grid.h"
#include "master.h"
#include "constants.h"
#include "tools.h"

namespace
{
    // TODO use interp2 functions instead of manual interpolation
    template<typename TF> __global__
    void calc_mom_2nd_g(TF* __restrict__ u, TF* __restrict__ v, TF* __restrict__ w,
                        TF* __restrict__ mom, TF* __restrict__ dz,
                        int istart, int jstart, int kstart,
                        int iend,   int jend,   int kend,
                        int jj,     int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;
        const int ii = 1;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            mom[ijk] = (0.5*(u[ijk]+u[ijk+ii]) + 0.5*(v[ijk]+v[ijk+jj]) + 0.5*(w[ijk]+w[ijk+kk]))*dz[k];
        }
    }

    template<typename TF> __global__
    void calc_tke_2nd_g(TF* __restrict__ u, TF* __restrict__ v, TF* __restrict__ w,
                        TF* __restrict__ tke, TF* __restrict__ dz,
                        int istart, int jstart, int kstart,
                        int iend,   int jend,   int kend,
                        int jj,     int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;
        const int ii = 1;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            tke[ijk] = ( 0.5*(pow(u[ijk],2)+pow(u[ijk+ii],2))
                       + 0.5*(pow(v[ijk],2)+pow(v[ijk+jj],2))
                       + 0.5*(pow(w[ijk],2)+pow(w[ijk+kk],2)))*dz[k];
        }
    }

    template<typename TF> __global__
    void calc_mass_2nd_g(TF* __restrict__ s, TF* __restrict__ mass, TF* __restrict__ dz,
                         int istart, int jstart, int kstart,
                         int iend,   int jend,   int kend,
                         int jj,     int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            mass[ijk] = s[ijk]*dz[k];
        }
    }
}

#ifdef USECUDA
template<typename TF>
TF Fields<TF>::check_momentum()
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    auto tmp1 = get_tmp_g();

    calc_mom_2nd_g<<<gridGPU, blockGPU>>>(
        mp["u"]->fld_g, mp["v"]->fld_g, mp["w"]->fld_g, 
        tmp1->fld_g, gd.dz_g,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    TF mom = field3d_operators.calc_mean(tmp1->fld_g);

    release_tmp_g(tmp1);

    return mom;
}
#endif

#ifdef USECUDA
template<typename TF>
TF Fields<TF>::check_tke()
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    auto tmp1 = get_tmp_g();

    calc_tke_2nd_g<<<gridGPU, blockGPU>>>(
        mp["u"]->fld_g, mp["v"]->fld_g, mp["w"]->fld_g, 
        tmp1->fld_g, gd.dz_g,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    TF tke = field3d_operators.calc_mean(tmp1->fld_g);
    tke *= 0.5;

    release_tmp_g(tmp1);

    return tke;
}
#endif

#ifdef USECUDA
template<typename TF>
TF Fields<TF>::check_mass()
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    TF mass;

    // CvH for now, do the mass check on the first scalar... Do we want to change this?
    auto it = sp.begin();
    if (sp.begin() != sp.end())
    {
        auto tmp1 = get_tmp_g();

        calc_mass_2nd_g<<<gridGPU, blockGPU>>>(
            it->second->fld_g, tmp1->fld_g, gd.dz_g,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.ijcells);
        cuda_check_error();

        mass = field3d_operators.calc_mean(tmp1->fld_g);

        release_tmp_g(tmp1);
    }
    else
        mass = 0; 

    return mass;
}
#endif

/**
 * This function allocates all field3d instances and fields at device
 */
template<typename TF>
void Fields<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();
    const int nmemsize   = gd.ncells*sizeof(TF);
    const int nmemsize1d = gd.kcells*sizeof(TF);

    // Prognostic fields
    for (auto& it : a)
        it.second->init_device();

    // Tendencies
    for (auto& it : at)
        cuda_safe_call(cudaMalloc(&it.second->fld_g, nmemsize));

    // Reference profiles
    cuda_safe_call(cudaMalloc(&rhoref_g,  nmemsize1d));
    cuda_safe_call(cudaMalloc(&rhorefh_g, nmemsize1d));

    // copy all the data to the GPU
    forward_device();
}

/**
 * This function deallocates all field3d instances and fields at device
 */
template<typename TF>
void Fields<TF>::clear_device()
{
    for (auto& it : a)
        it.second->clear_device();

    for (auto& it : at)
        cuda_safe_call(cudaFree(it.second->fld_g));

    cuda_safe_call(cudaFree(rhoref_g));
    cuda_safe_call(cudaFree(rhorefh_g));

    // Free the tmp fields
    for (auto& it : atmp_g)
        it->clear_device();
}

/**
 * This function copies all fields from host to device
 */
template<typename TF>
void Fields<TF>::forward_device()
{
    auto& gd = grid.get_grid_data();
    for (auto& it : a)
        forward_field3d_device(it.second.get());

    for (auto& it : at)
        forward_field_device_3d(it.second->fld_g, it.second->fld.data());

    forward_field_device_1d(rhoref_g,  rhoref.data() , gd.kcells);
    forward_field_device_1d(rhorefh_g, rhorefh.data(), gd.kcells);
}

/**
 * This function copies all fields required for statistics and output from device to host
 */
template<typename TF>
void Fields<TF>::backward_device()
{
    for (auto& it : a)
        backward_field3d_device(it.second.get());
}

/* BvS: it would make more sense to put this routine in field3d.cu, but how to solve this with the calls to fields.cu? */
/**
 * This function copies a field3d instance from host to device
 * @param fld Pointer to field3d instance
 */
template<typename TF>
void Fields<TF>::forward_field3d_device(Field3d<TF>* fld)
{
    auto& gd = grid.get_grid_data();
    forward_field_device_3d(fld->fld_g,      fld->fld.data());
    forward_field_device_2d(fld->fld_bot_g,  fld->fld_bot.data());
    forward_field_device_2d(fld->fld_top_g,  fld->fld_top.data());
    forward_field_device_2d(fld->grad_bot_g, fld->grad_bot.data());
    forward_field_device_2d(fld->grad_top_g, fld->grad_top.data());
    forward_field_device_2d(fld->flux_bot_g, fld->flux_bot.data());
    forward_field_device_2d(fld->flux_top_g, fld->flux_top.data());
    forward_field_device_1d(fld->fld_mean_g, fld->fld_mean.data(), gd.kcells);
}

/* BvS: it would make more sense to put this routine in field3d.cu, but how to solve this with the calls to fields.cu? */
/**
 * This function copies a field3d instance from device to host
 * @param fld Pointer to field3d instance
 */
template<typename TF>
void Fields<TF>::backward_field3d_device(Field3d<TF>* fld)
{
    auto& gd = grid.get_grid_data();
    backward_field_device_3d(fld->fld.data(),      fld->fld_g     );
    backward_field_device_2d(fld->fld_bot.data(),  fld->fld_bot_g );
    backward_field_device_2d(fld->fld_top.data(),  fld->fld_top_g );
    backward_field_device_2d(fld->grad_bot.data(), fld->grad_bot_g);
    backward_field_device_2d(fld->grad_top.data(), fld->grad_top_g);
    backward_field_device_2d(fld->flux_bot.data(), fld->flux_bot_g);
    backward_field_device_2d(fld->flux_top.data(), fld->flux_top_g);
    backward_field_device_1d(fld->fld_mean.data(), fld->fld_mean_g, gd.kcells);
}

/**
 * This function copies a single 3d field from host to device
 * @param field_g Pointer to 3d field at device
 * @param field Pointer to 3d field at host
 * @param sw Switch to align the host field to device memory
 */
template<typename TF>
void Fields<TF>::forward_field_device_3d(TF* field_g, TF* field)
{
    auto& gd = grid.get_grid_data();
    cuda_safe_call(cudaMemcpy(field_g, field, gd.ncells*sizeof(TF), cudaMemcpyHostToDevice));
}

/**
 * This function copies a single 2d field from host to device
 * @param field_g Pointer to 2d field at device
 * @param field Pointer to 2d field at host
 * @param sw Switch to align the host field to device memory
 */
template<typename TF>
void Fields<TF>::forward_field_device_2d(TF* field_g, TF* field)
{
    auto& gd = grid.get_grid_data();
    cuda_safe_call(cudaMemcpy(field_g, field, gd.ijcells*sizeof(TF), cudaMemcpyHostToDevice));
}

/**
 * This function copies an array from host to device
 * @param field_g Pointer array at device
 * @param field Pointer to array at host
 * @param ncells Number of (TF precision) values to copy
 */
template<typename TF>
void Fields<TF>::forward_field_device_1d(TF* field_g, TF* field, int ncells)
{
    cuda_safe_call(cudaMemcpy(field_g, field, ncells*sizeof(TF), cudaMemcpyHostToDevice));
}

/**
 * This function copies a single 3d field from device to host
 * @param field Pointer to 3d field at host
 * @param field_g Pointer to 3d field at device
 * @param sw Switch to align the host field to device memory
 */
template<typename TF>
void Fields<TF>::backward_field_device_3d(TF* field, TF* field_g)
{
    auto& gd = grid.get_grid_data();
    cuda_safe_call(cudaMemcpy(field, field_g, gd.ncells*sizeof(TF), cudaMemcpyDeviceToHost));
}

/**
 * This function copies a single 2d field from device to host
 * @param field Pointer to 2d field at host
 * @param field_g Pointer to 2d field at device
 * @param sw Switch to align the host field to device memory
 */
template<typename TF>
void Fields<TF>::backward_field_device_2d(TF* field, TF* field_g)
{
    auto& gd = grid.get_grid_data();
    cuda_safe_call(cudaMemcpy(field, field_g, gd.ijcells*sizeof(TF), cudaMemcpyDeviceToHost));
}

/**
 * This function copies an array from device to host
 * @param field Pointer to array at host
 * @param field_g Pointer array at device
 * @param ncells Number of (TF precision) values to copy
 */
template<typename TF>
void Fields<TF>::backward_field_device_1d(TF* field, TF* field_g, int ncells)
{
    cuda_safe_call(cudaMemcpy(field, field_g, ncells*sizeof(TF), cudaMemcpyDeviceToHost));
}

template class Fields<double>;
template class Fields<float>;
