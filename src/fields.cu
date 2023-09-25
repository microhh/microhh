/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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
#include "fields.h"
#include "grid.h"
#include "soil_grid.h"
#include "master.h"
#include "column.h"
#include "constants.h"
#include "tools.h"
#include "fast_math.h"
#include "soil_field3d.h"

namespace
{
    namespace fm = Fast_math;

    // TODO use interp2 functions instead of manual interpolation
    template<typename TF> __global__
    void calc_mom_2nd_g(const TF* __restrict__ u, const TF* __restrict__ v, const TF* __restrict__ w,
                        TF* __restrict__ mom, const TF* __restrict__ dz,
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
            mom[ijk] = (TF(0.5)*(u[ijk]+u[ijk+ii])
                      + TF(0.5)*(v[ijk]+v[ijk+jj])
                      + TF(0.5)*(w[ijk]+w[ijk+kk]))*dz[k];
        }
    }

    template<typename TF> __global__
    void calc_tke_2nd_g(const TF* __restrict__ u, const TF* __restrict__ v, const TF* __restrict__ w,
                        TF* __restrict__ tke, const TF* __restrict__ dz,
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
            tke[ijk] = ( TF(0.5)*(fm::pow2(u[ijk])+fm::pow2(u[ijk+ii]))
                       + TF(0.5)*(fm::pow2(v[ijk])+fm::pow2(v[ijk+jj]))
                       + TF(0.5)*(fm::pow2(w[ijk])+fm::pow2(w[ijk+kk])))*dz[k];
        }
    }

    template<typename TF> __global__
    void set_to_val_g(
            TF* const __restrict__ fld, const TF value,
            const int icells, const int jcells, const int kcells,
            const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;
        const int k = blockIdx.z;

        if (i < icells && j < jcells && k < kcells)
        {
            const int ijk = i + j*icells + k*ijcells;
            fld[ijk] = value;
        }
    }

    template<typename TF> __global__
    void set_to_val_g(
            TF* const __restrict__ fld, const TF value,
            const int icells, const int jcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij = i + j*icells;
            fld[ij] = value;
        }
    }

    template<typename TF> __global__
    void set_to_val_g(
            TF* const __restrict__ fld, const TF value, const int kcells)
    {
        const int k = blockIdx.z;

        if (k < kcells)
            fld[k] = value;
    }

    template<typename TF>
    void set_field3d_to_value(
            Field3d<TF>* fld, const TF value, const Grid_data<TF>& gd)
    {
        const int blocki = gd.ithread_block;
        const int blockj = gd.jthread_block;
        const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
        const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

        dim3 gridGPU2(gridi, gridj, 1);
        dim3 gridGPU3(gridi, gridj, gd.kcells);

        dim3 blockGPU(blocki, blockj, 1);

        //dim3 gridGPU1(1, 1, 1);
        //dim3 blockGPU1(1, 1, gd.kcells);

        set_to_val_g<TF><<<gridGPU3, blockGPU>>>(
                fld->fld_g, value,
                gd.icells, gd.jcells, gd.kcells, gd.ijcells);
        cuda_check_error();

        set_to_val_g<TF><<<gridGPU2, blockGPU>>>(fld->fld_bot_g, value, gd.icells, gd.jcells);
        cuda_check_error();
        set_to_val_g<TF><<<gridGPU2, blockGPU>>>(fld->fld_top_g, value, gd.icells, gd.jcells);
        cuda_check_error();

        set_to_val_g<TF><<<gridGPU2, blockGPU>>>(fld->flux_bot_g, value, gd.icells, gd.jcells);
        cuda_check_error();
        set_to_val_g<TF><<<gridGPU2, blockGPU>>>(fld->flux_top_g, value, gd.icells, gd.jcells);
        cuda_check_error();

        set_to_val_g<TF><<<gridGPU2, blockGPU>>>(fld->grad_bot_g, value, gd.icells, gd.jcells);
        cuda_check_error();
        set_to_val_g<TF><<<gridGPU2, blockGPU>>>(fld->grad_top_g, value, gd.icells, gd.jcells);
        cuda_check_error();

        //set_to_val_g<<<gridGPU1, blockGPU1>>>(fld->fld_mean_g, value, gd.kcells);
        //cuda_check_error();
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

    calc_mom_2nd_g<TF><<<gridGPU, blockGPU>>>(
        mp.at("u")->fld_g, mp.at("v")->fld_g, mp.at("w")->fld_g,
        tmp1->fld_g, gd.dz_g,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    TF mom = field3d_operators.calc_mean_g(tmp1->fld_g);

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

    calc_tke_2nd_g<TF><<<gridGPU, blockGPU>>>(
        mp.at("u")->fld_g, mp.at("v")->fld_g, mp.at("w")->fld_g,
        tmp1->fld_g, gd.dz_g,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    TF tke = field3d_operators.calc_mean_g(tmp1->fld_g);
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
        mass = field3d_operators.calc_mean_g(it->second->fld_g);
    else
        mass = 0;

    return mass;
}

template<typename TF>
void Fields<TF>::exec()
{
    // calculate the means for the prognostic scalars
    if (calc_mean_profs)
    {
        for (auto& it : ap)
            field3d_operators.calc_mean_profile_g(it.second->fld_mean_g, it.second->fld_g);
    }
}
#endif

#ifdef USECUDA
template<typename TF>
std::shared_ptr<Field3d<TF>> Fields<TF>::get_tmp_g()
{
    auto& gd = grid.get_grid_data();
    std::shared_ptr<Field3d<TF>> tmp;

    #pragma omp critical
    {
        cudaDeviceSynchronize();

        // In case of insufficient tmp fields, allocate a new one.
        if (atmp_g.empty())
        {
            init_tmp_field_g();
            tmp = atmp_g.back();
            tmp->init_device();

            // Initialise all values at zero
            set_field3d_to_value(tmp.get(), TF(0), gd);
        }
        else
            tmp = atmp_g.back();

        atmp_g.pop_back();
    }

    // Assign to a huge negative number in case of debug mode.
    #ifdef __CUDACC_DEBUG__
    set_field3d_to_value(tmp.get(), TF(-1e30), gd);
    #endif

    return tmp;
}

template<typename TF>
void Fields<TF>::release_tmp_g(std::shared_ptr<Field3d<TF>>& tmp)
{
    #pragma omp critical
    {
        cudaDeviceSynchronize();
        atmp_g.push_back(std::move(tmp));
    }
}
#endif

/**
 * This function allocates all field3d instances and fields at device
 */
template<typename TF>
void Fields<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    const int ijmemsize  = gd.ijcells*sizeof(TF);
    const int nmemsize   = gd.ncells*sizeof(TF);
    const int nmemsize1d = gd.kcells*sizeof(TF);
    const int nmemsize_soil = sgd.ncells*sizeof(TF);

    // Prognostic fields atmosphere
    for (auto& it : a)
        it.second->init_device();

    // Tendencies atmosphere
    for (auto& it : at)
        it.second->fld_g.resize(gd.ncells);

    // Prognostic 2D fields
    for (auto& it : ap2d)
        cuda_safe_call(cudaMalloc(&it.second->fld_g, ijmemsize));

    // Tendencies 2D fields
    for (auto& it : at2d)
        cuda_safe_call(cudaMalloc(&it.second->fld_g, ijmemsize));

    // Prognostic fields soil
    for (auto& it : sps)
        it.second->init_device();

    // Tendencies soil
    for (auto& it : sts)
        cuda_safe_call(cudaMalloc(&it.second->fld_g, nmemsize_soil));

    // Reference profiles
    rhoref_g.resize(gd.kcells);
    rhorefh_g.resize(gd.kcells);
    rhorefi_g.resize(gd.kcells);
    rhorefhi_g.resize(gd.kcells);

    // copy all the data to the GPU
    forward_device();
}

/**
 * This function deallocates all field3d instances and fields at device
 */
template<typename TF>
void Fields<TF>::clear_device()
{
    // Prognostic fields atmosphere
    for (auto& it : a)
        it.second->clear_device();

    // Tendencies atmosphere
    for (auto& it : at)
        it.second->fld_g.free();

    // Prognostic 2D fields
    for (auto& it : ap2d)
        cuda_safe_call(cudaFree(it.second->fld_g));

    // Tendencies 2D fields
    for (auto& it : at2d)
        cuda_safe_call(cudaFree(it.second->fld_g));

    // Prognostic fields soil
    for (auto& it : sps)
        it.second->clear_device();

    // Tendencies soil
    for (auto& it : sts)
        cuda_safe_call(cudaFree(it.second->fld_g));

    rhoref_g.free();
    rhorefh_g.free();
    rhorefi_g.free();
    rhorefhi_g.free();

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
    auto& sgd = soil_grid.get_grid_data();

    // Prognostic fields atmosphere
    for (auto& it : a)
        forward_field3d_device(it.second.get());

    // Tendencies atmosphere
    for (auto& it : at)
        forward_field_device_3d(it.second->fld_g, it.second->fld.data());

    // Prognostic 2D fields
    for (auto& it : ap2d)
        forward_field_device_2d(it.second->fld_g, it.second->fld.data());

    // Tendencies 2D fields
    for (auto& it : at2d)
        forward_field_device_2d(it.second->fld_g, it.second->fld.data());

    // Prognostic fields soil
    for (auto& it : sps)
        forward_soil_field3d_device(it.second.get());

    // Tendencies soil
    for (auto& it : sts)
        forward_field_device(it.second->fld_g, it.second->fld.data(), sgd.ncells);

    forward_field_device(rhoref_g,   rhoref.data() , gd.kcells);
    forward_field_device(rhorefh_g,  rhorefh.data(), gd.kcells);

    // Calculate reciprocal of rho
    std::vector<TF> rhorefi(gd.kcells);
    std::vector<TF> rhorefhi(gd.kcells);

    for (int k = 0; k < gd.kcells; k++) {
        rhorefi[k] = 1.0 / rhoref[k];
        rhorefhi[k] = 1.0 / rhorefh[k];
    }

    forward_field_device(rhorefi_g,  rhorefi.data() , gd.kcells);
    forward_field_device(rhorefhi_g, rhorefhi.data(), gd.kcells);
}

/**
 * This function copies all fields required for statistics and output from device to host
 */
template<typename TF>
void Fields<TF>::backward_device()
{
    auto& gd = grid.get_grid_data();

    // Prognostic fields atmosphere
    for (auto& it : a)
        backward_field3d_device(it.second.get());

    // Prognostic 2D fields
    for (auto& it : ap2d)
        backward_field_device_2d(it.second->fld.data(), it.second->fld_g);

    // Prognostic fields soil.
    for (auto& it : sps)
        backward_soil_field3d_device(it.second.get());
}

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
    forward_field_device(fld->fld_mean_g, fld->fld_mean.data(), gd.kcells);
}

/**
 * This function copies a Soil_field3d instance from host to device
 * @param fld Pointer to Soil_field3d instance
 */
template<typename TF>
void Fields<TF>::forward_soil_field3d_device(Soil_field3d<TF>* fld)
{
    auto& agd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    const int nmemsize  = sgd.ncells * sizeof(TF);
    const int ijmemsize = agd.ijcells * sizeof(TF);

    cuda_safe_call(cudaMemcpy(fld->fld_g,      fld->fld.data(),      nmemsize,  cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(fld->fld_bot_g,  fld->fld_bot.data(),  ijmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(fld->fld_top_g,  fld->fld_top.data(),  ijmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(fld->flux_bot_g, fld->flux_bot.data(), ijmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(fld->flux_top_g, fld->flux_top.data(), ijmemsize, cudaMemcpyHostToDevice));
}

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
    backward_field_device(fld->fld_mean.data(), fld->fld_mean_g, gd.kcells);
}

/**
 * This function copies a Soil_field3d instance from device to host
 * @param fld Pointer to Soil_field3d instance
 */
template<typename TF>
void Fields<TF>::backward_soil_field3d_device(Soil_field3d<TF>* fld)
{
    auto& agd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    const int nmemsize = sgd.ncells * sizeof(TF);
    const int ijmemsize = agd.ijcells * sizeof(TF);

    cuda_safe_call(cudaMemcpy(fld->fld.data(),      fld->fld_g,      nmemsize, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(fld->fld_bot.data(),  fld->fld_bot_g,  ijmemsize, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(fld->fld_top.data(),  fld->fld_top_g,  ijmemsize, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(fld->flux_bot.data(), fld->flux_bot_g, ijmemsize, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(fld->flux_top.data(), fld->flux_top_g, ijmemsize, cudaMemcpyDeviceToHost));
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
void Fields<TF>::forward_field_device(TF* field_g, TF* field, int ncells)
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
void Fields<TF>::backward_field_device(TF* field, TF* field_g, int ncells)
{
    cuda_safe_call(cudaMemcpy(field, field_g, ncells*sizeof(TF), cudaMemcpyDeviceToHost));
}

#ifdef USECUDA
template<typename TF>
void Fields<TF>::exec_column(Column<TF>& column)
{
    auto& gd = grid.get_grid_data();

    const TF no_offset = 0.;

    column.calc_column("u",mp.at("u")->fld_g, gd.utrans);
    column.calc_column("v",mp.at("v")->fld_g, gd.vtrans);
    column.calc_column("w",mp.at("w")->fld_g, no_offset);

    for (auto& it : sp)
    {
        column.calc_column(it.first, it.second->fld_g, no_offset);
    }

    column.calc_column("p", sd.at("p")->fld_g, no_offset);
}
#endif


#ifdef FLOAT_SINGLE
template class Fields<float>;
#else
template class Fields<double>;
#endif
