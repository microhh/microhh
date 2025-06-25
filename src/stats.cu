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

#include <cstdio>
#include <cmath>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <utility>
#include <netcdf.h>

#include "master.h"
#include "grid.h"
#include "soil_grid.h"
// #include "background_profs.h"
#include "fields.h"
#include "stats.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "timeloop.h"
#include "advec.h"
#include "diff.h"
#include "netcdf_interface.h"
#include "stats_functions.h"
#include "tools.h"
#include "grid_kernels.h"


using namespace Stats_functions;


#ifdef USECUDA
template<typename TF>
void Stats<TF>::calc_stats_mean(
        const std::string& varname, const Field3d<TF>& fld, const TF offset)
{
    auto& gd = grid.get_grid_data();

    unsigned int flag;
    const int* nmask;
    std::string name;
        cuda_check_error();
    auto masked = fields.get_tmp_g();
    // Calc mean of atmospheric variables
    if (std::find(varlist.begin(), varlist.end(), varname) != varlist.end())
    {
        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, fld.loc[2]);
            auto threads = grid.get_dim_gpu(gd.icells, gd.jcells, gd.kcells);
            apply_mask_g<<<threads.first, threads.second>>>(masked->fld_g.data(),  fld.fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);
            field3d_operators.calc_mean_profile_g(masked->fld_mean_g, masked->fld_g);
            cuda_safe_call(cudaMemcpy(m.second.profs.at(varname).data.data(), masked->fld_mean_g.data(), gd.kcells * sizeof(TF), cudaMemcpyDeviceToHost));
            for (int k=gd.kstart; k<gd.kend+1; ++k)
            {
                if (nmask[k])
                    m.second.profs.at(varname).data[k] *= gd.itot * gd.jtot / nmask[k];
            }
            master.sum(m.second.profs.at(varname).data.data(), gd.kcells);

            // Add the offset.
            for (auto& value : m.second.profs.at(varname).data)
                value += offset;
            set_fillvalue_prof(m.second.profs.at(varname).data.data(), nmask, gd.kstart, gd.kcells);
        }
    }
    fields.release_tmp_g(masked);

    name = varname + "_bot";
    calc_stats_2d_g(name, fld.fld_bot_g, offset);
}

template<typename TF>
void Stats<TF>::calc_stats_2d_g(
        const std::string& varname, const cuda_vector<TF>& fld, const TF offset)
{
    auto& gd = grid.get_grid_data();

    unsigned int flag;
    const int* nmask;
    std::string name;
    auto masked = fields.get_tmp_xy_g();
    // Calc mean of atmospheric variables
    if (std::find(varlist.begin(), varlist.end(), varname) != varlist.end())
    {
        for (auto& m : masks)
        {
            if (m.second.nmask_bot > 0)
            {

                set_flag(flag, nmask, m.second, 1);
                auto threads = grid.get_dim_gpu(gd.icells, gd.jcells, 1);

                apply_mask_g<<<threads.first, threads.second>>>(masked->data(),  fld.data(), mfield_g, flag, gd.icells, gd.jcells, 1, gd.ijcells);
                m.second.tseries.at(varname).data = field3d_operators.calc_sum_2d_g(masked->data())/m.second.nmask_bot;
                master.sum(&m.second.tseries.at(varname).data, 1);
                m.second.tseries.at(varname).data += offset;
            }
            else
                m.second.tseries.at(varname).data = netcdf_fp_fillvalue<TF>();

        }
    }
    fields.release_tmp_xy_g(masked);

}

#endif


#ifdef USECUDA

    __global__
    void set_to_val_g(
            unsigned int* const __restrict__ fld, const unsigned int value,
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
    // template<typename TF> __global__
    // void add_val(TF* __restrict__ a, const TF* const __restrict__ fld, int nsize, TF val)
    // {
    //     const int n = blockIdx.x*blockDim.x + threadIdx.x;

    //     if (n < nsize)
    //         a[n] = fld[n] + val;
    // }

    // template<typename TF> __global__
    // void raise_to_pow(TF* __restrict__ a, const int nsize, const int exponent)
    // {
    //     const int n = blockIdx.x*blockDim.x + threadIdx.x;

    //     if (n < nsize)
    //         a[n] = pow(a[n], exponent);
    // }

    template<typename TF> __global__
    void add_profile(
            TF* __restrict__ fld,
            const TF* __restrict__ profile,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if ((i < iend) && (j < jend) && (k < kend))
        {
            const int ijk = i + j*icells + k*ijcells;
            fld[ijk] += profile[k];
        }
    }

template<typename TF>
void Stats<TF>::calc_stats_moments(
        const std::string& varname, const Field3d<TF>& fld, const TF offset)
{
    using namespace Tools_g;
    auto& gd = grid.get_grid_data();

    unsigned int flag;
    const int* nmask;
    std::string name;
            cuda_check_error();

        cuda_check_error();
    auto masked = fields.get_tmp_g();
    auto dev = fields.get_tmp_g();
    // // Calc moments

    for (int power=2; power<=4; power++)
    {
        name = varname + "_" + std::to_string(power);

        if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
        {
            for (auto& m : masks)
            {
                set_flag(flag, nmask, m.second, fld.loc[2]);
                auto threads_1d = grid.get_dim_gpu(gd.ijcells);
                auto threads = grid.get_dim_gpu(gd.icells, gd.jcells, gd.kcells);

                for (int k = gd.kstart; k < gd.kend+1; ++k)
                {
                    int kk = k * gd.ijcells;
                    add_val<<<threads_1d.first, threads_1d.second>>>(&(dev->fld_g[kk]),&(fld.fld_g[kk]), gd.ijcells, - m.second.profs.at(varname).data[k] + offset);
                }
                threads_1d = grid.get_dim_gpu(gd.ncells);
                raise_to_pow<<<threads_1d.first, threads_1d.second>>>(dev->fld_g.data(), gd.ncells, power);

                apply_mask_g<<<threads.first, threads.second>>>(masked->fld_g.data(),  dev->fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);
                field3d_operators.calc_mean_profile_g(masked->fld_mean_g, masked->fld_g);
                cuda_safe_call(cudaMemcpy(m.second.profs.at(name).data.data(), masked->fld_mean_g.data(), gd.kcells * sizeof(TF), cudaMemcpyDeviceToHost));
                for (int k=gd.kstart; k<gd.kend+1; ++k)
                {
                    if (nmask[k])
                        m.second.profs.at(name).data[k] *= gd.itot * gd.jtot / nmask[k];
                }
                master.sum(m.second.profs.at(name).data.data(), gd.kcells);

                set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);
            }
        }
    }
    fields.release_tmp_g(masked);
    fields.release_tmp_g(dev);
}


template<typename TF>
void Stats<TF>::calc_stats_w(
        const std::string& varname, const Field3d<TF>& fld, const TF offset)
{
    using namespace Tools_g;
    auto& gd = grid.get_grid_data();

    unsigned int flag;
    const int* nmask;
    std::string name;
        cuda_check_error();
    auto masked = fields.get_tmp_g();
    auto advec_flux = fields.get_tmp_g();

    // Calc Resolved Flux
    name = varname + "_w";
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        advec.get_advec_flux(*advec_flux, fld);

        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, fld.loc[2]);
            auto threads = grid.get_dim_gpu(gd.icells, gd.jcells, gd.kcells);
            apply_mask_g<<<threads.first, threads.second>>>(masked->fld_g.data(),  advec_flux->fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);
            field3d_operators.calc_mean_profile_g(masked->fld_mean_g, masked->fld_g);
            cuda_safe_call(cudaMemcpy(m.second.profs.at(name).data.data(), masked->fld_mean_g.data(), gd.kcells * sizeof(TF), cudaMemcpyDeviceToHost));
            for (int k=gd.kstart; k<gd.kend+1; ++k)
            {
                if (nmask[k])
                    m.second.profs.at(name).data[k] *= gd.itot * gd.jtot / nmask[k];
            }
            master.sum(m.second.profs.at(name).data.data(), gd.kcells);

            set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);
        }
    }
    fields.release_tmp_g(advec_flux);
    fields.release_tmp_g(masked);

}


template<typename TF>
void Stats<TF>::calc_stats_diff(
        const std::string& varname, const Field3d<TF>& fld, const TF offset)
{
    using namespace Tools_g;
    auto& gd = grid.get_grid_data();

    unsigned int flag;
    const int* nmask;
    std::string name;
        cuda_check_error();
    auto masked = fields.get_tmp_g();
    auto diff_flux = fields.get_tmp_g();


    // Calc Diffusive Flux
    name = varname + "_diff";
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        diff.get_diff_flux(*diff_flux, fld);

        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, !fld.loc[2]);
            auto threads = grid.get_dim_gpu(gd.icells, gd.jcells, gd.kcells);
            apply_mask_g<<<threads.first, threads.second>>>(masked->fld_g.data(), diff_flux->fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);
            field3d_operators.calc_mean_profile_g(masked->fld_mean_g, masked->fld_g);
            cuda_safe_call(cudaMemcpy(m.second.profs.at(name).data.data(), masked->fld_mean_g.data(), gd.kcells * sizeof(TF), cudaMemcpyDeviceToHost));
            for (int k=gd.kstart; k<gd.kend+1; ++k)
            {
                if (nmask[k])
                    m.second.profs.at(name).data[k] *= gd.itot * gd.jtot / nmask[k];
            }
            master.sum(m.second.profs.at(name).data.data(), gd.kcells);

            set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);
        }
    }
    fields.release_tmp_g(diff_flux);
    fields.release_tmp_g(masked);

}





template<typename TF>
void Stats<TF>::calc_stats_path(
        const std::string& varname, const Field3d<TF>& fld)
{
    using namespace Tools_g;
    auto& gd = grid.get_grid_data();

    unsigned int flag;
    const int* nmask;
    std::string name;

    // Calc Integrated Path
    name = varname + "_path";

        cuda_check_error();
    auto masked = fields.get_tmp_g();
    auto mask_proj = fields.get_tmp_xy_g();
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, fld.loc[2]);
            auto threads = grid.get_dim_gpu(gd.icells, gd.jcells, gd.kcells);
            auto threads_ncells = grid.get_dim_gpu(gd.ncells);
            auto threads_ijcells = grid.get_dim_gpu(gd.ncells);

            set_to_val<<<threads_ncells.first, threads_ncells.second>>>(masked->fld_g.data(), gd.ncells, TF(1.));
            apply_mask_g<<<threads.first, threads.second>>>(masked->fld_g.data(),  masked->fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);

            field3d_operators.calc_proj_sum_g(mask_proj->data(), masked->fld_g.data());
            sign_by_arr<TF><<<threads_ijcells.first, threads_ijcells.second>>>(mask_proj->data(), gd.ijcells);
            TF denominator = field3d_operators.calc_sum_2d_g(mask_proj->data());

            apply_mask_g<<<threads.first, threads.second>>>(masked->fld_g.data(),  fld.fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);

            for (int k = gd.kstart; k < gd.kend+1; ++k)
            {
                int kk = k * gd.ijcells;
                mult_by_val<TF><<<threads_ijcells.first, threads_ijcells.second>>>(&masked->fld_g[kk], gd.ijcells, TF(fields.rhoref[k]*gd.dz[k]));
            }

            m.second.tseries.at(name).data = field3d_operators.calc_sum_g(masked->fld_g)/denominator;
            master.sum(&m.second.tseries.at(name).data, 1);
        }
    }

    fields.release_tmp_g(masked);
    fields.release_tmp_xy_g(mask_proj);

}


template<typename TF>
void Stats<TF>::calc_stats_cover(
        const std::string& varname, const Field3d<TF>& fld, const TF offset, const TF threshold)

{
    using namespace Tools_g;
    auto& gd = grid.get_grid_data();

    unsigned int flag;
    const int* nmask;
    std::string name;

    // Calc Integrated Path
    name = varname + "_cover";

        cuda_check_error();
    auto masked = fields.get_tmp_g();
    auto mask_proj = fields.get_tmp_xy_g();
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, fld.loc[2]);
            auto threads = grid.get_dim_gpu(gd.icells, gd.jcells, gd.kcells);
            auto threads_ijcells = grid.get_dim_gpu(gd.ijcells);
            auto threads_ncells = grid.get_dim_gpu(gd.ncells);

            set_to_val<<<threads_ncells.first, threads_ncells.second>>>(masked->fld_g.data(), gd.ncells, TF(1.));
            apply_mask_g<<<threads.first, threads.second>>>(masked->fld_g.data(),  masked->fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);

            field3d_operators.calc_proj_sum_g(mask_proj->data(), masked->fld_g.data());
            sign_by_arr<TF><<<threads_ijcells.first, threads_ijcells.second>>>(mask_proj->data(), gd.ijcells);
            TF denominator = field3d_operators.calc_sum_2d_g(mask_proj->data());

            add_val<<<threads_ncells.first, threads_ncells.second>>>(masked->fld_g.data(), fld.fld_g.data(), gd.ncells, offset - threshold);
            apply_mask_g<<<threads.first, threads.second>>>(masked->fld_g.data(),  masked->fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);

            field3d_operators.calc_proj_sum_g(mask_proj->data(), masked->fld_g.data());
            sign_by_arr<TF><<<threads_ijcells.first, threads_ijcells.second>>>(mask_proj->data(), gd.ijcells);
            TF numerator = field3d_operators.calc_sum_2d_g(mask_proj->data());

            m.second.tseries.at(name).data = numerator / denominator;
            master.sum(&m.second.tseries.at(name).data, 1);
        }
    }

    fields.release_tmp_g(masked);
    fields.release_tmp_xy_g(mask_proj);

}


template<typename TF>
void Stats<TF>::calc_stats_frac(
        const std::string& varname, const Field3d<TF>& fld, const TF offset, const TF threshold)

{
    using namespace Tools_g;

    auto& gd = grid.get_grid_data();

    unsigned int flag;
    const int* nmask;
    std::string name;
    // Calc Fraction
    name = varname + "_frac";

        cuda_check_error();
    auto masked = fields.get_tmp_g();
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        for (auto& m : masks)
        {
            auto threads = grid.get_dim_gpu(gd.icells, gd.jcells, gd.kcells);
            auto threads_ncells = grid.get_dim_gpu(gd.ncells);

            add_val<<<threads_ncells.first, threads_ncells.second>>>(masked->fld_g.data(), fld.fld_g.data(), gd.ncells, offset - threshold);
            sign_by_arr<TF><<<threads_ncells.first, threads_ncells.second>>>(masked->fld_g.data(), gd.ncells);
            set_flag(flag, nmask, m.second, fld.loc[2]);
            apply_mask_g<<<threads.first, threads.second>>>(masked->fld_g.data(),  fld.fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);
            field3d_operators.calc_mean_profile_g(masked->fld_mean_g, masked->fld_g);
            cuda_safe_call(cudaMemcpy(m.second.profs.at(name).data.data(), masked->fld_mean_g.data(), gd.kcells * sizeof(TF), cudaMemcpyDeviceToHost));
            for (int k=gd.kstart; k<gd.kend+1; ++k)
            {
                if (nmask[k])
                    m.second.profs.at(name).data[k] *= gd.itot * gd.jtot / nmask[k];
            }
            master.sum(m.second.profs.at(name).data.data(), gd.kcells);

            set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);
        }
    }
    fields.release_tmp_g(masked);

}


template<typename TF>
void Stats<TF>::calc_tend(Field3d<TF>& fld, const std::string& tend_name)
{

    if (!doing_tendency)
        return;

    using namespace Tools_g;

    auto& gd = grid.get_grid_data();

    unsigned int flag;
    const int* nmask;

    std::string name = fld.name + "_" + tend_name;
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        tendency_order.at(fld.name).push_back(tend_name);
        cuda_check_error();
        auto masked = fields.get_tmp_g();

        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, fld.loc[2]);
            auto threads = grid.get_dim_gpu(gd.icells, gd.jcells, gd.kcells);
            apply_mask_g<<<threads.first, threads.second>>>(masked->fld_g.data(),  fld.fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);
            field3d_operators.calc_mean_profile_g(masked->fld_mean_g, masked->fld_g);
            cuda_safe_call(cudaMemcpy(m.second.profs.at(name).data.data(), masked->fld_mean_g.data(), gd.kcells * sizeof(TF), cudaMemcpyDeviceToHost));
            for (int k=gd.kstart; k<gd.kend+1; ++k)
            {
                if (nmask[k])
                    m.second.profs.at(name).data[k] *= gd.itot * gd.jtot / nmask[k];
            }
            master.sum(m.second.profs.at(name).data.data(), gd.kcells);

            set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);
        }
        fields.release_tmp_g(masked);
    }
}

#endif



template<typename TF>
void Stats<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();

    cuda_safe_call(cudaMalloc(&mfield_g, gd.ncells * sizeof(unsigned int)));
    cuda_safe_call(cudaMalloc(&mfield_bot_g, gd.ijcells * sizeof(unsigned int)));
    cuda_check_error();
}

#ifdef USECUDA
template<typename TF>
void Stats<TF>::initialize_masks()
{
    auto& gd = grid.get_grid_data();
    unsigned int flagmax = 0;
    using namespace Tools_g;

    for (auto& it : masks){
        flagmax += it.second.flag + it.second.flagh;
    }
    auto threads = grid.get_dim_gpu(gd.ncells);
    set_to_val<<<threads.first, threads.second>>>(mfield_g, gd.ncells, flagmax);
    cuda_check_error();


    auto threads_ijcells = grid.get_dim_gpu(gd.ijcells);

    set_to_val<<<threads_ijcells.first, threads_ijcells.second>>>(mfield_bot_g, gd.ijcells, flagmax);
    cuda_check_error();

}


template<typename TF>
void Stats<TF>::finalize_masks()
{
    using namespace Tools_g;
    auto& gd = grid.get_grid_data();

    boundary_cyclic.exec_g(mfield_g);
    boundary_cyclic.exec_2d_g(mfield_bot_g);

    auto masked = fields.get_tmp_g();
    auto ones = fields.get_tmp_g();
    const int nmemsize = gd.kcells*sizeof(int);

    std::vector<TF> nmask_TF(gd.kcells);
    auto threads1d = grid.get_dim_gpu(gd.ncells);

    set_to_val<TF><<<threads1d.first, threads1d.second>>>(
                ones->fld_g, gd.ncells, TF(1.0));
    cuda_check_error();

    for (auto& it : masks)
    {
        unsigned int flag = it.second.flag;
        unsigned int flagh = it.second.flagh;

        // Mask at the full level.

        auto threads = grid.get_dim_gpu(gd.icells, gd.jcells, gd.kcells);
        apply_mask_g<<<threads.first, threads.second>>>(masked->fld_g.data(), ones->fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);

        Grid_kernels::calc_mean_prof_kernel(
                masked->fld_mean_g.data(),
                masked->fld_g.data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                0, gd.kcells,
                gd.icells, gd.ijcells,
                gd.itot*gd.jtot,
                master);

        cuda_safe_call(cudaMemcpy(nmask_TF.data(), masked->fld_mean_g.data(), gd.kcells*sizeof(TF), cudaMemcpyDeviceToHost));

        auto it1 = std::find(varlist.begin(), varlist.end(), "area");
        if (it1 != varlist.end())
        {
            for (int k=0; k<gd.kcells; ++k)
                it.second.profs.at("area").data[k] = nmask_TF[k];
        }

        for (int k=0; k<gd.kcells; ++k)
            it.second.nmask[k] = static_cast<int>(gd.itot*gd.jtot*nmask_TF[k] + 0.5);


        // Mask at the half level.
        apply_mask_g<<<threads.first, threads.second>>>(masked->fld_g.data(), ones->fld_g.data(), mfield_g, flagh, gd.icells, gd.jcells, gd.kcells, gd.ijcells);
        Grid_kernels::calc_mean_prof_kernel(
                masked->fld_mean_g.data(),
                masked->fld_g.data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                0, gd.kcells,
                gd.icells, gd.ijcells,
                gd.itot*gd.jtot,
                master);

        cuda_safe_call(cudaMemcpy(nmask_TF.data(), masked->fld_mean_g.data(), gd.kcells*sizeof(int), cudaMemcpyDeviceToHost));

        it1 = std::find(varlist.begin(), varlist.end(), "areah");
        if (it1 != varlist.end())
        {
            for (int k=0; k<gd.kcells; ++k)
                it.second.profs.at("areah").data[k] = nmask_TF[k];
        }

        for (int k=0; k<gd.kcells; ++k)
            it.second.nmaskh[k] = static_cast<int>(gd.itot*gd.jtot*nmask_TF[k] + 0.5);


        // Set the surface mask.
        it.second.nmask_bot = it.second.nmaskh[gd.kstart];
    }

    fields.release_tmp_g(masked);
    fields.release_tmp_g(ones);
}


template<typename TF>
void Stats<TF>::set_mask_thres(
        std::string mask_name, Field3d<TF>& fld, Field3d<TF>& fldh, TF threshold, Stats_mask_type mode)
{
    auto& gd = grid.get_grid_data();

    unsigned int flag, flagh;
    bool found_mask = false;

    for (auto& it : masks)
    {
        if (it.second.name == mask_name)
        {
            found_mask = true;
            flag = it.second.flag;
            flagh = it.second.flagh;
        }
    }

    if (!found_mask)
        throw std::runtime_error("Invalid mask name in set_mask_thres()");

    auto threads = grid.get_dim_gpu(gd.imax, gd.jmax, gd.kcells);

    if (mode == Stats_mask_type::Plus)
        calc_mask_thres_g<TF, Stats_mask_type::Plus><<<threads.first, threads.second>>>(
                mfield_g, mfield_bot_g, flag, flagh,
                fld.fld_g.data(), fldh.fld_g.data(), fldh.fld_bot_g.data(), threshold,
                gd.istart, gd.jstart, 0,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.jcells, gd.ijcells, gd.kcells);
    else if (mode == Stats_mask_type::Min)
        calc_mask_thres_g<TF, Stats_mask_type::Min><<<threads.first, threads.second>>>(
                mfield_g, mfield_bot_g, flag, flagh,
                fld.fld_g.data(), fldh.fld_g.data(), fldh.fld_bot_g.data(), threshold,
                gd.istart, gd.jstart, 0,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.jcells, gd.ijcells, gd.kcells);
    else
        throw std::runtime_error("Invalid mask type in set_mask_thres()");
}
#endif


#ifdef FLOAT_SINGLE
template class Stats<float>;
#else
template class Stats<double>;
#endif
