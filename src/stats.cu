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

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 gridGPU3(gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    unsigned int flag;
    const int* nmask;
    std::string name;
    auto masked = fields.get_tmp_g();
    // Calc mean of atmospheric variables
    if (std::find(varlist.begin(), varlist.end(), varname) != varlist.end())
    {
        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, fld.loc[2]);
            apply_mask_g<<<gridGPU3, blockGPU>>>(masked->fld_g.data(),  fld.fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);
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

}
#endif


#ifdef USECUDA

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

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 gridGPU3(gridi, gridj, gd.kcells);
    dim3 gridGPU2(gridi, gridj, 1);
    dim3 blockGPU(blocki, blockj, 1);
    const int nblock = gd.ithread_block;
    int ngrid  = gd.ncells/nblock + (gd.ncells%nblock > 0);

    unsigned int flag;
    const int* nmask;
    std::string name;
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

                for (int k = gd.kstart; k < gd.kend+1; ++k)
                {
                    int kk = k * gd.ijcells;
                    add_val<<<gd.ijcells, nblock>>>(&(dev->fld_g[kk]),&(fld.fld_g[kk]), gd.ijcells, - m.second.profs.at(varname).data[k] + offset);
                }
                raise_to_pow<<<gd.ncells, nblock>>>(dev->fld_g.data(), gd.ncells, power);
                apply_mask_g<<<gridGPU3, blockGPU>>>(masked->fld_g.data(),  dev->fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);
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

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);
    const int ngrid  = gd.ncells/blocki + (gd.ncells%blocki > 0);
    const int nblock = gd.ithread_block;
    dim3 gridGPU3(gridi, gridj, gd.kcells);
    dim3 gridGPU2(gridi, gridj, 1);
    dim3 blockGPU(blocki, blockj, 1);

    unsigned int flag;
    const int* nmask;
    std::string name;
    auto masked = fields.get_tmp_g();
    auto advec_flux = fields.get_tmp_g();

    // Calc Resolved Flux
    name = varname + "_w";
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        advec.get_advec_flux(*advec_flux, fld);
        // set_to_val<<<ngrid, nblock>>>(advec_flux->fld_g.data(), gd.ncells, 1.);

        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, fld.loc[2]);
            apply_mask_g<<<gridGPU3, blockGPU>>>(masked->fld_g.data(),  advec_flux->fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);
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

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 gridGPU3(gridi, gridj, gd.kcells);
    dim3 gridGPU2(gridi, gridj, 1);
    dim3 blockGPU(blocki, blockj, 1);

    unsigned int flag;
    const int* nmask;
    std::string name;
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
            apply_mask_g<<<gridGPU3, blockGPU>>>(masked->fld_g.data(), diff_flux->fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);
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

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);
    const int ijgrid  = gd.ijcells/blocki + (gd.ijcells%blocki > 0);
    const int nblock = gd.ithread_block;
    const int ngrid  = gd.ncells/blocki + (gd.ncells%blocki > 0);

    dim3 gridGPU3(gridi, gridj, gd.kcells);
    dim3 gridGPU2(gridi, gridj, 1);
    dim3 blockGPU(blocki, blockj, 1);

    unsigned int flag;
    const int* nmask;
    std::string name;

    // Calc Integrated Path
    name = varname + "_path";

    auto masked = fields.get_tmp_g();
    auto mask_proj = fields.get_tmp_xy_g();
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, fld.loc[2]);
            set_to_val<<<ngrid, nblock>>>(masked->fld_g.data(), gd.ncells, TF(1.));
            apply_mask_g<<<gridGPU3, blockGPU>>>(masked->fld_g.data(),  masked->fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);

            field3d_operators.calc_proj_sum_g(mask_proj->data(), masked->fld_g.data());
            sign_by_arr<TF><<<ijgrid, blocki>>>(mask_proj->data(), gd.ijcells);
            TF denominator = field3d_operators.calc_sum_2d_g(mask_proj->data());

            apply_mask_g<<<gridGPU3, blockGPU>>>(masked->fld_g.data(),  fld.fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);

            for (int k = gd.kstart; k < gd.kend+1; ++k)
            {
                int kk = k * gd.ijcells;
                mult_by_val<TF><<<ijgrid, blocki>>>(&masked->fld_g[kk], gd.ijcells, TF(fields.rhoref[k]*gd.dz[k]));
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

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);
    const int ijgrid  = gd.ijcells/blocki + (gd.ijcells%blocki > 0);
    const int nblock = gd.ithread_block;
    const int ngrid  = gd.ncells/blocki + (gd.ncells%blocki > 0);

    dim3 gridGPU3(gridi, gridj, gd.kcells);
    dim3 gridGPU2(gridi, gridj, 1);
    dim3 blockGPU(blocki, blockj, 1);

    unsigned int flag;
    const int* nmask;
    std::string name;

    // Calc Integrated Path
    name = varname + "_cover";

    auto masked = fields.get_tmp_g();
    auto mask_proj = fields.get_tmp_xy_g();
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, fld.loc[2]);

            set_to_val<<<ngrid, nblock>>>(masked->fld_g.data(), gd.ncells, TF(1.));
            apply_mask_g<<<gridGPU3, blockGPU>>>(masked->fld_g.data(),  masked->fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);

            field3d_operators.calc_proj_sum_g(mask_proj->data(), masked->fld_g.data());
            sign_by_arr<TF><<<ijgrid, blocki>>>(mask_proj->data(), gd.ijcells);
            TF denominator = field3d_operators.calc_sum_2d_g(mask_proj->data());

            add_val<<<gd.ncells, nblock>>>(masked->fld_g.data(), fld.fld_g.data(), gd.ncells, offset - threshold);
            apply_mask_g<<<gridGPU3, blockGPU>>>(masked->fld_g.data(),  masked->fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);

            field3d_operators.calc_proj_sum_g(mask_proj->data(), masked->fld_g.data());
            sign_by_arr<TF><<<ijgrid, blocki>>>(mask_proj->data(), gd.ijcells);
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

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);
    const int nblock = gd.ithread_block;
    const int ngrid  = gd.ncells/blocki + (gd.ncells%blocki > 0);

    dim3 gridGPU3(gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    unsigned int flag;
    const int* nmask;
    std::string name;
    // Calc Fraction
    name = varname + "_frac";

    auto masked = fields.get_tmp_g();
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        for (auto& m : masks)
        {
            add_val<<<ngrid, nblock>>>(masked->fld_g.data(), fld.fld_g.data(), gd.ncells, offset - threshold);
            sign_by_arr<TF><<<ngrid, nblock>>>(masked->fld_g.data(), gd.ncells);
            set_flag(flag, nmask, m.second, fld.loc[2]);
            apply_mask_g<<<gridGPU3, blockGPU>>>(masked->fld_g.data(),  fld.fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);
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



#endif


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
    const int nblock = gd.ithread_block;
    int ngrid  = gd.ncells/nblock + (gd.ncells%nblock > 0);
    set_to_val<<<ngrid, nblock>>>(mfield_g.data(), gd.ncells, flagmax);

    ngrid  = gd.ijcells/nblock + (gd.ijcells%nblock > 0);
    set_to_val<<<ngrid, nblock>>>(mfield_bot_g.data(), gd.ijcells, flagmax);
}


template<typename TF>
void Stats<TF>::finalize_masks()
{
    using namespace Tools_g;
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 gridGPU3(gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    boundary_cyclic.exec_g(mfield_g.data());
    boundary_cyclic.exec_2d_g(mfield_bot_g.data());

    auto masked = fields.get_tmp_g();
    auto ones = fields.get_tmp_g();
    const int nmemsize = gd.kcells*sizeof(int);

    std::vector<TF> nmask_TF(gd.kcells);

    set_to_val_g<TF><<<gridGPU3, blockGPU>>>(
                ones->fld_g, TF(1.0),
                gd.icells, gd.jcells, gd.kcells, gd.ijcells);

    for (auto& it : masks)
    {
        unsigned int flag = it.second.flag;
        unsigned int flagh = it.second.flagh;

        // Mask at the full level.
        apply_mask_g<<<gridGPU3, blockGPU>>>(masked->fld_g.data(), ones->fld_g.data(), mfield_g, flag, gd.icells, gd.jcells, gd.kcells, gd.ijcells);

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
        apply_mask_g<<<gridGPU3, blockGPU>>>(masked->fld_g.data(), ones->fld_g.data(), mfield_g, flagh, gd.icells, gd.jcells, gd.kcells, gd.ijcells);
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
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 gridGPU(gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);


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

    if (mode == Stats_mask_type::Plus)
        calc_mask_thres_g<TF, Stats_mask_type::Plus><<<gridGPU, blockGPU>>>(
                mfield_g.data(), mfield_bot_g.data(), flag, flagh,
                fld.fld_g.data(), fldh.fld_g.data(), fldh.fld_bot_g.data(), threshold,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.jcells, gd.ijcells, gd.kcells);
    else if (mode == Stats_mask_type::Min)
        calc_mask_thres_g<TF, Stats_mask_type::Min><<<gridGPU, blockGPU>>>(
                mfield_g.data(), mfield_bot_g.data(), flag, flagh,
                fld.fld_g.data(), fldh.fld_g.data(), fldh.fld_bot_g.data(), threshold,
                gd.istart, gd.jstart, gd.kstart,
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
