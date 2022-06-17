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
#include "grid.h"
#include "fields.h"
#include "thermo_dry.h"
#include "defines.h"
#include "constants.h"
#include "master.h"
#include "column.h"
#include "tools.h"
#include "stats.h"
#include "finite_difference.h"

namespace
{
    using namespace Constants;
    using namespace Finite_difference::O2;

    template<typename TF> __global__
    void calc_buoyancy_tend_2nd_g(TF* __restrict__ wt,
                                  TF* __restrict__ th, TF* __restrict__ threfh,
                                  int istart, int jstart, int kstart,
                                  int iend,   int jend,   int kend,
                                  int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            wt[ijk] += grav<TF>/threfh[k] * (static_cast<TF>(0.5)*(th[ijk-kk]+th[ijk]) - threfh[k]);
        }
    }

    template<typename TF> __global__
    void calc_baroclinic_2nd_g(TF* __restrict__ tht, const TF* __restrict__ v,
                               const TF dthetady_ls,
                               int istart, int jstart, int kstart,
                               int iend,   int jend,   int kend,
                               int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            tht[ijk] -= dthetady_ls * interp2(v[ijk], v[ijk+jj]);
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_g(TF* __restrict__ b,
                         TF* __restrict__ th, TF* __restrict__ thref,
                         int istart, int jstart,
                         int iend,   int jend,   int kcells,
                         int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z;

        if (i < iend && j < jend && k < kcells)
        {
            const int ijk = i + j*jj + k*kk;
            b[ijk] = grav<TF>/thref[k] * (th[ijk] - thref[k]);
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_bot_g(TF* __restrict__ b,     TF* __restrict__ bbot,
                             TF* __restrict__ th,    TF* __restrict__ thbot,
                             TF* __restrict__ thref, TF* __restrict__ threfh,
                             int kstart, int icells, int jcells,
                             int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;

            bbot[ij] = grav<TF>/threfh[kstart] * (thbot[ij] - threfh[kstart]);
            b[ijk]   = grav<TF>/thref [kstart] * (th[ijk]   - thref [kstart]);
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_flux_bot_g(TF* __restrict__ bfluxbot, TF* __restrict__ thfluxbot,
                                  TF* __restrict__ threfh,
                                  int kstart, int icells, int jcells,
                                  int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;
            bfluxbot[ij] = grav<TF>/threfh[kstart]*thfluxbot[ij];
        }
    }

    template<typename TF> __global__
    void calc_N2_g(TF* __restrict__ N2,    TF* __restrict__ th,
                   TF* __restrict__ thref, TF* __restrict__ dzi,
                   int istart, int jstart, int kstart,
                   int iend,   int jend,   int kend,
                   int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            N2[ijk] = grav<TF>/thref[k]*static_cast<TF>(0.5)*(th[ijk+kk] - th[ijk-kk])*dzi[k];
        }
    }
} // end namespace

template<typename TF>
void Thermo_dry<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();

    const int nmemsize = gd.kcells*sizeof(TF);

    // Allocate fields for Boussinesq and anelastic solver
    cuda_safe_call(cudaMalloc(&bs.thref_g,   nmemsize));
    cuda_safe_call(cudaMalloc(&bs.threfh_g,  nmemsize));
    cuda_safe_call(cudaMalloc(&bs.pref_g,    nmemsize));
    cuda_safe_call(cudaMalloc(&bs.prefh_g,   nmemsize));
    cuda_safe_call(cudaMalloc(&bs.exnref_g,  nmemsize));
    cuda_safe_call(cudaMalloc(&bs.exnrefh_g, nmemsize));


    // Copy fields to device
    cuda_safe_call(cudaMemcpy(bs.thref_g,   bs.thref.data(),   nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.threfh_g,  bs.threfh.data(),  nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.pref_g,    bs.pref.data(),    nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.prefh_g,   bs.prefh.data(),   nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.exnref_g,  bs.exnref.data(),  nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.exnrefh_g, bs.exnrefh.data(), nmemsize, cudaMemcpyHostToDevice));

}

template<typename TF>
void Thermo_dry<TF>::clear_device()
{
    cuda_safe_call(cudaFree(bs.thref_g ));
    cuda_safe_call(cudaFree(bs.threfh_g));
    cuda_safe_call(cudaFree(bs.pref_g  ));
    cuda_safe_call(cudaFree(bs.prefh_g ));
    cuda_safe_call(cudaFree(bs.exnref_g ));
    cuda_safe_call(cudaFree(bs.exnrefh_g));
    tdep_pbot->clear_device();
}

template<typename TF>
void Thermo_dry<TF>::forward_device()
{
    auto& gd = grid.get_grid_data();

    const int nmemsize = gd.kcells*sizeof(TF);
    // Copy fields to device
    cuda_safe_call(cudaMemcpy(bs.pref_g,    bs.pref.data(),    nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.prefh_g,   bs.prefh.data(),   nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.exnref_g,  bs.exnref.data(),  nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.exnrefh_g, bs.exnrefh.data(), nmemsize, cudaMemcpyHostToDevice));
}
template<typename TF>
void Thermo_dry<TF>::backward_device()
{
    auto& gd = grid.get_grid_data();

    const int nmemsize = gd.kcells*sizeof(TF);
    cudaMemcpy(bs.pref_g,    bs.pref.data(),    nmemsize, cudaMemcpyHostToDevice);
    cudaMemcpy(bs.prefh_g,   bs.prefh.data(),   nmemsize, cudaMemcpyHostToDevice);
    cudaMemcpy(bs.exnref_g,  bs.exnref.data(),  nmemsize, cudaMemcpyHostToDevice);
    cudaMemcpy(bs.exnrefh_g, bs.exnrefh.data(), nmemsize, cudaMemcpyHostToDevice);

    bs_stats = bs;
}

#ifdef USECUDA
template<typename TF>
void Thermo_dry<TF>::exec(const double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    if (grid.get_spatial_order() == Grid_order::Second)
    {
        calc_buoyancy_tend_2nd_g<TF><<<gridGPU, blockGPU>>>(
            fields.mt.at("w")->fld_g, fields.sp.at("th")->fld_g, bs.threfh_g,
            gd.istart, gd.jstart, gd.kstart+1,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.ijcells);
        cuda_check_error();

        if (swbaroclinic)
        {
            calc_baroclinic_2nd_g<TF><<<gridGPU, blockGPU>>>(
                fields.st.at("th")->fld_g, fields.mp.at("v")->fld_g,
                dthetady_ls,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();
        }
    }
    else if (grid.get_spatial_order() == Grid_order::Fourth)
    {
        std::string msg = "4th order thermo_dry not (yet) implemented";
        throw std::runtime_error(msg);
    }
    cudaDeviceSynchronize();
    stats.calc_tend(*fields.mt.at("w"), tend_name);

}
#endif

#ifdef USECUDA
template<typename TF>
void Thermo_dry<TF>::get_thermo_field_g(
        Field3d<TF>& fld, const std::string& name, const bool cyclic)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    dim3 gridGPU2 (gridi, gridj, gd.kmax);
    dim3 blockGPU2(blocki, blockj, 1);

    if (name == "b")
    {
        calc_buoyancy_g<TF><<<gridGPU, blockGPU>>>(
            fld.fld_g, fields.sp.at("th")->fld_g, bs.thref_g,
            gd.istart, gd.jstart,
            gd.iend, gd.jend, gd.kcells,
            gd.icells, gd.ijcells);
        cuda_check_error();
    }
    else if (name == "N2")
    {
        calc_N2_g<TF><<<gridGPU2, blockGPU2>>>(
            fld.fld_g, fields.sp.at("th")->fld_g, bs.thref_g, gd.dzi_g,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.ijcells);
        cuda_check_error();
    }
    else
    {
        std::string msg = "get_thermo_field \"" + name + "\" not supported";
        throw std::runtime_error(msg);
    }

    if (cyclic)
        boundary_cyclic.exec_g(fld.fld_g);
}
#endif

#ifdef USECUDA
template<typename TF>
void Thermo_dry<TF>::get_buoyancy_fluxbot_g(Field3d<TF>& b)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 gridGPU (gridi, gridj, 1);
    dim3 blockGPU(blocki, blockj, 1);

    calc_buoyancy_flux_bot_g<TF><<<gridGPU, blockGPU>>>(
        b.flux_bot_g, fields.sp.at("th")->flux_bot_g,
        bs.threfh_g, gd.kstart, gd.icells, gd.jcells,
        gd.icells, gd.ijcells);
    cuda_check_error();
}
#endif

#ifdef USECUDA
template<typename TF>
void Thermo_dry<TF>::get_buoyancy_surf_g(Field3d<TF>& b)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 gridGPU (gridi, gridj, 1);
    dim3 blockGPU(blocki, blockj, 1);

    calc_buoyancy_bot_g<TF><<<gridGPU, blockGPU>>>(
        b.fld_g, b.flux_bot_g,
        fields.sp.at("th")->fld_g, fields.sp.at("th")->fld_bot_g,
        bs.thref_g, bs.threfh_g, gd.kstart, gd.icells, gd.jcells,
        gd.icells, gd.ijcells);
    cuda_check_error();

    calc_buoyancy_flux_bot_g<TF><<<gridGPU, blockGPU>>>(
        b.flux_bot_g, fields.sp.at("th")->flux_bot_g,
        bs.threfh_g, gd.kstart, gd.icells, gd.jcells,
        gd.icells, gd.ijcells);
    cuda_check_error();
}
#endif

#ifdef USECUDA
template<typename TF>
void Thermo_dry<TF>::exec_column(Column<TF>& column)
{
    const TF no_offset = 0.;
    auto output = fields.get_tmp_g();

    get_thermo_field_g(*output, "b", false);
    column.calc_column("b", output->fld_g, no_offset);

    fields.release_tmp_g(output);
}
#endif
template class Thermo_dry<double>;
template class Thermo_dry<float>;
