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
#include "thermo_vapor.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "master.h"
#include "tools.h"
#include "column.h"
#include "stats.h"

#include "thermo_moist_functions.h"

namespace
{
    using namespace Constants;
    using namespace Finite_difference::O2;
    using namespace Thermo_moist_functions;

    template<typename TF> __global__
    void calc_buoyancy_tend_2nd_g(TF* __restrict__ wt, TF* __restrict__ th, TF* __restrict__ qt,
                                  TF* __restrict__ thvrefh, TF* __restrict__ exnh, TF* __restrict__ ph,
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

            // Half level temperature and vaporure content
            const TF thh = TF(0.5) * (th[ijk-kk] + th[ijk]);         // Half level liq. water pot. temp.
            const TF qth = TF(0.5) * (qt[ijk-kk] + qt[ijk]);         // Half level specific hum.
            wt[ijk] += buoyancy_no_ql(thh, qth, thvrefh[k]);
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_g(TF* __restrict__ b,  TF* __restrict__ th,
                         TF* __restrict__ qt, TF* __restrict__ thvref,
                         TF* __restrict__ p,  TF* __restrict__ exn,
                         int istart, int jstart, int kstart,
                         int iend,   int jend,   int kcells,
                         int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z;

        if (i < iend && j < jend && k < kcells)
        {
            const int ijk   = i + j*jj + k*kk;
            b[ijk] = buoyancy_no_ql(th[ijk], qt[ijk], thvref[k]);
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_bot_g(TF* __restrict__ b,      TF* __restrict__ bbot,
                             TF* __restrict__ th,     TF* __restrict__ thbot,
                             TF* __restrict__ qt,     TF* __restrict__ qtbot,
                             TF* __restrict__ thvref, TF* __restrict__ thvrefh,
                             int kstart, int icells, int jcells,
                             int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;

            bbot[ij ] = buoyancy_no_ql(thbot[ij], qtbot[ij], thvrefh[kstart]);
            b   [ijk] = buoyancy_no_ql(th[ijk],   qt[ijk],   thvref[kstart]);
        }
    }

    template<typename TF> __global__
    void calc_buoyancy_flux_bot_g(TF* __restrict__ bfluxbot,
                                  TF* __restrict__ th, TF* __restrict__ thfluxbot,
                                  TF* __restrict__ qt, TF* __restrict__ qtfluxbot,
                                  TF* __restrict__ thvrefh,
                                  int kstart, int icells, int jcells,
                                  int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            bfluxbot[ij] = buoyancy_flux_no_ql(th[ijk], thfluxbot[ij], qt[ijk], qtfluxbot[ij], thvrefh[kstart]);
        }
    }

    template<typename TF> __global__
    void calc_N2_g(TF* __restrict__ N2, TF* __restrict__ th,
                   TF* __restrict__ thvref, TF* __restrict__ dzi,
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
            N2[ijk] = grav<TF>/thvref[k]*TF(0.5)*(th[ijk+kk] - th[ijk-kk])*dzi[k];
        }
    }

} // end name    space

template<typename TF>
void Thermo_vapor<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();
    const int nmemsize = gd.kcells*sizeof(TF);

    // Allocate fields for Boussinesq and anelastic solver
    cuda_safe_call(cudaMalloc(&bs.thvref_g,  nmemsize));
    cuda_safe_call(cudaMalloc(&bs.thvrefh_g, nmemsize));
    cuda_safe_call(cudaMalloc(&bs.pref_g,    nmemsize));
    cuda_safe_call(cudaMalloc(&bs.prefh_g,   nmemsize));
    cuda_safe_call(cudaMalloc(&bs.exnref_g,  nmemsize));
    cuda_safe_call(cudaMalloc(&bs.exnrefh_g, nmemsize));

    // Copy fields to device
    cuda_safe_call(cudaMemcpy(bs.thvref_g,  bs.thvref.data(),  nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.thvrefh_g, bs.thvrefh.data(), nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.pref_g,    bs.pref.data(),    nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.prefh_g,   bs.prefh.data(),   nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.exnref_g,  bs.exnref.data(),  nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.exnrefh_g, bs.exnrefh.data(), nmemsize, cudaMemcpyHostToDevice));
}

template<typename TF>
void Thermo_vapor<TF>::clear_device()
{
    cuda_safe_call(cudaFree(bs.thvref_g ));
    cuda_safe_call(cudaFree(bs.thvrefh_g));
    cuda_safe_call(cudaFree(bs.pref_g   ));
    cuda_safe_call(cudaFree(bs.prefh_g  ));
    cuda_safe_call(cudaFree(bs.exnref_g ));
    cuda_safe_call(cudaFree(bs.exnrefh_g));
    tdep_pbot->clear_device();
}

template<typename TF>
void Thermo_vapor<TF>::forward_device()
{
    // Copy fields to device
    auto& gd = grid.get_grid_data();
    const int nmemsize = gd.kcells*sizeof(TF);
    cuda_safe_call(cudaMemcpy(bs.pref_g,    bs.pref.data(),    nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.prefh_g,   bs.prefh.data(),   nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.exnref_g,  bs.exnref.data(),  nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bs.exnrefh_g, bs.exnrefh.data(), nmemsize, cudaMemcpyHostToDevice));
}

template<typename TF>
void Thermo_vapor<TF>::backward_device()
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
void Thermo_vapor<TF>::exec(const double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);


    // Re-calculate hydrostatic pressure and exner
    if (bs.swupdatebasestate)
    {
        //calc_hydrostatic_pressure<TF><<<1, 1>>>(bs.pref_g, bs.prefh_g, bs.exnref_g, bs.exnrefh_g,
        //                                        fields.sp.at("thl")->fld_mean_g, fields.sp.at("qt")->fld_mean_g,
        //                                        gd.z_g, gd.dz_g, gd.dzh_g, bs.pbot, gd.kstart, gd.kend);
        //cuda_check_error();

        // BvS: Calculating hydrostatic pressure on GPU is extremely slow. As temporary solution, copy back mean profiles to host,
        //      calculate pressure there and copy back the required profiles.
        cudaMemcpy(fields.sp.at("thl")->fld_mean.data(), fields.sp.at("thl")->fld_mean_g, gd.kcells*sizeof(TF), cudaMemcpyDeviceToHost);
        cudaMemcpy(fields.sp.at("qt")->fld_mean.data(),  fields.sp.at("qt")->fld_mean_g,  gd.kcells*sizeof(TF), cudaMemcpyDeviceToHost);

        auto tmp = fields.get_tmp();

        calc_base_state_no_ql(bs.pref.data(), bs.prefh.data(),
                        &tmp->fld[0*gd.kcells], &tmp->fld[1*gd.kcells], &tmp->fld[2*gd.kcells], &tmp->fld[3*gd.kcells],
                        bs.exnref.data(), bs.exnrefh.data(), fields.sp.at("thl")->fld_mean.data(), fields.sp.at("qt")->fld_mean.data(),
                        bs.pbot, gd.kstart, gd.kend, gd.z.data(), gd.dz.data(), gd.dzh.data());

        fields.release_tmp(tmp);

        // Only half level pressure and bs.exner needed for BuoyancyTend()
        cudaMemcpy(bs.prefh_g,   bs.prefh.data(),   gd.kcells*sizeof(TF), cudaMemcpyHostToDevice);
        cudaMemcpy(bs.exnrefh_g, bs.exnrefh.data(), gd.kcells*sizeof(TF), cudaMemcpyHostToDevice);
    }

    calc_buoyancy_tend_2nd_g<TF><<<gridGPU, blockGPU>>>(
        fields.mt.at("w")->fld_g, fields.sp.at("thl")->fld_g,
        fields.sp.at("qt")->fld_g, bs.thvrefh_g, bs.exnrefh_g, bs.prefh_g,
        gd.istart,  gd.jstart, gd.kstart+1,
        gd.iend,    gd.jend,   gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    cudaDeviceSynchronize();
    stats.calc_tend(*fields.mt.at("w"), tend_name);

}
#endif

#ifdef USECUDA
template<typename TF>
void Thermo_vapor<TF>::get_thermo_field_g(
        Field3d<TF>& fld, const std::string& name, const bool cyclic )
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

    // Re-calculate hydrostatic pressure and exner
    if (bs.swupdatebasestate)
    {
        //calc_hydrostatic_pressure<TF><<<1, 1>>>(bs.pref_g, bs.prefh_g, bs.exnref_g, bs.exnrefh_g,
        //                                        fields.sp.at("thl")->fld_mean_g, fields.sp.at("qt")->fld_mean_g,
        //                                        gd.z_g, gd.dz_g, gd.dzh_g, bs.pbot, gd.kstart, gd.kend);
        //cuda_check_error();

        // BvS: Calculating hydrostatic pressure on GPU is extremely slow. As temporary solution, copy back mean profiles to host,
        //      calculate pressure there and copy back the required profiles.
        cudaMemcpy(fields.sp.at("thl")->fld_mean.data(), fields.sp.at("thl")->fld_mean_g, gd.kcells*sizeof(TF), cudaMemcpyDeviceToHost);
        cudaMemcpy(fields.sp.at("qt")->fld_mean.data(),  fields.sp.at("qt")->fld_mean_g,  gd.kcells*sizeof(TF), cudaMemcpyDeviceToHost);

        auto tmp = fields.get_tmp();

        calc_base_state_no_ql(bs.pref.data(), bs.prefh.data(),
                        &tmp->fld[0*gd.kcells], &tmp->fld[1*gd.kcells], &tmp->fld[2*gd.kcells], &tmp->fld[3*gd.kcells],
                        bs.exnref.data(), bs.exnrefh.data(), fields.sp.at("thl")->fld_mean.data(), fields.sp.at("qt")->fld_mean.data(),
                        bs.pbot, gd.kstart, gd.kend, gd.z.data(), gd.dz.data(), gd.dzh.data());

        fields.release_tmp(tmp);

        // Only half level pressure and bs.exner needed for BuoyancyTend()
        cudaMemcpy(bs.prefh_g,   bs.prefh.data(),   gd.kcells*sizeof(TF), cudaMemcpyHostToDevice);
        cudaMemcpy(bs.exnrefh_g, bs.exnrefh.data(), gd.kcells*sizeof(TF), cudaMemcpyHostToDevice);
    }


    if (name == "b")
    {
        calc_buoyancy_g<TF><<<gridGPU, blockGPU>>>(
            fld.fld_g, fields.sp.at("thl")->fld_g, fields.sp.at("qt")->fld_g,
            bs.thvref_g, bs.pref_g, bs.exnref_g,
            gd.istart,  gd.jstart, gd.kstart,
            gd.iend, gd.jend, gd.kcells,
            gd.icells, gd.ijcells);
        cuda_check_error();
    }
    else if (name == "N2")
    {
        calc_N2_g<TF><<<gridGPU2, blockGPU2>>>(
            fld.fld_g, fields.sp.at("thl")->fld_g, bs.thvref_g, gd.dzi_g,
            gd.istart,  gd.jstart, gd.kstart,
            gd.iend,    gd.jend,   gd.kend,
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
TF* Thermo_vapor<TF>::get_basestate_fld_g(std::string name)
{
    // BvS TO-DO: change std::string to enum
    if (name == "pref")
        return bs.pref_g;
    else if (name == "prefh")
        return bs.prefh_g;
    else if (name == "exner")
        return bs.exnref_g;
    else if (name == "exnerh")
        return bs.exnrefh_g;
    else
    {
        std::string error_message = "Can not get basestate field \"" + name + "\" from thermo_moist";
        throw std::runtime_error(error_message);
    }
}
#endif

#ifdef USECUDA
template<typename TF>
void Thermo_vapor<TF>::get_buoyancy_fluxbot_g(Field3d<TF>& bfield)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 gridGPU (gridi, gridj, 1);
    dim3 blockGPU(blocki, blockj, 1);

    calc_buoyancy_flux_bot_g<TF><<<gridGPU, blockGPU>>>(
        bfield.flux_bot_g,
        fields.sp.at("thl")->fld_g, fields.sp.at("thl")->flux_bot_g,
        fields.sp.at("qt")->fld_g, fields.sp.at("qt")->flux_bot_g,
        bs.thvrefh_g, gd.kstart, gd.icells, gd.jcells,
        gd.icells, gd.ijcells);
    cuda_check_error();
}
#endif

#ifdef USECUDA
template<typename TF>
void Thermo_vapor<TF>::get_buoyancy_surf_g(Field3d<TF>& bfield)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 gridGPU (gridi, gridj, 1);
    dim3 blockGPU(blocki, blockj, 1);

    calc_buoyancy_bot_g<TF><<<gridGPU, blockGPU>>>(
        bfield.fld_g, bfield.fld_bot_g,
        fields.sp.at("thl")->fld_g, fields.sp.at("thl")->fld_bot_g,
        fields.sp.at("qt")->fld_g, fields.sp.at("qt")->fld_bot_g,
        bs.thvref_g, bs.thvrefh_g, gd.kstart, gd.icells, gd.jcells,
        gd.icells, gd.ijcells);
    cuda_check_error();

    calc_buoyancy_flux_bot_g<TF><<<gridGPU, blockGPU>>>(
        bfield.flux_bot_g,
        fields.sp.at("thl")->fld_g, fields.sp.at("thl")->flux_bot_g,
        fields.sp.at("qt")->fld_g, fields.sp.at("qt")->flux_bot_g,
        bs.thvrefh_g, gd.kstart, gd.icells, gd.jcells,
        gd.icells, gd.ijcells);
    cuda_check_error();
}
#endif

#ifdef USECUDA
template<typename TF>
void Thermo_vapor<TF>::exec_column(Column<TF>& column)
{
    const TF no_offset = 0.;
    auto output = fields.get_tmp_g();

    get_thermo_field_g(*output, "b", false);
    column.calc_column("b", output->fld_g, no_offset);

    fields.release_tmp_g(output);
}
#endif
template class Thermo_vapor<double>;
template class Thermo_vapor<float>;
