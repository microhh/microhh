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
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "constants.h"
#include "thermo.h"
#include "model.h"
#include "master.h"
#include "tools.h"
#include "timedep.h"
#include "monin_obukhov.h"
#include "boundary_surface.h"
#include "boundary_surface_kernels_gpu.h"
#include "fast_math.h"
#include "column.h"

namespace
{
    namespace most = Monin_obukhov;
    namespace fm = Fast_math;
    namespace bsk = Boundary_surface_kernels_g;

    const int nzL = 10000; // Size of the lookup table for MO iterations.

    template<typename TF> __device__
    TF find_obuk_g(
            const float* const __restrict__ zL,
            const float* const __restrict__ f,
            int &n,
            const TF Ri,
            const TF zsl)
    {
        // Determine search direction.
        if ((f[n]-Ri) > 0.f)
            while ( (f[n-1]-Ri) > 0.f && n > 0) { --n; }
        else
            while ( (f[n]-Ri) < 0.f && n < (nzL-1) ) { ++n; }

        const TF zL0 = (n == 0 || n == nzL-1) ? zL[n] : zL[n-1] + (Ri-f[n-1]) / (f[n]-f[n-1]) * (zL[n]-zL[n-1]);

        return zsl/zL0;
    }


    template<typename TF> __device__
    TF calc_obuk_noslip_flux_g(
            const float* const __restrict__ zL,
            const float* const __restrict__ f,
            int& n,
            const TF du,
            const TF bfluxbot,
            const TF zsl)
    {
        // Calculate the appropriate Richardson number.
        const TF Ri = -Constants::kappa<TF> * bfluxbot * zsl / fm::pow3(du);
        return find_obuk_g(zL, f, n, Ri, zsl);
    }

    template<typename TF> __device__
    TF calc_obuk_noslip_dirichlet_g(
            const float* const __restrict__ zL,
            const float* const __restrict__ f,
            int& n,
            const TF du,
            const TF db,
            const TF zsl)
    {
        // Calculate the appropriate Richardson number.
        const TF Ri = Constants::kappa<TF> * db * zsl / fm::pow2(du);
        return find_obuk_g(zL, f, n, Ri, zsl);
    }

    template<typename TF> __global__
    void stability_g(
            TF* const __restrict__ ustar,
            TF* const __restrict__ obuk,
            int* const __restrict__ nobuk_g,
            const TF* const __restrict__ b,
            const TF* const __restrict__ bbot,
            const TF* const __restrict__ bfluxbot,
            const TF* const __restrict__ dutot,
            const TF* const __restrict__ z0m,
            const float* const __restrict__ zL_sl_g,
            const float* const __restrict__ f_sl_g,
            const TF db_ref,
            const TF zsl,
            const int icells, const int jcells,
            const int kstart, const int jj, int kk,
            const Boundary_type mbcbot,
            const Boundary_type thermobc)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;

            // case 1: fixed buoyancy flux and fixed ustar
            if (mbcbot == Boundary_type::Ustar_type && thermobc == Boundary_type::Flux_type)
            {
                obuk[ij] = -fm::pow3(ustar[ij]) / (Constants::kappa<TF>*bfluxbot[ij]);
            }
            // case 2: fixed buoyancy flux and free ustar
            else if (mbcbot == Boundary_type::Dirichlet_type && thermobc == Boundary_type::Flux_type)
            {
                obuk [ij] = calc_obuk_noslip_flux_g(zL_sl_g, f_sl_g, nobuk_g[ij], dutot[ij], bfluxbot[ij], zsl);
                ustar[ij] = dutot[ij] * most::fm(zsl, z0m[ij], obuk[ij]);
            }
            // case 3: fixed buoyancy surface value and free ustar
            else if (mbcbot == Boundary_type::Dirichlet_type && thermobc == Boundary_type::Dirichlet_type)
            {
                TF db = b[ijk] - bbot[ij] + db_ref;
                obuk [ij] = calc_obuk_noslip_dirichlet_g(zL_sl_g, f_sl_g, nobuk_g[ij], dutot[ij], db, zsl);
                ustar[ij] = dutot[ij] * most::fm(zsl, z0m[ij], obuk[ij]);
            }
        }
    }

    template<typename TF> __global__
    void stability_neutral_g(
            TF* const __restrict__ ustar,
            TF* const __restrict__ obuk,
            const TF* const __restrict__ dutot,
            const TF* const __restrict__ z0m,
            const TF zsl,
            const int icells, const int jcells, const int jj,
            const Boundary_type mbcbot,
            const Boundary_type thermobc)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;

            // case 1: fixed buoyancy flux and fixed ustar
            if (mbcbot == Boundary_type::Ustar_type && thermobc == Boundary_type::Flux_type)
            {
                obuk[ij] = -Constants::dbig;
            }
            // case 2: fixed buoyancy flux and free ustar
            else if (mbcbot == Boundary_type::Dirichlet_type && thermobc == Boundary_type::Flux_type)
            {
                obuk [ij] = -Constants::dbig;
                ustar[ij] = dutot[ij] * most::fm(zsl, z0m[ij], obuk[ij]);
            }
            // case 3: fixed buoyancy surface value and free ustar
            else if (mbcbot == Boundary_type::Dirichlet_type && thermobc == Boundary_type::Dirichlet_type)
            {
                obuk [ij] = -Constants::dbig;
                ustar[ij] = dutot[ij] * most::fm(zsl, z0m[ij], obuk[ij]);
            }
        }
    }

    template<typename TF> __global__
    void surfm_flux_g(
            TF* const __restrict__ ufluxbot,
            TF* const __restrict__ vfluxbot,
            const TF* const __restrict__ u,
            const TF* const __restrict__ v,
            const TF* const __restrict__ ubot,
            const TF* const __restrict__ vbot,
            const TF* const __restrict__ ustar,
            const TF* const __restrict__ obuk,
            const TF* const __restrict__ z0m,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int jj, const int kk,
            const Boundary_type bcbot)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ii  = 1;
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;

            if (bcbot == Boundary_type::Dirichlet_type)
            {
                // interpolate the whole stability function rather than ustar or obuk
                ufluxbot[ij] = -(u[ijk]-ubot[ij])*TF(0.5)*(ustar[ij-ii]*most::fm(zsl, z0m[ij], obuk[ij-ii]) +
                                                           ustar[ij]   *most::fm(zsl, z0m[ij], obuk[ij]));
                vfluxbot[ij] = -(v[ijk]-vbot[ij])*TF(0.5)*(ustar[ij-jj]*most::fm(zsl, z0m[ij], obuk[ij-jj]) +
                                                           ustar[ij]   *most::fm(zsl, z0m[ij], obuk[ij]));
            }
            else if (bcbot == Boundary_type::Ustar_type)
            {
                const TF minval = 1.e-2;

                // minimize the wind at 0.01, thus the wind speed squared at 0.0001
                const TF vonu2 = fmax(minval, TF(0.25)*( fm::pow2(v[ijk-ii]-vbot[ij-ii]) + fm::pow2(v[ijk-ii+jj]-vbot[ij-ii+jj])
                                                       + fm::pow2(v[ijk   ]-vbot[ij   ]) + fm::pow2(v[ijk   +jj]-vbot[ij   +jj])) );
                const TF uonv2 = fmax(minval, TF(0.25)*( fm::pow2(u[ijk-jj]-ubot[ij-jj]) + fm::pow2(u[ijk+ii-jj]-ubot[ij+ii-jj])
                                                       + fm::pow2(u[ijk   ]-ubot[ij   ]) + fm::pow2(u[ijk+ii   ]-ubot[ij+ii   ])) );

                const TF u2 = fmax(minval, fm::pow2(u[ijk]-ubot[ij]));
                const TF v2 = fmax(minval, fm::pow2(v[ijk]-vbot[ij]));

                const TF ustaronu4 = TF(0.5)*(fm::pow4(ustar[ij-ii]) + fm::pow4(ustar[ij]));
                const TF ustaronv4 = TF(0.5)*(fm::pow4(ustar[ij-jj]) + fm::pow4(ustar[ij]));

                ufluxbot[ij] = -copysign(TF(1.), u[ijk]-ubot[ij]) * pow(ustaronu4 / (TF(1.) + vonu2 / u2), TF(0.5));
                vfluxbot[ij] = -copysign(TF(1.), v[ijk]-vbot[ij]) * pow(ustaronv4 / (TF(1.) + uonv2 / v2), TF(0.5));
            }
        }
    }

    template<typename TF> __global__
    void surfm_grad_g(
            TF* const __restrict__ ugradbot,
            TF* const __restrict__ vgradbot,
            const TF* const __restrict__ u,
            const TF* const __restrict__ v,
            const TF* const __restrict__ ubot,
            const TF* const __restrict__ vbot,
            const TF zsl,
            const int icells, const int jcells,
            const int kstart, const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;

            ugradbot[ij] = (u[ijk]-ubot[ij])/zsl;
            vgradbot[ij] = (v[ijk]-vbot[ij])/zsl;
        }
    }

    template<typename TF> __global__
    void surfs_g(
            TF* const __restrict__ varfluxbot,
            TF* const __restrict__ vargradbot,
            TF* const __restrict__ varbot,
            const TF* const __restrict__ var,
            const TF* const __restrict__ ustar,
            const TF* const __restrict__ obuk,
            const TF* const __restrict__ z0h,
            const TF zsl,
            const int icells, const int jcells,
            const int kstart, const int jj, const int kk,
            const Boundary_type bcbot)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;

            if (bcbot == Boundary_type::Dirichlet_type)
            {
                varfluxbot[ij] = -(var[ijk]-varbot[ij])*ustar[ij]*most::fh(zsl, z0h[ij], obuk[ij]);
                vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
            }
            else if (bcbot == Boundary_type::Flux_type)
            {
                varbot[ij]     = varfluxbot[ij] / (ustar[ij]*most::fh(zsl, z0h[ij], obuk[ij])) + var[ijk];
                vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
            }
        }
    }
}

template<typename TF>
void Boundary_surface<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();

    // Prepare base boundary, for inflow profiles.
    Boundary<TF>::prepare_device();

    const int dmemsize2d = gd.ijcells*sizeof(TF);
    const int imemsize2d = gd.ijcells*sizeof(int);
    const int dimemsize  = gd.icells*sizeof(TF);
    const int iimemsize  = gd.icells*sizeof(int);

    // 2D fields:
    cuda_safe_call(cudaMalloc(&obuk_g,    dmemsize2d));
    cuda_safe_call(cudaMalloc(&ustar_g,   dmemsize2d));
    cuda_safe_call(cudaMalloc(&z0m_g,     dmemsize2d));
    cuda_safe_call(cudaMalloc(&z0h_g,     dmemsize2d));
    cuda_safe_call(cudaMalloc(&dudz_mo_g, dmemsize2d));
    cuda_safe_call(cudaMalloc(&dvdz_mo_g, dmemsize2d));
    cuda_safe_call(cudaMalloc(&dbdz_mo_g, dmemsize2d));

    // Lookuk table:
    if (sw_constant_z0)
    {
        cuda_safe_call(cudaMalloc(&nobuk_g, imemsize2d));
        cuda_safe_call(cudaMalloc(&zL_sl_g, nzL*sizeof(float)));
        cuda_safe_call(cudaMalloc(&f_sl_g,  nzL*sizeof(float)));
    }

    // Copy data to GPU:
    forward_device();
}

template<typename TF>
void Boundary_surface<TF>::forward_device()
{
    auto& gd = grid.get_grid_data();

    const int dmemsize2d = gd.ijcells * sizeof(TF);
    const int imemsize2d = gd.ijcells * sizeof(int);
    const int flutmemsize = nzL * sizeof(float);

    // 2D fields:
    cuda_safe_call(cudaMemcpy(obuk_g,    obuk.data(),    dmemsize2d, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(ustar_g,   ustar.data(),   dmemsize2d, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(z0m_g,     z0m.data(),     dmemsize2d, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(z0h_g,     z0h.data(),     dmemsize2d, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(dudz_mo_g, dudz_mo.data(), dmemsize2d, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(dvdz_mo_g, dvdz_mo.data(), dmemsize2d, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(dbdz_mo_g, dbdz_mo.data(), dmemsize2d, cudaMemcpyHostToDevice));

    if (sw_constant_z0)
    {
        // Lookup table:
        cuda_safe_call(cudaMemcpy(nobuk_g, nobuk.data(), imemsize2d,  cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(zL_sl_g, zL_sl.data(), flutmemsize, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(f_sl_g,  f_sl.data(),  flutmemsize, cudaMemcpyHostToDevice));
    }
}

template<typename TF>
void Boundary_surface<TF>::backward_device()
{
    auto& gd = grid.get_grid_data();

    const int dmemsize2d = gd.ijcells * sizeof(TF);

    // 2D fields:
    cuda_safe_call(cudaMemcpy(obuk.data(),    obuk_g,    dmemsize2d, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(ustar.data(),   ustar_g,   dmemsize2d, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(dudz_mo.data(), dudz_mo_g, dmemsize2d, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(dvdz_mo.data(), dvdz_mo_g, dmemsize2d, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(dbdz_mo.data(), dbdz_mo_g, dmemsize2d, cudaMemcpyDeviceToHost));
}

template<typename TF>
void Boundary_surface<TF>::clear_device()
{
    cuda_safe_call(cudaFree(obuk_g ));
    cuda_safe_call(cudaFree(ustar_g));
    cuda_safe_call(cudaFree(z0m_g));
    cuda_safe_call(cudaFree(z0h_g));

    cuda_safe_call(cudaFree(dudz_mo_g));
    cuda_safe_call(cudaFree(dvdz_mo_g));
    cuda_safe_call(cudaFree(dbdz_mo_g));

    if (sw_constant_z0)
    {
        cuda_safe_call(cudaFree(nobuk_g));
        cuda_safe_call(cudaFree(zL_sl_g));
        cuda_safe_call(cudaFree(f_sl_g ));
    }
}

#ifdef USECUDA
template<typename TF>
void Boundary_surface<TF>::exec(
        Thermo<TF>& thermo, Radiation<TF>& radiation,
        Microphys<TF>& microphys, Timeloop<TF>& timeloop)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;

    // For 2D field excluding ghost cells
    int gridi = gd.imax/blocki + (gd.imax%blocki > 0);
    int gridj = gd.jmax/blockj + (gd.jmax%blockj > 0);
    dim3 gridGPU (gridi,  gridj,  1);
    dim3 blockGPU(blocki, blockj, 1);

    // For 2D field including ghost cells
    gridi = gd.icells/blocki + (gd.icells%blocki > 0);
    gridj = gd.jcells/blockj + (gd.jcells%blockj > 0);
    dim3 gridGPU2 (gridi,  gridj,  1);
    dim3 blockGPU2(blocki, blockj, 1);

    // Calculate dutot in tmp2
    auto dutot = fields.get_tmp_g();

    bsk::calc_dutot_g<<<gridGPU, blockGPU>>>(
        dutot->fld_g,
        fields.mp.at("u")->fld_g,
        fields.mp.at("v")->fld_g,
        fields.mp.at("u")->fld_bot_g,
        fields.mp.at("v")->fld_bot_g,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // 2D cyclic boundaries on dutot
    boundary_cyclic.exec_2d_g(dutot->fld_g);

    // start with retrieving the stability information
    if (thermo.get_switch() == "0")
    {
        // Calculate ustar and Obukhov length, including ghost cells
        stability_neutral_g<<<gridGPU2, blockGPU2>>>(
            ustar_g, obuk_g,
            dutot->fld_g, z0m_g, gd.z[gd.kstart],
            gd.icells, gd.jcells, gd.icells,
            mbcbot, thermobc);
        cuda_check_error();
    }
    else
    {
        auto buoy = fields.get_tmp_g();
        thermo.get_buoyancy_surf_g(*buoy);
        const TF db_ref = thermo.get_db_ref();

        // Calculate ustar and Obukhov length, including ghost cells
        stability_g<<<gridGPU2, blockGPU2>>>(
            ustar_g, obuk_g, nobuk_g,
            buoy->fld_g, buoy->fld_bot_g, buoy->flux_bot_g,
            dutot->fld_g, z0m_g,
            zL_sl_g, f_sl_g,
            db_ref, gd.z[gd.kstart],
            gd.icells, gd.jcells,
            gd.kstart, gd.icells,
            gd.ijcells,
            mbcbot, thermobc);
        cuda_check_error();

        fields.release_tmp_g(buoy);
    }

    fields.release_tmp_g(dutot);

    // Calculate the surface value, gradient and flux depending on the chosen boundary condition.
    // Momentum:
    surfm_flux_g<<<gridGPU, blockGPU>>>(
        fields.mp.at("u")->flux_bot_g,
        fields.mp.at("v")->flux_bot_g,
        fields.mp.at("u")->fld_g,
        fields.mp.at("v")->fld_g,
        fields.mp.at("u")->fld_bot_g,
        fields.mp.at("v")->fld_bot_g,
        ustar_g, obuk_g, z0m_g, gd.z[gd.kstart],
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart,
        gd.icells, gd.ijcells,
        mbcbot);
    cuda_check_error();

    // 2D cyclic boundaries on the surface fluxes
    boundary_cyclic.exec_2d_g(fields.mp.at("u")->flux_bot_g);
    boundary_cyclic.exec_2d_g(fields.mp.at("v")->flux_bot_g);

    // Calculate surface gradients, including ghost cells
    surfm_grad_g<<<gridGPU2, blockGPU2>>>(
        fields.mp.at("u")->grad_bot_g,
        fields.mp.at("v")->grad_bot_g,
        fields.mp.at("u")->fld_g,
        fields.mp.at("v")->fld_g,
        fields.mp.at("u")->fld_bot_g,
        fields.mp.at("v")->fld_bot_g,
        gd.z[gd.kstart], gd.icells, gd.jcells,
        gd.kstart, gd.icells, gd.ijcells);
    cuda_check_error();

    // Scalars:
    for (auto it : fields.sp)
        surfs_g<<<gridGPU2, blockGPU2>>>(
            it.second->flux_bot_g,
            it.second->grad_bot_g,
            it.second->fld_bot_g,
            it.second->fld_g,
            ustar_g, obuk_g, z0h_g, gd.z[gd.kstart],
            gd.icells,  gd.jcells, gd.kstart,
            gd.icells, gd.ijcells,
            sbc.at(it.first).bcbot);
    cuda_check_error();

    // Calc MO gradients, for subgrid scheme
    bsk::calc_duvdz_mo_g<<<gridGPU2, blockGPU2>>>(
            dudz_mo_g, dvdz_mo_g,
            fields.mp.at("u")->fld_g,
            fields.mp.at("v")->fld_g,
            fields.mp.at("u")->fld_bot_g,
            fields.mp.at("v")->fld_bot_g,
            fields.mp.at("u")->flux_bot_g,
            fields.mp.at("v")->flux_bot_g,
            ustar_g, obuk_g, z0m_g,
            gd.z[gd.kstart],
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart,
            gd.icells, gd.ijcells);
    cuda_check_error();

    auto buoy = fields.get_tmp_g();
    thermo.get_buoyancy_fluxbot_g(*buoy);

    bsk::calc_dbdz_mo_g<<<gridGPU2, blockGPU2>>>(
            dbdz_mo_g, buoy->flux_bot_g,
            ustar_g, obuk_g,
            gd.z[gd.kstart],
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
    cuda_check_error();

    fields.release_tmp_g(buoy);
}

template<typename TF>
void Boundary_surface<TF>::exec_column(Column<TF>& column)
{
    const TF no_offset = 0.;
    column.calc_time_series("obuk", obuk_g, no_offset);
    column.calc_time_series("ustar", ustar_g, no_offset);
}
#endif

template class Boundary_surface<double>;
template class Boundary_surface<float>;
