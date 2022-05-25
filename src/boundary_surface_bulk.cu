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

#include <cmath>
#include "fast_math.h"
#include "constants.h"
#include "tools.h"
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "boundary_surface_bulk.h"
#include "boundary_surface_kernels_gpu.h"
#include "monin_obukhov.h"

namespace
{
    namespace most = Monin_obukhov;
    namespace fm = Fast_math;
    namespace bsk = Boundary_surface_kernels_g;

    template<typename TF> __global__
    void momentum_fluxgrad_g(
            TF* const __restrict__ ufluxbot,
            TF* const __restrict__ vfluxbot,
            TF* const __restrict__ ugradbot,
            TF* const __restrict__ vgradbot,
            const TF* const __restrict__ u,
            const TF* const __restrict__ v,
            const TF* const __restrict__ ubot,
            const TF* const __restrict__ vbot,
            const TF* const __restrict__ dutot,
            const TF Cm,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;

            ufluxbot[ij] = -Cm * dutot[ij] * (u[ijk]-ubot[ij]);
            vfluxbot[ij] = -Cm * dutot[ij] * (v[ijk]-vbot[ij]);

            ugradbot[ij] = (u[ijk]-ubot[ij])/zsl;
            vgradbot[ij] = (v[ijk]-vbot[ij])/zsl;
        }
    }

    template<typename TF> __global__
    void scalar_fluxgrad_g(
            TF* const __restrict__ sfluxbot,
            TF* const __restrict__ sgradbot,
            const TF* const __restrict__ s,
            const TF* const __restrict__ sbot,
            const TF* const __restrict__ dutot,
            const TF Cs,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;

            sfluxbot[ij] = -Cs * dutot[ij] * (s[ijk]-sbot[ij]);
            sgradbot[ij] = (s[ijk]-sbot[ij])/zsl;
        }
    }

    template<typename TF> __global__
    void surface_scaling_g(
            TF* const __restrict__ ustar,
            TF* const __restrict__ obuk,
            const TF* const __restrict__ dutot,
            const TF* const __restrict__ bfluxbot,
            const TF Cm,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int jj)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij = i + j*jj;

            ustar[ij] = sqrt(Cm) * dutot[ij];
            obuk[ij] = - fm::pow3(ustar[ij]) / (Constants::kappa<TF> * bfluxbot[ij]);
        }
    }
}

template<typename TF>
void Boundary_surface_bulk<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();

    // Prepare base boundary, for inflow profiles.
    Boundary<TF>::prepare_device();

    const int dmemsize2d = gd.ijcells*sizeof(TF);
    const int dimemsize  = gd.icells*sizeof(TF);

    cuda_safe_call(cudaMalloc(&obuk_g,  dmemsize2d));
    cuda_safe_call(cudaMalloc(&ustar_g, dmemsize2d));

    cuda_safe_call(cudaMalloc(&z0m_g, dmemsize2d));

    cuda_safe_call(cudaMalloc(&dudz_mo_g, dmemsize2d));
    cuda_safe_call(cudaMalloc(&dvdz_mo_g, dmemsize2d));
    cuda_safe_call(cudaMalloc(&dbdz_mo_g, dmemsize2d));

    cuda_safe_call(cudaMemcpy2D(obuk_g,  dimemsize, obuk.data(),  dimemsize, dimemsize, gd.jcells, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy2D(ustar_g, dimemsize, ustar.data(), dimemsize, dimemsize, gd.jcells, cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy2D(z0m_g, dimemsize, z0m.data(), dimemsize, dimemsize, gd.jcells, cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy2D(dudz_mo_g, dimemsize, dudz_mo.data(), dimemsize, dimemsize, gd.jcells, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy2D(dvdz_mo_g, dimemsize, dvdz_mo.data(), dimemsize, dimemsize, gd.jcells, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy2D(dbdz_mo_g, dimemsize, dbdz_mo.data(), dimemsize, dimemsize, gd.jcells, cudaMemcpyHostToDevice));
}

// TMP BVS
template<typename TF>
void Boundary_surface_bulk<TF>::forward_device()
{
    auto& gd = grid.get_grid_data();

    const int dimemsize   = gd.icells  * sizeof(TF);

    cuda_safe_call(cudaMemcpy2D(obuk_g,  dimemsize, obuk.data(),  dimemsize, dimemsize, gd.jcells, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy2D(ustar_g, dimemsize, ustar.data(), dimemsize, dimemsize, gd.jcells, cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy2D(dudz_mo_g, dimemsize, dudz_mo.data(), dimemsize, dimemsize, gd.jcells, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy2D(dvdz_mo_g, dimemsize, dvdz_mo.data(), dimemsize, dimemsize, gd.jcells, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy2D(dbdz_mo_g, dimemsize, dbdz_mo.data(), dimemsize, dimemsize, gd.jcells, cudaMemcpyHostToDevice));
}

// TMP BVS
template<typename TF>
void Boundary_surface_bulk<TF>::backward_device()
{
    auto& gd = grid.get_grid_data();

    const int dimemsize = gd.icells * sizeof(TF);

    cuda_safe_call(cudaMemcpy2D(obuk.data(),  dimemsize, obuk_g,  dimemsize, dimemsize, gd.jcells, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy2D(ustar.data(), dimemsize, ustar_g, dimemsize, dimemsize, gd.jcells, cudaMemcpyDeviceToHost));

    cuda_safe_call(cudaMemcpy2D(dudz_mo.data(), dimemsize, dudz_mo_g, dimemsize, dimemsize, gd.jcells, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy2D(dvdz_mo.data(), dimemsize, dvdz_mo_g, dimemsize, dimemsize, gd.jcells, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy2D(dbdz_mo.data(), dimemsize, dbdz_mo_g, dimemsize, dimemsize, gd.jcells, cudaMemcpyDeviceToHost));
}

template<typename TF>
void Boundary_surface_bulk<TF>::clear_device()
{
    cuda_safe_call(cudaFree(obuk_g ));
    cuda_safe_call(cudaFree(ustar_g));

    cuda_safe_call(cudaFree(z0m_g));

    cuda_safe_call(cudaFree(dudz_mo_g));
    cuda_safe_call(cudaFree(dvdz_mo_g));
    cuda_safe_call(cudaFree(dbdz_mo_g));
}

#ifdef USECUDA
template<typename TF>
void Boundary_surface_bulk<TF>::exec(
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

    const TF zsl = gd.z[gd.kstart];

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

    // Calculate surface momentum fluxes, excluding ghost cells
    momentum_fluxgrad_g<<<gridGPU, blockGPU>>>(
        fields.mp.at("u")->flux_bot_g,
        fields.mp.at("v")->flux_bot_g,
        fields.mp.at("u")->grad_bot_g,
        fields.mp.at("v")->grad_bot_g,
        fields.mp.at("u")->fld_g,
        fields.mp.at("v")->fld_g,
        fields.mp.at("u")->fld_bot_g,
        fields.mp.at("v")->fld_bot_g,
        dutot->fld_g, bulk_cm, zsl,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // 2D cyclic boundaries on the surface fluxes
    boundary_cyclic.exec_2d_g(fields.mp.at("u")->flux_bot_g);
    boundary_cyclic.exec_2d_g(fields.mp.at("v")->flux_bot_g);
    boundary_cyclic.exec_2d_g(fields.mp.at("u")->grad_bot_g);
    boundary_cyclic.exec_2d_g(fields.mp.at("v")->grad_bot_g);

    // Calculate scalar fluxes, gradients and/or values, including ghost cells
    for (auto it : fields.sp)
    {
        scalar_fluxgrad_g<<<gridGPU2, blockGPU2>>>(
            it.second->flux_bot_g,
            it.second->grad_bot_g,
            it.second->fld_g,
            it.second->fld_bot_g,
            dutot->fld_g,
            bulk_cs.at(it.first), zsl,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart,
            gd.icells, gd.ijcells);
        cuda_check_error();

        boundary_cyclic.exec_2d_g(it.second->flux_bot_g);
        boundary_cyclic.exec_2d_g(it.second->grad_bot_g);
    }

    // Calculate ustar and Obukhov length
    auto b= fields.get_tmp_g();
    thermo.get_buoyancy_fluxbot_g(*b);

    surface_scaling_g<<<gridGPU2, blockGPU2>>>(
        ustar_g,
        obuk_g,
        dutot->fld_g,
        b->flux_bot_g,
        bulk_cm,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.icells);

    // Calculate MO gradients for diffusion scheme
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

    bsk::calc_dbdz_mo_g<<<gridGPU2, blockGPU2>>>(
        dbdz_mo_g, b->flux_bot_g,
        ustar_g, obuk_g,
        gd.z[gd.kstart],
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.icells);
    cuda_check_error();

    fields.release_tmp_g(b);
    fields.release_tmp_g(dutot);
}
#endif

template class Boundary_surface_bulk<double>;
template class Boundary_surface_bulk<float>;
