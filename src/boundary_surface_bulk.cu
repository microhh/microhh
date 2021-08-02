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
#include "monin_obukhov.h"

namespace
{
    namespace most = Monin_obukhov;
    namespace fm = Fast_math;

    /* Calculate absolute wind speed */
    template<typename TF> __global__
    void calculate_du_g(
            TF* __restrict__ dutot,
            TF* __restrict__ u,    TF* __restrict__ v,
            TF* __restrict__ ubot, TF* __restrict__ vbot,
            const int istart, const int iend, const int jstart, const int jend, const int kstart,
            const int jj, const int kk)

    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ii  = 1;
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            const TF minval = 1.e-1;

            const TF du2 = fm::pow2(TF(0.5)*(u[ijk] + u[ijk+ii]) - TF(0.5)*(ubot[ij] + ubot[ij+ii]))
                         + fm::pow2(TF(0.5)*(v[ijk] + v[ijk+jj]) - TF(0.5)*(vbot[ij] + vbot[ij+jj]));
            dutot[ij] = fmax(sqrt(du2), minval);
        }
    }

    template<typename TF> __global__
    void momentum_fluxgrad_g(
            TF* __restrict__ ufluxbot, TF* __restrict__ vfluxbot,
            TF* __restrict__ ugradbot, TF* __restrict__ vgradbot,
            TF* __restrict__ u,        TF* __restrict__ v,
            TF* __restrict__ ubot,     TF* __restrict__ vbot,
            TF* __restrict__ dutot,    TF Cm, TF zsl,
            const int istart, const int iend, const int jstart, const int jend, const int kstart,
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
            TF* __restrict__ sfluxbot,
            TF* __restrict__ sgradbot,
            TF* __restrict__ s,
            TF* __restrict__ sbot,
            TF* __restrict__ dutot,    TF Cs, TF zsl,
            const int istart, const int iend, const int jstart, const int jend, const int kstart,
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
            TF* __restrict__ ustar,
            TF* __restrict__ obuk,
            TF* __restrict__ dutot,
            TF* __restrict__ bfluxbot,
            TF Cm,
            const int istart, const int iend, const int jstart, const int jend,
            const int jj)

    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij = i + j*jj;
            const TF sqrt_Cm = sqrt(Cm);

            ustar[ij] = sqrt_Cm * dutot[ij];
            obuk[ij] = - fm::pow3(ustar[ij]) / (Constants::kappa<TF> * bfluxbot[ij]);
        }
    }

template<typename TF> __global__
    void calc_duvdz_g(
            TF* const __restrict__ dudz,
            TF* const __restrict__ dvdz,
            const TF* const __restrict__ u,
            const TF* const __restrict__ v,
            const TF* const __restrict__ ubot,
            const TF* const __restrict__ vbot,
            const TF* const __restrict__ ufluxbot,
            const TF* const __restrict__ vfluxbot,
            const TF* const __restrict__ ustar,
            const TF* const __restrict__ obuk,
            const TF* const __restrict__ z0m,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        const int ii = 1;
        const int jj = icells;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*icells;
            const int ijk = i + j*icells + kstart*ijcells;

            const TF du_c = TF(0.5)*((u[ijk] - ubot[ij]) + (u[ijk+ii] - ubot[ij+ii]));
            const TF dv_c = TF(0.5)*((v[ijk] - vbot[ij]) + (v[ijk+jj] - vbot[ij+jj]));

            const TF ufluxbot = -du_c * ustar[ij] * most::fm(zsl, z0m[ij], obuk[ij]);
            const TF vfluxbot = -dv_c * ustar[ij] * most::fm(zsl, z0m[ij], obuk[ij]);

            dudz[ij] = -ufluxbot / (Constants::kappa<TF> * zsl * ustar[ij]) * most::phim(zsl/obuk[ij]);
            dvdz[ij] = -vfluxbot / (Constants::kappa<TF> * zsl * ustar[ij]) * most::phim(zsl/obuk[ij]);
        }
    }

    template<typename TF> __global__
    void calc_dbdz_g(
            TF* const __restrict__ dbdz,
            const TF* const __restrict__ bfluxbot,
            const TF* const __restrict__ ustar,
            const TF* const __restrict__ obuk,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        const int ii = 1;
        const int jj = icells;

        if (i < iend && j < jend)
        {
            const int ij = i + j*icells;
            dbdz[ij] = -bfluxbot[ij] / (Constants::kappa<TF> * zsl * ustar[ij]) * most::phih(zsl/obuk[ij]);
        }
    }
}

template<typename TF>
void Boundary_surface_bulk<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();

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

    calculate_du_g<<<gridGPU, blockGPU>>>(
        dutot->fld_g,
        fields.mp.at("u")->fld_g,     fields.mp.at("v")->fld_g,
        fields.mp.at("u")->fld_bot_g, fields.mp.at("v")->fld_bot_g,
        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.icells, gd.ijcells);

    cuda_check_error();

    // 2D cyclic boundaries on dutot
    boundary_cyclic.exec_2d_g(dutot->fld_g);

    // Calculate surface momentum fluxes, excluding ghost cells
    momentum_fluxgrad_g<<<gridGPU, blockGPU>>>(
        fields.mp.at("u")->flux_bot_g, fields.mp.at("v")->flux_bot_g,
        fields.mp.at("u")->grad_bot_g, fields.mp.at("v")->grad_bot_g,
        fields.mp.at("u")->fld_g,      fields.mp.at("v")->fld_g,
        fields.mp.at("u")->fld_bot_g,  fields.mp.at("v")->fld_bot_g,
        dutot->fld_g, bulk_cm, zsl,
        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.icells, gd.ijcells);
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
            it.second->flux_bot_g, it.second->grad_bot_g,
            it.second->fld_g, it.second->fld_bot_g,
            dutot->fld_g, bulk_cs.at(it.first), zsl,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.icells, gd.ijcells);
        cuda_check_error();
        boundary_cyclic.exec_2d_g(it.second->flux_bot_g);
        boundary_cyclic.exec_2d_g(it.second->grad_bot_g);
    }

    // Calculate ustar and Obukhov length
    auto b= fields.get_tmp_g();
    thermo.get_buoyancy_fluxbot_g(*b);

    surface_scaling_g<<<gridGPU2, blockGPU2>>>(
        ustar_g, obuk_g, dutot->fld_g, b->flux_bot_g, bulk_cm,
        gd.istart, gd.iend, gd.jstart, gd.jend, gd.icells);

    // Calculate MO gradients for diffusion scheme
    calc_duvdz_g<<<gridGPU2, blockGPU2>>>(
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

    calc_dbdz_g<<<gridGPU2, blockGPU2>>>(
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
