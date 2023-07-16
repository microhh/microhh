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
#include <cmath>
#include <algorithm>
#include <iostream>

#include "grid.h"
#include "fields.h"
#include "master.h"
#include "diff_smag2.h"
#include "boundary.h"
#include "boundary_surface.h"
#include "defines.h"
#include "constants.h"
#include "thermo.h"
#include "tools.h"
#include "stats.h"
#include "monin_obukhov.h"
#include "fast_math.h"

#include "diff_kernels.cuh"

namespace
{
    namespace most = Monin_obukhov;
    namespace fm = Fast_math;
    namespace dk = Diff_kernels_g;

    template<typename TF, Surface_model surface_model> __global__
    void evisc_g(
            TF* __restrict__ evisc,
            TF* __restrict__ N2,
            TF* __restrict__ bgradbot,
            TF* __restrict__ mlen0,
            TF* __restrict__ z0m,
            TF* __restrict__ z,
            const TF tPri,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj,     const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const TF n_mason = TF(2);

        if (i < iend && j < jend && k < kend)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + k*kk;

            if (k == kstart && surface_model == Surface_model::Enabled)
            {
                // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
                TF RitPrratio = bgradbot[ij] / evisc[ijk] * tPri;
                RitPrratio = fmin(RitPrratio, TF(1.-Constants::dsmall));

                const TF mlen = std::pow(TF(1.)/(TF(1.)/mlen0[k] + TF(1.)/(std::pow(Constants::kappa<TF>*(z[kstart]+z0m[ij]), n_mason))), TF(1.)/n_mason);
                evisc[ijk] = fm::pow2(mlen) * sqrt(evisc[ijk] * (TF(1.)-RitPrratio));
            }
            else if (surface_model == Surface_model::Enabled)
            {
                // Add the buoyancy production to the TKE
                TF RitPrratio = N2[ijk] / evisc[ijk] * tPri;
                RitPrratio = fmin(RitPrratio, TF(1.-Constants::dsmall));

                // Mason mixing length
                const TF mlen = std::pow(TF(1.)/(TF(1.)/mlen0[k] + TF(1.)/(std::pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);
                evisc[ijk] = fm::pow2(mlen) * sqrt(evisc[ijk] * (TF(1.)-RitPrratio));
            }
            else
            {
                // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
                TF RitPrratio = N2[ijk] / evisc[ijk] * tPri;
                RitPrratio = fmin(RitPrratio, TF(1.-Constants::dsmall));
                evisc[ijk] = fm::pow2(mlen0[k]) * sqrt(evisc[ijk] * (TF(1.)-RitPrratio));
            }
        }
    }

    template<typename TF> __global__
    void evisc_neutral_g(
            TF* __restrict__ evisc,
            TF* __restrict__ z0m,
            TF* __restrict__ z,
            TF* __restrict__ mlen0,
            const int istart, const int jstart, const int kstart,
            const int iend, const int jend, const int kend,
            const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const TF n_mason = TF(2);

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            const int ij = i + j*jj;

            const TF mlen = std::pow(TF(1.)/(TF(1.)/mlen0[k] + TF(1.)/(std::pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);
            evisc[ijk] = fm::pow2(mlen) * sqrt(evisc[ijk]);
        }
    }

    template<typename TF> __global__
    void evisc_neutral_vandriest_g(
            TF* __restrict__ evisc,
            const TF* __restrict__ u, const TF* __restrict__ v,
            const TF* __restrict__ mlen_smag,
            const TF* __restrict__ z, const TF* __restrict__ dzhi,
            const TF zsize, const TF visc,
            const int istart, const int jstart, const int kstart,
            const int iend, const int jend, const int kend,
            const int jj, const int kk)

    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const TF A_vandriest = TF(26.);

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            const int ijk_bot = i + j*jj + kstart*kk;
            const int ijk_top = i + j*jj + kend*kk;

            const TF u_tau_bot = pow(
                    fm::pow2( visc*(u[ijk_bot] - u[ijk_bot-kk] )*dzhi[kstart] )
                  + fm::pow2( visc*(v[ijk_bot] - v[ijk_bot-kk] )*dzhi[kstart] ), TF(0.25) );
            const TF u_tau_top = pow(
                    fm::pow2( visc*(u[ijk_top] - u[ijk_top-kk] )*dzhi[kend] )
                  + fm::pow2( visc*(v[ijk_top] - v[ijk_top-kk] )*dzhi[kend] ), TF(0.25) );

            const TF fac_bot = TF(1.) - exp( -(       z[k] *u_tau_bot) / (A_vandriest*visc) );
            const TF fac_top = TF(1.) - exp( -((zsize-z[k])*u_tau_top) / (A_vandriest*visc) );

            const TF fac = min(fac_bot, fac_top);

            evisc[ijk] = fm::pow2(fac * mlen_smag[k]) * sqrt(evisc[ijk]);
        }
    }

    template<typename TF> __global__
    void calc_ghostcells_evisc(
            TF* __restrict__ evisc,
            const int icells, const int jcells,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int kb = kstart;
            const int kt = kend-1;

            const int ijkb = i + j*jj + kb*kk;
            const int ijkt = i + j*jj + kt*kk;

            evisc[ijkb-kk] = evisc[ijkb];
            evisc[ijkt+kk] = evisc[ijkt];
        }
    }
}

/* Calculate the mixing length (mlen) offline, and put on GPU */
#ifdef USECUDA
template<typename TF>
void Diff_smag2<TF>::prepare_device(Boundary<TF>& boundary)
{
    auto& gd = grid.get_grid_data();

    std::vector<TF> mlen(gd.kcells);

    if (boundary.get_switch() == "default")
    {
        for (int k=0; k<gd.kcells; ++k)
            mlen[k] = cs * pow(gd.dx*gd.dy*gd.dz[k], 1./3.);
    }
    else
    {
        const TF n_mason = TF(2);
        for (int k=0; k<gd.kcells; ++k)
            mlen[k] = std::pow(cs * std::pow(gd.dx*gd.dy*gd.dz[k], TF(1./3.)), n_mason);
    }

    const int nmemsize = gd.kcells*sizeof(TF);
    cuda_safe_call(cudaMalloc(&mlen_g, nmemsize));
    cuda_safe_call(cudaMemcpy(mlen_g, mlen.data(), nmemsize, cudaMemcpyHostToDevice));
}

template<typename TF>
void Diff_smag2<TF>::clear_device()
{
    cuda_safe_call(cudaFree(mlen_g));
}
#endif

#ifdef USECUDA
template<typename TF>
void Diff_smag2<TF>::exec_viscosity(Stats<TF>&, Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    // Contain the full icells and jcells in this grid.
    const int grid2di  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int grid2dj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 grid2dGPU (grid2di, grid2dj);
    dim3 block2dGPU(blocki, blockj);

    // Use surface model.
    if (boundary.get_switch() != "default")
    {
        TF* z0m_g   = boundary.get_z0m_g();

        // Get MO gradients velocity:
        TF* dudz_g  = boundary.get_dudz_g();
        TF* dvdz_g  = boundary.get_dvdz_g();

        // Calculate total strain rate
        dk::calc_strain2_g<TF, Surface_model::Enabled><<<gridGPU, blockGPU>>>(
            fields.sd.at("evisc")->fld_g,
            fields.mp.at("u")->fld_g,
            fields.mp.at("v")->fld_g,
            fields.mp.at("w")->fld_g,
            dudz_g, dvdz_g,
            gd.dzi_g, gd.dzhi_g,
            gd.dxi, gd.dyi,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);
        cuda_check_error();

        if (thermo.get_switch() == "0")
        {
            // Start with retrieving the stability information
            evisc_neutral_g<TF><<<gridGPU, blockGPU>>>(
                fields.sd.at("evisc")->fld_g,
                z0m_g, gd.z_g, mlen_g,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();
        }
        else
        {
            // Assume buoyancy calculation is needed
            auto tmp1 = fields.get_tmp_g();
            thermo.get_thermo_field_g(*tmp1, "N2", false);

            // Get MO gradient buoyancy:
            TF* dbdz_g  = boundary.get_dbdz_g();

            // Calculate eddy viscosity
            TF tPri = 1./tPr;

            evisc_g<TF, Surface_model::Enabled><<<gridGPU, blockGPU>>>(
                fields.sd.at("evisc")->fld_g,
                tmp1->fld_g, dbdz_g,
                mlen_g, z0m_g, gd.z_g,
                tPri,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();

            fields.release_tmp_g(tmp1);
        }

        boundary_cyclic.exec_g(fields.sd.at("evisc")->fld_g);
    }
    // Do not use surface model.
    else
    {
        // Calculate total strain rate
        dk::calc_strain2_g<TF, Surface_model::Disabled><<<gridGPU, blockGPU>>>(
            fields.sd.at("evisc")->fld_g,
            fields.mp.at("u")->fld_g,
            fields.mp.at("v")->fld_g,
            fields.mp.at("w")->fld_g,
            nullptr, nullptr,
            gd.dzi_g, gd.dzhi_g,
            gd.dxi, gd.dyi,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);
        cuda_check_error();

        // start with retrieving the stability information
        if (thermo.get_switch() == "0")
        {
            evisc_neutral_vandriest_g<TF><<<gridGPU, blockGPU>>>(
                fields.sd.at("evisc")->fld_g,
                fields.mp.at("u")->fld_g,
                fields.mp.at("v")->fld_g,
                mlen_g, gd.z_g, gd.dzhi_g,
                gd.zsize, fields.visc,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend, gd.jend, gd.kend,
                gd.icells, gd.ijcells);
            cuda_check_error();
        }
        // assume buoyancy calculation is needed
        else
        {
            // store the buoyancyflux in datafluxbot of tmp1
            auto tmp1 = fields.get_tmp_g();
            thermo.get_buoyancy_fluxbot_g(*tmp1);
            // As we only use the fluxbot field of tmp1 we store the N2 in the interior.
            thermo.get_thermo_field_g(*tmp1, "N2", false);

            // Calculate eddy viscosity
            TF tPri = 1./tPr;

            evisc_g<TF, Surface_model::Enabled><<<gridGPU, blockGPU>>>(
                fields.sd.at("evisc")->fld_g,
                tmp1->fld_g, nullptr,
                mlen_g, nullptr, gd.z_g,
                tPri,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

            cuda_check_error();

            fields.release_tmp_g(tmp1);
        }

        boundary_cyclic.exec_g(fields.sd.at("evisc")->fld_g);
        calc_ghostcells_evisc<TF><<<grid2dGPU, block2dGPU>>>(
                fields.sd.at("evisc")->fld_g,
                gd.icells, gd.jcells,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }
}
#endif

#ifdef USECUDA
template<typename TF>
void Diff_smag2<TF>::exec(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    const TF dxidxi = TF(1)/(gd.dx * gd.dx);
    const TF dyidyi = TF(1)/(gd.dy * gd.dy);
    const TF tPri = TF(1)/tPr;

    // Do not use surface model.
    if (boundary.get_switch() == "default")
    {
        dk::diff_uvw_g<TF, Surface_model::Disabled><<<gridGPU, blockGPU>>>(
                fields.mt.at("u")->fld_g,
                fields.mt.at("v")->fld_g,
                fields.mt.at("w")->fld_g,
                fields.sd.at("evisc")->fld_g,
                fields.mp.at("u")->fld_g,
                fields.mp.at("v")->fld_g,
                fields.mp.at("w")->fld_g,
                fields.mp.at("u")->flux_bot_g,
                fields.mp.at("u")->flux_top_g,
                fields.mp.at("v")->flux_bot_g,
                fields.mp.at("v")->flux_top_g,
                gd.dzi_g,
                gd.dzhi_g,
                fields.rhoref_g,
                fields.rhorefh_g,
                gd.dxi,
                gd.dyi,
                fields.visc,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        cuda_check_error();

        for (auto it : fields.st)
        {
            dk::diff_c_g<TF, Surface_model::Disabled><<<gridGPU, blockGPU>>>(
                    it.second->fld_g,
                    fields.sp.at(it.first)->fld_g,
                    fields.sd.at("evisc")->fld_g,
                    fields.sp.at(it.first)->flux_bot_g,
                    fields.sp.at(it.first)->flux_top_g,
                    gd.dzi_g,
                    gd.dzhi_g,
                    fields.rhoref_g,
                    fields.rhorefh_g,
                    dxidxi,
                    dyidyi,
                    tPri,
                    fields.sp.at(it.first)->visc,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
        }
        cuda_check_error();
    }
    // Use surface model.
    else
    {
        dk::diff_uvw_g<TF, Surface_model::Enabled><<<gridGPU, blockGPU>>>(
                fields.mt.at("u")->fld_g,
                fields.mt.at("v")->fld_g,
                fields.mt.at("w")->fld_g,
                fields.sd.at("evisc")->fld_g,
                fields.mp.at("u")->fld_g,
                fields.mp.at("v")->fld_g,
                fields.mp.at("w")->fld_g,
                fields.mp.at("u")->flux_bot_g,
                fields.mp.at("u")->flux_top_g,
                fields.mp.at("v")->flux_bot_g,
                fields.mp.at("v")->flux_top_g,
                gd.dzi_g,
                gd.dzhi_g,
                fields.rhoref_g,
                fields.rhorefh_g,
                gd.dxi,
                gd.dyi,
                fields.visc,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
        cuda_check_error();

        for (auto it : fields.st)
            dk::diff_c_g<TF, Surface_model::Enabled><<<gridGPU, blockGPU>>>(
                    it.second->fld_g,
                    fields.sp.at(it.first)->fld_g,
                    fields.sd.at("evisc")->fld_g,
                    fields.sp.at(it.first)->flux_bot_g,
                    fields.sp.at(it.first)->flux_top_g,
                    gd.dzi_g,
                    gd.dzhi_g,
                    fields.rhoref_g,
                    fields.rhorefh_g,
                    dxidxi,
                    dyidyi,
                    tPri,
                    fields.sp.at(it.first)->visc,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
        cuda_check_error();
    }

    cudaDeviceSynchronize();
    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);
    for (auto it : fields.st)
        stats.calc_tend(*it.second, tend_name);

}
#endif

#ifdef USECUDA
template<typename TF>
unsigned long Diff_smag2<TF>::get_time_limit(unsigned long idt, double dt)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    const TF dxidxi = TF(1)/(gd.dx * gd.dx);
    const TF dyidyi = TF(1)/(gd.dy * gd.dy);
    const TF tPrfac_i = TF(1)/std::min(TF(1.), tPr);

    auto tmp1 = fields.get_tmp_g();

    // Calculate dnmul in tmp1 field
    dk::calc_dnmul_g<<<gridGPU, blockGPU>>>(
            tmp1->fld_g,
            fields.sd.at("evisc")->fld_g,
            gd.dzi_g,
            tPrfac_i,
            dxidxi, dyidyi,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.icells, gd.ijcells);
    cuda_check_error();

    // Get maximum from tmp1 field
    double dnmul = field3d_operators.calc_max_g(tmp1->fld_g);
    dnmul = std::max(Constants::dsmall, dnmul);

    const unsigned long idtlim = idt * dnmax/(dnmul*dt);

    fields.release_tmp_g(tmp1);

    return idtlim;
}
#endif

#ifdef USECUDA
template<typename TF>
double Diff_smag2<TF>::get_dn(double dt)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    const TF dxidxi = TF(1)/(gd.dx * gd.dx);
    const TF dyidyi = TF(1)/(gd.dy * gd.dy);
    const TF tPrfac_i = TF(1)/std::min(TF(1.), tPr);

    // Calculate dnmul in tmp1 field
    auto dnmul_tmp = fields.get_tmp_g();

    dk::calc_dnmul_g<<<gridGPU, blockGPU>>>(
        dnmul_tmp->fld_g,
        fields.sd.at("evisc")->fld_g,
        gd.dzi_g,
        tPrfac_i,
        dxidxi, dyidyi,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart, gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // Get maximum from tmp1 field
    // CvH This is odd, because there might be need for calc_max in CPU version.
    double dnmul = field3d_operators.calc_max_g(dnmul_tmp->fld_g);

    fields.release_tmp_g(dnmul_tmp);

    return dnmul*dt;
}
#endif

template class Diff_smag2<double>;
template class Diff_smag2<float>;
