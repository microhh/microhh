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

#include <algorithm>
#include <cmath>
#include <iostream>

#include "grid.h"
#include "fields.h"
#include "master.h"
#include "thermo.h"
#include "boundary.h"
#include "stats.h"
#include "monin_obukhov.h"
#include "tools.h"

#include "diff_tke2.h"
#include "diff_kernels.cuh"

namespace
{
    namespace most = Monin_obukhov;
    namespace fm = Fast_math;
    namespace dk = Diff_kernels_g;

    template<typename TF, Surface_model surface_model, bool sw_mason> __global__
    void calc_evisc_neutral_g(
            TF* const restrict evisc,
            const TF* const restrict sgstke,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict z,
            const TF* const restrict dz,
            const TF* const restrict z0m,
            const TF* const restrict mlen0,
            const TF dx, const TF dy,
            const TF cn, const TF cm,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int jj = icells;
        const int kk = ijcells;

        // Wall damping constant.
        constexpr TF n_mason = TF(2.);
        constexpr TF A_vandriest = TF(26.);

        TF fac;
        if (i < iend && j < jend && k < kend)
        {
            const int ij = i + j*jj;
            const int ijk = i + j*jj + k*kk;

            if (surface_model == Surface_model::Disabled)
                asm("trap;");
            else
            {
                if (sw_mason) // Apply Mason's wall correction
                    fac = pow(TF(1.)/(TF(1.)/pow(mlen0[k], n_mason) + TF(1.)/
                                (pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);
                else
                    fac = mlen0[k];

                // Calculate eddy diffusivity for momentum.
                evisc[ijk] = cm * fac * sqrt(sgstke[ijk]);
            }
        }
    }

    template<typename TF, Surface_model surface_model, bool sw_mason> __global__
    void calc_evisc_g(
            TF* const restrict evisc,
            const TF* const restrict sgstke,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict N2,
            const TF* const restrict bgradbot,
            const TF* const restrict z,
            const TF* const restrict dz,
            const TF* const restrict z0m,
            const TF* const restrict mlen0,
            const TF dx, const TF dy,
            const TF cn, const TF cm,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int jj = icells;
        const int kk = ijcells;

        // Variables for the wall damping and length scales
        const TF n_mason = TF(2.);
        TF mlen;
        TF fac;

        if (i < iend && j < jend && k < kend)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + k*kk;

            if (surface_model == Surface_model::Disabled)
                asm("trap;");
            else
            {
                mlen = mlen0[k];

                if (k == kstart)
                {
                    if ( bgradbot[ij] > 0 ) // Only if stably stratified, adapt length scale
                        mlen = cn * sqrt(sgstke[ijk]) / sqrt(bgradbot[ij]);
                }
                else
                {
                    if ( N2[ijk] > 0 ) // Only if stably stratified, adapt length scale
                        mlen = cn * sqrt(sgstke[ijk]) / sqrt(N2[ijk]);
                }

                fac  = min(mlen0[k], mlen);

                if (sw_mason) // Apply Mason's wall correction here
                    fac = pow(TF(1.)/(TF(1.)/pow(fac, n_mason) + TF(1.)/
                            (pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);

                // Calculate eddy diffusivity for momentum.
                evisc[ijk] = cm * fac * sqrt(sgstke[ijk]);
            }
        }
    }

    template<typename TF, Surface_model surface_model, bool sw_mason> __global__
    void calc_evisc_heat_g(
            TF* const restrict evisch,
            const TF* __restrict__ evisc,
            const TF* __restrict__ sgstke,
            const TF* __restrict__ N2,
            const TF* __restrict__ bgradbot,
            const TF* __restrict__ z,
            const TF* __restrict__ dz,
            const TF* __restrict__ z0m,
            const TF* __restrict__ mlen0,
            const TF dx, const TF dy,
            const TF cn, const TF ch1, const TF ch2,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int jj = icells;
        const int kk = ijcells;

        // Variables for the wall damping and length scales
        const TF n_mason = TF(2.);
        TF mlen;
        TF fac;

        if (i < iend && j < jend && k < kend)
        {

            if (surface_model == Surface_model::Disabled)
                asm("trap;");
            else
            {
                const int ij = i + j*jj;
                const int ijk = i + j*jj + k*kk;

                mlen = mlen0[k];

                if (k == kstart)
                {
                    if ( bgradbot[ij] > 0 ) // Only if stably stratified, adapt length scale
                        mlen = cn * sqrt(sgstke[ijk]) / sqrt(bgradbot[ij]);
                }
                else
                {
                    if ( N2[ijk] > 0 ) // Only if stably stratified, adapt length scale
                        mlen = cn * sqrt(sgstke[ijk]) / sqrt(N2[ijk]);
                }

                fac  = min(mlen0[k], mlen);

                if (sw_mason) // Apply Mason's wall correction here
                    fac = pow(TF(1.)/(TF(1.)/pow(fac, n_mason) + TF(1.)/
                                (pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);

                // Calculate eddy diffusivity for momentum.
                evisch[ijk] = (ch1 + ch2 * fac / mlen0[k]) * evisc[ijk];
            }
        }
    }

    template<typename TF> __global__
    void sgstke_shear_tend_g(
            TF* const __restrict__ at,
            const TF* __restrict__ a,
            const TF* __restrict__ evisc,
            const TF* __restrict__ strain2,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int jj = icells;
        const int kk = ijcells;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;

            // Calculate shear production of SGS TKE based on Deardorff (1980)
            // NOTE: `strain2` is defined/calculated as:
            // S^2 = 0.5 * (dui/dxj + duj/dxi)^2 = dui/dxj * (dui/dxj + duj/dxi)
            at[ijk] += evisc[ijk] * strain2[ijk];
        }
    }

    template<typename TF> __global__
    void sgstke_buoy_tend_g(
            TF* const __restrict__ at,
            const TF* __restrict__ a,
            const TF* __restrict__ evisch,
            const TF* __restrict__ N2,
            const TF* __restrict__ bgradbot,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int jj = icells;
        const int kk = ijcells;

        if (i < iend && j < jend && k < kend)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + k*kk;

            // Calculate buoyancy destruction of SGS TKE based on Deardorff (1980)
            if (k == kstart)
                at[ijk] -= evisch[ijk] * bgradbot[ij];
            else
                at[ijk] -= evisch[ijk] * N2[ijk];

        }
    }

    template<typename TF, bool sw_mason> __global__
    void sgstke_diss_tend_g(
            TF* const __restrict__ at,
            const TF* const __restrict__ a,
            const TF* const __restrict__ N2,
            const TF* const __restrict__ bgradbot,
            const TF* const __restrict__ z,
            const TF* const __restrict__ dz,
            const TF* const __restrict__ z0m,
            const TF* const __restrict__ mlen0,
            const TF dx, const TF dy,
            const TF cn, const TF ce1, const TF ce2,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int jj = icells;
        const int kk = ijcells;

        const TF n_mason = TF(2.);
        TF mlen;
        TF fac;

        if (i < iend && j < jend && k < kend)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + k*kk;

            // Calculate geometric filter width, based on Deardorff (1980)
            mlen = mlen0[k];

            // Only if stably stratified, adapt length scale
            if (k == kstart)
            {
                if (bgradbot[ij] > 0)
                    mlen = cn * sqrt(a[ijk]) / sqrt(bgradbot[ij]);
            }
            else
            {
                if (N2[ijk] > 0)
                    mlen = cn * sqrt(a[ijk]) / sqrt(N2[ijk]);
            }

            fac  = min(mlen0[k], mlen);

            // Apply Mason's wall correction here
            if (sw_mason)
                fac = pow(TF(1.)/(TF(1.)/pow(fac, n_mason) + TF(1.)/
                        (pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);

            // Calculate dissipation of SGS TKE based on Deardorff (1980)
            at[ijk] -= (ce1 + ce2 * fac / mlen0[k]) * pow(a[ijk], TF(3./2.)) / fac;
        }
    }

    template<typename TF,  bool sw_mason> __global__
    void sgstke_diss_tend_neutral_g(
            TF* const __restrict__ at,
            const TF* const __restrict__ a,
            const TF* const __restrict__ z,
            const TF* const __restrict__ dz,
            const TF* const __restrict__ z0m,
            const TF* const __restrict__ mlen0,
            const TF dx, const TF dy,
            const TF ce1, const TF ce2,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int jj = icells;
        const int kk = ijcells;

        const TF n_mason = TF(2.);
        TF fac;

        if (i < iend && j < jend && k < kend)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + k*kk;

            fac = mlen0[k];

            if (sw_mason) // Apply Mason's wall correction here
                fac = pow(TF(1.)/(TF(1.)/pow(mlen0[k], n_mason) + TF(1.)/
                        (pow(Constants::kappa<TF>*(z[k]+z0m[ij]), n_mason))), TF(1.)/n_mason);

            // Calculate dissipation of SGS TKE based on Deardorff (1980)
            at[ijk] -= (ce1 + ce2 * fac / mlen0[k]) * pow(a[ijk], TF(3./2.)) / fac ;
        }
    }
}

#ifdef USECUDA
template<typename TF>
unsigned long Diff_tke2<TF>::get_time_limit(const unsigned long idt, const double dt)
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
    const TF tPr_i_dummy = 1;

    auto tmp1 = fields.get_tmp_g();

    // When no buoyancy, use eddy viscosity for momentum.
    TF* evisc_g = !sw_buoy
        ? fields.sd.at("evisc")->fld_g
        : fields.sd.at("eviscs")->fld_g;

    // Calculate dnmul in tmp1 field
    dk::calc_dnmul_g<TF><<<gridGPU, blockGPU>>>(
            tmp1->fld_g,
            evisc_g,
            gd.dzi_g,
            tPr_i_dummy,
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

template<typename TF>
double Diff_tke2<TF>::get_dn(const double dt)
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
    const TF tPr_i_dummy = 1;

    // Calculate dnmul in tmp1 field
    auto dnmul_tmp = fields.get_tmp_g();

    // When no buoyancy, use eddy viscosity for momentum.
    TF* evisc_g = !sw_buoy
        ? fields.sd.at("evisc")->fld_g
        : fields.sd.at("eviscs")->fld_g;

    dk::calc_dnmul_g<TF><<<gridGPU, blockGPU>>>(
            dnmul_tmp->fld_g,
            evisc_g,
            gd.dzi_g,
            tPr_i_dummy,
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

template<typename TF>
void Diff_tke2<TF>::exec(Stats<TF>& stats)
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

    // Dummy tPr value for `diff_c`.
    const TF tPr_i_dummy = 1;

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
    {
        TF* evisc_ptr;

        if (it.first == "sgstke")  // sgstke diffuses with eddy viscosity for momentum
            evisc_ptr = fields.sd.at("evisc")->fld_g;
        else  // all other scalars, normally diffuse with eddy viscosity for heat/scalars
        {
            if (!sw_buoy) // but not if there is no buoyancy (then eviscs not defined)
                evisc_ptr = fields.sd.at("evisc")->fld_g;
            else
                evisc_ptr = fields.sd.at("eviscs")->fld_g;
        }

        dk::diff_c_g<TF, Surface_model::Enabled><<<gridGPU, blockGPU>>>(
                it.second->fld_g,
                fields.sp.at(it.first)->fld_g,
                evisc_ptr,
                fields.sp.at(it.first)->flux_bot_g,
                fields.sp.at(it.first)->flux_top_g,
                gd.dzi_g,
                gd.dzhi_g,
                fields.rhoref_g,
                fields.rhorefh_g,
                dxidxi,
                dyidyi,
                tPr_i_dummy,
                fields.sp.at(it.first)->visc,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }
    cuda_check_error();

    cudaDeviceSynchronize();
    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);
    for (auto it : fields.st)
        stats.calc_tend(*it.second, tend_name);
}

template<typename TF>
void Diff_tke2<TF>::exec_viscosity(Stats<TF>& stats, Thermo<TF>& thermo)
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

    if (boundary.get_switch() == "default")
        throw std::runtime_error("Resolved wall not supported in Deardorff SGSm.");

    // Get MO gradients velocity and roughness length:
    TF* dudz_g  = boundary.get_dudz_g();
    TF* dvdz_g  = boundary.get_dvdz_g();
    TF* z0m_g   = boundary.get_z0m_g();

    auto str2_tmp = fields.get_tmp_g();

    // Calculate total strain rate
    dk::calc_strain2_g<TF, Surface_model::Enabled><<<gridGPU, blockGPU>>>(
            str2_tmp->fld_g,
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

    // Start with retrieving the stability information
    if (!sw_buoy)
    {
        if (sw_mason)
            calc_evisc_neutral_g<TF, Surface_model::Enabled, true><<<gridGPU, blockGPU>>>(
                    fields.sd.at("evisc")->fld_g,
                    fields.sp.at("sgstke")->fld_g,
                    fields.mp.at("u")->fld_g,
                    fields.mp.at("v")->fld_g,
                    fields.mp.at("w")->fld_g,
                    gd.z_g, gd.dz_g,
                    z0m_g, mlen0_g,
                    gd.dx, gd.dy,
                    this->cn, this->cm,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells,
                    gd.ijcells);
        else
            calc_evisc_neutral_g<TF, Surface_model::Enabled, false><<<gridGPU, blockGPU>>>(
                    fields.sd.at("evisc")->fld_g,
                    fields.sp.at("sgstke")->fld_g,
                    fields.mp.at("u")->fld_g,
                    fields.mp.at("v")->fld_g,
                    fields.mp.at("w")->fld_g,
                    gd.z_g, gd.dz_g,
                    z0m_g, mlen0_g,
                    gd.dx, gd.dy,
                    this->cn, this->cm,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells,
                    gd.ijcells);

        if (sw_mason)
            sgstke_diss_tend_neutral_g<TF, true><<<gridGPU, blockGPU>>>(
                    fields.st.at("sgstke")->fld_g,
                    fields.sp.at("sgstke")->fld_g,
                    gd.z_g, gd.dz_g,
                    z0m_g, mlen0_g,
                    gd.dx, gd.dy,
                    this->ce1, this->ce2,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
        else
            sgstke_diss_tend_neutral_g<TF, false><<<gridGPU, blockGPU>>>(
                    fields.st.at("sgstke")->fld_g,
                    fields.sp.at("sgstke")->fld_g,
                    gd.z_g, gd.dz_g,
                    z0m_g, mlen0_g,
                    gd.dx, gd.dy,
                    this->ce1, this->ce2,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
    }
    else
    {
        // Assume buoyancy calculation is needed
        auto buoy_tmp = fields.get_tmp_g();
        thermo.get_thermo_field_g(*buoy_tmp, "N2", false);

        // Get MO gradient buoyancy:
        TF* dbdz_g  = boundary.get_dbdz_g();

        // Note BvS: templated lambda functions are not (yet?) allowed by NVCC :-(
        if (sw_mason)
            calc_evisc_g<TF, Surface_model::Enabled, true><<<gridGPU, blockGPU>>>(
                    fields.sd.at("evisc")->fld_g,
                    fields.sp.at("sgstke")->fld_g,
                    fields.mp.at("u")->fld_g,
                    fields.mp.at("v")->fld_g,
                    fields.mp.at("w")->fld_g,
                    buoy_tmp->fld_g, dbdz_g,
                    gd.z_g, gd.dz_g,
                    z0m_g, mlen0_g,
                    gd.dx, gd.dy,
                    this->cn, this->cm,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells,
                    gd.ijcells);
        else
            calc_evisc_g<TF, Surface_model::Enabled, false><<<gridGPU, blockGPU>>>(
                    fields.sd.at("evisc")->fld_g,
                    fields.sp.at("sgstke")->fld_g,
                    fields.mp.at("u")->fld_g,
                    fields.mp.at("v")->fld_g,
                    fields.mp.at("w")->fld_g,
                    buoy_tmp->fld_g, dbdz_g,
                    gd.z_g, gd.dz_g,
                    z0m_g, mlen0_g,
                    gd.dx, gd.dy,
                    this->cn, this->cm,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells,
                    gd.ijcells);

        if (sw_mason)
            calc_evisc_heat_g<TF, Surface_model::Enabled, true><<<gridGPU, blockGPU>>>(
                    fields.sd.at("eviscs")->fld_g,
                    fields.sd.at("evisc")->fld_g,
                    fields.sp.at("sgstke")->fld_g,
                    buoy_tmp->fld_g, dbdz_g,
                    gd.z_g, gd.dz_g,
                    z0m_g, mlen0_g,
                    gd.dx, gd.dy,
                    this->cn, this->ch1, this->ch2,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
        else
            calc_evisc_heat_g<TF, Surface_model::Enabled, false><<<gridGPU, blockGPU>>>(
                    fields.sd.at("eviscs")->fld_g,
                    fields.sd.at("evisc")->fld_g,
                    fields.sp.at("sgstke")->fld_g,
                    buoy_tmp->fld_g, dbdz_g,
                    gd.z_g, gd.dz_g,
                    z0m_g, mlen0_g,
                    gd.dx, gd.dy,
                    this->cn, this->ch1, this->ch2,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

        boundary_cyclic.exec_g(fields.sd.at("evisc")->fld_g);
        boundary_cyclic.exec_g(fields.sd.at("eviscs")->fld_g);

        // Calculate tendencies here, to prevent having to
        // re-calculate `strain2` in diffusion->exec()`.
        sgstke_buoy_tend_g<TF><<<gridGPU, blockGPU>>>(
               fields.st.at("sgstke")->fld_g,
               fields.sp.at("sgstke")->fld_g,
               fields.sd.at("eviscs")->fld_g,
               buoy_tmp->fld_g,
               dbdz_g,
               gd.istart, gd.iend,
               gd.jstart, gd.jend,
               gd.kstart, gd.kend,
               gd.icells, gd.ijcells);

        cudaDeviceSynchronize();
        stats.calc_tend(*fields.st.at("sgstke"), tend_name_buoy);

        if (sw_mason)
            sgstke_diss_tend_g<TF, true><<<gridGPU, blockGPU>>>(
                   fields.st.at("sgstke")->fld_g,
                   fields.sp.at("sgstke")->fld_g,
                   buoy_tmp->fld_g,
                   dbdz_g,
                   gd.z_g,
                   gd.dz_g,
                   z0m_g,
                   mlen0_g,
                   gd.dx, gd.dy,
                   this->cn,
                   this->ce1,
                   this->ce2,
                   gd.istart, gd.iend,
                   gd.jstart, gd.jend,
                   gd.kstart, gd.kend,
                   gd.icells, gd.ijcells);
        else
            sgstke_diss_tend_g<TF, false><<<gridGPU, blockGPU>>>(
                    fields.st.at("sgstke")->fld_g,
                    fields.sp.at("sgstke")->fld_g,
                    buoy_tmp->fld_g,
                    dbdz_g,
                    gd.z_g,
                    gd.dz_g,
                    z0m_g,
                    mlen0_g,
                    gd.dx, gd.dy,
                    this->cn,
                    this->ce1,
                    this->ce2,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

        cudaDeviceSynchronize();
        stats.calc_tend(*fields.st.at("sgstke"), tend_name_diss);

        fields.release_tmp_g(buoy_tmp);
    }

    sgstke_shear_tend_g<TF><<<gridGPU, blockGPU>>>(
           fields.st.at("sgstke")->fld_g,
           fields.sp.at("sgstke")->fld_g,
           fields.sd.at("evisc")->fld_g,
           str2_tmp->fld_g,
           gd.istart, gd.iend,
           gd.jstart, gd.jend,
           gd.kstart, gd.kend,
           gd.icells, gd.ijcells);

    cudaDeviceSynchronize();
    stats.calc_tend(*fields.st.at("sgstke"), tend_name_shear);

    // Release temporary fields
    fields.release_tmp_g(str2_tmp);
}

template<typename TF>
void Diff_tke2<TF>::prepare_device(Boundary<TF>& boundary)
{
    auto& gd = grid.get_grid_data();

    std::vector<TF> mlen0(gd.kcells);

    for (int k=0; k<gd.kcells; ++k)
        mlen0[k] = std::pow(gd.dx*gd.dy*gd.dz[k], TF(1./3.));

    const int nmemsize = gd.kcells*sizeof(TF);
    cuda_safe_call(cudaMalloc(&mlen0_g, nmemsize));
    cuda_safe_call(cudaMemcpy(mlen0_g, mlen0.data(), nmemsize, cudaMemcpyHostToDevice));
}

template<typename TF>
void Diff_tke2<TF>::clear_device()
{
    cuda_safe_call(cudaFree(mlen0_g));
}
#endif

#ifdef FLOAT_SINGLE
template class Diff_tke2<float>;
#else
template class Diff_tke2<double>;
#endif
