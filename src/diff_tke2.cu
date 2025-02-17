/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
 * Copyright (c) 2021-2024 Steven van der Linden
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

// Kernel Launcher
#include "cuda_launcher.h"
#include "diff_kl_kernels.cuh"
#include "diff_tke2_kl_kernels.cuh"
#include "cuda_buffer.h"

namespace
{
    namespace most = Monin_obukhov;
    namespace fm = Fast_math;
    namespace dk = Diff_kernels_g;

    template<typename TF, Surface_model surface_model, bool sw_mason> __global__
    void calc_evisc_neutral_g(
            TF* const restrict evisc,
            const TF* const restrict sgstke,
            const TF* const restrict z,
            const TF* const restrict z0m,
            const TF* const restrict mlen0,
            const TF cm,
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
                if constexpr (sw_mason) // Apply Mason's wall correction here
                {
                    if constexpr (n_mason == 2)
                        fac = std::sqrt(TF(1.) / ( TF(1.)/fm::pow2(mlen0[k]) + TF(1.)/(fm::pow2(Constants::kappa<TF>*(z[k]+z0m[ij]))) ) );
                    else
                        fac = std::pow(TF(1.) / (TF(1.)/std::pow(mlen0[k], TF(n_mason)) + TF(1.)/
                                    (std::pow(Constants::kappa<TF>*(z[k]+z0m[ij]), TF(n_mason)))), TF(1.)/TF(n_mason));
                }
                else
                    fac = mlen0[k];

                // Calculate eddy diffusivity for momentum.
                evisc[ijk] = cm * fac * sqrt(sgstke[ijk]);
            }
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

            if constexpr (sw_mason) // Apply Mason's wall correction here
            {
                if constexpr (n_mason == 2)
                    fac = sqrt(TF(1.) / ( TF(1.)/fm::pow2(mlen0[k]) + TF(1.)/(fm::pow2(Constants::kappa<TF>*(z[k]+z0m[ij]))) ) );
                else
                    fac = pow(TF(1.) / (TF(1.)/std::pow(mlen0[k], TF(n_mason)) + TF(1.)/
                                (pow(Constants::kappa<TF>*(z[k]+z0m[ij]), TF(n_mason)))), TF(1.)/TF(n_mason));
            }

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

    // Grid layout struct for cuda launcher.
    Grid_layout grid_layout = {
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.istride,
            gd.jstride,
            gd.kstride};

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

    launch_grid_kernel<Diff_les_kernels::diff_uvw_g<TF, true>>(
            grid_layout,
            fields.mt.at("u")->fld_g.view(),
            fields.mt.at("v")->fld_g.view(),
            fields.mt.at("w")->fld_g.view(),
            fields.sd.at("evisc")->fld_g,
            fields.mp.at("u")->fld_g,
            fields.mp.at("v")->fld_g,
            fields.mp.at("w")->fld_g,
            fields.mp.at("u")->flux_bot_g,
            fields.mp.at("u")->flux_top_g,
            fields.mp.at("v")->flux_bot_g,
            fields.mp.at("v")->flux_top_g,
            gd.dzi_g, gd.dzhi_g, gd.dxi, gd.dyi,
            fields.rhoref_g, fields.rhorefh_g,
            fields.rhorefi_g, fields.rhorefhi_g,
            fields.visc);

    cuda_check_error();

    for (auto it : fields.st)
    {
        cuda_vector<TF>* evisc_ptr;

        if (it.first == "sgstke")  // sgstke diffuses with eddy viscosity for momentum
            evisc_ptr = &fields.sd.at("evisc")->fld_g;
        else  // all other scalars, normally diffuse with eddy viscosity for heat/scalars
        {
            if (!sw_buoy) // but not if there is no buoyancy (then eviscs not defined)
                evisc_ptr = &fields.sd.at("evisc")->fld_g;
            else
                evisc_ptr = &fields.sd.at("eviscs")->fld_g;
        }

        launch_grid_kernel<Diff_les_kernels::diff_c_g<TF, true>>(
                grid_layout,
                it.second->fld_g.view(),
                fields.sp.at(it.first)->fld_g,
                *evisc_ptr,
                fields.sp.at(it.first)->flux_bot_g,
                fields.sp.at(it.first)->flux_top_g,
                gd.dzi_g, gd.dzhi_g,
                dxidxi, dyidyi,
                fields.rhorefi_g,
                fields.rhorefh_g,
                tPr_i_dummy,
                fields.sp.at(it.first)->visc);
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

    // Grid layout struct for cuda launcher.
    Grid_layout grid_layout = {
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.istride,
            gd.jstride,
            gd.kstride};

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
    auto& dudz_g  = boundary.get_dudz_g();
    auto& dvdz_g  = boundary.get_dvdz_g();
    auto& z0m_g   = boundary.get_z0m_g();

    auto str2_tmp = fields.get_tmp_g();

    // Calculate total strain rate
    launch_grid_kernel<Diff_les_kernels::calc_strain2_g<TF, true>>(
            grid_layout,
            str2_tmp->fld_g.view(),
            fields.mp.at("u")->fld_g,
            fields.mp.at("v")->fld_g,
            fields.mp.at("w")->fld_g,
            dudz_g, dvdz_g,
            gd.dzi_g, gd.dzhi_g,
            gd.dxi, gd.dyi);

    // Start with retrieving the stability information
    if (!sw_buoy)
    {
        if (sw_mason)
            calc_evisc_neutral_g<TF, Surface_model::Enabled, true><<<gridGPU, blockGPU>>>(
                    fields.sd.at("evisc")->fld_g,
                    fields.sp.at("sgstke")->fld_g,
                    gd.z_g,
                    z0m_g, mlen0_g,
                    this->cm,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells,
                    gd.ijcells);
        else
            calc_evisc_neutral_g<TF, Surface_model::Enabled, false><<<gridGPU, blockGPU>>>(
                    fields.sd.at("evisc")->fld_g,
                    fields.sp.at("sgstke")->fld_g,
                    gd.z_g,
                    z0m_g, mlen0_g,
                    this->cm,
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
        auto& dbdz_g = boundary.get_dbdz_g();

        // Note BvS: templated lambda functions are not (yet?) allowed by NVCC :-(
        if (sw_mason)
            launch_grid_kernel<Diff_tke2_kernels::evisc_g<TF, true, true>>(
                    grid_layout,
                    fields.sd.at("evisc")->fld_g.view(),
                    fields.sp.at("sgstke")->fld_g,
                    buoy_tmp->fld_g,
                    dbdz_g,
                    gd.z_g,
                    z0m_g,
                    mlen0_g,
                    this->cn,
                    this->cm);
        else
            launch_grid_kernel<Diff_tke2_kernels::evisc_g<TF, true, false>>(
                    grid_layout,
                    fields.sd.at("evisc")->fld_g.view(),
                    fields.sp.at("sgstke")->fld_g,
                    buoy_tmp->fld_g,
                    dbdz_g,
                    gd.z_g,
                    z0m_g,
                    mlen0_g,
                    this->cn,
                    this->cm);

        if (sw_mason)
            launch_grid_kernel<Diff_tke2_kernels::evisc_heat_g<TF, true, true>>(
                    grid_layout,
                    fields.sd.at("eviscs")->fld_g.view(),
                    fields.sd.at("evisc")->fld_g,
                    fields.sp.at("sgstke")->fld_g,
                    buoy_tmp->fld_g, dbdz_g,
                    gd.z_g,
                    z0m_g, mlen0_g,
                    this->cn, this->ch1, this->ch2);
        else
            launch_grid_kernel<Diff_tke2_kernels::evisc_heat_g<TF, true, false>>(
                    grid_layout,
                    fields.sd.at("eviscs")->fld_g.view(),
                    fields.sd.at("evisc")->fld_g,
                    fields.sp.at("sgstke")->fld_g,
                    buoy_tmp->fld_g, dbdz_g,
                    gd.z_g,
                    z0m_g, mlen0_g,
                    this->cn, this->ch1, this->ch2);

        boundary_cyclic.exec_g(fields.sd.at("evisc")->fld_g);
        boundary_cyclic.exec_g(fields.sd.at("eviscs")->fld_g);

        // Calculate tendencies here, to prevent having to
        // re-calculate `strain2` in diffusion->exec()`.
        launch_grid_kernel<Diff_tke2_kernels::sgstke_buoy_tend_g<TF>>(
                grid_layout,
                fields.st.at("sgstke")->fld_g.view(),
                fields.sp.at("sgstke")->fld_g,
                fields.sd.at("eviscs")->fld_g,
                buoy_tmp->fld_g,
                dbdz_g);

        cudaDeviceSynchronize();
        stats.calc_tend(*fields.st.at("sgstke"), tend_name_buoy);

        if (sw_mason)
            launch_grid_kernel<Diff_tke2_kernels::sgstke_diss_tend_g<TF, true>>(
                grid_layout,
                fields.st.at("sgstke")->fld_g.view(),
                fields.sp.at("sgstke")->fld_g,
                buoy_tmp->fld_g,
                dbdz_g,
                gd.z_g,
                z0m_g,
                mlen0_g,
                this->cn,
                this->ce1,
                this->ce2);
        else
            launch_grid_kernel<Diff_tke2_kernels::sgstke_diss_tend_g<TF, false>>(
                grid_layout,
                fields.st.at("sgstke")->fld_g.view(),
                fields.sp.at("sgstke")->fld_g,
                buoy_tmp->fld_g,
                dbdz_g,
                gd.z_g,
                z0m_g,
                mlen0_g,
                this->cn,
                this->ce1,
                this->ce2);

        cudaDeviceSynchronize();
        stats.calc_tend(*fields.st.at("sgstke"), tend_name_diss);

        fields.release_tmp_g(buoy_tmp);
    }

    launch_grid_kernel<Diff_tke2_kernels::sgstke_shear_tend_g<TF>>(
            grid_layout,
            fields.st.at("sgstke")->fld_g.view(),
            fields.sp.at("sgstke")->fld_g,
            fields.sd.at("evisc")->fld_g,
            str2_tmp->fld_g);

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

    mlen0_g.allocate(gd.kcells);
    cuda_safe_call(cudaMemcpy(mlen0_g, mlen0.data(), mlen0_g.size_in_bytes(), cudaMemcpyHostToDevice));
}

template<typename TF>
void Diff_tke2<TF>::clear_device()
{
}
#endif

#ifdef FLOAT_SINGLE
template class Diff_tke2<float>;
#else
template class Diff_tke2<double>;
#endif
