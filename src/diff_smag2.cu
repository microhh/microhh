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
#include "cuda_launcher.h"
#include "diff_smag2_kernels.cuh"

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

    mlen_g.allocate(gd.kcells);
    cuda_safe_call(cudaMemcpy(mlen_g, mlen.data(), mlen_g.size_in_bytes(), cudaMemcpyHostToDevice));
}

template<typename TF>
void Diff_smag2<TF>::clear_device()
{
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
        auto& z0m_g   = boundary.get_z0m_g();

        // Get MO gradients velocity:
        auto& dudz_g  = boundary.get_dudz_g();
        auto& dvdz_g  = boundary.get_dvdz_g();

        // Calculate total strain rate
        launch_grid_kernel<diff_smag2::calc_strain2_g<TF, true>>(
            gd,
            fields.sd.at("evisc")->fld_g.view(),
            fields.mp.at("u")->fld_g,
            fields.mp.at("v")->fld_g,
            fields.mp.at("w")->fld_g,
            dudz_g, dvdz_g,
            gd.dzi_g, gd.dzhi_g,
            gd.dxi, gd.dyi);

        if (thermo.get_switch() == Thermo_type::Disabled)
        {
            // Start with retrieving the stability information
            diff_smag2::evisc_neutral_g<TF><<<gridGPU, blockGPU>>>(
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
            auto& dbdz_g  = boundary.get_dbdz_g();

            // Calculate eddy viscosity
            TF tPri = 1./tPr;

            launch_grid_kernel<diff_smag2::evisc_g<TF, true>>(
                gd,
                fields.sd.at("evisc")->fld_g.view(),
                tmp1->fld_g, dbdz_g,
                mlen_g, z0m_g, gd.z_g,
                tPri);

            fields.release_tmp_g(tmp1);
        }

        boundary_cyclic.exec_g(fields.sd.at("evisc")->fld_g);
    }
    // Do not use surface model.
    else
    {
        // Calculate total strain rate
        launch_grid_kernel<diff_smag2::calc_strain2_g<TF, false>>(
            gd,
            fields.sd.at("evisc")->fld_g.view(),
            fields.mp.at("u")->fld_g,
            fields.mp.at("v")->fld_g,
            fields.mp.at("w")->fld_g,
            nullptr, nullptr,
            gd.dzi_g, gd.dzhi_g,
            gd.dxi, gd.dyi);

        // start with retrieving the stability information
        if (thermo.get_switch() == Thermo_type::Disabled)
        {
            diff_smag2::evisc_neutral_vandriest_g<TF><<<gridGPU, blockGPU>>>(
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

            launch_grid_kernel<diff_smag2::evisc_g<TF, true>>(
                gd,
                fields.sd.at("evisc")->fld_g.view(),
                tmp1->fld_g, nullptr,
                mlen_g, nullptr, gd.z_g,
                tPri);

            fields.release_tmp_g(tmp1);
        }

        boundary_cyclic.exec_g(fields.sd.at("evisc")->fld_g);
        diff_smag2::calc_ghostcells_evisc<TF><<<grid2dGPU, block2dGPU>>>(
                fields.sd.at("evisc")->fld_g,
                gd.icells, gd.jcells,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }

    cuda_check_error();
}
#endif

#ifdef USECUDA
template<typename TF>
void Diff_smag2<TF>::exec(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    const TF dxidxi = TF(1)/(gd.dx * gd.dx);
    const TF dyidyi = TF(1)/(gd.dy * gd.dy);
    const TF tPri = TF(1)/tPr;

    // Do not use surface model.
    if (boundary.get_switch() == "default")
    {
        launch_grid_kernel<diff_smag2::diff_uvw_g<TF, false>>(
                gd,
                fields.mt.at("u")->fld_g.view(), fields.mt.at("v")->fld_g.view(), fields.mt.at("w")->fld_g.view(),
                fields.sd.at("evisc")->fld_g,
                fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
                fields.mp.at("u")->flux_bot_g, fields.mp.at("u")->flux_top_g,
                fields.mp.at("v")->flux_bot_g, fields.mp.at("v")->flux_top_g,
                gd.dzi_g, gd.dzhi_g, gd.dxi, gd.dyi,
                fields.rhoref_g, fields.rhorefh_g,
                fields.rhorefi_g, fields.rhorefhi_g,
                fields.visc);

        for (auto it : fields.st)
        {
            launch_grid_kernel<diff_smag2::diff_c_g<TF, false>>(
                    gd,
                    it.second->fld_g.view(), fields.sp.at(it.first)->fld_g, fields.sd.at("evisc")->fld_g,
                    fields.sp.at(it.first)->flux_bot_g, fields.sp.at(it.first)->flux_top_g,
                    gd.dzi_g, gd.dzhi_g, dxidxi, dyidyi,
                    fields.rhorefi_g, fields.rhorefh_g,
                    tPri, fields.sp.at(it.first)->visc);
        }
        cuda_check_error();
    }
    // Use surface model.
    else
    {
        launch_grid_kernel<diff_smag2::diff_uvw_g<TF, true>>(
                gd,
                fields.mt.at("u")->fld_g.view(), fields.mt.at("v")->fld_g.view(), fields.mt.at("w")->fld_g.view(),
                fields.sd.at("evisc")->fld_g,
                fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
                fields.mp.at("u")->flux_bot_g, fields.mp.at("u")->flux_top_g,
                fields.mp.at("v")->flux_bot_g, fields.mp.at("v")->flux_top_g,
                gd.dzi_g, gd.dzhi_g, gd.dxi, gd.dyi,
                fields.rhoref_g, fields.rhorefh_g,
                fields.rhorefi_g, fields.rhorefhi_g,
                fields.visc);

        for (auto it : fields.st)
            launch_grid_kernel<diff_smag2::diff_c_g<TF, true>>(
                    gd,
                    it.second->fld_g.view(), fields.sp.at(it.first)->fld_g, fields.sd.at("evisc")->fld_g,
                    fields.sp.at(it.first)->flux_bot_g, fields.sp.at(it.first)->flux_top_g,
                    gd.dzi_g, gd.dzhi_g, dxidxi, dyidyi,
                    fields.rhorefi_g, fields.rhorefh_g,
                    tPri, fields.sp.at(it.first)->visc);
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
    diff_smag2::calc_dnmul_g<TF><<<gridGPU, blockGPU>>>(
            tmp1->fld_g, fields.sd.at("evisc")->fld_g,
            gd.dzi_g, tPrfac_i, dxidxi, dyidyi,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend,
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

    diff_smag2::calc_dnmul_g<TF><<<gridGPU, blockGPU>>>(
        dnmul_tmp->fld_g, fields.sd.at("evisc")->fld_g,
        gd.dzi_g, tPrfac_i, dxidxi, dyidyi,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // Get maximum from tmp1 field
    // CvH This is odd, because there might be need for calc_max in CPU version.
    double dnmul = field3d_operators.calc_max_g(dnmul_tmp->fld_g);

    fields.release_tmp_g(dnmul_tmp);

    return dnmul*dt;
}
#endif


#ifdef FLOAT_SINGLE
template class Diff_smag2<float>;
#else
template class Diff_smag2<double>;
#endif
