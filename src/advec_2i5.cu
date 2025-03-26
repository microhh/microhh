/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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

#include "advec_2i5.h"
#include "advec_2i5_kernels.cuh"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "tools.h"
#include "constants.h"
#include "finite_difference.h"
#include "field3d_operators.h"
#include "cuda_launcher.h"


#ifdef USECUDA
template<typename TF>
unsigned long Advec_2i5<TF>::get_time_limit(unsigned long idt, double dt)
{
    // Calculate cfl and prevent zero divisons.
    double cfl = get_cfl(dt);
    cfl = std::max(cflmin, cfl);
    const unsigned long idtlim = idt * cflmax / cfl;

    return idtlim;
}


template<typename TF>
double Advec_2i5<TF>::get_cfl(const double dt)
{
    const Grid_data<TF>& gd = grid.get_grid_data();
    auto tmp1 = fields.get_tmp_g();

    Grid_layout grid_layout = {
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.istride,
            gd.jstride,
            gd.kstride};

    launch_grid_kernel<Advec_2i5_kernels::calc_cfl_g<TF>>(
            grid_layout,
            tmp1->fld_g.view(),
            fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
            gd.dzi_g, gd.dxi, gd.dyi);
    cuda_check_error();

    TF cfl = field3d_operators.calc_max_g(tmp1->fld_g);
    fields.release_tmp_g(tmp1);

    master.max(&cfl, 1);

    cfl = cfl*dt;

    return static_cast<double>(cfl);
}


template<typename TF>
void Advec_2i5<TF>::exec(Stats<TF>& stats)
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    Grid_layout grid_layout = {
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kend,
            gd.istride,
            gd.jstride,
            gd.kstride};

    launch_grid_kernel<Advec_2i5_kernels::advec_u_g<TF>>(
            grid_layout,
            fields.mt.at("u")->fld_g.view(),
            fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
            fields.rhorefi_g, fields.rhorefh_g, gd.dzi_g, gd.dxi, gd.dyi);

    launch_grid_kernel<Advec_2i5_kernels::advec_v_g<TF>>(
            grid_layout,
            fields.mt.at("v")->fld_g.view(),
            fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
            fields.rhorefi_g, fields.rhorefh_g, gd.dzi_g, gd.dxi, gd.dyi);

    launch_grid_kernel<Advec_2i5_kernels::advec_w_g<TF>>(
            grid_layout,
            fields.mt.at("w")->fld_g.view(),
            fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
            fields.rhoref_g, fields.rhorefhi_g, gd.dzhi_g, gd.dxi, gd.dyi);

    for (const std::string& s : sp_limit)
    {
        launch_grid_kernel<Advec_2i5_kernels::advec_s_lim_g<TF>>(
                grid_layout,
                fields.st.at(s)->fld_g.view(), fields.sp.at(s)->fld_g,
                fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
                fields.rhorefi_g, fields.rhorefh_g, gd.dzi_g, gd.dxi, gd.dyi);
    }

    for (const std::string& s : sp_no_limit)
    {
        launch_grid_kernel<Advec_2i5_kernels::advec_s_g<TF>>(
                grid_layout,
                fields.st.at(s)->fld_g.view(), fields.sp.at(s)->fld_g,
                fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
                fields.rhorefi_g, fields.rhorefh_g, gd.dzi_g, gd.dxi, gd.dyi);
    }

    cudaDeviceSynchronize();
    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);
    for (auto it : fields.st)
        stats.calc_tend(*it.second, tend_name);
}

template<typename TF>
void Advec_2i5<TF>::get_advec_flux(
        Field3d<TF>& advec_flux, const Field3d<TF>& fld)
{
    using namespace Tools_g;

    auto& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);
    const int ngrid  = gd.ncells/blocki + (gd.ncells%blocki > 0);
    const int nblock = gd.ithread_block;
    dim3 gridGPU(gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);
//     if (fld.loc == gd.uloc)
//     {
//         advec_flux_u<TF><<<gridGPU, blockGPU>>>(
//                 advec_flux.fld_g.data(), fld.fld_g.data(), fields.mp.at("w")->fld_g.data(), gd.icells, gd.jcells, gd.kstart, gd.kend, gd.ijcells);
//     }
//     else if (fld.loc == gd.vloc)
//     {
//         advec_flux_v<TF><<<gridGPU, blockGPU>>>(
//                 advec_flux.fld_g.data(), fld.fld_g.data(), fields.mp.at("w")->fld_g.data(), gd.icells, gd.jcells, gd.kstart, gd.kend, gd.ijcells);
//     }
//     else if (fld.loc == gd.wloc)
//     {
//         advec_flux_w<TF><<<gridGPU, blockGPU>>>(
//                 advec_flux.fld_g.data(), fld.fld_g.data(), fields.mp.at("w")->fld_g.data(), gd.icells, gd.jcells, gd.kstart, gd.kend, gd.ijcells);
//     }
//     else if (fld.loc == gd.sloc)
//     {
//         advec_flux_s_g<TF><<<gridGPU, blockGPU>>>(
//                 advec_flux.fld_g.data(), fld.fld_g.data(), fields.mp.at("w")->fld_g.data(), gd.icells, gd.jcells, gd.kstart, gd.kend, gd.ijcells);
//     }
//     else
        throw std::runtime_error("Advec_2i5 cannot deliver flux field at that location");
}

#endif


#ifdef FLOAT_SINGLE
template class Advec_2i5<float>;
#else
template class Advec_2i5<double>;
#endif
