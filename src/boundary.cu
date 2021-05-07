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
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "timeloop.h"
#include "timedep.h"
#include "defines.h"
#include "finite_difference.h"
#include "constants.h"
#include "tools.h"
#include "boundary.h"

using namespace Finite_difference::O4;

namespace
{
    template<typename TF> __global__
    void set_bc_value_g(TF* __restrict__ a, TF aval,
                  const int icells, const int jcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*icells;
            a[ij] = aval;
        }
    }


    template<typename TF> __global__
    void calc_ghost_cells_bot_2nd_g(TF* __restrict__ a, TF* __restrict__ dzh, Boundary_type sw,
                                    TF* __restrict__ abot, TF* __restrict__ agradbot,
                                    const int icells, const int jcells, const int kstart)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int kk  = icells*jcells;
        const int ij  = i + j*icells;
        const int ijk = i + j*icells + kstart*kk;

        if (i < icells && j < jcells)
        {
            if (sw == Boundary_type::Dirichlet_type)
                a[ijk-kk] = TF(2.)*abot[ij] - a[ijk];

            else if (sw == Boundary_type::Neumann_type || sw == Boundary_type::Flux_type)
                a[ijk-kk] = -agradbot[ij]*dzh[kstart] + a[ijk];
        }
    }

    template<typename TF> __global__
    void calc_ghost_cells_top_2nd_g(TF* __restrict__ a, TF* __restrict__ dzh, const Boundary_type sw,
                                    TF* __restrict__ atop, TF* __restrict__ agradtop,
                                    const int icells, const int jcells, const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int kk  = icells*jcells;
        const int ij  = i + j*icells;
        const int ijk = i + j*icells + (kend-1)*kk;

        if (i < icells && j < jcells)
        {
            if (sw == Boundary_type::Dirichlet_type)
                a[ijk+kk] = TF(2.)*atop[ij] - a[ijk];

            else if (sw == Boundary_type::Off_type)
            {
                atop[ij] = TF(3./2.)*a[ijk] - TF(1./2.)*a[ijk-kk];
                a[ijk+kk] = TF(2.)*atop[ij] - a[ijk];
            }

            else if (sw == Boundary_type::Neumann_type || sw == Boundary_type::Flux_type)
                a[ijk+kk] = agradtop[ij]*dzh[kend] + a[ijk];
        }
    }

    template<typename TF> __global__
    void calc_ghost_cells_bot_4th_g(TF* __restrict__ a, TF* __restrict__ z, const Boundary_type sw,
                                    TF* __restrict__ abot, TF* __restrict__ agradbot,
                                    const int icells, const int jcells, const int kstart)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int kk1 = 1*icells*jcells;
        const int kk2 = 2*icells*jcells;

        const int ij  = i + j*icells;
        const int ijk = i + j*icells + kstart*kk1;

        if (i < icells && j < jcells)
        {
            if (sw == Boundary_type::Dirichlet_type)
            {
                a[ijk-kk1] = TF(8./3.)*abot[ij] - TF(2.)*a[ijk] + TF(1./3.)*a[ijk+kk1];
                a[ijk-kk2] = TF(8.)*abot[ij] - TF(9.)*a[ijk] + TF(2.)*a[ijk+kk1];
            }

            else if (sw == Boundary_type::Neumann_type || sw == Boundary_type::Flux_type)
            {
                a[ijk-kk1] = TF(-1.)*grad4(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk    ];
                a[ijk-kk2] = TF(-3.)*grad4(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk+kk1];
            }
        }
    }

    template<typename TF> __global__
    void calc_ghost_cells_top_4th_g(TF* __restrict__ a, TF* __restrict__ z,const Boundary_type sw,
                                    TF* __restrict__ atop, TF* __restrict__ agradtop,
                                    const int icells, const int jcells, const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int kk1 = 1*icells*jcells;
        const int kk2 = 2*icells*jcells;

        const int ij  = i + j*icells;
        const int ijk = i + j*icells + (kend-1)*kk1;

        if( i < icells && j < jcells)
        {
            if (sw == Boundary_type::Dirichlet_type)
            {
                a[ijk+kk1] = TF(8./3.)*atop[ij] - TF(2.)*a[ijk] + TF(1./3.)*a[ijk-kk1];
                a[ijk+kk2] = TF(8.)*atop[ij] - TF(9.)*a[ijk] + TF(2.)*a[ijk-kk1];
            }

            else if (sw == Boundary_type::Neumann_type || sw == Boundary_type::Flux_type)
            {
                a[ijk+kk1] = TF(1.)*grad4(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk    ];
                a[ijk+kk2] = TF(3.)*grad4(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk-kk1];
            }
        }
    }

    template<typename TF> __global__
    void calc_ghost_cells_botw_4th_g(TF* __restrict__ w,
                                      const int icells, const int jcells, const int kstart)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int kk1 = 1*icells*jcells;
        const int kk2 = 2*icells*jcells;
        const int kk3 = 3*icells*jcells;

        const int ijk = i + j*icells + kstart*kk1;

        if (i < icells && j < jcells)
        {
            w[ijk-kk1] = TF(-6.)*w[ijk+kk1] + TF(4.)*w[ijk+kk2] - w[ijk+kk3];
        }
    }

    template<typename TF> __global__
    void calc_ghost_cells_topw_4th_g(TF* __restrict__ w,
                                      const int icells, const int jcells, const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int kk1 = 1*icells*jcells;
        const int kk2 = 2*icells*jcells;
        const int kk3 = 3*icells*jcells;

        const int ijk = i + j*icells + kend*kk1;

        if (i < icells && j < jcells)
        {
            w[ijk+kk1] = TF(-6.)*w[ijk-kk1] + TF(4.)*w[ijk-kk2] - w[ijk-kk3];
        }
    }

    template<typename TF> __global__
    void calc_ghost_cells_botw_cons_4th_g(TF* __restrict__ w,
                                      const int icells, const int jcells, const int kstart)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int kk1 = 1*icells*jcells;
        const int kk2 = 2*icells*jcells;

        const int ijk = i + j*icells + kstart*kk1;

        if (i < icells && j < jcells)
        {
            w[ijk-kk1] = -w[ijk+kk1];
            w[ijk-kk2] = -w[ijk+kk2];
        }
    }

    template<typename TF> __global__
    void calc_ghost_cells_topw_cons_4th_g(TF* __restrict__ w,
                                      const int icells, const int jcells, const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int kk1 = 1*icells*jcells;
        const int kk2 = 2*icells*jcells;

        const int ijk = i + j*icells + kend*kk1;

        if (i < icells && j < jcells)
        {
            w[ijk+kk1] = -w[ijk-kk1];
            w[ijk+kk2] = -w[ijk-kk2];
        }
    }
}

#ifdef USECUDA
template<typename TF>
void Boundary<TF>::set_ghost_cells()
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 grid2dGPU (gridi, gridj);
    dim3 block2dGPU(blocki, blockj);

    if (grid.get_spatial_order() == Grid_order::Second)
    {
        calc_ghost_cells_bot_2nd_g<TF><<<grid2dGPU, block2dGPU>>>(
            fields.mp.at("u")->fld_g, gd.dzh_g, mbcbot,
            fields.mp.at("u")->fld_bot_g, fields.mp.at("u")->grad_bot_g,
            gd.icells, gd.jcells, gd.kstart);
        cuda_check_error();

        calc_ghost_cells_top_2nd_g<TF><<<grid2dGPU, block2dGPU>>>(
            fields.mp.at("u")->fld_g, gd.dzh_g, mbctop,
            fields.mp.at("u")->fld_top_g, fields.mp.at("u")->grad_top_g,
            gd.icells, gd.jcells, gd.kend);
        cuda_check_error();

        calc_ghost_cells_bot_2nd_g<TF><<<grid2dGPU, block2dGPU>>>(
            fields.mp.at("v")->fld_g, gd.dzh_g, mbcbot,
            fields.mp.at("v")->fld_bot_g, fields.mp.at("v")->grad_bot_g,
            gd.icells, gd.jcells, gd.kstart);
        cuda_check_error();

        calc_ghost_cells_top_2nd_g<TF><<<grid2dGPU, block2dGPU>>>(
            fields.mp.at("v")->fld_g, gd.dzh_g, mbctop,
            fields.mp.at("v")->fld_top_g, fields.mp.at("v")->grad_top_g,
            gd.icells, gd.jcells, gd.kend);
        cuda_check_error();

        for (auto& it : fields.sp)
        {
            calc_ghost_cells_bot_2nd_g<TF><<<grid2dGPU, block2dGPU>>>(
                it.second->fld_g, gd.dzh_g, sbc.at(it.first).bcbot,
                it.second->fld_bot_g, it.second->grad_bot_g,
                gd.icells, gd.jcells, gd.kstart);
            cuda_check_error();

            calc_ghost_cells_top_2nd_g<TF><<<grid2dGPU, block2dGPU>>>(
                it.second->fld_g, gd.dzh_g, sbc.at(it.first).bctop,
                it.second->fld_top_g, it.second->grad_top_g,
                gd.icells, gd.jcells, gd.kend);
            cuda_check_error();
        }
    }
    else if (grid.get_spatial_order() == Grid_order::Fourth)
    {
        calc_ghost_cells_bot_4th_g<TF><<<grid2dGPU, block2dGPU>>>(
            fields.mp.at("u")->fld_g, gd.z_g, mbcbot,
            fields.mp.at("u")->fld_bot_g, fields.mp.at("u")->grad_bot_g,
            gd.icells, gd.jcells, gd.kstart);
        cuda_check_error();

        calc_ghost_cells_top_4th_g<TF><<<grid2dGPU, block2dGPU>>>(
            fields.mp.at("u")->fld_g, gd.z_g, mbctop,
            fields.mp.at("u")->fld_top_g, fields.mp.at("u")->grad_top_g,
            gd.icells, gd.jcells, gd.kend);
        cuda_check_error();

        calc_ghost_cells_bot_4th_g<TF><<<grid2dGPU, block2dGPU>>>(
            fields.mp.at("v")->fld_g, gd.z_g, mbcbot,
            fields.mp.at("v")->fld_bot_g, fields.mp.at("v")->grad_bot_g,
            gd.icells, gd.jcells, gd.kstart);
        cuda_check_error();

        calc_ghost_cells_top_4th_g<TF><<<grid2dGPU, block2dGPU>>>(
            fields.mp.at("v")->fld_g, gd.z_g, mbctop,
            fields.mp.at("v")->fld_top_g, fields.mp.at("v")->grad_top_g,
            gd.icells, gd.jcells, gd.kend);
        cuda_check_error();

        for (auto& it : fields.sp)
        {
            calc_ghost_cells_bot_4th_g<TF><<<grid2dGPU, block2dGPU>>>(
                it.second->fld_g, gd.z_g, sbc.at(it.first).bcbot,
                it.second->fld_bot_g, it.second->grad_bot_g,
                gd.icells, gd.jcells, gd.kstart);
            cuda_check_error();

            calc_ghost_cells_top_4th_g<TF><<<grid2dGPU, block2dGPU>>>(
                it.second->fld_g, gd.z_g, sbc.at(it.first).bctop,
                it.second->fld_top_g, it.second->grad_top_g,
                gd.icells, gd.jcells, gd.kend);
            cuda_check_error();
        }
    }
}

template<typename TF>
void Boundary<TF>::set_ghost_cells_w(const Boundary_w_type boundary_w_type)
{
    auto& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 grid2dGPU (gridi,  gridj );
    dim3 block2dGPU(blocki, blockj);

    if (grid.get_spatial_order() == Grid_order::Fourth)
    {
        if (boundary_w_type == Boundary_w_type::Normal_type)
        {
            calc_ghost_cells_botw_4th_g<TF><<<grid2dGPU, block2dGPU>>>(
                fields.mp.at("w")->fld_g,
                gd.icells, gd.jcells, gd.kstart);
            cuda_check_error();

            calc_ghost_cells_topw_4th_g<TF><<<grid2dGPU, block2dGPU>>>(
                fields.mp.at("w")->fld_g,
                gd.icells, gd.jcells, gd.kend);
            cuda_check_error();
        }
        else if (boundary_w_type == Boundary_w_type::Conservation_type)
        {
            calc_ghost_cells_botw_cons_4th_g<TF><<<grid2dGPU, block2dGPU>>>(
                fields.mp.at("w")->fld_g,
                gd.icells, gd.jcells, gd.kstart);
            cuda_check_error();

            calc_ghost_cells_topw_cons_4th_g<TF><<<grid2dGPU, block2dGPU>>>(
                fields.mp.at("w")->fld_g,
                gd.icells, gd.jcells, gd.kend);
            cuda_check_error();
        }
    }
}
#endif

template<typename TF>
void Boundary<TF>::clear_device()
{
    for(auto& it : tdep_bc)
        it.second->clear_device();
}

template<typename TF>
void Boundary<TF>::set_bc_g(TF* restrict a, TF* restrict agrad, TF* restrict aflux,
                        Boundary_type sw, TF aval, TF visc, TF offset)
{
    auto& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 grid2dGPU (gridi, gridj);
    dim3 block2dGPU(blocki, blockj);

    if (sw == Boundary_type::Dirichlet_type)
    {
        set_bc_value_g<TF><<<grid2dGPU, block2dGPU>>>(a, aval-offset, gd.icells, gd.jcells);
        cuda_check_error();
    }
    else if (sw == Boundary_type::Neumann_type)
    {
        set_bc_value_g<TF><<<grid2dGPU, block2dGPU>>>(agrad, aval, gd.icells, gd.jcells);
        set_bc_value_g<TF><<<grid2dGPU, block2dGPU>>>(aflux, -aval*visc, gd.icells, gd.jcells);
        cuda_check_error();
    }
    else if (sw == Boundary_type::Flux_type)
    {
        set_bc_value_g<TF><<<grid2dGPU, block2dGPU>>>(aflux, aval, gd.icells, gd.jcells);
        set_bc_value_g<TF><<<grid2dGPU, block2dGPU>>>(agrad, -aval*visc, gd.icells, gd.jcells);
        cuda_check_error();
    }
}

#ifdef USECUDA
template <typename TF>
void Boundary<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    const TF no_offset = 0.;

    for (auto& it : tdep_bc)
    {
        it.second->update_time_dependent(sbc.at(it.first).bot,timeloop);
        set_bc_g(
                fields.sp.at(it.first)->fld_bot_g,
                fields.sp.at(it.first)->grad_bot_g,
                fields.sp.at(it.first)->flux_bot_g,
                sbc.at(it.first).bcbot,
                sbc.at(it.first).bot,
                fields.sp.at(it.first)->visc,
                no_offset);
    }
}
#endif

#ifdef USECUDA
template<typename TF>
void Boundary<TF>::set_prognostic_cyclic_bcs()
{
    /* Set cyclic boundary conditions of the
       prognostic 3D fields */
    boundary_cyclic.exec_g(fields.mp.at("u")->fld_g);
    boundary_cyclic.exec_g(fields.mp.at("v")->fld_g);
    boundary_cyclic.exec_g(fields.mp.at("w")->fld_g);

    for (auto& it : fields.sp)
        boundary_cyclic.exec_g(it.second->fld_g);

    // Overwrite here the ghost cells for the scalars with outflow BCs
    for (auto& s : scalar_outflow)
        boundary_outflow.exec(fields.sp.at(s)->fld_g);
}
#endif

#ifdef USECUDA
template<typename TF>
TF* Boundary<TF>::get_z0m_g()
{
    throw std::runtime_error("Function get_z0m_g() not implemented in base boundary.");
}

template<typename TF>
TF* Boundary<TF>::get_ustar_g()
{
    throw std::runtime_error("Function get_z0m_g() not implemented in base boundary.");
}

template<typename TF>
TF* Boundary<TF>::get_obuk_g()
{
    throw std::runtime_error("Function get_z0m_g() not implemented in base boundary.");
}
#endif

template class Boundary<double>;
template class Boundary<float>;
