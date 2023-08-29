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
    void set_bc_value_g(
            TF* const __restrict__ a,
            const TF aval,
            const int icells,
            const int jcells)
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
    void set_bc_2d_value_g(
            TF* const __restrict__ a,
            const TF* const __restrict__ aval,
            const TF offset,
            const TF factor,
            const int icells,
            const int jcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*icells;
            a[ij] = (aval[ij]-offset)*factor;
        }
    }

    template<typename TF> __global__
    void copy_xy_g(
            TF* const __restrict__ out,
            const TF* const __restrict__ in,
            const int icells,
            const int jcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*icells;
            out[ij] = in[ij];
        }
    }

    template<typename TF> __global__
    void interp_sbot_time_g(
            TF* const __restrict__ fld_out,
            const TF* const __restrict__ fld_prev,
            const TF* const __restrict__ fld_next,
            const TF fac0,
            const TF fac1,
            const int icells,
            const int jcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*icells;
            fld_out[ij] = fac0*fld_prev[ij] + fac1*fld_next[ij];
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
void Boundary<TF>::set_bc_g(
        TF* const __restrict__ a,
        TF* const __restrict__ agrad,
        TF* const __restrict__ aflux,
        Boundary_type sw,
        const TF aval,
        const TF visc,
        const TF offset,
        bool set_flux_grad)
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
        if (set_flux_grad)
            set_bc_value_g<TF><<<grid2dGPU, block2dGPU>>>(aflux, -aval*visc, gd.icells, gd.jcells);
        cuda_check_error();
    }
    else if (sw == Boundary_type::Flux_type)
    {
        set_bc_value_g<TF><<<grid2dGPU, block2dGPU>>>(aflux, aval, gd.icells, gd.jcells);
        if (set_flux_grad)
            set_bc_value_g<TF><<<grid2dGPU, block2dGPU>>>(agrad, -aval*visc, gd.icells, gd.jcells);
        cuda_check_error();
    }
}

template<typename TF>
void Boundary<TF>::set_bc_2d_g(
        TF* const __restrict__ a,
        TF* const __restrict__ agrad,
        TF* const __restrict__ aflux,
        TF* const __restrict__ aval,
        Boundary_type sw,
        const TF visc,
        const TF offset,
        bool set_flux_grad)
{
    auto& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 grid2dGPU (gridi, gridj);
    dim3 block2dGPU(blocki, blockj);

    const TF no_offset = 0;
    const TF no_factor = 1;

    if (sw == Boundary_type::Dirichlet_type)
    {
        set_bc_2d_value_g<TF><<<grid2dGPU, block2dGPU>>>(a, aval, offset, no_factor, gd.icells, gd.jcells);
        cuda_check_error();
    }
    else if (sw == Boundary_type::Neumann_type)
    {
        set_bc_2d_value_g<TF><<<grid2dGPU, block2dGPU>>>(agrad, aval, no_offset, no_factor, gd.icells, gd.jcells);
        if (set_flux_grad)
            set_bc_2d_value_g<TF><<<grid2dGPU, block2dGPU>>>(aflux, aval, no_offset, -visc, gd.icells, gd.jcells);
        cuda_check_error();
    }
    else if (sw == Boundary_type::Flux_type)
    {
        set_bc_2d_value_g<TF><<<grid2dGPU, block2dGPU>>>(aflux, aval, no_offset, no_factor, gd.icells, gd.jcells);
        if (set_flux_grad)
            set_bc_2d_value_g<TF><<<grid2dGPU, block2dGPU>>>(agrad, aval, no_offset, -visc, gd.icells, gd.jcells);
        cuda_check_error();
    }
}

#ifdef USECUDA
template <typename TF>
void Boundary<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    const TF no_offset = 0.;
    const bool set_flux_grad = (swboundary == "default");

    if (swtimedep_sbot_2d)
    {
        auto tmp_cpu = fields.get_tmp();
        auto tmp_gpu = fields.get_tmp_g();

        unsigned long itime = timeloop.get_itime();

        const int blocki = gd.ithread_block;
        const int blockj = gd.jthread_block;
        const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
        const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

        dim3 grid2dGPU (gridi, gridj);
        dim3 block2dGPU(blocki, blockj);

        if (itime > itime_sbot_2d_next)
        {
            // Read new surface sbot fields
            const double ifactor = timeloop.get_ifactor();
            unsigned long iiotimeprec = timeloop.get_iiotimeprec();

            itime_sbot_2d_prev = itime_sbot_2d_next;
            itime_sbot_2d_next = itime_sbot_2d_prev + sbot_2d_loadtime*ifactor;
            const int iotime1 = int(itime_sbot_2d_next / iiotimeprec);

            int nerror = 0;

            for (auto& fld : sbot_2d_list)
            {
                // Swap 2D input fields at both CPU and GPU.
                // CPU is not strictly necessary, but it keeps the
                // CPU and GPU fields identical...
                sbot_2d_prev.at(fld) = sbot_2d_next.at(fld);

                copy_xy_g<TF><<<grid2dGPU, block2dGPU>>>(
                        sbot_2d_prev_g.at(fld),
                        sbot_2d_next_g.at(fld),
                        gd.icells, gd.jcells);

                // Read new time step
                char filename[256];
                std::string name = fld + "_bot_in";
                std::sprintf(filename, "%s.%07d", name.c_str(), iotime1);
                master.print_message("Loading \"%s\" ... ", filename);

                if (field3d_io.load_xy_slice(
                        sbot_2d_next.at(fld).data(), tmp_cpu->fld.data(), filename))
                {
                    master.print_message("FAILED\n");
                    nerror += 1;
                }
                else
                    master.print_message("OK\n");

                boundary_cyclic.exec_2d(sbot_2d_next.at(fld).data());

                // Copy new field to the GPU.
                cuda_safe_call(cudaMemcpy(
                        sbot_2d_next_g.at(fld),
                        sbot_2d_next.at(fld).data(),
                        gd.ijcells*sizeof(TF),
                        cudaMemcpyHostToDevice));
            }

            master.sum(&nerror, 1);
            if (nerror)
                throw std::runtime_error("Error loading time dependent sbot fields");
        }

        // Interpolate sbot to current time
        const TF fac1 = TF(itime - itime_sbot_2d_prev ) / TF(itime_sbot_2d_next - itime_sbot_2d_prev);
        const TF fac0 = TF(1) - fac1;

        for (auto& fld : sbot_2d_list)
        {
            // Interpolate to current time.
            interp_sbot_time_g<TF><<<grid2dGPU, block2dGPU>>>(
                    tmp_gpu->fld_bot_g,
                    sbot_2d_prev_g.at(fld),
                    sbot_2d_next_g.at(fld),
                    fac0, fac1,
                    gd.icells, gd.jcells);

            // Set new boundary conditions.
            set_bc_2d_g(
                    fields.sp.at(fld)->fld_bot_g,
                    fields.sp.at(fld)->grad_bot_g,
                    fields.sp.at(fld)->flux_bot_g,
                    tmp_gpu->fld_bot_g,
                    sbc.at(fld).bcbot,
                    fields.sp.at(fld)->visc,
                    no_offset,
                    set_flux_grad);
        }

        fields.release_tmp(tmp_cpu);
        fields.release_tmp_g(tmp_gpu);
    }
    else
    {
        for (auto& it : tdep_bc)
        {
            // Interpolate to current time.
            it.second->update_time_dependent(sbc.at(it.first).bot, timeloop);

            // Set new boundary conditions.
            set_bc_g(
                    fields.sp.at(it.first)->fld_bot_g,
                    fields.sp.at(it.first)->grad_bot_g,
                    fields.sp.at(it.first)->flux_bot_g,
                    sbc.at(it.first).bcbot,
                    sbc.at(it.first).bot,
                    fields.sp.at(it.first)->visc,
                    no_offset,
                    set_flux_grad);
        }
    }

    if (swtimedep_outflow)
        throw std::runtime_error("Time dependent outflow is not (yet) supported on the GPU...");
}

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
}

template<typename TF>
void Boundary<TF>::set_prognostic_outflow_bcs()
{
    // Overwrite here the ghost cells for the scalars with outflow BCs
    for (auto& s : scalar_outflow)
        boundary_outflow.exec(fields.sp.at(s)->fld_g, inflow_profiles_g.at(s));
}

template<typename TF>
TF* Boundary<TF>::get_z0m_g()
{
    throw std::runtime_error("Function get_z0m_g() not implemented in base boundary.");
}

template<typename TF>
TF* Boundary<TF>::get_dudz_g()
{
    throw std::runtime_error("Function get_dudz_g() not implemented in base boundary.");
}

template<typename TF>
TF* Boundary<TF>::get_dvdz_g()
{
    throw std::runtime_error("Function get_dvdz_g() not implemented in base boundary.");
}

template<typename TF>
TF* Boundary<TF>::get_dbdz_g()
{
    throw std::runtime_error("Function get_dbdz_g() not implemented in base boundary.");
}

template<typename TF>
void Boundary<TF>::prepare_device(Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();

    const int kmemsize = gd.kcells * sizeof(TF);
    const int ijmemsize = gd.ijcells * sizeof(TF);

    for (auto& scalar : scalar_outflow)
    {
        inflow_profiles_g.emplace(scalar, nullptr);
        cuda_safe_call(cudaMalloc(&inflow_profiles_g.at(scalar), kmemsize));
        cuda_safe_call(cudaMemcpy(inflow_profiles_g.at(scalar), inflow_profiles.at(scalar).data(), kmemsize, cudaMemcpyHostToDevice));
    }

    if (swtimedep_sbot_2d)
    {
        for (auto& scalar : sbot_2d_list)
        {
            sbot_2d_prev_g.emplace(scalar, nullptr);
            sbot_2d_next_g.emplace(scalar, nullptr);

            cuda_safe_call(cudaMalloc(&sbot_2d_prev_g.at(scalar), ijmemsize));
            cuda_safe_call(cudaMalloc(&sbot_2d_next_g.at(scalar), ijmemsize));

            cuda_safe_call(cudaMemcpy(sbot_2d_prev_g.at(scalar), sbot_2d_prev.at(scalar).data(), ijmemsize, cudaMemcpyHostToDevice));
            cuda_safe_call(cudaMemcpy(sbot_2d_next_g.at(scalar), sbot_2d_next.at(scalar).data(), ijmemsize, cudaMemcpyHostToDevice));
        }
    }
}

template<typename TF>
void Boundary<TF>::clear_device(Thermo<TF>& thermo)
{
    for(auto& it : tdep_bc)
        it.second->clear_device();

    for (auto& scalar : scalar_outflow)
        cuda_safe_call(cudaFree(inflow_profiles_g.at(scalar)));

    if (swtimedep_sbot_2d)
    {
        for (auto& scalar : sbot_2d_list)
        {
            cuda_safe_call(cudaFree(sbot_2d_prev_g.at(scalar)));
            cuda_safe_call(cudaFree(sbot_2d_next_g.at(scalar)));
        }
    }
}

template<typename TF>
void Boundary<TF>::forward_device(Thermo<TF>& thermo)
{
}

template<typename TF>
void Boundary<TF>::backward_device(Thermo<TF>& thermo)
{
}
#endif

template class Boundary<double>;
template class Boundary<float>;
