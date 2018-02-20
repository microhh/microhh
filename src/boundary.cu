/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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
#include "boundary.h"
#include "defines.h"
#include "model.h"
#include "timeloop.h"
#include "finite_difference.h"
#include "constants.h"
#include "tools.h"

using namespace Finite_difference::O4;

namespace
{
    template<typename TF> __global__
    void set_bc_value_g(TF* __restrict__ a, TF aval,
                  const int icells, const int icellsp, const int jcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*icellsp;
            a[ij] = aval;
        }
    }


    template<typename TF> __global__
    void calc_ghost_cells_bot_2nd_g(TF* __restrict__ a, TF* __restrict__ dzh, Boundary_type sw,
                                    TF* __restrict__ abot, TF* __restrict__ agradbot,
                                    const int icells, const int icellsp,
                                    const int jcells, const int kstart)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int kk  = icellsp*jcells;
        const int ij  = i + j*icellsp;
        const int ijk = i + j*icellsp + kstart*kk;

        if (i < icells && j < jcells)
        {
            if (sw == Boundary<TF>::Dirichlet_type)
                a[ijk-kk] = 2.*abot[ij] - a[ijk];

            else if (sw == Boundary<TF>::Neumann_type || sw == Boundary<TF>::Flux_type)
                a[ijk-kk] = -agradbot[ij]*dzh[kstart] + a[ijk];
        }
    }

    template<typename TF> __global__
    void calc_ghost_cells_top_2nd_g(TF* __restrict__ a, TF* __restrict__ dzh, const Boundary_type sw,
                                    TF* __restrict__ atop, TF* __restrict__ agradtop,
                                    const int icells, const int icellsp,
                                    const int jcells, const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int kk  = icellsp*jcells;
        const int ij  = i + j*icellsp;
        const int ijk = i + j*icellsp + (kend-1)*kk;

        if (i < icells && j < jcells)
        {
            if (sw == Boundary<TF>::Dirichlet_type)
                a[ijk+kk] = 2.*atop[ij] - a[ijk];

            else if (sw == Boundary<TF>::Neumann_type || sw == Boundary<TF>::Flux_type)
                a[ijk+kk] = agradtop[ij]*dzh[kend] + a[ijk];
        }
    }

    template<typename TF> __global__
    void calc_ghost_cells_bot_4th_g(TF* __restrict__ a, TF* __restrict__ z, const Boundary_type sw,
                                    TF* __restrict__ abot, TF* __restrict__ agradbot,                                    
                                    const int icells, const int icellsp,
                                    const int jcells, const int kstart)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int kk1 = 1*icellsp*jcells;
        const int kk2 = 2*icellsp*jcells;

        const int ij  = i + j*icellsp;
        const int ijk = i + j*icellsp + kstart*kk1;

        if (i < icells && j < jcells)
        {
            if (sw == Boundary<TF>::Dirichlet_type)
            {
                a[ijk-kk1] = (8./3.)*abot[ij] - 2.*a[ijk] + (1./3.)*a[ijk+kk1];
                a[ijk-kk2] = 8.*abot[ij] - 9.*a[ijk] + 2.*a[ijk+kk1];
            }

            else if (sw == Boundary<TF>::Neumann_type || sw == Boundary<TF>::Flux_type)
            {
                a[ijk-kk1] = -(1./24.)*grad4x(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk    ];
                a[ijk-kk2] = -(1./ 8.)*grad4x(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk+kk1];
            }
        }
    }

    template<typename TF> __global__
    void calc_ghost_cells_top_4th_g(TF* __restrict__ a, TF* __restrict__ z,const Boundary_type sw,
                                    TF* __restrict__ atop, TF* __restrict__ agradtop,                                    
                                    const int icells, const int icellsp,
                                    const int jcells, const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int kk1 = 1*icellsp*jcells;
        const int kk2 = 2*icellsp*jcells;

        const int ij  = i + j*icellsp;
        const int ijk = i + j*icellsp + (kend-1)*kk1;

        if( i < icells && j < jcells)
        {
            if (sw == Boundary<TF>::Dirichlet_type)
            {
                a[ijk+kk1] = (8./3.)*atop[ij] - 2.*a[ijk] + (1./3.)*a[ijk-kk1];
                a[ijk+kk2] = 8.*atop[ij] - 9.*a[ijk] + 2.*a[ijk-kk1];
            }

            else if (sw == Boundary<TF>::Neumann_type || sw == Boundary<TF>::Flux_type)
            {
                a[ijk+kk1] = (1./24.)*grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk    ];
                a[ijk+kk2] = (1./ 8.)*grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk-kk1];
            }
        }
    }

    template<typename TF> __global__
    void calc_ghost_cells_botw_4th_g(TF* __restrict__ w,
                                      const int icells, const int icellsp,
                                      const int jcells, const int kstart)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int kk1 = 1*icellsp*jcells;
        const int kk2 = 2*icellsp*jcells;
        const int kk3 = 3*icellsp*jcells;

        const int ijk = i + j*icellsp + kstart*kk1;

        if (i < icells && j < jcells)
        {
            w[ijk-kk1] = -6.*w[ijk+kk1] + 4.*w[ijk+kk2] - w[ijk+kk3];
        }
    }

    template<typename TF> __global__
    void calc_ghost_cells_topw_4th_g(TF* __restrict__ w,
                                      const int icells, const int icellsp,
                                      const int jcells, const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int kk1 = 1*icellsp*jcells;
        const int kk2 = 2*icellsp*jcells;
        const int kk3 = 3*icellsp*jcells;

        const int ijk = i + j*icellsp + kend*kk1;

        if (i < icells && j < jcells)
        {
            w[ijk+kk1] = -6.*w[ijk-kk1] + 4.*w[ijk-kk2] - w[ijk-kk3];
        }
    }

    template<typename TF> __global__
    void calc_ghost_cells_botw_cons_4th_g(TF* __restrict__ w,
                                      const int icells, const int icellsp,
                                      const int jcells, const int kstart)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int kk1 = 1*icellsp*jcells;
        const int kk2 = 2*icellsp*jcells;

        const int ijk = i + j*icellsp + kstart*kk1;

        if (i < icells && j < jcells)
        {
            w[ijk-kk1] = -w[ijk+kk1];
            w[ijk-kk2] = -w[ijk+kk2];
        }
    }

    template<typename TF> __global__
    void calc_ghost_cells_topw_cons_4th_g(TF* __restrict__ w,
                                      const int icells, const int icellsp,
                                      const int jcells, const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int kk1 = 1*icellsp*jcells;
        const int kk2 = 2*icellsp*jcells;

        const int ijk = i + j*icellsp + kend*kk1;

        if (i < icells && j < jcells)
        {
            w[ijk+kk1] = -w[ijk-kk1];
            w[ijk+kk2] = -w[ijk-kk2];
        }
    }
}

#ifdef USECUDA
template<typename TF>
void Boundary<TF>::exec()
{
    auto& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 grid2dGPU (gridi, gridj);
    dim3 block2dGPU(blocki, blockj);

    const int offs = gd.memoffset;

    // Cyclic boundary conditions, do this before the bottom BC's.
    grid.boundary_cyclic_g(&fields.mp.at("u")->fld_g[offs]);
    grid.boundary_cyclic_g(&fields.mp.at("v")->fld_g[offs]);
    grid.boundary_cyclic_g(&fields.mp.at("w")->fld_g[offs]);

    for (auto& it : fields.sp) 
        grid.boundary_cyclic_g(&it.second->fld_g[offs]);

    // Calculate the boundary values.
    update_bcs();

    if(grid.swspatialorder == "2")
    {
        calc_ghost_cells_bot_2nd_g<<<grid2dGPU, block2dGPU>>>(
            &fields.mp.at("u")->fld_g[offs], gd.dzh_g, mbcbot,
            &fields.mp.at("u")->fld_bot_g[offs], &fields.mp.at("u")->grad_bot_g[offs],
            gd.icells, gd.icellsp,
            gd.jcells, gd.kstart);
        cuda_check_error();

        calc_ghost_cells_top_2nd_g<<<grid2dGPU, block2dGPU>>>(
            &fields.mp.at("u")->fld_g[offs], gd.dzh_g, mbctop,
            &fields.mp.at("u")->fld_top_g[offs], &fields.mp.at("u")->grad_top_g[offs],
            gd.icells, gd.icellsp,
            gd.jcells, gd.kend);
        cuda_check_error();

        calc_ghost_cells_bot_2nd_g<<<grid2dGPU, block2dGPU>>>(
            &fields.mp.at("v")->fld_g[offs], gd.dzh_g, mbcbot,
            &fields.mp.at("v")->fld_bot_g[offs], &fields.mp.at("v")->grad_bot_g[offs],
            gd.icells, gd.icellsp,
            gd.jcells, gd.kstart);
        cuda_check_error();

        calc_ghost_cells_top_2nd_g<<<grid2dGPU, block2dGPU>>>(
            &fields.mp.at("v")->fld_g[offs], gd.dzh_g, mbctop,
            &fields.mp.at("v")->fld_top_g[offs], &fields.mp.at("v")->grad_top_g[offs],
            gd.icells, gd.icellsp,
            gd.jcells, gd.kend);
        cuda_check_error();

        for (auto& it : fields.sp)
        {
            calc_ghost_cells_bot_2nd_g<<<grid2dGPU, block2dGPU>>>(
                &it.second->fld_g[offs], gd.dzh_g, sbc.at(it.first).bcbot,
                &it.second->fld_bot_g[offs], &it.second->grad_bot_g[offs],
                gd.icells, gd.icellsp,
                gd.jcells, gd.kstart);
            cuda_check_error();

            calc_ghost_cells_top_2nd_g<<<grid2dGPU, block2dGPU>>>(
                &it.second->fld_g[offs], gd.dzh_g, sbc.at(it.first).bctop,
                &it.second->fld_top_g[offs], &it.second->grad_top_g[offs],
                gd.icells, gd.icellsp,
                gd.jcells, gd.kend);
            cuda_check_error();
        }
    }
    else if(grid.swspatialorder == "4")
    {
        calc_ghost_cells_bot_4th_g<<<grid2dGPU, block2dGPU>>>(
            &fields.mp.at("u")->fld_g[offs], gd.dzh_g, mbcbot,
            &fields.mp.at("u")->fld_bot_g[offs], &fields.mp.at("u")->grad_bot_g[offs],
            gd.icells, gd.icellsp,
            gd.jcells, gd.kstart);
        cuda_check_error();

        calc_ghost_cells_top_4th_g<<<grid2dGPU, block2dGPU>>>(
            &fields.mp.at("u")->fld_g[offs], gd.dzh_g, mbctop,
            &fields.mp.at("u")->fld_top_g[offs], &fields.mp.at("u")->grad_top_g[offs],
            gd.icells, gd.icellsp,
            gd.jcells, gd.kend);
        cuda_check_error();

        calc_ghost_cells_bot_4th_g<<<grid2dGPU, block2dGPU>>>(
            &fields.mp.at("v")->fld_g[offs], gd.dzh_g, mbcbot,
            &fields.mp.at("v")->fld_bot_g[offs], &fields.mp.at("v")->grad_bot_g[offs],
            gd.icells, gd.icellsp,
            gd.jcells, gd.kstart);
        cuda_check_error();

        calc_ghost_cells_top_4th_g<<<grid2dGPU, block2dGPU>>>(
            &fields.mp.at("v")->fld_g[offs], gd.dzh_g, mbctop,
            &fields.mp.at("v")->fld_top_g[offs], &fields.mp.at("v")->grad_top_g[offs],
            gd.icells, gd.icellsp,
            gd.jcells, gd.kend);
        cuda_check_error();

        for (auto& it : fields.sp)
        {
            calc_ghost_cells_bot_4th_g<<<grid2dGPU, block2dGPU>>>(
                &it.second->fld_g[offs], gd.dzh_g, sbc.at(it.first).bcbot,
                &it.second->fld_bot_g[offs], &it.second->grad_bot_g[offs],
                gd.icells, gd.icellsp,
                gd.jcells, gd.kstart);
            cuda_check_error();

            calc_ghost_cells_top_4th_g<<<grid2dGPU, block2dGPU>>>(
                &it.second->fld_g[offs], gd.dzh_g, sbc.at(it.first).bctop,
                &it.second->fld_top_g[offs], &it.second->grad_top_g[offs],
                gd.icells, gd.icellsp,
                gd.jcells, gd.kend);
            cuda_check_error();
        }
    }
}
#endif

#ifdef USECUDA
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

    const int offs = gd.memoffset;

    if (grid.swspatialorder == "4")
    {
        if (boundary_w_type == Boundary<TF>::Normal_type)
        {
            calc_ghost_cells_botw_4th_g<<<grid2dGPU, block2dGPU>>>(
                &fields.mp.at("w")->fld_g[offs],
                gd.icells, gd.icellsp,
                gd.jcells, gd.kstart);
            cuda_check_error();

            calc_ghost_cells_topw_4th_g<<<grid2dGPU, block2dGPU>>>(
                &fields.mp.at("w")->fld_g[offs],
                gd.icells, gd.icellsp,
                gd.jcells, gd.kend);
            cuda_check_error();
        }
        else if (boundary_w_type == Boundary<TF>::Conservation_type)
        {
            calc_ghost_cells_botw_cons_4th_g<<<grid2dGPU, block2dGPU>>>(
                &fields.mp.at("w")->fld_g[offs],
                gd.icells, gd.icellsp,
                gd.jcells, gd.kstart);
            cuda_check_error();

            calc_ghost_cells_topw_cons_4th_g<<<grid2dGPU, block2dGPU>>>(
                &fields.mp.at("w")->fld_g[offs],
                gd.icells, gd.icellsp,
                gd.jcells, gd.kend);
            cuda_check_error();
        }
    }
}
#endif

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

    const int offs = gd.memoffset;

    if (sw == Boundary<TF>::Dirichlet_type)
    {
        set_bc_value_g<<<grid2dGPU, block2dGPU>>>(&a[offs], aval-offset,    gd.icells, gd.icellsp, gd.jcells);
        cuda_check_error();
    }
    else if (sw == Boundary<TF>::Neumann_type)
    {
        set_bc_value_g<<<grid2dGPU, block2dGPU>>>(&agrad[offs], aval,       gd.icells, gd.icellsp, gd.jcells);
        set_bc_value_g<<<grid2dGPU, block2dGPU>>>(&aflux[offs], -aval*visc, gd.icells, gd.icellsp, gd.jcells);
        cuda_check_error();
    }
    else if (sw == Boundary<TF>::Flux_type)
    {
        set_bc_value_g<<<grid2dGPU, block2dGPU>>>(&aflux[offs], aval,       gd.icells, gd.icellsp, gd.jcells);
        set_bc_value_g<<<grid2dGPU, block2dGPU>>>(&agrad[offs], -aval*visc, gd.icells, gd.icellsp, gd.jcells);
        cuda_check_error();
    }
}

template class Boundary<double>;
template class Boundary<float>;
