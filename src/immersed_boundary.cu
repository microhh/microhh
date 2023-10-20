/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

#include <iostream>

#include "immersed_boundary.h"
#include "boundary.h"
#include "fields.h"
#include "tools.h"

namespace
{
    template<typename TF> __global__
    void set_ghost_cells_g(
            TF* const __restrict__ fld, const TF* const __restrict__ boundary_value,
            const TF* const __restrict__ c_idw, const TF* const __restrict__ c_idw_sum,
            const TF* const __restrict__ di,
            const int* const __restrict__ gi,  const int* const __restrict__ gj,  const int* const __restrict__ gk,
            const int* const __restrict__ ipi, const int* const __restrict__ ipj, const int* const __restrict__ ipk,
            Boundary_type bc, const TF visc, const int n_ghostcells, const int n_idw,
            const int icells, const int ijcells)
    {
        const int n  = blockIdx.x*blockDim.x + threadIdx.x;
        const int n_idw_loc = (bc == Boundary_type::Dirichlet_type) ? n_idw-1 : n_idw;

        if (n < n_ghostcells)
        {
            const int ijkg = gi[n] + gj[n]*icells + gk[n]*ijcells;

            // Sum the IDW coefficient times the value at the neighbouring grid points
            TF vI = TF(0);
            for (int i=0; i<n_idw_loc; ++i)
            {
                const int ii = i + n*n_idw;
                const int ijki = ipi[ii] + ipj[ii]*icells + ipk[ii]*ijcells;
                vI += c_idw[ii] * fld[ijki];
            }

            // For Dirichlet BCs, add the boundary value
            if (bc == Boundary_type::Dirichlet_type)
            {
                const int ii = n_idw-1 + n*n_idw;
                vI += c_idw[ii] * boundary_value[n];
            }

            vI /= c_idw_sum[n];

            // Set the ghost cells, depending on the IB boundary conditions
            if (bc == Boundary_type::Dirichlet_type)
                fld[ijkg] = 2*boundary_value[n] - vI;         // Image value reflected across IB
            else if (bc == Boundary_type::Neumann_type)
                fld[ijkg] = vI - boundary_value[n] * di[n];   // Image value minus gradient times distance
            else if (bc == Boundary_type::Flux_type)
            {
                const TF grad = -boundary_value[n] / visc;
                fld[ijkg] = vI - grad * di[n];             // Image value minus gradient times distance
            }
        }
    }
}

#ifdef USECUDA
template <typename TF>
void Immersed_boundary<TF>::exec_momentum()
{
    if (sw_ib == IB_type::Disabled)
        return;

    auto& gd = grid.get_grid_data();

    const int blocki = 256;

    const int nghost_u = ghost.at("u").nghost;
    const int nghost_v = ghost.at("v").nghost;
    const int nghost_w = ghost.at("w").nghost;

    const int gridu  = nghost_u / blocki + (nghost_u % blocki > 0);
    const int gridv  = nghost_v / blocki + (nghost_v % blocki > 0);
    const int gridw  = nghost_w / blocki + (nghost_w % blocki > 0);
    
    dim3 gridGPU_u(gridu);
    dim3 gridGPU_v(gridv);
    dim3 gridGPU_w(gridw);

    dim3 blockGPU(blocki);


    set_ghost_cells_g<TF><<<gridGPU_u, blockGPU>>>(
            fields.mp.at("u")->fld_g, ghost.at("u").mbot_g, 
            ghost.at("u").c_idw_g, ghost.at("u").c_idw_sum_g, ghost.at("u").di_g,
            ghost.at("u").i_g, ghost.at("u").j_g, ghost.at("u").k_g,
            ghost.at("u").ip_i_g, ghost.at("u").ip_j_g, ghost.at("u").ip_k_g,
            Boundary_type::Dirichlet_type, fields.visc, nghost_u, n_idw_points,
            gd.icells, gd.ijcells);
    cuda_check_error();

    set_ghost_cells_g<TF><<<gridGPU_v, blockGPU>>>(
            fields.mp.at("v")->fld_g, ghost.at("v").mbot_g,
            ghost.at("v").c_idw_g, ghost.at("v").c_idw_sum_g, ghost.at("v").di_g,
            ghost.at("v").i_g, ghost.at("v").j_g, ghost.at("v").k_g,
            ghost.at("v").ip_i_g, ghost.at("v").ip_j_g, ghost.at("v").ip_k_g,
            Boundary_type::Dirichlet_type, fields.visc, nghost_v, n_idw_points,
            gd.icells, gd.ijcells);
    cuda_check_error();

    set_ghost_cells_g<TF><<<gridGPU_w, blockGPU>>>(
            fields.mp.at("w")->fld_g, ghost.at("w").mbot_g,
            ghost.at("w").c_idw_g, ghost.at("w").c_idw_sum_g, ghost.at("w").di_g,
            ghost.at("w").i_g, ghost.at("w").j_g, ghost.at("w").k_g,
            ghost.at("w").ip_i_g, ghost.at("w").ip_j_g, ghost.at("w").ip_k_g,
            Boundary_type::Dirichlet_type, fields.visc, nghost_w, n_idw_points,
            gd.icells, gd.ijcells);
    cuda_check_error();

    boundary_cyclic.exec_g(fields.mp.at("u")->fld_g);
    boundary_cyclic.exec_g(fields.mp.at("v")->fld_g);
    boundary_cyclic.exec_g(fields.mp.at("w")->fld_g);
}

template <typename TF>
void Immersed_boundary<TF>::exec_scalars()
{
    if (sw_ib == IB_type::Disabled)
        return;

    if (fields.sp.size() > 0)
    {
        auto& gd = grid.get_grid_data();

        const int blocki = 256;
        const int nghost_s = ghost.at("s").nghost;
        const int grid = nghost_s / blocki + (nghost_s % blocki > 0);
        
        dim3 gridGPU(grid);
        dim3 blockGPU(blocki);

        for (auto& it : fields.sp)
        {
            set_ghost_cells_g<TF><<<gridGPU, blockGPU>>>(
                    it.second->fld_g, ghost.at("s").sbot_g.at(it.first),
                    ghost.at("s").c_idw_g, ghost.at("s").c_idw_sum_g, ghost.at("s").di_g,
                    ghost.at("s").i_g, ghost.at("s").j_g, ghost.at("s").k_g,
                    ghost.at("s").ip_i_g, ghost.at("s").ip_j_g, ghost.at("s").ip_k_g,
                    sbcbot, it.second->visc, nghost_s, n_idw_points,
                    gd.icells, gd.ijcells);

            boundary_cyclic.exec_g(it.second->fld_g);
        }
    }
}
#endif

template <typename TF>
void Immersed_boundary<TF>::prepare_device()
{
    if (sw_ib == IB_type::Disabled)
        return;

    auto& gd = grid.get_grid_data();

    // Allocate and copy data for all ghost cells (u, v, w, and optionally scalars)
    for (auto& g : ghost)
    {
        const int nghost = g.second.nghost;

        const int imemsize_1d = nghost*sizeof(int);
        const int fmemsize_1d = nghost*sizeof(TF);

        const int imemsize_2d = nghost*n_idw_points*sizeof(int);
        const int fmemsize_2d = nghost*n_idw_points*sizeof(TF);

        // Allocate
        cuda_safe_call(cudaMalloc(&g.second.i_g, imemsize_1d));
        cuda_safe_call(cudaMalloc(&g.second.j_g, imemsize_1d));
        cuda_safe_call(cudaMalloc(&g.second.k_g, imemsize_1d));

        cuda_safe_call(cudaMalloc(&g.second.xb_g, fmemsize_1d));
        cuda_safe_call(cudaMalloc(&g.second.yb_g, fmemsize_1d));
        cuda_safe_call(cudaMalloc(&g.second.zb_g, fmemsize_1d));

        cuda_safe_call(cudaMalloc(&g.second.xi_g, fmemsize_1d));
        cuda_safe_call(cudaMalloc(&g.second.yi_g, fmemsize_1d));
        cuda_safe_call(cudaMalloc(&g.second.zi_g, fmemsize_1d));

        cuda_safe_call(cudaMalloc(&g.second.di_g, fmemsize_1d));

        cuda_safe_call(cudaMalloc(&g.second.ip_i_g, imemsize_2d));
        cuda_safe_call(cudaMalloc(&g.second.ip_j_g, imemsize_2d));
        cuda_safe_call(cudaMalloc(&g.second.ip_k_g, imemsize_2d));
        cuda_safe_call(cudaMalloc(&g.second.ip_d_g, fmemsize_2d));

        cuda_safe_call(cudaMalloc(&g.second.c_idw_g, fmemsize_2d));
        cuda_safe_call(cudaMalloc(&g.second.c_idw_sum_g, fmemsize_1d));

        if (g.first == "u" || g.first == "v" || g.first == "w")
            cuda_safe_call(cudaMalloc(&g.second.mbot_g, fmemsize_1d));
        else
            for (auto& it : g.second.sbot)
                cuda_safe_call(cudaMalloc(&g.second.sbot_g[it.first], fmemsize_1d));

        // Forward copy
        cuda_safe_call(cudaMemcpy(g.second.i_g, g.second.i.data(), imemsize_1d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(g.second.j_g, g.second.j.data(), imemsize_1d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(g.second.k_g, g.second.k.data(), imemsize_1d, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(g.second.xb_g, g.second.xb.data(), fmemsize_1d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(g.second.yb_g, g.second.yb.data(), fmemsize_1d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(g.second.zb_g, g.second.zb.data(), fmemsize_1d, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(g.second.xi_g, g.second.xi.data(), fmemsize_1d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(g.second.yi_g, g.second.yi.data(), fmemsize_1d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(g.second.zi_g, g.second.zi.data(), fmemsize_1d, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(g.second.di_g, g.second.di.data(), fmemsize_1d, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(g.second.ip_i_g, g.second.ip_i.data(), imemsize_2d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(g.second.ip_j_g, g.second.ip_j.data(), imemsize_2d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(g.second.ip_k_g, g.second.ip_k.data(), imemsize_2d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(g.second.ip_d_g, g.second.ip_d.data(), fmemsize_2d, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(g.second.c_idw_g, g.second.c_idw.data(), fmemsize_2d, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(g.second.c_idw_sum_g, g.second.c_idw_sum.data(), fmemsize_1d, cudaMemcpyHostToDevice));

        if (g.first == "u" || g.first == "v" || g.first == "w")
            cuda_safe_call(cudaMemcpy(g.second.mbot_g, g.second.mbot.data(), fmemsize_1d, cudaMemcpyHostToDevice));
        else
            for (auto& it : g.second.sbot)
                cuda_safe_call(cudaMemcpy(g.second.sbot_g.at(it.first), g.second.sbot.at(it.first).data(), 
                               fmemsize_1d, cudaMemcpyHostToDevice));
    }    
}

template <typename TF>
void Immersed_boundary<TF>::clear_device()
{
    if (sw_ib == IB_type::Disabled)
        return;

    // De-allocate all ghost cell properties
    for (auto& g : ghost)
    {
        cuda_safe_call(cudaFree(g.second.i_g));
        cuda_safe_call(cudaFree(g.second.j_g));
        cuda_safe_call(cudaFree(g.second.k_g));

        cuda_safe_call(cudaFree(g.second.xb_g));
        cuda_safe_call(cudaFree(g.second.yb_g));
        cuda_safe_call(cudaFree(g.second.zb_g));

        cuda_safe_call(cudaFree(g.second.xi_g));
        cuda_safe_call(cudaFree(g.second.yi_g));
        cuda_safe_call(cudaFree(g.second.zi_g));

        cuda_safe_call(cudaFree(g.second.di_g));

        cuda_safe_call(cudaFree(g.second.ip_i_g));
        cuda_safe_call(cudaFree(g.second.ip_j_g));
        cuda_safe_call(cudaFree(g.second.ip_k_g));
        cuda_safe_call(cudaFree(g.second.ip_d_g));

        cuda_safe_call(cudaFree(g.second.c_idw_g));
        cuda_safe_call(cudaFree(g.second.c_idw_sum_g));
    }
}


#ifdef FLOAT_SINGLE
template class Immersed_boundary<float>;
#else
template class Immersed_boundary<double>;
#endif
