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

#include "master.h"
#include "grid.h"
#include "boundary_outflow.h"
#include "tools.h"

namespace
{
    template<typename TF, Edge_location location, Flow_direction direction> __global__
    void compute_inoutflow_2nd_g(
            TF* const restrict a, const TF* const restrict inflow_prof,
            const int istart, const int iend, const int igc,
            const int jstart, const int jend, const int jgc,
            const int icells, const int jcells, const int kcells,
            const int ijcells)
    {
        const int i  = blockIdx.x*blockDim.x + threadIdx.x;
        const int j  = blockIdx.y*blockDim.y + threadIdx.y;
        const int k  = blockIdx.z*blockDim.z + threadIdx.z;

        int ijk;
        int ijk_gc;
        int ijk_d;

        // Set the ghost cells using extrapolation.
        if (location == Edge_location::West || location == Edge_location::East)
        {
            if (k < kcells && j < jcells && i < igc)
            {
                if (location == Edge_location::West)
                {
                    ijk_d  = (istart    ) + j*icells + k*ijcells;
                    ijk    = (istart+i  ) + j*icells + k*ijcells;
                    ijk_gc = (istart-1-i) + j*icells + k*ijcells;
                }
                else if (location == Edge_location::East)
                {
                    ijk_d  = (iend-1  ) + j*icells + k*ijcells;
                    ijk    = (iend-1-i) + j*icells + k*ijcells;
                    ijk_gc = (iend+i  ) + j*icells + k*ijcells;
                }

                if (direction == Flow_direction::Inflow)
                    a[ijk_gc] = a[ijk_d] - (i+1)*TF(2)*(a[ijk_d] - inflow_prof[k]);
                else
                    a[ijk_gc] = a[ijk];
            }
        }
        else if (location == Edge_location::North || location == Edge_location::South)
        {
            if (k < kcells && j < jgc && i < icells)
            {
                if (location == Edge_location::South)
                {
                    ijk_d  = i + (jstart    )*icells + k*ijcells;
                    ijk    = i + (jstart+j  )*icells + k*ijcells;
                    ijk_gc = i + (jstart-1-j)*icells + k*ijcells;
                }
                else if (location == Edge_location::North)
                {
                    ijk_d  = i + (jend-1  )*icells + k*ijcells;
                    ijk    = i + (jend-1-j)*icells + k*ijcells;
                    ijk_gc = i + (jend+j  )*icells + k*ijcells;
                }

                if (direction == Flow_direction::Inflow)
                    a[ijk_gc] = a[ijk_d] - (j+1)*TF(2)*(a[ijk_d] - inflow_prof[k]);
                else
                    a[ijk_gc] = a[ijk];
            }
        }
    }

    template<typename TF> __global__
    void compute_outflow_4th(
            TF* const restrict a,
            const int iend, const int icells,
            const int jcells, const int kcells,
            const int ijcells)
    {
        const int j  = blockIdx.x*blockDim.x + threadIdx.x;
        const int k  = blockIdx.y*blockDim.y + threadIdx.y;

        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;

        if (j < jcells && k < kcells)
        {
            const int ijk = (iend-1) + j*icells + k*ijcells;
            a[ijk+ii1] = TF(2.)*a[ijk] - TF( 3./2.)*a[ijk-ii1] + TF(1./2.)*a[ijk-ii2];
            a[ijk+ii2] = TF(3.)*a[ijk] - TF( 7./2.)*a[ijk-ii1] + TF(3./2.)*a[ijk-ii2];
            a[ijk+ii3] = TF(5.)*a[ijk] - TF(15./2.)*a[ijk-ii1] + TF(7./2.)*a[ijk-ii2];
        }
    }

    template<typename TF> __global__
    void compute_inflow_4th(
            TF* const restrict a,
            const TF value,
            const int istart, const int icells,
            const int jcells, const int kcells,
            const int ijcells)
    {
        const int j  = blockIdx.x*blockDim.x + threadIdx.x;
        const int k  = blockIdx.y*blockDim.y + threadIdx.y;

        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;

        if (j < jcells && k < kcells)
        {
            const int ijk = istart + j*icells + k*ijcells;
            a[ijk-ii1] = value + TF( 9./8.)*a[ijk] - TF( 14./8.)*a[ijk+ii1] + TF( 5./8.)*a[ijk+ii2];
            a[ijk-ii2] = value + TF(33./8.)*a[ijk] - TF( 54./8.)*a[ijk+ii1] + TF(21./8.)*a[ijk+ii2];
            a[ijk-ii3] = value + TF(65./8.)*a[ijk] - TF(110./8.)*a[ijk+ii1] + TF(45./8.)*a[ijk+ii2];
        }
    }
}

#ifdef USECUDA
template<typename TF>
void Boundary_outflow<TF>::exec(
    TF* const restrict data,
    const TF* const restrict inflow_prof)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    if (grid.get_spatial_order() == Grid_order::Fourth)
    {
        const int blocki = gd.jthread_block;
        const int blockj = 64;

        const int gridi  = gd.jcells/blocki + (gd.jcells%blocki > 0);
        const int gridj  = gd.kcells/blockj + (gd.kcells%blockj > 0);

        dim3 grid2dGPU (gridi, gridj);
        dim3 block2dGPU(blocki, blockj);

        // Dirichlet BCs on west boundary, Neumann on east boundary,
        // cyclic BCs on north-south boundaries
        if (md.mpicoordx == 0)
            compute_inflow_4th<<<grid2dGPU, block2dGPU>>>(
                    data, TF(0),
                    gd.istart,
                    gd.icells, gd.jcells, gd.kcells,
                    gd.ijcells);

        if (md.mpicoordx == md.npx-1)
            compute_outflow_4th<<<grid2dGPU, block2dGPU>>>(
                    data,
                    gd.iend,
                    gd.icells, gd.jcells, gd.kcells,
                    gd.ijcells);

        cuda_check_error();
    }
    else if (grid.get_spatial_order() == Grid_order::Second)
    {
        // Vertical grid
        const int blockk = 4;
        const int gridk  = gd.kcells/blockk + (gd.kcells%blockk > 0);

        // Grid x-direction
        const int blocki_x = gd.igc;
        const int blockj_x = 64;

        const int gridi_x  = 1;
        const int gridj_x  = gd.jcells/blockj_x + (gd.jcells%blockj_x > 0);

        dim3 gridGPU_x (gridi_x, gridj_x, gridk);
        dim3 blockGPU_x(blocki_x, blockj_x, blockk);

        // Grid y-direction
        const int blocki_y = 64;
        const int blockj_y = gd.jgc;

        const int gridi_y  = gd.icells/blocki_y + (gd.icells%blocki_y > 0);
        const int gridj_y  = 1;

        dim3 gridGPU_y (gridi_y, gridj_y, gridk);
        dim3 blockGPU_y(blocki_y, blockj_y, blockk);

        if (md.mpicoordx == 0)
        {
            const Edge_location edge = Edge_location::West;

            if (flow_direction[edge] == Flow_direction::Inflow)
                compute_inoutflow_2nd_g<TF, edge, Flow_direction::Inflow><<<gridGPU_x, blockGPU_x>>>(
                        data, inflow_prof,
                        gd.istart, gd.iend, gd.igc,
                        gd.jstart, gd.jend, gd.kgc,
                        gd.icells, gd.jcells, gd.kcells,
                        gd.ijcells);
            else
                compute_inoutflow_2nd_g<TF, edge, Flow_direction::Outflow><<<gridGPU_x, blockGPU_x>>>(
                        data, inflow_prof,
                        gd.istart, gd.iend, gd.igc,
                        gd.jstart, gd.jend, gd.kgc,
                        gd.icells, gd.jcells, gd.kcells,
                        gd.ijcells);
        }

        if (md.mpicoordx == md.npx-1)
        {
            const Edge_location edge = Edge_location::East;

            if (flow_direction[edge] == Flow_direction::Inflow)
                compute_inoutflow_2nd_g<TF, edge, Flow_direction::Inflow><<<gridGPU_x, blockGPU_x>>>(
                        data, inflow_prof,
                        gd.istart, gd.iend, gd.igc,
                        gd.jstart, gd.jend, gd.kgc,
                        gd.icells, gd.jcells, gd.kcells,
                        gd.ijcells);
            else
                compute_inoutflow_2nd_g<TF, edge, Flow_direction::Outflow><<<gridGPU_x, blockGPU_x>>>(
                        data, inflow_prof,
                        gd.istart, gd.iend, gd.igc,
                        gd.jstart, gd.jend, gd.kgc,
                        gd.icells, gd.jcells, gd.kcells,
                        gd.ijcells);
        }

        if (md.mpicoordy == 0)
        {
            const Edge_location edge = Edge_location::South;

            if (flow_direction[edge] == Flow_direction::Inflow)
                compute_inoutflow_2nd_g<TF, edge, Flow_direction::Inflow><<<gridGPU_y, blockGPU_y>>>(
                        data, inflow_prof,
                        gd.istart, gd.iend, gd.igc,
                        gd.jstart, gd.jend, gd.kgc,
                        gd.icells, gd.jcells, gd.kcells,
                        gd.ijcells);
            else
                compute_inoutflow_2nd_g<TF, edge, Flow_direction::Outflow><<<gridGPU_y, blockGPU_y>>>(
                        data, inflow_prof,
                        gd.istart, gd.iend, gd.igc,
                        gd.jstart, gd.jend, gd.kgc,
                        gd.icells, gd.jcells, gd.kcells,
                        gd.ijcells);
        }

        if (md.mpicoordy == md.npy-1)
        {
            const Edge_location edge = Edge_location::North;

            if (flow_direction[edge] == Flow_direction::Inflow)
                compute_inoutflow_2nd_g<TF, edge, Flow_direction::Inflow><<<gridGPU_y, blockGPU_y>>>(
                        data, inflow_prof,
                        gd.istart, gd.iend, gd.igc,
                        gd.jstart, gd.jend, gd.kgc,
                        gd.icells, gd.jcells, gd.kcells,
                        gd.ijcells);
            else
                compute_inoutflow_2nd_g<TF, edge, Flow_direction::Outflow><<<gridGPU_y, blockGPU_y>>>(
                        data, inflow_prof,
                        gd.istart, gd.iend, gd.igc,
                        gd.jstart, gd.jend, gd.kgc,
                        gd.icells, gd.jcells, gd.kcells,
                        gd.ijcells);
        }

        cuda_check_error();
    }
}
#endif


#ifdef FLOAT_SINGLE
template class Boundary_outflow<float>;
#else
template class Boundary_outflow<double>;
#endif
