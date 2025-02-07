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

#include "master.h"
#include "grid.h"
#include "boundary_cyclic.h"
#include "boundary_cyclic_kernels.h"


template<typename TF>
Boundary_cyclic<TF>::Boundary_cyclic(Master& masterin, Grid<TF>& gridin) :
    master(masterin),
    grid(gridin)
{
}


template<typename TF>
Boundary_cyclic<TF>::~Boundary_cyclic()
{
}


template<typename TF>
void Boundary_cyclic<TF>::exec(TF* const restrict data, Edge edge)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    TF* buffer_send = grid.get_tmp_3d();
    TF* buffer_recv = grid.get_tmp_3d();

    const bool use_gpu = false;

    Boundary_cyclic_kernels::cyclic_kernel<TF, use_gpu>(
            data,
            buffer_send,
            buffer_recv,
            edge,
            gd.istart, gd.iend, gd.jstart, gd.jend,
            gd.icells, gd.jcells, gd.kcells,
            gd.icells, gd.ijcells,
            md);

    grid.release_tmp_3d(buffer_send);
    grid.release_tmp_3d(buffer_recv);
}


template<typename TF>
void Boundary_cyclic<TF>::exec_2d(TF* const restrict data)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    // Turn into 2D
    TF* buffer_send = grid.get_tmp_2d();
    TF* buffer_recv = grid.get_tmp_2d();

    const bool use_gpu = false;

    Boundary_cyclic_kernels::cyclic_kernel<TF, use_gpu>(
            data,
            buffer_send,
            buffer_recv,
            Edge::Both_edges,
            gd.istart, gd.iend, gd.jstart, gd.jend,
            gd.icells, gd.jcells, 1,
            gd.icells, gd.ijcells,
            md);

    grid.release_tmp_2d(buffer_send);
    grid.release_tmp_2d(buffer_recv);
}


template<typename TF>
void Boundary_cyclic<TF>::exec(unsigned int* const restrict data, Edge edge)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    TF* buffer_send = grid.get_tmp_3d();
    TF* buffer_recv = grid.get_tmp_3d();

    const bool use_gpu = false;

    Boundary_cyclic_kernels::cyclic_kernel<unsigned int, use_gpu>(
            data,
            reinterpret_cast<unsigned int*>(buffer_send),
            reinterpret_cast<unsigned int*>(buffer_recv),
            edge,
            gd.istart, gd.iend, gd.jstart, gd.jend,
            gd.icells, gd.jcells, gd.kcells,
            gd.icells, gd.ijcells,
            md);

    grid.release_tmp_3d(buffer_send);
    grid.release_tmp_3d(buffer_recv);
}


template<typename TF>
void Boundary_cyclic<TF>::exec_2d(unsigned int* const restrict data)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    TF* buffer_send = grid.get_tmp_2d();
    TF* buffer_recv = grid.get_tmp_2d();

    const bool use_gpu = false;

    Boundary_cyclic_kernels::cyclic_kernel<unsigned int, use_gpu>(
            data,
            reinterpret_cast<unsigned int*>(buffer_send),
            reinterpret_cast<unsigned int*>(buffer_recv),
            Edge::Both_edges,
            gd.istart, gd.iend, gd.jstart, gd.jend,
            gd.icells, gd.jcells, gd.kcells,
            gd.icells, gd.ijcells,
            md);

    grid.release_tmp_2d(buffer_send);
    grid.release_tmp_2d(buffer_recv);
}


#ifdef FLOAT_SINGLE
template class Boundary_cyclic<float>;
#else
template class Boundary_cyclic<double>;
#endif
