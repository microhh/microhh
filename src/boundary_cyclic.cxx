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

#include "master.h"
#include "grid.h"
#include "boundary_cyclic.h"

template<typename TF>
Boundary_cyclic<TF>::Boundary_cyclic(Master& masterin, Grid<TF>& gridin) :
    master(masterin),
    grid(gridin),
    mpi_types_allocated(false)
{
}

template<typename TF>
Boundary_cyclic<TF>::~Boundary_cyclic()
{
    exit_mpi();
}

template<typename TF>
void Boundary_cyclic<TF>::init()
{
    init_mpi();
}

#ifdef USEMPI
namespace
{
    template<typename TF> MPI_Datatype mpi_fp_type();
    template<> MPI_Datatype mpi_fp_type<double>() { return MPI_DOUBLE; }
    template<> MPI_Datatype mpi_fp_type<float>() { return MPI_FLOAT; }
}

template<typename TF>
void Boundary_cyclic<TF>::init_mpi()
{
    auto& gd = grid.get_grid_data();

    // create the MPI types for the cyclic boundary conditions
    int datacount, datablock, datastride;

    // east west
    datacount  = gd.jcells*gd.kcells;
    datablock  = gd.igc;
    datastride = gd.icells;
    MPI_Type_vector(datacount, datablock, datastride, mpi_fp_type<TF>(), &eastwestedge);
    MPI_Type_commit(&eastwestedge);
    MPI_Type_vector(datacount, datablock, datastride, MPI_UNSIGNED, &eastwestedge_uint);
    MPI_Type_commit(&eastwestedge_uint);

    // north south
    datacount  = gd.kcells;
    datablock  = gd.icells*gd.jgc;
    datastride = gd.icells*gd.jcells;
    MPI_Type_vector(datacount, datablock, datastride, mpi_fp_type<TF>(), &northsouthedge);
    MPI_Type_commit(&northsouthedge);
    MPI_Type_vector(datacount, datablock, datastride, MPI_UNSIGNED, &northsouthedge_uint);
    MPI_Type_commit(&northsouthedge_uint);

    // east west 2d
    datacount  = gd.jcells;
    datablock  = gd.igc;
    datastride = gd.icells;
    MPI_Type_vector(datacount, datablock, datastride, mpi_fp_type<TF>(), &eastwestedge2d);
    MPI_Type_commit(&eastwestedge2d);
    MPI_Type_vector(datacount, datablock, datastride, MPI_UNSIGNED, &eastwestedge2d_uint);
    MPI_Type_commit(&eastwestedge2d_uint);

    // north south 2d
    datacount  = 1;
    datablock  = gd.icells*gd.jgc;
    datastride = gd.icells*gd.jcells;
    MPI_Type_vector(datacount, datablock, datastride, mpi_fp_type<TF>(), &northsouthedge2d);
    MPI_Type_commit(&northsouthedge2d);
    MPI_Type_vector(datacount, datablock, datastride, MPI_UNSIGNED, &northsouthedge2d_uint);
    MPI_Type_commit(&northsouthedge2d_uint);

    mpi_types_allocated = true;
}

template<typename TF>
void Boundary_cyclic<TF>::exit_mpi()
{
    if (mpi_types_allocated)
    {
        MPI_Type_free(&eastwestedge);
        MPI_Type_free(&northsouthedge);

        MPI_Type_free(&eastwestedge2d);
        MPI_Type_free(&northsouthedge2d);
    }
}

template<typename TF>
void Boundary_cyclic<TF>::exec(TF* const restrict data, Edge edge)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    const int ncount = 1;

    if (edge == Edge::East_west_edge || edge == Edge::Both_edges)
    {
        // Communicate east-west edges.
        const int eastout = gd.iend-gd.igc;
        const int westin  = 0;
        const int westout = gd.istart;
        const int eastin  = gd.iend;

        // Send and receive the ghost cells in east-west direction.
        MPI_Isend(&data[eastout], ncount, eastwestedge, md.neast, 1, md.commxy, master.get_request_ptr());
        MPI_Irecv(&data[ westin], ncount, eastwestedge, md.nwest, 1, md.commxy, master.get_request_ptr());
        MPI_Isend(&data[westout], ncount, eastwestedge, md.nwest, 2, md.commxy, master.get_request_ptr());
        MPI_Irecv(&data[ eastin], ncount, eastwestedge, md.neast, 2, md.commxy, master.get_request_ptr());
        master.wait_all();
    }

    if (edge == Edge::North_south_edge || edge == Edge::Both_edges)
    {
        // If the run is 3D, perform the cyclic boundary routine for the north-south direction.
        if (gd.jtot > 1)
        {
            // Communicate north-south edges.
            const int northout = (gd.jend-gd.jgc)*gd.icells;
            const int southin  = 0;
            const int southout = gd.jstart*gd.icells;
            const int northin  = gd.jend  *gd.icells;

            // Send and receive the ghost cells in the north-south direction.
            MPI_Isend(&data[northout], ncount, northsouthedge, md.nnorth, 1, md.commxy, master.get_request_ptr());
            MPI_Irecv(&data[ southin], ncount, northsouthedge, md.nsouth, 1, md.commxy, master.get_request_ptr());
            MPI_Isend(&data[southout], ncount, northsouthedge, md.nsouth, 2, md.commxy, master.get_request_ptr());
            MPI_Irecv(&data[ northin], ncount, northsouthedge, md.nnorth, 2, md.commxy, master.get_request_ptr());
            master.wait_all();
        }
        // In case of 2D, fill all the ghost cells in the y-direction with the same value.
        else
        {
            const int jj = gd.icells;
            const int kk = gd.icells*gd.jcells;

            for (int k=gd.kstart; k<gd.kend; ++k)
                for (int j=0; j<gd.jgc; ++j)
                    #pragma ivdep
                    for (int i=0; i<gd.icells; ++i)
                    {
                        const int ijkref   = i + gd.jstart*jj   + k*kk;
                        const int ijknorth = i + j*jj           + k*kk;
                        const int ijksouth = i + (gd.jend+j)*jj + k*kk;
                        data[ijknorth] = data[ijkref];
                        data[ijksouth] = data[ijkref];
                    }
        }
    }
}

template<typename TF>
void Boundary_cyclic<TF>::exec_2d(TF* const restrict data)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    const int ncount = 1;

    // Communicate east-west edges.
    const int eastout = gd.iend-gd.igc;
    const int westin  = 0;
    const int westout = gd.istart;
    const int eastin  = gd.iend;

    // Communicate north-south edges.
    const int northout = (gd.jend-gd.jgc)*gd.icells;
    const int southin  = 0;
    const int southout = gd.jstart*gd.icells;
    const int northin  = gd.jend  *gd.icells;

    // First, send and receive the ghost cells in east-west direction.
    MPI_Isend(&data[eastout], ncount, eastwestedge2d, md.neast, 1, md.commxy, master.get_request_ptr());
    MPI_Irecv(&data[ westin], ncount, eastwestedge2d, md.nwest, 1, md.commxy, master.get_request_ptr());
    MPI_Isend(&data[westout], ncount, eastwestedge2d, md.nwest, 2, md.commxy, master.get_request_ptr());
    MPI_Irecv(&data[ eastin], ncount, eastwestedge2d, md.neast, 2, md.commxy, master.get_request_ptr());
    master.wait_all();

    // If the run is 3D, apply the BCs.
    if (gd.jtot > 1)
    {
        // Second, send and receive the ghost cells in the north-south direction.
        MPI_Isend(&data[northout], ncount, northsouthedge2d, md.nnorth, 1, md.commxy, master.get_request_ptr());
        MPI_Irecv(&data[ southin], ncount, northsouthedge2d, md.nsouth, 1, md.commxy, master.get_request_ptr());
        MPI_Isend(&data[southout], ncount, northsouthedge2d, md.nsouth, 2, md.commxy, master.get_request_ptr());
        MPI_Irecv(&data[ northin], ncount, northsouthedge2d, md.nnorth, 2, md.commxy, master.get_request_ptr());
        master.wait_all();
    }
    // In case of 2D, fill all the ghost cells with the current value.
    else
    {
        // Local copies for fast performance in loop.
        const int jj = gd.icells;
        const int jstart = gd.jstart;
        const int jend = gd.jend;

        for (int j=0; j<gd.jgc; ++j)
            #pragma ivdep
            for (int i=0; i<gd.icells; ++i)
            {
                const int ijref   = i + jstart*jj;
                const int ijnorth = i + j*jj;
                const int ijsouth = i + (jend+j)*jj;
                data[ijnorth] = data[ijref];
                data[ijsouth] = data[ijref];
            }
    }
}

template<typename TF>
void Boundary_cyclic<TF>::exec(unsigned int* const restrict data, Edge edge)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    const int ncount = 1;

    if (edge == Edge::East_west_edge || edge == Edge::Both_edges)
    {
        // Communicate east-west edges.
        const int eastout = gd.iend-gd.igc;
        const int westin  = 0;
        const int westout = gd.istart;
        const int eastin  = gd.iend;

        // Send and receive the ghost cells in east-west direction.
        MPI_Isend(&data[eastout], ncount, eastwestedge_uint, md.neast, 1, md.commxy, master.get_request_ptr());
        MPI_Irecv(&data[ westin], ncount, eastwestedge_uint, md.nwest, 1, md.commxy, master.get_request_ptr());
        MPI_Isend(&data[westout], ncount, eastwestedge_uint, md.nwest, 2, md.commxy, master.get_request_ptr());
        MPI_Irecv(&data[ eastin], ncount, eastwestedge_uint, md.neast, 2, md.commxy, master.get_request_ptr());
        master.wait_all();
    }

    if (edge == Edge::North_south_edge || edge == Edge::Both_edges)
    {
        // If the run is 3D, perform the cyclic boundary routine for the north-south direction.
        if (gd.jtot > 1)
        {
            // Communicate north-south edges.
            const int northout = (gd.jend-gd.jgc)*gd.icells;
            const int southin  = 0;
            const int southout = gd.jstart*gd.icells;
            const int northin  = gd.jend  *gd.icells;

            // Send and receive the ghost cells in the north-south direction.
            MPI_Isend(&data[northout], ncount, northsouthedge_uint, md.nnorth, 1, md.commxy, master.get_request_ptr());
            MPI_Irecv(&data[ southin], ncount, northsouthedge_uint, md.nsouth, 1, md.commxy, master.get_request_ptr());
            MPI_Isend(&data[southout], ncount, northsouthedge_uint, md.nsouth, 2, md.commxy, master.get_request_ptr());
            MPI_Irecv(&data[ northin], ncount, northsouthedge_uint, md.nnorth, 2, md.commxy, master.get_request_ptr());
            master.wait_all();
        }
        // In case of 2D, fill all the ghost cells in the y-direction with the same value.
        else
        {
            const int jj = gd.icells;
            const int kk = gd.icells*gd.jcells;

            for (int k=gd.kstart; k<gd.kend; ++k)
                for (int j=0; j<gd.jgc; ++j)
                    #pragma ivdep
                    for (int i=0; i<gd.icells; ++i)
                    {
                        const int ijkref   = i + gd.jstart*jj   + k*kk;
                        const int ijknorth = i + j*jj           + k*kk;
                        const int ijksouth = i + (gd.jend+j)*jj + k*kk;
                        data[ijknorth] = data[ijkref];
                        data[ijksouth] = data[ijkref];
                    }
        }
    }
}

template<typename TF>
void Boundary_cyclic<TF>::exec_2d(unsigned int* const restrict data)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    const int ncount = 1;

    // Communicate east-west edges.
    const int eastout = gd.iend-gd.igc;
    const int westin  = 0;
    const int westout = gd.istart;
    const int eastin  = gd.iend;

    // Communicate north-south edges.
    const int northout = (gd.jend-gd.jgc)*gd.icells;
    const int southin  = 0;
    const int southout = gd.jstart*gd.icells;
    const int northin  = gd.jend  *gd.icells;

    // First, send and receive the ghost cells in east-west direction.
    MPI_Isend(&data[eastout], ncount, eastwestedge2d_uint, md.neast, 1, md.commxy, master.get_request_ptr());
    MPI_Irecv(&data[ westin], ncount, eastwestedge2d_uint, md.nwest, 1, md.commxy, master.get_request_ptr());
    MPI_Isend(&data[westout], ncount, eastwestedge2d_uint, md.nwest, 2, md.commxy, master.get_request_ptr());
    MPI_Irecv(&data[ eastin], ncount, eastwestedge2d_uint, md.neast, 2, md.commxy, master.get_request_ptr());
    master.wait_all();

    // If the run is 3D, apply the BCs.
    if (gd.jtot > 1)
    {
        // Second, send and receive the ghost cells in the north-south direction.
        MPI_Isend(&data[northout], ncount, northsouthedge2d_uint, md.nnorth, 1, md.commxy, master.get_request_ptr());
        MPI_Irecv(&data[ southin], ncount, northsouthedge2d_uint, md.nsouth, 1, md.commxy, master.get_request_ptr());
        MPI_Isend(&data[southout], ncount, northsouthedge2d_uint, md.nsouth, 2, md.commxy, master.get_request_ptr());
        MPI_Irecv(&data[ northin], ncount, northsouthedge2d_uint, md.nnorth, 2, md.commxy, master.get_request_ptr());
        master.wait_all();
    }
    // In case of 2D, fill all the ghost cells with the current value.
    else
    {
        // Local copies for fast performance in loop.
        const int jj = gd.icells;
        const int jstart = gd.jstart;
        const int jend = gd.jend;

        for (int j=0; j<gd.jgc; ++j)
            #pragma ivdep
            for (int i=0; i<gd.icells; ++i)
            {
                const int ijref   = i + jstart*jj;
                const int ijnorth = i + j*jj;
                const int ijsouth = i + (jend+j)*jj;
                data[ijnorth] = data[ijref];
                data[ijsouth] = data[ijref];
            }
    }
}

#else

template<typename TF>
void Boundary_cyclic<TF>::init_mpi()
{
}

template<typename TF>
void Boundary_cyclic<TF>::exit_mpi()
{
}

template<typename TF>
void Boundary_cyclic<TF>::exec(TF* restrict data, Edge edge)
{
    auto& gd = grid.get_grid_data();

    const int jj = gd.icells;
    const int kk = gd.icells*gd.jcells;

    if (edge == Edge::East_west_edge || edge == Edge::Both_edges)
    {
        // first, east west boundaries
        for (int k=0; k<gd.kcells; ++k)
            for (int j=0; j<gd.jcells; ++j)
                #pragma ivdep
                for (int i=0; i<gd.igc; ++i)
                {
                    const int ijk0 = i          + j*jj + k*kk;
                    const int ijk1 = gd.iend-gd.igc+i + j*jj + k*kk;
                    data[ijk0] = data[ijk1];
                }

        for (int k=0; k<gd.kcells; ++k)
            for (int j=0; j<gd.jcells; ++j)
                #pragma ivdep
                for (int i=0; i<gd.igc; ++i)
                {
                    const int ijk0 = i+gd.iend   + j*jj + k*kk;
                    const int ijk1 = i+gd.istart + j*jj + k*kk;
                    data[ijk0] = data[ijk1];
                }
    }

    if (edge == Edge::North_south_edge || edge == Edge::Both_edges)
    {
        // if the run is 3D, apply the BCs
        if (gd.jtot > 1)
        {
            // second, send and receive the ghost cells in the north-south direction
            for (int k=0; k<gd.kcells; ++k)
                for (int j=0; j<gd.jgc; ++j)
                    #pragma ivdep
                    for (int i=0; i<gd.icells; ++i)
                    {
                        const int ijk0 = i + j                 *jj + k*kk;
                        const int ijk1 = i + (gd.jend-gd.jgc+j)*jj + k*kk;
                        data[ijk0] = data[ijk1];
                    }

            for (int k=0; k<gd.kcells; ++k)
                for (int j=0; j<gd.jgc; ++j)
                    #pragma ivdep
                    for (int i=0; i<gd.icells; ++i)
                    {
                        const int ijk0 = i + (j+gd.jend  )*jj + k*kk;
                        const int ijk1 = i + (j+gd.jstart)*jj + k*kk;
                        data[ijk0] = data[ijk1];
                    }
        }
        // in case of 2D, fill all the ghost cells with the current value
        else
        {
            for (int k=gd.kstart; k<gd.kend; ++k)
                for (int j=0; j<gd.jgc; ++j)
                    #pragma ivdep
                    for (int i=0; i<gd.icells; ++i)
                    {
                        const int ijkref   = i + gd.jstart*jj   + k*kk;
                        const int ijknorth = i + j*jj           + k*kk;
                        const int ijksouth = i + (gd.jend+j)*jj + k*kk;
                        data[ijknorth] = data[ijkref];
                        data[ijksouth] = data[ijkref];
                    }
        }
    }
}

template<typename TF>
void Boundary_cyclic<TF>::exec_2d(TF* restrict data)
{
    auto& gd = grid.get_grid_data();

    const int jj = gd.icells;

    // First, east west boundaries.
    for (int j=0; j<gd.jcells; ++j)
        #pragma ivdep
        for (int i=0; i<gd.igc; ++i)
        {
            const int ij0 = i                + j*jj;
            const int ij1 = gd.iend-gd.igc+i + j*jj;
            data[ij0] = data[ij1];
        }

    for (int j=0; j<gd.jcells; ++j)
        #pragma ivdep
        for (int i=0; i<gd.igc; ++i)
        {
            const int ij0 = i+gd.iend   + j*jj;
            const int ij1 = i+gd.istart + j*jj;
            data[ij0] = data[ij1];
        }

    // If the run is 3D, apply the BCs.
    if (gd.jtot > 1)
    {
        // Second, send and receive the ghost cells in the north-south direction.
        for (int j=0; j<gd.jgc; ++j)
            #pragma ivdep
            for (int i=0; i<gd.icells; ++i)
            {
                const int ij0 = i + j                 *jj;
                const int ij1 = i + (gd.jend-gd.jgc+j)*jj;
                data[ij0] = data[ij1];
            }

        for (int j=0; j<gd.jgc; ++j)
            #pragma ivdep
            for (int i=0; i<gd.icells; ++i)
            {
                const int ij0 = i + (j+gd.jend  )*jj;
                const int ij1 = i + (j+gd.jstart)*jj;
                data[ij0] = data[ij1];
            }
    }
    // In case of 2D, fill all the ghost cells with the current value.
    else
    {
        for (int j=0; j<gd.jgc; ++j)
            #pragma ivdep
            for (int i=0; i<gd.icells; ++i)
            {
                const int ijref   = i + gd.jstart*jj;
                const int ijnorth = i + j*jj;
                const int ijsouth = i + (gd.jend+j)*jj;
                data[ijnorth] = data[ijref];
                data[ijsouth] = data[ijref];
            }
    }
}

template<typename TF>
void Boundary_cyclic<TF>::exec(unsigned int* restrict data, Edge edge)
{
    auto& gd = grid.get_grid_data();

    const int jj = gd.icells;
    const int kk = gd.icells*gd.jcells;

    if (edge == Edge::East_west_edge || edge == Edge::Both_edges)
    {
        // first, east west boundaries
        for (int k=0; k<gd.kcells; ++k)
            for (int j=0; j<gd.jcells; ++j)
                #pragma ivdep
                for (int i=0; i<gd.igc; ++i)
                {
                    const int ijk0 = i          + j*jj + k*kk;
                    const int ijk1 = gd.iend-gd.igc+i + j*jj + k*kk;
                    data[ijk0] = data[ijk1];
                }

        for (int k=0; k<gd.kcells; ++k)
            for (int j=0; j<gd.jcells; ++j)
                #pragma ivdep
                for (int i=0; i<gd.igc; ++i)
                {
                    const int ijk0 = i+gd.iend   + j*jj + k*kk;
                    const int ijk1 = i+gd.istart + j*jj + k*kk;
                    data[ijk0] = data[ijk1];
                }
    }

    if (edge == Edge::North_south_edge || edge == Edge::Both_edges)
    {
        // if the run is 3D, apply the BCs
        if (gd.jtot > 1)
        {
            // second, send and receive the ghost cells in the north-south direction
            for (int k=0; k<gd.kcells; ++k)
                for (int j=0; j<gd.jgc; ++j)
                    #pragma ivdep
                    for (int i=0; i<gd.icells; ++i)
                    {
                        const int ijk0 = i + j                 *jj + k*kk;
                        const int ijk1 = i + (gd.jend-gd.jgc+j)*jj + k*kk;
                        data[ijk0] = data[ijk1];
                    }

            for (int k=0; k<gd.kcells; ++k)
                for (int j=0; j<gd.jgc; ++j)
                    #pragma ivdep
                    for (int i=0; i<gd.icells; ++i)
                    {
                        const int ijk0 = i + (j+gd.jend  )*jj + k*kk;
                        const int ijk1 = i + (j+gd.jstart)*jj + k*kk;
                        data[ijk0] = data[ijk1];
                    }
        }
        // in case of 2D, fill all the ghost cells with the current value
        else
        {
            for (int k=gd.kstart; k<gd.kend; ++k)
                for (int j=0; j<gd.jgc; ++j)
                    #pragma ivdep
                    for (int i=0; i<gd.icells; ++i)
                    {
                        const int ijkref   = i + gd.jstart*jj   + k*kk;
                        const int ijknorth = i + j*jj           + k*kk;
                        const int ijksouth = i + (gd.jend+j)*jj + k*kk;
                        data[ijknorth] = data[ijkref];
                        data[ijksouth] = data[ijkref];
                    }
        }
    }
}

template<typename TF>
void Boundary_cyclic<TF>::exec_2d(unsigned int* restrict data)
{
    auto& gd = grid.get_grid_data();

    const int jj = gd.icells;

    // First, east west boundaries.
    for (int j=0; j<gd.jcells; ++j)
        #pragma ivdep
        for (int i=0; i<gd.igc; ++i)
        {
            const int ij0 = i                + j*jj;
            const int ij1 = gd.iend-gd.igc+i + j*jj;
            data[ij0] = data[ij1];
        }

    for (int j=0; j<gd.jcells; ++j)
        #pragma ivdep
        for (int i=0; i<gd.igc; ++i)
        {
            const int ij0 = i+gd.iend   + j*jj;
            const int ij1 = i+gd.istart + j*jj;
            data[ij0] = data[ij1];
        }

    // If the run is 3D, apply the BCs.
    if (gd.jtot > 1)
    {
        // Second, send and receive the ghost cells in the north-south direction.
        for (int j=0; j<gd.jgc; ++j)
            #pragma ivdep
            for (int i=0; i<gd.icells; ++i)
            {
                const int ij0 = i + j                 *jj;
                const int ij1 = i + (gd.jend-gd.jgc+j)*jj;
                data[ij0] = data[ij1];
            }

        for (int j=0; j<gd.jgc; ++j)
            #pragma ivdep
            for (int i=0; i<gd.icells; ++i)
            {
                const int ij0 = i + (j+gd.jend  )*jj;
                const int ij1 = i + (j+gd.jstart)*jj;
                data[ij0] = data[ij1];
            }
    }
    // In case of 2D, fill all the ghost cells with the current value.
    else
    {
        for (int j=0; j<gd.jgc; ++j)
            #pragma ivdep
            for (int i=0; i<gd.icells; ++i)
            {
                const int ijref   = i + gd.jstart*jj;
                const int ijnorth = i + j*jj;
                const int ijsouth = i + (gd.jend+j)*jj;
                data[ijnorth] = data[ijref];
                data[ijsouth] = data[ijref];
            }
    }
}
#endif


#ifdef FLOAT_SINGLE
template class Boundary_cyclic<float>;
#else
template class Boundary_cyclic<double>;
#endif
