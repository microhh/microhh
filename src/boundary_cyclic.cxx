/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
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
    grid(gridin)
{
}

template<typename TF>
Boundary_cyclic<TF>::~Boundary_cyclic()
{
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

    // north south
    datacount  = gd.kcells;
    datablock  = gd.icells*gd.jgc;
    datastride = gd.icells*gd.jcells;
    MPI_Type_vector(datacount, datablock, datastride, mpi_fp_type<TF>(), &northsouthedge);
    MPI_Type_commit(&northsouthedge);
}

template<typename TF>
void Boundary_cyclic<TF>::exec(TF* const restrict data, Edge edge)
{
    auto& gd = grid.get_grid_data();

    const int ncount = 1;

    if (edge == Edge::East_west_edge || edge == Edge::Both_edges)
    {
        // Communicate east-west edges.
        const int eastout = gd.iend-gd.igc;
        const int westin  = 0;
        const int westout = gd.istart;
        const int eastin  = gd.iend;

        // Send and receive the ghost cells in east-west direction.
        MPI_Isend(&data[eastout], ncount, eastwestedge, master.neast, 1, master.commxy, &master.reqs[master.reqsn]);
        master.reqsn++;
        MPI_Irecv(&data[westin], ncount, eastwestedge, master.nwest, 1, master.commxy, &master.reqs[master.reqsn]);
        master.reqsn++;
        MPI_Isend(&data[westout], ncount, eastwestedge, master.nwest, 2, master.commxy, &master.reqs[master.reqsn]);
        master.reqsn++;
        MPI_Irecv(&data[eastin], ncount, eastwestedge, master.neast, 2, master.commxy, &master.reqs[master.reqsn]);
        master.reqsn++;
        // Wait here for the MPI to have correct values in the corners of the cells.
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
            MPI_Isend(&data[northout], ncount, northsouthedge, master.nnorth, 1, master.commxy, &master.reqs[master.reqsn]);
            master.reqsn++;
            MPI_Irecv(&data[southin], ncount, northsouthedge, master.nsouth, 1, master.commxy, &master.reqs[master.reqsn]);
            master.reqsn++;
            MPI_Isend(&data[southout], ncount, northsouthedge, master.nsouth, 2, master.commxy, &master.reqs[master.reqsn]);
            master.reqsn++;
            MPI_Irecv(&data[northin], ncount, northsouthedge, master.nnorth, 2, master.commxy, &master.reqs[master.reqsn]);
            master.reqsn++;
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

#else

template<typename TF>
void Boundary_cyclic<TF>::init_mpi()
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
        for (int k=0; k<gd.kcells; k++)
            for (int j=0; j<gd.jcells; j++)
                #pragma ivdep
                for (int i=0; i<gd.igc; i++)
                {
                    const int ijk0 = i          + j*jj + k*kk;
                    const int ijk1 = gd.iend-gd.igc+i + j*jj + k*kk;
                    data[ijk0] = data[ijk1];
                }

        for (int k=0; k<gd.kcells; k++)
            for (int j=0; j<gd.jcells; j++)
                #pragma ivdep
                for (int i=0; i<gd.igc; i++)
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
            for (int k=0; k<gd.kcells; k++)
                for (int j=0; j<gd.jgc; j++)
                    #pragma ivdep
                    for (int i=0; i<gd.icells; i++)
                    {
                        const int ijk0 = i + j                 *jj + k*kk;
                        const int ijk1 = i + (gd.jend-gd.jgc+j)*jj + k*kk;
                        data[ijk0] = data[ijk1];
                    }

            for (int k=0; k<gd.kcells; k++)
                for (int j=0; j<gd.jgc; j++)
                    #pragma ivdep
                    for (int i=0; i<gd.icells; i++)
                    {
                        const int ijk0 = i + (j+gd.jend  )*jj + k*kk;
                        const int ijk1 = i + (j+gd.jstart)*jj + k*kk;
                        data[ijk0] = data[ijk1];
                    }
        }
        // in case of 2D, fill all the ghost cells with the current value
        else
        {
            for (int k=gd.kstart; k<gd.kend; k++)
                for (int j=0; j<gd.jgc; j++)
                    #pragma ivdep
                    for (int i=0; i<gd.icells; i++)
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
#endif

template class Boundary_cyclic<double>;
template class Boundary_cyclic<float>;
