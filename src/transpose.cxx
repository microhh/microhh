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
#include "transpose.h"

template<typename TF, typename TF_data>
Transpose<TF, TF_data>::Transpose(Master& masterin, Grid<TF>& gridin) :
    master(masterin),
    grid(gridin),
    mpi_types_allocated(false)
{
}

template<typename TF, typename TF_data>
Transpose<TF, TF_data>::~Transpose()
{
    exit_mpi();
}

template<typename TF, typename TF_data>
void Transpose<TF, TF_data>::init()
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

template<typename TF, typename TF_data>
void Transpose<TF, TF_data>::init_mpi()
{
    auto& gd = grid.get_grid_data();

    int datacount, datablock, datastride;

    // transposez
    datacount = gd.imax*gd.jmax*gd.kblock;
    MPI_Type_contiguous(datacount, mpi_fp_type<TF_data>(), &transposez);
    MPI_Type_commit(&transposez);

    // transposez iblock/jblock/kblock
    datacount = gd.iblock*gd.jblock*gd.kblock;
    MPI_Type_contiguous(datacount, mpi_fp_type<TF_data>(), &transposez2);
    MPI_Type_commit(&transposez2);

    // transposex imax
    datacount  = gd.jmax*gd.kblock;
    datablock  = gd.imax;
    datastride = gd.itot;
    MPI_Type_vector(datacount, datablock, datastride, mpi_fp_type<TF_data>(), &transposex);
    MPI_Type_commit(&transposex);

    // transposex iblock
    datacount  = gd.jmax*gd.kblock;
    datablock  = gd.iblock;
    datastride = gd.itot;
    MPI_Type_vector(datacount, datablock, datastride, mpi_fp_type<TF_data>(), &transposex2);
    MPI_Type_commit(&transposex2);

    // transposey
    datacount  = gd.kblock;
    datablock  = gd.iblock*gd.jmax;
    datastride = gd.iblock*gd.jtot;
    MPI_Type_vector(datacount, datablock, datastride, mpi_fp_type<TF_data>(), &transposey);
    MPI_Type_commit(&transposey);

    // transposey2
    datacount  = gd.kblock;
    datablock  = gd.iblock*gd.jblock;
    datastride = gd.iblock*gd.jtot;
    MPI_Type_vector(datacount, datablock, datastride, mpi_fp_type<TF_data>(), &transposey2);
    MPI_Type_commit(&transposey2);

    mpi_types_allocated = true;
}

template<typename TF, typename TF_data>
void Transpose<TF, TF_data>::exit_mpi()
{
    if (mpi_types_allocated)
    {
        MPI_Type_free(&transposez);
        MPI_Type_free(&transposez2);
        MPI_Type_free(&transposex);
        MPI_Type_free(&transposex2);
        MPI_Type_free(&transposey);
        MPI_Type_free(&transposey2);
    }
}

template<typename TF, typename TF_data>
void Transpose<TF, TF_data>::exec_zx(TF_data* const restrict ar, TF_data* const restrict as)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    const int ncount = 1;
    const int tag = 1;

    const int jj = gd.imax;
    const int kk = gd.imax*gd.jmax;

    for (int n=0; n<md.npx; ++n)
    {
        // Determine where to fetch the data and where to store it.
        const int ijks = n*gd.kblock*kk;
        const int ijkr = n*jj;

        // Send and receive the data.
        MPI_Isend(&as[ijks], ncount, transposez, n, tag, md.commx, master.get_request_ptr());
        MPI_Irecv(&ar[ijkr], ncount, transposex, n, tag, md.commx, master.get_request_ptr());
    }

    master.wait_all();
}

template<typename TF, typename TF_data>
void Transpose<TF, TF_data>::exec_xz(TF_data* const restrict ar, TF_data* const restrict as)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    const int ncount = 1;
    const int tag = 1;

    const int jj = gd.imax;
    const int kk = gd.imax*gd.jmax;

    for (int n=0; n<md.npx; ++n)
    {
        // Determine where to fetch the data and where to store it.
        const int ijks = n*jj;
        const int ijkr = n*gd.kblock*kk;

        // Send and receive the data.
        MPI_Isend(&as[ijks], ncount, transposex, n, tag, md.commx, master.get_request_ptr());
        MPI_Irecv(&ar[ijkr], ncount, transposez, n, tag, md.commx, master.get_request_ptr());
    }

    master.wait_all();
}

template<typename TF, typename TF_data>
void Transpose<TF, TF_data>::exec_xy(TF_data* const restrict ar, TF_data* const restrict as)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    const int ncount = 1;
    const int tag = 1;

    const int jj = gd.iblock;
    const int kk = gd.iblock*gd.jmax;

    for (int n=0; n<md.npy; ++n)
    {
        // determine where to fetch the data and where to store it
        const int ijks = n*jj;
        const int ijkr = n*kk;

        // send and receive the data
        MPI_Isend(&as[ijks], ncount, transposex2, n, tag, md.commy, master.get_request_ptr());
        MPI_Irecv(&ar[ijkr], ncount, transposey , n, tag, md.commy, master.get_request_ptr());
    }

    master.wait_all();
}

template<typename TF, typename TF_data>
void Transpose<TF, TF_data>::exec_yx(TF_data* const restrict ar, TF_data* const restrict as)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    const int ncount = 1;
    const int tag = 1;

    const int jj = gd.iblock;
    const int kk = gd.iblock*gd.jmax;

    for (int n=0; n<md.npy; ++n)
    {
        // determine where to fetch the data and where to store it
        const int ijks = n*kk;
        const int ijkr = n*jj;

        // send and receive the data
        MPI_Isend(&as[ijks], ncount, transposey , n, tag, md.commy, master.get_request_ptr());
        MPI_Irecv(&ar[ijkr], ncount, transposex2, n, tag, md.commy, master.get_request_ptr());
    }

    master.wait_all();
}

template<typename TF, typename TF_data>
void Transpose<TF, TF_data>::exec_yz(TF_data* const restrict ar, TF_data* const restrict as)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    const int ncount = 1;
    const int tag = 1;

    const int jj = gd.iblock;
    const int kk = gd.iblock*gd.jblock;

    for (int n=0; n<md.npx; ++n)
    {
        // determine where to fetch the data and where to store it
        const int ijks = n*gd.jblock*jj;
        const int ijkr = n*gd.kblock*kk;

        // send and receive the data
        MPI_Isend(&as[ijks], ncount, transposey2, n, tag, md.commx, master.get_request_ptr());
        MPI_Irecv(&ar[ijkr], ncount, transposez2, n, tag, md.commx, master.get_request_ptr());
    }

    master.wait_all();
}

template<typename TF, typename TF_data>
void Transpose<TF, TF_data>::exec_zy(TF_data* const restrict ar, TF_data* const restrict as)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    const int ncount = 1;
    const int tag = 1;

    const int jj = gd.iblock;
    const int kk = gd.iblock*gd.jblock;

    for (int n=0; n<md.npx; ++n)
    {
        // determine where to fetch the data and where to store it
        const int ijks = n*gd.kblock*kk;
        const int ijkr = n*gd.jblock*jj;

        // send and receive the data
        MPI_Isend(&as[ijks], ncount, transposez2, n, tag, md.commx, master.get_request_ptr());
        MPI_Irecv(&ar[ijkr], ncount, transposey2, n, tag, md.commx, master.get_request_ptr());
    }

    master.wait_all();
}
#else

template<typename TF, typename TF_data>
void Transpose<TF, TF_data>::init_mpi()
{
}

template<typename TF, typename TF_data>
void Transpose<TF, TF_data>::exit_mpi()
{
}
#endif


#ifdef FLOAT_SINGLE
template class Transpose<float, float>;
#else
template class Transpose<double, double>;
#endif
