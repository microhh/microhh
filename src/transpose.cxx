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


template<typename TF>
Transpose<TF>::Transpose(Master& masterin, Grid<TF>& gridin) :
    master(masterin),
    grid(gridin),
    mpi_types_allocated(false)
{
}


template<typename TF>
Transpose<TF>::~Transpose()
{
    exit_mpi();
}


template<typename TF>
void Transpose<TF>::init()
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
void Transpose<TF>::init_mpi()
{
    auto& gd = grid.get_grid_data();

    int datacount, datablock, datastride;

    // transposez
    datacount = gd.imax*gd.jmax*gd.kblock;
    MPI_Type_contiguous(datacount, mpi_fp_type<TF>(), &transposez);
    MPI_Type_commit(&transposez);

    // transposez iblock/jblock/kblock
    datacount = gd.iblock*gd.jblock*gd.kblock;
    MPI_Type_contiguous(datacount, mpi_fp_type<TF>(), &transposez2);
    MPI_Type_commit(&transposez2);

    // transposex imax
    datacount  = gd.jmax*gd.kblock;
    datablock  = gd.imax;
    datastride = gd.itot;
    MPI_Type_vector(datacount, datablock, datastride, mpi_fp_type<TF>(), &transposex);
    MPI_Type_commit(&transposex);

    // transposex iblock
    datacount  = gd.jmax*gd.kblock;
    datablock  = gd.iblock;
    datastride = gd.itot;
    MPI_Type_vector(datacount, datablock, datastride, mpi_fp_type<TF>(), &transposex2);
    MPI_Type_commit(&transposex2);

    // transposey
    datacount  = gd.kblock;
    datablock  = gd.iblock*gd.jmax;
    datastride = gd.iblock*gd.jtot;
    MPI_Type_vector(datacount, datablock, datastride, mpi_fp_type<TF>(), &transposey);
    MPI_Type_commit(&transposey);

    // transposey2
    datacount  = gd.kblock;
    datablock  = gd.iblock*gd.jblock;
    datastride = gd.iblock*gd.jtot;
    MPI_Type_vector(datacount, datablock, datastride, mpi_fp_type<TF>(), &transposey2);
    MPI_Type_commit(&transposey2);

    mpi_types_allocated = true;
}


template<typename TF>
void Transpose<TF>::exit_mpi()
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


template<typename TF>
void Transpose<TF>::exec_zx(TF* const restrict data, TF* const restrict buffer_send, TF* const restrict buffer_recv)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    if (md.npx == 1)
        return;

    // Compute the appropriate strides.
    const int jj = gd.imax;
    const int kk = gd.imax*gd.jmax;
    const int nn = gd.imax*gd.jmax*gd.kblock;

    // Local copies to make OpenACC work.
    const int npx = md.npx;
    const int kblock = gd.kblock;
    const int jmax = gd.jmax;
    const int imax = gd.imax;

    // Pack the buffers.
    // #pragma acc parallel loop present(buffer_send, data) collapse(4)
    for (int n=0; n<npx; ++n)
        for (int k=0; k<kblock; ++k)
            for (int j=0; j<jmax; ++j)
                for (int i=0; i<imax; ++i)
                {
                    const int ijk = i + j*jj + k*kk + n*nn;
                    const int ijk_buf = i + j*jj + k*kk + n*nn;
                    buffer_send[ijk_buf] = data[ijk];
                }

    // Send and receive the buffers.
    const int tag = 1;
    const int nreqs = md.npx*2;
    MPI_Request reqs[nreqs];

    // #pragma acc host_data use_device(buffer_send, buffer_recv)
    {
        for (int n=0; n<md.npx; ++n)
        {
            MPI_Isend(&buffer_send[n*nn], nn, mpi_fp_type<TF>(), n, tag, md.commx, &reqs[2*n  ]);
            MPI_Irecv(&buffer_recv[n*nn], nn, mpi_fp_type<TF>(), n, tag, md.commx, &reqs[2*n+1]);
        }
        MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
    }

    // Unpack the buffer.
    const int jj_x = gd.itot;
    const int kk_x = gd.itot*gd.jmax;

    // #pragma acc parallel loop present(buffer_recv, data) collapse(4)
    for (int n=0; n<npx; ++n)
        for (int k=0; k<kblock; ++k)
            for (int j=0; j<jmax; ++j)
                for (int i=0; i<imax; ++i)
                {
                    const int ijk = (i + n*imax) + j*jj_x + k*kk_x;
                    const int ijk_buf = i + j*jj + k*kk + n*nn;
                    data[ijk] = buffer_recv[ijk_buf];
                }
}


template<typename TF>
void Transpose<TF>::exec_xz(TF* const restrict ar, TF* const restrict as)
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

template<typename TF>
void Transpose<TF>::exec_xy(TF* const restrict ar, TF* const restrict as)
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

template<typename TF>
void Transpose<TF>::exec_yx(TF* const restrict ar, TF* const restrict as)
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

template<typename TF>
void Transpose<TF>::exec_yz(TF* const restrict ar, TF* const restrict as)
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

template<typename TF>
void Transpose<TF>::exec_zy(TF* const restrict ar, TF* const restrict as)
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

template<typename TF>
void Transpose<TF>::init_mpi()
{
}

template<typename TF>
void Transpose<TF>::exit_mpi()
{
}
#endif


#ifdef FLOAT_SINGLE
template class Transpose<float>;
#else
template class Transpose<double>;
#endif
