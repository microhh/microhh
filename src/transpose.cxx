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
    grid(gridin)
{
}


template<typename TF>
Transpose<TF>::~Transpose()
{
}


#ifdef USEMPI
namespace
{
    template<typename TF> MPI_Datatype mpi_fp_type();
    template<> MPI_Datatype mpi_fp_type<double>() { return MPI_DOUBLE; }
    template<> MPI_Datatype mpi_fp_type<float>() { return MPI_FLOAT; }
}


template<typename TF>
void Transpose<TF>::exec_zx(TF* const restrict data)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    if (md.npx == 1)
        return;

    TF* restrict buffer_send = grid.get_tmp_3d();
    TF* restrict buffer_recv = grid.get_tmp_3d();

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
    #pragma acc parallel loop present(buffer_send, data) collapse(4)
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

    #pragma acc host_data use_device(buffer_send, buffer_recv)
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

    #pragma acc parallel loop present(buffer_recv, data) collapse(4)
    for (int n=0; n<npx; ++n)
        for (int k=0; k<kblock; ++k)
            for (int j=0; j<jmax; ++j)
                for (int i=0; i<imax; ++i)
                {
                    const int ijk = (i + n*imax) + j*jj_x + k*kk_x;
                    const int ijk_buf = i + j*jj + k*kk + n*nn;
                    data[ijk] = buffer_recv[ijk_buf];
                }

    grid.release_tmp_3d(buffer_send);
    grid.release_tmp_3d(buffer_recv);
}


template<typename TF>
void Transpose<TF>::exec_xz(TF* const restrict data)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    if (md.npx == 1)
        return;

    TF* restrict buffer_send = grid.get_tmp_3d();
    TF* restrict buffer_recv = grid.get_tmp_3d();

    const int jj = gd.imax;
    const int kk = gd.imax*gd.jmax;
    const int nn = gd.imax*gd.jmax*gd.kblock;

    // Pack the buffer.
    const int jj_x = gd.itot;
    const int kk_x = gd.itot*gd.jmax;

    // Local copies to make OpenACC work.
    const int npx = md.npx;
    const int kblock = gd.kblock;
    const int jmax = gd.jmax;
    const int imax = gd.imax;

    #pragma acc parallel loop present(buffer_send, data) collapse(4)
    for (int n=0; n<npx; ++n)
        for (int k=0; k<kblock; ++k)
            for (int j=0; j<jmax; ++j)
                for (int i=0; i<imax; ++i)
                {
                    const int ijk = (i + n*imax) + j*jj_x + k*kk_x;
                    const int ijk_buf = i + j*jj + k*kk + n*nn;
                    buffer_send[ijk_buf] = data[ijk];
                }

    // Send and receive the buffers.
    const int tag = 1;
    const int nreqs = md.npx*2;
    MPI_Request reqs[nreqs];
    #pragma acc host_data use_device(buffer_send, buffer_recv)
    {
        for (int n=0; n<md.npx; ++n)
        {
            MPI_Isend(&buffer_send[n*nn], nn, mpi_fp_type<TF>(), n, tag, md.commx, &reqs[2*n  ]);
            MPI_Irecv(&buffer_recv[n*nn], nn, mpi_fp_type<TF>(), n, tag, md.commx, &reqs[2*n+1]);
        }
        MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
    }

    // Unpack the buffers.
    #pragma acc parallel loop present(buffer_recv, data) collapse(4)
    for (int n=0; n<npx; ++n)
        for (int k=0; k<kblock; ++k)
            for (int j=0; j<jmax; ++j)
                for (int i=0; i<imax; ++i)
                {
                    const int ijk = i + j*jj + k*kk + n*nn;
                    const int ijk_buf = i + j*jj + k*kk + n*nn;
                    data[ijk] = buffer_recv[ijk_buf];
                }

    grid.release_tmp_3d(buffer_send);
    grid.release_tmp_3d(buffer_recv);
}


template<typename TF>
void Transpose<TF>::exec_xy(TF* const restrict data)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    if (md.npy == 1)
        return;

    TF* restrict buffer_send = grid.get_tmp_3d();
    TF* restrict buffer_recv = grid.get_tmp_3d();

    // Compute the appropriate strides.
    const int jj_x = gd.itot;
    const int kk_x = gd.itot*gd.jmax;

    const int jj_y = gd.iblock;
    const int kk_y = gd.iblock*gd.jtot;

    const int jj = gd.iblock;
    const int kk = gd.iblock*gd.jmax;
    const int nn = gd.iblock*gd.jmax*gd.kblock;

    // Local copies to make OpenACC work.
    const int npy = md.npy;
    const int kblock = gd.kblock;
    const int jmax = gd.jmax;
    const int iblock = gd.iblock;

    // Pack the buffers.
    #pragma acc parallel loop present(buffer_send, data) collapse(4)
    for (int n=0; n<npy; ++n)
        for (int k=0; k<kblock; ++k)
            for (int j=0; j<jmax; ++j)
                for (int i=0; i<iblock; ++i)
                {
                    const int ijk = (i + n*iblock) + j*jj_x + k*kk_x;
                    const int ijk_buf = i + j*jj + k*kk + n*nn;
                    buffer_send[ijk_buf] = data[ijk];
                }

    // Send and receive the buffers.
    const int tag = 1;
    const int nreqs = md.npy*2;
    MPI_Request reqs[nreqs];
    #pragma acc host_data use_device(buffer_send, buffer_recv)
    {
        for (int n=0; n<md.npy; ++n)
        {
            MPI_Isend(&buffer_send[n*nn], nn, mpi_fp_type<TF>(), n, tag, md.commy, &reqs[2*n  ]);
            MPI_Irecv(&buffer_recv[n*nn], nn, mpi_fp_type<TF>(), n, tag, md.commy, &reqs[2*n+1]);
        }
        MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
    }

    // Unpack the buffer.
    #pragma acc parallel loop present(buffer_recv, data) collapse(4)
    for (int n=0; n<npy; ++n)
        for (int k=0; k<kblock; ++k)
            for (int j=0; j<jmax; ++j)
                for (int i=0; i<iblock; ++i)
                {
                    const int ijk = i + (j + n*jmax)*jj_y + k*kk_y;
                    const int ijk_buf = i + j*jj + k*kk + n*nn;
                    data[ijk] = buffer_recv[ijk_buf];
                }

    grid.release_tmp_3d(buffer_send);
    grid.release_tmp_3d(buffer_recv);
}


template<typename TF>
void Transpose<TF>::exec_yx(TF* const restrict data)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    if (md.npy == 1)
        return;

    TF* restrict buffer_send = grid.get_tmp_3d();
    TF* restrict buffer_recv = grid.get_tmp_3d();

    // Compute the appropriate strides.
    const int jj_y = gd.iblock;
    const int kk_y = gd.iblock*gd.jtot;

    const int jj_x = gd.itot;
    const int kk_x = gd.itot*gd.jmax;

    const int jj = gd.iblock;
    const int kk = gd.iblock*gd.jmax;
    const int nn = gd.iblock*gd.jmax*gd.kblock;

    // Local copies to make OpenACC work.
    const int npy = md.npy;
    const int kblock = gd.kblock;
    const int jmax = gd.jmax;
    const int iblock = gd.iblock;

    // Pack the buffer.
    #pragma acc parallel loop present(buffer_send, data) collapse(4)
    for (int n=0; n<npy; ++n)
        for (int k=0; k<kblock; ++k)
            for (int j=0; j<jmax; ++j)
                for (int i=0; i<iblock; ++i)
                {
                    const int ijk = i + (j + n*jmax)*jj_y + k*kk_y;
                    const int ijk_buf = i + j*jj + k*kk + n*nn;
                    buffer_send[ijk_buf] = data[ijk];
                }

    // Send and receive the buffers.
    const int tag = 1;
    const int nreqs = md.npy*2;
    MPI_Request reqs[nreqs];
    #pragma acc host_data use_device(buffer_send, buffer_recv)
    {
        for (int n=0; n<md.npy; ++n)
        {
            MPI_Isend(&buffer_send[n*nn], nn, mpi_fp_type<TF>(), n, tag, md.commy, &reqs[2*n  ]);
            MPI_Irecv(&buffer_recv[n*nn], nn, mpi_fp_type<TF>(), n, tag, md.commy, &reqs[2*n+1]);
        }
        MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
    }

    // Unpack the buffers.
    #pragma acc parallel loop present(buffer_recv, data) collapse(4)
    for (int n=0; n<npy; ++n)
        for (int k=0; k<kblock; ++k)
            for (int j=0; j<jmax; ++j)
                for (int i=0; i<iblock; ++i)
                {
                    const int ijk = (i + n*iblock) + j*jj_x + k*kk_x;
                    const int ijk_buf = i + j*jj + k*kk + n*nn;
                    data[ijk] = buffer_recv[ijk_buf];
                }

    grid.release_tmp_3d(buffer_send);
    grid.release_tmp_3d(buffer_recv);
}


template<typename TF>
void Transpose<TF>::exec_yz(TF* const restrict data)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    if (md.npx == 1)
        return;

    TF* restrict buffer_send = grid.get_tmp_3d();
    TF* restrict buffer_recv = grid.get_tmp_3d();

    // Compute the appropriate strides.
    const int jj_y = gd.iblock;
    const int kk_y = gd.iblock*gd.jtot;

    const int jj_z = gd.iblock;
    const int kk_z = gd.iblock*gd.jblock;

    const int jj = gd.iblock;
    const int kk = gd.iblock*gd.jblock;
    const int nn = gd.iblock*gd.jblock*gd.kblock;

    // Local copies to make OpenACC work.
    const int npx = md.npx;
    const int kblock = gd.kblock;
    const int jblock = gd.jblock;
    const int iblock = gd.iblock;

    // Pack the buffers.
    #pragma acc parallel loop present(buffer_send, data) collapse(4)
    for (int n=0; n<npx; ++n)
        for (int k=0; k<kblock; ++k)
            for (int j=0; j<jblock; ++j)
                for (int i=0; i<iblock; ++i)
                {
                    const int ijk = i + (j + n*jblock)*jj_y + k*kk_y;
                    const int ijk_buf = i + j*jj + k*kk + n*nn;
                    buffer_send[ijk_buf] = data[ijk];
                }

    // Send and receive the buffers.
    const int tag = 1;
    const int nreqs = md.npx*2;
    MPI_Request reqs[nreqs];
    #pragma acc host_data use_device(buffer_send, buffer_recv)
    {
        for (int n=0; n<md.npx; ++n)
        {
            MPI_Isend(&buffer_send[n*nn], nn, mpi_fp_type<TF>(), n, tag, md.commx, &reqs[2*n  ]);
            MPI_Irecv(&buffer_recv[n*nn], nn, mpi_fp_type<TF>(), n, tag, md.commx, &reqs[2*n+1]);
        }
        MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
    }

    // Unpack the buffer.
    #pragma acc parallel loop present(buffer_recv, data) collapse(4)
    for (int n=0; n<npx; ++n)
        for (int k=0; k<kblock; ++k)
            for (int j=0; j<jblock; ++j)
                for (int i=0; i<iblock; ++i)
                {
                    const int ijk = i + j*jj_z + (k + n*kblock)*kk_z;
                    const int ijk_buf = i + j*jj + k*kk + n*nn;
                    data[ijk] = buffer_recv[ijk_buf];
                }

    grid.release_tmp_3d(buffer_send);
    grid.release_tmp_3d(buffer_recv);
}


template<typename TF>
void Transpose<TF>::exec_zy(TF* const restrict data)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    if (md.npx == 1)
        return;

    TF* restrict buffer_send = grid.get_tmp_3d();
    TF* restrict buffer_recv = grid.get_tmp_3d();

    // Compute the appropriate strides.
    const int jj_y = gd.iblock;
    const int kk_y = gd.iblock*gd.jtot;

    const int jj_z = gd.iblock;
    const int kk_z = gd.iblock*gd.jblock;

    const int jj = gd.iblock;
    const int kk = gd.iblock*gd.jblock;
    const int nn = gd.iblock*gd.jblock*gd.kblock;

    // Local copies to make OpenACC work.
    const int npx = md.npx;
    const int kblock = gd.kblock;
    const int jblock = gd.jblock;
    const int iblock = gd.iblock;

    // Pack the buffer.
    #pragma acc parallel loop present(buffer_send, data) collapse(4)
    for (int n=0; n<npx; ++n)
        for (int k=0; k<kblock; ++k)
            for (int j=0; j<jblock; ++j)
                for (int i=0; i<iblock; ++i)
                {
                    const int ijk = i + j*jj_z + (k + n*kblock)*kk_z;
                    const int ijk_buf = i + j*jj + k*kk + n*nn;
                    buffer_send[ijk_buf] = data[ijk];
                }

    // Send and receive the buffers.
    const int tag = 1;
    const int nreqs = md.npx*2;
    MPI_Request reqs[nreqs];
    #pragma acc host_data use_device(buffer_send, buffer_recv)
    {
        for (int n=0; n<md.npx; ++n)
        {
            MPI_Isend(&buffer_send[n*nn], nn, mpi_fp_type<TF>(), n, tag, md.commx, &reqs[2*n  ]);
            MPI_Irecv(&buffer_recv[n*nn], nn, mpi_fp_type<TF>(), n, tag, md.commx, &reqs[2*n+1]);
        }
        MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
    }

    // Unpack the buffers.
    #pragma acc parallel loop present(buffer_recv, data) collapse(4)
    for (int n=0; n<npx; ++n)
        for (int k=0; k<kblock; ++k)
            for (int j=0; j<jblock; ++j)
                for (int i=0; i<iblock; ++i)
                {
                    const int ijk = i + (j + n*jblock)*jj_y + k*kk_y;
                    const int ijk_buf = i + j*jj + k*kk + n*nn;
                    data[ijk] = buffer_recv[ijk_buf];
                }

    grid.release_tmp_3d(buffer_send);
    grid.release_tmp_3d(buffer_recv);
}

#else
template<typename TF> void Transpose<TF>::exec_zx(TF* const restrict data, TF* const restrict buffer_send, TF* const restrict buffer_recv) {}
template<typename TF> void Transpose<TF>::exec_xz(TF* const restrict data, TF* const restrict buffer_send, TF* const restrict buffer_recv) {}
template<typename TF> void Transpose<TF>::exec_xy(TF* const restrict data, TF* const restrict buffer_send, TF* const restrict buffer_recv) {}
template<typename TF> void Transpose<TF>::exec_yx(TF* const restrict data, TF* const restrict buffer_send, TF* const restrict buffer_recv) {}
template<typename TF> void Transpose<TF>::exec_yz(TF* const restrict data, TF* const restrict buffer_send, TF* const restrict buffer_recv) {}
template<typename TF> void Transpose<TF>::exec_zy(TF* const restrict data, TF* const restrict buffer_send, TF* const restrict buffer_recv) {}
#endif


#ifdef FLOAT_SINGLE
template class Transpose<float>;
#else
template class Transpose<double>;
#endif
