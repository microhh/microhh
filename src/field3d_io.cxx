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
#include <iostream>
#include <cmath>
#include <chrono>   // tmp
#include "master.h"
#include "grid.h"
#include "field3d.h"
#include "fields.h"
#include "defines.h"
#include "field3d_io.h"

template<typename TF>
Field3d_io<TF>::Field3d_io(Master& masterin, Grid<TF>& gridin) :
    master(masterin), grid(gridin)
{
}

template<typename TF>
Field3d_io<TF>::~Field3d_io()
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
int Field3d_io<TF>::save_field3d(
        TF* const restrict data,
        TF* const restrict tmp1, TF* const restrict tmp2,
        const char* filename, const TF offset,
        const int kstart, const int kend)
{
    // Save the data in transposed order to have large chunks of contiguous disk space.
    // MPI-IO is not stable on Juqueen and Supermuc otherwise
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();
    auto tp = Transpose<TF>(master, grid);
    tp.init();

    // Extract the data from the 3d field without the ghost cells
    const int jj    = gd.icells;
    const int kk    = gd.icells*gd.jcells;
    const int jjb   = gd.imax;
    const int kkb   = gd.imax*gd.jmax;
    const int kmax  = kend-kstart;
    const int count = gd.imax*gd.jmax*kmax;

    MPI_Datatype subarray;   // MPI datatype containing the dimensions of the total array that is contained in one process.

    // For full 3D fields, use the transposed save to increase IO performance
    bool sw_transpose = (kmax == gd.kmax) ? true : false;

    for (int k=0; k<kmax; ++k)
        for (int j=0; j<gd.jmax; ++j)
            #pragma ivdep
            for (int i=0; i<gd.imax; ++i)
            {
                const int ijk  = i+gd.igc + (j+gd.jgc)*jj + (k+kstart)*kk;
                const int ijkb = i + j*jjb + k*kkb;
                tmp1[ijkb] = data[ijk] + offset;
            }

    if (sw_transpose)
    {
        // Transpose the 3D field
        tp.exec_zx(tmp2, tmp1);

        // Create MPI datatype for writing transposed field
        int totsize [3] = {gd.kmax,   gd.jtot, gd.itot};
        int subsize [3] = {gd.kblock, gd.jmax, gd.itot};
        int substart[3] = {md.mpicoordx*gd.kblock, md.mpicoordy*gd.jmax, 0};
        MPI_Type_create_subarray(3, totsize, subsize, substart, MPI_ORDER_C, mpi_fp_type<TF>(), &subarray);
        MPI_Type_commit(&subarray);
    }
    else
    {
        // Create MPI datatype for writing non-transposed field
        int totsize [3] = {kmax, gd.jtot, gd.itot};
        int subsize [3] = {kmax, gd.jmax, gd.imax};
        int substart[3] = {0, md.mpicoordy*gd.jmax, md.mpicoordx*gd.imax};
        MPI_Type_create_subarray(3, totsize, subsize, substart, MPI_ORDER_C, mpi_fp_type<TF>(), &subarray);
        MPI_Type_commit(&subarray);
    }

    MPI_File fh;
    if (MPI_File_open(md.commxy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
        return 1;

    // select noncontiguous part of 3d array to store the selected data
    MPI_Offset fileoff = 0; // the offset within the file (header size)
    char name[] = "native";

    if (MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subarray, name, MPI_INFO_NULL))
        return 1;

    if (sw_transpose)
    {
        if (MPI_File_write_all(fh, tmp2, count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE))
            return 1;
    }
    else
    {
        if (MPI_File_write_all(fh, tmp1, count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE))
            return 1;
    }

    if (MPI_File_close(&fh))
        return 1;

    MPI_Type_free(&subarray);

    return 0;
}

template<typename TF>
int Field3d_io<TF>::load_field3d(
        TF* const restrict data,
        TF* const restrict tmp1, TF* const restrict tmp2,
        const char* filename, TF offset,
        const int kstart, const int kend)
{
    // Read the data (optionally) in transposed order to have large chunks of contiguous disk space.
    // MPI-IO is not stable on Juqueen and supermuc otherwise.
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();
    auto tp = Transpose<TF>(master, grid);
    tp.init();

    MPI_Datatype subarray;   // MPI datatype containing the dimensions of the total array that is contained in one process.

    const int kmax = kend-kstart;

    // For full 3D fields, use the transposed read to increase IO performance
    bool sw_transpose = (kmax == gd.kmax) ? true : false;

    if (sw_transpose)
    {
        // Create MPI datatype for reading transposed field
        int totsize [3] = {gd.kmax  , gd.jtot, gd.itot};
        int subsize [3] = {gd.kblock, gd.jmax, gd.itot};
        int substart[3] = {md.mpicoordx*gd.kblock, md.mpicoordy*gd.jmax, 0};
        MPI_Type_create_subarray(3, totsize, subsize, substart, MPI_ORDER_C, mpi_fp_type<TF>(), &subarray);
        MPI_Type_commit(&subarray);
    }
    else
    {
        // Create MPI datatype for reading non-transposed field
        int totsize [3] = {kmax, gd.jtot, gd.itot};
        int subsize [3] = {kmax, gd.jmax, gd.imax};
        int substart[3] = {0, md.mpicoordy*gd.jmax, md.mpicoordx*gd.imax};
        MPI_Type_create_subarray(3, totsize, subsize, substart, MPI_ORDER_C, mpi_fp_type<TF>(), &subarray);
        MPI_Type_commit(&subarray);
    }

    // Read the file
    MPI_File fh;
    if (MPI_File_open(md.commxy, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
        return 1;

    // select noncontiguous part of 3d array to store the selected data
    MPI_Offset fileoff = 0; // the offset within the file (header size)
    char name[] = "native";
    MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subarray, name, MPI_INFO_NULL);

    // extract the data from the 3d field without the ghost cells
    int count = gd.imax*gd.jmax*kmax;

    if (sw_transpose)
    {
        if (MPI_File_read_all(fh, tmp1, count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE))
            return 1;
    }
    else
    {
        if (MPI_File_read_all(fh, tmp2, count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE))
            return 1;
    }

    if (MPI_File_close(&fh))
        return 1;

    // Transpose the 3D field
    if (sw_transpose)
        tp.exec_xz(tmp2, tmp1);

    const int jj  = gd.icells;
    const int kk  = gd.icells*gd.jcells;
    const int jjb = gd.imax;
    const int kkb = gd.imax*gd.jmax;

    for (int k=0; k<kmax; ++k)
        for (int j=0; j<gd.jmax; ++j)
            #pragma ivdep
            for (int i=0; i<gd.imax; ++i)
            {
                const int ijk  = i+gd.igc + (j+gd.jgc)*jj + (k+kstart)*kk;
                const int ijkb = i + j*jjb + k*kkb;
                data[ijk] = tmp2[ijkb] - offset;
            }

    MPI_Type_free(&subarray);

    return 0;
}

template<typename TF>
int Field3d_io<TF>::save_xz_slice(
        TF* const restrict data, TF* const restrict tmp,
        const char* filename, const int jslice,
        const int kstart, const int kend)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    int nerror = 0;

    const int jj   = gd.icells;
    const int kk   = gd.icells*gd.jcells;
    const int kkb  = gd.imax;
    const int kmax = kend-kstart;

    int count = gd.imax*kmax;

    if (md.mpicoordy == jslice/gd.jmax)
    {
        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int i=0; i<gd.imax; i++)
            {
                // take the modulus of jslice and gd.jmax to have the right offset within proc
                const int ijk  = i+gd.igc + ((jslice%gd.jmax)+gd.jgc)*jj + k*kk;
                const int ijkb = i + (k-kstart)*kkb;
                tmp[ijkb] = data[ijk];
            }

#ifdef DISABLE_2D_MPIIO
//
// NOTE: on GPFS (IBM Spectrum Scale) file systems, saving small bits of data from
// each MPI task into a single binary using MPI-IO is **extremely** slow. With `DISABLE_2D_MPIIO`,
// the 2D slices are first gathered on MPI rank 0, and then written without MPI-IO.
//
        // MPI task which gathers/writes the slice:
        const int mpi_rank_recv = md.mpicoordy * md.npx;

        // Create send/receive MPI types.
        MPI_Datatype send_type;
        MPI_Type_vector(kmax, gd.imax, gd.imax, mpi_fp_type<TF>(), &send_type);
        MPI_Type_commit(&send_type);

        MPI_Datatype recv_type;
        int totxzsize_recv [2] = {kmax, gd.itot};
        int subxzsize_recv [2] = {kmax, gd.imax};
        int subxzstart_recv[2] = {0, md.mpicoordx*gd.imax};
        MPI_Type_create_subarray(2, totxzsize_recv, subxzsize_recv, subxzstart_recv, MPI_ORDER_C, mpi_fp_type<TF>(), &recv_type);
        MPI_Type_commit(&recv_type);

        MPI_Datatype recv_type_r;
        MPI_Type_create_resized(recv_type, 0, sizeof(TF), &recv_type_r);
        MPI_Type_commit(&recv_type_r);

        // TMP: benchmarking
        auto begin = std::chrono::high_resolution_clock::now();

        // Create size/offset arrays for MPI_Gatherv().
        std::vector<int> counts(md.npx);
        std::fill(counts.begin(), counts.end(), 1);

        std::vector<int> offset(md.npx);
        for (int i=0; i<md.npx; ++i)
            offset[i] = i*gd.imax;

        // Gather the data!
        std::vector<TF> recv;
        if (md.mpicoordx == mpi_rank_recv)
            recv.resize(gd.itot*gd.ktot);

        MPI_Gatherv(tmp, 1, send_type, recv.data(), counts.data(), offset.data(), recv_type_r, mpi_rank_recv, md.commx);

        // Only MPI rank 0 writes the data.
        if (md.mpicoordx == mpi_rank_recv)
        {
            FILE *pFile;
            pFile = fopen(filename, "wbx");
            if (pFile == NULL)
                return 1;

            fwrite(recv.data(), sizeof(TF), gd.itot*kmax, pFile);
            fclose(pFile);

            auto end = std::chrono::high_resolution_clock::now();
            auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
            const int bytes = gd.itot*kmax*sizeof(TF);
            const double Mbytes_s = double(bytes)/double(diff) * 1e9/1024/1024;
            std::cout << "Write XZ without MPI-IO = " << Mbytes_s << " MB/s" << std::endl;
        }
#else
//
// Use MPI-IO to write the data from different MPI tasks into a single file
//
        // Create MPI datatype for XZ-slice
        MPI_Datatype subxzslice;
        int totxzsize [2] = {kmax, gd.itot};
        int subxzsize [2] = {kmax, gd.imax};
        int subxzstart[2] = {0, md.mpicoordx*gd.imax};
        MPI_Type_create_subarray(2, totxzsize, subxzsize, subxzstart, MPI_ORDER_C, mpi_fp_type<TF>(), &subxzslice);
        MPI_Type_commit(&subxzslice);

        // TMP: benchmarking
        auto begin = std::chrono::high_resolution_clock::now();

        MPI_File fh;
        if (MPI_File_open(md.commx, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
            ++nerror;

        // select noncontiguous part of 3d array to store the selected data
        MPI_Offset fileoff = 0; // the offset within the file (header size)
        char name[] = "native";

        if (!nerror)
            if (MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subxzslice, name, MPI_INFO_NULL))
                ++nerror;

        // only write at the procs that contain the slice
        if (!nerror)
            if (MPI_File_write_all(fh, tmp, count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE))
                ++nerror;

        if (!nerror)
            MPI_File_sync(fh);

        if (!nerror)
            if (MPI_File_close(&fh))
                ++nerror;

        MPI_Type_free(&subxzslice);

        if (md.mpiid == 0)
        {
            auto end = std::chrono::high_resolution_clock::now();
            auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
            const int bytes = gd.itot*kmax*sizeof(TF);
            const double Mbytes_s = double(bytes)/double(diff) * 1e9/1024/1024;
            std::cout << "Write XZ with MPI-IO = " << Mbytes_s << " MB/s" << std::endl;
        }
#endif
    }

    // Gather errors from other processes
    master.sum(&nerror,1);
    MPI_Barrier(md.commxy);

    return nerror;
}

template<typename TF>
int Field3d_io<TF>::save_yz_slice(
        TF* restrict data, TF* restrict tmp,
        const char* filename, const int islice,
        const int kstart, const int kend)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    int nerror = 0;

    const int jj   = gd.icells;
    const int kk   = gd.ijcells;
    const int kkb  = gd.jmax;
    const int kmax = kend-kstart;

    int count = gd.jmax*kmax;

    if (md.mpicoordx == islice/gd.imax)
    {
        // Strip off the ghost cells
        for (int k=kstart; k<kend; k++)
            #pragma ivdep
            for (int j=0; j<gd.jmax; j++)
            {
                // take the modulus of jslice and jmax to have the right offset within proc
                const int ijk  = (islice%gd.imax)+gd.igc + (j+gd.jgc)*jj + k*kk;
                const int ijkb = j + (k-kstart)*kkb;
                tmp[ijkb] = data[ijk];
            }

#ifdef DISABLE_2D_MPIIO
//
// NOTE: on GPFS (IBM Spectrum Scale) file systems, saving small bits of data from
// each MPI task into a single binary using MPI-IO is **extremely** slow. With `DISABLE_2D_MPIIO`,
// the 2D slices are first gathered on MPI rank 0, and then written without MPI-IO.
//
        // MPI task which gathers/writes the slice:
        const int mpi_rank_recv = md.mpicoordx * md.npy;

        // Create send/receive MPI types.
        MPI_Datatype send_type;
        MPI_Type_vector(kmax, gd.jmax, gd.jmax, mpi_fp_type<TF>(), &send_type);
        MPI_Type_commit(&send_type);

        MPI_Datatype recv_type;
        int totyzsize_recv [2] = {kmax, gd.jtot};
        int subyzsize_recv [2] = {kmax, gd.jmax};
        int subyzstart_recv[2] = {0, md.mpicoordy*gd.jmax};
        MPI_Type_create_subarray(2, totyzsize_recv, subyzsize_recv, subyzstart_recv, MPI_ORDER_C, mpi_fp_type<TF>(), &recv_type);
        MPI_Type_commit(&recv_type);

        MPI_Datatype recv_type_r;
        MPI_Type_create_resized(recv_type, 0, sizeof(TF), &recv_type_r);
        MPI_Type_commit(&recv_type_r);

        // Create size/offset arrays for MPI_Gatherv().
        std::vector<int> counts(md.npy);
        std::fill(counts.begin(), counts.end(), 1);

        std::vector<int> offset(md.npy);
        for (int i=0; i<md.npy; ++i)
            offset[i] = i*gd.jmax;

        // Gather the data!
        std::vector<TF> recv;
        if (md.mpicoordy == mpi_rank_recv)
            recv.resize(gd.jtot*gd.ktot);

        MPI_Gatherv(tmp, 1, send_type, recv.data(), counts.data(), offset.data(), recv_type_r, mpi_rank_recv, md.commy);

        // Only MPI rank 0 writes the data.
        if (md.mpicoordy == mpi_rank_recv)
        {
            FILE *pFile;
            pFile = fopen(filename, "wbx");
            if (pFile == NULL)
                return 1;

            fwrite(recv.data(), sizeof(TF), gd.jtot*kmax, pFile);
            fclose(pFile);
        }
#else
//
// Use MPI-IO to write the data from different MPI tasks into a single file
//
        // Create MPI datatype for YZ-slice
        MPI_Datatype subyzslice;
        int totyzsize [2] = {kmax, gd.jtot};
        int subyzsize [2] = {kmax, gd.jmax};
        int subyzstart[2] = {0, md.mpicoordy*gd.jmax};
        MPI_Type_create_subarray(2, totyzsize, subyzsize, subyzstart, MPI_ORDER_C, mpi_fp_type<TF>(), &subyzslice);
        MPI_Type_commit(&subyzslice);

        MPI_File fh;
        if (MPI_File_open(md.commy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
            ++nerror;

        // select noncontiguous part of 3d array to store the selected data
        MPI_Offset fileoff = 0; // the offset within the file (header size)
        char name[] = "native";

        if (!nerror)
            if (MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subyzslice, name, MPI_INFO_NULL))
                ++nerror;

        // only write at the procs that contain the slice
        if (!nerror)
            if (MPI_File_write_all(fh, tmp, count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE))
                ++nerror;

        if (!nerror)
            MPI_File_sync(fh);

        if (!nerror)
            if (MPI_File_close(&fh))
                ++nerror;

        MPI_Type_free(&subyzslice);
#endif
    }

    // Gather errors from other processes
    master.sum(&nerror,1);
    MPI_Barrier(md.commxy);

    return nerror;
}

template<typename TF>
int Field3d_io<TF>::save_xy_slice(
        TF* const restrict data, TF* const restrict tmp,
        const char* filename, const int kslice)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    // Extract the data from the 3d field without the ghost cells
    const int jj  = gd.icells;
    const int kk  = gd.icells*gd.jcells;
    const int jjb = gd.imax;

    const int count = gd.imax*gd.jmax;

    for (int j=0; j<gd.jmax; j++)
        #pragma ivdep
        for (int i=0; i<gd.imax; i++)
        {
            // Take the modulus of jslice and jmax to have the right offset within proc
            const int ijk  = i+gd.igc + (j+gd.jgc)*jj + kslice*kk;
            const int ijkb = i + j*jjb;
            tmp[ijkb] = data[ijk];
        }

#ifdef DISABLE_2D_MPIIO
//
// NOTE: on GPFS (IBM Spectrum Scale) file systems, saving small bits of data from
// each MPI task into a single binary using MPI-IO is **extremely** slow. With `DISABLE_2D_MPIIO`,
// the 2D slices are first gathered on MPI rank 0, and then written without MPI-IO.
//
    // Create send/receive MPI types.
    MPI_Datatype send_type;
    MPI_Type_vector(gd.jmax, gd.imax, gd.imax, mpi_fp_type<TF>(), &send_type);
    MPI_Type_commit(&send_type);

    MPI_Datatype recv_type;
    int totxysize_recv [2] = {gd.jtot, gd.itot};
    int subxysize_recv [2] = {gd.jmax, gd.imax};
    int subxystart_recv[2] = {md.mpicoordy*gd.jmax, md.mpicoordx*gd.imax};
    MPI_Type_create_subarray(2, totxysize_recv, subxysize_recv, subxystart_recv, MPI_ORDER_C, mpi_fp_type<TF>(), &recv_type);
    MPI_Type_commit(&recv_type);

    MPI_Datatype recv_type_r;
    MPI_Type_create_resized(recv_type, 0, sizeof(TF), &recv_type_r);
    MPI_Type_commit(&recv_type_r);

    // Create size/offset arrays for MPI_Gatherv().
    std::vector<int> counts(md.nprocs);
    std::fill(counts.begin(), counts.end(), 1);

    std::vector<int> offset(md.nprocs);
    for (int i=0; i<md.npx; ++i)
        for (int j=0; j<md.npy; ++j)
        {
            const int ii = i+j*md.npx;
            offset[ii] = i*gd.imax + j*gd.jmax*gd.itot;
        }

    // TMP: benchmarking
    auto begin = std::chrono::high_resolution_clock::now();

    // Gather the data!
    std::vector<TF> recv = std::vector<TF>(gd.itot*gd.jtot);
    MPI_Gatherv(tmp, 1, send_type, recv.data(), counts.data(), offset.data(), recv_type_r, 0, md.commxy);

    // Only MPI rank 0 writes the data.
    if (md.mpiid == 0)
    {
        FILE *pFile;
        pFile = fopen(filename, "wbx");
        if (pFile == NULL)
            return 1;

        fwrite(recv.data(), sizeof(TF), gd.itot*gd.jtot, pFile);
        fclose(pFile);

        auto end = std::chrono::high_resolution_clock::now();
        auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
        const int bytes = gd.itot*gd.jtot*sizeof(TF);
        const double Mbytes_s = double(bytes)/double(diff) * 1e9/1024/1024;
        std::cout << "Write XY without MPI-IO = " << Mbytes_s << " MB/s" << std::endl;
    }

#else
//
// Use MPI-IO to write the data from different MPI tasks into a single file
//
    // Define MPI datatype for XY-slice
    MPI_Datatype subxyslice;
    int totxysize [2] = {gd.jtot, gd.itot};
    int subxysize [2] = {gd.jmax, gd.imax};
    int subxystart[2] = {md.mpicoordy*gd.jmax, md.mpicoordx*gd.imax};
    MPI_Type_create_subarray(2, totxysize, subxysize, subxystart, MPI_ORDER_C, mpi_fp_type<TF>(), &subxyslice);
    MPI_Type_commit(&subxyslice);

    // TMP: benchmarking
    auto begin = std::chrono::high_resolution_clock::now();

    MPI_File fh;
    if (MPI_File_open(md.commxy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
        return 1;

    // Select noncontiguous part of 3d array to store the selected data
    MPI_Offset fileoff = 0; // the offset within the file (header size)
    char name[] = "native";

    if (MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subxyslice, name, MPI_INFO_NULL))
        return 1;

    // Only write at the procs that contain the slice
    if (MPI_File_write_all(fh, tmp, count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE))
        return 1;

    MPI_File_sync(fh);

    if (MPI_File_close(&fh))
        return 1;

    MPI_Type_free(&subxyslice);

    MPI_Barrier(md.commxy);

    if (md.mpiid == 0)
    {
        auto end = std::chrono::high_resolution_clock::now();
        auto diff = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count();
        const int bytes = gd.itot*gd.jtot*sizeof(TF);
        const double Mbytes_s = double(bytes)/double(diff) * 1e9/1024/1024;
        std::cout << "Write XY with MPI-IO = " << Mbytes_s << " MB/s" << std::endl;
    }

#endif

    return 0;
}

template<typename TF>
int Field3d_io<TF>::load_xy_slice(
        TF* const restrict data, TF* const restrict tmp,
        const char* filename, int kslice)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    // Extract the data from the 3d field without the ghost cells
    const int jj = gd.icells;
    const int kk = gd.icells*gd.jcells;
    const int jjb = gd.imax;

    // Subtract the ghost cells in case of a pure 2d plane that does not have ghost cells.
    if (kslice == -1)
        kslice = -gd.kgc;

    int count = gd.imax*gd.jmax;

    // Create MPI datatype for XY-slice read
    MPI_Datatype subxyslice;
    int totxysize [2] = {gd.jtot, gd.itot};
    int subxysize [2] = {gd.jmax, gd.imax};
    int subxystart[2] = {md.mpicoordy*gd.jmax, md.mpicoordx*gd.imax};
    MPI_Type_create_subarray(2, totxysize, subxysize, subxystart, MPI_ORDER_C, mpi_fp_type<TF>(), &subxyslice);
    MPI_Type_commit(&subxyslice);

    MPI_File fh;
    if (MPI_File_open(md.commxy, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
        return 1;

    // Select noncontiguous part of 3d array to store the selected data
    MPI_Offset fileoff = 0; // the offset within the file (header size)
    char name[] = "native";

    if (MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subxyslice, name, MPI_INFO_NULL))
        return 1;

    // Only write at the procs that contain the slice
    if (MPI_File_read_all(fh, tmp, count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE))
        return 1;

    if (MPI_File_close(&fh))
        return 1;

    MPI_Type_free(&subxyslice);

    MPI_Barrier(md.commxy);

    for (int j=0; j<gd.jmax; j++)
        #pragma ivdep
        for (int i=0; i<gd.imax; i++)
        {
            // take the modulus of jslice and jmax to have the right offset within proc
            const int ijk  = i+gd.igc + (j+gd.jgc)*jj + (kslice+gd.kgc)*kk;
            const int ijkb = i + j*jjb;
            data[ijk] = tmp[ijkb];
        }

    return 0;
}

#else

template<typename TF>
int Field3d_io<TF>::save_field3d(
        TF* const restrict data,
        TF* const restrict tmp1, TF* const restrict tmp2,
        const char* filename, const TF offset,
        const int kstart, const int kend)
{
    auto& gd = grid.get_grid_data();

    FILE *pFile;
    pFile = fopen(filename, "wbx");

    if (pFile == NULL)
        return 1;

    const int jj = gd.icells;
    const int kk = gd.icells*gd.jcells;

    // First, add the offset to the data
    for (int k=kstart; k<kend; ++k)
        for (int j=gd.jstart; j<gd.jend; ++j)
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                tmp1[ijk] = data[ijk] + offset;
            }

    // Second, save the data to disk
    for (int k=kstart; k<kend; ++k)
        for (int j=gd.jstart; j<gd.jend; ++j)
        {
            const int ijk = gd.istart + j*jj + k*kk;
            if( fwrite(&tmp1[ijk], sizeof(TF), gd.imax, pFile) != (unsigned)gd.imax)
                return 1;
        }

    fclose(pFile);

    return 0;
}

template<typename TF>
int Field3d_io<TF>::load_field3d(
        TF* const restrict data,
        TF* const restrict tmp1, TF* const restrict tmp2,
        const char* filename, const TF offset,
        const int kstart, const int kend)
{
    auto& gd = grid.get_grid_data();

    FILE *pFile;
    pFile = fopen(filename, "rb");

    if (pFile == NULL)
        return 1;

    const int jj = gd.icells;
    const int kk = gd.icells*gd.jcells;

    // First, load the data from disk
    for (int k=kstart; k<kend; k++)
        for (int j=gd.jstart; j<gd.jend; j++)
        {
            const int ijk = gd.istart + j*jj + k*kk;
            if( fread(&tmp1[ijk], sizeof(TF), gd.imax, pFile) != (unsigned)gd.imax )
                return 1;
        }

    fclose(pFile);

    // Second, remove the offset
    for (int k=kstart; k<kend; k++)
        for (int j=gd.jstart; j<gd.jend; j++)
            for (int i=gd.istart; i<gd.iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                data[ijk] = tmp1[ijk] - offset;
            }

    return 0;
}

template<typename TF>
int Field3d_io<TF>::save_xz_slice(
        TF* const restrict data, TF* const restrict tmp,
        const char* filename, const int jslice,
        const int kstart, const int kend)
{
    auto& gd = grid.get_grid_data();

    // extract the data from the 3d field without the ghost cells
    const int jj   = gd.icells;
    const int kk   = gd.icells*gd.jcells;
    const int kkb  = gd.imax;
    const int kmax = kend-kstart;

    const int count = gd.imax*kmax;

    for (int k=kstart; k<kend; k++)
        #pragma ivdep
        for (int i=0; i<gd.imax; i++)
        {
            // take the modulus of jslice and gd.jmax to have the right offset within proc
            const int ijk  = i+gd.igc + (jslice+gd.jgc)*jj + k*kk;
            const int ijkb = i + (k-kstart)*kkb;
            tmp[ijkb] = data[ijk];
        }

    FILE *pFile;
    pFile = fopen(filename, "wbx");
    if (pFile == NULL)
        return 1;

    fwrite(tmp, sizeof(TF), count, pFile);
    fclose(pFile);

    return 0;
}

template<typename TF>
int Field3d_io<TF>::save_yz_slice(
        TF* const restrict data, TF* const restrict tmp,
        const char* filename, const int islice,
        const int kstart, const int kend)
{
    auto& gd = grid.get_grid_data();

    // Extract the data from the 3d field without the ghost cells
    const int jj   = gd.icells;
    const int kk   = gd.ijcells;
    const int kkb  = gd.jmax;
    const int kmax = kend-kstart;

    int count = gd.jmax*kmax;

    // Strip off the ghost cells
    for (int k=kstart; k<kend; k++)
        #pragma ivdep
        for (int j=0; j<gd.jmax; j++)
        {
            // take the modulus of jslice and jmax to have the right offset within proc
            const int ijk  = (islice%gd.imax)+gd.igc + (j+gd.jgc)*jj + k*kk;
            const int ijkb = j + (k-kstart)*kkb;
            tmp[ijkb] = data[ijk];
        }

    FILE *pFile;
    pFile = fopen(filename, "wbx");
    if (pFile == NULL)
        return 1;

    fwrite(tmp, sizeof(TF), count, pFile);
    fclose(pFile);

    return 0;
}

template<typename TF>
int Field3d_io<TF>::save_xy_slice(
        TF* const restrict data, TF* const restrict tmp,
        const char* filename, const int kslice)
{
    auto& gd = grid.get_grid_data();

    // extract the data from the 3d field without the ghost cells
    const int jj  = gd.icells;
    const int kk  = gd.icells*gd.jcells;
    const int jjb = gd.imax;

    const int count = gd.imax*gd.jmax;

    for (int j=0; j<gd.jmax; j++)
        #pragma ivdep
        for (int i=0; i<gd.imax; i++)
        {
            // Take the modulus of jslice and jmax to have the right offset within proc
            const int ijk  = i+gd.igc + (j+gd.jgc)*jj + kslice*kk;
            const int ijkb = i + j*jjb;
            tmp[ijkb] = data[ijk];
        }

    FILE *pFile;
    pFile = fopen(filename, "wbx");
    if (pFile == NULL)
        return 1;

    fwrite(tmp, sizeof(TF), count, pFile);
    fclose(pFile);

    return 0;
}

template<typename TF>
int Field3d_io<TF>::load_xy_slice(
        TF* const restrict data, TF* const restrict tmp,
        const char* filename, int kslice)
{
    auto& gd = grid.get_grid_data();

    const int count = gd.imax*gd.jmax;

    FILE *pFile;
    pFile = fopen(filename, "rb");
    if (pFile == NULL)
        return 1;

    fread(tmp, sizeof(TF), count, pFile);
    fclose(pFile);

    // Subtract the ghost cells in case of a pure 2d plane that does not have ghost cells.
    if (kslice == -1)
        kslice = -gd.kgc;

    // put the data back into a field with ghost cells
    const int jj  = gd.icells;
    const int kk  = gd.icells*gd.jcells;
    const int jjb = gd.imax;

    for (int j=0; j<gd.jmax; j++)
        #pragma ivdep
        for (int i=0; i<gd.imax; i++)
        {
            const int ijk  = i+gd.igc + (j+gd.jgc)*jj + (kslice+gd.kgc)*kk;
            const int ijkb = i + j*jjb;

            data[ijk] = tmp[ijkb];
        }

    return 0;
}
#endif

template class Field3d_io<double>;
template class Field3d_io<float>;
