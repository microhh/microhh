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

#ifdef USEMPI
#include <cstdio>
#include "master.h"
#include "grid.h"
#include "defines.h"

namespace
{
    template<typename TF> MPI_Datatype mpi_fp_type();
    template<> MPI_Datatype mpi_fp_type<double>() { return MPI_DOUBLE; }
    template<> MPI_Datatype mpi_fp_type<float>() { return MPI_FLOAT; }
}

// MPI functions
template<typename TF>
void Grid<TF>::init_mpi()
{
    // file saving and loading, take C-ordering into account
    int totsizei  = gd.itot;
    int subsizei  = gd.imax;
    int substarti = master.mpicoordx*gd.imax;
    MPI_Type_create_subarray(1, &totsizei, &subsizei, &substarti, MPI_ORDER_C, mpi_fp_type<TF>(), &subi);
    MPI_Type_commit(&subi);

    int totsizej  = gd.jtot;
    int subsizej  = gd.jmax;
    int substartj = master.mpicoordy*gd.jmax;
    MPI_Type_create_subarray(1, &totsizej, &subsizej, &substartj, MPI_ORDER_C, mpi_fp_type<TF>(), &subj);
    MPI_Type_commit(&subj);

    // the lines below describe the array in case transposes are not used before saving
    // int totsize [3] = {kmax, jtot, itot};
    // int subsize [3] = {kmax, jmax, imax};
    // int substart[3] = {0, master.mpicoordy*jmax, master.mpicoordx*imax};
    int totsize [3] = {gd.kmax  , gd.jtot, gd.itot};
    int subsize [3] = {gd.kblock, gd.jmax, gd.itot};
    int substart[3] = {master.mpicoordx*gd.kblock, master.mpicoordy*gd.jmax, 0};
    MPI_Type_create_subarray(3, totsize, subsize, substart, MPI_ORDER_C, mpi_fp_type<TF>(), &subarray);
    MPI_Type_commit(&subarray);

    // save mpitype for a xz-slice for cross section processing
    int totxzsize [2] = {gd.kmax, gd.itot};
    int subxzsize [2] = {gd.kmax, gd.imax};
    int subxzstart[2] = {0, master.mpicoordx*gd.imax};
    MPI_Type_create_subarray(2, totxzsize, subxzsize, subxzstart, MPI_ORDER_C, mpi_fp_type<TF>(), &subxzslice);
    MPI_Type_commit(&subxzslice);

    // save mpitype for a yz-slice for cross section processing
    int totyzsize [2] = {gd.kmax, gd.jtot};
    int subyzsize [2] = {gd.kmax, gd.jmax};
    int subyzstart[2] = {0, master.mpicoordy*gd.jmax};
    MPI_Type_create_subarray(2, totyzsize, subyzsize, subyzstart, MPI_ORDER_C, mpi_fp_type<TF>(), &subyzslice);
    MPI_Type_commit(&subyzslice);

    // save mpitype for a xy-slice for cross section processing
    int totxysize [2] = {gd.jtot, gd.itot};
    int subxysize [2] = {gd.jmax, gd.imax};
    int subxystart[2] = {master.mpicoordy*gd.jmax, master.mpicoordx*gd.imax};
    MPI_Type_create_subarray(2, totxysize, subxysize, subxystart, MPI_ORDER_C, mpi_fp_type<TF>(), &subxyslice);
    MPI_Type_commit(&subxyslice);

    mpitypes = true;
}

template<typename TF>
void Grid<TF>::exit_mpi()
{
    if (mpitypes)
    {
        MPI_Type_free(&subi);
        MPI_Type_free(&subj);
        MPI_Type_free(&subarray);
        MPI_Type_free(&subxzslice);
        MPI_Type_free(&subyzslice);
        MPI_Type_free(&subxyslice);
    }
}

template<typename TF>
void Grid<TF>::save_grid()
{
    char filename[256];
    std::sprintf(filename, "%s.%07d", "grid", 0);
    master.print_message("Saving \"%s\" ... ", filename);

    MPI_File fh;
    if (MPI_File_open(master.commxy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
    {
        master.print_message("FAILED\n");
        throw 1;
    }

    // select noncontiguous part of 3d array to store the selected data
    MPI_Offset fileoff = 0; // the offset within the file (header size)
    char name[] = "native";

    MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subi, name, MPI_INFO_NULL);
    if (master.mpicoordy == 0)
        MPI_File_write(fh, &gd.x[gd.istart], gd.imax, mpi_fp_type<TF>(), MPI_STATUS_IGNORE);
    MPI_Barrier(master.commxy);
    fileoff += gd.itot*sizeof(TF);

    MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subi, name, MPI_INFO_NULL);
    if (master.mpicoordy == 0)
        MPI_File_write(fh, &gd.xh[gd.istart], gd.imax, mpi_fp_type<TF>(), MPI_STATUS_IGNORE);
    MPI_Barrier(master.commxy);
    fileoff += gd.itot*sizeof(TF);

    MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subj, name, MPI_INFO_NULL);
    if (master.mpicoordx == 0)
        MPI_File_write(fh, &gd.y[gd.jstart], gd.jmax, mpi_fp_type<TF>(), MPI_STATUS_IGNORE);
    MPI_Barrier(master.commxy);
    fileoff += gd.jtot*sizeof(TF);

    MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subj, name, MPI_INFO_NULL);
    if (master.mpicoordx == 0)
        MPI_File_write(fh, &gd.yh[gd.jstart], gd.jmax, mpi_fp_type<TF>(), MPI_STATUS_IGNORE);
    MPI_Barrier(master.commxy);

    MPI_File_sync(fh);
    if (MPI_File_close(&fh))
        throw 1;

    if (master.mpiid == 0)
    {
        FILE *pFile;
        pFile = fopen(filename, "ab");
        fwrite(&gd.z [gd.kstart], sizeof(TF), gd.kmax, pFile);
        fwrite(&gd.zh[gd.kstart], sizeof(TF), gd.kmax, pFile);
        fclose(pFile);
    }

    // the saving procedure is a success
    master.print_message("OK\n");
}

template<typename TF>
void Grid<TF>::load_grid()
{
    int nerror = 0;

    // LOAD THE GRID
    char filename[256];
    std::sprintf(filename, "%s.%07d", "grid", 0);
    if (master.mpiid == 0) std::printf("Loading \"%s\" ... ", filename);

    FILE *pFile;
    if (master.mpiid == 0)
    {
        pFile = fopen(filename, "rb");
        if (pFile == NULL)
        {
            ++nerror;
        }
        else
        {
            int n = (2*gd.itot+2*gd.jtot)*sizeof(TF);
            fseek(pFile, n, SEEK_SET);
            fread(&gd.z [gd.kstart], sizeof(TF), gd.kmax, pFile);
            fread(&gd.zh[gd.kstart], sizeof(TF), gd.kmax, pFile);
            fclose(pFile);
        }
    }

    // Communicate the file read error over all procs.
    master.broadcast(&nerror, 1);
    if (nerror)
    {
        master.print_message("FAILED\n");
        throw 1;
    }
    else
        master.print_message("OK\n");

    master.broadcast(&gd.z [gd.kstart], gd.kmax);
    master.broadcast(&gd.zh[gd.kstart], gd.kmax);

    // Calculate the missing coordinates.
    calculate();
}

template<typename TF>
int Grid<TF>::save_field3d(TF* restrict data, TF* restrict tmp1, TF* restrict tmp2, char* filename, TF offset)
{
    // save the data in transposed order to have large chunks of contiguous disk space
    // MPI-IO is not stable on Juqueen and supermuc otherwise

    // extract the data from the 3d field without the ghost cells
    const int jj  = gd.icells;
    const int kk  = gd.icells*gd.jcells;
    const int jjb = gd.imax;
    const int kkb = gd.imax*gd.jmax;

    int count = gd.imax*gd.jmax*gd.kmax;

    for (int k=0; k<gd.kmax; ++k)
        for (int j=0; j<gd.jmax; ++j)
            #pragma ivdep
            for (int i=0; i<gd.imax; ++i)
            {
                const int ijk  = i+gd.igc + (j+gd.jgc)*jj + (k+gd.kgc)*kk;
                const int ijkb = i + j*jjb + k*kkb;
                tmp1[ijkb] = data[ijk] + offset;
            }

    transpose.exec_zx(tmp2, tmp1);

    MPI_File fh;
    if (MPI_File_open(master.commxy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
        return 1;

    // select noncontiguous part of 3d array to store the selected data
    MPI_Offset fileoff = 0; // the offset within the file (header size)
    char name[] = "native";

    if (MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subarray, name, MPI_INFO_NULL))
        return 1;

    if (MPI_File_write_all(fh, tmp2, count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE))
        return 1;

    if (MPI_File_close(&fh))
        return 1;

    return 0;
}

template<typename TF>
int Grid<TF>::load_field3d(TF* const restrict data, TF* const restrict tmp1, TF* const restrict tmp2, char* filename, TF offset)
{
    // Save the data in transposed order to have large chunks of contiguous disk space.
    // MPI-IO is not stable on Juqueen and supermuc otherwise.

    // read the file
    MPI_File fh;
    if (MPI_File_open(master.commxy, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
        return 1;

    // select noncontiguous part of 3d array to store the selected data
    MPI_Offset fileoff = 0; // the offset within the file (header size)
    char name[] = "native";
    MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subarray, name, MPI_INFO_NULL);

    // extract the data from the 3d field without the ghost cells
    int count = gd.imax*gd.jmax*gd.kmax;

    if (MPI_File_read_all(fh, tmp1, count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE))
        return 1;

    if (MPI_File_close(&fh))
        return 1;

    // transpose the data back
    transpose.exec_xz(tmp2, tmp1);

    const int jj  = gd.icells;
    const int kk  = gd.icells*gd.jcells;
    const int jjb = gd.imax;
    const int kkb = gd.imax*gd.jmax;

    for (int k=0; k<gd.kmax; ++k)
        for (int j=0; j<gd.jmax; ++j)
            #pragma ivdep
            for (int i=0; i<gd.imax; ++i)
            {
                const int ijk  = i+gd.igc + (j+gd.jgc)*jj + (k+gd.kgc)*kk;
                const int ijkb = i + j*jjb + k*kkb;
                data[ijk] = tmp2[ijkb] - offset;
            }

    return 0;
}

template<typename TF>
int Grid<TF>::save_xz_slice(TF* restrict data, TF* restrict tmp, char* filename, int jslice)
{
    // extract the data from the 3d field without the ghost cells
    int nerror=0;

    const int jj  = gd.icells;
    const int kk  = gd.icells*gd.jcells;
    const int kkb = gd.imax;

    int count = gd.imax*gd.kmax;

    for (int k=0; k<gd.kmax; k++)
#pragma ivdep
        for (int i=0; i<gd.imax; i++)
        {
            // take the modulus of jslice and gd.jmax to have the right offset within proc
            const int ijk  = i+gd.igc + ((jslice%gd.jmax)+gd.jgc)*jj + (k+gd.kgc)*kk;
            const int ijkb = i + k*kkb;
            tmp[ijkb] = data[ijk];
        }

    if (master.mpicoordy == jslice/gd.jmax)
    {
        MPI_File fh;
        if (MPI_File_open(master.commx, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
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
    }

    // Gather errors from other processes
    master.sum(&nerror,1);

    MPI_Barrier(master.commxy);

    return nerror;
}

template<typename TF>
int Grid<TF>::save_yz_slice(TF* restrict data, TF* restrict tmp, char* filename, int islice)
{
    // extract the data from the 3d field without the ghost cells
    int nerror=0;

    const int jj = gd.icells;
    const int kk = gd.ijcells;

    const int kkb = gd.jmax;

    int count = gd.jmax*gd.kmax;

    // Strip off the ghost cells
    for (int k=0; k<gd.kmax; k++)
        #pragma ivdep
        for (int j=0; j<gd.jmax; j++)
        {
            // take the modulus of jslice and jmax to have the right offset within proc
            const int ijk  = (islice%gd.imax)+gd.igc + (j+gd.jgc)*jj + (k+gd.kgc)*kk;
            const int ijkb = j + k*kkb;
            tmp[ijkb] = data[ijk];
        }

    if (master.mpicoordx == islice/gd.imax)
    {
        MPI_File fh;
        if (MPI_File_open(master.commy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
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
    }

    // Gather errors from other processes
    master.sum(&nerror,1);

    MPI_Barrier(master.commxy);

    return nerror;
}

template<typename TF>
int Grid<TF>::save_xy_slice(TF* restrict data, TF* restrict tmp, char* filename, int kslice)
{
    // extract the data from the 3d field without the ghost cells
    const int jj  = gd.icells;
    const int kk  = gd.icells*gd.jcells;
    const int jjb = gd.imax;

    // Subtract the ghost cells in case of a pure 2d plane that does not have ghost cells.
    if (kslice == -1)
        kslice = -gd.kgc;

    int count = gd.imax*gd.jmax;

    for (int j=0; j<gd.jmax; j++)
#pragma ivdep
        for (int i=0; i<gd.imax; i++)
        {
            // take the modulus of jslice and jmax to have the right offset within proc
            const int ijk  = i+gd.igc + (j+gd.jgc)*jj + (kslice+gd.kgc)*kk;
            const int ijkb = i + j*jjb;
            tmp[ijkb] = data[ijk];
        }

    MPI_File fh;
    if (MPI_File_open(master.commxy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
        return 1;

    // select noncontiguous part of 3d array to store the selected data
    MPI_Offset fileoff = 0; // the offset within the file (header size)
    char name[] = "native";

    if (MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subxyslice, name, MPI_INFO_NULL))
        return 1;

    // only write at the procs that contain the slice
    if (MPI_File_write_all(fh, tmp, count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE))
        return 1;

    MPI_File_sync(fh);

    if (MPI_File_close(&fh))
        return 1;

    MPI_Barrier(master.commxy);

    return 0;
}

template<typename TF>
int Grid<TF>::load_xy_slice(TF* restrict data, TF* restrict tmp, char* filename, int kslice)
{
    // extract the data from the 3d field without the ghost cells
    const int jj  = gd.icells;
    const int kk  = gd.icells*gd.jcells;
    const int jjb = gd.imax;

    // Subtract the ghost cells in case of a pure 2d plane that does not have ghost cells.
    if (kslice == -1)
        kslice = -gd.kgc;

    int count = gd.imax*gd.jmax;

    MPI_File fh;
    if (MPI_File_open(master.commxy, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
        return 1;

    // select noncontiguous part of 3d array to store the selected data
    MPI_Offset fileoff = 0; // the offset within the file (header size)
    char name[] = "native";

    if (MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subxyslice, name, MPI_INFO_NULL))
        return 1;

    // only write at the procs that contain the slice
    if (MPI_File_read_all(fh, tmp, count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE))
        return 1;

    if (MPI_File_close(&fh))
        return 1;

    MPI_Barrier(master.commxy);

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

template class Grid<double>;
template class Grid<float>;
#endif
