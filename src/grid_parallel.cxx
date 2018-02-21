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
#include <fftw3.h>
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

    // east west 2d
    datacount  = gd.jcells;
    datablock  = gd.igc;
    datastride = gd.icells;
    MPI_Type_vector(datacount, datablock, datastride, mpi_fp_type<TF>(), &eastwestedge2d);
    MPI_Type_commit(&eastwestedge2d);

    // north south 2d
    datacount  = 1;
    datablock  = gd.icells*gd.jgc;
    datastride = gd.icells*gd.jcells;
    MPI_Type_vector(datacount, datablock, datastride, mpi_fp_type<TF>(), &northsouthedge2d);
    MPI_Type_commit(&northsouthedge2d);

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

    // allocate the array for the profiles
    // profl = new double[kcells];

    mpitypes = true;
}

template<typename TF>
void Grid<TF>::exit_mpi()
{
    if (mpitypes)
    {
        MPI_Type_free(&eastwestedge);
        MPI_Type_free(&northsouthedge);
        MPI_Type_free(&eastwestedge2d);
        MPI_Type_free(&northsouthedge2d);
        MPI_Type_free(&transposez);
        MPI_Type_free(&transposez2);
        MPI_Type_free(&transposex);
        MPI_Type_free(&transposex2);
        MPI_Type_free(&transposey);
        MPI_Type_free(&transposey2);
        MPI_Type_free(&subi);
        MPI_Type_free(&subj);
        MPI_Type_free(&subarray);
        MPI_Type_free(&subxzslice);
        MPI_Type_free(&subyzslice);
        MPI_Type_free(&subxyslice);

        // delete[] profl;
    }
}

template<typename TF>
void Grid<TF>::boundary_cyclic(TF* const restrict data, Edge edge)
{
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

template<typename TF>
void Grid<TF>::boundary_cyclic_2d(TF* const restrict data)
{
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
    MPI_Isend(&data[eastout], ncount, eastwestedge2d, master.neast, 1, master.commxy, &master.reqs[master.reqsn]);
    master.reqsn++;
    MPI_Irecv(&data[westin], ncount, eastwestedge2d, master.nwest, 1, master.commxy, &master.reqs[master.reqsn]);
    master.reqsn++;
    MPI_Isend(&data[westout], ncount, eastwestedge2d, master.nwest, 2, master.commxy, &master.reqs[master.reqsn]);
    master.reqsn++;
    MPI_Irecv(&data[eastin], ncount, eastwestedge2d, master.neast, 2, master.commxy, &master.reqs[master.reqsn]);
    master.reqsn++;
    // Wait here for the mpi to have correct values in the corners of the cells.
    master.wait_all();

    // If the run is 3D, apply the BCs.
    if (gd.jtot > 1)
    {
        // Second, send and receive the ghost cells in the north-south direction.
        MPI_Isend(&data[northout], ncount, northsouthedge2d, master.nnorth, 1, master.commxy, &master.reqs[master.reqsn]);
        master.reqsn++;
        MPI_Irecv(&data[southin], ncount, northsouthedge2d, master.nsouth, 1, master.commxy, &master.reqs[master.reqsn]);
        master.reqsn++;
        MPI_Isend(&data[southout], ncount, northsouthedge2d, master.nsouth, 2, master.commxy, &master.reqs[master.reqsn]);
        master.reqsn++;
        MPI_Irecv(&data[northin], ncount, northsouthedge2d, master.nnorth, 2, master.commxy, &master.reqs[master.reqsn]);
        master.reqsn++;
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
void Grid<TF>::transpose_zx(TF* const restrict ar, TF* const restrict as)
{
    const int ncount = 1;
    const int tag = 1;

    const int jj = gd.imax;
    const int kk = gd.imax*gd.jmax;

    for (int n=0; n<master.npx; ++n)
    {
        // Determine where to fetch the data and where to store it.
        const int ijks = n*gd.kblock*kk;
        const int ijkr = n*jj;

        // Send and receive the data.
        MPI_Isend(&as[ijks], ncount, transposez, n, tag, master.commx, &master.reqs[master.reqsn]);
        master.reqsn++;
        MPI_Irecv(&ar[ijkr], ncount, transposex, n, tag, master.commx, &master.reqs[master.reqsn]);
        master.reqsn++;
    }

    master.wait_all();
}

template<typename TF>
void Grid<TF>::transpose_xz(TF* const restrict ar, TF* const restrict as)
{
    const int ncount = 1;
    const int tag = 1;

    const int jj = gd.imax;
    const int kk = gd.imax*gd.jmax;

    for (int n=0; n<master.npx; ++n)
    {
        // Determine where to fetch the data and where to store it.
        const int ijks = n*jj;
        const int ijkr = n*gd.kblock*kk;

        // Send and receive the data.
        MPI_Isend(&as[ijks], ncount, transposex, n, tag, master.commx, &master.reqs[master.reqsn]);
        master.reqsn++;
        MPI_Irecv(&ar[ijkr], ncount, transposez, n, tag, master.commx, &master.reqs[master.reqsn]);
        master.reqsn++;
    }

    master.wait_all();
}

template<typename TF>
void Grid<TF>::transpose_xy(TF* const restrict ar, TF* const restrict as)
{
    const int ncount = 1;
    const int tag = 1;

    const int jj = gd.iblock;
    const int kk = gd.iblock*gd.jmax;

    for (int n=0; n<master.npy; ++n)
    {
        // determine where to fetch the data and where to store it
        const int ijks = n*jj;
        const int ijkr = n*kk;

        // send and receive the data
        MPI_Isend(&as[ijks], ncount, transposex2, n, tag, master.commy, &master.reqs[master.reqsn]);
        master.reqsn++;
        MPI_Irecv(&ar[ijkr], ncount, transposey , n, tag, master.commy, &master.reqs[master.reqsn]);
        master.reqsn++;
    }

    master.wait_all();
}

template<typename TF>
void Grid<TF>::transpose_yx(TF* const restrict ar, TF* const restrict as)
{
    const int ncount = 1;
    const int tag = 1;

    const int jj = gd.iblock;
    const int kk = gd.iblock*gd.jmax;

    for (int n=0; n<master.npy; ++n)
    {
        // determine where to fetch the data and where to store it
        const int ijks = n*kk;
        const int ijkr = n*jj;

        // send and receive the data
        MPI_Isend(&as[ijks], ncount, transposey , n, tag, master.commy, &master.reqs[master.reqsn]);
        master.reqsn++;
        MPI_Irecv(&ar[ijkr], ncount, transposex2, n, tag, master.commy, &master.reqs[master.reqsn]);
        master.reqsn++;
    }

    master.wait_all();
}

template<typename TF>
void Grid<TF>::transpose_yz(TF* const restrict ar, TF* const restrict as)
{
    const int ncount = 1;
    const int tag = 1;

    const int jj = gd.iblock;
    const int kk = gd.iblock*gd.jblock;

    for (int n=0; n<master.npx; ++n)
    {
        // determine where to fetch the data and where to store it
        const int ijks = n*gd.jblock*jj;
        const int ijkr = n*gd.kblock*kk;

        // send and receive the data
        MPI_Isend(&as[ijks], ncount, transposey2, n, tag, master.commx, &master.reqs[master.reqsn]);
        master.reqsn++;
        MPI_Irecv(&ar[ijkr], ncount, transposez2, n, tag, master.commx, &master.reqs[master.reqsn]);
        master.reqsn++;
    }

    master.wait_all();
}

template<typename TF>
void Grid<TF>::transpose_zy(TF* const restrict ar, TF* const restrict as)
{
    const int ncount = 1;
    const int tag = 1;

    const int jj = gd.iblock;
    const int kk = gd.iblock*gd.jblock;

    for (int n=0; n<master.npx; ++n)
    {
        // determine where to fetch the data and where to store it
        const int ijks = n*gd.kblock*kk;
        const int ijkr = n*gd.jblock*jj;

        // send and receive the data
        MPI_Isend(&as[ijks], ncount, transposez2, n, tag, master.commx, &master.reqs[master.reqsn]);
        master.reqsn++;
        MPI_Irecv(&ar[ijkr], ncount, transposey2, n, tag, master.commx, &master.reqs[master.reqsn]);
        master.reqsn++;
    }

    master.wait_all();
}

// template<typename TF>
// void Grid<TF>::get_max(double *var)
// {
//     double varl = *var;
//     MPI_Allreduce(&varl, var, 1, mpi_fp_type<TF>(), MPI_MAX, master.commxy);
// }
//
// template<typename TF>
// void Grid<TF>::get_max(int *var)
// {
//     int varl = *var;
//     MPI_Allreduce(&varl, var, 1, MPI_INT, MPI_MAX, master.commxy);
// }
//
// template<typename TF>
// void Grid<TF>::get_sum(double *var)
// {
//     double varl = *var;
//     MPI_Allreduce(&varl, var, 1, mpi_fp_type<TF>(), MPI_SUM, master.commxy);
// }
//
// template<typename TF>
// void Grid<TF>::get_prof(double *prof, int kcellsin)
// {
//     for (int k=0; k<kcellsin; k++)
//         profl[k] = prof[k] / master.nprocs;
//
//     MPI_Allreduce(profl, prof, kcellsin, mpi_fp_type<TF>(), MPI_SUM, master.commxy);
// }

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

    transpose_zx(tmp2, tmp1);

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
    transpose_xz(tmp2, tmp1);

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

namespace
{
    template<typename> void fftw_execute_wrapper(const fftw_plan&, const fftwf_plan&);

    template<>
    void fftw_execute_wrapper<double>(const fftw_plan& p, const fftwf_plan& pf)
    {
        fftw_execute(p);
    }

    template<>
    void fftw_execute_wrapper<float>(const fftw_plan& p, const fftwf_plan& pf)
    {
        fftwf_execute(pf);
    }
}

template<typename TF>
void Grid<TF>::fft_forward(TF* const restrict data,   TF* const restrict tmp1,
                           TF* const restrict fftini, TF* const restrict fftouti,
                           TF* const restrict fftinj, TF* const restrict fftoutj)
{
    // Transpose the pressure field.
    transpose_zx(tmp1, data);

    int kk = gd.itot*gd.jmax;

    // Process the fourier transforms slice by slice.
    for (int k=0; k<gd.kblock; ++k)
    {
        #pragma ivdep
        for (int n=0; n<gd.itot*gd.jmax; ++n)
        {
            const int ij = n;
            const int ijk = n + k*kk;
            fftini[ij] = tmp1[ijk];
        }

        fftw_execute_wrapper<TF>(iplanf, iplanff);

        #pragma ivdep
        for (int n=0; n<gd.itot*gd.jmax; ++n)
        {
            const int ij = n;
            const int ijk = n + k*kk;
            tmp1[ijk] = fftouti[ij];
        }
    }

    // Transpose again.
    transpose_xy(data, tmp1);

    kk = gd.iblock*gd.jtot;

    // Do the second fourier transform.
    for (int k=0; k<gd.kblock; ++k)
    {
        #pragma ivdep
        for (int n=0; n<gd.iblock*gd.jtot; ++n)
        {
            const int ij = n;
            const int ijk = n + k*kk;
            fftinj[ij] = data[ijk];
        }

        fftw_execute_wrapper<TF>(jplanf, jplanff);

        #pragma ivdep
        for (int n=0; n<gd.iblock*gd.jtot; ++n)
        {
            const int ij = n;
            const int ijk = n + k*kk;
            // Shift to use p in pressure solver.
            tmp1[ijk] = fftoutj[ij];
        }
    }

    // Transpose back to original orientation.
    transpose_yz(data, tmp1);
}

template<typename TF>
void Grid<TF>::fft_backward(TF* const restrict data,   TF* const restrict tmp1,
                            TF* const restrict fftini, TF* const restrict fftouti,
                            TF* const restrict fftinj, TF* const restrict fftoutj)
{
    // Transpose back to y.
    transpose_zy(tmp1, data);

    int kk = gd.iblock*gd.jtot;

    // Transform the second transform back.
    for (int k=0; k<gd.kblock; ++k)
    {
        #pragma ivdep
        for (int n=0; n<gd.iblock*gd.jtot; ++n)
        {
            const int ij = n;
            const int ijk = n + k*kk;
            fftinj[ij] = tmp1[ijk];
        }

        fftw_execute_wrapper<TF>(jplanb, jplanbf);

        #pragma ivdep
        for (int n=0; n<gd.iblock*gd.jtot; ++n)
        {
            const int ij = n;
            const int ijk = n + k*kk;
            data[ijk] = fftoutj[ij] / gd.jtot;
        }
    }

    // Transpose back to x.
    transpose_yx(tmp1, data);

    kk = gd.itot*gd.jmax;

    // Transform the first transform back.
    for (int k=0; k<gd.kblock; ++k)
    {
        #pragma ivdep
        for (int n=0; n<gd.itot*gd.jmax; ++n)
        {
            const int ij = n;
            const int ijk = n + k*kk;
            fftini[ij] = tmp1[ijk];
        }

        fftw_execute_wrapper<TF>(iplanb, iplanbf);

        #pragma ivdep
        for (int n=0; n<gd.itot*gd.jmax; ++n)
        {
            const int ij = n;
            const int ijk = n + k*kk;
            // swap array here to avoid unnecessary 3d loop
            data[ijk] = fftouti[ij] / gd.itot;
        }
    }

    // And transpose back...
    transpose_xz(tmp1, data);
}

//
// template<typename TF>
// int Grid<TF>::save_xz_slice(double* restrict data, double* restrict tmp, char* filename, int jslice)
// {
//     // extract the data from the 3d field without the ghost cells
//     int nerror=0;
//
//     const int jj  = icells;
//     const int kk  = icells*jcells;
//     const int kkb = imax;
//
//     int count = imax*kmax;
//
//     for (int k=0; k<kmax; k++)
// #pragma ivdep
//         for (int i=0; i<imax; i++)
//         {
//             // take the modulus of jslice and jmax to have the right offset within proc
//             const int ijk  = i+igc + ((jslice%jmax)+jgc)*jj + (k+kgc)*kk;
//             const int ijkb = i + k*kkb;
//             tmp[ijkb] = data[ijk];
//         }
//
//     if (master.mpicoordy == jslice/jmax)
//     {
//         MPI_File fh;
//         if (MPI_File_open(master.commx, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
//             ++nerror;
//
//         // select noncontiguous part of 3d array to store the selected data
//         MPI_Offset fileoff = 0; // the offset within the file (header size)
//         char name[] = "native";
//
//         if (!nerror)
//             if (MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subxzslice, name, MPI_INFO_NULL))
//                 ++nerror;
//
//         // only write at the procs that contain the slice
//         if (!nerror)
//             if (MPI_File_write_all(fh, tmp, count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE))
//                 ++nerror;
//
//         if (!nerror)
//             MPI_File_sync(fh);
//
//         if (!nerror)
//             if (MPI_File_close(&fh))
//                 ++nerror;
//     }
//
//     // Gather errors from other processes
//     master.sum(&nerror,1);
//
//     MPI_Barrier(master.commxy);
//
//     return nerror;
// }
//
// template<typename TF>
// int Grid<TF>::save_yz_slice(double* restrict data, double* restrict tmp, char* filename, int islice)
// {
//     // extract the data from the 3d field without the ghost cells
//     int nerror=0;
//
//     const int jj = icells;
//     const int kk = ijcells;
//
//     const int kkb = jmax;
//
//     int count = jmax*kmax;
//
//     // Strip off the ghost cells
//     for (int k=0; k<kmax; k++)
//         #pragma ivdep
//         for (int j=0; j<jmax; j++)
//         {
//             // take the modulus of jslice and jmax to have the right offset within proc
//             const int ijk  = (islice%imax)+igc + (j+jgc)*jj + (k+kgc)*kk;
//             const int ijkb = j + k*kkb;
//             tmp[ijkb] = data[ijk];
//         }
//
//     if (master.mpicoordx == islice/imax)
//     {
//         MPI_File fh;
//         if (MPI_File_open(master.commy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
//             ++nerror;
//
//         // select noncontiguous part of 3d array to store the selected data
//         MPI_Offset fileoff = 0; // the offset within the file (header size)
//         char name[] = "native";
//
//         if (!nerror)
//             if (MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subyzslice, name, MPI_INFO_NULL))
//                 ++nerror;
//
//         // only write at the procs that contain the slice
//         if (!nerror)
//             if (MPI_File_write_all(fh, tmp, count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE))
//                 ++nerror;
//
//         if (!nerror)
//             MPI_File_sync(fh);
//
//         if (!nerror)
//             if (MPI_File_close(&fh))
//                 ++nerror;
//     }
//
//     // Gather errors from other processes
//     master.sum(&nerror,1);
//
//     MPI_Barrier(master.commxy);
//
//     return nerror;
// }
//
// template<typename TF>
// int Grid<TF>::save_xy_slice(double* restrict data, double* restrict tmp, char* filename, int kslice)
// {
//     // extract the data from the 3d field without the ghost cells
//     const int jj  = icells;
//     const int kk  = icells*jcells;
//     const int jjb = imax;
//
//     // Subtract the ghost cells in case of a pure 2d plane that does not have ghost cells.
//     if (kslice == -1)
//         kslice = -kgc;
//
//     int count = imax*jmax;
//
//     for (int j=0; j<jmax; j++)
// #pragma ivdep
//         for (int i=0; i<imax; i++)
//         {
//             // take the modulus of jslice and jmax to have the right offset within proc
//             const int ijk  = i+igc + (j+jgc)*jj + (kslice+kgc)*kk;
//             const int ijkb = i + j*jjb;
//             tmp[ijkb] = data[ijk];
//         }
//
//     MPI_File fh;
//     if (MPI_File_open(master.commxy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
//         return 1;
//
//     // select noncontiguous part of 3d array to store the selected data
//     MPI_Offset fileoff = 0; // the offset within the file (header size)
//     char name[] = "native";
//
//     if (MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subxyslice, name, MPI_INFO_NULL))
//         return 1;
//
//     // only write at the procs that contain the slice
//     if (MPI_File_write_all(fh, tmp, count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE))
//         return 1;
//
//     MPI_File_sync(fh);
//
//     if (MPI_File_close(&fh))
//         return 1;
//
//     MPI_Barrier(master.commxy);
//
//     return 0;
// }
//
// template<typename TF>
// int Grid<TF>::load_xy_slice(double* restrict data, double* restrict tmp, char* filename, int kslice)
// {
//     // extract the data from the 3d field without the ghost cells
//     const int jj  = icells;
//     const int kk  = icells*jcells;
//     const int jjb = imax;
//
//     // Subtract the ghost cells in case of a pure 2d plane that does not have ghost cells.
//     if (kslice == -1)
//         kslice = -kgc;
//
//     int count = imax*jmax;
//
//     MPI_File fh;
//     if (MPI_File_open(master.commxy, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
//         return 1;
//
//     // select noncontiguous part of 3d array to store the selected data
//     MPI_Offset fileoff = 0; // the offset within the file (header size)
//     char name[] = "native";
//
//     if (MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subxyslice, name, MPI_INFO_NULL))
//         return 1;
//
//     // only write at the procs that contain the slice
//     if (MPI_File_read_all(fh, tmp, count, mpi_fp_type<TF>(), MPI_STATUS_IGNORE))
//         return 1;
//
//     if (MPI_File_close(&fh))
//         return 1;
//
//     MPI_Barrier(master.commxy);
//
//     for (int j=0; j<jmax; j++)
// #pragma ivdep
//         for (int i=0; i<imax; i++)
//         {
//             // take the modulus of jslice and jmax to have the right offset within proc
//             const int ijk  = i+igc + (j+jgc)*jj + (kslice+kgc)*kk;
//             const int ijkb = i + j*jjb;
//             data[ijk] = tmp[ijkb];
//         }
//
//     return 0;
// }

template class Grid<double>;
template class Grid<float>;
#endif
