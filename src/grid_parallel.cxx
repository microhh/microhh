/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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

// MPI functions
void Grid::init_mpi()
{
    // create the MPI types for the cyclic boundary conditions
    int datacount, datablock, datastride;

    // east west
    datacount  = jcells*kcells;
    datablock  = igc;
    datastride = icells;
    MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &eastwestedge);
    MPI_Type_commit(&eastwestedge);

    // north south
    datacount  = kcells;
    datablock  = icells*jgc;
    datastride = icells*jcells;
    MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &northsouthedge);
    MPI_Type_commit(&northsouthedge);

    // east west 2d
    datacount  = jcells;
    datablock  = igc;
    datastride = icells;
    MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &eastwestedge2d);
    MPI_Type_commit(&eastwestedge2d);

    // north south 2d
    datacount  = 1;
    datablock  = icells*jgc;
    datastride = icells*jcells;
    MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &northsouthedge2d);
    MPI_Type_commit(&northsouthedge2d);

    // transposez
    datacount = imax*jmax*kblock;
    MPI_Type_contiguous(datacount, MPI_DOUBLE, &transposez);
    MPI_Type_commit(&transposez);

    // transposez iblock/jblock/kblock
    datacount = iblock*jblock*kblock;
    MPI_Type_contiguous(datacount, MPI_DOUBLE, &transposez2);
    MPI_Type_commit(&transposez2);

    // transposex imax
    datacount  = jmax*kblock;
    datablock  = imax;
    datastride = itot;
    MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &transposex);
    MPI_Type_commit(&transposex);

    // transposex iblock
    datacount  = jmax*kblock;
    datablock  = iblock;
    datastride = itot;
    MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &transposex2);
    MPI_Type_commit(&transposex2);

    // transposey
    datacount  = kblock;
    datablock  = iblock*jmax;
    datastride = iblock*jtot;
    MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &transposey);
    MPI_Type_commit(&transposey);

    // transposey2
    datacount  = kblock;
    datablock  = iblock*jblock;
    datastride = iblock*jtot;
    MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE, &transposey2);
    MPI_Type_commit(&transposey2);

    // file saving and loading, take C-ordering into account
    int totsizei  = itot;
    int subsizei  = imax;
    int substarti = master->mpicoordx*imax;
    MPI_Type_create_subarray(1, &totsizei, &subsizei, &substarti, MPI_ORDER_C, MPI_DOUBLE, &subi);
    MPI_Type_commit(&subi);

    int totsizej  = jtot;
    int subsizej  = jmax;
    int substartj = master->mpicoordy*jmax;
    MPI_Type_create_subarray(1, &totsizej, &subsizej, &substartj, MPI_ORDER_C, MPI_DOUBLE, &subj);
    MPI_Type_commit(&subj);

    // the lines below describe the array in case transposes are not used before saving
    // int totsize [3] = {kmax, jtot, itot};
    // int subsize [3] = {kmax, jmax, imax};
    // int substart[3] = {0, master->mpicoordy*jmax, master->mpicoordx*imax};
    int totsize [3] = {kmax  , jtot, itot};
    int subsize [3] = {kblock, jmax, itot};
    int substart[3] = {master->mpicoordx*kblock, master->mpicoordy*jmax, 0};
    MPI_Type_create_subarray(3, totsize, subsize, substart, MPI_ORDER_C, MPI_DOUBLE, &subarray);
    MPI_Type_commit(&subarray);

    // save mpitype for a xz-slice for cross section processing
    int totxzsize [2] = {kmax, itot};
    int subxzsize [2] = {kmax, imax};
    int subxzstart[2] = {0, master->mpicoordx*imax};
    MPI_Type_create_subarray(2, totxzsize, subxzsize, subxzstart, MPI_ORDER_C, MPI_DOUBLE, &subxzslice);
    MPI_Type_commit(&subxzslice);
    
    // save mpitype for a yz-slice for cross section processing
    int totyzsize [2] = {kmax, jtot};
    int subyzsize [2] = {kmax, jmax};
    int subyzstart[2] = {0, master->mpicoordy*jmax};
    MPI_Type_create_subarray(2, totyzsize, subyzsize, subyzstart, MPI_ORDER_C, MPI_DOUBLE, &subyzslice);
    MPI_Type_commit(&subyzslice);

    // save mpitype for a xy-slice for cross section processing
    int totxysize [2] = {jtot, itot};
    int subxysize [2] = {jmax, imax};
    int subxystart[2] = {master->mpicoordy*jmax, master->mpicoordx*imax};
    MPI_Type_create_subarray(2, totxysize, subxysize, subxystart, MPI_ORDER_C, MPI_DOUBLE, &subxyslice);
    MPI_Type_commit(&subxyslice);

    // allocate the array for the profiles
    profl = new double[kcells];

    mpitypes = true;
} 

void Grid::exit_mpi()
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

        delete[] profl;
    }
}

void Grid::boundary_cyclic(double* restrict data, Edge edge)
{
    const int ncount = 1;

    if (edge == East_west_edge || edge == Both_edges)
    {
        // Communicate east-west edges.
        const int eastout = iend-igc;
        const int westin  = 0;
        const int westout = istart;
        const int eastin  = iend;

        // Send and receive the ghost cells in east-west direction.
        MPI_Isend(&data[eastout], ncount, eastwestedge, master->neast, 1, master->commxy, &master->reqs[master->reqsn]);
        master->reqsn++;
        MPI_Irecv(&data[westin], ncount, eastwestedge, master->nwest, 1, master->commxy, &master->reqs[master->reqsn]);
        master->reqsn++;
        MPI_Isend(&data[westout], ncount, eastwestedge, master->nwest, 2, master->commxy, &master->reqs[master->reqsn]);
        master->reqsn++;
        MPI_Irecv(&data[eastin], ncount, eastwestedge, master->neast, 2, master->commxy, &master->reqs[master->reqsn]);
        master->reqsn++;
        // Wait here for the MPI to have correct values in the corners of the cells.
        master->wait_all();
    }

    if (edge == North_south_edge || edge == Both_edges)
    {
        // If the run is 3D, perform the cyclic boundary routine for the north-south direction.
        if (jtot > 1)
        {
            // Communicate north-south edges.
            const int northout = (jend-jgc)*icells;
            const int southin  = 0;
            const int southout = jstart*icells;
            const int northin  = jend  *icells;

            // Send and receive the ghost cells in the north-south direction.
            MPI_Isend(&data[northout], ncount, northsouthedge, master->nnorth, 1, master->commxy, &master->reqs[master->reqsn]);
            master->reqsn++;
            MPI_Irecv(&data[southin], ncount, northsouthedge, master->nsouth, 1, master->commxy, &master->reqs[master->reqsn]);
            master->reqsn++;
            MPI_Isend(&data[southout], ncount, northsouthedge, master->nsouth, 2, master->commxy, &master->reqs[master->reqsn]);
            master->reqsn++;
            MPI_Irecv(&data[northin], ncount, northsouthedge, master->nnorth, 2, master->commxy, &master->reqs[master->reqsn]);
            master->reqsn++;
            master->wait_all();
        }
        // In case of 2D, fill all the ghost cells in the y-direction with the same value.
        else
        {
            const int jj = icells;
            const int kk = icells*jcells;

            for (int k=kstart; k<kend; k++)
                for (int j=0; j<jgc; j++)
#pragma ivdep
                    for (int i=0; i<icells; i++)
                    {
                        const int ijkref   = i + jstart*jj   + k*kk;
                        const int ijknorth = i + j*jj        + k*kk;
                        const int ijksouth = i + (jend+j)*jj + k*kk;
                        data[ijknorth] = data[ijkref];
                        data[ijksouth] = data[ijkref];
                    }
        }
    }
}

void Grid::boundary_cyclic_2d(double* restrict data)
{
    int ncount = 1;

    // communicate east-west edges
    int eastout = iend-igc;
    int westin  = 0;
    int westout = istart;
    int eastin  = iend;

    // communicate north-south edges
    int northout = (jend-jgc)*icells;
    int southin  = 0;
    int southout = jstart*icells;
    int northin  = jend  *icells;

    // first, send and receive the ghost cells in east-west direction
    MPI_Isend(&data[eastout], ncount, eastwestedge2d, master->neast, 1, master->commxy, &master->reqs[master->reqsn]);
    master->reqsn++;
    MPI_Irecv(&data[westin], ncount, eastwestedge2d, master->nwest, 1, master->commxy, &master->reqs[master->reqsn]);
    master->reqsn++;
    MPI_Isend(&data[westout], ncount, eastwestedge2d, master->nwest, 2, master->commxy, &master->reqs[master->reqsn]);
    master->reqsn++;
    MPI_Irecv(&data[eastin], ncount, eastwestedge2d, master->neast, 2, master->commxy, &master->reqs[master->reqsn]);
    master->reqsn++;
    // wait here for the mpi to have correct values in the corners of the cells
    master->wait_all();

    // if the run is 3D, apply the BCs
    if (jtot > 1)
    {
        // second, send and receive the ghost cells in the north-south direction
        MPI_Isend(&data[northout], ncount, northsouthedge2d, master->nnorth, 1, master->commxy, &master->reqs[master->reqsn]);
        master->reqsn++;
        MPI_Irecv(&data[southin], ncount, northsouthedge2d, master->nsouth, 1, master->commxy, &master->reqs[master->reqsn]);
        master->reqsn++;
        MPI_Isend(&data[southout], ncount, northsouthedge2d, master->nsouth, 2, master->commxy, &master->reqs[master->reqsn]);
        master->reqsn++;
        MPI_Irecv(&data[northin], ncount, northsouthedge2d, master->nnorth, 2, master->commxy, &master->reqs[master->reqsn]);
        master->reqsn++;
        master->wait_all();
    }
    // in case of 2D, fill all the ghost cells with the current value
    else
    {
        // 2d essential variables
        const int jj = icells;

        for (int j=0; j<jgc; j++)
#pragma ivdep
            for (int i=0; i<icells; i++)
            {
                const int ijref   = i + jstart*jj;
                const int ijnorth = i + j*jj;
                const int ijsouth = i + (jend+j)*jj;
                data[ijnorth] = data[ijref];
                data[ijsouth] = data[ijref];
            }
    }
}

void Grid::transpose_zx(double* restrict ar, double* restrict as)
{
    const int ncount = 1;
    const int tag = 1;

    const int jj = imax;
    const int kk = imax*jmax;

    for (int n=0; n<master->npx; n++)
    {
        // determine where to fetch the data and where to store it
        const int ijks = n*kblock*kk;
        const int ijkr = n*jj;

        // send and receive the data
        MPI_Isend(&as[ijks], ncount, transposez, n, tag, master->commx, &master->reqs[master->reqsn]);
        master->reqsn++;
        MPI_Irecv(&ar[ijkr], ncount, transposex, n, tag, master->commx, &master->reqs[master->reqsn]);
        master->reqsn++;
    }
    master->wait_all();
}

void Grid::transpose_xz(double* restrict ar, double* restrict as)
{
    const int ncount = 1;
    const int tag = 1;

    const int jj = imax;
    const int kk = imax*jmax;

    for (int n=0; n<master->npx; n++)
    {
        // determine where to fetch the data and where to store it
        const int ijks = n*jj;
        const int ijkr = n*kblock*kk;

        // send and receive the data
        MPI_Isend(&as[ijks], ncount, transposex, n, tag, master->commx, &master->reqs[master->reqsn]);
        master->reqsn++;
        MPI_Irecv(&ar[ijkr], ncount, transposez, n, tag, master->commx, &master->reqs[master->reqsn]);
        master->reqsn++;
    }
    master->wait_all();
}

void Grid::transpose_xy(double* restrict ar, double* restrict as)
{
    const int ncount = 1;
    const int tag = 1;

    const int jj = iblock;
    const int kk = iblock*jmax;

    for (int n=0; n<master->npy; n++)
    {
        // determine where to fetch the data and where to store it
        const int ijks = n*jj;
        const int ijkr = n*kk;

        // send and receive the data
        MPI_Isend(&as[ijks], ncount, transposex2, n, tag, master->commy, &master->reqs[master->reqsn]);
        master->reqsn++;
        MPI_Irecv(&ar[ijkr], ncount, transposey , n, tag, master->commy, &master->reqs[master->reqsn]);
        master->reqsn++;
    }
    master->wait_all();
}

void Grid::transpose_yx(double* restrict ar, double* restrict as)
{
    const int ncount = 1;
    const int tag = 1;

    const int jj = iblock;
    const int kk = iblock*jmax;

    for (int n=0; n<master->npy; n++)
    {
        // determine where to fetch the data and where to store it
        const int ijks = n*kk;
        const int ijkr = n*jj;

        // send and receive the data
        MPI_Isend(&as[ijks], ncount, transposey , n, tag, master->commy, &master->reqs[master->reqsn]);
        master->reqsn++;
        MPI_Irecv(&ar[ijkr], ncount, transposex2, n, tag, master->commy, &master->reqs[master->reqsn]);
        master->reqsn++;
    }
    master->wait_all();
}

void Grid::transpose_yz(double* restrict ar, double* restrict as)
{
    const int ncount = 1;
    const int tag = 1;

    const int jj = iblock;
    const int kk = iblock*jblock;

    for (int n=0; n<master->npx; n++)
    {
        // determine where to fetch the data and where to store it
        const int ijks = n*jblock*jj;
        const int ijkr = n*kblock*kk;

        // send and receive the data
        MPI_Isend(&as[ijks], ncount, transposey2, n, tag, master->commx, &master->reqs[master->reqsn]);
        master->reqsn++;
        MPI_Irecv(&ar[ijkr], ncount, transposez2, n, tag, master->commx, &master->reqs[master->reqsn]);
        master->reqsn++;
    }
    master->wait_all();
}

void Grid::transpose_zy(double* restrict ar, double* restrict as)
{
    const int ncount = 1;
    const int tag = 1;

    const int jj = iblock;
    const int kk = iblock*jblock;

    for (int n=0; n<master->npx; n++)
    {
        // determine where to fetch the data and where to store it
        const int ijks = n*kblock*kk;
        const int ijkr = n*jblock*jj;

        // send and receive the data
        MPI_Isend(&as[ijks], ncount, transposez2, n, tag, master->commx, &master->reqs[master->reqsn]);
        master->reqsn++;
        MPI_Irecv(&ar[ijkr], ncount, transposey2, n, tag, master->commx, &master->reqs[master->reqsn]);
        master->reqsn++;
    }
    master->wait_all();
}

void Grid::get_max(double *var)
{
    double varl = *var;
    MPI_Allreduce(&varl, var, 1, MPI_DOUBLE, MPI_MAX, master->commxy);
}

void Grid::get_sum(double *var)
{
    double varl = *var;
    MPI_Allreduce(&varl, var, 1, MPI_DOUBLE, MPI_SUM, master->commxy);
}

void Grid::get_prof(double *prof, int kcellsin)
{
    for (int k=0; k<kcellsin; k++)
        profl[k] = prof[k] / master->nprocs;

    MPI_Allreduce(profl, prof, kcellsin, MPI_DOUBLE, MPI_SUM, master->commxy);
}

// IO functions
void Grid::save()
{
    char filename[256];
    std::sprintf(filename, "%s.%07d", "grid", 0);
    master->print_message("Saving \"%s\" ... ", filename);

    MPI_File fh;
    if (MPI_File_open(master->commxy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
    {
        master->print_message("FAILED\n");
        throw 1;
    }

    // select noncontiguous part of 3d array to store the selected data
    MPI_Offset fileoff = 0; // the offset within the file (header size)
    char name[] = "native";

    MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subi, name, MPI_INFO_NULL);
    if (master->mpicoordy == 0)
        MPI_File_write(fh, &x[istart], imax, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_Barrier(master->commxy);
    fileoff += itot*sizeof(double);

    MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subi, name, MPI_INFO_NULL);
    if (master->mpicoordy == 0)
        MPI_File_write(fh, &xh[istart], imax, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_Barrier(master->commxy);
    fileoff += itot*sizeof(double);

    MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subj, name, MPI_INFO_NULL);
    if (master->mpicoordx == 0)
        MPI_File_write(fh, &y[jstart], jmax, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_Barrier(master->commxy);
    fileoff += jtot*sizeof(double);

    MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subj, name, MPI_INFO_NULL);
    if (master->mpicoordx == 0)
        MPI_File_write(fh, &yh[jstart], jmax, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_Barrier(master->commxy);

    MPI_File_sync(fh);
    if (MPI_File_close(&fh))
        throw 1;

    if (master->mpiid == 0)
    {
        FILE *pFile;
        pFile = fopen(filename, "ab");
        fwrite(&z [kstart], sizeof(double), kmax, pFile);
        fwrite(&zh[kstart], sizeof(double), kmax, pFile);
        fclose(pFile);
    }

    // the saving procedure is a success
    master->print_message("OK\n");

    // SAVE THE FFTW PLAN IN ORDER TO ENSURE BITWISE IDENTICAL RESTARTS
    // use the FFTW3 many interface in order to reduce function call overhead
    int rank = 1;
    int ni[] = {itot};
    int nj[] = {jtot};
    int istride = 1;
    int jstride = iblock;
    int idist = itot;
    int jdist = 1;

    fftw_r2r_kind kindf[] = {FFTW_R2HC};
    fftw_r2r_kind kindb[] = {FFTW_HC2R};

    iplanf = fftw_plan_many_r2r(rank, ni, jmax, fftini, ni, istride, idist,
                                fftouti, ni, istride, idist, kindf, FFTW_EXHAUSTIVE);
    iplanb = fftw_plan_many_r2r(rank, ni, jmax, fftini, ni, istride, idist,
                                fftouti, ni, istride, idist, kindb, FFTW_EXHAUSTIVE);
    jplanf = fftw_plan_many_r2r(rank, nj, iblock, fftinj, nj, jstride, jdist,
                                fftoutj, nj, jstride, jdist, kindf, FFTW_EXHAUSTIVE);
    jplanb = fftw_plan_many_r2r(rank, nj, iblock, fftinj, nj, jstride, jdist,
                                fftoutj, nj, jstride, jdist, kindb, FFTW_EXHAUSTIVE);

    fftwplan = true;

    if (master->mpiid == 0)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", "fftwplan", 0);

        master->print_message("Saving \"%s\" ... ", filename);

        int n = fftw_export_wisdom_to_filename(filename);
        if (n == 0)
        {
            master->print_error("\"%s\" cannot be saved\n", filename);
            throw 1;
        }
        else
            master->print_message("OK\n");
    }
}

void Grid::load()
{
    int nerror = 0;

    // LOAD THE GRID
    char filename[256];
    std::sprintf(filename, "%s.%07d", "grid", 0);
    if (master->mpiid == 0) std::printf("Loading \"%s\" ... ", filename);

    FILE *pFile;
    if (master->mpiid == 0)
    {
        pFile = fopen(filename, "rb");
        if (pFile == NULL)
        {
            ++nerror;
        }
        else
        {
            int n = (2*itot+2*jtot)*sizeof(double);
            fseek(pFile, n, SEEK_SET);
            fread(&z [kstart], sizeof(double), kmax, pFile);
            fread(&zh[kstart], sizeof(double), kmax, pFile);
            fclose(pFile);
        }
    }

    // communicate the file read error over all procs
    master->broadcast(&nerror, 1);
    if (nerror)
    {
        master->print_message("FAILED\n");
        throw 1;
    }
    else
        master->print_message("OK\n");

    master->broadcast(&z [kstart], kmax);
    master->broadcast(&zh[kstart], kmax);

    // calculate the missing coordinates
    calculate();

    // LOAD THE FFTW PLAN
    std::sprintf(filename, "%s.%07d", "fftwplan", 0);

    master->print_message("Loading \"%s\" ... ", filename);

    int n = fftw_import_wisdom_from_filename(filename);
    if (n == 0)
    {
        master->print_message("FAILED\n");
        throw 1;
    }
    else
        master->print_message("OK\n");

    // use the FFTW3 many interface in order to reduce function call overhead
    int rank = 1;
    int ni[] = {itot};
    int nj[] = {jtot};
    int istride = 1;
    int jstride = iblock;
    int idist = itot;
    int jdist = 1;
    fftw_r2r_kind kindf[] = {FFTW_R2HC};
    fftw_r2r_kind kindb[] = {FFTW_HC2R};
    iplanf = fftw_plan_many_r2r(rank, ni, jmax, fftini, ni, istride, idist,
                                fftouti, ni, istride, idist, kindf, FFTW_EXHAUSTIVE);
    iplanb = fftw_plan_many_r2r(rank, ni, jmax, fftini, ni, istride, idist,
                                fftouti, ni, istride, idist, kindb, FFTW_EXHAUSTIVE);
    jplanf = fftw_plan_many_r2r(rank, nj, iblock, fftinj, nj, jstride, jdist,
                                fftoutj, nj, jstride, jdist, kindf, FFTW_EXHAUSTIVE);
    jplanb = fftw_plan_many_r2r(rank, nj, iblock, fftinj, nj, jstride, jdist,
                                fftoutj, nj, jstride, jdist, kindb, FFTW_EXHAUSTIVE);

    fftwplan = true;

    fftw_forget_wisdom();
}

int Grid::save_field3d(double* restrict data, double* restrict tmp1, double* restrict tmp2, char* filename, double offset)
{
    // save the data in transposed order to have large chunks of contiguous disk space
    // MPI-IO is not stable on Juqueen and supermuc otherwise

    // extract the data from the 3d field without the ghost cells
    const int jj  = icells;
    const int kk  = icells*jcells;
    const int jjb = imax;
    const int kkb = imax*jmax;

    int count = imax*jmax*kmax;

    for (int k=0; k<kmax; k++)
        for (int j=0; j<jmax; j++)
#pragma ivdep
            for (int i=0; i<imax; i++)
            {
                const int ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;
                const int ijkb = i + j*jjb + k*kkb;
                tmp1[ijkb] = data[ijk] + offset;
            }

    transpose_zx(tmp2, tmp1);

    MPI_File fh;
    if (MPI_File_open(master->commxy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
        return 1;

    // select noncontiguous part of 3d array to store the selected data
    MPI_Offset fileoff = 0; // the offset within the file (header size)
    char name[] = "native";

    if (MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subarray, name, MPI_INFO_NULL))
        return 1;

    if (MPI_File_write_all(fh, tmp2, count, MPI_DOUBLE, MPI_STATUS_IGNORE))
        return 1;

    if (MPI_File_close(&fh))
        return 1;

    return 0;
}

int Grid::load_field3d(double* restrict data, double* restrict tmp1, double* restrict tmp2, char* filename, double offset)
{
    // save the data in transposed order to have large chunks of contiguous disk space
    // MPI-IO is not stable on Juqueen and supermuc otherwise

    // read the file
    MPI_File fh;
    if (MPI_File_open(master->commxy, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
        return 1;

    // select noncontiguous part of 3d array to store the selected data
    MPI_Offset fileoff = 0; // the offset within the file (header size)
    char name[] = "native";
    MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subarray, name, MPI_INFO_NULL);

    // extract the data from the 3d field without the ghost cells
    int count = imax*jmax*kmax;

    if (MPI_File_read_all(fh, tmp1, count, MPI_DOUBLE, MPI_STATUS_IGNORE))
        return 1;

    if (MPI_File_close(&fh))
        return 1;

    // transpose the data back
    transpose_xz(tmp2, tmp1);

    const int jj  = icells;
    const int kk  = icells*jcells;
    const int jjb = imax;
    const int kkb = imax*jmax;

    for (int k=0; k<kmax; k++)
        for (int j=0; j<jmax; j++)
#pragma ivdep
            for (int i=0; i<imax; i++)
            {
                const int ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;
                const int ijkb = i + j*jjb + k*kkb;
                data[ijk] = tmp2[ijkb] - offset;
            }

    return 0;
}

void Grid::fft_forward(double* restrict data,   double* restrict tmp1,
                       double* restrict fftini, double* restrict fftouti,
                       double* restrict fftinj, double* restrict fftoutj)
{
    // transpose the pressure field
    transpose_zx(tmp1,data);

    int kk = itot*jmax;

    // process the fourier transforms slice by slice
    for (int k=0; k<kblock; k++)
    {
#pragma ivdep
        for (int n=0; n<itot*jmax; n++)
        {
            const int ij  = n;
            const int ijk = n + k*kk;
            fftini[ij] = tmp1[ijk];
        }

        fftw_execute(iplanf);

#pragma ivdep
        for (int n=0; n<itot*jmax; n++)
        {
            const int ij  = n;
            const int ijk = n + k*kk;
            tmp1[ijk] = fftouti[ij];
        }
    }

    // transpose again
    transpose_xy(data,tmp1);

    kk = iblock*jtot;

    // do the second fourier transform
    for (int k=0; k<kblock; k++)
    {
#pragma ivdep
        for (int n=0; n<iblock*jtot; n++)
        {
            const int ij  = n;
            const int ijk = n + k*kk;
            fftinj[ij] = data[ijk];
        }

        fftw_execute(jplanf);

#pragma ivdep
        for (int n=0; n<iblock*jtot; n++)
        {
            const int ij  = n;
            const int ijk = n + k*kk;
            // shift to use p in pressure solver
            tmp1[ijk] = fftoutj[ij];
        }
    }

    // transpose back to original orientation
    transpose_yz(data,tmp1);
}

void Grid::fft_backward(double* restrict data,   double* restrict tmp1,
                        double* restrict fftini, double* restrict fftouti,
                        double* restrict fftinj, double* restrict fftoutj)
{
    // transpose back to y
    transpose_zy(tmp1, data);

    int kk = iblock*jtot;

    // transform the second transform back
    for (int k=0; k<kblock; k++)
    {
#pragma ivdep
        for (int n=0; n<iblock*jtot; n++)
        {
            const int ij  = n;
            const int ijk = n + k*kk;
            fftinj[ij] = tmp1[ijk];
        }

        fftw_execute(jplanb);

#pragma ivdep
        for (int n=0; n<iblock*jtot; n++)
        {
            const int ij  = n;
            const int ijk = n + k*kk;
            data[ijk] = fftoutj[ij] / jtot;
        }
    }

    // transpose back to x
    transpose_yx(tmp1, data);

    kk = itot*jmax;

    // transform the first transform back
    for (int k=0; k<kblock; k++)
    {
#pragma ivdep
        for (int n=0; n<itot*jmax; n++)
        {
            const int ij  = n;
            const int ijk = n + k*kk;
            fftini[ij] = tmp1[ijk];
        }

        fftw_execute(iplanb);

#pragma ivdep
        for (int n=0; n<itot*jmax; n++)
        {
            const int ij  = n;
            const int ijk = n + k*kk;
            // swap array here to avoid unnecessary 3d loop
            data[ijk] = fftouti[ij] / itot;
        }
    }

    // and transpose back...
    transpose_xz(tmp1, data);
}

int Grid::save_xz_slice(double* restrict data, double* restrict tmp, char* filename, int jslice)
{
    // extract the data from the 3d field without the ghost cells
    int nerror=0;

    const int jj  = icells;
    const int kk  = icells*jcells;
    const int kkb = imax;

    int count = imax*kmax;

    for (int k=0; k<kmax; k++)
#pragma ivdep
        for (int i=0; i<imax; i++)
        {
            // take the modulus of jslice and jmax to have the right offset within proc
            const int ijk  = i+igc + ((jslice%jmax)+jgc)*jj + (k+kgc)*kk;
            const int ijkb = i + k*kkb;
            tmp[ijkb] = data[ijk];
        }

    if (master->mpicoordy == jslice/jmax)
    {
        MPI_File fh;
        if (MPI_File_open(master->commx, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
            ++nerror;

        // select noncontiguous part of 3d array to store the selected data
        MPI_Offset fileoff = 0; // the offset within the file (header size)
        char name[] = "native";

        if (!nerror)
            if (MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subxzslice, name, MPI_INFO_NULL))
                ++nerror;

        // only write at the procs that contain the slice
        if (!nerror)
            if (MPI_File_write_all(fh, tmp, count, MPI_DOUBLE, MPI_STATUS_IGNORE))
                ++nerror;

        if (!nerror)
            MPI_File_sync(fh);

        if (!nerror)
            if (MPI_File_close(&fh))
                ++nerror;
    }

    // Gather errors from other processes
    master->sum(&nerror,1);

    MPI_Barrier(master->commxy);

    return nerror;
}

int Grid::save_yz_slice(double* restrict data, double* restrict tmp, char* filename, int islice)
{
    // extract the data from the 3d field without the ghost cells
    int nerror=0;

    const int ii  = jcells;
    const int kk  = icells*jcells;
    const int kkb = jmax;

    int count = jmax*kmax;

    for (int k=0; k<kmax; k++)
#pragma ivdep
        for (int j=0; j<jmax; j++)
        {
            // take the modulus of jslice and jmax to have the right offset within proc
            const int ijk  = j+jgc + ((islice%imax)+igc)*ii + (k+kgc)*kk;
            const int ijkb = j + k*kkb;
            tmp[ijkb] = data[ijk];
        }

    if (master->mpicoordx == islice/imax)
    {
        MPI_File fh;
        if (MPI_File_open(master->commy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
            ++nerror;

        // select noncontiguous part of 3d array to store the selected data
        MPI_Offset fileoff = 0; // the offset within the file (header size)
        char name[] = "native";

        if (!nerror)
            if (MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subyzslice, name, MPI_INFO_NULL))
                ++nerror;

        // only write at the procs that contain the slice
        if (!nerror)
            if (MPI_File_write_all(fh, tmp, count, MPI_DOUBLE, MPI_STATUS_IGNORE))
                ++nerror;

        if (!nerror)
            MPI_File_sync(fh);

        if (!nerror)
            if (MPI_File_close(&fh))
                ++nerror;
    }

    // Gather errors from other processes
    master->sum(&nerror,1);

    MPI_Barrier(master->commxy);

    return nerror;
}

int Grid::save_xy_slice(double* restrict data, double* restrict tmp, char* filename, int kslice)
{
    // extract the data from the 3d field without the ghost cells
    const int jj  = icells;
    const int kk  = icells*jcells;
    const int jjb = imax;

    // Subtract the ghost cells in case of a pure 2d plane that does not have ghost cells.
    if (kslice == -1)
        kslice = -kgc;

    int count = imax*jmax;

    for (int j=0; j<jmax; j++)
#pragma ivdep
        for (int i=0; i<imax; i++)
        {
            // take the modulus of jslice and jmax to have the right offset within proc
            const int ijk  = i+igc + (j+jgc)*jj + (kslice+kgc)*kk;
            const int ijkb = i + j*jjb;
            tmp[ijkb] = data[ijk];
        }

    MPI_File fh;
    if (MPI_File_open(master->commxy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
        return 1;

    // select noncontiguous part of 3d array to store the selected data
    MPI_Offset fileoff = 0; // the offset within the file (header size)
    char name[] = "native";

    if (MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subxyslice, name, MPI_INFO_NULL))
        return 1;

    // only write at the procs that contain the slice
    if (MPI_File_write_all(fh, tmp, count, MPI_DOUBLE, MPI_STATUS_IGNORE))
        return 1;

    MPI_File_sync(fh);

    if (MPI_File_close(&fh))
        return 1;

    MPI_Barrier(master->commxy);

    return 0;
}

int Grid::load_xy_slice(double* restrict data, double* restrict tmp, char* filename, int kslice)
{
    // extract the data from the 3d field without the ghost cells
    const int jj  = icells;
    const int kk  = icells*jcells;
    const int jjb = imax;

    // Subtract the ghost cells in case of a pure 2d plane that does not have ghost cells.
    if (kslice == -1)
        kslice = -kgc;

    int count = imax*jmax;

    MPI_File fh;
    if (MPI_File_open(master->commxy, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh))
        return 1;

    // select noncontiguous part of 3d array to store the selected data
    MPI_Offset fileoff = 0; // the offset within the file (header size)
    char name[] = "native";

    if (MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subxyslice, name, MPI_INFO_NULL))
        return 1;

    // only write at the procs that contain the slice
    if (MPI_File_read_all(fh, tmp, count, MPI_DOUBLE, MPI_STATUS_IGNORE))
        return 1;

    if (MPI_File_close(&fh))
        return 1;

    MPI_Barrier(master->commxy);

    for (int j=0; j<jmax; j++)
#pragma ivdep
        for (int i=0; i<imax; i++)
        {
            // take the modulus of jslice and jmax to have the right offset within proc
            const int ijk  = i+igc + (j+jgc)*jj + (kslice+kgc)*kk;
            const int ijkb = i + j*jjb;
            data[ijk] = tmp[ijkb];
        }

    return 0;
}
#endif
