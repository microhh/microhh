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

#ifndef USEMPI
#include <fftw3.h>
#include <cstdio>
#include "master.h"
#include "grid.h"
#include "defines.h"

// MPI functions
void Grid::init_mpi()
{
    mpitypes = true;
} 

void Grid::exit_mpi()
{
}

void Grid::boundary_cyclic(double* restrict data, Edge edge)
{
    const int jj = icells;
    const int kk = icells*jcells;

    if (edge == East_west_edge || edge == Both_edges)
    {
        // first, east west boundaries
        for (int k=0; k<kcells; k++)
            for (int j=0; j<jcells; j++)
#pragma ivdep
                for (int i=0; i<igc; i++)
                {
                    const int ijk0 = i          + j*jj + k*kk;
                    const int ijk1 = iend-igc+i + j*jj + k*kk;
                    data[ijk0] = data[ijk1];
                }

        for (int k=0; k<kcells; k++)
            for (int j=0; j<jcells; j++)
#pragma ivdep
                for (int i=0; i<igc; i++)
                {
                    const int ijk0 = i+iend   + j*jj + k*kk;
                    const int ijk1 = i+istart + j*jj + k*kk;
                    data[ijk0] = data[ijk1];
                }
    }

    if (edge == North_south_edge || edge == Both_edges)
    {
        // if the run is 3D, apply the BCs
        if (jtot > 1)
        {
            // second, send and receive the ghost cells in the north-south direction
            for (int k=0; k<kcells; k++)
                for (int j=0; j<jgc; j++)
#pragma ivdep
                    for (int i=0; i<icells; i++)
                    {
                        const int ijk0 = i + j           *jj + k*kk;
                        const int ijk1 = i + (jend-jgc+j)*jj + k*kk;
                        data[ijk0] = data[ijk1];
                    }

            for (int k=0; k<kcells; k++)
                for (int j=0; j<jgc; j++)
#pragma ivdep
                    for (int i=0; i<icells; i++)
                    {
                        const int ijk0 = i + (j+jend  )*jj + k*kk;
                        const int ijk1 = i + (j+jstart)*jj + k*kk;
                        data[ijk0] = data[ijk1];
                    }
        }
        // in case of 2D, fill all the ghost cells with the current value
        else
        {
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
    const int jj = icells;

    // first, east west boundaries
    for (int j=0; j<jcells; j++)
#pragma ivdep
        for (int i=0; i<igc; i++)
        {
            const int ij0 = i          + j*jj;
            const int ij1 = iend-igc+i + j*jj;
            data[ij0] = data[ij1];
        }

    for (int j=0; j<jcells; j++)
#pragma ivdep
        for (int i=0; i<igc; i++)
        {
            const int ij0 = i+iend   + j*jj;
            const int ij1 = i+istart + j*jj;
            data[ij0] = data[ij1];
        }

    // if the run is 3D, apply the BCs
    if (jtot > 1)
    {
        // second, send and receive the ghost cells in the north-south direction
        for (int j=0; j<jgc; j++)
#pragma ivdep
            for (int i=0; i<icells; i++)
            {
                const int ij0 = i + j           *jj;
                const int ij1 = i + (jend-jgc+j)*jj;
                data[ij0] = data[ij1];
            }

        for (int j=0; j<jgc; j++)
#pragma ivdep
            for (int i=0; i<icells; i++)
            {
                const int ij0 = i + (j+jend  )*jj;
                const int ij1 = i + (j+jstart)*jj;
                data[ij0] = data[ij1];
            }
    }
    // in case of 2D, fill all the ghost cells with the current value
    else
    {
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
    const int jj = imax;
    const int kk = imax*jmax;

    for (int k=0; k<kmax; k++)
        for (int j=0; j<jmax; j++)
#pragma ivdep
            for (int i=0; i<imax; i++)
            {
                const int ijk = i + j*jj + k*kk;
                ar[ijk] = as[ijk];
            }
}

void Grid::transpose_xz(double* restrict ar, double* restrict as)
{
    const int jj = imax;
    const int kk = imax*jmax;

    for (int k=0; k<kmax; k++)
        for (int j=0; j<jmax; j++)
#pragma ivdep
            for (int i=0; i<imax; i++)
            {
                const int ijk = i + j*jj + k*kk;
                ar[ijk] = as[ijk];
            }
}

void Grid::transpose_xy(double* restrict ar, double* restrict as)
{
    const int jj = imax;
    const int kk = imax*jmax;

    for (int k=0; k<kmax; k++)
        for (int j=0; j<jmax; j++)
#pragma ivdep
            for (int i=0; i<imax; i++)
            {
                const int ijk = i + j*jj + k*kk;
                ar[ijk] = as[ijk];
            }
}

void Grid::transpose_yx(double* restrict ar, double* restrict as)
{
    const int jj = imax;
    const int kk = imax*jmax;

    for (int k=0; k<kmax; k++)
        for (int j=0; j<jmax; j++)
#pragma ivdep
            for (int i=0; i<imax; i++)
            {
                const int ijk = i + j*jj + k*kk;
                ar[ijk] = as[ijk];
            }
}

void Grid::transpose_yz(double* restrict ar, double* restrict as)
{
    const int jj = imax;
    const int kk = imax*jmax;

    for (int k=0; k<kmax; k++)
        for (int j=0; j<jmax; j++)
#pragma ivdep
            for (int i=0; i<imax; i++)
            {
                const int ijk = i + j*jj + k*kk;
                ar[ijk] = as[ijk];
            }
}

void Grid::transpose_zy(double* restrict ar, double* restrict as)
{
    const int jj = imax;
    const int kk = imax*jmax;

    for (int k=0; k<kmax; k++)
        for (int j=0; j<jmax; j++)
#pragma ivdep
            for (int i=0; i<imax; i++)
            {
                const int ijk = i + j*jj + k*kk;
                ar[ijk] = as[ijk];
            }
}

void Grid::get_max(double *var)
{
}

void Grid::get_max(int *var)
{
}

void Grid::get_sum(double *var)
{
}

void Grid::get_prof(double *prof, int kcellsin)
{
}

// IO functions
void Grid::save()
{
    // SAVE THE GRID
    FILE *pFile;
    char filename[256];
    std::sprintf(filename, "%s.%07d", "grid", 0);
    pFile = fopen(filename, "wbx");
    master->print_message("Saving \"%s\" ... ", filename);

    if (pFile == NULL)
    {
        master->print_message("FAILED\n");
        throw 1;
    }
    else
        master->print_message("OK\n");

    fwrite(&x [istart], sizeof(double), itot, pFile);
    fwrite(&xh[istart], sizeof(double), itot, pFile);
    fwrite(&y [jstart], sizeof(double), jtot, pFile);
    fwrite(&yh[jstart], sizeof(double), jtot, pFile);
    fwrite(&z [kstart], sizeof(double), ktot, pFile);
    fwrite(&zh[kstart], sizeof(double), ktot, pFile);
    fclose(pFile);

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
            master->print_message("FAILED\n");
            throw 1;
        }
        else
            master->print_message("OK\n");
    }
}

void Grid::load()
{
    // LOAD THE GRID
    FILE *pFile;
    char filename[256];
    std::sprintf(filename, "%s.%07d", "grid", 0);
    pFile = fopen(filename, "rb");
    master->print_message("Loading \"%s\" ... ", filename);

    if (pFile == NULL)
    {
        master->print_message("FAILED\n");
        throw 1;
    }
    else
        master->print_message("OK\n");

    fread(&x [istart], sizeof(double), itot, pFile);
    fread(&xh[istart], sizeof(double), itot, pFile);
    fread(&y [jstart], sizeof(double), jtot, pFile);
    fread(&yh[jstart], sizeof(double), jtot, pFile);
    fread(&z [kstart], sizeof(double), ktot, pFile);
    fread(&zh[kstart], sizeof(double), ktot, pFile);
    fclose(pFile);

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
    FILE *pFile;
    pFile = fopen(filename, "wbx");

    if (pFile == NULL)
        return 1;

    const int jj = icells;
    const int kk = icells*jcells;

    // first, add the offset to the data
    for (int k=kstart; k<kend; k++)
        for (int j=jstart; j<jend; j++)
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                tmp1[ijk] = data[ijk] + offset;
            }

    // second, save the data to disk
    for (int k=kstart; k<kend; k++)
        for (int j=jstart; j<jend; j++)
        {
            const int ijk = istart + j*jj + k*kk;
            fwrite(&tmp1[ijk], sizeof(double), imax, pFile);
        }

    fclose(pFile);

    return 0;
}

int Grid::load_field3d(double* restrict data, double* restrict tmp1, double* restrict tmp2, char* filename, double offset)
{
    FILE *pFile;
    pFile = fopen(filename, "rb");

    if (pFile == NULL)
        return 1;

    const int jj = icells;
    const int kk = icells*jcells;

    // first, load the data from disk
    for (int k=kstart; k<kend; k++)
        for (int j=jstart; j<jend; j++)
        {
            const int ijk = istart + j*jj + k*kk;
            fread(&tmp1[ijk], sizeof(double), imax, pFile);
        }

    fclose(pFile);

    // second, remove the offset
    for (int k=kstart; k<kend; k++)
        for (int j=jstart; j<jend; j++)
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                data[ijk] = tmp1[ijk] - offset;
            }

    return 0;
}

void Grid::fft_forward(double* restrict data,   double* restrict tmp1,
                       double* restrict fftini, double* restrict fftouti,
                       double* restrict fftinj, double* restrict fftoutj)
{
    int kk = itot*jmax;

    // process the fourier transforms slice by slice
    for (int k=0; k<kblock; k++)
    {
#pragma ivdep
        for (int n=0; n<itot*jmax; n++)
        {
            const int ij  = n;
            const int ijk = n + k*kk;
            fftini[ij] = data[ijk];
        }

        fftw_execute(iplanf);

#pragma ivdep
        for (int n=0; n<itot*jmax; n++)
        {
            const int ij  = n;
            const int ijk = n + k*kk;
            data[ijk] = fftouti[ij];
        }
    }

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
            data[ijk] = fftoutj[ij];
        }
    }
}

void Grid::fft_backward(double* restrict data,   double* restrict tmp1,
                        double* restrict fftini, double* restrict fftouti,
                        double* restrict fftinj, double* restrict fftoutj)
{
    int kk = iblock*jtot;

    // transform the second transform back
    for (int k=0; k<kblock; k++)
    {
#pragma ivdep
        for (int n=0; n<iblock*jtot; n++)
        {
            const int ij  = n;
            const int ijk = n + k*kk;
            fftinj[ij] = data[ijk];
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

    kk = itot*jmax;

    // transform the first transform back
    for (int k=0; k<kblock; k++)
    {
#pragma ivdep
        for (int n=0; n<itot*jmax; n++)
        {
            const int ij  = n;
            const int ijk = n + k*kk;
            fftini[ij] = data[ijk];
        }

        fftw_execute(iplanb);

#pragma ivdep
        for (int n=0; n<itot*jmax; n++)
        {
            const int ij  = n;
            const int ijk = n + k*kk;
            // swap array here to avoid unnecessary 3d loop
            tmp1[ijk] = fftouti[ij] / itot;
        }
    }
}

int Grid::save_xz_slice(double* restrict data, double* restrict tmp, char* filename, int jslice)
{
    // extract the data from the 3d field without the ghost cells
    const int jj  = icells;
    const int kk  = icells*jcells;
    const int kkb = imax;

    const int count = imax*kmax;

    for (int k=0; k<kmax; k++)
#pragma ivdep
        for (int i=0; i<imax; i++)
        {
            // take the modulus of jslice and jmax to have the right offset within proc
            const int ijk  = i+igc + (jslice+jgc)*jj + (k+kgc)*kk;
            const int ijkb = i + k*kkb;
            tmp[ijkb] = data[ijk];
        }

    FILE *pFile;
    pFile = fopen(filename, "wbx");
    if (pFile == NULL)
        return 1;

    fwrite(tmp, sizeof(double), count, pFile);
    fclose(pFile);

    return 0;
}

int Grid::save_xy_slice(double* restrict data, double* restrict tmp, char* filename, int kslice)
{
    // extract the data from the 3d field without the ghost cells
    const int jj  = icells;
    const int kk  = icells*jcells;
    const int jjb = imax;

    const int count = imax*jmax;

    // Subtract the ghost cells in case of a pure 2d plane that does not have ghost cells.
    if (kslice == -1)
        kslice = -kgc;

    for (int j=0; j<jmax; j++)
#pragma ivdep
        for (int i=0; i<imax; i++)
        {
            // take the modulus of jslice and jmax to have the right offset within proc
            const int ijk  = i+igc + (j+jgc)*jj + (kslice+kgc)*kk;
            const int ijkb = i + j*jjb;
            tmp[ijkb] = data[ijk];
        }

    FILE *pFile;
    pFile = fopen(filename, "wbx");
    if (pFile == NULL)
        return 1;

    fwrite(tmp, sizeof(double), count, pFile);
    fclose(pFile);

    return 0;
}

int Grid::load_xy_slice(double* restrict data, double* restrict tmp, char* filename, int kslice)
{
    const int count = imax*jmax;

    FILE *pFile;
    pFile = fopen(filename, "rb");
    if (pFile == NULL)
        return 1;

    fread(tmp, sizeof(double), count, pFile);
    fclose(pFile);

    // Subtract the ghost cells in case of a pure 2d plane that does not have ghost cells.
    if (kslice == -1)
        kslice = -kgc;

    // put the data back into a field with ghost cells
    const int jj  = icells;
    const int kk  = icells*jcells;
    const int jjb = imax;

    for (int j=0; j<jmax; j++)
#pragma ivdep
        for (int i=0; i<imax; i++)
        {
            const int ijk  = i+igc + (j+jgc)*jj + (kslice+kgc)*kk;
            const int ijkb = i + j*jjb;
            data[ijk] = tmp[ijkb];
        }

    return 0;
}
#endif
