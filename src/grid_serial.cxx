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

#ifndef USEMPI
#include <fftw3.h>
#include <cstdio>
#include "master.h"
#include "grid.h"
#include "defines.h"

// MPI functions
template<typename TF>
void Grid<TF>::init_mpi()
{
    mpitypes = true;
} 

template<typename TF>
void Grid<TF>::exit_mpi()
{
}

template<typename TF>
void Grid<TF>::boundary_cyclic(TF* restrict data, Edge edge)
{
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

template<typename TF>
void Grid<TF>::boundary_cyclic_2d(TF* restrict data)
{
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

// 
// template<typename TF>
// void Grid<TF>::transpose_zx(double* restrict ar, double* restrict as)
// {
//     const int jj = imax;
//     const int kk = imax*jmax;
// 
//     for (int k=0; k<kmax; k++)
//         for (int j=0; j<jmax; j++)
// #pragma ivdep
//             for (int i=0; i<imax; i++)
//             {
//                 const int ijk = i + j*jj + k*kk;
//                 ar[ijk] = as[ijk];
//             }
// }
// 
// template<typename TF>
// void Grid<TF>::transpose_xz(double* restrict ar, double* restrict as)
// {
//     const int jj = imax;
//     const int kk = imax*jmax;
// 
//     for (int k=0; k<kmax; k++)
//         for (int j=0; j<jmax; j++)
// #pragma ivdep
//             for (int i=0; i<imax; i++)
//             {
//                 const int ijk = i + j*jj + k*kk;
//                 ar[ijk] = as[ijk];
//             }
// }
// 
// template<typename TF>
// void Grid<TF>::transpose_xy(double* restrict ar, double* restrict as)
// {
//     const int jj = imax;
//     const int kk = imax*jmax;
// 
//     for (int k=0; k<kmax; k++)
//         for (int j=0; j<jmax; j++)
// #pragma ivdep
//             for (int i=0; i<imax; i++)
//             {
//                 const int ijk = i + j*jj + k*kk;
//                 ar[ijk] = as[ijk];
//             }
// }
// 
// template<typename TF>
// void Grid<TF>::transpose_yx(double* restrict ar, double* restrict as)
// {
//     const int jj = imax;
//     const int kk = imax*jmax;
// 
//     for (int k=0; k<kmax; k++)
//         for (int j=0; j<jmax; j++)
// #pragma ivdep
//             for (int i=0; i<imax; i++)
//             {
//                 const int ijk = i + j*jj + k*kk;
//                 ar[ijk] = as[ijk];
//             }
// }
// 
// template<typename TF>
// void Grid<TF>::transpose_yz(double* restrict ar, double* restrict as)
// {
//     const int jj = imax;
//     const int kk = imax*jmax;
// 
//     for (int k=0; k<kmax; k++)
//         for (int j=0; j<jmax; j++)
// #pragma ivdep
//             for (int i=0; i<imax; i++)
//             {
//                 const int ijk = i + j*jj + k*kk;
//                 ar[ijk] = as[ijk];
//             }
// }
// 
// template<typename TF>
// void Grid<TF>::transpose_zy(double* restrict ar, double* restrict as)
// {
//     const int jj = imax;
//     const int kk = imax*jmax;
// 
//     for (int k=0; k<kmax; k++)
//         for (int j=0; j<jmax; j++)
// #pragma ivdep
//             for (int i=0; i<imax; i++)
//             {
//                 const int ijk = i + j*jj + k*kk;
//                 ar[ijk] = as[ijk];
//             }
// }
// 
// template<typename TF>
// void Grid<TF>::get_max(double *var)
// {
// }
// 
// template<typename TF>
// void Grid<TF>::get_max(int *var)
// {
// }
// 
// template<typename TF>
// void Grid<TF>::get_sum(double *var)
// {
// }
// 
// template<typename TF>
// void Grid<TF>::get_prof(double *prof, int kcellsin)
// {
// }

// IO functions
template<typename TF>
void Grid<TF>::save_grid()
{
    // SAVE THE GRID
    FILE *pFile;
    char filename[256];
    std::sprintf(filename, "%s.%07d", "grid", 0);
    pFile = fopen(filename, "wbx");
    master.print_message("Saving \"%s\" ... ", filename);

    if (pFile == NULL)
    {
        master.print_message("FAILED\n");
        throw 1;
    }
    else
        master.print_message("OK\n");

    fwrite(&gd.x [gd.istart], sizeof(TF), gd.itot, pFile);
    fwrite(&gd.xh[gd.istart], sizeof(TF), gd.itot, pFile);
    fwrite(&gd.y [gd.jstart], sizeof(TF), gd.jtot, pFile);
    fwrite(&gd.yh[gd.jstart], sizeof(TF), gd.jtot, pFile);
    fwrite(&gd.z [gd.kstart], sizeof(TF), gd.ktot, pFile);
    fwrite(&gd.zh[gd.kstart], sizeof(TF), gd.ktot, pFile);
    fclose(pFile);
}

template<typename TF>
void Grid<TF>::load_grid()
{
    // LOAD THE GRID
    FILE *pFile;
    char filename[256];
    std::sprintf(filename, "%s.%07d", "grid", 0);
    pFile = fopen(filename, "rb");
    master.print_message("Loading \"%s\" ... ", filename);

    if (pFile == NULL)
    {
        master.print_message("FAILED\n");
        throw 1;
    }
    else
        master.print_message("OK\n");

    fread(&gd.x [gd.istart], sizeof(TF), gd.itot, pFile);
    fread(&gd.xh[gd.istart], sizeof(TF), gd.itot, pFile);
    fread(&gd.y [gd.jstart], sizeof(TF), gd.jtot, pFile);
    fread(&gd.yh[gd.jstart], sizeof(TF), gd.jtot, pFile);
    fread(&gd.z [gd.kstart], sizeof(TF), gd.ktot, pFile);
    fread(&gd.zh[gd.kstart], sizeof(TF), gd.ktot, pFile);
    fclose(pFile);

    // calculate the missing coordinates
    calculate();
}

template<>
void Grid<double>::load_fftw()
{
    // LOAD THE FFTW PLAN
    char filename[256];
    std::sprintf(filename, "%s.%07d", "fftwplan", 0);

    master.print_message("Loading \"%s\" ... ", filename);

    int n = fftw_import_wisdom_from_filename(filename);
    if (n == 0)
    {
        master.print_message("FAILED\n");
        throw 1;
    }
    else
        master.print_message("OK\n");

    // use the FFTW3 many interface in order to reduce function call overhead
    int rank = 1;
    int ni[] = {gd.itot};
    int nj[] = {gd.jtot};
    int istride = 1;
    int jstride = gd.iblock;
    int idist = gd.itot;
    int jdist = 1;
    fftw_r2r_kind kindf[] = {FFTW_R2HC};
    fftw_r2r_kind kindb[] = {FFTW_HC2R};
    iplanf = fftw_plan_many_r2r(rank, ni, gd.jmax, fftini, ni, istride, idist,
            fftouti, ni, istride, idist, kindf, FFTW_EXHAUSTIVE);
    iplanb = fftw_plan_many_r2r(rank, ni, gd.jmax, fftini, ni, istride, idist,
            fftouti, ni, istride, idist, kindb, FFTW_EXHAUSTIVE);
    jplanf = fftw_plan_many_r2r(rank, nj, gd.iblock, fftinj, nj, jstride, jdist,
            fftoutj, nj, jstride, jdist, kindf, FFTW_EXHAUSTIVE);
    jplanb = fftw_plan_many_r2r(rank, nj, gd.iblock, fftinj, nj, jstride, jdist,
            fftoutj, nj, jstride, jdist, kindb, FFTW_EXHAUSTIVE);

    fftwplan = true;

    fftw_forget_wisdom();
}


template<>
void Grid<float>::load_fftw()
{
    // LOAD THE FFTW PLAN
    char filename[256];
    std::sprintf(filename, "%s.%07d", "fftwplan", 0);

    master.print_message("Loading \"%s\" ... ", filename);

    int n = fftwf_import_wisdom_from_filename(filename);
    if (n == 0)
    {
        master.print_message("FAILED\n");
        throw 1;
    }
    else
        master.print_message("OK\n");

    // use the FFTW3 many interface in order to reduce function call overhead
    int rank = 1;
    int ni[] = {gd.itot};
    int nj[] = {gd.jtot};
    int istride = 1;
    int jstride = gd.iblock;
    int idist = gd.itot;
    int jdist = 1;
    fftwf_r2r_kind kindf[] = {FFTW_R2HC};
    fftwf_r2r_kind kindb[] = {FFTW_HC2R};
    iplanff = fftwf_plan_many_r2r(rank, ni, gd.jmax, fftini, ni, istride, idist,
            fftouti, ni, istride, idist, kindf, FFTW_EXHAUSTIVE);
    iplanbf = fftwf_plan_many_r2r(rank, ni, gd.jmax, fftini, ni, istride, idist,
            fftouti, ni, istride, idist, kindb, FFTW_EXHAUSTIVE);
    jplanff = fftwf_plan_many_r2r(rank, nj, gd.iblock, fftinj, nj, jstride, jdist,
            fftoutj, nj, jstride, jdist, kindf, FFTW_EXHAUSTIVE);
    jplanbf = fftwf_plan_many_r2r(rank, nj, gd.iblock, fftinj, nj, jstride, jdist,
            fftoutj, nj, jstride, jdist, kindb, FFTW_EXHAUSTIVE);

    fftwplan = true;

    fftwf_forget_wisdom();
}

template<typename TF>
void Grid<TF>::load()
{
    load_grid();
    load_fftw();
}

template<typename TF>
int Grid<TF>::save_field3d(TF* restrict data, TF* restrict tmp1, TF* restrict tmp2,
        char* filename, const TF offset)
{
    FILE *pFile;
    pFile = fopen(filename, "wbx");

    if (pFile == NULL)
        return 1;

    const int jj = gd.icells;
    const int kk = gd.icells*gd.jcells;

    // first, add the offset to the data
    for (int k=gd.kstart; k<gd.kend; ++k)
        for (int j=gd.jstart; j<gd.jend; ++j)
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                tmp1[ijk] = data[ijk] + offset;
            }

    // second, save the data to disk
    for (int k=gd.kstart; k<gd.kend; ++k)
        for (int j=gd.jstart; j<gd.jend; ++j)
        {
            const int ijk = gd.istart + j*jj + k*kk;
            if( fwrite(&tmp1[ijk], sizeof(TF), gd.imax, pFile) != gd.imax)
                return 1;
        }

    fclose(pFile);

    return 0;
}

template<typename TF>
int Grid<TF>::load_field3d(TF* restrict data, TF* restrict tmp1, TF* restrict tmp2,
        char* filename, const TF offset)
{
    FILE *pFile;
    pFile = fopen(filename, "rb");

    if (pFile == NULL)
        return 1;

    const int jj = gd.icells;
    const int kk = gd.icells*gd.jcells;

    // first, load the data from disk
    for (int k=gd.kstart; k<gd.kend; k++)
        for (int j=gd.jstart; j<gd.jend; j++)
        {
            const int ijk = gd.istart + j*jj + k*kk;
            if( fread(&tmp1[ijk], sizeof(TF), gd.imax, pFile) != gd.imax )
                return 1;
        }

    fclose(pFile);

    // second, remove the offset
    for (int k=gd.kstart; k<gd.kend; k++)
        for (int j=gd.jstart; j<gd.jend; j++)
            for (int i=gd.istart; i<gd.iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                data[ijk] = tmp1[ijk] - offset;
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
    int kk = gd.itot*gd.jmax;

    // Process the fourier transforms slice by slice.
    for (int k=0; k<gd.kblock; ++k)
    {
#pragma ivdep
        for (int n=0; n<gd.itot*gd.jmax; ++n)
        {
            const int ij = n;
            const int ijk = n + k*kk;
            fftini[ij] = data[ijk];
        }

        fftw_execute_wrapper<TF>(iplanf, iplanff);

#pragma ivdep
        for (int n=0; n<gd.itot*gd.jmax; ++n)
        {
            const int ij = n;
            const int ijk = n + k*kk;
            data[ijk] = fftouti[ij];
        }
    }

    kk = gd.iblock*gd.jtot;

    // do the second fourier transform
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
            // shift to use p in pressure solver
            data[ijk] = fftoutj[ij];
        }
    }
}

template<typename TF>
void Grid<TF>::fft_backward(TF* const restrict data,   TF* const restrict tmp1,
                            TF* const restrict fftini, TF* const restrict fftouti,
                            TF* const restrict fftinj, TF* const restrict fftoutj)
{
    int kk = gd.iblock*gd.jtot;

    // transform the second transform back
    for (int k=0; k<gd.kblock; ++k)
    {
        #pragma ivdep
        for (int n=0; n<gd.iblock*gd.jtot; ++n)
        {
            const int ij = n;
            const int ijk = n + k*kk;
            fftinj[ij] = data[ijk];
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

    kk = gd.itot*gd.jmax;

    // transform the first transform back
    for (int k=0; k<gd.kblock; k++)
    {
#pragma ivdep
        for (int n=0; n<gd.itot*gd.jmax; ++n)
        {
            const int ij = n;
            const int ijk = n + k*kk;
            fftini[ij] = data[ijk];
        }

        fftw_execute_wrapper<TF>(iplanb, iplanbf);

#pragma ivdep
        for (int n=0; n<gd.itot*gd.jmax; ++n)
        {
            const int ij = n;
            const int ijk = n + k*kk;
            // swap array here to avoid unnecessary 3d loop
            tmp1[ijk] = fftouti[ij] / gd.itot;
        }
    }
}
// 
// template<typename TF>
// int Grid<TF>::save_xz_slice(double* restrict data, double* restrict tmp, char* filename, int jslice)
// {
//     // extract the data from the 3d field without the ghost cells
//     const int jj  = icells;
//     const int kk  = icells*jcells;
//     const int kkb = imax;
// 
//     const int count = imax*kmax;
// 
//     for (int k=0; k<kmax; k++)
// #pragma ivdep
//         for (int i=0; i<imax; i++)
//         {
//             // take the modulus of jslice and jmax to have the right offset within proc
//             const int ijk  = i+igc + (jslice+jgc)*jj + (k+kgc)*kk;
//             const int ijkb = i + k*kkb;
//             tmp[ijkb] = data[ijk];
//         }
// 
//     FILE *pFile;
//     pFile = fopen(filename, "wbx");
//     if (pFile == NULL)
//         return 1;
// 
//     fwrite(tmp, sizeof(double), count, pFile);
//     fclose(pFile);
// 
//     return 0;
// }
// 
// template<typename TF>
// int Grid<TF>::save_yz_slice(double* restrict data, double* restrict tmp, char* filename, int islice)
// {
//     // Extract the data from the 3d field without the ghost cells
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
//     FILE *pFile;
//     pFile = fopen(filename, "wbx");
//     if (pFile == NULL)
//         return 1;
// 
//     fwrite(tmp, sizeof(double), count, pFile);
//     fclose(pFile);
// 
//     return 0;
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
//     const int count = imax*jmax;
// 
//     // Subtract the ghost cells in case of a pure 2d plane that does not have ghost cells.
//     if (kslice == -1)
//         kslice = -kgc;
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
//     FILE *pFile;
//     pFile = fopen(filename, "wbx");
//     if (pFile == NULL)
//         return 1;
// 
//     fwrite(tmp, sizeof(double), count, pFile);
//     fclose(pFile);
// 
//     return 0;
// }
// 
// template<typename TF>
// int Grid<TF>::load_xy_slice(double* restrict data, double* restrict tmp, char* filename, int kslice)
// {
//     const int count = imax*jmax;
// 
//     FILE *pFile;
//     pFile = fopen(filename, "rb");
//     if (pFile == NULL)
//         return 1;
// 
//     fread(tmp, sizeof(double), count, pFile);
//     fclose(pFile);
// 
//     // Subtract the ghost cells in case of a pure 2d plane that does not have ghost cells.
//     if (kslice == -1)
//         kslice = -kgc;
// 
//     // put the data back into a field with ghost cells
//     const int jj  = icells;
//     const int kk  = icells*jcells;
//     const int jjb = imax;
// 
//     for (int j=0; j<jmax; j++)
// #pragma ivdep
//         for (int i=0; i<imax; i++)
//         {
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
