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

template<typename TF>
int Grid<TF>::save_xz_slice(TF* restrict data, TF* restrict tmp, char* filename, int jslice)
{
    // extract the data from the 3d field without the ghost cells
    const int jj  = gd.icells;
    const int kk  = gd.icells*gd.jcells;
    const int kkb = gd.imax;

    const int count = gd.imax*gd.kmax;

    for (int k=0; k<gd.kmax; k++)
#pragma ivdep
        for (int i=0; i<gd.imax; i++)
        {
            // take the modulus of jslice and gd.jmax to have the right offset within proc
            const int ijk  = i+gd.igc + (jslice+gd.jgc)*jj + (k+gd.kgc)*kk;
            const int ijkb = i + k*kkb;
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
int Grid<TF>::save_yz_slice(TF* restrict data, TF* restrict tmp, char* filename, int islice)
{
    // Extract the data from the 3d field without the ghost cells
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

    FILE *pFile;
    pFile = fopen(filename, "wbx");
    if (pFile == NULL)
        return 1;

    fwrite(tmp, sizeof(TF), count, pFile);
    fclose(pFile);

    return 0;
}

template<typename TF>
int Grid<TF>::save_xy_slice(TF* restrict data, TF* restrict tmp, char* filename, int kslice)
{
    // extract the data from the 3d field without the ghost cells
    const int jj  = gd.icells;
    const int kk  = gd.icells*gd.jcells;
    const int jjb = gd.imax;

    const int count = gd.imax*gd.jmax;

    // Subtract the ghost cells in case of a pure 2d plane that does not have ghost cells.
    if (kslice == -1)
        kslice = -gd.kgc;

    for (int j=0; j<gd.jmax; j++)
#pragma ivdep
        for (int i=0; i<gd.imax; i++)
        {
            // take the modulus of jslice and jmax to have the right offset within proc
            const int ijk  = i+gd.igc + (j+gd.jgc)*jj + (kslice+gd.kgc)*kk;
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
int Grid<TF>::load_xy_slice(TF* restrict data, TF* restrict tmp, char* filename, int kslice)
{
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

template class Grid<double>;
template class Grid<float>;
#endif
