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

#include <cstdio>
#include <iostream>
#include <cmath>
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

#ifndef USEMPI
template<typename TF>
int Field3d_io<TF>::save_field3d(TF* restrict data, TF* restrict tmp1, TF* restrict tmp2,
        char* filename, const TF offset)
{
    auto& gd = grid.get_grid_data();

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
int Field3d_io<TF>::load_field3d(TF* restrict data, TF* restrict tmp1, TF* restrict tmp2,
        char* filename, const TF offset)
{
    auto& gd = grid.get_grid_data();

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
int Field3d_io<TF>::save_xz_slice(TF* restrict data, TF* restrict tmp, char* filename, int jslice)
{
    auto& gd = grid.get_grid_data();

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
int Field3d_io<TF>::save_yz_slice(TF* restrict data, TF* restrict tmp, char* filename, int islice)
{
    auto& gd = grid.get_grid_data();

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
int Field3d_io<TF>::save_xy_slice(TF* restrict data, TF* restrict tmp, char* filename, int kslice)
{
    auto& gd = grid.get_grid_data();

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
int Field3d_io<TF>::load_xy_slice(TF* restrict data, TF* restrict tmp, char* filename, int kslice)
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
