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

#ifndef USEMPI

#include <cstdio>
#include <stdexcept>

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

    int nerror = 0;
    if (pFile == NULL)
    {
        master.print_message("FAILED\n");
        nerror++;
    }
    else
        master.print_message("OK\n");
    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("Error in grid");

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

    int nerror = 0;
    if (pFile == NULL)
        nerror++;

    if(fread(&gd.x [gd.istart], sizeof(TF), gd.itot, pFile) != (unsigned)gd.itot )
        ++nerror;
    if(fread(&gd.xh[gd.istart], sizeof(TF), gd.itot, pFile) != (unsigned)gd.itot )
        ++nerror;
    if(fread(&gd.y [gd.jstart], sizeof(TF), gd.jtot, pFile) != (unsigned)gd.jtot )
        ++nerror;
    if(fread(&gd.yh[gd.jstart], sizeof(TF), gd.jtot, pFile) != (unsigned)gd.jtot )
        ++nerror;
    if(fread(&gd.z [gd.kstart], sizeof(TF), gd.ktot, pFile) != (unsigned)gd.ktot )
        ++nerror;
    if(fread(&gd.zh[gd.kstart], sizeof(TF), gd.ktot, pFile) != (unsigned)gd.ktot )
        ++nerror;
    fclose(pFile);

    master.sum(&nerror, 1);

    if (nerror)
    {
        master.print_message("FAILED\n");
        throw std::runtime_error("Error in grid");
    }
    else
        master.print_message("OK\n");

    // calculate the missing coordinates
    calculate();
}

template class Grid<double>;
template class Grid<float>;
#endif
