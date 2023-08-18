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

#ifdef USEMPI

#include <cstdio>
#include <stdexcept>

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
    auto& md = master.get_MPI_data();

    // file saving and loading, take C-ordering into account
    int totsizei  = gd.itot;
    int subsizei  = gd.imax;
    int substarti = md.mpicoordx*gd.imax;
    MPI_Type_create_subarray(1, &totsizei, &subsizei, &substarti, MPI_ORDER_C, mpi_fp_type<TF>(), &subi);
    MPI_Type_commit(&subi);

    int totsizej  = gd.jtot;
    int subsizej  = gd.jmax;
    int substartj = md.mpicoordy*gd.jmax;
    MPI_Type_create_subarray(1, &totsizej, &subsizej, &substartj, MPI_ORDER_C, mpi_fp_type<TF>(), &subj);
    MPI_Type_commit(&subj);

    mpitypes = true;
}

template<typename TF>
void Grid<TF>::exit_mpi()
{
    if (mpitypes)
    {
        MPI_Type_free(&subi);
        MPI_Type_free(&subj);
    }
}

template<typename TF>
void Grid<TF>::save_grid()
{
    auto& md = master.get_MPI_data();
    int nerror = 0;
    char filename[256];
    std::sprintf(filename, "%s.%07d", "grid", 0);
    master.print_message("Saving \"%s\" ... ", filename);

    MPI_File fh;
    if (MPI_File_open(md.commxy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh))
    {
        master.print_message("FAILED\n");
        nerror++;
    }

    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("Error in grid");

    // select noncontiguous part of 3d array to store the selected data
    MPI_Offset fileoff = 0; // the offset within the file (header size)
    char name[] = "native";

    MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subi, name, MPI_INFO_NULL);
    if (md.mpicoordy == 0)
        MPI_File_write(fh, &gd.x[gd.istart], gd.imax, mpi_fp_type<TF>(), MPI_STATUS_IGNORE);
    MPI_Barrier(md.commxy);
    fileoff += gd.itot*sizeof(TF);

    MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subi, name, MPI_INFO_NULL);
    if (md.mpicoordy == 0)
        MPI_File_write(fh, &gd.xh[gd.istart], gd.imax, mpi_fp_type<TF>(), MPI_STATUS_IGNORE);
    MPI_Barrier(md.commxy);
    fileoff += gd.itot*sizeof(TF);

    MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subj, name, MPI_INFO_NULL);
    if (md.mpicoordx == 0)
        MPI_File_write(fh, &gd.y[gd.jstart], gd.jmax, mpi_fp_type<TF>(), MPI_STATUS_IGNORE);
    MPI_Barrier(md.commxy);
    fileoff += gd.jtot*sizeof(TF);

    MPI_File_set_view(fh, fileoff, mpi_fp_type<TF>(), subj, name, MPI_INFO_NULL);
    if (md.mpicoordx == 0)
        MPI_File_write(fh, &gd.yh[gd.jstart], gd.jmax, mpi_fp_type<TF>(), MPI_STATUS_IGNORE);
    MPI_Barrier(md.commxy);

    MPI_File_sync(fh);
    if (MPI_File_close(&fh))
        nerror++;

    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("Grid save FAILED");

    if (master.get_mpiid() == 0)
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
    // auto& md = master.get_MPI_data();

    int nerror = 0;

    // LOAD THE GRID
    char filename[256];
    std::sprintf(filename, "%s.%07d", "grid", 0);
    if (master.get_mpiid() == 0)
        std::printf("Loading \"%s\" ... ", filename);

    FILE *pFile;
    if (master.get_mpiid() == 0)
    {
        pFile = fopen(filename, "rb");
        if (pFile == NULL)
            ++nerror;
        else
        {
            int n = (2*gd.itot+2*gd.jtot)*sizeof(TF);
            if (fseek(pFile, n, SEEK_SET) != 0);
                ++nerror;
            if(fread(&gd.z [gd.kstart], sizeof(TF), gd.kmax, pFile) != (unsigned)gd.kmax )
                ++nerror;
            if(fread(&gd.zh[gd.kstart], sizeof(TF), gd.kmax, pFile) != (unsigned)gd.kmax )
                ++nerror;
            fclose(pFile);
        }
    }

    // Communicate the file read error over all procs.

    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("Error in grid");
    else
        master.print_message("OK\n");

    master.broadcast(&gd.z [gd.kstart], gd.kmax);
    master.broadcast(&gd.zh[gd.kstart], gd.kmax);

    // Calculate the missing coordinates.
    calculate();
}

template class Grid<double>;
template class Grid<float>;
#endif
