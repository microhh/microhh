/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
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

#include "master.h"
#include "grid.h"
#include "fft.h"

template<typename TF>
FFT<TF>::FFT(Master& masterin, Grid<TF>& gridin) :
    master(masterin), grid(gridin),
    transpose(master, grid)
{
    fftw_plan = false;

    // Initialize the pointers to zero.
    fftini  = nullptr;
    fftouti = nullptr;
    fftinj  = nullptr;
    fftoutj = nullptr;
}

template<typename TF>
void FFT<TF>::init()
{
    allocate_fftw();
}

template<> void FFT<double>::allocate_fftw()
{
    auto& gd = grid.get_grid_data();

    fftini  = fftw_alloc_real(gd.itot*gd.jmax);
    fftouti = fftw_alloc_real(gd.itot*gd.jmax);
    fftinj  = fftw_alloc_real(gd.jtot*gd.iblock);
    fftoutj = fftw_alloc_real(gd.jtot*gd.iblock);
}

template<> void FFT<float>::allocate_fftw()
{
    auto& gd = grid.get_grid_data();

    fftini  = fftwf_alloc_real(gd.itot*gd.jmax);
    fftouti = fftwf_alloc_real(gd.itot*gd.jmax);
    fftinj  = fftwf_alloc_real(gd.jtot*gd.iblock);
    fftoutj = fftwf_alloc_real(gd.jtot*gd.iblock);
}

template<>
void FFT<double>::fftw_finish()
{
    if (fftw_plan)
    {
        fftw_destroy_plan(iplanf);
        fftw_destroy_plan(iplanb);
        fftw_destroy_plan(jplanf);
        fftw_destroy_plan(jplanb);
    }

    fftw_free(fftini);
    fftw_free(fftouti);
    fftw_free(fftinj);
    fftw_free(fftoutj);

    fftw_cleanup();
}

template<>
void FFT<float>::fftw_finish()
{
    if (fftw_plan)
    {
        fftwf_destroy_plan(iplanff);
        fftwf_destroy_plan(iplanbf);
        fftwf_destroy_plan(jplanff);
        fftwf_destroy_plan(jplanbf);
    }

    fftwf_free(fftini);
    fftwf_free(fftouti);
    fftwf_free(fftinj);
    fftwf_free(fftoutj);

    fftwf_cleanup();
}

template<typename TF>
FFT<TF>::~FFT()
{
    fftw_finish();
}

template<>
void FFT<double>::load()
{
    // LOAD THE FFTW PLAN
    auto& gd = grid.get_grid_data();

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

    fftw_plan = true;

    fftw_forget_wisdom();
}


template<>
void FFT<float>::load()
{
    // LOAD THE FFTW PLAN
    auto& gd = grid.get_grid_data();

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

    fftw_plan = true;

    fftwf_forget_wisdom();
}


template<>
void FFT<double>::save()
{
    // SAVE THE FFTW PLAN IN ORDER TO ENSURE BITWISE IDENTICAL RESTARTS
    // Use the FFTW3 many interface in order to reduce function call overhead.
    auto& gd = grid.get_grid_data();

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

    fftw_plan = true;

    if (master.mpiid == 0)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", "fftwplan", 0);

        master.print_message("Saving \"%s\" ... ", filename);

        int n = fftw_export_wisdom_to_filename(filename);
        if (n == 0)
        {
            master.print_message("FAILED\n");
            throw 1;
        }
        else
            master.print_message("OK\n");
    }
}

template<>
void FFT<float>::save()
{
    // SAVE THE FFTW PLAN IN ORDER TO ENSURE BITWISE IDENTICAL RESTARTS
    // Use the FFTW3 many interface in order to reduce function call overhead.
    auto& gd = grid.get_grid_data();

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

    fftw_plan = true;

    if (master.mpiid == 0)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", "fftwplan", 0);

        master.print_message("Saving \"%s\" ... ", filename);

        int n = fftwf_export_wisdom_to_filename(filename);
        if (n == 0)
        {
            master.print_message("FAILED\n");
            throw 1;
        }
        else
            master.print_message("OK\n");
    }
}

template class FFT<double>;
template class FFT<float>;
