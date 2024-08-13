/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

#include <stdexcept>

#include "master.h"
#include "grid.h"
#include "fft.h"


template<typename TF>
FFT<TF>::FFT(Master& masterin, Grid<TF>& gridin) :
    master(masterin), grid(gridin),
    transpose(master, grid)
{
    has_fftw_plan = false;

    // Initialize the pointers to zero.
    fftini  = nullptr;
    fftouti = nullptr;
    fftinj  = nullptr;
    fftoutj = nullptr;
}


#ifdef FLOAT_SINGLE
template<>
void FFT<float>::init()
{
    auto& gd = grid.get_grid_data();

    fftini  = fftwf_alloc_real(gd.itot*gd.jmax);
    fftouti = fftwf_alloc_real(gd.itot*gd.jmax);
    fftinj  = fftwf_alloc_real(gd.jtot*gd.iblock);
    fftoutj = fftwf_alloc_real(gd.jtot*gd.iblock);

    transpose.init();
}
#else
template<>
void FFT<double>::init()
{
    auto& gd = grid.get_grid_data();

    fftini  = fftw_alloc_real(gd.itot*gd.jmax);
    fftouti = fftw_alloc_real(gd.itot*gd.jmax);
    fftinj  = fftw_alloc_real(gd.jtot*gd.iblock);
    fftoutj = fftw_alloc_real(gd.jtot*gd.iblock);

    transpose.init();
}
#endif


#ifdef FLOAT_SINGLE
template<>
FFT<float>::~FFT()
{
    if (has_fftw_plan)
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
#else
template<>
FFT<double>::~FFT()
{
    if (has_fftw_plan)
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
#endif


#ifdef FLOAT_SINGLE
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
        throw std::runtime_error("Error loading FFTW Plan");
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
            fftouti, ni, istride, idist, kindf, FFTW_ESTIMATE);
    iplanbf = fftwf_plan_many_r2r(rank, ni, gd.jmax, fftini, ni, istride, idist,
            fftouti, ni, istride, idist, kindb, FFTW_ESTIMATE);
    jplanff = fftwf_plan_many_r2r(rank, nj, gd.iblock, fftinj, nj, jstride, jdist,
            fftoutj, nj, jstride, jdist, kindf, FFTW_ESTIMATE);
    jplanbf = fftwf_plan_many_r2r(rank, nj, gd.iblock, fftinj, nj, jstride, jdist,
            fftoutj, nj, jstride, jdist, kindb, FFTW_ESTIMATE);

    has_fftw_plan = true;

    fftwf_forget_wisdom();
}
#else
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
        throw std::runtime_error("Error loading FFTW Plan");
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
            fftouti, ni, istride, idist, kindf, FFTW_ESTIMATE);
    iplanb = fftw_plan_many_r2r(rank, ni, gd.jmax, fftini, ni, istride, idist,
            fftouti, ni, istride, idist, kindb, FFTW_ESTIMATE);
    jplanf = fftw_plan_many_r2r(rank, nj, gd.iblock, fftinj, nj, jstride, jdist,
            fftoutj, nj, jstride, jdist, kindf, FFTW_ESTIMATE);
    jplanb = fftw_plan_many_r2r(rank, nj, gd.iblock, fftinj, nj, jstride, jdist,
            fftoutj, nj, jstride, jdist, kindb, FFTW_ESTIMATE);

    has_fftw_plan = true;

    fftw_forget_wisdom();
}
#endif


#ifdef FLOAT_SINGLE
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
                                  fftouti, ni, istride, idist, kindf, FFTW_ESTIMATE);
    iplanbf = fftwf_plan_many_r2r(rank, ni, gd.jmax, fftini, ni, istride, idist,
                                  fftouti, ni, istride, idist, kindb, FFTW_ESTIMATE);
    jplanff = fftwf_plan_many_r2r(rank, nj, gd.iblock, fftinj, nj, jstride, jdist,
                                  fftoutj, nj, jstride, jdist, kindf, FFTW_ESTIMATE);
    jplanbf = fftwf_plan_many_r2r(rank, nj, gd.iblock, fftinj, nj, jstride, jdist,
                                  fftoutj, nj, jstride, jdist, kindb, FFTW_ESTIMATE);

    has_fftw_plan = true;

    int nerror = 0;
    if (master.get_mpiid() == 0)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", "fftwplan", 0);

        master.print_message("Saving \"%s\" ... ", filename);

        int n = fftwf_export_wisdom_to_filename(filename);
        if (n == 0)
        {
            master.print_message("FAILED\n");
            nerror++;
        }
        else
            master.print_message("OK\n");
    }

    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("Error saving FFTW plan");
}
#else
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
                                fftouti, ni, istride, idist, kindf, FFTW_ESTIMATE);
    iplanb = fftw_plan_many_r2r(rank, ni, gd.jmax, fftini, ni, istride, idist,
                                fftouti, ni, istride, idist, kindb, FFTW_ESTIMATE);
    jplanf = fftw_plan_many_r2r(rank, nj, gd.iblock, fftinj, nj, jstride, jdist,
                                fftoutj, nj, jstride, jdist, kindf, FFTW_ESTIMATE);
    jplanb = fftw_plan_many_r2r(rank, nj, gd.iblock, fftinj, nj, jstride, jdist,
                                fftoutj, nj, jstride, jdist, kindb, FFTW_ESTIMATE);

    has_fftw_plan = true;

    int nerror = 0;
    if (master.get_mpiid() == 0)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", "fftwplan", 0);

        master.print_message("Saving \"%s\" ... ", filename);

        int n = fftw_export_wisdom_to_filename(filename);
        if (n == 0)
        {
            master.print_message("FAILED\n");
            nerror++;
        }
        else
            master.print_message("OK\n");
    }

    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("Error saving FFTW plan");
}
#endif


namespace
{
    template<typename> void fftw_execute_wrapper(const fftw_plan&, const fftwf_plan&);

    #ifdef FLOAT_SINGLE
    template<>
    void fftw_execute_wrapper<float>(const fftw_plan& p, const fftwf_plan& pf)
    {
        fftwf_execute(pf);
    }
    #else
    template<>
    void fftw_execute_wrapper<double>(const fftw_plan& p, const fftwf_plan& pf)
    {
        fftw_execute(p);
    }
    #endif

    #ifndef USEMPI
    template<typename TF>
    void fft_forward(TF* const restrict data,   TF* const restrict tmp1,
                     TF* const restrict fftini, TF* const restrict fftouti,
                     TF* const restrict fftinj, TF* const restrict fftoutj,
                     fftw_plan& iplanf, fftwf_plan& iplanff,
                     fftw_plan& jplanf, fftwf_plan& jplanff,
                     const Grid_data<TF>& gd, Transpose<TF>& transpose)
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
    void fft_backward(TF* const restrict data,   TF* const restrict tmp1,
                      TF* const restrict fftini, TF* const restrict fftouti,
                      TF* const restrict fftinj, TF* const restrict fftoutj,
                      fftw_plan& iplanb, fftwf_plan& iplanbf,
                      fftw_plan& jplanb, fftwf_plan& jplanbf,
                      const Grid_data<TF>& gd, Transpose<TF>& transpose)
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

    #else
    template<typename TF>
    void fft_forward(TF* const restrict data,   TF* const restrict tmp1,
                     TF* const restrict fftini, TF* const restrict fftouti,
                     TF* const restrict fftinj, TF* const restrict fftoutj,
                     fftw_plan& iplanf, fftwf_plan& iplanff,
                     fftw_plan& jplanf, fftwf_plan& jplanff,
                     const Grid_data<TF>& gd, Transpose<TF>& transpose)
    {
        // Transpose the pressure field.
        transpose.exec_zx(tmp1, data);

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
        transpose.exec_xy(data, tmp1);

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
        transpose.exec_yz(data, tmp1);
    }

    template<typename TF>
    void fft_backward(TF* const restrict data,   TF* const restrict tmp1,
                      TF* const restrict fftini, TF* const restrict fftouti,
                      TF* const restrict fftinj, TF* const restrict fftoutj,
                      fftw_plan& iplanb, fftwf_plan& iplanbf,
                      fftw_plan& jplanb, fftwf_plan& jplanbf,
                      const Grid_data<TF>& gd, Transpose<TF>& transpose)
    {
        // Transpose back to y.
        transpose.exec_zy(tmp1, data);

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
        transpose.exec_yx(tmp1, data);

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
        transpose.exec_xz(tmp1, data);
    }
    #endif
}

template<typename TF>
void FFT<TF>::exec_forward(TF* const restrict data, TF* const restrict tmp1)
{
    fft_forward(data, tmp1, fftini, fftouti, fftinj, fftoutj,
            iplanf, iplanff, jplanf, jplanff, grid.get_grid_data(), transpose);
}

template<typename TF>
void FFT<TF>::exec_backward(TF* const restrict data, TF* const restrict tmp1)
{
    fft_backward(data, tmp1, fftini, fftouti, fftinj, fftoutj,
            iplanb, iplanbf, jplanb, jplanbf, grid.get_grid_data(), transpose);
}


#ifdef FLOAT_SINGLE
template class FFT<float>;
#else
template class FFT<double>;
#endif
