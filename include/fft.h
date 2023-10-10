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

#ifndef FFT_H
#define FFT_H

#include <fftw3.h>
#include "transpose.h"

class Master;
template<typename> class Grid;

template<typename TF>
class FFT
{
    public:
        FFT(Master&, Grid<TF>&);
        ~FFT();

        void exec_forward (TF* const restrict, TF* const restrict);
        void exec_backward(TF* const restrict, TF* const restrict);

        void init();
        void load();
        void save();

    private:
        Master& master; // Reference to master class.
        Grid<TF>& grid; // Reference to grid class.
        Transpose<TF> transpose; // Reference to grid class.

        TF *fftini, *fftouti; // Help arrays for fast-fourier transforms in x-direction.
        TF *fftinj, *fftoutj; // Help arrays for fast-fourier transforms in y-direction.
        fftw_plan iplanf, iplanb; // FFTW3 plans for forward and backward transforms in x-direction.
        fftw_plan jplanf, jplanb; // FFTW3 plans for forward and backward transforms in y-direction.
        fftwf_plan iplanff, iplanbf; // FFTW3 plans for forward and backward transforms in x-direction.
        fftwf_plan jplanff, jplanbf; // FFTW3 plans for forward and backward transforms in y-direction.

        bool has_fftw_plan;
};
#endif
