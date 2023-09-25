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

#ifndef PRES_4_H
#define PRES_4_H

#include "pres.h"
#include "defines.h"
#include "boundary_cyclic.h"

#ifdef USECUDA
#include <cufft.h>
#endif

template<typename TF>
class Pres_4 : public Pres<TF>
{
    public:
        Pres_4(Master&, Grid<TF>&, Fields<TF>&, FFT<TF>&, Input&);
        ~Pres_4();

        void init();
        void set_values();
        void create(Stats<TF>&);

        void exec(const double, Stats<TF>&);
        TF check_divergence();

        #ifdef USECUDA
        void prepare_device();
        void clear_device();
        #endif

    private:
        using Pres<TF>::master;
        using Pres<TF>::grid;
        using Pres<TF>::fields;
        using Pres<TF>::field3d_operators;
        using Pres<TF>::fft;
        Boundary_cyclic<TF> boundary_cyclic;

        std::vector<TF> bmati;
        std::vector<TF> bmatj;
        std::vector<TF> m1;
        std::vector<TF> m2;
        std::vector<TF> m3;
        std::vector<TF> m4;
        std::vector<TF> m5;
        std::vector<TF> m6;
        std::vector<TF> m7;

        #ifdef USECUDA
        using Pres<TF>::make_cufft_plan;
        using Pres<TF>::fft_forward;
        using Pres<TF>::fft_backward;

        TF* bmati_g;
        TF* bmatj_g;
        TF* m1_g;
        TF* m2_g;
        TF* m3_g;
        TF* m4_g;
        TF* m5_g;
        TF* m6_g;
        TF* m7_g;

        cufftDoubleComplex* ffti_complex_g;
        cufftDoubleComplex* fftj_complex_g;
        cufftHandle iplanf;
        cufftHandle jplanf;
        cufftHandle iplanb;
        cufftHandle jplanb;
        #endif

        template<bool>
        void input(TF* restrict,
                   const TF* restrict, const TF* restrict, const TF* restrict,
                   TF* restrict, TF* restrict, TF* restrict,
                   const TF* restrict, const TF);

        void solve(TF* restrict, TF* restrict, const TF* restrict,
                   const TF* restrict, const TF* restrict, const TF* restrict, const TF* restrict,
                   const TF* restrict, const TF* restrict, const TF* restrict,
                   TF* restrict, TF* restrict, TF* restrict, TF* restrict,
                   TF* restrict, TF* restrict, TF* restrict, TF* restrict,
                   TF* restrict, TF* restrict,
                   const int);

        template<bool>
        void output(TF* restrict, TF* restrict, TF* restrict,
                    const TF* restrict, const TF* restrict);

        void hdma(TF* restrict, TF* restrict, TF* restrict, TF* restrict,
                  TF* restrict, TF* restrict, TF* restrict, TF* restrict,
                  int);

        TF calc_divergence(const TF* restrict, const TF* restrict, const TF* restrict, const TF* restrict);

        const std::string tend_name = "pres";
        const std::string tend_longname = "Pressure";
};
#endif
