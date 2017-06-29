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

#ifndef PRES_4
#define PRES_4

#include "pres.h"
#include "defines.h"

#ifdef USECUDA
#include <cufft.h>
#endif

class Model;

class Pres_4 : public Pres
{
    public:
        Pres_4(Model*, Input*);
        ~Pres_4();

        void init();
        void set_values();

        void exec(double);
        double check_divergence();

#ifdef USECUDA
        void prepare_device();
        void clear_device();
#endif

    private:
        double* bmati;
        double* bmatj;
        double* m1;
        double* m2;
        double* m3;
        double* m4;
        double* m5;
        double* m6;
        double* m7;

#ifdef USECUDA
        double* bmati_g;
        double* bmatj_g;
        double* m1_g;
        double* m2_g;
        double* m3_g;
        double* m4_g;
        double* m5_g;
        double* m6_g;
        double* m7_g;

        cufftDoubleComplex* ffti_complex_g;
        cufftDoubleComplex* fftj_complex_g;
        cufftHandle iplanf;
        cufftHandle jplanf; 
        cufftHandle iplanb;
        cufftHandle jplanb; 
#endif

        template<bool>
        void input(double* restrict, 
                   double* restrict, double* restrict, double* restrict,
                   double* restrict, double* restrict, double* restrict,
                   double* restrict, double);

        void solve(double* restrict, double* restrict, double* restrict,
                   double* restrict, double* restrict, double* restrict, double* restrict,
                   double* restrict, double* restrict, double* restrict,
                   double* restrict, double* restrict, double* restrict, double* restrict,
                   double* restrict, double* restrict, double* restrict, double* restrict,
                   double* restrict, double* restrict,
                   int);

        template<bool>
        void output(double* restrict, double* restrict, double* restrict,
                    double* restrict, double* restrict);

        void hdma(double* restrict, double* restrict, double* restrict, double* restrict,
                  double* restrict, double* restrict, double* restrict, double* restrict,
                  int);

        double calc_divergence(double* restrict, double* restrict, double* restrict, double* restrict);
};
#endif
