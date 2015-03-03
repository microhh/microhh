/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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

class Pres4 : public Pres
{
  public:
    Pres4(Model *, Input *);
    ~Pres4();

    void init();
    void set_values();

    void exec(double);
    double checkDivergence();

    #ifdef USECUDA
    void prepare_device();
    void clear_device();
    #endif

  private:
    double *bmati, *bmatj;
    double *m1,*m2,*m3,*m4,*m5,*m6,*m7;

    #ifdef USECUDA
    double *bmati_g, *bmatj_g;
    double *m1_g,*m2_g,*m3_g,*m4_g,*m5_g,*m6_g,*m7_g;

    cufftDoubleComplex *ffti_complex_g, *fftj_complex_g; 
    cufftHandle iplanf, jplanf; 
    cufftHandle iplanb, jplanb; 
    #endif

    template<bool>
    void input(double * restrict, 
               double * restrict, double * restrict, double * restrict,
               double * restrict, double * restrict, double * restrict,
               double * restrict, double);

    void solve(double * restrict, double * restrict, double * restrict,
               double * restrict, double * restrict, double * restrict, double * restrict,
               double * restrict, double * restrict, double * restrict,
               double * restrict, double * restrict, double * restrict, double * restrict,
               double * restrict, double * restrict, double * restrict, double * restrict,
               double * restrict, double * restrict,
               int);

    template<bool>
    void output(double * restrict, double * restrict, double * restrict,
                double * restrict, double * restrict);

    void hdma(double * restrict, double * restrict, double * restrict, double * restrict,
              double * restrict, double * restrict, double * restrict, double * restrict,
              int);

    double calcDivergence(double * restrict, double * restrict, double * restrict, double * restrict);
};
#endif
