/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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
    void setValues();

    void exec(double);
    double checkDivergence();

    #ifdef USECUDA
    void prepareDevice();
    void clearDevice();
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

    void input(double *, 
               double *, double *, double *,
               double *, double *, double *,
               double *, double);

    void solve(double *, double *, double *,
               double *, double *, double *, double *,
               double *, double *, double *,
               double *, double *, double *, double *,
               double *, double *, double *, double *,
               double *, double *,
               int);

    void output(double *, double *, double *,
                double *, double *);

    void hdma(double *, double *, double *, double *,
              double *, double *, double *, double *,
              int);

    double calcDivergence(double *, double *, double *, double *);
};
#endif
