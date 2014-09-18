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

#ifndef PRES_2
#define PRES_2

#include "pres.h"

#ifdef USECUDA
#include <cufft.h>
#endif

// forward declaration
class cmodel;

class cpres_2 : public cpres
{
  public:
    cpres_2(cmodel *, cinput *);
    ~cpres_2();

    void init();
    void setvalues();

    void exec(double);
    double check();

#ifdef USECUDA
    int prepareGPU();
#endif

  private:
    double *bmati, *bmatj;
    double *a, *c;
    double *work2d;

    // GPU
#ifdef USECUDA
    double *bmati_g, *bmatj_g;
    double *a_g, *c_g;
    double *work2d_g;

    cufftDoubleComplex *ffti_complex_g, *fftj_complex_g; 
    cufftHandle iplanf, jplanf; 
    cufftHandle iplanb, jplanb; 
#endif

    void pres_in(double *, 
                 double *, double *, double *,
                 double *, double *, double *,
                 double *, double *, double *,
                 double);
    void pres_solve(double *, double *, double *,
                    double *, double *,
                    double *, double *, double *, double *);
    void pres_out(double *, double *, double *,
                  double *, double *);
    double calcdivergence(double *, double *, double *, double *, double *, double *);

    void tdma(double *, double *, double *, double *, 
              double *, double *);
};
#endif
