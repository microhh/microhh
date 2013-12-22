/*
 * MicroHH
 * Copyright (c) 2011-2013 Chiel van Heerwaarden
 * Copyright (c) 2011-2013 Thijs Heus
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

#ifndef PRES_G2
#define PRES_G2

#include <fftw3.h>
#include "grid.h"
#include "fields.h"
#include "pres.h"
#include "mpiinterface.h"

// forward declaration
class cmodel;

class cpres_g2 : public cpres
{
  public:
    cpres_g2(cmodel *);
    ~cpres_g2();

    int init();
    int setvalues();
    int exec(double);
    double check();

  private:
    bool allocated;

    double *bmati, *bmatj;
    double *a, *c;
    double *work2d;

    int pres_in(double *, 
                double *, double *, double *,
                double *, double *, double *,
                double *, double);
    int pres_solve(double *, double *, double *, double *,
                   double *, double *, 
                   double *, double *);
    int pres_out(double *, double *, double *,
                 double *, double *);
    double calcdivergence(double *, double *, double *, double *);

    // functions
    int tdma(double *, double *, double *, double *, 
             double *, double *);
};
#endif
