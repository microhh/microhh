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

#ifndef PRES_4
#define PRES_4

#include "pres.h"

// forward declaration
class cmodel;

class cpres_4 : public cpres
{
  public:
    cpres_4(cmodel *);
    ~cpres_4();

    int init();
    int setvalues();
    int exec(double);
    double check();

  private:
    bool allocated;

    double *bmati, *bmatj;
    double *m0,*m1,*m2,*m3,*m4,*m5,*m6,*m7,*m8;
    double *work2d;

    // CvH remove later
    double *m0temp,*m1temp,*m2temp,*m3temp,*m4temp,*m5temp,*m6temp,*m7temp,*m8temp,*ptemp;

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
    int hdma(double *, double *, double *, double *,
             double *, double *, double *, double *);
};
#endif
