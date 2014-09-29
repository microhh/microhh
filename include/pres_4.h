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

// forward declaration
class cmodel;

class cpres_4 : public cpres
{
  public:
    cpres_4(cmodel *, cinput *);
    ~cpres_4();

    void init();
    void setvalues();

    void exec(double);
    double check();

  private:
    int jslice;

    double *bmati, *bmatj;
    double *m1,*m2,*m3,*m4,*m5,*m6,*m7;

    void pres_in(double *, 
                 double *, double *, double *,
                 double *, double *, double *,
                 double *, double);
    void pres_solve(double *, double *, double *,
                    double *, double *, double *, double *,
                    double *, double *, double *,
                    double *, double *, double *, double *,
                    double *, double *, double *, double *,
                    double *, double *);
    void pres_out(double *, double *, double *,
                  double *, double *);
    double calcdivergence(double *, double *, double *, double *);

    // functions
    void hdma(double *, double *, double *, double *,
              double *, double *, double *, double *);
};
#endif
