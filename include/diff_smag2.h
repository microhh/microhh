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

#ifndef DIFF_SMAG_2
#define DIFF_SMAG_2

#include "diff.h"

class DiffSmag2 : public Diff
{
  public:
    DiffSmag2(Model *, Input *);
    ~DiffSmag2();

    void exec();
    void execViscosity();

    unsigned long getTimeLimit(unsigned long, double);
    double getdn(double);

    double tPr;

    // GPU functions and variables
    void prepareDevice();
    void clearDevice();

    double *mlen_g; 

  private:
    void strain2(double *,
                 double *, double *, double *,
                 double *, double *,
                 double *, double *,
                 double *, double *, double *);

    void evisc(double *,
               double *, double *, double *, double *,
               double *, double *, double *,
               double *, double *,
               double *, double *, double *,
               double);
    void eviscNeutral(double *,
                      double *, double *, double *,
                      double *, double *,
                      double *, double *,
                      double);
    void diffu(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
    void diffv(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
    void diffw(double *, double *, double *, double *, double *, double *, double *, double *, double *);
    void diffc(double *, double *, double *, double *, double *, double *, double *, double *, double *, double);

    double calc_dnmul(double *, double *, double);

    inline double phim(double);
    inline double phih(double);

    double cs;
};
#endif
