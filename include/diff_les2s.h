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

#ifndef DIFF_LES2S
#define DIFF_LES2S

#include "diff.h"

class Diff_les2s : public Diff
{
  public:
    Diff_les2s(Model *, Input *);
    ~Diff_les2s();

    int exec();
    int execvisc();

    unsigned long gettimelim(unsigned long, double);
    double getdn(double);

    double tPr;

    // GPU functions and variables
    int prepareDevice();
    int clearDevice();
    double *mlen_g; 

  private:
    int strain2(double *,
                double *, double *, double *,
                double *, double *,
                double *, double *,
                double *, double *, double *);

    int evisc(double *,
              double *, double *, double *, double *,
              double *, double *, double *,
              double *, double *,
              double *, double *, double *,
              double);
    int evisc_neutral(double *,
                      double *, double *, double *,
                      double *, double *,
                      double *, double *,
                      double);
    int diffu(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
    int diffv(double *, double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
    int diffw(double *, double *, double *, double *, double *, double *, double *, double *, double *);
    int diffc(double *, double *, double *, double *, double *, double *, double *, double *, double *, double);

    double calcdnmul(double *, double *, double);

    inline double phim(double);
    inline double phih(double);

    double cs;
};
#endif
