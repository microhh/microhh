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

#ifndef ADVEC_4M
#define ADVEC_4M

#include "advec.h"

/**
 * Derived class for fully conservative 4th order advection scheme.
 * Fully mass, momentum and energy conserving advection scheme based on the paper
 * of Morinishi et al., (1998).
 */
class cadvec_4m : public cadvec
{
  public:
    cadvec_4m(cmodel *, cinput *); ///< Constructor of the advection class.
    ~cadvec_4m();                  ///< Destructor of the advection class.

    unsigned long gettimelim(long unsigned int, double); ///< Get the limit on the time step imposed by the advection scheme.
    double getcfl(double);                               ///< Get the CFL number.
    int exec();                                          ///< Execute the advection scheme.

  private:
    double calccfl(double *, double *, double *, double *, double);         ///< Calculate the CFL number.
    int advecu(double *, double *, double *, double *, double *);           ///< Calculate longitudinal velocity advection.
    int advecv(double *, double *, double *, double *, double *);           ///< Calculate latitudinal velocity advection.
    int advecw(double *, double *, double *, double *, double *);           ///< Calculate vertical velocity advection.
    int advecs(double *, double *, double *, double *, double *, double *); ///< Calculate scalar advection.

    inline double interp2(const double, const double);                                           ///< 2nd order interpolation.
    inline double interp4(const double, const double, const double, const double);               ///< 4th order interpolation.
    inline double grad4  (const double, const double, const double, const double, const double); ///< 4th order gradient.
    inline double grad4x (const double, const double, const double, const double);               ///< 4th order gradient (only numerator).

    inline double interp4biasbot(const double, const double, const double, const double); ///< 4th order interpolation (bottom boundary).
    inline double interp4biastop(const double, const double, const double, const double); ///< 4th order interpolation (top boundary).
    inline double grad4xbiasbot (const double, const double, const double, const double); ///< 4th order interpolation (bottom boundary, only numerator).
    inline double grad4xbiastop (const double, const double, const double, const double); ///< 4th order interpolation (top boundary, only numerator).
};
#endif
