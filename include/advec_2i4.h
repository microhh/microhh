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

#ifndef ADVEC_2I4
#define ADVEC_2I4

#include "advec.h"

/**
 * Derived class for 2nd order advection scheme with 4th order interpolation.
 */
class Advec2i4 : public Advec
{
  public:
    Advec2i4(Model *, Input *); ///< Constructor of the advection class.
    ~Advec2i4();                ///< Destructor of the advection class.

    void exec(); ///< Execute the advection scheme.

    unsigned long getTimeLimit(long unsigned int, double); ///< Get the limit on the time step imposed by the advection scheme.
    double get_cfl(double);                                ///< Get the CFL number.

  private:
    double calc_cfl(double *, double *, double *, double *, double); ///< Calculate the CFL number.

    void advecu(double *, double *, double *, double *, double *, double *, double *);           ///< Calculate longitudinal velocity advection.
    void advecv(double *, double *, double *, double *, double *, double *, double *);           ///< Calculate latitudinal velocity advection.
    void advecw(double *, double *, double *, double *, double *, double *, double *);           ///< Calculate vertical velocity advection.
    void advecs(double *, double *, double *, double *, double *, double *, double *, double *); ///< Calculate scalar advection.
};
#endif
