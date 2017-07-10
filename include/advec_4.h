/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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

#ifndef ADVEC_4
#define ADVEC_4

#include "advec.h"
#include "defines.h"

class Model;
class Input;

/**
 * Derived class for 4th order advection scheme.
 */
class Advec_4 : public Advec
{
    public:
        Advec_4(Model*, Input*); ///< Constructor of the advection class.
        ~Advec_4();              ///< Destructor of the advection class.

        void exec(); ///< Execute the advection scheme.
        unsigned long get_time_limit(long unsigned int, double); ///< Get the limit on the time step imposed by the advection scheme.
        double get_cfl(double); ///< Get the CFL number.

    private:
        double calc_cfl(double*, double*, double*, double*, double); ///< Calculate the CFL number.

        template<bool>
        void advec_u(double* restrict, double* restrict, double* restrict, double* restrict, double* restrict); ///< Calculate longitudinal velocity advection.
        template<bool>
        void advec_v(double* restrict, double* restrict, double* restrict, double* restrict, double* restrict); ///< Calculate latitudinal velocity advection.
        template<bool>
        void advec_w(double* restrict, double* restrict, double* restrict, double* restrict, double* restrict); ///< Calculate vertical velocity advection.
        template<bool>
        void advec_s(double* restrict, double* restrict, double* restrict, double* restrict, double* restrict, double* restrict); ///< Calculate scalar advection.
};
#endif
