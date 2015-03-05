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

#ifndef ADVEC_4M
#define ADVEC_4M

#include "advec.h"

/**
 * Derived class for fully conservative 4th order advection scheme.
 * Fully mass, momentum and energy conserving advection scheme based on the paper
 * of Morinishi et al., (1998).
 */
class Advec_4m : public Advec
{
    public:
        Advec_4m(Model*, Input*); ///< Constructor of the advection class.
        ~Advec_4m();              ///< Destructor of the advection class.

        void exec(); ///< Execute the advection scheme.

        unsigned long get_time_limit(long unsigned int, double); ///< Get the limit on the time step imposed by the advection scheme.

        double get_cfl(double); ///< Get the CFL number.

    private:
        double calc_cfl(double*, double*, double*, double*, double); ///< Calculate the CFL number.

        void advec_u(double*, double*, double*, double*, double*);          ///< Calculate longitudinal velocity advection.
        void advec_v(double*, double*, double*, double*, double*);          ///< Calculate latitudinal velocity advection.
        void advec_w(double*, double*, double*, double*, double*);          ///< Calculate vertical velocity advection.
        void advec_s(double*, double*, double*, double*, double*, double*); ///< Calculate scalar advection.
};
#endif
