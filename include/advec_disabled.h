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

#ifndef ADVEC_DISABLED
#define ADVEC_DISABLED

#include "advec.h"

class Model;
class Input;

/**
 * Derived class for a disabled advection scheme
 */
class Advec_disabled : public Advec
{
    public:
        Advec_disabled(Model*, Input*); ///< Constructor of the advection class.
        ~Advec_disabled();              ///< Destructor of the advection class.

        void exec(); ///< Execute the advection scheme.

        unsigned long get_time_limit(unsigned long, double); ///< Get the maximum time step imposed by advection scheme

        double get_cfl(double); ///< Retrieve the CFL number.
};
#endif
