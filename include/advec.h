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

#ifndef ADVEC
#define ADVEC

#include <string>

class Master;
class Input;
class Model;
class Grid;
class Fields;
class Input;

/**
 * Base class for the advection scheme.
 * This class handles the case when advection is turned off. Derived classes are
 * implemented that handle different advection schemes.
 */
class Advec
{
    public:
        Advec(Model*, Input*); ///< Constructor of the advection class.
        virtual ~Advec();      ///< Destructor of the advection class.

        static Advec* factory(Master*, Input*, Model*, const std::string); ///< Factory function for advection class generation.

        // Pure virtual functions that have to be implemented in derived class.
        virtual void exec() = 0; ///< Execute the advection scheme.
        virtual unsigned long get_time_limit(unsigned long, double) = 0; ///< Get the maximum time step imposed by advection scheme
        virtual double get_cfl(double) = 0; ///< Retrieve the CFL number.

    protected:
        Master* master; ///< Pointer to master class.
        Model*  model;  ///< Pointer to model class.
        Grid*   grid;   ///< Pointer to grid class.
        Fields* fields; ///< Pointer to fields class.

        double cflmax; ///< Maximum allowed value for the CFL criterion.
        static const double cflmin; ///< Minimum value for CFL used to avoid overflows.
};
#endif
