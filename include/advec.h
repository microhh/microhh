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

#ifndef ADVEC
#define ADVEC

#include <string>

// forward declarations to speed up build time
class cmaster;
class cinput;
class cmodel;
class cgrid;
class cfields;

/**
 * Base class for the advection scheme.
 * This class handles the case when advection is turned off. Derived classes are
 * implemented that handle different advection schemes.
 */
class cadvec
{
  public:
    cadvec(cmodel *, cinput *); ///< Constructor of the advection class.
    virtual ~cadvec();          ///< Destructor of the advection class.

    static cadvec* factory(cmaster *, cinput *, cmodel *, const std::string); ///< Factory function for advection class generation.

    virtual int exec(); ///< Execute the advection scheme.

    virtual double getcfl(double); ///< Retrieve the CFL number.

    virtual unsigned long gettimelim(unsigned long, double); ///< Get the maximum time step imposed by advection scheme

  protected:
    cmaster *master; ///< Pointer to master class.
    cmodel  *model;  ///< Pointer to model class.
    cgrid   *grid;   ///< Pointer to grid class.
    cfields *fields; ///< Pointer to fields class.

    double cflmax; ///< Maximum allowed value for the CFL criterion.
};
#endif
