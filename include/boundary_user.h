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

#ifndef BOUNDARY_USER
#define BOUNDARY_USER

#include "boundary.h"

// forward declaration to reduce compilation time
class cmodel;

class cboundary_user : public cboundary
{
  public:
    cboundary_user(cmodel *, cinput *);

    void init(cinput *);

    int setvalues();

  private:
    int setbc_patch(double *, double *, double *, int, double, double, double,
                    double *, double, double); ///< Set the values for the boundary fields.

    // patch type
    int    patch_dim;
    double patch_xh;
    double patch_xr;
    double patch_xi;
    double patch_facr;
    double patch_facl;
};
#endif
