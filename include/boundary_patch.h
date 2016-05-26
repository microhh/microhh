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

#ifndef BOUNDARY_PATCH
#define BOUNDARY_PATCH

#include "boundary_surface.h"

class Model;
class Input;
class Fields;
struct Mask;

class Boundary_patch : public Boundary_surface
{
    public:
        Boundary_patch(Model*, Input*);

        void init(Input*);
        void set_values();
        void get_mask(Field3d*, Field3d*, Mask*);
        void get_surface_mask(Field3d*);

    private:
        void calc_patch(double*, const double*, const double*, int, double, double, double, double, double, double, double, double);  ///< Calculate the patches
        void set_bc_patch(double*, double*, double*, double*, const double, const double, const int, const double, const double, const double);       ///< Set the values for the boundary fields.

        // Patch properties.
        int    patch_dim;
        double patch_xh;
        double patch_xr;
        double patch_xi;
        double patch_xoffs;
        double patch_yh;
        double patch_yr;
        double patch_yi;
        double patch_yoffs;
        double patch_facr;
        double patch_facl;

        std::map<std::string, double> patch_facr_map;
        std::map<std::string, double> patch_facl_map;

};
#endif
