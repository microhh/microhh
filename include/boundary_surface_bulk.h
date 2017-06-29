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

#ifndef BOUNDARY_SURFACE_BULK
#define BOUNDARY_SURFACE_BULK

#include "boundary_surface.h"
#include "stats.h"

class Model;
class Input;
class Stats;
struct Mask;

class Boundary_surface_bulk : public Boundary_surface
{
    public:
        Boundary_surface_bulk(Model*, Input*);
        ~Boundary_surface_bulk();

        void init(Input*);
        void set_values();

    private:
        // surface scheme
        void update_bcs();
        void calculate_du(double*, double*, double*, double*, double*);
        void momentum_fluxgrad(double*, double*, double*, double*, double*, double*, double*, double*, double*, double, double);
        void scalar_fluxgrad(double*, double*, double*, double*, double*, double, double);
        void surface_scaling(double*, double*, double*, double*, double);

        // transfer coefficients
        double bulk_cm;
        std::map<std::string, double> bulk_cs;
};
#endif
