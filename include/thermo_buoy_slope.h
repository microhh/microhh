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

#ifndef THERMO_BUOY_SLOPE
#define THERMO_BUOY_SLOPE

#include "thermo.h"

class Master;
class Grid;
class Fields;

/**
 * Class for the dry thermodynamics.
 * This class is responsible for the computation of the right hand side term related to
 * the acceleration by buoyancy. In the dry thermodynamics temperature and buoyancy are
 * equivalent and no complex buoyancy function is required.
 */
class Thermo_buoy_slope : public Thermo
{
    public:
        Thermo_buoy_slope(Model *, Input *); ///< Constructor of the dry thermodynamics class.
        virtual ~Thermo_buoy_slope();        ///< Destructor of the dry thermodynamics class.

        void exec(); ///< Add the tendencies belonging to the buoyancy.

        bool check_thermo_field(std::string name);
        void get_buoyancy_surf(Field3d *);             ///< Compute the near-surface and bottom buoyancy for usage in another routine.
        void get_buoyancy_fluxbot(Field3d*);           ///< Compute the bottom buoyancy flux for usage in another routine.
        void get_prog_vars(std::vector<std::string>*); ///< Retrieve a list of prognostic variables.
        void get_thermo_field(Field3d*, Field3d*, std::string name);

        // Empty functions that are allowed to pass.
        void init() {}
        void create(Input*) {}
        void exec_stats(Mask*) {}
        void exec_cross() {}
        void exec_dump() {}
        void get_mask(Field3d*, Field3d*, Mask*) {}

    private:
        void calc_buoyancy(double *, double *); ///< Calculation of the buoyancy.
        void calc_buoyancy_bot(double *, double *,
                               double *, double *); ///< Calculation of the near-surface and surface buoyancy.
        void calc_buoyancy_fluxbot(double *, double *); ///< Calculation of the buoyancy flux at the bottom.
        void calc_buoyancy_tend_u_4th(double *, double *); ///< Calculation of the buoyancy tendency with 4th order accuracy.
        void calc_buoyancy_tend_w_4th(double *, double *); ///< Calculation of the buoyancy tendency with 4th order accuracy.
        void calc_buoyancy_tend_b_4th(double *, double *, double *); ///< Calculation of the buoyancy tendency with 4th order accuracy.

        double alpha; ///< Slope angle in radians.
        double n2;    ///< Background stratification.
};
#endif
