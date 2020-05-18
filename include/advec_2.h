/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
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

#ifndef ADVEC_2_H
#define ADVEC_2_H

#include "advec.h"

class Master;
template<typename> class Advec;
template<typename> class Grid;
template<typename> class Fields;
class Input;

/**
 * Derived class for 2nd order advection scheme.
 */
template<typename TF>
class Advec_2 : public Advec<TF>
{
    public:
        Advec_2(Master&, Grid<TF>&, Fields<TF>&, Input&); ///< Constructor of the advection class.
        ~Advec_2(); ///< Destructor of the advection class.

        void create(Stats<TF>&);
        void exec(Stats<TF>&); ///< Execute the advection scheme.
        unsigned long get_time_limit(long unsigned int, double); ///< Get the limit on the time step imposed by the advection scheme.
        double get_cfl(double); ///< Get the CFL number.

        void get_advec_flux(Field3d<TF>&, const Field3d<TF>&);
        Advection_type get_switch() const { return Advection_type::Advec_2; }

    private:
        using Advec<TF>::master;
        using Advec<TF>::grid;
        using Advec<TF>::fields;
        using Advec<TF>::field3d_operators;

        using Advec<TF>::cflmax;
        using Advec<TF>::cflmin;

        const std::string tend_name = "advec";
        const std::string tend_longname = "Advection";
};
#endif
