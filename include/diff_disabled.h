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

#ifndef DIFF_DISABLED
#define DIFF_DISABLED

#include "diff.h"

template<typename TF>
class Diff_disabled : public Diff<TF>
{
    public:
        Diff_disabled(Master&, Grid<TF>&, Fields<TF>&, Input&);  ///< Constructor of the diffusion class
        ~Diff_disabled();                                ///< Destructor of the diffusion class

        Diffusion_type get_switch() const;
        unsigned long get_time_limit(unsigned long, double);
        double get_dn(double);

        // Empty functions which simply pass for disabled diffusion
        void set_values() {}
        void exec_viscosity(Thermo<TF>&) {}
        void exec() {}

        //#ifdef USECUDA
        //// GPU functions and variables
        //void prepare_device() {}
        //#endif
    private:

        const Diffusion_type swdiff = Diffusion_type::Disabled;
};
#endif
