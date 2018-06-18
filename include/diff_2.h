/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
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

#ifndef DIFF_2_H
#define DIFF_2_H

#include "diff.h"

template<typename TF>
class Diff_2 : public Diff<TF>
{
    public:
        Diff_2(Master&, Grid<TF>&, Fields<TF>&, Input&);  ///< Constructor of the diffusion class
        ~Diff_2();                                ///< Destructor of the diffusion class

        Diffusion_type get_switch() const;
        unsigned long get_time_limit(unsigned long, double);
        double get_dn(double);

        void set_values();
        void init() {};
        void exec(Boundary<TF>&);

        // Empty functions, these are allowed to pass.
        void exec_viscosity(Boundary<TF>&, Thermo<TF>&) {}

        #ifdef USECUDA
        void prepare_device() {};
        #endif

    private:
        using Diff<TF>::master;
        using Diff<TF>::grid;
        using Diff<TF>::fields;

        const Diffusion_type swdiff = Diffusion_type::Diff_2;

        double dnmax;
        double dnmul;
};
#endif
