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

#ifndef DIFF_4_H
#define DIFF_4_H

#include "diff.h"

template<typename TF>
class Diff_4 : public Diff<TF>
{
    public:
        Diff_4(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Diff_4();

        Diffusion_type get_switch() const;
        unsigned long get_time_limit(unsigned long, double);
        double get_dn(double);

        void create(Stats<TF>&);
        void init() {};
        void exec(Boundary<TF>&);
        void exec_stats(Stats<TF>&) {};

        // Empty functions, these are allowed to pass.
        void exec_viscosity(Boundary<TF>&, Thermo<TF>&) {}
        void diff_flux(Field3d<TF>&, const Field3d<TF>& , const int*);

        #ifdef USECUDA
        void prepare_device(Boundary<TF>&) {};
        #endif

    private:
        using Diff<TF>::master;
        using Diff<TF>::grid;
        using Diff<TF>::fields;

        const Diffusion_type swdiff = Diffusion_type::Diff_4;

        double dnmax;
        double dnmul;
};
#endif
