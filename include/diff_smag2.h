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

#ifndef DIFF_SMAG2_H
#define DIFF_SMAG2_H

#include "diff.h"
#include "boundary_cyclic.h"

template<typename TF>
class Diff_smag2 : public Diff<TF>
{
    public:
        Diff_smag2(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Diff_smag2();

        Diffusion_type get_switch() const;
        unsigned long get_time_limit(unsigned long, double);
        double get_dn(double);

        void set_values();
        void init();
        void exec();
        void exec_viscosity(Boundary<TF>&, Thermo<TF>&);

        //#ifdef USECUDA
        //void prepare_device() {};
        //#endif

    private:
        using Diff<TF>::master;
        using Diff<TF>::grid;
        using Diff<TF>::fields;
        Boundary_cyclic<TF> boundary_cyclic;

        using Diff<TF>::tPr;

        const Diffusion_type swdiff = Diffusion_type::Diff_smag2;

        double dnmax;
        double dnmul;

        double cs;
};
#endif
