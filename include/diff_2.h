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

#ifndef DIFF_2_H
#define DIFF_2_H

#include "diff.h"

template<typename TF>
class Diff_2 : public Diff<TF>
{
    public:
        Diff_2(Master&, Grid<TF>&, Fields<TF>&, Boundary<TF>&, Input&); ///< Constructor of the diffusion class
        ~Diff_2(); ///< Destructor of the diffusion class

        Diffusion_type get_switch() const;
        unsigned long get_time_limit(unsigned long, double);
        double get_dn(double);

        void create(Stats<TF>&);
        void init() {};
        void exec(Stats<TF>&);

        // Empty functions, these are allowed to pass.
        void exec_viscosity(Stats<TF>&, Thermo<TF>&) {}
        void diff_flux(Field3d<TF>&, const Field3d<TF>&);
        void exec_stats(Stats<TF>&, Thermo<TF>&) {};

        #ifdef USECUDA
        void prepare_device(Boundary<TF>&) {};
        #endif

    private:
        using Diff<TF>::master;
        using Diff<TF>::grid;
        using Diff<TF>::fields;

        const Diffusion_type swdiff = Diffusion_type::Diff_2;

        double dnmax;
        double dnmul;
        const std::string tend_name = "diff";
        const std::string tend_longname = "Diffusion";
};
#endif
