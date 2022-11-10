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

#ifndef DIFF_DEARDORFF_H
#define DIFF_DEARDORFF_H

#include "diff.h"
#include "boundary_cyclic.h"
#include "field3d_operators.h"

template<typename> class Stats;


template<typename TF>
class Diff_deardorff : public Diff<TF>
{
    public:
        Diff_deardorff(Master&, Grid<TF>&, Fields<TF>&, Boundary<TF>&, Input&);
        ~Diff_deardorff();

        Diffusion_type get_switch() const;
        unsigned long get_time_limit(unsigned long, double);
        double get_dn(double);

        void create(Stats<TF>&);
        void init();
        void exec(Stats<TF>&);
        void exec_viscosity(Stats<TF>&, Thermo<TF>&); 
        void diff_flux(Field3d<TF>&, const Field3d<TF>&);
        void exec_stats(Stats<TF>&, Thermo<TF>&); // Pass Thermo<TF>& also in other diffusion classes

        #ifdef USECUDA
        void prepare_device(Boundary<TF>&);
        void clear_device();
        #endif

    private:
        using Diff<TF>::master;
        using Diff<TF>::grid;
        using Diff<TF>::fields;
        using Diff<TF>::boundary;
        Boundary_cyclic<TF> boundary_cyclic;
        Field3d_operators<TF> field3d_operators;

        using Diff<TF>::tPr;

        const Diffusion_type swdiff = Diffusion_type::Diff_deardorff;

        void create_stats(Stats<TF>&);

        TF* mlen_g;

        double dnmax;
        double dnmul;

        TF cs;
        TF ap;
        TF cf;
        TF ce1;
        TF ce2;
        TF cm;
        TF ch1;
        TF ch2;
        TF cn;

        bool sw_buoy;
        bool swmason;  ///< Switch for use of Mason's wall correction

        const std::string tend_name = "diff";
        const std::string tend_longname = "Diffusion";

        // Specific tendency names for sgs tke-scheme
        const std::string tend_name_shear = "shear";
        const std::string tend_longname_shear = "SGS TKE Shear Term";

        const std::string tend_name_buoy = "buoy";
        const std::string tend_longname_buoy = "SGS TKE Buoyancy Term";

        const std::string tend_name_diss = "diss";
        const std::string tend_longname_diss = "SGS TKE Dissipation";
};
#endif
