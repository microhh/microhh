/*
 * MicroHH
 * Copyright (c) 2011-2019 Chiel van Heerwaarden
 * Copyright (c) 2011-2019 Thijs Heus
 * Copyright (c) 2014-2019 Bart van Stratum
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

#include <cstdio>
#include <cmath>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "finite_difference.h"
#include "model.h"
#include "thermo.h"
#include "diff.h"
#include "advec.h"
#include "force.h"
#include "stats.h"

#include "budget.h"
#include "budget_4.h"

template<typename TF>
Budget_4<TF>::Budget_4(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin,
        Thermo<TF>& thermoin, Diff<TF>& diffin, Advec<TF>& advecin, Force<TF>& forcein, Input& inputin) :
    Budget<TF>(masterin, gridin, fieldsin, thermoin, diffin, advecin, forcein, inputin)
{
}

template<typename TF>
Budget_4<TF>::~Budget_4()
{
}

template<typename TF>
void Budget_4<TF>::init()
{
}

template<typename TF>
void Budget_4<TF>::create(Stats<TF>& stats)
{
    // Add the profiles for the kinetic energy to the statistics.
    stats.add_prof("ke" , "Kinetic energy" , "m2 s-2", "z");
    stats.add_prof("tke", "Turbulent kinetic energy" , "m2 s-2", "z");

    // Add the profiles for the kinetic energy budget to the statistics.
    stats.add_prof("u2_shear" , "Shear production term in U2 budget" , "m2 s-3", "z" );
    stats.add_prof("v2_shear" , "Shear production term in V2 budget" , "m2 s-3", "z" );
    stats.add_prof("tke_shear", "Shear production term in TKE budget", "m2 s-3", "z" );
    stats.add_prof("uw_shear" , "Shear production term in UW budget" , "m2 s-3", "zh");

    stats.add_prof("u2_turb" , "Turbulent transport term in U2 budget" , "m2 s-3", "z" );
    stats.add_prof("v2_turb" , "Turbulent transport term in V2 budget" , "m2 s-3", "z" );
    stats.add_prof("w2_turb" , "Turbulent transport term in W2 budget" , "m2 s-3", "zh");
    stats.add_prof("tke_turb", "Turbulent transport term in TKE budget", "m2 s-3", "z" );
    stats.add_prof("uw_turb" , "Turbulent transport term in UW budget" , "m2 s-3", "zh");

    stats.add_prof("u2_visc" , "Viscous transport term in U2 budget" , "m2 s-3", "z" );
    stats.add_prof("v2_visc" , "Viscous transport term in V2 budget" , "m2 s-3", "z" );
    stats.add_prof("w2_visc" , "Viscous transport term in W2 budget" , "m2 s-3", "zh");
    stats.add_prof("tke_visc", "Viscous transport term in TKE budget", "m2 s-3", "z" );
    stats.add_prof("uw_visc" , "Viscous transport term in UW budget" , "m2 s-3", "zh");

    stats.add_prof("u2_diss" , "Dissipation term in U2 budget" , "m2 s-3", "z" );
    stats.add_prof("v2_diss" , "Dissipation term in V2 budget" , "m2 s-3", "z" );
    stats.add_prof("w2_diss" , "Dissipation term in W2 budget" , "m2 s-3", "zh");
    stats.add_prof("tke_diss", "Dissipation term in TKE budget", "m2 s-3", "z" );
    stats.add_prof("uw_diss" , "Dissipation term in UW budget" , "m2 s-3", "zh");

    stats.add_prof("w2_pres" , "Pressure transport term in W2 budget" , "m2 s-3", "zh");
    stats.add_prof("tke_pres", "Pressure transport term in TKE budget", "m2 s-3", "z" );
    stats.add_prof("uw_pres" , "Pressure transport term in UW budget" , "m2 s-3", "zh");

    stats.add_prof("u2_rdstr", "Pressure redistribution term in U2 budget", "m2 s-3", "z" );
    stats.add_prof("v2_rdstr", "Pressure redistribution term in V2 budget", "m2 s-3", "z" );
    stats.add_prof("w2_rdstr", "Pressure redistribution term in W2 budget", "m2 s-3", "zh");
    stats.add_prof("uw_rdstr", "Pressure redistribution term in UW budget", "m2 s-3", "zh");

    if (thermo.get_switch() != "0")
    {
        stats.add_prof("w2_buoy" , "Buoyancy production/destruction term in W2 budget" , "m2 s-3", "zh");
        stats.add_prof("tke_buoy", "Buoyancy production/destruction term in TKE budget", "m2 s-3", "z" );
        stats.add_prof("uw_buoy" , "Buoyancy production/destruction term in UW budget" , "m2 s-3", "zh");

        stats.add_prof("b2_shear", "Shear production term in B2 budget"   , "m2 s-5", "z");
        stats.add_prof("b2_turb" , "Turbulent transport term in B2 budget", "m2 s-5", "z");
        stats.add_prof("b2_visc" , "Viscous transport term in B2 budget"  , "m2 s-5", "z");
        stats.add_prof("b2_diss" , "Dissipation term in B2 budget"        , "m2 s-5", "z");

        stats.add_prof("bw_shear", "Shear production term in BW budget"   , "m2 s-4", "zh");
        stats.add_prof("bw_turb" , "Turbulent transport term in BW budget", "m2 s-4", "zh");
        stats.add_prof("bw_visc" , "Viscous transport term in BW budget"  , "m2 s-4", "zh");
        stats.add_prof("bw_rdstr", "Redistribution term in BW budget"     , "m2 s-4", "zh");
        stats.add_prof("bw_buoy" , "Buoyancy term in BW budget"           , "m2 s-4", "zh");
        stats.add_prof("bw_diss" , "Dissipation term in BW budget"        , "m2 s-4", "zh");
        stats.add_prof("bw_pres" , "Pressure transport term in BW budget" , "m2 s-4", "zh");
    }

    if (thermo.get_switch() != "0")
    {
        // Add the profiles for the potential energy budget to the statistics.
        stats.add_prof("bsort", "Sorted buoyancy", "m s-2", "z");
        stats.add_prof("zsort", "Height diff buoyancy and sorted buoyancy", "m", "z");
        stats.add_prof("pe"   , "Total potential energy", "m2 s-2", "z");
        stats.add_prof("ape"  , "Available potential energy", "m2 s-2", "z");
        stats.add_prof("bpe"  , "Background potential energy", "m2 s-2", "z");

        // Add the budget terms for the potential energy.
        stats.add_prof("pe_turb", "Turbulent transport term in potential energy budget", "m2 s-3", "z");
        stats.add_prof("pe_visc", "Viscous transport term in potential energy budget", "m2 s-3", "z");
        stats.add_prof("pe_bous", "Boussinesq term in potential energy budget", "m2 s-3", "z");

        // add the budget terms for the background potential energy
        // stats.add_prof("bpe_turb", "Turbulent transport term in background potential energy budget", "m2 s-3", "z");
        // stats.add_prof("bpe_visc", "Viscous transport term in background potential energy budget", "m2 s-3", "z");
        // stats.add_prof("bpe_diss", "Dissipation term in background potential energy budget", "m2 s-3", "z");
    }
}

template<typename TF>
void Budget_4<TF>::exec_stats(Stats<TF>& stats)
{
}

template class Budget_4<double>;
template class Budget_4<float>;
