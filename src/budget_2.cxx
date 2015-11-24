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

#include <cstdio>
#include <cmath>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "finite_difference.h"
#include "model.h"
#include "thermo.h"
#include "stats.h"
#include <netcdfcpp.h>

#include "budget.h"
#include "budget_2.h"

using namespace Finite_difference::O2;

Budget_2::Budget_2(Input* inputin, Master* masterin, Grid* gridin, Fields* fieldsin, Thermo* thermoin, Stats* statsin) :
    Budget(inputin, masterin, gridin, fieldsin, thermoin, statsin)
{
    umodel = 0;
    vmodel = 0;
}

Budget_2::~Budget_2()
{
    delete[] umodel;
    delete[] vmodel;
}

void Budget_2::init()
{
    umodel = new double[grid.kcells];
    vmodel = new double[grid.kcells];

    for (int k=0; k<grid.kcells; ++k)
    {
        umodel[k] = 0.;
        vmodel[k] = 0.;
    }
}

void Budget_2::create()
{
    // add the profiles for the kinetic energy to the statistics
    stats.add_prof("ke" , "Kinetic energy" , "m2 s-2", "z");
    stats.add_prof("tke", "Turbulent kinetic energy" , "m2 s-2", "z");

    // add the profiles for the kinetic energy budget to the statistics
    stats.add_prof("u2_shear" , "Shear production term in U2 budget" , "m2 s-3", "z" );
    stats.add_prof("v2_shear" , "Shear production term in V2 budget" , "m2 s-3", "z" );
    stats.add_prof("tke_shear", "Shear production term in TKE budget", "m2 s-3", "z" );

    stats.add_prof("u2_turb" , "Turbulent transport term in U2 budget" , "m2 s-3", "z" );
    stats.add_prof("v2_turb" , "Turbulent transport term in V2 budget" , "m2 s-3", "z" );
    stats.add_prof("w2_turb" , "Turbulent transport term in W2 budget" , "m2 s-3", "zh");
    stats.add_prof("tke_turb", "Turbulent transport term in TKE budget", "m2 s-3", "z" );

    stats.add_prof("u2_visc" , "Viscous transport term in U2 budget" , "m2 s-3", "z" );
    stats.add_prof("v2_visc" , "Viscous transport term in V2 budget" , "m2 s-3", "z" );
    stats.add_prof("w2_visc" , "Viscous transport term in W2 budget" , "m2 s-3", "zh");
    stats.add_prof("tke_visc", "Viscous transport term in TKE budget", "m2 s-3", "z" );

    stats.add_prof("u2_diss" , "Dissipation term in U2 budget" , "m2 s-3", "z" );
    stats.add_prof("v2_diss" , "Dissipation term in V2 budget" , "m2 s-3", "z" );
    stats.add_prof("w2_diss" , "Dissipation term in W2 budget" , "m2 s-3", "zh");
    stats.add_prof("tke_diss", "Dissipation term in TKE budget", "m2 s-3", "z" );

    stats.add_prof("w2_pres" , "Pressure transport term in W2 budget" , "m2 s-3", "zh");
    stats.add_prof("tke_pres", "Pressure transport term in TKE budget", "m2 s-3", "z" );

    stats.add_prof("u2_rdstr", "Pressure redistribution term in U2 budget", "m2 s-3", "z" );
    stats.add_prof("v2_rdstr", "Pressure redistribution term in V2 budget", "m2 s-3", "z" );
    stats.add_prof("w2_rdstr", "Pressure redistribution term in W2 budget", "m2 s-3", "zh");

    if (thermo.get_switch() != "0")
    {
        stats.add_prof("w2_buoy" , "Buoyancy production/destruction term in W2 budget" , "m2 s-3", "zh");
        stats.add_prof("tke_buoy", "Buoyancy production/destruction term in TKE budget", "m2 s-3", "z" );
    }
}

void Budget_2::exec_stats(Mask* m)
{
    // calculate the mean of the fields
    grid.calc_mean(umodel, fields.u->data, grid.kcells);
    grid.calc_mean(vmodel, fields.v->data, grid.kcells);

    // Calculate kinetic and turbulent kinetic energy
    calc_ke(m->profs["ke"].data, m->profs["tke"].data,
            fields.u->data, fields.v->data, fields.w->data, umodel, vmodel, grid.utrans, grid.vtrans);
}

void Budget_2::calc_ke(double* const restrict ke, double* const restrict tke,
                       const double* const restrict u, const double* const restrict v, const double* const restrict w,
                       const double* const restrict umodel, const double* const restrict vmodel,
                       const double utrans, const double vtrans)
{
    const int ii = 1;
    const int jj = grid.icells;
    const int kk = grid.ijcells;

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        ke[k]  = 0;
        tke[k] = 0;
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                const double u2 = interp2(pow(u[ijk]+utrans, 2), pow(u[ijk+ii]+utrans, 2));
                const double v2 = interp2(pow(v[ijk]+vtrans, 2), pow(v[ijk+jj]+vtrans, 2));
                const double w2 = interp2(pow(w[ijk],        2), pow(w[ijk+ii],        2));

                ke[k] += 0.5 * (u2 + v2 + w2);
            }

        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                const double u2 = interp2(pow(u[ijk]-umodel[k], 2), pow(u[ijk+ii]-umodel[k], 2));
                const double v2 = interp2(pow(v[ijk]-vmodel[k], 2), pow(v[ijk+jj]-vmodel[k], 2));
                const double w2 = interp2(pow(w[ijk],           2), pow(w[ijk+ii],           2));

                tke[k] += 0.5 * (u2 + v2 + w2);
            }
    }

    master.sum(ke , grid.kcells);
    master.sum(tke, grid.kcells);

    const int n = grid.itot*grid.jtot;
    for (int k=grid.kstart; k<grid.kend; ++k)
    {
       ke[k]  /= n;
       tke[k] /= n;
    }
}
