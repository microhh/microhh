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

#include "budget.h"
#include "budget_4.h"

using namespace Finite_difference::O4;

Budget_4::Budget_4(Input* inputin, Master* masterin, Grid* gridin, Fields* fieldsin, Thermo* thermoin, Diff* diffin, Advec* advecin, Force* forcein, Stats* statsin) :
    Budget(inputin, masterin, gridin, fieldsin, thermoin, diffin, advecin, forcein, statsin)
{
    umodel = 0;
    vmodel = 0;
}

Budget_4::~Budget_4()
{
    delete[] umodel;
    delete[] vmodel;
}

void Budget_4::init()
{
    umodel = new double[grid.kcells];
    vmodel = new double[grid.kcells];

    for (int k=0; k<grid.kcells; ++k)
    {
        umodel[k] = 0.;
        vmodel[k] = 0.;
    }
}

void Budget_4::create()
{
    // add the profiles for the kinetic energy to the statistics
    stats.add_prof("ke" , "Kinetic energy" , "m2 s-2", "z");
    stats.add_prof("tke", "Turbulent kinetic energy" , "m2 s-2", "z");

    // add the profiles for the kinetic energy budget to the statistics
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
        // add the profiles for the potential energy budget to the statistics
        stats.add_prof("bsort", "Sorted buoyancy", "m s-2", "z");
        stats.add_prof("zsort", "Height diff buoyancy and sorted buoyancy", "m", "z");
        stats.add_prof("pe"   , "Total potential energy", "m2 s-2", "z");
        stats.add_prof("ape"  , "Available potential energy", "m2 s-2", "z");
        stats.add_prof("bpe"  , "Background potential energy", "m2 s-2", "z");

        // add the budget terms for the potential energy
        stats.add_prof("pe_turb", "Turbulent transport term in potential energy budget", "m2 s-3", "z");
        stats.add_prof("pe_visc", "Viscous transport term in potential energy budget", "m2 s-3", "z");
        stats.add_prof("pe_bous", "Boussinesq term in potential energy budget", "m2 s-3", "z");

        // add the budget terms for the background potential energy
        // stats.add_prof("bpe_turb", "Turbulent transport term in background potential energy budget", "m2 s-3", "z");
        // stats.add_prof("bpe_visc", "Viscous transport term in background potential energy budget", "m2 s-3", "z");
        // stats.add_prof("bpe_diss", "Dissipation term in background potential energy budget", "m2 s-3", "z");
    }
}

void Budget_4::exec_stats(Mask* m)
{
    // calculate the mean of the fields
    grid.calc_mean(umodel, fields.u->data, grid.kcells);
    grid.calc_mean(vmodel, fields.v->data, grid.kcells);

    if (grid.swspatialorder == "4")
    {
        // calculate the TKE budget
        calc_ke(fields.u->data, fields.v->data, fields.w->data,
                umodel, vmodel,
                grid.utrans, grid.vtrans,
                m->profs["ke"].data, m->profs["tke"].data);

        calc_tke_budget_shear_turb(fields.u->data, fields.v->data, fields.w->data,
                                   fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data,
                                   umodel, vmodel,
                                   m->profs["u2_shear"].data, m->profs["v2_shear"].data, m->profs["tke_shear"].data, m->profs["uw_shear"].data,
                                   m->profs["u2_turb"].data, m->profs["v2_turb"].data, m->profs["w2_turb"].data, m->profs["tke_turb"].data, m->profs["uw_turb"].data,
                                   grid.dzi4, grid.dzhi4);

        calc_tke_budget(fields.u->data, fields.v->data, fields.w->data, fields.sd["p"]->data,
                        fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data,
                        umodel, vmodel,
                        m->profs["u2_visc"].data, m->profs["v2_visc"].data, m->profs["w2_visc"].data, m->profs["tke_visc"].data, m->profs["uw_visc"].data,
                        m->profs["u2_diss"].data, m->profs["v2_diss"].data, m->profs["w2_diss"].data, m->profs["tke_diss"].data, m->profs["uw_diss"].data,
                        m->profs["w2_pres"].data, m->profs["tke_pres"].data, m->profs["uw_pres"].data,
                        m->profs["u2_rdstr"].data, m->profs["v2_rdstr"].data, m->profs["w2_rdstr"].data, m->profs["uw_rdstr"].data,
                        grid.dzi4, grid.dzhi4, fields.visc);

        // calculate the buoyancy term of the TKE budget
        if (thermo.get_switch() != "0")
        {
            // store the buoyancy in the tmp1 field
            thermo.get_thermo_field(fields.atmp["tmp1"], fields.atmp["tmp2"], "b", true);

            grid.calc_mean(fields.atmp["tmp1"]->datamean, fields.atmp["tmp1"]->data, grid.kcells);
            grid.calc_mean(fields.sd["p"]->datamean, fields.sd["p"]->data, grid.kcells);

            calc_tke_budget_buoy(fields.u->data, fields.w->data, fields.atmp["tmp1"]->data,
                                 umodel, fields.atmp["tmp1"]->datamean,
                                 m->profs["w2_buoy"].data, m->profs["tke_buoy"].data, m->profs["uw_buoy"].data);

            calc_b2_budget(fields.w->data, fields.atmp["tmp1"]->data,
                           fields.atmp["tmp1"]->datamean,
                           m->profs["b2_shear"].data, m->profs["b2_turb"].data, m->profs["b2_visc"].data, m->profs["b2_diss"].data,
                           grid.dzi4, grid.dzhi4,
                           fields.visc);

            calc_bw_budget(fields.w->data, fields.sd["p"]->data, fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data,
                           fields.sd["p"]->datamean, fields.atmp["tmp1"]->datamean,
                           m->profs["bw_shear"].data, m->profs["bw_turb"].data, m->profs["bw_visc"].data,
                           m->profs["bw_buoy"].data, m->profs["bw_rdstr"].data, m->profs["bw_diss"].data, m->profs["bw_pres"].data,
                           grid.dzi4, grid.dzhi4,
                           fields.visc);
        }

        // calculate the potential energy budget
        if (thermo.get_switch() != "0")
        {
            // calculate the sorted buoyancy profile, tmp1 still contains the buoyancy
            stats.calc_sorted_prof(fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, m->profs["bsort"].data);

            // calculate the potential energy back, tmp1 contains the buoyancy, tmp2 will contain height that the local buoyancy
            // will reach in the sorted profile
            calc_pe(fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, fields.atmp["tmp2"]->databot, fields.atmp["tmp2"]->datatop,
                    grid.z,
                    m->profs["bsort"].data,
                    m->profs["pe"].data, m->profs["ape"].data, m->profs["bpe"].data,
                    m->profs["zsort"].data);


            // calculate the budget of background potential energy, start with this one, because tmp2 contains the needed height
            // which will be overwritten inside of the routine
            // calcBpeBudget(fields.w->data, fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, fields.atmp["tmp2"]->databot, fields.atmp["tmp2"]->datatop,
            //               m->profs["bpe_turb"].data, m->profs["bpe_visc"].data, m->profs["bpe_diss"].data,
            //               // TODO put the correct value for visc here!!!!!
            //               m->profs["bsort"].data,
            //               grid.z, grid.dzi4, grid.dzhi4,
            //               fields.visc);

            // calculate the budget of potential energy
            calc_pe_budget(fields.w->data, fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, fields.atmp["tmp2"]->datatop,
                           m->profs["pe_turb"].data, m->profs["pe_visc"].data, m->profs["pe_bous"].data,
                           // TODO put the correct value for visc here!!!!!
                           grid.z, grid.zh, grid.dzi4, grid.dzhi4,
                           fields.visc);
        }
    }
}

void Budget_4::calc_ke(double* restrict u, double* restrict v, double* restrict w,
                       double* restrict umodel, double* restrict vmodel,
                       double utrans, double vtrans,
                       double* restrict ke, double* restrict tke)
{
    double u2,v2,w2;

    const int ii1 = 1;
    const int ii2 = 2;
    const int jj1 = 1*grid.icells;
    const int jj2 = 2*grid.icells;
    const int kk1 = 1*grid.ijcells;
    const int kk2 = 2*grid.ijcells;

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        ke [k] = 0;
        tke[k] = 0;
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                u2 = ci0*std::pow(u[ijk-ii1] + utrans, 2) + ci1*std::pow(u[ijk    ] + utrans, 2)
                   + ci2*std::pow(u[ijk+ii1] + utrans, 2) + ci3*std::pow(u[ijk+ii2] + utrans, 2);
                v2 = ci0*std::pow(v[ijk-jj1] + vtrans, 2) + ci1*std::pow(v[ijk    ] + vtrans, 2)
                   + ci2*std::pow(v[ijk+jj1] + vtrans, 2) + ci3*std::pow(v[ijk+jj2] + vtrans, 2);
                w2 = ci0*std::pow(w[ijk-kk1], 2) + ci1*std::pow(w[ijk], 2) + ci2*std::pow(w[ijk+kk1], 2) + ci3*std::pow(w[ijk+kk2], 2);
                ke[k] += 0.5*(u2 + v2 + w2);
            }

        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                u2 = ci0*std::pow(u[ijk-ii1] - umodel[k], 2) + ci1*std::pow(u[ijk    ] - umodel[k], 2)
                   + ci2*std::pow(u[ijk+ii1] - umodel[k], 2) + ci3*std::pow(u[ijk+ii2] - umodel[k], 2);
                v2 = ci0*std::pow(v[ijk-jj1] - vmodel[k], 2) + ci1*std::pow(v[ijk    ] - vmodel[k], 2)
                   + ci2*std::pow(v[ijk+jj1] - vmodel[k], 2) + ci3*std::pow(v[ijk+jj2] - vmodel[k], 2);
                w2 = ci0*std::pow(w[ijk-kk1], 2) + ci1*std::pow(w[ijk], 2) + ci2*std::pow(w[ijk+kk1], 2) + ci3*std::pow(w[ijk+kk2], 2);
                tke[k] += 0.5*(u2 + v2 + w2);
            }
    }

    master.sum(ke , grid.kcells);
    master.sum(tke, grid.kcells);

    int n = grid.itot*grid.jtot;
    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        ke [k] /= n;
        tke[k] /= n;
    }
}

void Budget_4::calc_tke_budget_shear_turb(double* restrict u, double* restrict v, double* restrict w,
                                          double* restrict wx, double* restrict wy,
                                          double* restrict umean, double* restrict vmean,
                                          double* restrict u2_shear, double* restrict v2_shear, double* restrict tke_shear, double* restrict uw_shear,
                                          double* restrict u2_turb, double* restrict v2_turb, double* restrict w2_turb, double* restrict tke_turb, double* restrict uw_turb,
                                          double* restrict dzi4, double* restrict dzhi4)
{
    // 1. INTERPOLATE THE VERTICAL VELOCITY TO U AND V LOCATION
    const int wloc [3] = {0,0,1};
    const int wxloc[3] = {1,0,1};
    const int wyloc[3] = {0,1,1};

    grid.interpolate_4th(wx, w, wloc, wxloc);
    grid.interpolate_4th(wy, w, wloc, wyloc);

    const int jj1 = 1*grid.icells;
    const int kk1 = 1*grid.ijcells;
    const int kk2 = 2*grid.ijcells;
    const int kk3 = 3*grid.ijcells;

    double n = grid.itot*grid.jtot;

    // 2. CALCULATE THE SHEAR TERM u'w*dumean/dz
    // bottom boundary
    int k = grid.kstart;
    u2_shear [k] = 0.;
    v2_shear [k] = 0.;
    tke_shear[k] = 0.;

    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            u2_shear[k] -= 2.*(u[ijk]-umean[k])*(ci0*wx[ijk-kk1] + ci1*wx[ijk] + ci2*wx[ijk+kk1] + ci3*wx[ijk+kk2])
                         * ( cg0*(bi0*umean[k-2] + bi1*umean[k-1] + bi2*umean[k  ] + bi3*umean[k+1])
                           + cg1*(ci0*umean[k-2] + ci1*umean[k-1] + ci2*umean[k  ] + ci3*umean[k+1])
                           + cg2*(ci0*umean[k-1] + ci1*umean[k  ] + ci2*umean[k+1] + ci3*umean[k+2])
                           + cg3*(ci0*umean[k  ] + ci1*umean[k+1] + ci2*umean[k+2] + ci3*umean[k+3])) * dzi4[k];

            v2_shear[k] -= 2.*(v[ijk]-vmean[k])*(ci0*wy[ijk-kk1] + ci1*wy[ijk] + ci2*wy[ijk+kk1] + ci3*wy[ijk+kk2])
                         * ( cg0*(bi0*vmean[k-2] + bi1*vmean[k-1] + bi2*vmean[k  ] + bi3*vmean[k+1])
                           + cg1*(ci0*vmean[k-2] + ci1*vmean[k-1] + ci2*vmean[k  ] + ci3*vmean[k+1])
                           + cg2*(ci0*vmean[k-1] + ci1*vmean[k  ] + ci2*vmean[k+1] + ci3*vmean[k+2])
                           + cg3*(ci0*vmean[k  ] + ci1*vmean[k+1] + ci2*vmean[k+2] + ci3*vmean[k+3])) * dzi4[k];
        }
    tke_shear[k] += 0.5*(u2_shear[k] + v2_shear[k]);

    // interior
    for (int k=grid.kstart+1; k<grid.kend-1; ++k)
    {
        u2_shear [k] = 0.;
        v2_shear [k] = 0.;
        tke_shear[k] = 0.;
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                u2_shear[k] -= 2.*(u[ijk]-umean[k])*(ci0*wx[ijk-kk1] + ci1*wx[ijk] + ci2*wx[ijk+kk1] + ci3*wx[ijk+kk2])
                             * ( cg0*(ci0*umean[k-3] + ci1*umean[k-2] + ci2*umean[k-1] + ci3*umean[k  ])
                               + cg1*(ci0*umean[k-2] + ci1*umean[k-1] + ci2*umean[k  ] + ci3*umean[k+1])
                               + cg2*(ci0*umean[k-1] + ci1*umean[k  ] + ci2*umean[k+1] + ci3*umean[k+2])
                               + cg3*(ci0*umean[k  ] + ci1*umean[k+1] + ci2*umean[k+2] + ci3*umean[k+3])) * dzi4[k];

                v2_shear[k] -= 2.*(v[ijk]-vmean[k])*(ci0*wy[ijk-kk1] + ci1*wy[ijk] + ci2*wy[ijk+kk1] + ci3*wy[ijk+kk2])
                             * ( cg0*(ci0*vmean[k-3] + ci1*vmean[k-2] + ci2*vmean[k-1] + ci3*vmean[k  ])
                               + cg1*(ci0*vmean[k-2] + ci1*vmean[k-1] + ci2*vmean[k  ] + ci3*vmean[k+1])
                               + cg2*(ci0*vmean[k-1] + ci1*vmean[k  ] + ci2*vmean[k+1] + ci3*vmean[k+2])
                               + cg3*(ci0*vmean[k  ] + ci1*vmean[k+1] + ci2*vmean[k+2] + ci3*vmean[k+3])) * dzi4[k];
            }
        tke_shear[k] += 0.5*(u2_shear[k] + v2_shear[k]);
    }

    // top boundary
    k = grid.kend-1;

    u2_shear [k] = 0.;
    v2_shear [k] = 0.;
    tke_shear[k] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            u2_shear[k] -= 2.*(u[ijk]-umean[k])*(ci0*wx[ijk-kk1] + ci1*wx[ijk] + ci2*wx[ijk+kk1] + ci3*wx[ijk+kk2])
                         * ( cg0*(ci0*umean[k-3] + ci1*umean[k-2] + ci2*umean[k-1] + ci3*umean[k  ])
                           + cg1*(ci0*umean[k-2] + ci1*umean[k-1] + ci2*umean[k  ] + ci3*umean[k+1])
                           + cg2*(ci0*umean[k-1] + ci1*umean[k  ] + ci2*umean[k+1] + ci3*umean[k+2])
                           + cg3*(ti0*umean[k  ] + ti1*umean[k+1] + ti2*umean[k+2] + ti3*umean[k+3])) * dzi4[k];

            v2_shear[k] -= 2.*(v[ijk]-vmean[k])*(ci0*wy[ijk-kk1] + ci1*wy[ijk] + ci2*wy[ijk+kk1] + ci3*wy[ijk+kk2])
                         * ( cg0*(ci0*vmean[k-3] + ci1*vmean[k-2] + ci2*vmean[k-1] + ci3*vmean[k  ])
                           + cg1*(ci0*vmean[k-2] + ci1*vmean[k-1] + ci2*vmean[k  ] + ci3*vmean[k+1])
                           + cg2*(ci0*vmean[k-1] + ci1*vmean[k  ] + ci2*vmean[k+1] + ci3*vmean[k+2])
                           + cg3*(ti0*vmean[k-1] + ti1*vmean[k  ] + ti2*vmean[k+1] + ti3*vmean[k+2])) * dzi4[k];
        }
    tke_shear[k] += 0.5*(u2_shear[k] + v2_shear[k]);

    // create the profiles
    master.sum(u2_shear, grid.kcells);
    master.sum(v2_shear, grid.kcells);
    master.sum(tke_shear, grid.kcells);

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        u2_shear [k] /= n;
        v2_shear [k] /= n;
        tke_shear[k] /= n;
    }

    // Reynolds stresses
    uw_shear [k] = 0.;

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        uw_shear [k] = 0.;
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                uw_shear[k] -= ( std::pow(wx[ijk],2) * ( cg0*umean[k-2] + cg1*umean[k-1] + cg2*umean[k] + cg3*umean[k+1] ) ) * dzhi4[k];
            }
    }

    master.sum(uw_shear, grid.kcells);

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        uw_shear [k] /= n;
    }


    // 3. CALCULATE TURBULENT FLUXES
    // bottom boundary
    k = grid.kstart;

    u2_turb [k] = 0.;
    v2_turb [k] = 0.;
    tke_turb[k] = 0.;
    uw_turb [k] = 0.;

    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            u2_turb[k]  -= ( cg0*((bi0*std::pow(u[ijk-kk2]-umean[k-2],2) + bi1*std::pow(u[ijk-kk1]-umean[k-1],2) + bi2*std::pow(u[ijk    ]-umean[k  ],2) + bi3*std::pow(u[ijk+kk1]-umean[k+1],2))*wx[ijk-kk1])
                           + cg1*((ci0*std::pow(u[ijk-kk2]-umean[k-2],2) + ci1*std::pow(u[ijk-kk1]-umean[k-1],2) + ci2*std::pow(u[ijk    ]-umean[k  ],2) + ci3*std::pow(u[ijk+kk1]-umean[k+1],2))*wx[ijk    ])
                           + cg2*((ci0*std::pow(u[ijk-kk1]-umean[k-1],2) + ci1*std::pow(u[ijk    ]-umean[k  ],2) + ci2*std::pow(u[ijk+kk1]-umean[k+1],2) + ci3*std::pow(u[ijk+kk2]-umean[k+2],2))*wx[ijk+kk1])
                           + cg3*((ci0*std::pow(u[ijk    ]-umean[k  ],2) + ci1*std::pow(u[ijk+kk1]-umean[k+1],2) + ci2*std::pow(u[ijk+kk2]-umean[k+2],2) + ci3*std::pow(u[ijk+kk3]-umean[k+3],2))*wx[ijk+kk2]) ) * dzi4[k];

            v2_turb[k]  -= ( cg0*((bi0*std::pow(v[ijk-kk2]-vmean[k-2],2) + bi1*std::pow(v[ijk-kk1]-vmean[k-1],2) + bi2*std::pow(v[ijk    ]-vmean[k  ],2) + bi3*std::pow(v[ijk+kk1]-vmean[k+1],2))*wy[ijk-kk1])
                           + cg1*((ci0*std::pow(v[ijk-kk2]-vmean[k-2],2) + ci1*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci2*std::pow(v[ijk    ]-vmean[k  ],2) + ci3*std::pow(v[ijk+kk1]-vmean[k+1],2))*wy[ijk    ])
                           + cg2*((ci0*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci1*std::pow(v[ijk    ]-vmean[k  ],2) + ci2*std::pow(v[ijk+kk1]-vmean[k+1],2) + ci3*std::pow(v[ijk+kk2]-vmean[k+2],2))*wy[ijk+kk1])
                           + cg3*((ci0*std::pow(v[ijk    ]-vmean[k  ],2) + ci1*std::pow(v[ijk+kk1]-vmean[k+1],2) + ci2*std::pow(v[ijk+kk2]-vmean[k+2],2) + ci3*std::pow(v[ijk+kk3]-vmean[k+3],2))*wy[ijk+kk2]) ) * dzi4[k];

            tke_turb[k] -= 0.5*( cg0*std::pow(w[ijk-kk1], 3) + cg1*std::pow(w[ijk], 3) + cg2*std::pow(w[ijk+kk1], 3) + cg3*std::pow(w[ijk+kk2], 3)) * dzi4[k];
        }
    tke_turb[k] += 0.5*(u2_turb[k] + v2_turb[k]);

    // interior
    for (int k=grid.kstart+1; k<grid.kend-1; ++k)
    {
        u2_turb [k] = 0.;
        v2_turb [k] = 0.;
        tke_turb[k] = 0.;

        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                u2_turb[k]  -= ( cg0*((ci0*std::pow(u[ijk-kk3]-umean[k-3],2) + ci1*std::pow(u[ijk-kk2]-umean[k-2],2) + ci2*std::pow(u[ijk-kk1]-umean[k-1],2) + ci3*std::pow(u[ijk    ]-umean[k  ],2))*wx[ijk-kk1])
                               + cg1*((ci0*std::pow(u[ijk-kk2]-umean[k-2],2) + ci1*std::pow(u[ijk-kk1]-umean[k-1],2) + ci2*std::pow(u[ijk    ]-umean[k  ],2) + ci3*std::pow(u[ijk+kk1]-umean[k+1],2))*wx[ijk    ])
                               + cg2*((ci0*std::pow(u[ijk-kk1]-umean[k-1],2) + ci1*std::pow(u[ijk    ]-umean[k  ],2) + ci2*std::pow(u[ijk+kk1]-umean[k+1],2) + ci3*std::pow(u[ijk+kk2]-umean[k+2],2))*wx[ijk+kk1])
                               + cg3*((ci0*std::pow(u[ijk    ]-umean[k  ],2) + ci1*std::pow(u[ijk+kk1]-umean[k+1],2) + ci2*std::pow(u[ijk+kk2]-umean[k+2],2) + ci3*std::pow(u[ijk+kk3]-umean[k+3],2))*wx[ijk+kk2]) ) * dzi4[k];

                v2_turb[k]  -= ( cg0*((ci0*std::pow(v[ijk-kk3]-vmean[k-3],2) + ci1*std::pow(v[ijk-kk2]-vmean[k-2],2) + ci2*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci3*std::pow(v[ijk    ]-vmean[k  ],2))*wy[ijk-kk1])
                               + cg1*((ci0*std::pow(v[ijk-kk2]-vmean[k-2],2) + ci1*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci2*std::pow(v[ijk    ]-vmean[k  ],2) + ci3*std::pow(v[ijk+kk1]-vmean[k+1],2))*wy[ijk    ])
                               + cg2*((ci0*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci1*std::pow(v[ijk    ]-vmean[k  ],2) + ci2*std::pow(v[ijk+kk1]-vmean[k+1],2) + ci3*std::pow(v[ijk+kk2]-vmean[k+2],2))*wy[ijk+kk1])
                               + cg3*((ci0*std::pow(v[ijk    ]-vmean[k  ],2) + ci1*std::pow(v[ijk+kk1]-vmean[k+1],2) + ci2*std::pow(v[ijk+kk2]-vmean[k+2],2) + ci3*std::pow(v[ijk+kk3]-vmean[k+3],2))*wy[ijk+kk2]) ) * dzi4[k];

                tke_turb[k] -= 0.5*( cg0*std::pow(w[ijk-kk1], 3) + cg1*std::pow(w[ijk], 3) + cg2*std::pow(w[ijk+kk1], 3) + cg3*std::pow(w[ijk+kk2], 3)) * dzi4[k];
            }
        tke_turb[k] += 0.5*(u2_turb[k] + v2_turb[k]);
    }

    // top boundary
    k = grid.kend-1;

    u2_turb [k] = 0.;
    v2_turb [k] = 0.;
    tke_turb[k] = 0.;

    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            u2_turb[k]  -= ( cg0*((ci0*std::pow(u[ijk-kk3]-umean[k-3],2) + ci1*std::pow(u[ijk-kk2]-umean[k-2],2) + ci2*std::pow(u[ijk-kk1]-umean[k-1],2) + ci3*std::pow(u[ijk    ]-umean[k  ],2))*wx[ijk-kk1])
                           + cg1*((ci0*std::pow(u[ijk-kk2]-umean[k-2],2) + ci1*std::pow(u[ijk-kk1]-umean[k-1],2) + ci2*std::pow(u[ijk    ]-umean[k  ],2) + ci3*std::pow(u[ijk+kk1]-umean[k+1],2))*wx[ijk    ])
                           + cg2*((ci0*std::pow(u[ijk-kk1]-umean[k-1],2) + ci1*std::pow(u[ijk    ]-umean[k  ],2) + ci2*std::pow(u[ijk+kk1]-umean[k+1],2) + ci3*std::pow(u[ijk+kk2]-umean[k+2],2))*wx[ijk+kk1])
                           + cg3*((ti0*std::pow(u[ijk-kk1]-umean[k-1],2) + ti1*std::pow(u[ijk    ]-umean[k  ],2) + ti2*std::pow(u[ijk+kk1]-umean[k+1],2) + ti3*std::pow(u[ijk+kk2]-umean[k+2],2))*wx[ijk+kk1]) ) * dzi4[k];

            v2_turb[k]  -= ( cg0*((ci0*std::pow(v[ijk-kk3]-vmean[k-3],2) + ci1*std::pow(v[ijk-kk2]-vmean[k-2],2) + ci2*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci3*std::pow(v[ijk    ]-vmean[k  ],2))*wy[ijk-kk1])
                           + cg1*((ci0*std::pow(v[ijk-kk2]-vmean[k-2],2) + ci1*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci2*std::pow(v[ijk    ]-vmean[k  ],2) + ci3*std::pow(v[ijk+kk1]-vmean[k+1],2))*wy[ijk    ])
                           + cg2*((ci0*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci1*std::pow(v[ijk    ]-vmean[k  ],2) + ci2*std::pow(v[ijk+kk1]-vmean[k+1],2) + ci3*std::pow(v[ijk+kk2]-vmean[k+2],2))*wy[ijk+kk1])
                           + cg3*((ti0*std::pow(v[ijk-kk1]-vmean[k-1],2) + ti1*std::pow(v[ijk    ]-vmean[k  ],2) + ti2*std::pow(v[ijk+kk1]-vmean[k+1],2) + ti3*std::pow(v[ijk+kk2]-vmean[k+2],2))*wy[ijk+kk1]) ) * dzi4[k];

            tke_turb[k] -= 0.5*( cg0*std::pow(w[ijk-kk1], 3) + cg1*std::pow(w[ijk], 3) + cg2*std::pow(w[ijk+kk1], 3) + cg3*std::pow(w[ijk+kk2], 3)) * dzi4[k];
        }
    tke_turb[k] += 0.5*(u2_turb[k] + v2_turb[k]);

    // calculate the vertical velocity term and the vertical reynold stresses
    k = grid.kstart;
    w2_turb[k] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            w2_turb[k] -= ( cg0*(bi0*std::pow(w[ijk-kk2],3) + bi1*std::pow(w[ijk-kk1],3) + bi2*std::pow(w[ijk    ],3) + bi3*std::pow(w[ijk+kk1],3))
                          + cg1*(ci0*std::pow(w[ijk-kk2],3) + ci1*std::pow(w[ijk-kk1],3) + ci2*std::pow(w[ijk    ],3) + ci3*std::pow(w[ijk+kk1],3))
                          + cg2*(ci0*std::pow(w[ijk-kk1],3) + ci1*std::pow(w[ijk    ],3) + ci2*std::pow(w[ijk+kk1],3) + ci3*std::pow(w[ijk+kk2],3))
                          + cg3*(ci0*std::pow(w[ijk    ],3) + ci1*std::pow(w[ijk+kk1],3) + ci2*std::pow(w[ijk+kk2],3) + ci3*std::pow(w[ijk+kk3],3)) ) * dzhi4[k];

            uw_turb[k] -= ( ( cg0*( std::pow(bi0*wx[ijk-kk2] + bi1*wx[ijk-kk1] + bi2*wx[ijk    ] + bi3*wx[ijk+kk1], 2) * (u[ijk-kk2]-umean[k-2]) )
                            + cg1*( std::pow(ci0*wx[ijk-kk2] + ci1*wx[ijk-kk1] + ci2*wx[ijk    ] + ci3*wx[ijk+kk1], 2) * (u[ijk-kk1]-umean[k-1]) )
                            + cg2*( std::pow(ci0*wx[ijk-kk1] + ci1*wx[ijk    ] + ci2*wx[ijk+kk1] + ci3*wx[ijk+kk2], 2) * (u[ijk    ]-umean[k  ]) )
                            + cg3*( std::pow(ci0*wx[ijk    ] + ci1*wx[ijk+kk1] + ci2*wx[ijk+kk2] + ci3*wx[ijk+kk3], 2) * (u[ijk+kk1]-umean[k+1]) ) )
                          * dzhi4[k  ] );
        }

    for (int k=grid.kstart+1; k<grid.kend; ++k)
    {
        w2_turb[k] = 0.;
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                w2_turb[k] -= ( cg0*(ci0*std::pow(w[ijk-kk3],3) + ci1*std::pow(w[ijk-kk2],3) + ci2*std::pow(w[ijk-kk1],3) + ci3*std::pow(w[ijk    ],3))
                              + cg1*(ci0*std::pow(w[ijk-kk2],3) + ci1*std::pow(w[ijk-kk1],3) + ci2*std::pow(w[ijk    ],3) + ci3*std::pow(w[ijk+kk1],3))
                              + cg2*(ci0*std::pow(w[ijk-kk1],3) + ci1*std::pow(w[ijk    ],3) + ci2*std::pow(w[ijk+kk1],3) + ci3*std::pow(w[ijk+kk2],3))
                              + cg3*(ci0*std::pow(w[ijk    ],3) + ci1*std::pow(w[ijk+kk1],3) + ci2*std::pow(w[ijk+kk2],3) + ci3*std::pow(w[ijk+kk3],3)) ) * dzhi4[k];

                uw_turb[k] -= ( ( cg0*( std::pow(ci0*wx[ijk-kk3] + ci1*wx[ijk-kk2] + ci2*wx[ijk-kk1] + ci3*wx[ijk    ], 2) * (u[ijk-kk2]-umean[k-2]) )
                                + cg1*( std::pow(ci0*wx[ijk-kk2] + ci1*wx[ijk-kk1] + ci2*wx[ijk    ] + ci3*wx[ijk+kk1], 2) * (u[ijk-kk1]-umean[k-1]) )
                                + cg2*( std::pow(ci0*wx[ijk-kk1] + ci1*wx[ijk    ] + ci2*wx[ijk+kk1] + ci3*wx[ijk+kk2], 2) * (u[ijk    ]-umean[k  ]) )
                                + cg3*( std::pow(ci0*wx[ijk    ] + ci1*wx[ijk+kk1] + ci2*wx[ijk+kk2] + ci3*wx[ijk+kk3], 2) * (u[ijk+kk1]-umean[k+1]) ) )
                              * dzhi4[k  ] );
            }
    }

    k = grid.kend;
    w2_turb[k] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            w2_turb[k] -= ( cg0*(ci0*std::pow(w[ijk-kk3],3) + ci1*std::pow(w[ijk-kk2],3) + ci2*std::pow(w[ijk-kk1],3) + ci3*std::pow(w[ijk    ],3))
                          + cg1*(ci0*std::pow(w[ijk-kk2],3) + ci1*std::pow(w[ijk-kk1],3) + ci2*std::pow(w[ijk    ],3) + ci3*std::pow(w[ijk+kk1],3))
                          + cg2*(ci0*std::pow(w[ijk-kk1],3) + ci1*std::pow(w[ijk    ],3) + ci2*std::pow(w[ijk+kk1],3) + ci3*std::pow(w[ijk+kk2],3))
                          + cg3*(ti0*std::pow(w[ijk-kk1],3) + ti1*std::pow(w[ijk    ],3) + ti2*std::pow(w[ijk+kk1],3) + ti3*std::pow(w[ijk+kk2],3)) ) * dzhi4[k];

            uw_turb[k] -= ( ( cg0*( ( ci0*wx[ijk-kk3] + ci1*wx[ijk-kk2] + ci2*wx[ijk-kk1] + ci3*wx[ijk    ] ) * ( u[ijk-kk2] - umean[k-2] ) )
                            + cg1*( ( ci0*wx[ijk-kk2] + ci1*wx[ijk-kk1] + ci2*wx[ijk    ] + ci3*wx[ijk+kk1] ) * ( u[ijk-kk1] - umean[k-1] ) )
                            + cg2*( ( ci0*wx[ijk-kk1] + ci1*wx[ijk    ] + ci2*wx[ijk+kk1] + ci3*wx[ijk+kk2] ) * ( u[ijk    ] - umean[k  ] ) )
                            + cg3*( ( ti0*wx[ijk-kk1] + ti1*wx[ijk    ] + ti2*wx[ijk+kk1] + ti3*wx[ijk+kk2] ) * ( u[ijk+kk1] - umean[k+1] ) ) )
                          * dzhi4[k  ] );
        }

    // calculate the profiles
    master.sum(u2_turb , grid.kcells);
    master.sum(v2_turb , grid.kcells);
    master.sum(w2_turb , grid.kcells);
    master.sum(tke_turb, grid.kcells);
    master.sum(uw_turb , grid.kcells);

    for (k=grid.kstart; k<grid.kend; ++k)
    {
        u2_turb [k] /= n;
        v2_turb [k] /= n;
        tke_turb[k] /= n;
    }

    for (k=grid.kstart; k<grid.kend+1; ++k)
    {
        w2_turb [k] /= n;
        uw_turb [k] /= n;
    }
}

void Budget_4::calc_tke_budget(double* restrict u, double* restrict v, double* restrict w, double* restrict p,
                               double* restrict wz, double* restrict uz,
                               double* restrict umean, double* restrict vmean,
                               double* restrict u2_visc, double* restrict v2_visc, double* restrict w2_visc, double* restrict tke_visc, double* restrict uw_visc,
                               double* restrict u2_diss, double* restrict v2_diss, double* restrict w2_diss, double* restrict tke_diss, double* restrict uw_diss,
                               double* restrict w2_pres, double* restrict tke_pres, double* restrict uw_pres,
                               double* restrict u2_rdstr, double* restrict v2_rdstr, double* restrict w2_rdstr, double* restrict uw_rdstr,
                               double* restrict dzi4, double* restrict dzhi4, double visc)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int ii3 = 3;
    const int jj1 = 1*grid.icells;
    const int jj2 = 2*grid.icells;
    const int jj3 = 3*grid.icells;
    const int kk1 = 1*grid.ijcells;
    const int kk2 = 2*grid.ijcells;
    const int kk3 = 3*grid.ijcells;
    const int kk4 = 4*grid.ijcells;

    const int kstart = grid.kstart;
    const int kend   = grid.kend;

    const double dxi = 1./grid.dx;
    const double dyi = 1./grid.dy;

    const double dzhi4bot = grid.dzhi4bot;
    const double dzhi4top = grid.dzhi4top;

    const double n = grid.itot*grid.jtot;

    // 4. CALCULATE THE PRESSURE TRANSPORT TERM
    // bottom boundary
    int k = grid.kstart;
    tke_pres[k] = 0.;

    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            tke_pres[k] -= ( cg0*((bi0*p[ijk-kk2] + bi1*p[ijk-kk1] + bi2*p[ijk    ] + bi3*p[ijk+kk1])*w[ijk-kk1])
                           + cg1*((ci0*p[ijk-kk2] + ci1*p[ijk-kk1] + ci2*p[ijk    ] + ci3*p[ijk+kk1])*w[ijk    ])
                           + cg2*((ci0*p[ijk-kk1] + ci1*p[ijk    ] + ci2*p[ijk+kk1] + ci3*p[ijk+kk2])*w[ijk+kk1])
                           + cg3*((ci0*p[ijk    ] + ci1*p[ijk+kk1] + ci2*p[ijk+kk2] + ci3*p[ijk+kk3])*w[ijk+kk2]) ) * dzi4[k];
        }

    // interior
    for (int k=grid.kstart+1; k<grid.kend-1; ++k)
    {
        tke_pres[k] = 0.;

        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                tke_pres[k] -= ( cg0*((ci0*p[ijk-kk3] + ci1*p[ijk-kk2] + ci2*p[ijk-kk1] + ci3*p[ijk    ])*w[ijk-kk1])
                               + cg1*((ci0*p[ijk-kk2] + ci1*p[ijk-kk1] + ci2*p[ijk    ] + ci3*p[ijk+kk1])*w[ijk    ])
                               + cg2*((ci0*p[ijk-kk1] + ci1*p[ijk    ] + ci2*p[ijk+kk1] + ci3*p[ijk+kk2])*w[ijk+kk1])
                               + cg3*((ci0*p[ijk    ] + ci1*p[ijk+kk1] + ci2*p[ijk+kk2] + ci3*p[ijk+kk3])*w[ijk+kk2]) ) * dzi4[k];
            }
    }

    // top boundary
    k = grid.kend-1;
    tke_pres[k] = 0.;

    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            tke_pres[k] -= ( cg0*((ci0*p[ijk-kk3] + ci1*p[ijk-kk2] + ci2*p[ijk-kk1] + ci3*p[ijk    ])*w[ijk-kk1])
                           + cg1*((ci0*p[ijk-kk2] + ci1*p[ijk-kk1] + ci2*p[ijk    ] + ci3*p[ijk+kk1])*w[ijk    ])
                           + cg2*((ci0*p[ijk-kk1] + ci1*p[ijk    ] + ci2*p[ijk+kk1] + ci3*p[ijk+kk2])*w[ijk+kk1])
                           + cg3*((ti0*p[ijk-kk1] + ti1*p[ijk    ] + ti2*p[ijk+kk1] + ti3*p[ijk+kk2])*w[ijk+kk2]) ) * dzi4[k];
        }

    // calculate the vertical velocity pressure transport term
    // \TODO implement the proper BC as soon as the full BC's for pressure are added
    // bottom boundary
    k = grid.kstart;
    w2_pres[k] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            w2_pres[k] -= 0.*( cg0*((bi0*w[ijk-kk2] + bi1*w[ijk-kk1] + bi2*w[ijk    ] + bi3*w[ijk+kk1])*p[ijk-kk2])
                             + cg1*((ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1])*p[ijk-kk1])
                             + cg2*((ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])*p[ijk    ])
                             + cg3*((ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3])*p[ijk+kk1]) ) * dzhi4[k];
        }

    // interior
    for (int k=grid.kstart+1; k<grid.kend; ++k)
    {
        w2_pres[k] = 0.;
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                w2_pres[k] -= 2.*( cg0*((ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ])*p[ijk-kk2])
                                 + cg1*((ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1])*p[ijk-kk1])
                                 + cg2*((ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])*p[ijk    ])
                                 + cg3*((ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3])*p[ijk+kk1]) ) * dzhi4[k];
            }
    }

    // top boundary
    k = grid.kend;
    w2_pres[k] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            w2_pres[k] -= 0.*( cg0*((ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ])*p[ijk-kk2])
                             + cg1*((ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1])*p[ijk-kk1])
                             + cg2*((ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])*p[ijk    ])
                             + cg3*((ti0*w[ijk-kk1] + ti1*w[ijk    ] + ti2*w[ijk+kk1] + ti3*w[ijk+kk2])*p[ijk+kk1]) ) * dzhi4[k];
        }

    // UW term
    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        uw_pres[k] = 0.;
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                uw_pres[k] -= ( ( ( cg0*( ( u[ijk        -kk2] - umean[k-2] ) * ( ci0*p[ijk-ii2    -kk2] + ci1*p[ijk-ii1    -kk2] + ci2*p[ijk        -kk2] + ci3*p[ijk+ii1    -kk2] ) )
                                  + cg1*( ( u[ijk        -kk1] - umean[k-1] ) * ( ci0*p[ijk-ii2    -kk1] + ci1*p[ijk-ii1    -kk1] + ci2*p[ijk        -kk1] + ci3*p[ijk+ii1    -kk1] ) )
                                  + cg2*( ( u[ijk            ] - umean[k  ] ) * ( ci0*p[ijk-ii2        ] + ci1*p[ijk-ii1        ] + ci2*p[ijk            ] + ci3*p[ijk+ii1        ] ) )
                                  + cg3*( ( u[ijk        +kk1] - umean[k+1] ) * ( ci0*p[ijk-ii2    +kk1] + ci1*p[ijk-ii1    +kk1] + ci2*p[ijk        +kk1] + ci3*p[ijk+ii1    +kk1] ) ) )

                                * dzhi4[k  ] )

                              + ( ( cg0*( w[ijk-ii2        ] * ( ci0*p[ijk-ii2    -kk2] + ci1*p[ijk-ii2    -kk1] + ci2*p[ijk-ii2        ] + ci3*p[ijk-ii2    +kk1] ) )
                                  + cg1*( w[ijk-ii1        ] * ( ci0*p[ijk-ii1    -kk2] + ci1*p[ijk-ii1    -kk1] + ci2*p[ijk-ii1        ] + ci3*p[ijk-ii1    +kk1] ) )
                                  + cg2*( w[ijk            ] * ( ci0*p[ijk        -kk2] + ci1*p[ijk        -kk1] + ci2*p[ijk            ] + ci3*p[ijk        +kk1] ) )
                                  + cg3*( w[ijk+ii1        ] * ( ci0*p[ijk+ii1    -kk2] + ci1*p[ijk+ii1    -kk1] + ci2*p[ijk+ii1        ] + ci3*p[ijk+ii1    +kk1] ) ) )

                                * cgi*dxi ) );
            }
    }

    master.sum(w2_pres , grid.kcells);
    master.sum(tke_pres, grid.kcells);
    master.sum(uw_pres , grid.kcells);

    for (int k=grid.kstart; k<grid.kend; ++k)
        tke_pres[k] /= n;

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        w2_pres[k] /= n;
        uw_pres[k] /= n;
    }

    // 5. CALCULATE THE VISCOUS TRANSPORT TERM
    // first, interpolate the vertical velocity to the scalar levels using temporary array wz
    for (int k=grid.kstart; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                wz[ijk] = ci0*w[ijk-kk1] + ci1*w[ijk] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2];
            }

    // calculate the ghost cells at the bottom
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + kstart*kk1;
            wz[ijk-kk1] = - 2.*wz[ijk] + (1./3.)*wz[ijk+kk1];
            wz[ijk-kk2] = - 9.*wz[ijk] + 2.*wz[ijk+kk1];
        }

    // calculate the ghost cells at the top
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + (kend-1)*kk1;
            wz[ijk+kk1] = - 2.*wz[ijk] + (1./3.)*wz[ijk-kk1];
            wz[ijk+kk2] = - 9.*wz[ijk] + 2.*wz[ijk-kk1];
        }

    // first, interpolate the horizontal velocity to the flux levels using temporary array uz
    k = kstart-1;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            uz[ijk] = bi0*u[ijk-kk1] + bi1*u[ijk] + bi2*u[ijk+kk1] + bi3*u[ijk+kk2];
        }

    for (int k=grid.kstart; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                uz[ijk] = ci0*u[ijk-kk2] + ci1*u[ijk-kk1] + ci2*u[ijk] + ci3*u[ijk+kk1];
            }

    k = kend;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            uz[ijk] = ti0*u[ijk-kk2] + ti1*u[ijk-kk1] + ti2*u[ijk] + ti3*u[ijk+kk1];
        }


    // bottom boundary
    k = grid.kstart;

    u2_visc [k] = 0.;
    v2_visc [k] = 0.;
    tke_visc[k] = 0.;

    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk  = i + j*jj1 + k*kk1;
            u2_visc[k]  += visc * ( cg0*((bg0*std::pow(u[ijk-kk2]-umean[k-2],2) + bg1*std::pow(u[ijk-kk1]-umean[k-1],2) + bg2*std::pow(u[ijk    ]-umean[k  ],2) + bg3*std::pow(u[ijk+kk1]-umean[k+1],2)) * dzhi4[k-1])
                                  + cg1*((cg0*std::pow(u[ijk-kk2]-umean[k-2],2) + cg1*std::pow(u[ijk-kk1]-umean[k-1],2) + cg2*std::pow(u[ijk    ]-umean[k  ],2) + cg3*std::pow(u[ijk+kk1]-umean[k+1],2)) * dzhi4[k  ])
                                  + cg2*((cg0*std::pow(u[ijk-kk1]-umean[k-1],2) + cg1*std::pow(u[ijk    ]-umean[k  ],2) + cg2*std::pow(u[ijk+kk1]-umean[k+1],2) + cg3*std::pow(u[ijk+kk2]-umean[k+2],2)) * dzhi4[k+1])
                                  + cg3*((cg0*std::pow(u[ijk    ]-umean[k  ],2) + cg1*std::pow(u[ijk+kk1]-umean[k+1],2) + cg2*std::pow(u[ijk+kk2]-umean[k+2],2) + cg3*std::pow(u[ijk+kk3]-umean[k+3],2)) * dzhi4[k+2]) ) * dzi4[k];

            v2_visc[k]  += visc * ( cg0*((bg0*std::pow(v[ijk-kk2]-vmean[k-2],2) + bg1*std::pow(v[ijk-kk1]-vmean[k-1],2) + bg2*std::pow(v[ijk    ]-vmean[k  ],2) + bg3*std::pow(v[ijk+kk1]-vmean[k+1],2)) * dzhi4[k-1])
                                  + cg1*((cg0*std::pow(v[ijk-kk2]-vmean[k-2],2) + cg1*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg2*std::pow(v[ijk    ]-vmean[k  ],2) + cg3*std::pow(v[ijk+kk1]-vmean[k+1],2)) * dzhi4[k  ])
                                  + cg2*((cg0*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg1*std::pow(v[ijk    ]-vmean[k  ],2) + cg2*std::pow(v[ijk+kk1]-vmean[k+1],2) + cg3*std::pow(v[ijk+kk2]-vmean[k+2],2)) * dzhi4[k+1])
                                  + cg3*((cg0*std::pow(v[ijk    ]-vmean[k  ],2) + cg1*std::pow(v[ijk+kk1]-vmean[k+1],2) + cg2*std::pow(v[ijk+kk2]-vmean[k+2],2) + cg3*std::pow(v[ijk+kk3]-vmean[k+3],2)) * dzhi4[k+2]) ) * dzi4[k];

            tke_visc[k] += 0.5 * visc * ( cg0*((bg0*std::pow(wz[ijk-kk2],2) + bg1*std::pow(wz[ijk-kk1],2) + bg2*std::pow(wz[ijk    ],2) + bg3*std::pow(wz[ijk+kk1],2)) * dzhi4[k-1])
                                        + cg1*((cg0*std::pow(wz[ijk-kk2],2) + cg1*std::pow(wz[ijk-kk1],2) + cg2*std::pow(wz[ijk    ],2) + cg3*std::pow(wz[ijk+kk1],2)) * dzhi4[k  ])
                                        + cg2*((cg0*std::pow(wz[ijk-kk1],2) + cg1*std::pow(wz[ijk    ],2) + cg2*std::pow(wz[ijk+kk1],2) + cg3*std::pow(wz[ijk+kk2],2)) * dzhi4[k+1])
                                        + cg3*((cg0*std::pow(wz[ijk    ],2) + cg1*std::pow(wz[ijk+kk1],2) + cg2*std::pow(wz[ijk+kk2],2) + cg3*std::pow(wz[ijk+kk3],2)) * dzhi4[k+2]) ) * dzi4[k];
        }
    tke_visc[k] += 0.5*(u2_visc[k] + v2_visc[k]);

    // interior
    for (int k=grid.kstart+1; k<grid.kend-1; ++k)
    {
        u2_visc [k] = 0.;
        v2_visc [k] = 0.;
        tke_visc[k] = 0.;

        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                u2_visc[k]  += visc * ( cg0*((cg0*std::pow(u[ijk-kk3]-umean[k-3],2) + cg1*std::pow(u[ijk-kk2]-umean[k-2],2) + cg2*std::pow(u[ijk-kk1]-umean[k-1],2) + cg3*std::pow(u[ijk    ]-umean[k  ],2)) * dzhi4[k-1])
                                      + cg1*((cg0*std::pow(u[ijk-kk2]-umean[k-2],2) + cg1*std::pow(u[ijk-kk1]-umean[k-1],2) + cg2*std::pow(u[ijk    ]-umean[k  ],2) + cg3*std::pow(u[ijk+kk1]-umean[k+1],2)) * dzhi4[k  ])
                                      + cg2*((cg0*std::pow(u[ijk-kk1]-umean[k-1],2) + cg1*std::pow(u[ijk    ]-umean[k  ],2) + cg2*std::pow(u[ijk+kk1]-umean[k+1],2) + cg3*std::pow(u[ijk+kk2]-umean[k+2],2)) * dzhi4[k+1])
                                      + cg3*((cg0*std::pow(u[ijk    ]-umean[k  ],2) + cg1*std::pow(u[ijk+kk1]-umean[k+1],2) + cg2*std::pow(u[ijk+kk2]-umean[k+2],2) + cg3*std::pow(u[ijk+kk3]-umean[k+3],2)) * dzhi4[k+2]) ) * dzi4[k];

                v2_visc[k]  += visc * ( cg0*((cg0*std::pow(v[ijk-kk3]-vmean[k-3],2) + cg1*std::pow(v[ijk-kk2]-vmean[k-2],2) + cg2*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg3*std::pow(v[ijk    ]-vmean[k  ],2)) * dzhi4[k-1])
                                      + cg1*((cg0*std::pow(v[ijk-kk2]-vmean[k-2],2) + cg1*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg2*std::pow(v[ijk    ]-vmean[k  ],2) + cg3*std::pow(v[ijk+kk1]-vmean[k+1],2)) * dzhi4[k  ])
                                      + cg2*((cg0*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg1*std::pow(v[ijk    ]-vmean[k  ],2) + cg2*std::pow(v[ijk+kk1]-vmean[k+1],2) + cg3*std::pow(v[ijk+kk2]-vmean[k+2],2)) * dzhi4[k+1])
                                      + cg3*((cg0*std::pow(v[ijk    ]-vmean[k  ],2) + cg1*std::pow(v[ijk+kk1]-vmean[k+1],2) + cg2*std::pow(v[ijk+kk2]-vmean[k+2],2) + cg3*std::pow(v[ijk+kk3]-vmean[k+3],2)) * dzhi4[k+2]) ) * dzi4[k];

                tke_visc[k] += 0.5 * visc * ( cg0*((cg0*std::pow(wz[ijk-kk3],2) + cg1*std::pow(wz[ijk-kk2],2) + cg2*std::pow(wz[ijk-kk1],2) + cg3*std::pow(wz[ijk    ],2)) * dzhi4[k-1])
                                            + cg1*((cg0*std::pow(wz[ijk-kk2],2) + cg1*std::pow(wz[ijk-kk1],2) + cg2*std::pow(wz[ijk    ],2) + cg3*std::pow(wz[ijk+kk1],2)) * dzhi4[k  ])
                                            + cg2*((cg0*std::pow(wz[ijk-kk1],2) + cg1*std::pow(wz[ijk    ],2) + cg2*std::pow(wz[ijk+kk1],2) + cg3*std::pow(wz[ijk+kk2],2)) * dzhi4[k+1])
                                            + cg3*((cg0*std::pow(wz[ijk    ],2) + cg1*std::pow(wz[ijk+kk1],2) + cg2*std::pow(wz[ijk+kk2],2) + cg3*std::pow(wz[ijk+kk3],2)) * dzhi4[k+2]) ) * dzi4[k];
            }
        tke_visc[k] += 0.5*(u2_visc[k] + v2_visc[k]);
    }

    // top boundary
    k = grid.kend-1;
    u2_visc [k] = 0.;
    v2_visc [k] = 0.;
    tke_visc[k] = 0.;

    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            u2_visc[k]  += visc * ( cg0*((cg0*std::pow(u[ijk-kk3]-umean[k-3],2) + cg1*std::pow(u[ijk-kk2]-umean[k-2],2) + cg2*std::pow(u[ijk-kk1]-umean[k-1],2) + cg3*std::pow(u[ijk    ]-umean[k  ],2)) * dzhi4[k-1])
                                  + cg1*((cg0*std::pow(u[ijk-kk2]-umean[k-2],2) + cg1*std::pow(u[ijk-kk1]-umean[k-1],2) + cg2*std::pow(u[ijk    ]-umean[k  ],2) + cg3*std::pow(u[ijk+kk1]-umean[k+1],2)) * dzhi4[k  ])
                                  + cg2*((cg0*std::pow(u[ijk-kk1]-umean[k-1],2) + cg1*std::pow(u[ijk    ]-umean[k  ],2) + cg2*std::pow(u[ijk+kk1]-umean[k+1],2) + cg3*std::pow(u[ijk+kk2]-umean[k+2],2)) * dzhi4[k+1])
                                  + cg3*((tg0*std::pow(u[ijk-kk1]-umean[k-1],2) + tg1*std::pow(u[ijk    ]-umean[k  ],2) + tg2*std::pow(u[ijk+kk1]-umean[k+1],2) + tg3*std::pow(u[ijk+kk2]-umean[k+2],2)) * dzhi4[k+2]) ) * dzi4[k];

            v2_visc[k]  += visc * ( cg0*((cg0*std::pow(v[ijk-kk3]-vmean[k-3],2) + cg1*std::pow(v[ijk-kk2]-vmean[k-2],2) + cg2*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg3*std::pow(v[ijk    ]-vmean[k  ],2)) * dzhi4[k-1])
                                  + cg1*((cg0*std::pow(v[ijk-kk2]-vmean[k-2],2) + cg1*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg2*std::pow(v[ijk    ]-vmean[k  ],2) + cg3*std::pow(v[ijk+kk1]-vmean[k+1],2)) * dzhi4[k  ])
                                  + cg2*((cg0*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg1*std::pow(v[ijk    ]-vmean[k  ],2) + cg2*std::pow(v[ijk+kk1]-vmean[k+1],2) + cg3*std::pow(v[ijk+kk2]-vmean[k+2],2)) * dzhi4[k+1])
                                  + cg3*((tg0*std::pow(v[ijk-kk1]-vmean[k-1],2) + tg1*std::pow(v[ijk    ]-vmean[k  ],2) + tg2*std::pow(v[ijk+kk1]-vmean[k+1],2) + tg3*std::pow(v[ijk+kk2]-vmean[k+2],2)) * dzhi4[k+2]) ) * dzi4[k];

            tke_visc[k] += 0.5 * visc * ( cg0*((cg0*std::pow(wz[ijk-kk3],2) + cg1*std::pow(wz[ijk-kk2],2) + cg2*std::pow(wz[ijk-kk1],2) + cg3*std::pow(wz[ijk    ],2)) * dzhi4[k-1])
                                        + cg1*((cg0*std::pow(wz[ijk-kk2],2) + cg1*std::pow(wz[ijk-kk1],2) + cg2*std::pow(wz[ijk    ],2) + cg3*std::pow(wz[ijk+kk1],2)) * dzhi4[k  ])
                                        + cg2*((cg0*std::pow(wz[ijk-kk1],2) + cg1*std::pow(wz[ijk    ],2) + cg2*std::pow(wz[ijk+kk1],2) + cg3*std::pow(wz[ijk+kk2],2)) * dzhi4[k+1])
                                        + cg3*((tg0*std::pow(wz[ijk-kk1],2) + tg1*std::pow(wz[ijk    ],2) + tg2*std::pow(wz[ijk+kk1],2) + tg3*std::pow(wz[ijk+kk2],2)) * dzhi4[k+2]) ) * dzi4[k];
        }
    tke_visc[k] += 0.5*(u2_visc[k] + v2_visc[k]);

    // Calculate the viscous transport of vertical velocity variance and the fluxes.
    // Interpolate

    // Bottom boundary.
    k = grid.kstart;
    w2_visc[k] = 0.;
    uw_visc[k] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            w2_visc[k] += visc * ( bg0*((bg0*std::pow(w[ijk-kk1],2) + bg1*std::pow(w[ijk    ],2) + bg2*std::pow(w[ijk+kk1],2) + bg3*std::pow(w[ijk+kk2],2)) * dzi4[k-1])
                                 + bg1*((cg0*std::pow(w[ijk-kk1],2) + cg1*std::pow(w[ijk    ],2) + cg2*std::pow(w[ijk+kk1],2) + cg3*std::pow(w[ijk+kk2],2)) * dzi4[k  ])
                                 + bg2*((cg0*std::pow(w[ijk    ],2) + cg1*std::pow(w[ijk+kk1],2) + cg2*std::pow(w[ijk+kk2],2) + cg3*std::pow(w[ijk+kk3],2)) * dzi4[k+1])
                                 + bg3*((cg0*std::pow(w[ijk+kk1],2) + cg1*std::pow(w[ijk+kk2],2) + cg2*std::pow(w[ijk+kk3],2) + cg3*std::pow(w[ijk+kk4],2)) * dzi4[k+2]) ) * dzhi4bot;


            uw_visc[k] += ( ( visc
                            * ( bg0*( ( bg0*( uz[ijk        -kk1] * ( ci0*w[ijk-ii2    -kk1] + ci1*w[ijk-ii1    -kk1] + ci2*w[ijk        -kk1] + ci3*w[ijk+ii1    -kk1] ) )
                                      + bg1*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) )
                                      + bg2*( uz[ijk        +kk1] * ( ci0*w[ijk-ii2    +kk1] + ci1*w[ijk-ii1    +kk1] + ci2*w[ijk        +kk1] + ci3*w[ijk+ii1    +kk1] ) )
                                      + bg3*( uz[ijk        +kk2] * ( ci0*w[ijk-ii2    +kk2] + ci1*w[ijk-ii1    +kk2] + ci2*w[ijk        +kk2] + ci3*w[ijk+ii1    +kk2] ) ) )

                                    * dzi4[k-1] )

                              + bg1*( ( cg0*( uz[ijk        -kk1] * ( ci0*w[ijk-ii2    -kk1] + ci1*w[ijk-ii1    -kk1] + ci2*w[ijk        -kk1] + ci3*w[ijk+ii1    -kk1] ) )
                                      + cg1*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) )
                                      + cg2*( uz[ijk        +kk1] * ( ci0*w[ijk-ii2    +kk1] + ci1*w[ijk-ii1    +kk1] + ci2*w[ijk        +kk1] + ci3*w[ijk+ii1    +kk1] ) )
                                      + cg3*( uz[ijk        +kk2] * ( ci0*w[ijk-ii2    +kk2] + ci1*w[ijk-ii1    +kk2] + ci2*w[ijk        +kk2] + ci3*w[ijk+ii1    +kk2] ) ) )

                                    * dzi4[k  ] )

                              + bg2*( ( cg0*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) )
                                      + cg1*( uz[ijk        +kk1] * ( ci0*w[ijk-ii2    +kk1] + ci1*w[ijk-ii1    +kk1] + ci2*w[ijk        +kk1] + ci3*w[ijk+ii1    +kk1] ) )
                                      + cg2*( uz[ijk        +kk2] * ( ci0*w[ijk-ii2    +kk2] + ci1*w[ijk-ii1    +kk2] + ci2*w[ijk        +kk2] + ci3*w[ijk+ii1    +kk2] ) )
                                      + cg3*( uz[ijk        +kk3] * ( ci0*w[ijk-ii2    +kk3] + ci1*w[ijk-ii1    +kk3] + ci2*w[ijk        +kk3] + ci3*w[ijk+ii1    +kk3] ) ) )

                                    * dzi4[k+1] )

                              + bg3*( ( cg0*( uz[ijk        +kk1] * ( ci0*w[ijk-ii2    +kk1] + ci1*w[ijk-ii1    +kk1] + ci2*w[ijk        +kk1] + ci3*w[ijk+ii1    +kk1] ) )
                                      + cg1*( uz[ijk        +kk2] * ( ci0*w[ijk-ii2    +kk2] + ci1*w[ijk-ii1    +kk2] + ci2*w[ijk        +kk2] + ci3*w[ijk+ii1    +kk2] ) )
                                      + cg2*( uz[ijk        +kk3] * ( ci0*w[ijk-ii2    +kk3] + ci1*w[ijk-ii1    +kk3] + ci2*w[ijk        +kk3] + ci3*w[ijk+ii1    +kk3] ) )
                                      + cg3*( uz[ijk        +kk4] * ( ci0*w[ijk-ii2    +kk4] + ci1*w[ijk-ii1    +kk4] + ci2*w[ijk        +kk4] + ci3*w[ijk+ii1    +kk4] ) ) )

                                    * dzi4[k+2] ) ) )


                          * dzhi4bot );

        }

    // Bottom boundary + 1.
    k = grid.kstart+1;
    w2_visc[k] = 0.;
    uw_visc[k] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            w2_visc[k] += visc * ( cg0*((bg0*std::pow(w[ijk-kk2],2) + bg1*std::pow(w[ijk-kk1],2) + bg2*std::pow(w[ijk    ],2) + bg3*std::pow(w[ijk+kk1],2)) * dzi4[k-2])
                                 + cg1*((cg0*std::pow(w[ijk-kk2],2) + cg1*std::pow(w[ijk-kk1],2) + cg2*std::pow(w[ijk    ],2) + cg3*std::pow(w[ijk+kk1],2)) * dzi4[k-1])
                                 + cg2*((cg0*std::pow(w[ijk-kk1],2) + cg1*std::pow(w[ijk    ],2) + cg2*std::pow(w[ijk+kk1],2) + cg3*std::pow(w[ijk+kk2],2)) * dzi4[k  ])
                                 + cg3*((cg0*std::pow(w[ijk    ],2) + cg1*std::pow(w[ijk+kk1],2) + cg2*std::pow(w[ijk+kk2],2) + cg3*std::pow(w[ijk+kk3],2)) * dzi4[k+1]) ) * dzhi4[k];


            uw_visc[k] += ( ( visc


                            * ( cg0*( ( bg0*( uz[ijk        -kk2] * ( ci0*w[ijk-ii2    -kk2] + ci1*w[ijk-ii1    -kk2] + ci2*w[ijk        -kk2] + ci3*w[ijk+ii1    -kk2] ) )
                                      + bg1*( uz[ijk        -kk1] * ( ci0*w[ijk-ii2    -kk1] + ci1*w[ijk-ii1    -kk1] + ci2*w[ijk        -kk1] + ci3*w[ijk+ii1    -kk1] ) )
                                      + bg2*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) )
                                      + bg3*( uz[ijk        +kk1] * ( ci0*w[ijk-ii2    +kk1] + ci1*w[ijk-ii1    +kk1] + ci2*w[ijk        +kk1] + ci3*w[ijk+ii1    +kk1] ) ) )

                                    * dzi4[k-2] )

                              + cg1*( ( cg0*( uz[ijk        -kk2] * ( ci0*w[ijk-ii2    -kk2] + ci1*w[ijk-ii1    -kk2] + ci2*w[ijk        -kk2] + ci3*w[ijk+ii1    -kk2] ) )
                                      + cg1*( uz[ijk        -kk1] * ( ci0*w[ijk-ii2    -kk1] + ci1*w[ijk-ii1    -kk1] + ci2*w[ijk        -kk1] + ci3*w[ijk+ii1    -kk1] ) )
                                      + cg2*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) )
                                      + cg3*( uz[ijk        +kk1] * ( ci0*w[ijk-ii2    +kk1] + ci1*w[ijk-ii1    +kk1] + ci2*w[ijk        +kk1] + ci3*w[ijk+ii1    +kk1] ) ) )

                                    * dzi4[k-1] )

                              + cg2*( ( cg0*( uz[ijk        -kk1] * ( ci0*w[ijk-ii2    -kk1] + ci1*w[ijk-ii1    -kk1] + ci2*w[ijk        -kk1] + ci3*w[ijk+ii1    -kk1] ) )
                                      + cg1*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) )
                                      + cg2*( uz[ijk        +kk1] * ( ci0*w[ijk-ii2    +kk1] + ci1*w[ijk-ii1    +kk1] + ci2*w[ijk        +kk1] + ci3*w[ijk+ii1    +kk1] ) )
                                      + cg3*( uz[ijk        +kk2] * ( ci0*w[ijk-ii2    +kk2] + ci1*w[ijk-ii1    +kk2] + ci2*w[ijk        +kk2] + ci3*w[ijk+ii1    +kk2] ) ) )

                                    * dzi4[k  ] )

                              + cg3*( ( cg0*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) )
                                      + cg1*( uz[ijk        +kk1] * ( ci0*w[ijk-ii2    +kk1] + ci1*w[ijk-ii1    +kk1] + ci2*w[ijk        +kk1] + ci3*w[ijk+ii1    +kk1] ) )
                                      + cg2*( uz[ijk        +kk2] * ( ci0*w[ijk-ii2    +kk2] + ci1*w[ijk-ii1    +kk2] + ci2*w[ijk        +kk2] + ci3*w[ijk+ii1    +kk2] ) )
                                      + cg3*( uz[ijk        +kk3] * ( ci0*w[ijk-ii2    +kk3] + ci1*w[ijk-ii1    +kk3] + ci2*w[ijk        +kk3] + ci3*w[ijk+ii1    +kk3] ) ) )

                                    * dzi4[k+1] ) ) )


                          * dzhi4[k  ] );

        }

    // Interior.
    for (int k=grid.kstart+2; k<grid.kend-1; ++k)
    {
        w2_visc[k] = 0.;
        uw_visc[k] = 0.;
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                w2_visc[k] += visc * ( cg0*((cg0*std::pow(w[ijk-kk3],2) + cg1*std::pow(w[ijk-kk2],2) + cg2*std::pow(w[ijk-kk1],2) + cg3*std::pow(w[ijk    ],2)) * dzi4[k-2])
                                     + cg1*((cg0*std::pow(w[ijk-kk2],2) + cg1*std::pow(w[ijk-kk1],2) + cg2*std::pow(w[ijk    ],2) + cg3*std::pow(w[ijk+kk1],2)) * dzi4[k-1])
                                     + cg2*((cg0*std::pow(w[ijk-kk1],2) + cg1*std::pow(w[ijk    ],2) + cg2*std::pow(w[ijk+kk1],2) + cg3*std::pow(w[ijk+kk2],2)) * dzi4[k  ])
                                     + cg3*((cg0*std::pow(w[ijk    ],2) + cg1*std::pow(w[ijk+kk1],2) + cg2*std::pow(w[ijk+kk2],2) + cg3*std::pow(w[ijk+kk3],2)) * dzi4[k+1]) ) * dzhi4[k];


                uw_visc[k] += ( ( visc


                                * ( cg0*( ( cg0*( uz[ijk        -kk3] * ( ci0*w[ijk-ii2    -kk3] + ci1*w[ijk-ii1    -kk3] + ci2*w[ijk        -kk3] + ci3*w[ijk+ii1    -kk3] ) )
                                          + cg1*( uz[ijk        -kk2] * ( ci0*w[ijk-ii2    -kk2] + ci1*w[ijk-ii1    -kk2] + ci2*w[ijk        -kk2] + ci3*w[ijk+ii1    -kk2] ) )
                                          + cg2*( uz[ijk        -kk1] * ( ci0*w[ijk-ii2    -kk1] + ci1*w[ijk-ii1    -kk1] + ci2*w[ijk        -kk1] + ci3*w[ijk+ii1    -kk1] ) )
                                          + cg3*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) ) )

                                        * dzi4[k-2] )

                                  + cg1*( ( cg0*( uz[ijk        -kk2] * ( ci0*w[ijk-ii2    -kk2] + ci1*w[ijk-ii1    -kk2] + ci2*w[ijk        -kk2] + ci3*w[ijk+ii1    -kk2] ) )
                                          + cg1*( uz[ijk        -kk1] * ( ci0*w[ijk-ii2    -kk1] + ci1*w[ijk-ii1    -kk1] + ci2*w[ijk        -kk1] + ci3*w[ijk+ii1    -kk1] ) )
                                          + cg2*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) )
                                          + cg3*( uz[ijk        +kk1] * ( ci0*w[ijk-ii2    +kk1] + ci1*w[ijk-ii1    +kk1] + ci2*w[ijk        +kk1] + ci3*w[ijk+ii1    +kk1] ) ) )

                                        * dzi4[k-1] )

                                  + cg2*( ( cg0*( uz[ijk        -kk1] * ( ci0*w[ijk-ii2    -kk1] + ci1*w[ijk-ii1    -kk1] + ci2*w[ijk        -kk1] + ci3*w[ijk+ii1    -kk1] ) )
                                          + cg1*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) )
                                          + cg2*( uz[ijk        +kk1] * ( ci0*w[ijk-ii2    +kk1] + ci1*w[ijk-ii1    +kk1] + ci2*w[ijk        +kk1] + ci3*w[ijk+ii1    +kk1] ) )
                                          + cg3*( uz[ijk        +kk2] * ( ci0*w[ijk-ii2    +kk2] + ci1*w[ijk-ii1    +kk2] + ci2*w[ijk        +kk2] + ci3*w[ijk+ii1    +kk2] ) ) )

                                        * dzi4[k  ] )

                                  + cg3*( ( cg0*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) )
                                          + cg1*( uz[ijk        +kk1] * ( ci0*w[ijk-ii2    +kk1] + ci1*w[ijk-ii1    +kk1] + ci2*w[ijk        +kk1] + ci3*w[ijk+ii1    +kk1] ) )
                                          + cg2*( uz[ijk        +kk2] * ( ci0*w[ijk-ii2    +kk2] + ci1*w[ijk-ii1    +kk2] + ci2*w[ijk        +kk2] + ci3*w[ijk+ii1    +kk2] ) )
                                          + cg3*( uz[ijk        +kk3] * ( ci0*w[ijk-ii2    +kk3] + ci1*w[ijk-ii1    +kk3] + ci2*w[ijk        +kk3] + ci3*w[ijk+ii1    +kk3] ) ) )

                                        * dzi4[k+1] ) ) )


                              * dzhi4[k  ] );
            }
    }

    // Top boundary - 1.
    k = grid.kend-1;
    w2_visc[k] = 0.;
    uw_visc[k] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            w2_visc[k] += visc * ( cg0*((cg0*std::pow(w[ijk-kk3],2) + cg1*std::pow(w[ijk-kk2],2) + cg2*std::pow(w[ijk-kk1],2) + cg3*std::pow(w[ijk    ],2)) * dzi4[k-2])
                                 + cg1*((cg0*std::pow(w[ijk-kk2],2) + cg1*std::pow(w[ijk-kk1],2) + cg2*std::pow(w[ijk    ],2) + cg3*std::pow(w[ijk+kk1],2)) * dzi4[k-1])
                                 + cg2*((cg0*std::pow(w[ijk-kk1],2) + cg1*std::pow(w[ijk    ],2) + cg2*std::pow(w[ijk+kk1],2) + cg3*std::pow(w[ijk+kk2],2)) * dzi4[k  ])
                                 + cg3*((tg0*std::pow(w[ijk-kk1],2) + tg1*std::pow(w[ijk    ],2) + tg2*std::pow(w[ijk+kk1],2) + tg3*std::pow(w[ijk+kk2],2)) * dzi4[k+1]) ) * dzhi4[k];


            uw_visc[k] += ( ( visc


                            * ( cg0*( ( cg0*( uz[ijk        -kk3] * ( ci0*w[ijk-ii2    -kk3] + ci1*w[ijk-ii1    -kk3] + ci2*w[ijk        -kk3] + ci3*w[ijk+ii1    -kk3] ) )
                                      + cg1*( uz[ijk        -kk2] * ( ci0*w[ijk-ii2    -kk2] + ci1*w[ijk-ii1    -kk2] + ci2*w[ijk        -kk2] + ci3*w[ijk+ii1    -kk2] ) )
                                      + cg2*( uz[ijk        -kk1] * ( ci0*w[ijk-ii2    -kk1] + ci1*w[ijk-ii1    -kk1] + ci2*w[ijk        -kk1] + ci3*w[ijk+ii1    -kk1] ) )
                                      + cg3*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) ) )

                                    * dzi4[k-2] )

                              + cg1*( ( cg0*( uz[ijk        -kk2] * ( ci0*w[ijk-ii2    -kk2] + ci1*w[ijk-ii1    -kk2] + ci2*w[ijk        -kk2] + ci3*w[ijk+ii1    -kk2] ) )
                                      + cg1*( uz[ijk        -kk1] * ( ci0*w[ijk-ii2    -kk1] + ci1*w[ijk-ii1    -kk1] + ci2*w[ijk        -kk1] + ci3*w[ijk+ii1    -kk1] ) )
                                      + cg2*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) )
                                      + cg3*( uz[ijk        +kk1] * ( ci0*w[ijk-ii2    +kk1] + ci1*w[ijk-ii1    +kk1] + ci2*w[ijk        +kk1] + ci3*w[ijk+ii1    +kk1] ) ) )

                                    * dzi4[k-1] )

                              + cg2*( ( cg0*( uz[ijk        -kk1] * ( ci0*w[ijk-ii2    -kk1] + ci1*w[ijk-ii1    -kk1] + ci2*w[ijk        -kk1] + ci3*w[ijk+ii1    -kk1] ) )
                                      + cg1*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) )
                                      + cg2*( uz[ijk        +kk1] * ( ci0*w[ijk-ii2    +kk1] + ci1*w[ijk-ii1    +kk1] + ci2*w[ijk        +kk1] + ci3*w[ijk+ii1    +kk1] ) )
                                      + cg3*( uz[ijk        +kk2] * ( ci0*w[ijk-ii2    +kk2] + ci1*w[ijk-ii1    +kk2] + ci2*w[ijk        +kk2] + ci3*w[ijk+ii1    +kk2] ) ) )

                                    * dzi4[k  ] )

                              + cg3*( ( tg0*( uz[ijk        -kk1] * ( ci0*w[ijk-ii2    -kk1] + ci1*w[ijk-ii1    -kk1] + ci2*w[ijk        -kk1] + ci3*w[ijk+ii1    -kk1] ) )
                                      + tg1*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) )
                                      + tg2*( uz[ijk        +kk1] * ( ci0*w[ijk-ii2    +kk1] + ci1*w[ijk-ii1    +kk1] + ci2*w[ijk        +kk1] + ci3*w[ijk+ii1    +kk1] ) )
                                      + tg3*( uz[ijk        +kk2] * ( ci0*w[ijk-ii2    +kk2] + ci1*w[ijk-ii1    +kk2] + ci2*w[ijk        +kk2] + ci3*w[ijk+ii1    +kk2] ) ) )

                                    * dzi4[k+1] ) ) )


                          * dzhi4[k  ] );
        }

    // Top boundary.
    k = grid.kend;
    w2_visc[k] = 0.;
    uw_visc[k] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            w2_visc[k] += visc * ( tg0*((cg0*std::pow(w[ijk-kk4],2) + cg1*std::pow(w[ijk-kk3],2) + cg2*std::pow(w[ijk-kk2],2) + cg3*std::pow(w[ijk-kk1],2)) * dzi4[k-3])
                                 + tg1*((cg0*std::pow(w[ijk-kk3],2) + cg1*std::pow(w[ijk-kk2],2) + cg2*std::pow(w[ijk-kk1],2) + cg3*std::pow(w[ijk    ],2)) * dzi4[k-2])
                                 + tg2*((cg0*std::pow(w[ijk-kk2],2) + cg1*std::pow(w[ijk-kk1],2) + cg2*std::pow(w[ijk    ],2) + cg3*std::pow(w[ijk+kk1],2)) * dzi4[k-1])
                                 + tg3*((tg0*std::pow(w[ijk-kk2],2) + tg1*std::pow(w[ijk-kk1],2) + tg2*std::pow(w[ijk    ],2) + tg3*std::pow(w[ijk+kk1],2)) * dzi4[k  ]) ) * dzhi4top;


            uw_visc[k] += ( ( visc


                            * ( tg0*( ( cg0*( uz[ijk        -kk4] * ( ci0*w[ijk-ii2    -kk4] + ci1*w[ijk-ii1    -kk4] + ci2*w[ijk        -kk4] + ci3*w[ijk+ii1    -kk4] ) )
                                      + cg1*( uz[ijk        -kk3] * ( ci0*w[ijk-ii2    -kk3] + ci1*w[ijk-ii1    -kk3] + ci2*w[ijk        -kk3] + ci3*w[ijk+ii1    -kk3] ) )
                                      + cg2*( uz[ijk        -kk2] * ( ci0*w[ijk-ii2    -kk2] + ci1*w[ijk-ii1    -kk2] + ci2*w[ijk        -kk2] + ci3*w[ijk+ii1    -kk2] ) )
                                      + cg3*( uz[ijk        -kk1] * ( ci0*w[ijk-ii2    -kk1] + ci1*w[ijk-ii1    -kk1] + ci2*w[ijk        -kk1] + ci3*w[ijk+ii1    -kk1] ) ) )

                                    * dzi4[k-3] )

                              + tg1*( ( cg0*( uz[ijk        -kk3] * ( ci0*w[ijk-ii2    -kk3] + ci1*w[ijk-ii1    -kk3] + ci2*w[ijk        -kk3] + ci3*w[ijk+ii1    -kk3] ) )
                                      + cg1*( uz[ijk        -kk2] * ( ci0*w[ijk-ii2    -kk2] + ci1*w[ijk-ii1    -kk2] + ci2*w[ijk        -kk2] + ci3*w[ijk+ii1    -kk2] ) )
                                      + cg2*( uz[ijk        -kk1] * ( ci0*w[ijk-ii2    -kk1] + ci1*w[ijk-ii1    -kk1] + ci2*w[ijk        -kk1] + ci3*w[ijk+ii1    -kk1] ) )
                                      + cg3*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) ) )

                                    * dzi4[k-2] )

                              + tg2*( ( cg0*( uz[ijk        -kk2] * ( ci0*w[ijk-ii2    -kk2] + ci1*w[ijk-ii1    -kk2] + ci2*w[ijk        -kk2] + ci3*w[ijk+ii1    -kk2] ) )
                                      + cg1*( uz[ijk        -kk1] * ( ci0*w[ijk-ii2    -kk1] + ci1*w[ijk-ii1    -kk1] + ci2*w[ijk        -kk1] + ci3*w[ijk+ii1    -kk1] ) )
                                      + cg2*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) )
                                      + cg3*( uz[ijk        +kk1] * ( ci0*w[ijk-ii2    +kk1] + ci1*w[ijk-ii1    +kk1] + ci2*w[ijk        +kk1] + ci3*w[ijk+ii1    +kk1] ) ) )

                                    * dzi4[k-1] )

                              + tg3*( ( tg0*( uz[ijk        -kk2] * ( ci0*w[ijk-ii2    -kk2] + ci1*w[ijk-ii1    -kk2] + ci2*w[ijk        -kk2] + ci3*w[ijk+ii1    -kk2] ) )
                                      + tg1*( uz[ijk        -kk1] * ( ci0*w[ijk-ii2    -kk1] + ci1*w[ijk-ii1    -kk1] + ci2*w[ijk        -kk1] + ci3*w[ijk+ii1    -kk1] ) )
                                      + tg2*( uz[ijk            ] * ( ci0*w[ijk-ii2        ] + ci1*w[ijk-ii1        ] + ci2*w[ijk            ] + ci3*w[ijk+ii1        ] ) )
                                      + tg3*( uz[ijk        +kk1] * ( ci0*w[ijk-ii2    +kk1] + ci1*w[ijk-ii1    +kk1] + ci2*w[ijk        +kk1] + ci3*w[ijk+ii1    +kk1] ) ) )

                                    * dzi4[k  ] ) ) )


                          * dzhi4top );

        }

    master.sum(u2_visc , grid.kcells);
    master.sum(v2_visc , grid.kcells);
    master.sum(w2_visc , grid.kcells);
    master.sum(tke_visc, grid.kcells);
    master.sum(uw_visc , grid.kcells);

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        u2_visc [k] /= n;
        v2_visc [k] /= n;
        tke_visc[k] /= n;
    }
    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        w2_visc [k] /= n;
        uw_visc [k] /= n;
    }


    // 6. CALCULATE THE DISSIPATION TERM

    // bottom boundary
    k = grid.kstart;

    u2_diss [k] = 0.;
    v2_diss [k] = 0.;
    tke_diss[k] = 0.;

    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            u2_diss[k]  -= 2.*visc * (
                             std::pow( ( cg0*((ci0*(u[ijk-ii3]-umean[k]) + ci1*(u[ijk-ii2]-umean[k]) + ci2*(u[ijk-ii1]-umean[k]) + ci3*(u[ijk    ]-umean[k])))
                                       + cg1*((ci0*(u[ijk-ii2]-umean[k]) + ci1*(u[ijk-ii1]-umean[k]) + ci2*(u[ijk    ]-umean[k]) + ci3*(u[ijk+ii1]-umean[k])))
                                       + cg2*((ci0*(u[ijk-ii1]-umean[k]) + ci1*(u[ijk    ]-umean[k]) + ci2*(u[ijk+ii1]-umean[k]) + ci3*(u[ijk+ii2]-umean[k])))
                                       + cg3*((ci0*(u[ijk    ]-umean[k]) + ci1*(u[ijk+ii1]-umean[k]) + ci2*(u[ijk+ii2]-umean[k]) + ci3*(u[ijk+ii3]-umean[k]))) ) * cgi*dxi, 2)

                           + std::pow( ( cg0*((ci0*(u[ijk-jj3]-umean[k]) + ci1*(u[ijk-jj2]-umean[k]) + ci2*(u[ijk-jj1]-umean[k]) + ci3*(u[ijk    ]-umean[k])))
                                       + cg1*((ci0*(u[ijk-jj2]-umean[k]) + ci1*(u[ijk-jj1]-umean[k]) + ci2*(u[ijk    ]-umean[k]) + ci3*(u[ijk+jj1]-umean[k])))
                                       + cg2*((ci0*(u[ijk-jj1]-umean[k]) + ci1*(u[ijk    ]-umean[k]) + ci2*(u[ijk+jj1]-umean[k]) + ci3*(u[ijk+jj2]-umean[k])))
                                       + cg3*((ci0*(u[ijk    ]-umean[k]) + ci1*(u[ijk+jj1]-umean[k]) + ci2*(u[ijk+jj2]-umean[k]) + ci3*(u[ijk+jj3]-umean[k]))) ) * cgi*dyi, 2)

                           + std::pow( ( cg0*((bi0*(u[ijk-kk2]-umean[k-2]) + bi1*(u[ijk-kk1]-umean[k-1]) + bi2*(u[ijk    ]-umean[k  ]) + bi3*(u[ijk+kk1]-umean[k+1])))
                                       + cg1*((ci0*(u[ijk-kk2]-umean[k-2]) + ci1*(u[ijk-kk1]-umean[k-1]) + ci2*(u[ijk    ]-umean[k  ]) + ci3*(u[ijk+kk1]-umean[k+1])))
                                       + cg2*((ci0*(u[ijk-kk1]-umean[k-1]) + ci1*(u[ijk    ]-umean[k  ]) + ci2*(u[ijk+kk1]-umean[k+1]) + ci3*(u[ijk+kk2]-umean[k+2])))
                                       + cg3*((ci0*(u[ijk    ]-umean[k  ]) + ci1*(u[ijk+kk1]-umean[k+1]) + ci2*(u[ijk+kk2]-umean[k+2]) + ci3*(u[ijk+kk3]-umean[k+3]))) ) * dzi4[k], 2) );

            v2_diss[k]  -= 2.*visc * (
                             std::pow( ( cg0*((ci0*(v[ijk-ii3]-vmean[k]) + ci1*(v[ijk-ii2]-vmean[k]) + ci2*(v[ijk-ii1]-vmean[k]) + ci3*(v[ijk    ]-vmean[k])))
                                       + cg1*((ci0*(v[ijk-ii2]-vmean[k]) + ci1*(v[ijk-ii1]-vmean[k]) + ci2*(v[ijk    ]-vmean[k]) + ci3*(v[ijk+ii1]-vmean[k])))
                                       + cg2*((ci0*(v[ijk-ii1]-vmean[k]) + ci1*(v[ijk    ]-vmean[k]) + ci2*(v[ijk+ii1]-vmean[k]) + ci3*(v[ijk+ii2]-vmean[k])))
                                       + cg3*((ci0*(v[ijk    ]-vmean[k]) + ci1*(v[ijk+ii1]-vmean[k]) + ci2*(v[ijk+ii2]-vmean[k]) + ci3*(v[ijk+ii3]-vmean[k]))) ) * cgi*dxi, 2)

                           + std::pow( ( cg0*((ci0*(v[ijk-jj3]-vmean[k]) + ci1*(v[ijk-jj2]-vmean[k]) + ci2*(v[ijk-jj1]-vmean[k]) + ci3*(v[ijk    ]-vmean[k])))
                                       + cg1*((ci0*(v[ijk-jj2]-vmean[k]) + ci1*(v[ijk-jj1]-vmean[k]) + ci2*(v[ijk    ]-vmean[k]) + ci3*(v[ijk+jj1]-vmean[k])))
                                       + cg2*((ci0*(v[ijk-jj1]-vmean[k]) + ci1*(v[ijk    ]-vmean[k]) + ci2*(v[ijk+jj1]-vmean[k]) + ci3*(v[ijk+jj2]-vmean[k])))
                                       + cg3*((ci0*(v[ijk    ]-vmean[k]) + ci1*(v[ijk+jj1]-vmean[k]) + ci2*(v[ijk+jj2]-vmean[k]) + ci3*(v[ijk+jj3]-vmean[k]))) ) * cgi*dyi, 2)

                           + std::pow( ( cg0*((bi0*(v[ijk-kk2]-vmean[k-2]) + bi1*(v[ijk-kk1]-vmean[k-1]) + bi2*(v[ijk    ]-vmean[k  ]) + bi3*(v[ijk+kk1]-vmean[k+1])))
                                       + cg1*((ci0*(v[ijk-kk2]-vmean[k-2]) + ci1*(v[ijk-kk1]-vmean[k-1]) + ci2*(v[ijk    ]-vmean[k  ]) + ci3*(v[ijk+kk1]-vmean[k+1])))
                                       + cg2*((ci0*(v[ijk-kk1]-vmean[k-1]) + ci1*(v[ijk    ]-vmean[k  ]) + ci2*(v[ijk+kk1]-vmean[k+1]) + ci3*(v[ijk+kk2]-vmean[k+2])))
                                       + cg3*((ci0*(v[ijk    ]-vmean[k  ]) + ci1*(v[ijk+kk1]-vmean[k+1]) + ci2*(v[ijk+kk2]-vmean[k+2]) + ci3*(v[ijk+kk3]-vmean[k+3]))) ) * dzi4[k], 2) );

            tke_diss[k] -= visc * (
                                    std::pow( (cg0*w[ijk-ii1] + cg1*w[ijk] + cg2*w[ijk+ii1] + cg3*w[ijk+ii2]) * cgi*dxi, 2)
                                  + std::pow( (cg0*w[ijk-jj1] + cg1*w[ijk] + cg2*w[ijk+jj1] + cg3*w[ijk+jj2]) * cgi*dyi, 2)
                                  + std::pow( (cg0*w[ijk-kk1] + cg1*w[ijk] + cg2*w[ijk+kk1] + cg3*w[ijk+kk2]) * dzi4[k], 2) );
        }
    tke_diss[k] += 0.5*(u2_diss[k] + v2_diss[k]);

    // interior
    for (int k=grid.kstart+1; k<grid.kend-1; ++k)
    {
        u2_diss [k] = 0.;
        v2_diss [k] = 0.;
        tke_diss[k] = 0.;

        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                u2_diss[k]  -= 2.*visc * (
                                 std::pow( ( cg0*((ci0*(u[ijk-ii3]-umean[k]) + ci1*(u[ijk-ii2]-umean[k]) + ci2*(u[ijk-ii1]-umean[k]) + ci3*(u[ijk    ]-umean[k])))
                                           + cg1*((ci0*(u[ijk-ii2]-umean[k]) + ci1*(u[ijk-ii1]-umean[k]) + ci2*(u[ijk    ]-umean[k]) + ci3*(u[ijk+ii1]-umean[k])))
                                           + cg2*((ci0*(u[ijk-ii1]-umean[k]) + ci1*(u[ijk    ]-umean[k]) + ci2*(u[ijk+ii1]-umean[k]) + ci3*(u[ijk+ii2]-umean[k])))
                                           + cg3*((ci0*(u[ijk    ]-umean[k]) + ci1*(u[ijk+ii1]-umean[k]) + ci2*(u[ijk+ii2]-umean[k]) + ci3*(u[ijk+ii3]-umean[k]))) ) * cgi*dxi, 2)

                               + std::pow( ( cg0*((ci0*(u[ijk-jj3]-umean[k]) + ci1*(u[ijk-jj2]-umean[k]) + ci2*(u[ijk-jj1]-umean[k]) + ci3*(u[ijk    ]-umean[k])))
                                           + cg1*((ci0*(u[ijk-jj2]-umean[k]) + ci1*(u[ijk-jj1]-umean[k]) + ci2*(u[ijk    ]-umean[k]) + ci3*(u[ijk+jj1]-umean[k])))
                                           + cg2*((ci0*(u[ijk-jj1]-umean[k]) + ci1*(u[ijk    ]-umean[k]) + ci2*(u[ijk+jj1]-umean[k]) + ci3*(u[ijk+jj2]-umean[k])))
                                           + cg3*((ci0*(u[ijk    ]-umean[k]) + ci1*(u[ijk+jj1]-umean[k]) + ci2*(u[ijk+jj2]-umean[k]) + ci3*(u[ijk+jj3]-umean[k]))) ) * cgi*dyi, 2)

                               + std::pow( ( cg0*((ci0*(u[ijk-kk3]-umean[k-3]) + ci1*(u[ijk-kk2]-umean[k-2]) + ci2*(u[ijk-kk1]-umean[k-1]) + ci3*(u[ijk    ]-umean[k  ])))
                                           + cg1*((ci0*(u[ijk-kk2]-umean[k-2]) + ci1*(u[ijk-kk1]-umean[k-1]) + ci2*(u[ijk    ]-umean[k  ]) + ci3*(u[ijk+kk1]-umean[k+1])))
                                           + cg2*((ci0*(u[ijk-kk1]-umean[k-1]) + ci1*(u[ijk    ]-umean[k  ]) + ci2*(u[ijk+kk1]-umean[k+1]) + ci3*(u[ijk+kk2]-umean[k+2])))
                                           + cg3*((ci0*(u[ijk    ]-umean[k  ]) + ci1*(u[ijk+kk1]-umean[k+1]) + ci2*(u[ijk+kk2]-umean[k+2]) + ci3*(u[ijk+kk3]-umean[k+3]))) ) * dzi4[k], 2) );

                v2_diss[k]  -= 2.*visc * (
                                 std::pow( ( cg0*((ci0*(v[ijk-ii3]-vmean[k]) + ci1*(v[ijk-ii2]-vmean[k]) + ci2*(v[ijk-ii1]-vmean[k]) + ci3*(v[ijk    ]-vmean[k])))
                                           + cg1*((ci0*(v[ijk-ii2]-vmean[k]) + ci1*(v[ijk-ii1]-vmean[k]) + ci2*(v[ijk    ]-vmean[k]) + ci3*(v[ijk+ii1]-vmean[k])))
                                           + cg2*((ci0*(v[ijk-ii1]-vmean[k]) + ci1*(v[ijk    ]-vmean[k]) + ci2*(v[ijk+ii1]-vmean[k]) + ci3*(v[ijk+ii2]-vmean[k])))
                                           + cg3*((ci0*(v[ijk    ]-vmean[k]) + ci1*(v[ijk+ii1]-vmean[k]) + ci2*(v[ijk+ii2]-vmean[k]) + ci3*(v[ijk+ii3]-vmean[k]))) ) * cgi*dxi, 2)

                               + std::pow( ( cg0*((ci0*(v[ijk-jj3]-vmean[k]) + ci1*(v[ijk-jj2]-vmean[k]) + ci2*(v[ijk-jj1]-vmean[k]) + ci3*(v[ijk    ]-vmean[k])))
                                           + cg1*((ci0*(v[ijk-jj2]-vmean[k]) + ci1*(v[ijk-jj1]-vmean[k]) + ci2*(v[ijk    ]-vmean[k]) + ci3*(v[ijk+jj1]-vmean[k])))
                                           + cg2*((ci0*(v[ijk-jj1]-vmean[k]) + ci1*(v[ijk    ]-vmean[k]) + ci2*(v[ijk+jj1]-vmean[k]) + ci3*(v[ijk+jj2]-vmean[k])))
                                           + cg3*((ci0*(v[ijk    ]-vmean[k]) + ci1*(v[ijk+jj1]-vmean[k]) + ci2*(v[ijk+jj2]-vmean[k]) + ci3*(v[ijk+jj3]-vmean[k]))) ) * cgi*dyi, 2)

                               + std::pow( ( cg0*((ci0*(v[ijk-kk3]-vmean[k-3]) + ci1*(v[ijk-kk2]-vmean[k-2]) + ci2*(v[ijk-kk1]-vmean[k-1]) + ci3*(v[ijk    ]-vmean[k  ])))
                                           + cg1*((ci0*(v[ijk-kk2]-vmean[k-2]) + ci1*(v[ijk-kk1]-vmean[k-1]) + ci2*(v[ijk    ]-vmean[k  ]) + ci3*(v[ijk+kk1]-vmean[k+1])))
                                           + cg2*((ci0*(v[ijk-kk1]-vmean[k-1]) + ci1*(v[ijk    ]-vmean[k  ]) + ci2*(v[ijk+kk1]-vmean[k+1]) + ci3*(v[ijk+kk2]-vmean[k+2])))
                                           + cg3*((ci0*(v[ijk    ]-vmean[k  ]) + ci1*(v[ijk+kk1]-vmean[k+1]) + ci2*(v[ijk+kk2]-vmean[k+2]) + ci3*(v[ijk+kk3]-vmean[k+3]))) ) * dzi4[k], 2) );

                tke_diss[k] -= visc * (
                                   std::pow( (cg0*w[ijk-ii1] + cg1*w[ijk] + cg2*w[ijk+ii1] + cg3*w[ijk+ii2]) * cgi*dxi, 2)
                                 + std::pow( (cg0*w[ijk-jj1] + cg1*w[ijk] + cg2*w[ijk+jj1] + cg3*w[ijk+jj2]) * cgi*dyi, 2)
                                 + std::pow( (cg0*w[ijk-kk1] + cg1*w[ijk] + cg2*w[ijk+kk1] + cg3*w[ijk+kk2]) * dzi4[k], 2) );
            }
        tke_diss[k] += 0.5*(u2_diss[k] + v2_diss[k]);
    }

    // top boundary
    k = grid.kend-1;

    u2_diss [k] = 0.;
    v2_diss [k] = 0.;
    tke_diss[k] = 0.;

    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            u2_diss[k]  -= 2.*visc * (
                             std::pow( ( cg0*((ci0*(u[ijk-ii3]-umean[k]) + ci1*(u[ijk-ii2]-umean[k]) + ci2*(u[ijk-ii1]-umean[k]) + ci3*(u[ijk    ]-umean[k])))
                                       + cg1*((ci0*(u[ijk-ii2]-umean[k]) + ci1*(u[ijk-ii1]-umean[k]) + ci2*(u[ijk    ]-umean[k]) + ci3*(u[ijk+ii1]-umean[k])))
                                       + cg2*((ci0*(u[ijk-ii1]-umean[k]) + ci1*(u[ijk    ]-umean[k]) + ci2*(u[ijk+ii1]-umean[k]) + ci3*(u[ijk+ii2]-umean[k])))
                                       + cg3*((ci0*(u[ijk    ]-umean[k]) + ci1*(u[ijk+ii1]-umean[k]) + ci2*(u[ijk+ii2]-umean[k]) + ci3*(u[ijk+ii3]-umean[k]))) ) * cgi*dxi, 2)

                           + std::pow( ( cg0*((ci0*(u[ijk-jj3]-umean[k]) + ci1*(u[ijk-jj2]-umean[k]) + ci2*(u[ijk-jj1]-umean[k]) + ci3*(u[ijk    ]-umean[k])))
                                       + cg1*((ci0*(u[ijk-jj2]-umean[k]) + ci1*(u[ijk-jj1]-umean[k]) + ci2*(u[ijk    ]-umean[k]) + ci3*(u[ijk+jj1]-umean[k])))
                                       + cg2*((ci0*(u[ijk-jj1]-umean[k]) + ci1*(u[ijk    ]-umean[k]) + ci2*(u[ijk+jj1]-umean[k]) + ci3*(u[ijk+jj2]-umean[k])))
                                       + cg3*((ci0*(u[ijk    ]-umean[k]) + ci1*(u[ijk+jj1]-umean[k]) + ci2*(u[ijk+jj2]-umean[k]) + ci3*(u[ijk+jj3]-umean[k]))) ) * cgi*dyi, 2)

                           + std::pow( ( cg0*((ci0*(u[ijk-kk3]-umean[k-3]) + ci1*(u[ijk-kk2]-umean[k-2]) + ci2*(u[ijk-kk1]-umean[k-1]) + ci3*(u[ijk    ]-umean[k  ])))
                                       + cg1*((ci0*(u[ijk-kk2]-umean[k-2]) + ci1*(u[ijk-kk1]-umean[k-1]) + ci2*(u[ijk    ]-umean[k  ]) + ci3*(u[ijk+kk1]-umean[k+1])))
                                       + cg2*((ci0*(u[ijk-kk1]-umean[k-1]) + ci1*(u[ijk    ]-umean[k  ]) + ci2*(u[ijk+kk1]-umean[k+1]) + ci3*(u[ijk+kk2]-umean[k+2])))
                                       + cg3*((ti0*(u[ijk-kk1]-umean[k-1]) + ti1*(u[ijk    ]-umean[k  ]) + ti2*(u[ijk+kk1]-umean[k+1]) + ti3*(u[ijk+kk2]-umean[k+2]))) ) * dzi4[k], 2) );

            v2_diss[k]  -= 2.*visc * (
                             std::pow( ( cg0*((ci0*(v[ijk-ii3]-vmean[k]) + ci1*(v[ijk-ii2]-vmean[k]) + ci2*(v[ijk-ii1]-vmean[k]) + ci3*(v[ijk    ]-vmean[k])))
                                       + cg1*((ci0*(v[ijk-ii2]-vmean[k]) + ci1*(v[ijk-ii1]-vmean[k]) + ci2*(v[ijk    ]-vmean[k]) + ci3*(v[ijk+ii1]-vmean[k])))
                                       + cg2*((ci0*(v[ijk-ii1]-vmean[k]) + ci1*(v[ijk    ]-vmean[k]) + ci2*(v[ijk+ii1]-vmean[k]) + ci3*(v[ijk+ii2]-vmean[k])))
                                       + cg3*((ci0*(v[ijk    ]-vmean[k]) + ci1*(v[ijk+ii1]-vmean[k]) + ci2*(v[ijk+ii2]-vmean[k]) + ci3*(v[ijk+ii3]-vmean[k]))) ) * cgi*dxi, 2)

                           + std::pow( ( cg0*((ci0*(v[ijk-jj3]-vmean[k]) + ci1*(v[ijk-jj2]-vmean[k]) + ci2*(v[ijk-jj1]-vmean[k]) + ci3*(v[ijk    ]-vmean[k])))
                                       + cg1*((ci0*(v[ijk-jj2]-vmean[k]) + ci1*(v[ijk-jj1]-vmean[k]) + ci2*(v[ijk    ]-vmean[k]) + ci3*(v[ijk+jj1]-vmean[k])))
                                       + cg2*((ci0*(v[ijk-jj1]-vmean[k]) + ci1*(v[ijk    ]-vmean[k]) + ci2*(v[ijk+jj1]-vmean[k]) + ci3*(v[ijk+jj2]-vmean[k])))
                                       + cg3*((ci0*(v[ijk    ]-vmean[k]) + ci1*(v[ijk+jj1]-vmean[k]) + ci2*(v[ijk+jj2]-vmean[k]) + ci3*(v[ijk+jj3]-vmean[k]))) ) * cgi*dyi, 2)

                           + std::pow( ( cg0*((ci0*(v[ijk-kk3]-vmean[k-3]) + ci1*(v[ijk-kk2]-vmean[k-2]) + ci2*(v[ijk-kk1]-vmean[k-1]) + ci3*(v[ijk    ]-vmean[k  ])))
                                       + cg1*((ci0*(v[ijk-kk2]-vmean[k-2]) + ci1*(v[ijk-kk1]-vmean[k-1]) + ci2*(v[ijk    ]-vmean[k  ]) + ci3*(v[ijk+kk1]-vmean[k+1])))
                                       + cg2*((ci0*(v[ijk-kk1]-vmean[k-1]) + ci1*(v[ijk    ]-vmean[k  ]) + ci2*(v[ijk+kk1]-vmean[k+1]) + ci3*(v[ijk+kk2]-vmean[k+2])))
                                       + cg3*((ti0*(v[ijk-kk1]-vmean[k-1]) + ti1*(v[ijk    ]-vmean[k  ]) + ti2*(v[ijk+kk1]-vmean[k+1]) + ti3*(v[ijk+kk2]-vmean[k+2]))) ) * dzi4[k], 2) );

            tke_diss[k] -= visc * (
                               std::pow( (cg0*w[ijk-ii1] + cg1*w[ijk] + cg2*w[ijk+ii1] + cg3*w[ijk+ii2]) * cgi*dxi, 2)
                             + std::pow( (cg0*w[ijk-jj1] + cg1*w[ijk] + cg2*w[ijk+jj1] + cg3*w[ijk+jj2]) * cgi*dyi, 2)
                             + std::pow( (cg0*w[ijk-kk1] + cg1*w[ijk] + cg2*w[ijk+kk1] + cg3*w[ijk+kk2]) * dzi4[k], 2) );
        }
    tke_diss[k] += 0.5*(u2_diss[k] + v2_diss[k]);

    // calculate the w2 budget term
    for (int k=grid.kstart+1; k<grid.kend; ++k)
    {
        w2_diss[k] = 0.;
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                w2_diss[k]  -= 2.*visc * (
                                 std::pow( ( cg0*(ci0*w[ijk-ii3] + ci1*w[ijk-ii2] + ci2*w[ijk-ii1] + ci3*w[ijk    ])
                                           + cg1*(ci0*w[ijk-ii2] + ci1*w[ijk-ii1] + ci2*w[ijk    ] + ci3*w[ijk+ii1])
                                           + cg2*(ci0*w[ijk-ii1] + ci1*w[ijk    ] + ci2*w[ijk+ii1] + ci3*w[ijk+ii2])
                                           + cg3*(ci0*w[ijk    ] + ci1*w[ijk+ii1] + ci2*w[ijk+ii2] + ci3*w[ijk+ii3]) ) * cgi*dxi, 2)

                               + std::pow( ( cg0*(ci0*w[ijk-jj3] + ci1*w[ijk-jj2] + ci2*w[ijk-jj1] + ci3*w[ijk    ])
                                           + cg1*(ci0*w[ijk-jj2] + ci1*w[ijk-jj1] + ci2*w[ijk    ] + ci3*w[ijk+jj1])
                                           + cg2*(ci0*w[ijk-jj1] + ci1*w[ijk    ] + ci2*w[ijk+jj1] + ci3*w[ijk+jj2])
                                           + cg3*(ci0*w[ijk    ] + ci1*w[ijk+jj1] + ci2*w[ijk+jj2] + ci3*w[ijk+jj3]) ) * cgi*dyi, 2)

                               + std::pow( ( cg0*(ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ])
                                           + cg1*(ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1])
                                           + cg2*(ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])
                                           + cg3*(ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3]) ) * dzhi4[k], 2) );
            }
    }

    k = grid.kstart;
    uw_diss[k] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;

            uw_diss[k] -= ( ( 2 * visc )
            
            
                          * ( ( ( ( cg0*( ci0*( ci0*( u[ijk-ii3-kk2] - umean[k-2] ) + ci1*( u[ijk-ii2-kk2] - umean[k-2] ) + ci2*( u[ijk-ii1-kk2] - umean[k-2] ) + ci3*( u[ijk    -kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-ii3-kk1] - umean[k-1] ) + ci1*( u[ijk-ii2-kk1] - umean[k-1] ) + ci2*( u[ijk-ii1-kk1] - umean[k-1] ) + ci3*( u[ijk    -kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-ii3    ] - umean[k  ] ) + ci1*( u[ijk-ii2    ] - umean[k  ] ) + ci2*( u[ijk-ii1    ] - umean[k  ] ) + ci3*( u[ijk        ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-ii3+kk1] - umean[k+1] ) + ci1*( u[ijk-ii2+kk1] - umean[k+1] ) + ci2*( u[ijk-ii1+kk1] - umean[k+1] ) + ci3*( u[ijk    +kk1] - umean[k+1] ) ) )
            
                                  + cg1*( ci0*( ci0*( u[ijk-ii2-kk2] - umean[k-2] ) + ci1*( u[ijk-ii1-kk2] - umean[k-2] ) + ci2*( u[ijk    -kk2] - umean[k-2] ) + ci3*( u[ijk+ii1-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-ii2-kk1] - umean[k-1] ) + ci1*( u[ijk-ii1-kk1] - umean[k-1] ) + ci2*( u[ijk    -kk1] - umean[k-1] ) + ci3*( u[ijk+ii1-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-ii2    ] - umean[k  ] ) + ci1*( u[ijk-ii1    ] - umean[k  ] ) + ci2*( u[ijk        ] - umean[k  ] ) + ci3*( u[ijk+ii1    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-ii2+kk1] - umean[k+1] ) + ci1*( u[ijk-ii1+kk1] - umean[k+1] ) + ci2*( u[ijk    +kk1] - umean[k+1] ) + ci3*( u[ijk+ii1+kk1] - umean[k+1] ) ) )
            
                                  + cg2*( ci0*( ci0*( u[ijk-ii1-kk2] - umean[k-2] ) + ci1*( u[ijk    -kk2] - umean[k-2] ) + ci2*( u[ijk+ii1-kk2] - umean[k-2] ) + ci3*( u[ijk+ii2-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-ii1-kk1] - umean[k-1] ) + ci1*( u[ijk    -kk1] - umean[k-1] ) + ci2*( u[ijk+ii1-kk1] - umean[k-1] ) + ci3*( u[ijk+ii2-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-ii1    ] - umean[k  ] ) + ci1*( u[ijk        ] - umean[k  ] ) + ci2*( u[ijk+ii1    ] - umean[k  ] ) + ci3*( u[ijk+ii2    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-ii1+kk1] - umean[k+1] ) + ci1*( u[ijk    +kk1] - umean[k+1] ) + ci2*( u[ijk+ii1+kk1] - umean[k+1] ) + ci3*( u[ijk+ii2+kk1] - umean[k+1] ) ) )
            
                                  + cg3*( ci0*( ci0*( u[ijk    -kk2] - umean[k-2] ) + ci1*( u[ijk+ii1-kk2] - umean[k-2] ) + ci2*( u[ijk+ii2-kk2] - umean[k-2] ) + ci3*( u[ijk+ii3-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk    -kk1] - umean[k-1] ) + ci1*( u[ijk+ii1-kk1] - umean[k-1] ) + ci2*( u[ijk+ii2-kk1] - umean[k-1] ) + ci3*( u[ijk+ii3-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk        ] - umean[k  ] ) + ci1*( u[ijk+ii1    ] - umean[k  ] ) + ci2*( u[ijk+ii2    ] - umean[k  ] ) + ci3*( u[ijk+ii3    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk    +kk1] - umean[k+1] ) + ci1*( u[ijk+ii1+kk1] - umean[k+1] ) + ci2*( u[ijk+ii2+kk1] - umean[k+1] ) + ci3*( u[ijk+ii3+kk1] - umean[k+1] ) ) ) )
            
            
                                * cgi*dxi )
            
            
                              * ( cg0*w[ijk-ii2] + cg1*w[ijk-ii1] + cg2*w[ijk    ] + cg3*w[ijk+ii1] ) )
            
            
                            * cgi*dxi ) );
            
            uw_diss[k] -= ( ( 2 * visc )
            
                          * ( ( ( ( cg0*( ci0*( ci0*( u[ijk-jj3-kk2] - umean[k-2] ) + ci1*( u[ijk-jj2-kk2] - umean[k-2] ) + ci2*( u[ijk-jj1-kk2] - umean[k-2] ) + ci3*( u[ijk    -kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-jj3-kk1] - umean[k-1] ) + ci1*( u[ijk-jj2-kk1] - umean[k-1] ) + ci2*( u[ijk-jj1-kk1] - umean[k-1] ) + ci3*( u[ijk    -kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-jj3    ] - umean[k  ] ) + ci1*( u[ijk-jj2    ] - umean[k  ] ) + ci2*( u[ijk-jj1    ] - umean[k  ] ) + ci3*( u[ijk        ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-jj3+kk1] - umean[k+1] ) + ci1*( u[ijk-jj2+kk1] - umean[k+1] ) + ci2*( u[ijk-jj1+kk1] - umean[k+1] ) + ci3*( u[ijk    +kk1] - umean[k+1] ) ) )
            
                                  + cg1*( ci0*( ci0*( u[ijk-jj2-kk2] - umean[k-2] ) + ci1*( u[ijk-jj1-kk2] - umean[k-2] ) + ci2*( u[ijk    -kk2] - umean[k-2] ) + ci3*( u[ijk+jj1-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-jj2-kk1] - umean[k-1] ) + ci1*( u[ijk-jj1-kk1] - umean[k-1] ) + ci2*( u[ijk    -kk1] - umean[k-1] ) + ci3*( u[ijk+jj1-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-jj2    ] - umean[k  ] ) + ci1*( u[ijk-jj1    ] - umean[k  ] ) + ci2*( u[ijk        ] - umean[k  ] ) + ci3*( u[ijk+jj1    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-jj2+kk1] - umean[k+1] ) + ci1*( u[ijk-jj1+kk1] - umean[k+1] ) + ci2*( u[ijk    +kk1] - umean[k+1] ) + ci3*( u[ijk+jj1+kk1] - umean[k+1] ) ) )
            
                                  + cg2*( ci0*( ci0*( u[ijk-jj1-kk2] - umean[k-2] ) + ci1*( u[ijk    -kk2] - umean[k-2] ) + ci2*( u[ijk+jj1-kk2] - umean[k-2] ) + ci3*( u[ijk+jj2-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-jj1-kk1] - umean[k-1] ) + ci1*( u[ijk    -kk1] - umean[k-1] ) + ci2*( u[ijk+jj1-kk1] - umean[k-1] ) + ci3*( u[ijk+jj2-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-jj1    ] - umean[k  ] ) + ci1*( u[ijk        ] - umean[k  ] ) + ci2*( u[ijk+jj1    ] - umean[k  ] ) + ci3*( u[ijk+jj2    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-jj1+kk1] - umean[k+1] ) + ci1*( u[ijk    +kk1] - umean[k+1] ) + ci2*( u[ijk+jj1+kk1] - umean[k+1] ) + ci3*( u[ijk+jj2+kk1] - umean[k+1] ) ) )
            
                                  + cg3*( ci0*( ci0*( u[ijk    -kk2] - umean[k-2] ) + ci1*( u[ijk+jj1-kk2] - umean[k-2] ) + ci2*( u[ijk+jj2-kk2] - umean[k-2] ) + ci3*( u[ijk+jj3-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk    -kk1] - umean[k-1] ) + ci1*( u[ijk+jj1-kk1] - umean[k-1] ) + ci2*( u[ijk+jj2-kk1] - umean[k-1] ) + ci3*( u[ijk+jj3-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk        ] - umean[k  ] ) + ci1*( u[ijk+jj1    ] - umean[k  ] ) + ci2*( u[ijk+jj2    ] - umean[k  ] ) + ci3*( u[ijk+jj3    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk    +kk1] - umean[k+1] ) + ci1*( u[ijk+jj1+kk1] - umean[k+1] ) + ci2*( u[ijk+jj2+kk1] - umean[k+1] ) + ci3*( u[ijk+jj3+kk1] - umean[k+1] ) ) ) )
            
                                * cgi*dyi )
            
            
                              * ( cg0*( ci0*( ci0*w[ijk-ii2-jj3] + ci1*w[ijk-ii1-jj3] + ci2*w[ijk    -jj3] + ci3*w[ijk+ii1-jj3] )
                                      + ci1*( ci0*w[ijk-ii2-jj2] + ci1*w[ijk-ii1-jj2] + ci2*w[ijk    -jj2] + ci3*w[ijk+ii1-jj2] )
                                      + ci2*( ci0*w[ijk-ii2-jj1] + ci1*w[ijk-ii1-jj1] + ci2*w[ijk    -jj1] + ci3*w[ijk+ii1-jj1] )
                                      + ci3*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] ) )
            
                                + cg1*( ci0*( ci0*w[ijk-ii2-jj2] + ci1*w[ijk-ii1-jj2] + ci2*w[ijk    -jj2] + ci3*w[ijk+ii1-jj2] )
                                      + ci1*( ci0*w[ijk-ii2-jj1] + ci1*w[ijk-ii1-jj1] + ci2*w[ijk    -jj1] + ci3*w[ijk+ii1-jj1] )
                                      + ci2*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci3*( ci0*w[ijk-ii2+jj1] + ci1*w[ijk-ii1+jj1] + ci2*w[ijk    +jj1] + ci3*w[ijk+ii1+jj1] ) )
            
                                + cg2*( ci0*( ci0*w[ijk-ii2-jj1] + ci1*w[ijk-ii1-jj1] + ci2*w[ijk    -jj1] + ci3*w[ijk+ii1-jj1] )
                                      + ci1*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci2*( ci0*w[ijk-ii2+jj1] + ci1*w[ijk-ii1+jj1] + ci2*w[ijk    +jj1] + ci3*w[ijk+ii1+jj1] )
                                      + ci3*( ci0*w[ijk-ii2+jj2] + ci1*w[ijk-ii1+jj2] + ci2*w[ijk    +jj2] + ci3*w[ijk+ii1+jj2] ) )
            
                                + cg3*( ci0*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci1*( ci0*w[ijk-ii2+jj1] + ci1*w[ijk-ii1+jj1] + ci2*w[ijk    +jj1] + ci3*w[ijk+ii1+jj1] )
                                      + ci2*( ci0*w[ijk-ii2+jj2] + ci1*w[ijk-ii1+jj2] + ci2*w[ijk    +jj2] + ci3*w[ijk+ii1+jj2] )
                                      + ci3*( ci0*w[ijk-ii2+jj3] + ci1*w[ijk-ii1+jj3] + ci2*w[ijk    +jj3] + ci3*w[ijk+ii1+jj3] ) ) ) )
            
            
                            * cgi*dyi ) );
            
            uw_diss[k] -= ( ( 2 * visc )
            
            
                          * ( ( ( ( cg0*( u[ijk-kk2] - umean[k-2] ) + cg1*( u[ijk-kk1] - umean[k-1] ) + cg2*( u[ijk    ] - umean[k  ] ) + cg3*( u[ijk+kk1] - umean[k+1] ) ) * dzhi4[k] )
            
            
                              * ( bg0*( bi0*( ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk    -kk1] + ci3*w[ijk+ii1-kk1] )
                                      + bi1*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + bi2*( ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk    +kk1] + ci3*w[ijk+ii1+kk1] )
                                      + bi3*( ci0*w[ijk-ii2+kk2] + ci1*w[ijk-ii1+kk2] + ci2*w[ijk    +kk2] + ci3*w[ijk+ii1+kk2] ) )
            
                                + bg1*( ci0*( ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk    -kk1] + ci3*w[ijk+ii1-kk1] )
                                      + ci1*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci2*( ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk    +kk1] + ci3*w[ijk+ii1+kk1] )
                                      + ci3*( ci0*w[ijk-ii2+kk2] + ci1*w[ijk-ii1+kk2] + ci2*w[ijk    +kk2] + ci3*w[ijk+ii1+kk2] ) )
            
                                + bg2*( ci0*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci1*( ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk    +kk1] + ci3*w[ijk+ii1+kk1] )
                                      + ci2*( ci0*w[ijk-ii2+kk2] + ci1*w[ijk-ii1+kk2] + ci2*w[ijk    +kk2] + ci3*w[ijk+ii1+kk2] )
                                      + ci3*( ci0*w[ijk-ii2+kk3] + ci1*w[ijk-ii1+kk3] + ci2*w[ijk    +kk3] + ci3*w[ijk+ii1+kk3] ) )
            
                                + bg3*( ci0*( ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk    +kk1] + ci3*w[ijk+ii1+kk1] )
                                      + ci1*( ci0*w[ijk-ii2+kk2] + ci1*w[ijk-ii1+kk2] + ci2*w[ijk    +kk2] + ci3*w[ijk+ii1+kk2] )
                                      + ci2*( ci0*w[ijk-ii2+kk3] + ci1*w[ijk-ii1+kk3] + ci2*w[ijk    +kk3] + ci3*w[ijk+ii1+kk3] )
                                      + ci3*( ci0*w[ijk-ii2+kk4] + ci1*w[ijk-ii1+kk4] + ci2*w[ijk    +kk4] + ci3*w[ijk+ii1+kk4] ) ) ) )
            
            
                            * dzhi4bot ) );
        }

    k = grid.kstart+1;
    uw_diss[k] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            
            uw_diss[k] -= ( ( 2 * visc )
            
            
                          * ( ( ( ( cg0*( ci0*( ci0*( u[ijk-ii3-kk2] - umean[k-2] ) + ci1*( u[ijk-ii2-kk2] - umean[k-2] ) + ci2*( u[ijk-ii1-kk2] - umean[k-2] ) + ci3*( u[ijk    -kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-ii3-kk1] - umean[k-1] ) + ci1*( u[ijk-ii2-kk1] - umean[k-1] ) + ci2*( u[ijk-ii1-kk1] - umean[k-1] ) + ci3*( u[ijk    -kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-ii3    ] - umean[k  ] ) + ci1*( u[ijk-ii2    ] - umean[k  ] ) + ci2*( u[ijk-ii1    ] - umean[k  ] ) + ci3*( u[ijk        ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-ii3+kk1] - umean[k+1] ) + ci1*( u[ijk-ii2+kk1] - umean[k+1] ) + ci2*( u[ijk-ii1+kk1] - umean[k+1] ) + ci3*( u[ijk    +kk1] - umean[k+1] ) ) )
            
                                  + cg1*( ci0*( ci0*( u[ijk-ii2-kk2] - umean[k-2] ) + ci1*( u[ijk-ii1-kk2] - umean[k-2] ) + ci2*( u[ijk    -kk2] - umean[k-2] ) + ci3*( u[ijk+ii1-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-ii2-kk1] - umean[k-1] ) + ci1*( u[ijk-ii1-kk1] - umean[k-1] ) + ci2*( u[ijk    -kk1] - umean[k-1] ) + ci3*( u[ijk+ii1-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-ii2    ] - umean[k  ] ) + ci1*( u[ijk-ii1    ] - umean[k  ] ) + ci2*( u[ijk        ] - umean[k  ] ) + ci3*( u[ijk+ii1    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-ii2+kk1] - umean[k+1] ) + ci1*( u[ijk-ii1+kk1] - umean[k+1] ) + ci2*( u[ijk    +kk1] - umean[k+1] ) + ci3*( u[ijk+ii1+kk1] - umean[k+1] ) ) )
            
                                  + cg2*( ci0*( ci0*( u[ijk-ii1-kk2] - umean[k-2] ) + ci1*( u[ijk    -kk2] - umean[k-2] ) + ci2*( u[ijk+ii1-kk2] - umean[k-2] ) + ci3*( u[ijk+ii2-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-ii1-kk1] - umean[k-1] ) + ci1*( u[ijk    -kk1] - umean[k-1] ) + ci2*( u[ijk+ii1-kk1] - umean[k-1] ) + ci3*( u[ijk+ii2-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-ii1    ] - umean[k  ] ) + ci1*( u[ijk        ] - umean[k  ] ) + ci2*( u[ijk+ii1    ] - umean[k  ] ) + ci3*( u[ijk+ii2    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-ii1+kk1] - umean[k+1] ) + ci1*( u[ijk    +kk1] - umean[k+1] ) + ci2*( u[ijk+ii1+kk1] - umean[k+1] ) + ci3*( u[ijk+ii2+kk1] - umean[k+1] ) ) )
            
                                  + cg3*( ci0*( ci0*( u[ijk    -kk2] - umean[k-2] ) + ci1*( u[ijk+ii1-kk2] - umean[k-2] ) + ci2*( u[ijk+ii2-kk2] - umean[k-2] ) + ci3*( u[ijk+ii3-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk    -kk1] - umean[k-1] ) + ci1*( u[ijk+ii1-kk1] - umean[k-1] ) + ci2*( u[ijk+ii2-kk1] - umean[k-1] ) + ci3*( u[ijk+ii3-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk        ] - umean[k  ] ) + ci1*( u[ijk+ii1    ] - umean[k  ] ) + ci2*( u[ijk+ii2    ] - umean[k  ] ) + ci3*( u[ijk+ii3    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk    +kk1] - umean[k+1] ) + ci1*( u[ijk+ii1+kk1] - umean[k+1] ) + ci2*( u[ijk+ii2+kk1] - umean[k+1] ) + ci3*( u[ijk+ii3+kk1] - umean[k+1] ) ) ) )
            
            
                                * cgi*dxi )
            
            
                              * ( cg0*w[ijk-ii2] + cg1*w[ijk-ii1] + cg2*w[ijk    ] + cg3*w[ijk+ii1] ) )
            
            
                            * cgi*dxi ) );
            
            uw_diss[k] -= ( ( 2 * visc )
            
            
                          * ( ( ( ( cg0*( ci0*( ci0*( u[ijk-jj3-kk2] - umean[k-2] ) + ci1*( u[ijk-jj2-kk2] - umean[k-2] ) + ci2*( u[ijk-jj1-kk2] - umean[k-2] ) + ci3*( u[ijk    -kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-jj3-kk1] - umean[k-1] ) + ci1*( u[ijk-jj2-kk1] - umean[k-1] ) + ci2*( u[ijk-jj1-kk1] - umean[k-1] ) + ci3*( u[ijk    -kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-jj3    ] - umean[k  ] ) + ci1*( u[ijk-jj2    ] - umean[k  ] ) + ci2*( u[ijk-jj1    ] - umean[k  ] ) + ci3*( u[ijk        ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-jj3+kk1] - umean[k+1] ) + ci1*( u[ijk-jj2+kk1] - umean[k+1] ) + ci2*( u[ijk-jj1+kk1] - umean[k+1] ) + ci3*( u[ijk    +kk1] - umean[k+1] ) ) )
            
                                  + cg1*( ci0*( ci0*( u[ijk-jj2-kk2] - umean[k-2] ) + ci1*( u[ijk-jj1-kk2] - umean[k-2] ) + ci2*( u[ijk    -kk2] - umean[k-2] ) + ci3*( u[ijk+jj1-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-jj2-kk1] - umean[k-1] ) + ci1*( u[ijk-jj1-kk1] - umean[k-1] ) + ci2*( u[ijk    -kk1] - umean[k-1] ) + ci3*( u[ijk+jj1-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-jj2    ] - umean[k  ] ) + ci1*( u[ijk-jj1    ] - umean[k  ] ) + ci2*( u[ijk        ] - umean[k  ] ) + ci3*( u[ijk+jj1    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-jj2+kk1] - umean[k+1] ) + ci1*( u[ijk-jj1+kk1] - umean[k+1] ) + ci2*( u[ijk    +kk1] - umean[k+1] ) + ci3*( u[ijk+jj1+kk1] - umean[k+1] ) ) )
            
                                  + cg2*( ci0*( ci0*( u[ijk-jj1-kk2] - umean[k-2] ) + ci1*( u[ijk    -kk2] - umean[k-2] ) + ci2*( u[ijk+jj1-kk2] - umean[k-2] ) + ci3*( u[ijk+jj2-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-jj1-kk1] - umean[k-1] ) + ci1*( u[ijk    -kk1] - umean[k-1] ) + ci2*( u[ijk+jj1-kk1] - umean[k-1] ) + ci3*( u[ijk+jj2-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-jj1    ] - umean[k  ] ) + ci1*( u[ijk        ] - umean[k  ] ) + ci2*( u[ijk+jj1    ] - umean[k  ] ) + ci3*( u[ijk+jj2    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-jj1+kk1] - umean[k+1] ) + ci1*( u[ijk    +kk1] - umean[k+1] ) + ci2*( u[ijk+jj1+kk1] - umean[k+1] ) + ci3*( u[ijk+jj2+kk1] - umean[k+1] ) ) )
            
                                  + cg3*( ci0*( ci0*( u[ijk    -kk2] - umean[k-2] ) + ci1*( u[ijk+jj1-kk2] - umean[k-2] ) + ci2*( u[ijk+jj2-kk2] - umean[k-2] ) + ci3*( u[ijk+jj3-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk    -kk1] - umean[k-1] ) + ci1*( u[ijk+jj1-kk1] - umean[k-1] ) + ci2*( u[ijk+jj2-kk1] - umean[k-1] ) + ci3*( u[ijk+jj3-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk        ] - umean[k  ] ) + ci1*( u[ijk+jj1    ] - umean[k  ] ) + ci2*( u[ijk+jj2    ] - umean[k  ] ) + ci3*( u[ijk+jj3    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk    +kk1] - umean[k+1] ) + ci1*( u[ijk+jj1+kk1] - umean[k+1] ) + ci2*( u[ijk+jj2+kk1] - umean[k+1] ) + ci3*( u[ijk+jj3+kk1] - umean[k+1] ) ) ) )
            
            
                                * cgi*dyi )
            
            
                              * ( cg0*( ci0*( ci0*w[ijk-ii2-jj3] + ci1*w[ijk-ii1-jj3] + ci2*w[ijk    -jj3] + ci3*w[ijk+ii1-jj3] )
                                      + ci1*( ci0*w[ijk-ii2-jj2] + ci1*w[ijk-ii1-jj2] + ci2*w[ijk    -jj2] + ci3*w[ijk+ii1-jj2] )
                                      + ci2*( ci0*w[ijk-ii2-jj1] + ci1*w[ijk-ii1-jj1] + ci2*w[ijk    -jj1] + ci3*w[ijk+ii1-jj1] )
                                      + ci3*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] ) )
            
                                + cg1*( ci0*( ci0*w[ijk-ii2-jj2] + ci1*w[ijk-ii1-jj2] + ci2*w[ijk    -jj2] + ci3*w[ijk+ii1-jj2] )
                                      + ci1*( ci0*w[ijk-ii2-jj1] + ci1*w[ijk-ii1-jj1] + ci2*w[ijk    -jj1] + ci3*w[ijk+ii1-jj1] )
                                      + ci2*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci3*( ci0*w[ijk-ii2+jj1] + ci1*w[ijk-ii1+jj1] + ci2*w[ijk    +jj1] + ci3*w[ijk+ii1+jj1] ) )
            
                                + cg2*( ci0*( ci0*w[ijk-ii2-jj1] + ci1*w[ijk-ii1-jj1] + ci2*w[ijk    -jj1] + ci3*w[ijk+ii1-jj1] )
                                      + ci1*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci2*( ci0*w[ijk-ii2+jj1] + ci1*w[ijk-ii1+jj1] + ci2*w[ijk    +jj1] + ci3*w[ijk+ii1+jj1] )
                                      + ci3*( ci0*w[ijk-ii2+jj2] + ci1*w[ijk-ii1+jj2] + ci2*w[ijk    +jj2] + ci3*w[ijk+ii1+jj2] ) )
            
                                + cg3*( ci0*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci1*( ci0*w[ijk-ii2+jj1] + ci1*w[ijk-ii1+jj1] + ci2*w[ijk    +jj1] + ci3*w[ijk+ii1+jj1] )
                                      + ci2*( ci0*w[ijk-ii2+jj2] + ci1*w[ijk-ii1+jj2] + ci2*w[ijk    +jj2] + ci3*w[ijk+ii1+jj2] )
                                      + ci3*( ci0*w[ijk-ii2+jj3] + ci1*w[ijk-ii1+jj3] + ci2*w[ijk    +jj3] + ci3*w[ijk+ii1+jj3] ) ) ) )
            
            
                            * cgi*dyi ) );
            
            uw_diss[k] -= ( ( 2 * visc )
            
            
                          * ( ( ( ( cg0*( u[ijk-kk2] - umean[k-2] ) + cg1*( u[ijk-kk1] - umean[k-1] ) + cg2*( u[ijk    ] - umean[k  ] ) + cg3*( u[ijk+kk1] - umean[k+1] ) ) * dzhi4[k] )
            
            
                              * ( cg0*( bi0*( ci0*w[ijk-ii2-kk2] + ci1*w[ijk-ii1-kk2] + ci2*w[ijk    -kk2] + ci3*w[ijk+ii1-kk2] )
                                      + bi1*( ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk    -kk1] + ci3*w[ijk+ii1-kk1] )
                                      + bi2*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + bi3*( ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk    +kk1] + ci3*w[ijk+ii1+kk1] ) )
            
                                + cg1*( ci0*( ci0*w[ijk-ii2-kk2] + ci1*w[ijk-ii1-kk2] + ci2*w[ijk    -kk2] + ci3*w[ijk+ii1-kk2] )
                                      + ci1*( ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk    -kk1] + ci3*w[ijk+ii1-kk1] )
                                      + ci2*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci3*( ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk    +kk1] + ci3*w[ijk+ii1+kk1] ) )
            
                                + cg2*( ci0*( ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk    -kk1] + ci3*w[ijk+ii1-kk1] )
                                      + ci1*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci2*( ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk    +kk1] + ci3*w[ijk+ii1+kk1] )
                                      + ci3*( ci0*w[ijk-ii2+kk2] + ci1*w[ijk-ii1+kk2] + ci2*w[ijk    +kk2] + ci3*w[ijk+ii1+kk2] ) )
            
                                + cg3*( ci0*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci1*( ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk    +kk1] + ci3*w[ijk+ii1+kk1] )
                                      + ci2*( ci0*w[ijk-ii2+kk2] + ci1*w[ijk-ii1+kk2] + ci2*w[ijk    +kk2] + ci3*w[ijk+ii1+kk2] )
                                      + ci3*( ci0*w[ijk-ii2+kk3] + ci1*w[ijk-ii1+kk3] + ci2*w[ijk    +kk3] + ci3*w[ijk+ii1+kk3] ) ) ) )
            
            
                            * dzhi4[k] ) ); 
        }

    for (int k=grid.kstart+2; k<grid.kend-1; ++k)
    {
        uw_diss[k] = 0.;
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
               
                uw_diss[k] -= ( ( 2 * visc )
                
                
                              * ( ( ( ( cg0*( ci0*( ci0*( u[ijk-ii3-kk2] - umean[k-2] ) + ci1*( u[ijk-ii2-kk2] - umean[k-2] ) + ci2*( u[ijk-ii1-kk2] - umean[k-2] ) + ci3*( u[ijk    -kk2] - umean[k-2] ) )
                                            + ci1*( ci0*( u[ijk-ii3-kk1] - umean[k-1] ) + ci1*( u[ijk-ii2-kk1] - umean[k-1] ) + ci2*( u[ijk-ii1-kk1] - umean[k-1] ) + ci3*( u[ijk    -kk1] - umean[k-1] ) )
                                            + ci2*( ci0*( u[ijk-ii3    ] - umean[k  ] ) + ci1*( u[ijk-ii2    ] - umean[k  ] ) + ci2*( u[ijk-ii1    ] - umean[k  ] ) + ci3*( u[ijk        ] - umean[k  ] ) )
                                            + ci3*( ci0*( u[ijk-ii3+kk1] - umean[k+1] ) + ci1*( u[ijk-ii2+kk1] - umean[k+1] ) + ci2*( u[ijk-ii1+kk1] - umean[k+1] ) + ci3*( u[ijk    +kk1] - umean[k+1] ) ) )
                
                                      + cg1*( ci0*( ci0*( u[ijk-ii2-kk2] - umean[k-2] ) + ci1*( u[ijk-ii1-kk2] - umean[k-2] ) + ci2*( u[ijk    -kk2] - umean[k-2] ) + ci3*( u[ijk+ii1-kk2] - umean[k-2] ) )
                                            + ci1*( ci0*( u[ijk-ii2-kk1] - umean[k-1] ) + ci1*( u[ijk-ii1-kk1] - umean[k-1] ) + ci2*( u[ijk    -kk1] - umean[k-1] ) + ci3*( u[ijk+ii1-kk1] - umean[k-1] ) )
                                            + ci2*( ci0*( u[ijk-ii2    ] - umean[k  ] ) + ci1*( u[ijk-ii1    ] - umean[k  ] ) + ci2*( u[ijk        ] - umean[k  ] ) + ci3*( u[ijk+ii1    ] - umean[k  ] ) )
                                            + ci3*( ci0*( u[ijk-ii2+kk1] - umean[k+1] ) + ci1*( u[ijk-ii1+kk1] - umean[k+1] ) + ci2*( u[ijk    +kk1] - umean[k+1] ) + ci3*( u[ijk+ii1+kk1] - umean[k+1] ) ) )
                
                                      + cg2*( ci0*( ci0*( u[ijk-ii1-kk2] - umean[k-2] ) + ci1*( u[ijk    -kk2] - umean[k-2] ) + ci2*( u[ijk+ii1-kk2] - umean[k-2] ) + ci3*( u[ijk+ii2-kk2] - umean[k-2] ) )
                                            + ci1*( ci0*( u[ijk-ii1-kk1] - umean[k-1] ) + ci1*( u[ijk    -kk1] - umean[k-1] ) + ci2*( u[ijk+ii1-kk1] - umean[k-1] ) + ci3*( u[ijk+ii2-kk1] - umean[k-1] ) )
                                            + ci2*( ci0*( u[ijk-ii1    ] - umean[k  ] ) + ci1*( u[ijk        ] - umean[k  ] ) + ci2*( u[ijk+ii1    ] - umean[k  ] ) + ci3*( u[ijk+ii2    ] - umean[k  ] ) )
                                            + ci3*( ci0*( u[ijk-ii1+kk1] - umean[k+1] ) + ci1*( u[ijk    +kk1] - umean[k+1] ) + ci2*( u[ijk+ii1+kk1] - umean[k+1] ) + ci3*( u[ijk+ii2+kk1] - umean[k+1] ) ) )
                
                                      + cg3*( ci0*( ci0*( u[ijk    -kk2] - umean[k-2] ) + ci1*( u[ijk+ii1-kk2] - umean[k-2] ) + ci2*( u[ijk+ii2-kk2] - umean[k-2] ) + ci3*( u[ijk+ii3-kk2] - umean[k-2] ) )
                                            + ci1*( ci0*( u[ijk    -kk1] - umean[k-1] ) + ci1*( u[ijk+ii1-kk1] - umean[k-1] ) + ci2*( u[ijk+ii2-kk1] - umean[k-1] ) + ci3*( u[ijk+ii3-kk1] - umean[k-1] ) )
                                            + ci2*( ci0*( u[ijk        ] - umean[k  ] ) + ci1*( u[ijk+ii1    ] - umean[k  ] ) + ci2*( u[ijk+ii2    ] - umean[k  ] ) + ci3*( u[ijk+ii3    ] - umean[k  ] ) )
                                            + ci3*( ci0*( u[ijk    +kk1] - umean[k+1] ) + ci1*( u[ijk+ii1+kk1] - umean[k+1] ) + ci2*( u[ijk+ii2+kk1] - umean[k+1] ) + ci3*( u[ijk+ii3+kk1] - umean[k+1] ) ) ) )
                
                
                                    * cgi*dxi )
                
                
                                  * ( cg0*w[ijk-ii2] + cg1*w[ijk-ii1] + cg2*w[ijk    ] + cg3*w[ijk+ii1] ) )
                
                
                                * cgi*dxi ) );
                
                uw_diss[k] -= ( ( 2 * visc )
                
                
                              * ( ( ( ( cg0*( ci0*( ci0*( u[ijk-jj3-kk2] - umean[k-2] ) + ci1*( u[ijk-jj2-kk2] - umean[k-2] ) + ci2*( u[ijk-jj1-kk2] - umean[k-2] ) + ci3*( u[ijk    -kk2] - umean[k-2] ) )
                                            + ci1*( ci0*( u[ijk-jj3-kk1] - umean[k-1] ) + ci1*( u[ijk-jj2-kk1] - umean[k-1] ) + ci2*( u[ijk-jj1-kk1] - umean[k-1] ) + ci3*( u[ijk    -kk1] - umean[k-1] ) )
                                            + ci2*( ci0*( u[ijk-jj3    ] - umean[k  ] ) + ci1*( u[ijk-jj2    ] - umean[k  ] ) + ci2*( u[ijk-jj1    ] - umean[k  ] ) + ci3*( u[ijk        ] - umean[k  ] ) )
                                            + ci3*( ci0*( u[ijk-jj3+kk1] - umean[k+1] ) + ci1*( u[ijk-jj2+kk1] - umean[k+1] ) + ci2*( u[ijk-jj1+kk1] - umean[k+1] ) + ci3*( u[ijk    +kk1] - umean[k+1] ) ) )
                
                                      + cg1*( ci0*( ci0*( u[ijk-jj2-kk2] - umean[k-2] ) + ci1*( u[ijk-jj1-kk2] - umean[k-2] ) + ci2*( u[ijk    -kk2] - umean[k-2] ) + ci3*( u[ijk+jj1-kk2] - umean[k-2] ) )
                                            + ci1*( ci0*( u[ijk-jj2-kk1] - umean[k-1] ) + ci1*( u[ijk-jj1-kk1] - umean[k-1] ) + ci2*( u[ijk    -kk1] - umean[k-1] ) + ci3*( u[ijk+jj1-kk1] - umean[k-1] ) )
                                            + ci2*( ci0*( u[ijk-jj2    ] - umean[k  ] ) + ci1*( u[ijk-jj1    ] - umean[k  ] ) + ci2*( u[ijk        ] - umean[k  ] ) + ci3*( u[ijk+jj1    ] - umean[k  ] ) )
                                            + ci3*( ci0*( u[ijk-jj2+kk1] - umean[k+1] ) + ci1*( u[ijk-jj1+kk1] - umean[k+1] ) + ci2*( u[ijk    +kk1] - umean[k+1] ) + ci3*( u[ijk+jj1+kk1] - umean[k+1] ) ) )
                
                                      + cg2*( ci0*( ci0*( u[ijk-jj1-kk2] - umean[k-2] ) + ci1*( u[ijk    -kk2] - umean[k-2] ) + ci2*( u[ijk+jj1-kk2] - umean[k-2] ) + ci3*( u[ijk+jj2-kk2] - umean[k-2] ) )
                                            + ci1*( ci0*( u[ijk-jj1-kk1] - umean[k-1] ) + ci1*( u[ijk    -kk1] - umean[k-1] ) + ci2*( u[ijk+jj1-kk1] - umean[k-1] ) + ci3*( u[ijk+jj2-kk1] - umean[k-1] ) )
                                            + ci2*( ci0*( u[ijk-jj1    ] - umean[k  ] ) + ci1*( u[ijk        ] - umean[k  ] ) + ci2*( u[ijk+jj1    ] - umean[k  ] ) + ci3*( u[ijk+jj2    ] - umean[k  ] ) )
                                            + ci3*( ci0*( u[ijk-jj1+kk1] - umean[k+1] ) + ci1*( u[ijk    +kk1] - umean[k+1] ) + ci2*( u[ijk+jj1+kk1] - umean[k+1] ) + ci3*( u[ijk+jj2+kk1] - umean[k+1] ) ) )
                
                                      + cg3*( ci0*( ci0*( u[ijk    -kk2] - umean[k-2] ) + ci1*( u[ijk+jj1-kk2] - umean[k-2] ) + ci2*( u[ijk+jj2-kk2] - umean[k-2] ) + ci3*( u[ijk+jj3-kk2] - umean[k-2] ) )
                                            + ci1*( ci0*( u[ijk    -kk1] - umean[k-1] ) + ci1*( u[ijk+jj1-kk1] - umean[k-1] ) + ci2*( u[ijk+jj2-kk1] - umean[k-1] ) + ci3*( u[ijk+jj3-kk1] - umean[k-1] ) )
                                            + ci2*( ci0*( u[ijk        ] - umean[k  ] ) + ci1*( u[ijk+jj1    ] - umean[k  ] ) + ci2*( u[ijk+jj2    ] - umean[k  ] ) + ci3*( u[ijk+jj3    ] - umean[k  ] ) )
                                            + ci3*( ci0*( u[ijk    +kk1] - umean[k+1] ) + ci1*( u[ijk+jj1+kk1] - umean[k+1] ) + ci2*( u[ijk+jj2+kk1] - umean[k+1] ) + ci3*( u[ijk+jj3+kk1] - umean[k+1] ) ) ) )
                
                
                                    * cgi*dyi )
                
                
                                  * ( cg0*( ci0*( ci0*w[ijk-ii2-jj3] + ci1*w[ijk-ii1-jj3] + ci2*w[ijk    -jj3] + ci3*w[ijk+ii1-jj3] )
                                          + ci1*( ci0*w[ijk-ii2-jj2] + ci1*w[ijk-ii1-jj2] + ci2*w[ijk    -jj2] + ci3*w[ijk+ii1-jj2] )
                                          + ci2*( ci0*w[ijk-ii2-jj1] + ci1*w[ijk-ii1-jj1] + ci2*w[ijk    -jj1] + ci3*w[ijk+ii1-jj1] )
                                          + ci3*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] ) )
                
                                    + cg1*( ci0*( ci0*w[ijk-ii2-jj2] + ci1*w[ijk-ii1-jj2] + ci2*w[ijk    -jj2] + ci3*w[ijk+ii1-jj2] )
                                          + ci1*( ci0*w[ijk-ii2-jj1] + ci1*w[ijk-ii1-jj1] + ci2*w[ijk    -jj1] + ci3*w[ijk+ii1-jj1] )
                                          + ci2*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                          + ci3*( ci0*w[ijk-ii2+jj1] + ci1*w[ijk-ii1+jj1] + ci2*w[ijk    +jj1] + ci3*w[ijk+ii1+jj1] ) )
                
                                    + cg2*( ci0*( ci0*w[ijk-ii2-jj1] + ci1*w[ijk-ii1-jj1] + ci2*w[ijk    -jj1] + ci3*w[ijk+ii1-jj1] )
                                          + ci1*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                          + ci2*( ci0*w[ijk-ii2+jj1] + ci1*w[ijk-ii1+jj1] + ci2*w[ijk    +jj1] + ci3*w[ijk+ii1+jj1] )
                                          + ci3*( ci0*w[ijk-ii2+jj2] + ci1*w[ijk-ii1+jj2] + ci2*w[ijk    +jj2] + ci3*w[ijk+ii1+jj2] ) )
                
                                    + cg3*( ci0*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                          + ci1*( ci0*w[ijk-ii2+jj1] + ci1*w[ijk-ii1+jj1] + ci2*w[ijk    +jj1] + ci3*w[ijk+ii1+jj1] )
                                          + ci2*( ci0*w[ijk-ii2+jj2] + ci1*w[ijk-ii1+jj2] + ci2*w[ijk    +jj2] + ci3*w[ijk+ii1+jj2] )
                                          + ci3*( ci0*w[ijk-ii2+jj3] + ci1*w[ijk-ii1+jj3] + ci2*w[ijk    +jj3] + ci3*w[ijk+ii1+jj3] ) ) ) )
                
                
                                * cgi*dyi ) );
                
                uw_diss[k] -= ( ( 2 * visc )
                
                
                              * ( ( ( ( cg0*( u[ijk-kk2] - umean[k-2] ) + cg1*( u[ijk-kk1] - umean[k-1] ) + cg2*( u[ijk    ] - umean[k  ] ) + cg3*( u[ijk+kk1] - umean[k+1] ) ) * dzhi4[k] )
                
                
                                  * ( cg0*( ci0*( ci0*w[ijk-ii2-kk3] + ci1*w[ijk-ii1-kk3] + ci2*w[ijk    -kk3] + ci3*w[ijk+ii1-kk3] )
                                          + ci1*( ci0*w[ijk-ii2-kk2] + ci1*w[ijk-ii1-kk2] + ci2*w[ijk    -kk2] + ci3*w[ijk+ii1-kk2] )
                                          + ci2*( ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk    -kk1] + ci3*w[ijk+ii1-kk1] )
                                          + ci3*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] ) )
                
                                    + cg1*( ci0*( ci0*w[ijk-ii2-kk2] + ci1*w[ijk-ii1-kk2] + ci2*w[ijk    -kk2] + ci3*w[ijk+ii1-kk2] )
                                          + ci1*( ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk    -kk1] + ci3*w[ijk+ii1-kk1] )
                                          + ci2*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                          + ci3*( ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk    +kk1] + ci3*w[ijk+ii1+kk1] ) )
                
                                    + cg2*( ci0*( ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk    -kk1] + ci3*w[ijk+ii1-kk1] )
                                          + ci1*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                          + ci2*( ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk    +kk1] + ci3*w[ijk+ii1+kk1] )
                                          + ci3*( ci0*w[ijk-ii2+kk2] + ci1*w[ijk-ii1+kk2] + ci2*w[ijk    +kk2] + ci3*w[ijk+ii1+kk2] ) )
                
                                    + cg3*( ci0*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                          + ci1*( ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk    +kk1] + ci3*w[ijk+ii1+kk1] )
                                          + ci2*( ci0*w[ijk-ii2+kk2] + ci1*w[ijk-ii1+kk2] + ci2*w[ijk    +kk2] + ci3*w[ijk+ii1+kk2] )
                                          + ci3*( ci0*w[ijk-ii2+kk3] + ci1*w[ijk-ii1+kk3] + ci2*w[ijk    +kk3] + ci3*w[ijk+ii1+kk3] ) ) ) )
                
                
                                * dzhi4[k] ) );
               
            }
    }

    k = grid.kend-1;
    uw_diss[k] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;

            uw_diss[k] -= ( ( 2 * visc )
            
            
                          * ( ( ( ( cg0*( ci0*( ci0*( u[ijk-ii3-kk2] - umean[k-2] ) + ci1*( u[ijk-ii2-kk2] - umean[k-2] ) + ci2*( u[ijk-ii1-kk2] - umean[k-2] ) + ci3*( u[ijk    -kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-ii3-kk1] - umean[k-1] ) + ci1*( u[ijk-ii2-kk1] - umean[k-1] ) + ci2*( u[ijk-ii1-kk1] - umean[k-1] ) + ci3*( u[ijk    -kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-ii3    ] - umean[k  ] ) + ci1*( u[ijk-ii2    ] - umean[k  ] ) + ci2*( u[ijk-ii1    ] - umean[k  ] ) + ci3*( u[ijk        ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-ii3+kk1] - umean[k+1] ) + ci1*( u[ijk-ii2+kk1] - umean[k+1] ) + ci2*( u[ijk-ii1+kk1] - umean[k+1] ) + ci3*( u[ijk    +kk1] - umean[k+1] ) ) )
            
                                  + cg1*( ci0*( ci0*( u[ijk-ii2-kk2] - umean[k-2] ) + ci1*( u[ijk-ii1-kk2] - umean[k-2] ) + ci2*( u[ijk    -kk2] - umean[k-2] ) + ci3*( u[ijk+ii1-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-ii2-kk1] - umean[k-1] ) + ci1*( u[ijk-ii1-kk1] - umean[k-1] ) + ci2*( u[ijk    -kk1] - umean[k-1] ) + ci3*( u[ijk+ii1-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-ii2    ] - umean[k  ] ) + ci1*( u[ijk-ii1    ] - umean[k  ] ) + ci2*( u[ijk        ] - umean[k  ] ) + ci3*( u[ijk+ii1    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-ii2+kk1] - umean[k+1] ) + ci1*( u[ijk-ii1+kk1] - umean[k+1] ) + ci2*( u[ijk    +kk1] - umean[k+1] ) + ci3*( u[ijk+ii1+kk1] - umean[k+1] ) ) )
            
                                  + cg2*( ci0*( ci0*( u[ijk-ii1-kk2] - umean[k-2] ) + ci1*( u[ijk    -kk2] - umean[k-2] ) + ci2*( u[ijk+ii1-kk2] - umean[k-2] ) + ci3*( u[ijk+ii2-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-ii1-kk1] - umean[k-1] ) + ci1*( u[ijk    -kk1] - umean[k-1] ) + ci2*( u[ijk+ii1-kk1] - umean[k-1] ) + ci3*( u[ijk+ii2-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-ii1    ] - umean[k  ] ) + ci1*( u[ijk        ] - umean[k  ] ) + ci2*( u[ijk+ii1    ] - umean[k  ] ) + ci3*( u[ijk+ii2    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-ii1+kk1] - umean[k+1] ) + ci1*( u[ijk    +kk1] - umean[k+1] ) + ci2*( u[ijk+ii1+kk1] - umean[k+1] ) + ci3*( u[ijk+ii2+kk1] - umean[k+1] ) ) )
            
                                  + cg3*( ci0*( ci0*( u[ijk    -kk2] - umean[k-2] ) + ci1*( u[ijk+ii1-kk2] - umean[k-2] ) + ci2*( u[ijk+ii2-kk2] - umean[k-2] ) + ci3*( u[ijk+ii3-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk    -kk1] - umean[k-1] ) + ci1*( u[ijk+ii1-kk1] - umean[k-1] ) + ci2*( u[ijk+ii2-kk1] - umean[k-1] ) + ci3*( u[ijk+ii3-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk        ] - umean[k  ] ) + ci1*( u[ijk+ii1    ] - umean[k  ] ) + ci2*( u[ijk+ii2    ] - umean[k  ] ) + ci3*( u[ijk+ii3    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk    +kk1] - umean[k+1] ) + ci1*( u[ijk+ii1+kk1] - umean[k+1] ) + ci2*( u[ijk+ii2+kk1] - umean[k+1] ) + ci3*( u[ijk+ii3+kk1] - umean[k+1] ) ) ) )
            
            
                                * cgi*dxi )
            
            
                              * ( cg0*w[ijk-ii2] + cg1*w[ijk-ii1] + cg2*w[ijk    ] + cg3*w[ijk+ii1] ) )
            
            
                            * cgi*dxi ) );
            
            uw_diss[k] -= ( ( 2 * visc )
            
            
                          * ( ( ( ( cg0*( ci0*( ci0*( u[ijk-jj3-kk2] - umean[k-2] ) + ci1*( u[ijk-jj2-kk2] - umean[k-2] ) + ci2*( u[ijk-jj1-kk2] - umean[k-2] ) + ci3*( u[ijk    -kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-jj3-kk1] - umean[k-1] ) + ci1*( u[ijk-jj2-kk1] - umean[k-1] ) + ci2*( u[ijk-jj1-kk1] - umean[k-1] ) + ci3*( u[ijk    -kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-jj3    ] - umean[k  ] ) + ci1*( u[ijk-jj2    ] - umean[k  ] ) + ci2*( u[ijk-jj1    ] - umean[k  ] ) + ci3*( u[ijk        ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-jj3+kk1] - umean[k+1] ) + ci1*( u[ijk-jj2+kk1] - umean[k+1] ) + ci2*( u[ijk-jj1+kk1] - umean[k+1] ) + ci3*( u[ijk    +kk1] - umean[k+1] ) ) )
            
                                  + cg1*( ci0*( ci0*( u[ijk-jj2-kk2] - umean[k-2] ) + ci1*( u[ijk-jj1-kk2] - umean[k-2] ) + ci2*( u[ijk    -kk2] - umean[k-2] ) + ci3*( u[ijk+jj1-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-jj2-kk1] - umean[k-1] ) + ci1*( u[ijk-jj1-kk1] - umean[k-1] ) + ci2*( u[ijk    -kk1] - umean[k-1] ) + ci3*( u[ijk+jj1-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-jj2    ] - umean[k  ] ) + ci1*( u[ijk-jj1    ] - umean[k  ] ) + ci2*( u[ijk        ] - umean[k  ] ) + ci3*( u[ijk+jj1    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-jj2+kk1] - umean[k+1] ) + ci1*( u[ijk-jj1+kk1] - umean[k+1] ) + ci2*( u[ijk    +kk1] - umean[k+1] ) + ci3*( u[ijk+jj1+kk1] - umean[k+1] ) ) )
            
                                  + cg2*( ci0*( ci0*( u[ijk-jj1-kk2] - umean[k-2] ) + ci1*( u[ijk    -kk2] - umean[k-2] ) + ci2*( u[ijk+jj1-kk2] - umean[k-2] ) + ci3*( u[ijk+jj2-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-jj1-kk1] - umean[k-1] ) + ci1*( u[ijk    -kk1] - umean[k-1] ) + ci2*( u[ijk+jj1-kk1] - umean[k-1] ) + ci3*( u[ijk+jj2-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-jj1    ] - umean[k  ] ) + ci1*( u[ijk        ] - umean[k  ] ) + ci2*( u[ijk+jj1    ] - umean[k  ] ) + ci3*( u[ijk+jj2    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-jj1+kk1] - umean[k+1] ) + ci1*( u[ijk    +kk1] - umean[k+1] ) + ci2*( u[ijk+jj1+kk1] - umean[k+1] ) + ci3*( u[ijk+jj2+kk1] - umean[k+1] ) ) )
            
                                  + cg3*( ci0*( ci0*( u[ijk    -kk2] - umean[k-2] ) + ci1*( u[ijk+jj1-kk2] - umean[k-2] ) + ci2*( u[ijk+jj2-kk2] - umean[k-2] ) + ci3*( u[ijk+jj3-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk    -kk1] - umean[k-1] ) + ci1*( u[ijk+jj1-kk1] - umean[k-1] ) + ci2*( u[ijk+jj2-kk1] - umean[k-1] ) + ci3*( u[ijk+jj3-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk        ] - umean[k  ] ) + ci1*( u[ijk+jj1    ] - umean[k  ] ) + ci2*( u[ijk+jj2    ] - umean[k  ] ) + ci3*( u[ijk+jj3    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk    +kk1] - umean[k+1] ) + ci1*( u[ijk+jj1+kk1] - umean[k+1] ) + ci2*( u[ijk+jj2+kk1] - umean[k+1] ) + ci3*( u[ijk+jj3+kk1] - umean[k+1] ) ) ) )
            
            
                                * cgi*dyi )
            
            
                              * ( cg0*( ci0*( ci0*w[ijk-ii2-jj3] + ci1*w[ijk-ii1-jj3] + ci2*w[ijk    -jj3] + ci3*w[ijk+ii1-jj3] )
                                      + ci1*( ci0*w[ijk-ii2-jj2] + ci1*w[ijk-ii1-jj2] + ci2*w[ijk    -jj2] + ci3*w[ijk+ii1-jj2] )
                                      + ci2*( ci0*w[ijk-ii2-jj1] + ci1*w[ijk-ii1-jj1] + ci2*w[ijk    -jj1] + ci3*w[ijk+ii1-jj1] )
                                      + ci3*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] ) )
            
                                + cg1*( ci0*( ci0*w[ijk-ii2-jj2] + ci1*w[ijk-ii1-jj2] + ci2*w[ijk    -jj2] + ci3*w[ijk+ii1-jj2] )
                                      + ci1*( ci0*w[ijk-ii2-jj1] + ci1*w[ijk-ii1-jj1] + ci2*w[ijk    -jj1] + ci3*w[ijk+ii1-jj1] )
                                      + ci2*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci3*( ci0*w[ijk-ii2+jj1] + ci1*w[ijk-ii1+jj1] + ci2*w[ijk    +jj1] + ci3*w[ijk+ii1+jj1] ) )
            
                                + cg2*( ci0*( ci0*w[ijk-ii2-jj1] + ci1*w[ijk-ii1-jj1] + ci2*w[ijk    -jj1] + ci3*w[ijk+ii1-jj1] )
                                      + ci1*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci2*( ci0*w[ijk-ii2+jj1] + ci1*w[ijk-ii1+jj1] + ci2*w[ijk    +jj1] + ci3*w[ijk+ii1+jj1] )
                                      + ci3*( ci0*w[ijk-ii2+jj2] + ci1*w[ijk-ii1+jj2] + ci2*w[ijk    +jj2] + ci3*w[ijk+ii1+jj2] ) )
            
                                + cg3*( ci0*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci1*( ci0*w[ijk-ii2+jj1] + ci1*w[ijk-ii1+jj1] + ci2*w[ijk    +jj1] + ci3*w[ijk+ii1+jj1] )
                                      + ci2*( ci0*w[ijk-ii2+jj2] + ci1*w[ijk-ii1+jj2] + ci2*w[ijk    +jj2] + ci3*w[ijk+ii1+jj2] )
                                      + ci3*( ci0*w[ijk-ii2+jj3] + ci1*w[ijk-ii1+jj3] + ci2*w[ijk    +jj3] + ci3*w[ijk+ii1+jj3] ) ) ) )
            
            
                            * cgi*dyi ) );
            
            uw_diss[k] -= ( ( 2 * visc )
            
            
                          * ( ( ( ( cg0*( u[ijk-kk2] - umean[k-2] ) + cg1*( u[ijk-kk1] - umean[k-1] ) + cg2*( u[ijk    ] - umean[k  ] ) + cg3*( u[ijk+kk1] - umean[k+1] ) ) * dzhi4[k] )
            
            
                              * ( cg0*( ci0*( ci0*w[ijk-ii2-kk3] + ci1*w[ijk-ii1-kk3] + ci2*w[ijk    -kk3] + ci3*w[ijk+ii1-kk3] )
                                      + ci1*( ci0*w[ijk-ii2-kk2] + ci1*w[ijk-ii1-kk2] + ci2*w[ijk    -kk2] + ci3*w[ijk+ii1-kk2] )
                                      + ci2*( ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk    -kk1] + ci3*w[ijk+ii1-kk1] )
                                      + ci3*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] ) )
            
                                + cg1*( ci0*( ci0*w[ijk-ii2-kk2] + ci1*w[ijk-ii1-kk2] + ci2*w[ijk    -kk2] + ci3*w[ijk+ii1-kk2] )
                                      + ci1*( ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk    -kk1] + ci3*w[ijk+ii1-kk1] )
                                      + ci2*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci3*( ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk    +kk1] + ci3*w[ijk+ii1+kk1] ) )
            
                                + cg2*( ci0*( ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk    -kk1] + ci3*w[ijk+ii1-kk1] )
                                      + ci1*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci2*( ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk    +kk1] + ci3*w[ijk+ii1+kk1] )
                                      + ci3*( ci0*w[ijk-ii2+kk2] + ci1*w[ijk-ii1+kk2] + ci2*w[ijk    +kk2] + ci3*w[ijk+ii1+kk2] ) )
            
                                + cg3*( ti0*( ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk    -kk1] + ci3*w[ijk+ii1-kk1] )
                                      + ti1*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ti2*( ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk    +kk1] + ci3*w[ijk+ii1+kk1] )
                                      + ti3*( ci0*w[ijk-ii2+kk2] + ci1*w[ijk-ii1+kk2] + ci2*w[ijk    +kk2] + ci3*w[ijk+ii1+kk2] ) ) ) )
            
            
                            * dzhi4[k] ) );
        }

    k = grid.kend;
    uw_diss[k] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
 
            uw_diss[k] -= ( ( 2 * visc )
            
            
                          * ( ( ( ( cg0*( ci0*( ci0*( u[ijk-ii3-kk2] - umean[k-2] ) + ci1*( u[ijk-ii2-kk2] - umean[k-2] ) + ci2*( u[ijk-ii1-kk2] - umean[k-2] ) + ci3*( u[ijk    -kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-ii3-kk1] - umean[k-1] ) + ci1*( u[ijk-ii2-kk1] - umean[k-1] ) + ci2*( u[ijk-ii1-kk1] - umean[k-1] ) + ci3*( u[ijk    -kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-ii3    ] - umean[k  ] ) + ci1*( u[ijk-ii2    ] - umean[k  ] ) + ci2*( u[ijk-ii1    ] - umean[k  ] ) + ci3*( u[ijk        ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-ii3+kk1] - umean[k+1] ) + ci1*( u[ijk-ii2+kk1] - umean[k+1] ) + ci2*( u[ijk-ii1+kk1] - umean[k+1] ) + ci3*( u[ijk    +kk1] - umean[k+1] ) ) )
            
                                  + cg1*( ci0*( ci0*( u[ijk-ii2-kk2] - umean[k-2] ) + ci1*( u[ijk-ii1-kk2] - umean[k-2] ) + ci2*( u[ijk    -kk2] - umean[k-2] ) + ci3*( u[ijk+ii1-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-ii2-kk1] - umean[k-1] ) + ci1*( u[ijk-ii1-kk1] - umean[k-1] ) + ci2*( u[ijk    -kk1] - umean[k-1] ) + ci3*( u[ijk+ii1-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-ii2    ] - umean[k  ] ) + ci1*( u[ijk-ii1    ] - umean[k  ] ) + ci2*( u[ijk        ] - umean[k  ] ) + ci3*( u[ijk+ii1    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-ii2+kk1] - umean[k+1] ) + ci1*( u[ijk-ii1+kk1] - umean[k+1] ) + ci2*( u[ijk    +kk1] - umean[k+1] ) + ci3*( u[ijk+ii1+kk1] - umean[k+1] ) ) )
            
                                  + cg2*( ci0*( ci0*( u[ijk-ii1-kk2] - umean[k-2] ) + ci1*( u[ijk    -kk2] - umean[k-2] ) + ci2*( u[ijk+ii1-kk2] - umean[k-2] ) + ci3*( u[ijk+ii2-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-ii1-kk1] - umean[k-1] ) + ci1*( u[ijk    -kk1] - umean[k-1] ) + ci2*( u[ijk+ii1-kk1] - umean[k-1] ) + ci3*( u[ijk+ii2-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-ii1    ] - umean[k  ] ) + ci1*( u[ijk        ] - umean[k  ] ) + ci2*( u[ijk+ii1    ] - umean[k  ] ) + ci3*( u[ijk+ii2    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-ii1+kk1] - umean[k+1] ) + ci1*( u[ijk    +kk1] - umean[k+1] ) + ci2*( u[ijk+ii1+kk1] - umean[k+1] ) + ci3*( u[ijk+ii2+kk1] - umean[k+1] ) ) )
            
                                  + cg3*( ci0*( ci0*( u[ijk    -kk2] - umean[k-2] ) + ci1*( u[ijk+ii1-kk2] - umean[k-2] ) + ci2*( u[ijk+ii2-kk2] - umean[k-2] ) + ci3*( u[ijk+ii3-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk    -kk1] - umean[k-1] ) + ci1*( u[ijk+ii1-kk1] - umean[k-1] ) + ci2*( u[ijk+ii2-kk1] - umean[k-1] ) + ci3*( u[ijk+ii3-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk        ] - umean[k  ] ) + ci1*( u[ijk+ii1    ] - umean[k  ] ) + ci2*( u[ijk+ii2    ] - umean[k  ] ) + ci3*( u[ijk+ii3    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk    +kk1] - umean[k+1] ) + ci1*( u[ijk+ii1+kk1] - umean[k+1] ) + ci2*( u[ijk+ii2+kk1] - umean[k+1] ) + ci3*( u[ijk+ii3+kk1] - umean[k+1] ) ) ) )
            
            
                                * cgi*dxi )
            
            
                              * ( cg0*w[ijk-ii2] + cg1*w[ijk-ii1] + cg2*w[ijk    ] + cg3*w[ijk+ii1] ) )
            
            
                            * cgi*dxi ) );
            
            uw_diss[k] -= ( ( 2 * visc )
            
            
                          * ( ( ( ( cg0*( ci0*( ci0*( u[ijk-jj3-kk2] - umean[k-2] ) + ci1*( u[ijk-jj2-kk2] - umean[k-2] ) + ci2*( u[ijk-jj1-kk2] - umean[k-2] ) + ci3*( u[ijk    -kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-jj3-kk1] - umean[k-1] ) + ci1*( u[ijk-jj2-kk1] - umean[k-1] ) + ci2*( u[ijk-jj1-kk1] - umean[k-1] ) + ci3*( u[ijk    -kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-jj3    ] - umean[k  ] ) + ci1*( u[ijk-jj2    ] - umean[k  ] ) + ci2*( u[ijk-jj1    ] - umean[k  ] ) + ci3*( u[ijk        ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-jj3+kk1] - umean[k+1] ) + ci1*( u[ijk-jj2+kk1] - umean[k+1] ) + ci2*( u[ijk-jj1+kk1] - umean[k+1] ) + ci3*( u[ijk    +kk1] - umean[k+1] ) ) )
            
                                  + cg1*( ci0*( ci0*( u[ijk-jj2-kk2] - umean[k-2] ) + ci1*( u[ijk-jj1-kk2] - umean[k-2] ) + ci2*( u[ijk    -kk2] - umean[k-2] ) + ci3*( u[ijk+jj1-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-jj2-kk1] - umean[k-1] ) + ci1*( u[ijk-jj1-kk1] - umean[k-1] ) + ci2*( u[ijk    -kk1] - umean[k-1] ) + ci3*( u[ijk+jj1-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-jj2    ] - umean[k  ] ) + ci1*( u[ijk-jj1    ] - umean[k  ] ) + ci2*( u[ijk        ] - umean[k  ] ) + ci3*( u[ijk+jj1    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-jj2+kk1] - umean[k+1] ) + ci1*( u[ijk-jj1+kk1] - umean[k+1] ) + ci2*( u[ijk    +kk1] - umean[k+1] ) + ci3*( u[ijk+jj1+kk1] - umean[k+1] ) ) )
            
                                  + cg2*( ci0*( ci0*( u[ijk-jj1-kk2] - umean[k-2] ) + ci1*( u[ijk    -kk2] - umean[k-2] ) + ci2*( u[ijk+jj1-kk2] - umean[k-2] ) + ci3*( u[ijk+jj2-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk-jj1-kk1] - umean[k-1] ) + ci1*( u[ijk    -kk1] - umean[k-1] ) + ci2*( u[ijk+jj1-kk1] - umean[k-1] ) + ci3*( u[ijk+jj2-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk-jj1    ] - umean[k  ] ) + ci1*( u[ijk        ] - umean[k  ] ) + ci2*( u[ijk+jj1    ] - umean[k  ] ) + ci3*( u[ijk+jj2    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk-jj1+kk1] - umean[k+1] ) + ci1*( u[ijk    +kk1] - umean[k+1] ) + ci2*( u[ijk+jj1+kk1] - umean[k+1] ) + ci3*( u[ijk+jj2+kk1] - umean[k+1] ) ) )
            
                                  + cg3*( ci0*( ci0*( u[ijk    -kk2] - umean[k-2] ) + ci1*( u[ijk+jj1-kk2] - umean[k-2] ) + ci2*( u[ijk+jj2-kk2] - umean[k-2] ) + ci3*( u[ijk+jj3-kk2] - umean[k-2] ) )
                                        + ci1*( ci0*( u[ijk    -kk1] - umean[k-1] ) + ci1*( u[ijk+jj1-kk1] - umean[k-1] ) + ci2*( u[ijk+jj2-kk1] - umean[k-1] ) + ci3*( u[ijk+jj3-kk1] - umean[k-1] ) )
                                        + ci2*( ci0*( u[ijk        ] - umean[k  ] ) + ci1*( u[ijk+jj1    ] - umean[k  ] ) + ci2*( u[ijk+jj2    ] - umean[k  ] ) + ci3*( u[ijk+jj3    ] - umean[k  ] ) )
                                        + ci3*( ci0*( u[ijk    +kk1] - umean[k+1] ) + ci1*( u[ijk+jj1+kk1] - umean[k+1] ) + ci2*( u[ijk+jj2+kk1] - umean[k+1] ) + ci3*( u[ijk+jj3+kk1] - umean[k+1] ) ) ) )
            
            
                                * cgi*dyi )
            
            
                              * ( cg0*( ci0*( ci0*w[ijk-ii2-jj3] + ci1*w[ijk-ii1-jj3] + ci2*w[ijk    -jj3] + ci3*w[ijk+ii1-jj3] )
                                      + ci1*( ci0*w[ijk-ii2-jj2] + ci1*w[ijk-ii1-jj2] + ci2*w[ijk    -jj2] + ci3*w[ijk+ii1-jj2] )
                                      + ci2*( ci0*w[ijk-ii2-jj1] + ci1*w[ijk-ii1-jj1] + ci2*w[ijk    -jj1] + ci3*w[ijk+ii1-jj1] )
                                      + ci3*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] ) )
            
                                + cg1*( ci0*( ci0*w[ijk-ii2-jj2] + ci1*w[ijk-ii1-jj2] + ci2*w[ijk    -jj2] + ci3*w[ijk+ii1-jj2] )
                                      + ci1*( ci0*w[ijk-ii2-jj1] + ci1*w[ijk-ii1-jj1] + ci2*w[ijk    -jj1] + ci3*w[ijk+ii1-jj1] )
                                      + ci2*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci3*( ci0*w[ijk-ii2+jj1] + ci1*w[ijk-ii1+jj1] + ci2*w[ijk    +jj1] + ci3*w[ijk+ii1+jj1] ) )
            
                                + cg2*( ci0*( ci0*w[ijk-ii2-jj1] + ci1*w[ijk-ii1-jj1] + ci2*w[ijk    -jj1] + ci3*w[ijk+ii1-jj1] )
                                      + ci1*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci2*( ci0*w[ijk-ii2+jj1] + ci1*w[ijk-ii1+jj1] + ci2*w[ijk    +jj1] + ci3*w[ijk+ii1+jj1] )
                                      + ci3*( ci0*w[ijk-ii2+jj2] + ci1*w[ijk-ii1+jj2] + ci2*w[ijk    +jj2] + ci3*w[ijk+ii1+jj2] ) )
            
                                + cg3*( ci0*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci1*( ci0*w[ijk-ii2+jj1] + ci1*w[ijk-ii1+jj1] + ci2*w[ijk    +jj1] + ci3*w[ijk+ii1+jj1] )
                                      + ci2*( ci0*w[ijk-ii2+jj2] + ci1*w[ijk-ii1+jj2] + ci2*w[ijk    +jj2] + ci3*w[ijk+ii1+jj2] )
                                      + ci3*( ci0*w[ijk-ii2+jj3] + ci1*w[ijk-ii1+jj3] + ci2*w[ijk    +jj3] + ci3*w[ijk+ii1+jj3] ) ) ) )
            
            
                            * cgi*dyi ) );
            
            uw_diss[k] -= ( ( 2 * visc )
            
            
                          * ( ( ( ( cg0*( u[ijk-kk2] - umean[k-2] ) + cg1*( u[ijk-kk1] - umean[k-1] ) + cg2*( u[ijk    ] - umean[k  ] ) + cg3*( u[ijk+kk1] - umean[k+1] ) ) * dzhi4[k] )
            
            
                              * ( tg0*( ci0*( ci0*w[ijk-ii2-kk4] + ci1*w[ijk-ii1-kk4] + ci2*w[ijk    -kk4] + ci3*w[ijk+ii1-kk4] )
                                      + ci1*( ci0*w[ijk-ii2-kk3] + ci1*w[ijk-ii1-kk3] + ci2*w[ijk    -kk3] + ci3*w[ijk+ii1-kk3] )
                                      + ci2*( ci0*w[ijk-ii2-kk2] + ci1*w[ijk-ii1-kk2] + ci2*w[ijk    -kk2] + ci3*w[ijk+ii1-kk2] )
                                      + ci3*( ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk    -kk1] + ci3*w[ijk+ii1-kk1] ) )
            
                                + tg1*( ci0*( ci0*w[ijk-ii2-kk3] + ci1*w[ijk-ii1-kk3] + ci2*w[ijk    -kk3] + ci3*w[ijk+ii1-kk3] )
                                      + ci1*( ci0*w[ijk-ii2-kk2] + ci1*w[ijk-ii1-kk2] + ci2*w[ijk    -kk2] + ci3*w[ijk+ii1-kk2] )
                                      + ci2*( ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk    -kk1] + ci3*w[ijk+ii1-kk1] )
                                      + ci3*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] ) )
            
                                + tg2*( ci0*( ci0*w[ijk-ii2-kk2] + ci1*w[ijk-ii1-kk2] + ci2*w[ijk    -kk2] + ci3*w[ijk+ii1-kk2] )
                                      + ci1*( ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk    -kk1] + ci3*w[ijk+ii1-kk1] )
                                      + ci2*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ci3*( ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk    +kk1] + ci3*w[ijk+ii1+kk1] ) )
            
                                + tg3*( ti0*( ci0*w[ijk-ii2-kk2] + ci1*w[ijk-ii1-kk2] + ci2*w[ijk    -kk2] + ci3*w[ijk+ii1-kk2] )
                                      + ti1*( ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk    -kk1] + ci3*w[ijk+ii1-kk1] )
                                      + ti2*( ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk        ] + ci3*w[ijk+ii1    ] )
                                      + ti3*( ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk    +kk1] + ci3*w[ijk+ii1+kk1] ) ) ) )
            
            
                            * dzhi4top ) ); 
        }

    master.sum(u2_diss , grid.kcells);
    master.sum(v2_diss , grid.kcells);
    master.sum(w2_diss , grid.kcells);
    master.sum(tke_diss, grid.kcells);
    master.sum(uw_diss , grid.kcells);

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        u2_diss [k] /= n;
        v2_diss [k] /= n;
        tke_diss[k] /= n;
    }

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        w2_diss [k] /= n;
        uw_diss [k] /= n;
    }

    // 7. CALCULATE THE PRESSURE REDISTRIBUTION TERM
    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        u2_rdstr [k] = 0.;
        v2_rdstr [k] = 0.;

        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                u2_rdstr [k] += 2.*(ci0*p[ijk-ii2] + ci1*p[ijk-ii1] + ci2*p[ijk] + ci3*p[ijk+ii1])
                              * ( cg0*((ci0*(u[ijk-ii3]-umean[k]) + ci1*(u[ijk-ii2]-umean[k]) + ci2*(u[ijk-ii1]-umean[k]) + ci3*(u[ijk    ]-umean[k])))
                                + cg1*((ci0*(u[ijk-ii2]-umean[k]) + ci1*(u[ijk-ii1]-umean[k]) + ci2*(u[ijk    ]-umean[k]) + ci3*(u[ijk+ii1]-umean[k])))
                                + cg2*((ci0*(u[ijk-ii1]-umean[k]) + ci1*(u[ijk    ]-umean[k]) + ci2*(u[ijk+ii1]-umean[k]) + ci3*(u[ijk+ii2]-umean[k])))
                                + cg3*((ci0*(u[ijk    ]-umean[k]) + ci1*(u[ijk+ii1]-umean[k]) + ci2*(u[ijk+ii2]-umean[k]) + ci3*(u[ijk+ii3]-umean[k]))) ) * cgi*dxi;
                v2_rdstr [k] += 2.*(ci0*p[ijk-jj2] + ci1*p[ijk-jj1] + ci2*p[ijk] + ci3*p[ijk+jj1])
                              * ( cg0*((ci0*(v[ijk-jj3]-vmean[k]) + ci1*(v[ijk-jj2]-vmean[k]) + ci2*(v[ijk-jj1]-vmean[k]) + ci3*(v[ijk    ]-vmean[k])))
                                + cg1*((ci0*(v[ijk-jj2]-vmean[k]) + ci1*(v[ijk-jj1]-vmean[k]) + ci2*(v[ijk    ]-vmean[k]) + ci3*(v[ijk+jj1]-vmean[k])))
                                + cg2*((ci0*(v[ijk-jj1]-vmean[k]) + ci1*(v[ijk    ]-vmean[k]) + ci2*(v[ijk+jj1]-vmean[k]) + ci3*(v[ijk+jj2]-vmean[k])))
                                + cg3*((ci0*(v[ijk    ]-vmean[k]) + ci1*(v[ijk+jj1]-vmean[k]) + ci2*(v[ijk+jj2]-vmean[k]) + ci3*(v[ijk+jj3]-vmean[k]))) ) * cgi*dyi;
            }
    }

    for (int k=grid.kstart+1; k<grid.kend; ++k)
    {
        w2_rdstr[k] = 0.;
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                w2_rdstr[k] += 2.*(ci0*p[ijk-kk2] + ci1*p[ijk-kk1] + ci2*p[ijk] + ci3*p[ijk+kk1])
                             * ( cg0*(ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ])
                               + cg1*(ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1])
                               + cg2*(ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])
                               + cg3*(ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3]) ) * dzhi4[k];
            }
    }

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        uw_rdstr[k] = 0.;
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                uw_rdstr[k] += ( ( ci0*( ci0*p[ijk-ii2-kk2] + ci1*p[ijk-ii1-kk2] + ci2*p[ijk-kk2] + ci3*p[ijk+ii1-kk2] )
                                 + ci1*( ci0*p[ijk-ii2-kk1] + ci1*p[ijk-ii1-kk1] + ci2*p[ijk-kk1] + ci3*p[ijk+ii1-kk1] )
                                 + ci2*( ci0*p[ijk-ii2    ] + ci1*p[ijk-ii1    ] + ci2*p[ijk    ] + ci3*p[ijk+ii1    ] )
                                 + ci3*( ci0*p[ijk-ii2+kk1] + ci1*p[ijk-ii1+kk1] + ci2*p[ijk+kk1] + ci3*p[ijk+ii1+kk1] ) )

                               * ( ( ( cg0*( u[ijk-kk2] - umean[k-2] ) + cg1*( u[ijk-kk1] - umean[k-1] ) + cg2*( u[ijk] - umean[k] ) + cg3*( u[ijk+kk1] - umean[k+1] ) ) * dzhi4[k] ) + ( ( cg0*w[ijk-ii2] + cg1*w[ijk-ii1] + cg2*w[ijk] + cg3*w[ijk+ii1] ) * cgi*dxi ) ) );
            }
    }

    master.sum(u2_rdstr, grid.kcells);
    master.sum(v2_rdstr, grid.kcells);
    master.sum(w2_rdstr, grid.kcells);
    master.sum(uw_rdstr, grid.kcells);

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        u2_rdstr[k] /= n;
        v2_rdstr[k] /= n;
    }

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        w2_rdstr[k] /= n;
        uw_rdstr[k] /= n;
    }
}

void Budget_4::calc_tke_budget_buoy(double* restrict u, double* restrict w, double* restrict b,
                                    double* restrict umean, double* restrict bmean,
                                    double* restrict w2_buoy, double* restrict tke_buoy, double* restrict uw_buoy)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int jj1 = 1*grid.icells;
    const int kk1 = 1*grid.ijcells;
    const int kk2 = 2*grid.ijcells;

    const double n = grid.itot*grid.jtot;

    // calculate the buoyancy term
    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        tke_buoy[k] = 0.;

        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                tke_buoy[k] += (ci0*w[ijk-kk1] + ci1*w[ijk] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])*( b[ijk] - bmean[k] );
            }
    }
    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        w2_buoy[k] = 0.;
        uw_buoy[k] = 0.;
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                w2_buoy[k] += 2.*(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk] + ci3*b[ijk+kk1])*w[ijk];

                uw_buoy[k] += ( ( ci0*( u[ijk-kk2] - umean[k-2] ) + ci1*( u[ijk-kk1] - umean[k-1] ) + ci2*( u[ijk    ] - umean[k  ] ) + ci3*( u[ijk+kk1] - umean[k+1] ) )

                              * ( ci0*( ci0*( b[ijk-ii2-kk2] - bmean[k-2] ) + ci1*( b[ijk-ii1-kk2] - bmean[k-2] ) + ci2*( b[ijk    -kk2] - bmean[k-2] ) + ci3*( b[ijk+ii1-kk2] - bmean[k-2] ) )
                                + ci1*( ci0*( b[ijk-ii2-kk1] - bmean[k-1] ) + ci1*( b[ijk-ii1-kk1] - bmean[k-1] ) + ci2*( b[ijk    -kk1] - bmean[k-1] ) + ci3*( b[ijk+ii1-kk1] - bmean[k-1] ) )
                                + ci2*( ci0*( b[ijk-ii2    ] - bmean[k  ] ) + ci1*( b[ijk-ii1    ] - bmean[k  ] ) + ci2*( b[ijk        ] - bmean[k  ] ) + ci3*( b[ijk+ii1    ] - bmean[k  ] ) )
                                + ci3*( ci0*( b[ijk-ii2+kk1] - bmean[k+1] ) + ci1*( b[ijk-ii1+kk1] - bmean[k+1] ) + ci2*( b[ijk    +kk1] - bmean[k+1] ) + ci3*( b[ijk+ii1+kk1] - bmean[k+1] ) ) ) );
            }
    }

    master.sum(tke_buoy, grid.kcells);
    master.sum(w2_buoy , grid.kcells);
    master.sum(uw_buoy , grid.kcells);

    for (int k=grid.kstart; k<grid.kend; ++k)
        tke_buoy[k] /= n;

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        w2_buoy[k] /= n;
        uw_buoy[k] /= n;
    }
}

void Budget_4::calc_b2_budget(double* restrict w, double* restrict b,
                              double* restrict bmean,
                              double* restrict b2_shear, double* restrict b2_turb, double* restrict b2_visc, double* restrict b2_diss,
                              double* restrict dzi4, double* restrict dzhi4,
                              const double visc)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int ii3 = 3;
    const int jj1 = 1*grid.icells;
    const int jj2 = 2*grid.icells;
    const int jj3 = 3*grid.icells;
    const int kk1 = 1*grid.ijcells;
    const int kk2 = 2*grid.ijcells;
    const int kk3 = 3*grid.ijcells;

    const double n = grid.itot*grid.jtot;

    const double dxi = grid.dxi;
    const double dyi = grid.dyi;

    // 1. CALCULATE THE GRADIENT PRODUCTION TERM
    int k = grid.kstart;
    b2_shear[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            b2_shear[k] -= ( ( ( ( 2.0 * ( b[ijk] - bmean[k] ) ) * ( ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2] ) )

                             * ( cg0*( bi0*bmean[k-2] + bi1*bmean[k-1] + bi2*bmean[k  ] + bi3*bmean[k+1] )
                               + cg1*( ci0*bmean[k-2] + ci1*bmean[k-1] + ci2*bmean[k  ] + ci3*bmean[k+1] )
                               + cg2*( ci0*bmean[k-1] + ci1*bmean[k  ] + ci2*bmean[k+1] + ci3*bmean[k+2] )
                               + cg3*( ci0*bmean[k  ] + ci1*bmean[k+1] + ci2*bmean[k+2] + ci3*bmean[k+3] ) ) )

                           * dzi4[k] );
        }

    for (int k=grid.kstart+1; k<grid.kend-1; ++k)
    {
        b2_shear[k] = 0;
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                b2_shear[k] -= ( ( ( ( 2.0 * ( b[ijk] - bmean[k] ) ) * ( ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2] ) )

                                 * ( cg0*( ci0*bmean[k-3] + ci1*bmean[k-2] + ci2*bmean[k-1] + ci3*bmean[k  ] )
                                   + cg1*( ci0*bmean[k-2] + ci1*bmean[k-1] + ci2*bmean[k  ] + ci3*bmean[k+1] )
                                   + cg2*( ci0*bmean[k-1] + ci1*bmean[k  ] + ci2*bmean[k+1] + ci3*bmean[k+2] )
                                   + cg3*( ci0*bmean[k  ] + ci1*bmean[k+1] + ci2*bmean[k+2] + ci3*bmean[k+3] ) ) )

                               * dzi4[k] );
            }
    }

    k = grid.kend-1;
    b2_shear[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            b2_shear[k] -= ( ( ( ( 2.0 * ( b[ijk] - bmean[k] ) ) * ( ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2] ) )

                             * ( cg0*( ci0*bmean[k-3] + ci1*bmean[k-2] + ci2*bmean[k-1] + ci3*bmean[k  ] )
                               + cg1*( ci0*bmean[k-2] + ci1*bmean[k-1] + ci2*bmean[k  ] + ci3*bmean[k+1] )
                               + cg2*( ci0*bmean[k-1] + ci1*bmean[k  ] + ci2*bmean[k+1] + ci3*bmean[k+2] )
                               + cg3*( ti0*bmean[k-1] + ti1*bmean[k  ] + ti2*bmean[k+1] + ti3*bmean[k+2] ) ) )

                           * dzi4[k] );
        }


    // 2. CALCULATE THE TURBULENT TRANSPORT TERM
    k = grid.kstart;
    b2_turb[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            b2_turb[k] -= ( ( cg0*( std::pow( ( bi0*( b[ijk-kk2] - bmean[k-2] ) + bi1*( b[ijk-kk1] - bmean[k-1] ) + bi2*( b[ijk    ] - bmean[k  ] ) + bi3*( b[ijk+kk1] - bmean[k+1] ) ), 2 ) * w[ijk-kk1] )
                            + cg1*( std::pow( ( ci0*( b[ijk-kk2] - bmean[k-2] ) + ci1*( b[ijk-kk1] - bmean[k-1] ) + ci2*( b[ijk    ] - bmean[k  ] ) + ci3*( b[ijk+kk1] - bmean[k+1] ) ), 2 ) * w[ijk    ] )
                            + cg2*( std::pow( ( ci0*( b[ijk-kk1] - bmean[k-1] ) + ci1*( b[ijk    ] - bmean[k  ] ) + ci2*( b[ijk+kk1] - bmean[k+1] ) + ci3*( b[ijk+kk2] - bmean[k+2] ) ), 2 ) * w[ijk+kk1] )
                            + cg3*( std::pow( ( ci0*( b[ijk    ] - bmean[k  ] ) + ci1*( b[ijk+kk1] - bmean[k+1] ) + ci2*( b[ijk+kk2] - bmean[k+2] ) + ci3*( b[ijk+kk3] - bmean[k+3] ) ), 2 ) * w[ijk+kk2] ) )

                          * dzi4[k] );
        }

    for (int k=grid.kstart+1; k<grid.kend-1; ++k)
    {
        b2_turb[k] = 0;
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                b2_turb[k] -= ( ( cg0*( std::pow( ( ci0*( b[ijk-kk3] - bmean[k-3] ) + ci1*( b[ijk-kk2] - bmean[k-2] ) + ci2*( b[ijk-kk1] - bmean[k-1] ) + ci3*( b[ijk    ] - bmean[k  ] ) ), 2 ) * w[ijk-kk1] )
                                + cg1*( std::pow( ( ci0*( b[ijk-kk2] - bmean[k-2] ) + ci1*( b[ijk-kk1] - bmean[k-1] ) + ci2*( b[ijk    ] - bmean[k  ] ) + ci3*( b[ijk+kk1] - bmean[k+1] ) ), 2 ) * w[ijk    ] )
                                + cg2*( std::pow( ( ci0*( b[ijk-kk1] - bmean[k-1] ) + ci1*( b[ijk    ] - bmean[k  ] ) + ci2*( b[ijk+kk1] - bmean[k+1] ) + ci3*( b[ijk+kk2] - bmean[k+2] ) ), 2 ) * w[ijk+kk1] )
                                + cg3*( std::pow( ( ci0*( b[ijk    ] - bmean[k  ] ) + ci1*( b[ijk+kk1] - bmean[k+1] ) + ci2*( b[ijk+kk2] - bmean[k+2] ) + ci3*( b[ijk+kk3] - bmean[k+3] ) ), 2 ) * w[ijk+kk2] ) )

                              * dzi4[k] );
            }
    }

    k = grid.kend-1;
    b2_turb[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            b2_turb[k] -= ( ( cg0*( std::pow( ( ci0*( b[ijk-kk3] - bmean[k-3] ) + ci1*( b[ijk-kk2] - bmean[k-2] ) + ci2*( b[ijk-kk1] - bmean[k-1] ) + ci3*( b[ijk    ] - bmean[k  ] ) ), 2 ) * w[ijk-kk1] )
                            + cg1*( std::pow( ( ci0*( b[ijk-kk2] - bmean[k-2] ) + ci1*( b[ijk-kk1] - bmean[k-1] ) + ci2*( b[ijk    ] - bmean[k  ] ) + ci3*( b[ijk+kk1] - bmean[k+1] ) ), 2 ) * w[ijk    ] )
                            + cg2*( std::pow( ( ci0*( b[ijk-kk1] - bmean[k-1] ) + ci1*( b[ijk    ] - bmean[k  ] ) + ci2*( b[ijk+kk1] - bmean[k+1] ) + ci3*( b[ijk+kk2] - bmean[k+2] ) ), 2 ) * w[ijk+kk1] )
                            + cg3*( std::pow( ( ti0*( b[ijk-kk1] - bmean[k-1] ) + ti1*( b[ijk    ] - bmean[k  ] ) + ti2*( b[ijk+kk1] - bmean[k+1] ) + ti3*( b[ijk+kk2] - bmean[k+2] ) ), 2 ) * w[ijk+kk2] ) )

                          * dzi4[k] );
        }


    // 3. CALCULATE THE VISCOUS TRANSPORT TERM
    k = grid.kstart;
    b2_visc[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            b2_visc[k] += ( ( visc

                            * ( cg0*( ( bg0*std::pow( ( b[ijk-kk2] - bmean[k-2] ), 2 ) + bg1*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + bg2*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + bg3*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) ) * dzhi4[k-1] )
                              + cg1*( ( cg0*std::pow( ( b[ijk-kk2] - bmean[k-2] ), 2 ) + cg1*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + cg2*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + cg3*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) ) * dzhi4[k  ] )
                              + cg2*( ( cg0*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + cg1*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + cg2*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) + cg3*std::pow( ( b[ijk+kk2] - bmean[k+2] ), 2 ) ) * dzhi4[k+1] )
                              + cg3*( ( cg0*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + cg1*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) + cg2*std::pow( ( b[ijk+kk2] - bmean[k+2] ), 2 ) + cg3*std::pow( ( b[ijk+kk3] - bmean[k+3] ), 2 ) ) * dzhi4[k+2] ) ) )

                          * dzi4[k] );
        }

    for (int k=grid.kstart+1; k<grid.kend-1; ++k)
    {
        b2_visc[k] = 0;
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                b2_visc[k] += ( ( visc

                                * ( cg0*( ( cg0*std::pow( ( b[ijk-kk3] - bmean[k-3] ), 2 ) + cg1*std::pow( ( b[ijk-kk2] - bmean[k-2] ), 2 ) + cg2*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + cg3*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) ) * dzhi4[k-1] )
                                  + cg1*( ( cg0*std::pow( ( b[ijk-kk2] - bmean[k-2] ), 2 ) + cg1*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + cg2*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + cg3*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) ) * dzhi4[k  ] )
                                  + cg2*( ( cg0*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + cg1*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + cg2*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) + cg3*std::pow( ( b[ijk+kk2] - bmean[k+2] ), 2 ) ) * dzhi4[k+1] )
                                  + cg3*( ( cg0*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + cg1*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) + cg2*std::pow( ( b[ijk+kk2] - bmean[k+2] ), 2 ) + cg3*std::pow( ( b[ijk+kk3] - bmean[k+3] ), 2 ) ) * dzhi4[k+2] ) ) )

                              * dzi4[k] );
            }
    }

    k = grid.kend-1;
    b2_visc[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            b2_visc[k] += ( ( visc

                            * ( cg0*( ( cg0*std::pow( ( b[ijk-kk3] - bmean[k-3] ), 2 ) + cg1*std::pow( ( b[ijk-kk2] - bmean[k-2] ), 2 ) + cg2*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + cg3*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) ) * dzhi4[k-1] )
                              + cg1*( ( cg0*std::pow( ( b[ijk-kk2] - bmean[k-2] ), 2 ) + cg1*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + cg2*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + cg3*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) ) * dzhi4[k  ] )
                              + cg2*( ( cg0*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + cg1*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + cg2*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) + cg3*std::pow( ( b[ijk+kk2] - bmean[k+2] ), 2 ) ) * dzhi4[k+1] )
                              + cg3*( ( tg0*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + tg1*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + tg2*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) + tg3*std::pow( ( b[ijk+kk2] - bmean[k+2] ), 2 ) ) * dzhi4[k+2] ) ) )

                          * dzi4[k] );
        }


    // 4. CALCULATE THE DISSIPATION TERM
    k = grid.kstart;
    b2_diss[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            b2_diss[k] -= ( ( 2.0 * visc )

                          * ( ( std::pow( ( ( cg0*( ci0*( b[ijk-ii3] - bmean[k] ) + ci1*( b[ijk-ii2] - bmean[k] ) + ci2*( b[ijk-ii1] - bmean[k] ) + ci3*( b[ijk    ] - bmean[k] ) )
                                            + cg1*( ci0*( b[ijk-ii2] - bmean[k] ) + ci1*( b[ijk-ii1] - bmean[k] ) + ci2*( b[ijk    ] - bmean[k] ) + ci3*( b[ijk+ii1] - bmean[k] ) )
                                            + cg2*( ci0*( b[ijk-ii1] - bmean[k] ) + ci1*( b[ijk    ] - bmean[k] ) + ci2*( b[ijk+ii1] - bmean[k] ) + ci3*( b[ijk+ii2] - bmean[k] ) )
                                            + cg3*( ci0*( b[ijk    ] - bmean[k] ) + ci1*( b[ijk+ii1] - bmean[k] ) + ci2*( b[ijk+ii2] - bmean[k] ) + ci3*( b[ijk+ii3] - bmean[k] ) ) )

                                          * cgi*dxi )

                                , 2 )

                              + std::pow( ( ( cg0*( ci0*( b[ijk-jj3] - bmean[k] ) + ci1*( b[ijk-jj2] - bmean[k] ) + ci2*( b[ijk-jj1] - bmean[k] ) + ci3*( b[ijk    ] - bmean[k] ) )
                                            + cg1*( ci0*( b[ijk-jj2] - bmean[k] ) + ci1*( b[ijk-jj1] - bmean[k] ) + ci2*( b[ijk    ] - bmean[k] ) + ci3*( b[ijk+jj1] - bmean[k] ) )
                                            + cg2*( ci0*( b[ijk-jj1] - bmean[k] ) + ci1*( b[ijk    ] - bmean[k] ) + ci2*( b[ijk+jj1] - bmean[k] ) + ci3*( b[ijk+jj2] - bmean[k] ) )
                                            + cg3*( ci0*( b[ijk    ] - bmean[k] ) + ci1*( b[ijk+jj1] - bmean[k] ) + ci2*( b[ijk+jj2] - bmean[k] ) + ci3*( b[ijk+jj3] - bmean[k] ) ) )

                                          * cgi*dyi )

                                , 2 ) )

                            + std::pow( ( ( cg0*( bi0*( b[ijk-kk2] - bmean[k-2] ) + bi1*( b[ijk-kk1] - bmean[k-1] ) + bi2*( b[ijk    ] - bmean[k  ] ) + bi3*( b[ijk+kk1] - bmean[k+1] ) )
                                          + cg1*( ci0*( b[ijk-kk2] - bmean[k-2] ) + ci1*( b[ijk-kk1] - bmean[k-1] ) + ci2*( b[ijk    ] - bmean[k  ] ) + ci3*( b[ijk+kk1] - bmean[k+1] ) )
                                          + cg2*( ci0*( b[ijk-kk1] - bmean[k-1] ) + ci1*( b[ijk    ] - bmean[k  ] ) + ci2*( b[ijk+kk1] - bmean[k+1] ) + ci3*( b[ijk+kk2] - bmean[k+2] ) )
                                          + cg3*( ci0*( b[ijk    ] - bmean[k  ] ) + ci1*( b[ijk+kk1] - bmean[k+1] ) + ci2*( b[ijk+kk2] - bmean[k+2] ) + ci3*( b[ijk+kk3] - bmean[k+3] ) ) )

                                        * dzi4[k] )

                              , 2 ) ) );
        }

    for (int k=grid.kstart+1; k<grid.kend-1; ++k)
    {
        b2_diss[k] = 0;
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                b2_diss[k] -= ( ( 2.0 * visc )

                              * ( ( std::pow( ( ( cg0*( ci0*( b[ijk-ii3] - bmean[k] ) + ci1*( b[ijk-ii2] - bmean[k] ) + ci2*( b[ijk-ii1] - bmean[k] ) + ci3*( b[ijk    ] - bmean[k] ) )
                                                + cg1*( ci0*( b[ijk-ii2] - bmean[k] ) + ci1*( b[ijk-ii1] - bmean[k] ) + ci2*( b[ijk    ] - bmean[k] ) + ci3*( b[ijk+ii1] - bmean[k] ) )
                                                + cg2*( ci0*( b[ijk-ii1] - bmean[k] ) + ci1*( b[ijk    ] - bmean[k] ) + ci2*( b[ijk+ii1] - bmean[k] ) + ci3*( b[ijk+ii2] - bmean[k] ) )
                                                + cg3*( ci0*( b[ijk    ] - bmean[k] ) + ci1*( b[ijk+ii1] - bmean[k] ) + ci2*( b[ijk+ii2] - bmean[k] ) + ci3*( b[ijk+ii3] - bmean[k] ) ) )

                                              * cgi*dxi )

                                    , 2 )

                                  + std::pow( ( ( cg0*( ci0*( b[ijk-jj3] - bmean[k] ) + ci1*( b[ijk-jj2] - bmean[k] ) + ci2*( b[ijk-jj1] - bmean[k] ) + ci3*( b[ijk    ] - bmean[k] ) )
                                                + cg1*( ci0*( b[ijk-jj2] - bmean[k] ) + ci1*( b[ijk-jj1] - bmean[k] ) + ci2*( b[ijk    ] - bmean[k] ) + ci3*( b[ijk+jj1] - bmean[k] ) )
                                                + cg2*( ci0*( b[ijk-jj1] - bmean[k] ) + ci1*( b[ijk    ] - bmean[k] ) + ci2*( b[ijk+jj1] - bmean[k] ) + ci3*( b[ijk+jj2] - bmean[k] ) )
                                                + cg3*( ci0*( b[ijk    ] - bmean[k] ) + ci1*( b[ijk+jj1] - bmean[k] ) + ci2*( b[ijk+jj2] - bmean[k] ) + ci3*( b[ijk+jj3] - bmean[k] ) ) )

                                              * cgi*dyi )

                                    , 2 ) )

                                + std::pow( ( ( cg0*( ci0*( b[ijk-kk3] - bmean[k-3] ) + ci1*( b[ijk-kk2] - bmean[k-2] ) + ci2*( b[ijk-kk1] - bmean[k-1] ) + ci3*( b[ijk    ] - bmean[k  ] ) )
                                              + cg1*( ci0*( b[ijk-kk2] - bmean[k-2] ) + ci1*( b[ijk-kk1] - bmean[k-1] ) + ci2*( b[ijk    ] - bmean[k  ] ) + ci3*( b[ijk+kk1] - bmean[k+1] ) )
                                              + cg2*( ci0*( b[ijk-kk1] - bmean[k-1] ) + ci1*( b[ijk    ] - bmean[k  ] ) + ci2*( b[ijk+kk1] - bmean[k+1] ) + ci3*( b[ijk+kk2] - bmean[k+2] ) )
                                              + cg3*( ci0*( b[ijk    ] - bmean[k  ] ) + ci1*( b[ijk+kk1] - bmean[k+1] ) + ci2*( b[ijk+kk2] - bmean[k+2] ) + ci3*( b[ijk+kk3] - bmean[k+3] ) ) )

                                            * dzi4[k] )

                                  , 2 ) ) );
            }
    }

    k = grid.kend-1;
    b2_diss[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            b2_diss[k] -= ( ( 2.0 * visc )

                          * ( ( std::pow( ( ( cg0*( ci0*( b[ijk-ii3] - bmean[k] ) + ci1*( b[ijk-ii2] - bmean[k] ) + ci2*( b[ijk-ii1] - bmean[k] ) + ci3*( b[ijk    ] - bmean[k] ) )
                                            + cg1*( ci0*( b[ijk-ii2] - bmean[k] ) + ci1*( b[ijk-ii1] - bmean[k] ) + ci2*( b[ijk    ] - bmean[k] ) + ci3*( b[ijk+ii1] - bmean[k] ) )
                                            + cg2*( ci0*( b[ijk-ii1] - bmean[k] ) + ci1*( b[ijk    ] - bmean[k] ) + ci2*( b[ijk+ii1] - bmean[k] ) + ci3*( b[ijk+ii2] - bmean[k] ) )
                                            + cg3*( ci0*( b[ijk    ] - bmean[k] ) + ci1*( b[ijk+ii1] - bmean[k] ) + ci2*( b[ijk+ii2] - bmean[k] ) + ci3*( b[ijk+ii3] - bmean[k] ) ) )

                                          * cgi*dxi )

                                , 2 )

                              + std::pow( ( ( cg0*( ci0*( b[ijk-jj3] - bmean[k] ) + ci1*( b[ijk-jj2] - bmean[k] ) + ci2*( b[ijk-jj1] - bmean[k] ) + ci3*( b[ijk    ] - bmean[k] ) )
                                            + cg1*( ci0*( b[ijk-jj2] - bmean[k] ) + ci1*( b[ijk-jj1] - bmean[k] ) + ci2*( b[ijk    ] - bmean[k] ) + ci3*( b[ijk+jj1] - bmean[k] ) )
                                            + cg2*( ci0*( b[ijk-jj1] - bmean[k] ) + ci1*( b[ijk    ] - bmean[k] ) + ci2*( b[ijk+jj1] - bmean[k] ) + ci3*( b[ijk+jj2] - bmean[k] ) )
                                            + cg3*( ci0*( b[ijk    ] - bmean[k] ) + ci1*( b[ijk+jj1] - bmean[k] ) + ci2*( b[ijk+jj2] - bmean[k] ) + ci3*( b[ijk+jj3] - bmean[k] ) ) )

                                          * cgi*dyi )

                                , 2 ) )

                            + std::pow( ( ( cg0*( ci0*( b[ijk-kk3] - bmean[k-3] ) + ci1*( b[ijk-kk2] - bmean[k-2] ) + ci2*( b[ijk-kk1] - bmean[k-1] ) + ci3*( b[ijk    ] - bmean[k  ] ) )
                                          + cg1*( ci0*( b[ijk-kk2] - bmean[k-2] ) + ci1*( b[ijk-kk1] - bmean[k-1] ) + ci2*( b[ijk    ] - bmean[k  ] ) + ci3*( b[ijk+kk1] - bmean[k+1] ) )
                                          + cg2*( ci0*( b[ijk-kk1] - bmean[k-1] ) + ci1*( b[ijk    ] - bmean[k  ] ) + ci2*( b[ijk+kk1] - bmean[k+1] ) + ci3*( b[ijk+kk2] - bmean[k+2] ) )
                                          + cg3*( ti0*( b[ijk-kk1] - bmean[k-1] ) + ti1*( b[ijk    ] - bmean[k  ] ) + ti2*( b[ijk+kk1] - bmean[k+1] ) + ti3*( b[ijk+kk2] - bmean[k+2] ) ) )

                                        * dzi4[k] )

                              , 2 ) ) );
        }

    master.sum(b2_shear, grid.kcells);
    master.sum(b2_turb , grid.kcells);
    master.sum(b2_visc , grid.kcells);
    master.sum(b2_diss , grid.kcells);

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        b2_shear[k] /= n;
        b2_turb [k] /= n;
        b2_visc [k] /= n;
        b2_diss [k] /= n;
    }

}

void Budget_4::calc_bw_budget(double* restrict w, double* restrict p, double* restrict b, double* restrict bz,
                              double* restrict pmean, double* restrict bmean,
                              double* restrict bw_shear, double* restrict bw_turb, double* restrict bw_visc,
                              double* restrict bw_buoy, double* restrict bw_rdstr, double* restrict bw_diss, double* restrict bw_pres,
                              double* restrict dzi4, double* restrict dzhi4,
                              const double visc)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int ii3 = 3;
    const int jj1 = 1*grid.icells;
    const int jj2 = 2*grid.icells;
    const int jj3 = 3*grid.icells;
    const int kk1 = 1*grid.ijcells;
    const int kk2 = 2*grid.ijcells;
    const int kk3 = 3*grid.ijcells;
    const int kk4 = 4*grid.ijcells;

    const double n = grid.itot*grid.jtot;

    const double dxi = grid.dxi;
    const double dyi = grid.dyi;

    const double dzhi4bot = grid.dzhi4bot;
    const double dzhi4top = grid.dzhi4top;

    // 0. Create an interpolated field for b on the cell face.
    int k = grid.kstart-1;
    for (int j=0; j<grid.jcells; ++j)
        #pragma ivdep
        for (int i=0; i<grid.icells; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            bz[ijk] = ( bi0*( b[ijk-kk1] - bmean[k-1] ) + bi1*( b[ijk    ] - bmean[k  ] ) + bi2*( b[ijk+kk1] - bmean[k+1] ) + bi3*( b[ijk+kk2] - bmean[k+2] ) );
        }

    for (int k=grid.kstart; k<grid.kend+1; ++k)
        for (int j=0; j<grid.jcells; ++j)
            #pragma ivdep
            for (int i=0; i<grid.icells; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bz[ijk] = ( ci0*( b[ijk-kk2] - bmean[k-2] ) + ci1*( b[ijk-kk1] - bmean[k-1] ) + ci2*( b[ijk    ] - bmean[k  ] ) + ci3*( b[ijk+kk1] - bmean[k+1] ) );
            }

    k = grid.kend+1;
    for (int j=0; j<grid.jcells; ++j)
        #pragma ivdep
        for (int i=0; i<grid.icells; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            bz[ijk] = ( ti0*( b[ijk-kk3] - bmean[k-3] ) + ti1*( b[ijk-kk2] - bmean[k-2] ) + ti2*( b[ijk-kk1] - bmean[k-1] ) + ti3*( b[ijk    ] - bmean[k  ] ) );
        }


    // 1. CALCULATE THE GRADIENT PRODUCTION TERM
    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        bw_shear[k] = 0;
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_shear[k] -= ( ( std::pow( w[ijk], 2 ) * ( cg0*bmean[k-2] + cg1*bmean[k-1] + cg2*bmean[k  ] + cg3*bmean[k+1] ) ) * dzhi4[k] );
            }
    }


    // 2. CALCULATE THE TURBULENT TRANSPORT TERM
    k = grid.kstart;
    bw_turb[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            bw_turb[k] -= ( ( bg0*( std::pow( ( bi0*w[ijk-kk1] + bi1*w[ijk    ] + bi2*w[ijk+kk1] + bi3*w[ijk+kk2] ), 2 ) * ( b[ijk-kk1] - bmean[k-1] ) )
                            + bg1*( std::pow( ( ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2] ), 2 ) * ( b[ijk    ] - bmean[k  ] ) )
                            + bg2*( std::pow( ( ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3] ), 2 ) * ( b[ijk+kk1] - bmean[k+1] ) )
                            + bg3*( std::pow( ( ci0*w[ijk+kk1] + ci1*w[ijk+kk2] + ci2*w[ijk+kk3] + ci3*w[ijk+kk4] ), 2 ) * ( b[ijk+kk2] - bmean[k+2] ) ) )

                          * dzhi4bot );
        }

    k = grid.kstart+1;
    bw_turb[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            bw_turb[k] -= ( ( cg0*( std::pow( ( bi0*w[ijk-kk2] + bi1*w[ijk-kk1] + bi2*w[ijk    ] + bi3*w[ijk+kk1] ), 2 ) * ( b[ijk-kk2] - bmean[k-2] ) )
                            + cg1*( std::pow( ( ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1] ), 2 ) * ( b[ijk-kk1] - bmean[k-1] ) )
                            + cg2*( std::pow( ( ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2] ), 2 ) * ( b[ijk    ] - bmean[k  ] ) )
                            + cg3*( std::pow( ( ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3] ), 2 ) * ( b[ijk+kk1] - bmean[k+1] ) ) )

                          * dzhi4[k] );
        }

    for (int k=grid.kstart+2; k<grid.kend-1; ++k)
    {
        bw_turb[k] = 0;
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_turb[k] -= ( ( cg0*( std::pow( ( ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ] ), 2 ) * ( b[ijk-kk2] - bmean[k-2] ) )
                                + cg1*( std::pow( ( ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1] ), 2 ) * ( b[ijk-kk1] - bmean[k-1] ) )
                                + cg2*( std::pow( ( ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2] ), 2 ) * ( b[ijk    ] - bmean[k  ] ) )
                                + cg3*( std::pow( ( ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3] ), 2 ) * ( b[ijk+kk1] - bmean[k+1] ) ) )

                              * dzhi4[k] );
            }
    }

    k = grid.kend-1;
    bw_turb[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            bw_turb[k] -= ( ( cg0*( std::pow( ( ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ] ), 2 ) * ( b[ijk-kk2] - bmean[k-2] ) )
                            + cg1*( std::pow( ( ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1] ), 2 ) * ( b[ijk-kk1] - bmean[k-1] ) )
                            + cg2*( std::pow( ( ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2] ), 2 ) * ( b[ijk    ] - bmean[k  ] ) )
                            + cg3*( std::pow( ( ti0*w[ijk-kk1] + ti1*w[ijk    ] + ti2*w[ijk+kk1] + ti3*w[ijk+kk2] ), 2 ) * ( b[ijk+kk1] - bmean[k+1] ) ) )

                          * dzhi4[k] );
        }

    k = grid.kend;
    bw_turb[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            bw_turb[k] -= ( ( tg0*( std::pow( ( ci0*w[ijk-kk4] + ci1*w[ijk-kk3] + ci2*w[ijk-kk2] + ci3*w[ijk-kk1] ), 2 ) * ( b[ijk-kk3] - bmean[k-3] ) )
                            + tg1*( std::pow( ( ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ] ), 2 ) * ( b[ijk-kk2] - bmean[k-2] ) )
                            + tg2*( std::pow( ( ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1] ), 2 ) * ( b[ijk-kk1] - bmean[k-1] ) )
                            + tg3*( std::pow( ( ti0*w[ijk-kk2] + ti1*w[ijk-kk1] + ti2*w[ijk    ] + ti3*w[ijk+kk1] ), 2 ) * ( b[ijk    ] - bmean[k  ] ) ) )

                          * dzhi4top );
        }


    // 3. CALCULATE THE VISCOUS TRANSPORT TERM
    k = grid.kstart;
    bw_visc[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            bw_visc[k] += ( ( visc
            
                            * ( bg0*( ( bg0*( w[ijk-kk1] * bz[ijk-kk1] ) + bg1*( w[ijk    ] * bz[ijk    ] ) + bg2*( w[ijk+kk1] * bz[ijk+kk1] ) + bg3*( w[ijk+kk2] * bz[ijk+kk2] ) ) * dzi4[k-1] )
                              + bg1*( ( cg0*( w[ijk-kk1] * bz[ijk-kk1] ) + cg1*( w[ijk    ] * bz[ijk    ] ) + cg2*( w[ijk+kk1] * bz[ijk+kk1] ) + cg3*( w[ijk+kk2] * bz[ijk+kk2] ) ) * dzi4[k  ] )
                              + bg2*( ( cg0*( w[ijk    ] * bz[ijk    ] ) + cg1*( w[ijk+kk1] * bz[ijk+kk1] ) + cg2*( w[ijk+kk2] * bz[ijk+kk2] ) + cg3*( w[ijk+kk3] * bz[ijk+kk3] ) ) * dzi4[k+1] )
                              + bg3*( ( cg0*( w[ijk+kk1] * bz[ijk+kk1] ) + cg1*( w[ijk+kk2] * bz[ijk+kk2] ) + cg2*( w[ijk+kk3] * bz[ijk+kk3] ) + cg3*( w[ijk+kk4] * bz[ijk+kk4] ) ) * dzi4[k+2] ) ) )
            
                          * dzhi4bot );
        }

    k = grid.kstart+1;
    bw_visc[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            bw_visc[k] += ( ( visc
            
                            * ( cg0*( ( bg0*( w[ijk-kk2] * bz[ijk-kk2] ) + bg1*( w[ijk-kk1] * bz[ijk-kk1] ) + bg2*( w[ijk    ] * bz[ijk    ] ) + bg3*( w[ijk+kk1] * bz[ijk+kk1] ) ) * dzi4[k-2] )
                              + cg1*( ( cg0*( w[ijk-kk2] * bz[ijk-kk2] ) + cg1*( w[ijk-kk1] * bz[ijk-kk1] ) + cg2*( w[ijk    ] * bz[ijk    ] ) + cg3*( w[ijk+kk1] * bz[ijk+kk1] ) ) * dzi4[k-1] )
                              + cg2*( ( cg0*( w[ijk-kk1] * bz[ijk-kk1] ) + cg1*( w[ijk    ] * bz[ijk    ] ) + cg2*( w[ijk+kk1] * bz[ijk+kk1] ) + cg3*( w[ijk+kk2] * bz[ijk+kk2] ) ) * dzi4[k  ] )
                              + cg3*( ( cg0*( w[ijk    ] * bz[ijk    ] ) + cg1*( w[ijk+kk1] * bz[ijk+kk1] ) + cg2*( w[ijk+kk2] * bz[ijk+kk2] ) + cg3*( w[ijk+kk3] * bz[ijk+kk3] ) ) * dzi4[k+1] ) ) )
            
                          * dzhi4[k] );
        }

    for (int k=grid.kstart+2; k<grid.kend-1; ++k)
    {
        bw_visc[k] = 0;
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_visc[k] += ( ( visc
                
                                * ( cg0*( ( cg0*( w[ijk-kk3] * bz[ijk-kk3] ) + cg1*( w[ijk-kk2] * bz[ijk-kk2] ) + cg2*( w[ijk-kk1] * bz[ijk-kk1] ) + cg3*( w[ijk    ] * bz[ijk    ] ) ) * dzi4[k-2] )
                                  + cg1*( ( cg0*( w[ijk-kk2] * bz[ijk-kk2] ) + cg1*( w[ijk-kk1] * bz[ijk-kk1] ) + cg2*( w[ijk    ] * bz[ijk    ] ) + cg3*( w[ijk+kk1] * bz[ijk+kk1] ) ) * dzi4[k-1] )
                                  + cg2*( ( cg0*( w[ijk-kk1] * bz[ijk-kk1] ) + cg1*( w[ijk    ] * bz[ijk    ] ) + cg2*( w[ijk+kk1] * bz[ijk+kk1] ) + cg3*( w[ijk+kk2] * bz[ijk+kk2] ) ) * dzi4[k  ] )
                                  + cg3*( ( cg0*( w[ijk    ] * bz[ijk    ] ) + cg1*( w[ijk+kk1] * bz[ijk+kk1] ) + cg2*( w[ijk+kk2] * bz[ijk+kk2] ) + cg3*( w[ijk+kk3] * bz[ijk+kk3] ) ) * dzi4[k+1] ) ) )
                
                              * dzhi4[k] );
            }
    }

    k = grid.kend-1;
    bw_visc[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            bw_visc[k] += ( ( visc
            
                            * ( cg0*( ( cg0*( w[ijk-kk3] * bz[ijk-kk3] ) + cg1*( w[ijk-kk2] * bz[ijk-kk2] ) + cg2*( w[ijk-kk1] * bz[ijk-kk1] ) + cg3*( w[ijk    ] * bz[ijk    ] ) ) * dzi4[k-2] )
                              + cg1*( ( cg0*( w[ijk-kk2] * bz[ijk-kk2] ) + cg1*( w[ijk-kk1] * bz[ijk-kk1] ) + cg2*( w[ijk    ] * bz[ijk    ] ) + cg3*( w[ijk+kk1] * bz[ijk+kk1] ) ) * dzi4[k-1] )
                              + cg2*( ( cg0*( w[ijk-kk1] * bz[ijk-kk1] ) + cg1*( w[ijk    ] * bz[ijk    ] ) + cg2*( w[ijk+kk1] * bz[ijk+kk1] ) + cg3*( w[ijk+kk2] * bz[ijk+kk2] ) ) * dzi4[k  ] )
                              + cg3*( ( tg0*( w[ijk-kk1] * bz[ijk-kk1] ) + tg1*( w[ijk    ] * bz[ijk    ] ) + tg2*( w[ijk+kk1] * bz[ijk+kk1] ) + tg3*( w[ijk+kk2] * bz[ijk+kk2] ) ) * dzi4[k+1] ) ) )
            
                          * dzhi4[k] );
        }

    k = grid.kend;
    bw_visc[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            bw_visc[k] += ( ( visc
            
                            * ( tg0*( ( cg0*( w[ijk-kk4] * bz[ijk-kk4] ) + cg1*( w[ijk-kk3] * bz[ijk-kk3] ) + cg2*( w[ijk-kk2] * bz[ijk-kk2] ) + cg3*( w[ijk-kk1] * bz[ijk-kk1] ) ) * dzi4[k-3] )
                              + tg1*( ( cg0*( w[ijk-kk3] * bz[ijk-kk3] ) + cg1*( w[ijk-kk2] * bz[ijk-kk2] ) + cg2*( w[ijk-kk1] * bz[ijk-kk1] ) + cg3*( w[ijk    ] * bz[ijk    ] ) ) * dzi4[k-2] )
                              + tg2*( ( cg0*( w[ijk-kk2] * bz[ijk-kk2] ) + cg1*( w[ijk-kk1] * bz[ijk-kk1] ) + cg2*( w[ijk    ] * bz[ijk    ] ) + cg3*( w[ijk+kk1] * bz[ijk+kk1] ) ) * dzi4[k-1] )
                              + tg3*( ( tg0*( w[ijk-kk2] * bz[ijk-kk2] ) + tg1*( w[ijk-kk1] * bz[ijk-kk1] ) + tg2*( w[ijk    ] * bz[ijk    ] ) + tg3*( w[ijk+kk1] * bz[ijk+kk1] ) ) * dzi4[k  ] ) ) )
            
                          * dzhi4top );
        }


    // 4. CALCULATE THE BUOYANCY TERM
    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        bw_buoy[k] = 0;
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_buoy[k] += std::pow( bz[ijk], 2 );
            }
    }

    // 5. CALCULATE THE REDISTRIBUTION TERM
    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        bw_rdstr[k] = 0;
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_rdstr[k] += ( ( ( ci0*( p[ijk-kk2] - pmean[k-2] ) + ci1*( p[ijk-kk1] - pmean[k-1] ) + ci2*( p[ijk    ] - pmean[k  ] ) + ci3*( p[ijk+kk1] - pmean[k+1] ) ) * ( cg0*( b[ijk-kk2] - bmean[k-2] ) + cg1*( b[ijk-kk1] - bmean[k-1] ) + cg2*( b[ijk    ] - bmean[k  ] ) + cg3*( b[ijk+kk1] - bmean[k+1] ) ) ) * dzhi4[k] );
            }
    }


    // 6. CALCULATE THE DISSIPATION TERM
    k = grid.kstart;
    bw_diss[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            bw_diss[k] -= ( ( 2.0 * visc )
            
                          * ( ( ( ( ( ( cg0*( ci0*w[ijk-ii3] + ci1*w[ijk-ii2] + ci2*w[ijk-ii1] + ci3*w[ijk    ] )
                                      + cg1*( ci0*w[ijk-ii2] + ci1*w[ijk-ii1] + ci2*w[ijk    ] + ci3*w[ijk+ii1] )
                                      + cg2*( ci0*w[ijk-ii1] + ci1*w[ijk    ] + ci2*w[ijk+ii1] + ci3*w[ijk+ii2] )
                                      + cg3*( ci0*w[ijk    ] + ci1*w[ijk+ii1] + ci2*w[ijk+ii2] + ci3*w[ijk+ii3] ) )
            
                                    * cgi*dxi )
            
                                  * ( cg0*( ci0*bz[ijk-ii3] + ci1*bz[ijk-ii2] + ci2*bz[ijk-ii1] + ci3*bz[ijk    ] )
                                    + cg1*( ci0*bz[ijk-ii2] + ci1*bz[ijk-ii1] + ci2*bz[ijk    ] + ci3*bz[ijk+ii1] )
                                    + cg2*( ci0*bz[ijk-ii1] + ci1*bz[ijk    ] + ci2*bz[ijk+ii1] + ci3*bz[ijk+ii2] )
                                    + cg3*( ci0*bz[ijk    ] + ci1*bz[ijk+ii1] + ci2*bz[ijk+ii2] + ci3*bz[ijk+ii3] ) ) )
            
                                * cgi*dxi )
            
                              + ( ( ( ( cg0*( ci0*w[ijk-jj3] + ci1*w[ijk-jj2] + ci2*w[ijk-jj1] + ci3*w[ijk    ] )
                                      + cg1*( ci0*w[ijk-jj2] + ci1*w[ijk-jj1] + ci2*w[ijk    ] + ci3*w[ijk+jj1] )
                                      + cg2*( ci0*w[ijk-jj1] + ci1*w[ijk    ] + ci2*w[ijk+jj1] + ci3*w[ijk+jj2] )
                                      + cg3*( ci0*w[ijk    ] + ci1*w[ijk+jj1] + ci2*w[ijk+jj2] + ci3*w[ijk+jj3] ) )
            
                                    * cgi*dyi )
            
                                  * ( cg0*( ci0*bz[ijk-jj3] + ci1*bz[ijk-jj2] + ci2*bz[ijk-jj1] + ci3*bz[ijk    ] )
                                    + cg1*( ci0*bz[ijk-jj2] + ci1*bz[ijk-jj1] + ci2*bz[ijk    ] + ci3*bz[ijk+jj1] )
                                    + cg2*( ci0*bz[ijk-jj1] + ci1*bz[ijk    ] + ci2*bz[ijk+jj1] + ci3*bz[ijk+jj2] )
                                    + cg3*( ci0*bz[ijk    ] + ci1*bz[ijk+jj1] + ci2*bz[ijk+jj2] + ci3*bz[ijk+jj3] ) ) )
            
                                * cgi*dyi ) )
            
                            + ( ( ( ( bg0*( bi0*w[ijk-kk1] + bi1*w[ijk    ] + bi2*w[ijk+kk1] + bi3*w[ijk+kk2] )
                                    + bg1*( ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2] )
                                    + bg2*( ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3] )
                                    + bg3*( ci0*w[ijk+kk1] + ci1*w[ijk+kk2] + ci2*w[ijk+kk3] + ci3*w[ijk+kk4] ) )
            
                                  * dzhi4bot )
            
                                * ( cg0*( b[ijk-kk2] - bmean[k-2] ) + cg1*( b[ijk-kk1] - bmean[k-1] ) + cg2*( b[ijk    ] - bmean[k  ] ) + cg3*( b[ijk+kk1] - bmean[k+1] ) ) )
            
                              * dzhi4bot ) ) );
        }

    k = grid.kstart+1;
    bw_diss[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            bw_diss[k] -= ( ( 2.0 * visc )
            
                          * ( ( ( ( ( ( cg0*( ci0*w[ijk-ii3] + ci1*w[ijk-ii2] + ci2*w[ijk-ii1] + ci3*w[ijk    ] )
                                      + cg1*( ci0*w[ijk-ii2] + ci1*w[ijk-ii1] + ci2*w[ijk    ] + ci3*w[ijk+ii1] )
                                      + cg2*( ci0*w[ijk-ii1] + ci1*w[ijk    ] + ci2*w[ijk+ii1] + ci3*w[ijk+ii2] )
                                      + cg3*( ci0*w[ijk    ] + ci1*w[ijk+ii1] + ci2*w[ijk+ii2] + ci3*w[ijk+ii3] ) )
            
                                    * cgi*dxi )
            
                                  * ( cg0*( ci0*bz[ijk-ii3] + ci1*bz[ijk-ii2] + ci2*bz[ijk-ii1] + ci3*bz[ijk    ] )
                                    + cg1*( ci0*bz[ijk-ii2] + ci1*bz[ijk-ii1] + ci2*bz[ijk    ] + ci3*bz[ijk+ii1] )
                                    + cg2*( ci0*bz[ijk-ii1] + ci1*bz[ijk    ] + ci2*bz[ijk+ii1] + ci3*bz[ijk+ii2] )
                                    + cg3*( ci0*bz[ijk    ] + ci1*bz[ijk+ii1] + ci2*bz[ijk+ii2] + ci3*bz[ijk+ii3] ) ) )
            
                                * cgi*dxi )
            
                              + ( ( ( ( cg0*( ci0*w[ijk-jj3] + ci1*w[ijk-jj2] + ci2*w[ijk-jj1] + ci3*w[ijk    ] )
                                      + cg1*( ci0*w[ijk-jj2] + ci1*w[ijk-jj1] + ci2*w[ijk    ] + ci3*w[ijk+jj1] )
                                      + cg2*( ci0*w[ijk-jj1] + ci1*w[ijk    ] + ci2*w[ijk+jj1] + ci3*w[ijk+jj2] )
                                      + cg3*( ci0*w[ijk    ] + ci1*w[ijk+jj1] + ci2*w[ijk+jj2] + ci3*w[ijk+jj3] ) )
            
                                    * cgi*dyi )
            
                                  * ( cg0*( ci0*bz[ijk-jj3] + ci1*bz[ijk-jj2] + ci2*bz[ijk-jj1] + ci3*bz[ijk    ] )
                                    + cg1*( ci0*bz[ijk-jj2] + ci1*bz[ijk-jj1] + ci2*bz[ijk    ] + ci3*bz[ijk+jj1] )
                                    + cg2*( ci0*bz[ijk-jj1] + ci1*bz[ijk    ] + ci2*bz[ijk+jj1] + ci3*bz[ijk+jj2] )
                                    + cg3*( ci0*bz[ijk    ] + ci1*bz[ijk+jj1] + ci2*bz[ijk+jj2] + ci3*bz[ijk+jj3] ) ) )
            
                                * cgi*dyi ) )
            
                            + ( ( ( ( cg0*( bi0*w[ijk-kk2] + bi1*w[ijk-kk1] + bi2*w[ijk    ] + bi3*w[ijk+kk1] )
                                    + cg1*( ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1] )
                                    + cg2*( ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2] )
                                    + cg3*( ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3] ) )
            
                                  * dzhi4[k] )
            
                                * ( cg0*( b[ijk-kk2] - bmean[k-2] ) + cg1*( b[ijk-kk1] - bmean[k-1] ) + cg2*( b[ijk    ] - bmean[k  ] ) + cg3*( b[ijk+kk1] - bmean[k+1] ) ) )
            
                              * dzhi4[k] ) ) );
        }

    for (int k=grid.kstart+2; k<grid.kend-1; ++k)
    {
        bw_diss[k] = 0;
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_diss[k] -= ( ( 2.0 * visc )
                
                              * ( ( ( ( ( ( cg0*( ci0*w[ijk-ii3] + ci1*w[ijk-ii2] + ci2*w[ijk-ii1] + ci3*w[ijk    ] )
                                          + cg1*( ci0*w[ijk-ii2] + ci1*w[ijk-ii1] + ci2*w[ijk    ] + ci3*w[ijk+ii1] )
                                          + cg2*( ci0*w[ijk-ii1] + ci1*w[ijk    ] + ci2*w[ijk+ii1] + ci3*w[ijk+ii2] )
                                          + cg3*( ci0*w[ijk    ] + ci1*w[ijk+ii1] + ci2*w[ijk+ii2] + ci3*w[ijk+ii3] ) )
                
                                        * cgi*dxi )
                
                                      * ( cg0*( ci0*bz[ijk-ii3] + ci1*bz[ijk-ii2] + ci2*bz[ijk-ii1] + ci3*bz[ijk    ] )
                                        + cg1*( ci0*bz[ijk-ii2] + ci1*bz[ijk-ii1] + ci2*bz[ijk    ] + ci3*bz[ijk+ii1] )
                                        + cg2*( ci0*bz[ijk-ii1] + ci1*bz[ijk    ] + ci2*bz[ijk+ii1] + ci3*bz[ijk+ii2] )
                                        + cg3*( ci0*bz[ijk    ] + ci1*bz[ijk+ii1] + ci2*bz[ijk+ii2] + ci3*bz[ijk+ii3] ) ) )
                
                                    * cgi*dxi )
                
                                  + ( ( ( ( cg0*( ci0*w[ijk-jj3] + ci1*w[ijk-jj2] + ci2*w[ijk-jj1] + ci3*w[ijk    ] )
                                          + cg1*( ci0*w[ijk-jj2] + ci1*w[ijk-jj1] + ci2*w[ijk    ] + ci3*w[ijk+jj1] )
                                          + cg2*( ci0*w[ijk-jj1] + ci1*w[ijk    ] + ci2*w[ijk+jj1] + ci3*w[ijk+jj2] )
                                          + cg3*( ci0*w[ijk    ] + ci1*w[ijk+jj1] + ci2*w[ijk+jj2] + ci3*w[ijk+jj3] ) )
                
                                        * cgi*dyi )
                
                                      * ( cg0*( ci0*bz[ijk-jj3] + ci1*bz[ijk-jj2] + ci2*bz[ijk-jj1] + ci3*bz[ijk    ] )
                                        + cg1*( ci0*bz[ijk-jj2] + ci1*bz[ijk-jj1] + ci2*bz[ijk    ] + ci3*bz[ijk+jj1] )
                                        + cg2*( ci0*bz[ijk-jj1] + ci1*bz[ijk    ] + ci2*bz[ijk+jj1] + ci3*bz[ijk+jj2] )
                                        + cg3*( ci0*bz[ijk    ] + ci1*bz[ijk+jj1] + ci2*bz[ijk+jj2] + ci3*bz[ijk+jj3] ) ) )
                
                                    * cgi*dyi ) )
                
                                + ( ( ( ( cg0*( ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ] )
                                        + cg1*( ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1] )
                                        + cg2*( ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2] )
                                        + cg3*( ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3] ) )
                
                                      * dzhi4[k] )
                
                                    * ( cg0*( b[ijk-kk2] - bmean[k-2] ) + cg1*( b[ijk-kk1] - bmean[k-1] ) + cg2*( b[ijk    ] - bmean[k  ] ) + cg3*( b[ijk+kk1] - bmean[k+1] ) ) )
                
                                  * dzhi4[k] ) ) );
            }
    }

    k = grid.kend-1;
    bw_diss[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            bw_diss[k] -= ( ( 2.0 * visc )
            
                          * ( ( ( ( ( ( cg0*( ci0*w[ijk-ii3] + ci1*w[ijk-ii2] + ci2*w[ijk-ii1] + ci3*w[ijk    ] )
                                      + cg1*( ci0*w[ijk-ii2] + ci1*w[ijk-ii1] + ci2*w[ijk    ] + ci3*w[ijk+ii1] )
                                      + cg2*( ci0*w[ijk-ii1] + ci1*w[ijk    ] + ci2*w[ijk+ii1] + ci3*w[ijk+ii2] )
                                      + cg3*( ci0*w[ijk    ] + ci1*w[ijk+ii1] + ci2*w[ijk+ii2] + ci3*w[ijk+ii3] ) )
            
                                    * cgi*dxi )
            
                                  * ( cg0*( ci0*bz[ijk-ii3] + ci1*bz[ijk-ii2] + ci2*bz[ijk-ii1] + ci3*bz[ijk    ] )
                                    + cg1*( ci0*bz[ijk-ii2] + ci1*bz[ijk-ii1] + ci2*bz[ijk    ] + ci3*bz[ijk+ii1] )
                                    + cg2*( ci0*bz[ijk-ii1] + ci1*bz[ijk    ] + ci2*bz[ijk+ii1] + ci3*bz[ijk+ii2] )
                                    + cg3*( ci0*bz[ijk    ] + ci1*bz[ijk+ii1] + ci2*bz[ijk+ii2] + ci3*bz[ijk+ii3] ) ) )
            
                                * cgi*dxi )
            
                              + ( ( ( ( cg0*( ci0*w[ijk-jj3] + ci1*w[ijk-jj2] + ci2*w[ijk-jj1] + ci3*w[ijk    ] )
                                      + cg1*( ci0*w[ijk-jj2] + ci1*w[ijk-jj1] + ci2*w[ijk    ] + ci3*w[ijk+jj1] )
                                      + cg2*( ci0*w[ijk-jj1] + ci1*w[ijk    ] + ci2*w[ijk+jj1] + ci3*w[ijk+jj2] )
                                      + cg3*( ci0*w[ijk    ] + ci1*w[ijk+jj1] + ci2*w[ijk+jj2] + ci3*w[ijk+jj3] ) )
            
                                    * cgi*dyi )
            
                                  * ( cg0*( ci0*bz[ijk-jj3] + ci1*bz[ijk-jj2] + ci2*bz[ijk-jj1] + ci3*bz[ijk    ] )
                                    + cg1*( ci0*bz[ijk-jj2] + ci1*bz[ijk-jj1] + ci2*bz[ijk    ] + ci3*bz[ijk+jj1] )
                                    + cg2*( ci0*bz[ijk-jj1] + ci1*bz[ijk    ] + ci2*bz[ijk+jj1] + ci3*bz[ijk+jj2] )
                                    + cg3*( ci0*bz[ijk    ] + ci1*bz[ijk+jj1] + ci2*bz[ijk+jj2] + ci3*bz[ijk+jj3] ) ) )
            
                                * cgi*dyi ) )
            
                            + ( ( ( ( cg0*( ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ] )
                                    + cg1*( ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1] )
                                    + cg2*( ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2] )
                                    + cg3*( ti0*w[ijk-kk1] + ti1*w[ijk    ] + ti2*w[ijk+kk1] + ti3*w[ijk+kk2] ) )
            
                                  * dzhi4[k] )
            
                                * ( cg0*( b[ijk-kk2] - bmean[k-2] ) + cg1*( b[ijk-kk1] - bmean[k-1] ) + cg2*( b[ijk    ] - bmean[k  ] ) + cg3*( b[ijk+kk1] - bmean[k+1] ) ) )
            
                              * dzhi4[k] ) ) );


        }

    k = grid.kend;
    bw_diss[k] = 0;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            bw_diss[k] -= ( ( 2.0 * visc )
            
                          * ( ( ( ( ( ( cg0*( ci0*w[ijk-ii3] + ci1*w[ijk-ii2] + ci2*w[ijk-ii1] + ci3*w[ijk    ] )
                                      + cg1*( ci0*w[ijk-ii2] + ci1*w[ijk-ii1] + ci2*w[ijk    ] + ci3*w[ijk+ii1] )
                                      + cg2*( ci0*w[ijk-ii1] + ci1*w[ijk    ] + ci2*w[ijk+ii1] + ci3*w[ijk+ii2] )
                                      + cg3*( ci0*w[ijk    ] + ci1*w[ijk+ii1] + ci2*w[ijk+ii2] + ci3*w[ijk+ii3] ) )
            
                                    * cgi*dxi )
            
                                  * ( cg0*( ci0*bz[ijk-ii3] + ci1*bz[ijk-ii2] + ci2*bz[ijk-ii1] + ci3*bz[ijk    ] )
                                    + cg1*( ci0*bz[ijk-ii2] + ci1*bz[ijk-ii1] + ci2*bz[ijk    ] + ci3*bz[ijk+ii1] )
                                    + cg2*( ci0*bz[ijk-ii1] + ci1*bz[ijk    ] + ci2*bz[ijk+ii1] + ci3*bz[ijk+ii2] )
                                    + cg3*( ci0*bz[ijk    ] + ci1*bz[ijk+ii1] + ci2*bz[ijk+ii2] + ci3*bz[ijk+ii3] ) ) )
            
                                * cgi*dxi )
            
                              + ( ( ( ( cg0*( ci0*w[ijk-jj3] + ci1*w[ijk-jj2] + ci2*w[ijk-jj1] + ci3*w[ijk    ] )
                                      + cg1*( ci0*w[ijk-jj2] + ci1*w[ijk-jj1] + ci2*w[ijk    ] + ci3*w[ijk+jj1] )
                                      + cg2*( ci0*w[ijk-jj1] + ci1*w[ijk    ] + ci2*w[ijk+jj1] + ci3*w[ijk+jj2] )
                                      + cg3*( ci0*w[ijk    ] + ci1*w[ijk+jj1] + ci2*w[ijk+jj2] + ci3*w[ijk+jj3] ) )
            
                                    * cgi*dyi )
            
                                  * ( cg0*( ci0*bz[ijk-jj3] + ci1*bz[ijk-jj2] + ci2*bz[ijk-jj1] + ci3*bz[ijk    ] )
                                    + cg1*( ci0*bz[ijk-jj2] + ci1*bz[ijk-jj1] + ci2*bz[ijk    ] + ci3*bz[ijk+jj1] )
                                    + cg2*( ci0*bz[ijk-jj1] + ci1*bz[ijk    ] + ci2*bz[ijk+jj1] + ci3*bz[ijk+jj2] )
                                    + cg3*( ci0*bz[ijk    ] + ci1*bz[ijk+jj1] + ci2*bz[ijk+jj2] + ci3*bz[ijk+jj3] ) ) )
            
                                * cgi*dyi ) )
            
                            + ( ( ( ( tg0*( ci0*w[ijk-kk4] + ci1*w[ijk-kk3] + ci2*w[ijk-kk2] + ci3*w[ijk-kk1] )
                                    + tg1*( ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ] )
                                    + tg2*( ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1] )
                                    + tg3*( ti0*w[ijk-kk2] + ti1*w[ijk-kk1] + ti2*w[ijk    ] + ti3*w[ijk+kk1] ) )
            
                                  * dzhi4top )
            
                                * ( cg0*( b[ijk-kk2] - bmean[k-2] ) + cg1*( b[ijk-kk1] - bmean[k-1] ) + cg2*( b[ijk    ] - bmean[k  ] ) + cg3*( b[ijk+kk1] - bmean[k+1] ) ) )
            
                              * dzhi4top ) ) );
        }

    // 7. CALCULATE THE PRESSURE TRANSPORT TERM
    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        bw_pres[k] = 0;
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_pres[k] -= ( ( cg0*( ( p[ijk-kk2] - pmean[k-2] ) * ( b[ijk-kk2] - bmean[k-2] ) ) + cg1*( ( p[ijk-kk1] - pmean[k-1] ) * ( b[ijk-kk1] - bmean[k-1] ) ) + cg2*( ( p[ijk    ] - pmean[k  ] ) * ( b[ijk    ] - bmean[k  ] ) ) + cg3*( ( p[ijk+kk1] - pmean[k+1] ) * ( b[ijk+kk1] - bmean[k+1] ) ) ) * dzhi4[k] );
            }
    }


    master.sum(bw_shear, grid.kcells);
    master.sum(bw_turb , grid.kcells);
    master.sum(bw_visc , grid.kcells);
    master.sum(bw_rdstr, grid.kcells);
    master.sum(bw_buoy , grid.kcells);
    master.sum(bw_diss , grid.kcells);
    master.sum(bw_pres , grid.kcells);

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        bw_shear[k] /= n;
        bw_turb [k] /= n;
        bw_visc [k] /= n;
        bw_rdstr[k] /= n;
        bw_buoy [k] /= n;
        bw_diss [k] /= n;
        bw_pres [k] /= n;
    }
}

void Budget_4::calc_pe(double* restrict b, double* restrict zsort, double* restrict zsortbot, double* restrict zsorttop,
                       double* restrict z,
                       double* restrict bsort,
                       double* restrict pe_total, double* restrict pe_avail, double* restrict pe_bg,
                       double* restrict zsortprof)
{
    const int jj = grid.icells;
    const int kk1 = 1*grid.ijcells;
    const int kk2 = 2*grid.ijcells;
    const int kstart = grid.kstart;
    const int kend = grid.kend;

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        pe_total[k] = 0;
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk1;
                pe_total[k] -= b[ijk] * z[k];
            }
    }

    master.sum(pe_total, grid.kcells);

    int n = grid.itot*grid.jtot;
    for (int k=grid.kstart; k<grid.kend; ++k)
        pe_total[k] /= n;

    // now find out the available potential energy
    // int ks;
    double zsortval;
    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        zsortprof[k] = 0.;
        pe_bg    [k] = 0.;
        pe_avail [k] = 0.;
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk1;
                /*
                   ks  = k;
                   if (b[ijk] > bsort[k])
                   {
                   while (b[ijk] > bsort[ks] && ks < grid.kend-1)
                   ++ks;

                // linearly interpolate the height
                zsortval = z[ks-1] + (b[ijk]-bsort[ks-1])/(bsort[ks]-bsort[ks-1]) * (z[ks]-z[ks-1]);

                }
                else if (b[ijk] < bsort[k])
                {
                while (b[ijk] < bsort[ks] && ks > grid.kstart)
                --ks;

                // linearly interpolate the height
                zsortval = z[ks] + (b[ijk]-bsort[ks])/(bsort[ks+1]-bsort[ks]) * (z[ks+1]-z[ks]);
                }
                else
                zsortval = z[ks];
                */
                zsortval = calc_zsort(b[ijk], bsort, z, k);
                zsort[ijk] = zsortval;

                zsortprof[k] += zsortval;
                pe_bg    [k] -= zsortval*b[ijk];
                pe_avail [k] -= (z[k]-zsortval)*b[ijk];
            }
    }

    master.sum(zsortprof, grid.kcells);
    master.sum(pe_bg    , grid.kcells);
    master.sum(pe_avail , grid.kcells);

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        zsortprof[k] /= n;
        pe_bg    [k] /= n;
        pe_avail [k] /= n;
    }

    // now, calculate the boundary conditions for zsort
    // bottom bc
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk1;
            zsortbot[ij] = calc_zsort(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk] + ci3*b[ijk+kk1], bsort, z, kstart);
        }

    // top bc
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + (kend-1)*kk1;
            zsorttop[ij] = calc_zsort(ci0*b[ijk-kk1] + ci1*b[ijk] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2], bsort, z, kend-1);
        }

    // calculate the ghost cells at the bottom
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk1;
            zsort[ijk-kk1] = (8./3.)*zsortbot[ij] - 2.*zsort[ijk] + (1./3.)*zsort[ijk+kk1];
            zsort[ijk-kk2] = 8.*zsortbot[ij] - 9.*zsort[ijk] + 2.*zsort[ijk+kk1];
        }

    // calculate the ghost cells at the top
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + (kend-1)*kk1;
            zsort[ijk+kk1] = (8./3.)*zsorttop[ij] - 2.*zsort[ijk] + (1./3.)*zsort[ijk-kk1];
            zsort[ijk+kk2] = 8.*zsorttop[ij] - 9.*zsort[ijk] + 2.*zsort[ijk-kk1];
        }
}

double Budget_4::calc_zsort(double b, double* restrict bsort, double* restrict z, int k)
{
    double zsortval;
    int ks = k;

    if (b > bsort[k])
    {
        while (b > bsort[ks] && ks < grid.kend-1)
            ++ks;

        // linearly interpolate the height
        zsortval = z[ks-1] + (b-bsort[ks-1])/(bsort[ks]-bsort[ks-1]) * (z[ks]-z[ks-1]);

    }
    else if (b < bsort[k])
    {
        while (b < bsort[ks] && ks > grid.kstart)
            --ks;

        // linearly interpolate the height
        zsortval = z[ks] + (b-bsort[ks])/(bsort[ks+1]-bsort[ks]) * (z[ks+1]-z[ks]);
    }
    else
        zsortval = z[ks];

    return zsortval;
}

double Budget_4::calc_dzstardb(double b, double* restrict bsort, double* restrict z)
{
    // start the iteration below the grid to make sure not to miss values below the first full level
    int k = grid.kstart-1;
    while (bsort[k+1] < b && k < grid.kend)
        ++k;

    // our required value is in between bsort[k] and bsort[k+1]
    // calculate a spline of the form zstar(b) = a*zstar(b_k) + b*zstar(b_(k+1)) + c * zstar''(b_k) + d * zstar''(b_(k+1))
    double ca, cb;//, cc, cd;
    ca = (bsort[k+1]-b) / (bsort[k+1]-bsort[k]);
    cb = 1.-ca;
    // cc = (1./6.) * (ca*ca*ca-ca) * std::pow(bsort[k+1]-bsort[k],2);
    // cd = (1./6.) * (cb*cb*cb-cb) * std::pow(bsort[k+1]-bsort[k],2);

    // calculate the second derivatives using second order accuracy since the profile is very smooth
    double d2zstarb2k, d2zstarb2kp;
    d2zstarb2k  = 2.*((z[k+1]-z[k  ])/(bsort[k+1]-bsort[k  ]) - (z[k  ]-z[k-1])/(bsort[k  ]-bsort[k-1])) / (bsort[k+1]-bsort[k-1]);
    d2zstarb2kp = 2.*((z[k+2]-z[k+1])/(bsort[k+2]-bsort[k+1]) - (z[k+1]-z[k  ])/(bsort[k+1]-bsort[k  ])) / (bsort[k+2]-bsort[k  ]);

    // std::printf("CvH %E, %E, %E, %E, %E, %E\n", bsort[k], b, bsort[k+1], z[k], ca*z[k] + cb*z[k+1], z[k+1]);

    // the derivative is computed according to:
    // dzstar/db = (zstar[k+1]-zstar[k])/(b[k+1]-b[k])
    //           - (3*a^2-1)/6 * (b[k+1]-b[k])*zstar''(b_k)
    //           + (3*b^2-1)/6 * (b[k+1]-b[k])*zstar''(b_(k+1))

    double dzstardb = (z[k+1]-z[k]) / (bsort[k+1]-bsort[k])
        - (3.*ca*ca-1.)/6. * (bsort[k+1]-bsort[k])*d2zstarb2k
        + (3.*cb*cb-1.)/6. * (bsort[k+1]-bsort[k])*d2zstarb2kp;

    return dzstardb;
}


void Budget_4::calc_pe_budget(double* restrict w, double* restrict b, double* restrict bz, double* restrict bztop,
                              double* restrict pe_turb, double* restrict pe_visc, double* restrict pe_bous,
                              double* restrict z, double* restrict zh, double* restrict dzi4, double* restrict dzhi4,
                              double visc)
{
    const int jj1 = 1*grid.icells;
    const int kk1 = 1*grid.ijcells;
    const int kk2 = 2*grid.ijcells;
    const int kk3 = 3*grid.ijcells;
    const int kstart = grid.kstart;
    const int kend = grid.kend;

    const double zsize = grid.zsize;

    // first, calculate the Boussinesq term (kappa*db/dz). Here bz contains the buoyancy and not the PE yet
    // bottom boundary
    pe_bous[kstart] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + kstart*kk1;
            pe_bous[kstart] += visc * ( cg0*(bi0*b[ijk-kk2] + bi1*b[ijk-kk1] + bi2*b[ijk    ] + bi3*b[ijk+kk1])
                                      + cg1*(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1])
                                      + cg2*(ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2])
                                      + cg3*(ci0*b[ijk    ] + ci1*b[ijk+kk1] + ci2*b[ijk+kk2] + ci3*b[ijk+kk3]) )
                                    * dzi4[kstart];
        }

    for (int k=grid.kstart+1; k<grid.kend-1; ++k)
    {
        pe_bous[k] = 0.;
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                pe_bous[k] += visc * ( cg0*(ci0*b[ijk-kk3] + ci1*b[ijk-kk2] + ci2*b[ijk-kk1] + ci3*b[ijk    ])
                                     + cg1*(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1])
                                     + cg2*(ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2])
                                     + cg3*(ci0*b[ijk    ] + ci1*b[ijk+kk1] + ci2*b[ijk+kk2] + ci3*b[ijk+kk3]) )
                                   * dzi4[k];
            }
    }

    // top boundary
    pe_bous[kend-1] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + (kend-1)*kk1;
            pe_bous[kend-1] += visc * ( cg0*(ci0*b[ijk-kk3] + ci1*b[ijk-kk2] + ci2*b[ijk-kk1] + ci3*b[ijk    ])
                                      + cg1*(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1])
                                      + cg2*(ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2])
                                      + cg3*(ti0*b[ijk-kk1] + ti1*b[ijk    ] + ti2*b[ijk+kk1] + ti3*b[ijk+kk2]) )
                                    * dzi4[kend-1];
        }

    // now, convert the buoyancy field into a potential energy field
    // first, before destroying the field, calculate the potential energy at the top
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ij  = i + j*jj1;
            const int ijk = i + j*jj1 + (kend-1)*kk1;
            bztop[ij] = -zsize*(ci0*b[ijk-kk1] + ci1*b[ijk] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2]);
        }

    // calculate the potential energy
    for (int k=grid.kstart; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bz[ijk] = -b[ijk] * z[k];
            }

    // calculate the ghost cells at the bottom, making use of the fact that bz = 0
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + kstart*kk1;
            bz[ijk-kk1] = - 2.*bz[ijk] + (1./3.)*bz[ijk+kk1];
            bz[ijk-kk2] = - 9.*bz[ijk] + 2.*bz[ijk+kk1];
        }

    // calculate the ghost cells at the top
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ij  = i + j*jj1;
            const int ijk = i + j*jj1 + (kend-1)*kk1;
            bz[ijk+kk1] = (8./3.)*bztop[ij] - 2.*bz[ijk] + (1./3.)*bz[ijk-kk1];
            bz[ijk+kk2] = 8.*bztop[ij] - 9.*bz[ijk] + 2.*bz[ijk-kk1];
        }

    // calculate the advective transport term
    // bottom boundary
    pe_turb[kstart] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + kstart*kk1;
            pe_turb[kstart] -= ( cg0*(w[ijk-kk1] * (bi0*bz[ijk-kk2] + bi1*bz[ijk-kk1] + bi2*bz[ijk    ] + bi3*bz[ijk+kk1]))
                               + cg1*(w[ijk    ] * (ci0*bz[ijk-kk2] + ci1*bz[ijk-kk1] + ci2*bz[ijk    ] + ci3*bz[ijk+kk1]))
                               + cg2*(w[ijk+kk1] * (ci0*bz[ijk-kk1] + ci1*bz[ijk    ] + ci2*bz[ijk+kk1] + ci3*bz[ijk+kk2]))
                               + cg3*(w[ijk+kk2] * (ci0*bz[ijk    ] + ci1*bz[ijk+kk1] + ci2*bz[ijk+kk2] + ci3*bz[ijk+kk3])) )
                             * dzi4[kstart];
        }

    for (int k=grid.kstart+1; k<grid.kend-1; ++k)
    {
        pe_turb[k] = 0.;
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                pe_turb[k] -= ( cg0*(w[ijk-kk1] * (ci0*bz[ijk-kk3] + ci1*bz[ijk-kk2] + ci2*bz[ijk-kk1] + ci3*bz[ijk    ]))
                              + cg1*(w[ijk    ] * (ci0*bz[ijk-kk2] + ci1*bz[ijk-kk1] + ci2*bz[ijk    ] + ci3*bz[ijk+kk1]))
                              + cg2*(w[ijk+kk1] * (ci0*bz[ijk-kk1] + ci1*bz[ijk    ] + ci2*bz[ijk+kk1] + ci3*bz[ijk+kk2]))
                              + cg3*(w[ijk+kk2] * (ci0*bz[ijk    ] + ci1*bz[ijk+kk1] + ci2*bz[ijk+kk2] + ci3*bz[ijk+kk3])) )
                            * dzi4[k];
            }
    }

    // top boundary
    pe_turb[kend-1] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + (kend-1)*kk1;
            pe_turb[kend-1] -= ( cg0*(w[ijk-kk1] * (ci0*bz[ijk-kk3] + ci1*bz[ijk-kk2] + ci2*bz[ijk-kk1] + ci3*bz[ijk    ]))
                               + cg1*(w[ijk    ] * (ci0*bz[ijk-kk2] + ci1*bz[ijk-kk1] + ci2*bz[ijk    ] + ci3*bz[ijk+kk1]))
                               + cg2*(w[ijk+kk1] * (ci0*bz[ijk-kk1] + ci1*bz[ijk    ] + ci2*bz[ijk+kk1] + ci3*bz[ijk+kk2]))
                               + cg3*(w[ijk+kk2] * (ti0*bz[ijk-kk1] + ti1*bz[ijk    ] + ti2*bz[ijk+kk1] + ti3*bz[ijk+kk2])) )
                             * dzi4[kend-1];
        }

    // calculate the diffusion of potential energy
    // bottom boundary
    pe_visc[kstart] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + kstart*kk1;
            pe_visc[kstart] -= visc * ( cg0*zh[kstart-1]*(bg0*b[ijk-kk2] + bg1*b[ijk-kk1] + bg2*b[ijk    ] + bg3*b[ijk+kk1]) * dzhi4[kstart-1]
                                      + cg1*zh[kstart  ]*(cg0*b[ijk-kk2] + cg1*b[ijk-kk1] + cg2*b[ijk    ] + cg3*b[ijk+kk1]) * dzhi4[kstart  ]
                                      + cg2*zh[kstart+1]*(cg0*b[ijk-kk1] + cg1*b[ijk    ] + cg2*b[ijk+kk1] + cg3*b[ijk+kk2]) * dzhi4[kstart+1]
                                      + cg3*zh[kstart+2]*(cg0*b[ijk    ] + cg1*b[ijk+kk1] + cg2*b[ijk+kk2] + cg3*b[ijk+kk3]) * dzhi4[kstart+2] )
                             * dzi4[kstart];
        }

    for (int k=grid.kstart+1; k<grid.kend-1; ++k)
    {
        pe_visc[k] = 0.;
        for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                pe_visc[k] -= visc * ( cg0*zh[k-1]*(cg0*b[ijk-kk3] + cg1*b[ijk-kk2] + cg2*b[ijk-kk1] + cg3*b[ijk    ]) * dzhi4[k-1]
                                     + cg1*zh[k  ]*(cg0*b[ijk-kk2] + cg1*b[ijk-kk1] + cg2*b[ijk    ] + cg3*b[ijk+kk1]) * dzhi4[k  ]
                                     + cg2*zh[k+1]*(cg0*b[ijk-kk1] + cg1*b[ijk    ] + cg2*b[ijk+kk1] + cg3*b[ijk+kk2]) * dzhi4[k+1]
                                     + cg3*zh[k+2]*(cg0*b[ijk    ] + cg1*b[ijk+kk1] + cg2*b[ijk+kk2] + cg3*b[ijk+kk3]) * dzhi4[k+2] )
                            * dzi4[k];
            }
    }

    // top boundary
    pe_visc[kend-1] = 0.;
    for (int j=grid.jstart; j<grid.jend; ++j)
#pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj1 + (kend-1)*kk1;
            pe_visc[kend-1] -= visc * ( cg0*zh[kend-2]*(cg0*b[ijk-kk3] + cg1*b[ijk-kk2] + cg2*b[ijk-kk1] + cg3*b[ijk    ]) * dzhi4[kend-2]
                                      + cg1*zh[kend-1]*(cg0*b[ijk-kk2] + cg1*b[ijk-kk1] + cg2*b[ijk    ] + cg3*b[ijk+kk1]) * dzhi4[kend-1]
                                      + cg2*zh[kend  ]*(cg0*b[ijk-kk1] + cg1*b[ijk    ] + cg2*b[ijk+kk1] + cg3*b[ijk+kk2]) * dzhi4[kend  ]
                                      + cg3*zh[kend+1]*(tg0*b[ijk-kk1] + tg1*b[ijk    ] + tg2*b[ijk+kk1] + tg3*b[ijk+kk2]) * dzhi4[kend+1] )
                                    * dzi4[kend-1];
        }

    master.sum(pe_turb, grid.kcells);
    master.sum(pe_visc, grid.kcells);
    master.sum(pe_bous, grid.kcells);

    int n = grid.itot*grid.jtot;
    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        pe_turb[k] /= n;
        pe_visc[k] /= n;
        pe_bous[k] /= n;
    }
}

// void Budget_4::calcBpeBudget(double* restrict w, double* restrict b,
//                            double* restrict bz, double* restrict bzbot, double* restrict bztop,
//                            double* restrict bpe_turb, double* restrict bpe_visc, double* restrict bpe_diss,
//                            double* restrict bsort,
//                            double* restrict z, double* restrict dzi4, double* restrict dzhi4,
//                            double visc)
// {
//   int ijk,ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
//   int kstart,kend;
//   double dzstardb;
//
//   ii1 = 1;
//   ii2 = 2;
//   ii3 = 3;
//   jj1 = 1*grid.icells;
//   jj2 = 2*grid.icells;
//   jj3 = 3*grid.icells;
//   kk1 = 1*grid.ijcells;
//   kk2 = 2*grid.ijcells;
//   kk3 = 3*grid.ijcells;
//   kstart = grid.kstart;
//   kend = grid.kend;
//
//   // calculate the diffusion of potential energy
//   // bottom boundary
//   bpe_visc[kstart] = 0.;
//   for (int j=grid.jstart; j<grid.jend; ++j)
// #pragma ivdep
//     for (int i=grid.istart; i<grid.iend; ++i)
//     {
//       const int ijk = i + j*jj1 + kstart*kk1;
//       bpe_visc[kstart] -= visc *
//                             ( cg0*(bg0*b [ijk-kk2] + bg1*b [ijk-kk1] + bg2*b [ijk    ] + bg3*b [ijk+kk1]) * dzhi4[kstart-1]
//                                  *(bi0*bz[ijk-kk2] + bi1*bz[ijk-kk1] + bi2*bz[ijk    ] + bi3*bz[ijk+kk1])
//                             + cg1*(cg0*b [ijk-kk2] + cg1*b [ijk-kk1] + cg2*b [ijk    ] + cg3*b [ijk+kk1]) * dzhi4[kstart  ]
//                                  *(ci0*bz[ijk-kk2] + ci1*bz[ijk-kk1] + ci2*bz[ijk    ] + ci3*bz[ijk+kk1])
//                             + cg2*(cg0*b [ijk-kk1] + cg1*b [ijk    ] + cg2*b [ijk+kk1] + cg3*b [ijk+kk2]) * dzhi4[kstart+1]
//                                  *(ci0*bz[ijk-kk1] + ci1*bz[ijk    ] + ci2*bz[ijk+kk1] + ci3*bz[ijk+kk2])
//                             + cg3*(cg0*b [ijk    ] + cg1*b [ijk+kk1] + cg2*b [ijk+kk2] + cg3*b [ijk+kk3]) * dzhi4[kstart+2]
//                                  *(ci0*bz[ijk    ] + ci1*bz[ijk+kk1] + ci2*bz[ijk+kk2] + ci3*bz[ijk+kk3]) )
//                             * dzi4[kstart];
//     }
//
//   for (int k=grid.kstart+1; k<grid.kend-1; ++k)
//   {
//     bpe_visc[k] = 0.;
//     for (int j=grid.jstart; j<grid.jend; ++j)
// #pragma ivdep
//       for (int i=grid.istart; i<grid.iend; ++i)
//       {
//         const int ijk = i + j*jj1 + k*kk1;
//         bpe_visc[k] -= visc *
//                          ( cg0*(cg0*b [ijk-kk3] + cg1*b [ijk-kk2] + cg2*b [ijk-kk1] + cg3*b [ijk    ]) * dzhi4[k-1]
//                               *(ci0*bz[ijk-kk3] + ci1*bz[ijk-kk2] + ci2*bz[ijk-kk1] + ci3*bz[ijk    ])
//                          + cg1*(cg0*b [ijk-kk2] + cg1*b [ijk-kk1] + cg2*b [ijk    ] + cg3*b [ijk+kk1]) * dzhi4[k  ]
//                               *(ci0*bz[ijk-kk2] + ci1*bz[ijk-kk1] + ci2*bz[ijk    ] + ci3*bz[ijk+kk1])
//                          + cg2*(cg0*b [ijk-kk1] + cg1*b [ijk    ] + cg2*b [ijk+kk1] + cg3*b [ijk+kk2]) * dzhi4[k+1]
//                               *(ci0*bz[ijk-kk1] + ci1*bz[ijk    ] + ci2*bz[ijk+kk1] + ci3*bz[ijk+kk2])
//                          + cg3*(cg0*b [ijk    ] + cg1*b [ijk+kk1] + cg2*b [ijk+kk2] + cg3*b [ijk+kk3]) * dzhi4[k+2]
//                               *(ci0*bz[ijk    ] + ci1*bz[ijk+kk1] + ci2*bz[ijk+kk2] + ci3*bz[ijk+kk3]) )
//                          * dzi4[k];
//       }
//   }
//
//   // top boundary
//   bpe_visc[kend-1] = 0.;
//   for (int j=grid.jstart; j<grid.jend; ++j)
// #pragma ivdep
//     for (int i=grid.istart; i<grid.iend; ++i)
//     {
//       const int ijk = i + j*jj1 + (kend-1)*kk1;
//       bpe_visc[kend-1] -= visc *
//                             ( cg0*(cg0*b [ijk-kk3] + cg1*b [ijk-kk2] + cg2*b [ijk-kk1] + cg3*b [ijk    ]) * dzhi4[kend-2]
//                                  *(ci0*bz[ijk-kk3] + ci1*bz[ijk-kk2] + ci2*bz[ijk-kk1] + ci3*bz[ijk    ])
//                             + cg1*(cg0*b [ijk-kk2] + cg1*b [ijk-kk1] + cg2*b [ijk    ] + cg3*b [ijk+kk1]) * dzhi4[kend-1]
//                                  *(ci0*bz[ijk-kk2] + ci1*bz[ijk-kk1] + ci2*bz[ijk    ] + ci3*bz[ijk+kk1])
//                             + cg2*(cg0*b [ijk-kk1] + cg1*b [ijk    ] + cg2*b [ijk+kk1] + cg3*b [ijk+kk2]) * dzhi4[kend  ]
//                                  *(ci0*bz[ijk-kk1] + ci1*bz[ijk    ] + ci2*bz[ijk+kk1] + ci3*bz[ijk+kk2])
//                             + cg3*(tg0*b [ijk-kk1] + tg1*b [ijk    ] + tg2*b [ijk+kk1] + tg3*b [ijk+kk2]) * dzhi4[kend+1]
//                                  *(ti0*bz[ijk-kk1] + ti1*bz[ijk    ] + ti2*bz[ijk+kk1] + ti3*bz[ijk+kk2]) )
//                             * dzi4[kend-1];
//     }
//
//   // calculate the dissipation term
//   double dxi,dyi;
//   dxi = 1./grid.dx;
//   dyi = 1./grid.dy;
//
//   bpe_diss[kstart] = 0.;
//   for (int j=grid.jstart; j<grid.jend; j++)
// #pragma ivdep
//     for (int i=grid.istart; i<grid.iend; i++)
//     {
//       ijk  = i + j*jj1 + kstart*kk1;
//       dzstardb = calc_dzstardb(b[ijk], bsort, z);
//       bpe_diss[kstart] += visc * dzstardb * (
//                          std::pow( ( cg0*(ci0*b[ijk-ii3] + ci1*b[ijk-ii2] + ci2*b[ijk-ii1] + ci3*b[ijk    ])
//                                    + cg1*(ci0*b[ijk-ii2] + ci1*b[ijk-ii1] + ci2*b[ijk    ] + ci3*b[ijk+ii1])
//                                    + cg2*(ci0*b[ijk-ii1] + ci1*b[ijk    ] + ci2*b[ijk+ii1] + ci3*b[ijk+ii2])
//                                    + cg3*(ci0*b[ijk    ] + ci1*b[ijk+ii1] + ci2*b[ijk+ii2] + ci3*b[ijk+ii3]) ) * cgi*dxi, 2)
//
//                        + std::pow( ( cg0*(ci0*b[ijk-jj3] + ci1*b[ijk-jj2] + ci2*b[ijk-jj1] + ci3*b[ijk    ])
//                                    + cg1*(ci0*b[ijk-jj2] + ci1*b[ijk-jj1] + ci2*b[ijk    ] + ci3*b[ijk+jj1])
//                                    + cg2*(ci0*b[ijk-jj1] + ci1*b[ijk    ] + ci2*b[ijk+jj1] + ci3*b[ijk+jj2])
//                                    + cg3*(ci0*b[ijk    ] + ci1*b[ijk+jj1] + ci2*b[ijk+jj2] + ci3*b[ijk+jj3]) ) * cgi*dyi, 2)
//
//                        + std::pow( ( cg0*(bi0*b[ijk-kk2] + bi1*b[ijk-kk1] + bi2*b[ijk    ] + bi3*b[ijk+kk1])
//                                    + cg1*(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1])
//                                    + cg2*(ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2])
//                                    + cg3*(ci0*b[ijk    ] + ci1*b[ijk+kk1] + ci2*b[ijk+kk2] + ci3*b[ijk+kk3]) ) * dzi4[kstart], 2) );
//     }
//
//   // interior
//   for (int k=grid.kstart+1; k<grid.kend-1; k++)
//   {
//     bpe_diss[k] = 0.;
//     for (int j=grid.jstart; j<grid.jend; j++)
// #pragma ivdep
//       for (int i=grid.istart; i<grid.iend; i++)
//       {
//         ijk  = i + j*jj1 + k*kk1;
//         dzstardb = calc_dzstardb(b[ijk], bsort, z);
//         bpe_diss[k] += visc * dzstardb * (
//                         std::pow( ( cg0*(ci0*b[ijk-ii3] + ci1*b[ijk-ii2] + ci2*b[ijk-ii1] + ci3*b[ijk    ])
//                                   + cg1*(ci0*b[ijk-ii2] + ci1*b[ijk-ii1] + ci2*b[ijk    ] + ci3*b[ijk+ii1])
//                                   + cg2*(ci0*b[ijk-ii1] + ci1*b[ijk    ] + ci2*b[ijk+ii1] + ci3*b[ijk+ii2])
//                                   + cg3*(ci0*b[ijk    ] + ci1*b[ijk+ii1] + ci2*b[ijk+ii2] + ci3*b[ijk+ii3]) ) * cgi*dxi, 2)
//
//                       + std::pow( ( cg0*(ci0*b[ijk-jj3] + ci1*b[ijk-jj2] + ci2*b[ijk-jj1] + ci3*b[ijk    ])
//                                   + cg1*(ci0*b[ijk-jj2] + ci1*b[ijk-jj1] + ci2*b[ijk    ] + ci3*b[ijk+jj1])
//                                   + cg2*(ci0*b[ijk-jj1] + ci1*b[ijk    ] + ci2*b[ijk+jj1] + ci3*b[ijk+jj2])
//                                   + cg3*(ci0*b[ijk    ] + ci1*b[ijk+jj1] + ci2*b[ijk+jj2] + ci3*b[ijk+jj3]) ) * cgi*dyi, 2)
//
//                       + std::pow( ( cg0*(ci0*b[ijk-kk3] + ci1*b[ijk-kk2] + ci2*b[ijk-kk1] + ci3*b[ijk    ])
//                                   + cg1*(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1])
//                                   + cg2*(ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2])
//                                   + cg3*(ci0*b[ijk    ] + ci1*b[ijk+kk1] + ci2*b[ijk+kk2] + ci3*b[ijk+kk3]) ) * dzi4[k], 2) );
//       }
//   }
//
//   // top
//   bpe_diss[kend-1] = 0.;
//   for (int j=grid.jstart; j<grid.jend; j++)
// #pragma ivdep
//     for (int i=grid.istart; i<grid.iend; i++)
//     {
//       const int ijk = i + j*jj1 + (kend-1)*kk1;
//       dzstardb = calc_dzstardb(b[ijk], bsort, z);
//       bpe_diss[kend-1] += visc * dzstardb * (
//                       std::pow( ( cg0*(ci0*b[ijk-ii3] + ci1*b[ijk-ii2] + ci2*b[ijk-ii1] + ci3*b[ijk    ])
//                                 + cg1*(ci0*b[ijk-ii2] + ci1*b[ijk-ii1] + ci2*b[ijk    ] + ci3*b[ijk+ii1])
//                                 + cg2*(ci0*b[ijk-ii1] + ci1*b[ijk    ] + ci2*b[ijk+ii1] + ci3*b[ijk+ii2])
//                                 + cg3*(ci0*b[ijk    ] + ci1*b[ijk+ii1] + ci2*b[ijk+ii2] + ci3*b[ijk+ii3]) ) * cgi*dxi, 2)
//
//                     + std::pow( ( cg0*(ci0*b[ijk-jj3] + ci1*b[ijk-jj2] + ci2*b[ijk-jj1] + ci3*b[ijk    ])
//                                 + cg1*(ci0*b[ijk-jj2] + ci1*b[ijk-jj1] + ci2*b[ijk    ] + ci3*b[ijk+jj1])
//                                 + cg2*(ci0*b[ijk-jj1] + ci1*b[ijk    ] + ci2*b[ijk+jj1] + ci3*b[ijk+jj2])
//                                 + cg3*(ci0*b[ijk    ] + ci1*b[ijk+jj1] + ci2*b[ijk+jj2] + ci3*b[ijk+jj3]) ) * cgi*dyi, 2)
//
//                     + std::pow( ( cg0*(ci0*b[ijk-kk3] + ci1*b[ijk-kk2] + ci2*b[ijk-kk1] + ci3*b[ijk    ])
//                                 + cg1*(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1])
//                                 + cg2*(ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2])
//                                 + cg3*(ti0*b[ijk-kk1] + ti1*b[ijk    ] + ti2*b[ijk+kk1] + ti3*b[ijk+kk2]) ) * dzi4[kend-1], 2) );
//     }
//
//   /*
//   // CONVERT THE BZ FIELD INTO BACKGROUND POTENTIAL ENERGY
//   // first, calculate the potential energy at the bottom, the bot field contains the zsort at the bottom boundary
//   for (int j=grid.jstart; j<grid.jend; ++j)
// #pragma ivdep
//     for (int i=grid.istart; i<grid.iend; ++i)
//     {
//       const int ij  = i + j*jj1;
//       const int ijk = i + j*jj1 + kstart*kk1;
//       bzbot[ij] *= -(ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk] + ci3*b[ijk+kk1]);
//     }
//
//
//   // calculate the potential energy at the top, the top field contains the zsort at the top boundary
//   for (int j=grid.jstart; j<grid.jend; ++j)
// #pragma ivdep
//     for (int i=grid.istart; i<grid.iend; ++i)
//     {
//       const int ij  = i + j*jj1;
//       const int ijk = i + j*jj1 + (kend-1)*kk1;
//       bztop[ij] *= -(ci0*b[ijk-kk1] + ci1*b[ijk] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2]);
//     }
//
//   // calculate the potential energy
//   for (int k=grid.kstart; k<grid.kend; ++k)
//     for (int j=grid.jstart; j<grid.jend; ++j)
// #pragma ivdep
//       for (int i=grid.istart; i<grid.iend; ++i)
//       {
//         const int ijk = i + j*jj1 + k*kk1;
//         bz[ijk] = -b[ijk] * bz[ijk];
//       }
//
//   // calculate the ghost cells at the bottom
//   for (int j=grid.jstart; j<grid.jend; ++j)
// #pragma ivdep
//     for (int i=grid.istart; i<grid.iend; ++i)
//     {
//       const int ij  = i + j*jj1;
//       const int ijk = i + j*jj1 + kstart*kk1;
//       bz[ijk-kk1] = (8./3.)*bzbot[ij] - 2.*bz[ijk] + (1./3.)*bz[ijk+kk1];
//       bz[ijk-kk2] = 8.*bzbot[ij] - 9.*bz[ijk] + 2.*bz[ijk+kk1];
//     }
//
//   // calculate the ghost cells at the top
//   for (int j=grid.jstart; j<grid.jend; ++j)
// #pragma ivdep
//     for (int i=grid.istart; i<grid.iend; ++i)
//     {
//       const int ij  = i + j*jj1;
//       const int ijk = i + j*jj1 + (kend-1)*kk1;
//       bz[ijk+kk1] = (8./3.)*bztop[ij] - 2.*bz[ijk] + (1./3.)*bz[ijk-kk1];
//       bz[ijk+kk2] = 8.*bztop[ij] - 9.*bz[ijk] + 2.*bz[ijk-kk1];
//     }
//   */
//
//   // calculate the advective transport term
//   // bottom boundary
//   bpe_turb[kstart] = 0.;
//   for (int j=grid.jstart; j<grid.jend; ++j)
// #pragma ivdep
//     for (int i=grid.istart; i<grid.iend; ++i)
//     {
//       const int ijk = i + j*jj1 + kstart*kk1;
//       bpe_turb[kstart] += bz[ijk]*
//                           ( cg0*(w[ijk-kk1] * (bi0*b[ijk-kk2] + bi1*b[ijk-kk1] + bi2*b[ijk    ] + bi3*b[ijk+kk1]))
//                           + cg1*(w[ijk    ] * (ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1]))
//                           + cg2*(w[ijk+kk1] * (ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2]))
//                           + cg3*(w[ijk+kk2] * (ci0*b[ijk    ] + ci1*b[ijk+kk1] + ci2*b[ijk+kk2] + ci3*b[ijk+kk3])) )
//                           * dzi4[kstart];
//     }
//
//   for (int k=grid.kstart+1; k<grid.kend-1; ++k)
//   {
//     bpe_turb[k] = 0.;
//     for (int j=grid.jstart; j<grid.jend; ++j)
// #pragma ivdep
//       for (int i=grid.istart; i<grid.iend; ++i)
//       {
//         const int ijk = i + j*jj1 + k*kk1;
//         bpe_turb[k] += bz[ijk]*
//                        ( cg0*(w[ijk-kk1] * (ci0*b[ijk-kk3] + ci1*b[ijk-kk2] + ci2*b[ijk-kk1] + ci3*b[ijk    ]))
//                        + cg1*(w[ijk    ] * (ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1]))
//                        + cg2*(w[ijk+kk1] * (ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2]))
//                        + cg3*(w[ijk+kk2] * (ci0*b[ijk    ] + ci1*b[ijk+kk1] + ci2*b[ijk+kk2] + ci3*b[ijk+kk3])) )
//                        * dzi4[k];
//       }
//   }
//
//   // top boundary
//   bpe_turb[kend-1] = 0.;
//   for (int j=grid.jstart; j<grid.jend; ++j)
// #pragma ivdep
//     for (int i=grid.istart; i<grid.iend; ++i)
//     {
//       const int ijk = i + j*jj1 + (kend-1)*kk1;
//       bpe_turb[kend-1] += bz[ijk]*
//                           ( cg0*(w[ijk-kk1] * (ci0*b[ijk-kk3] + ci1*b[ijk-kk2] + ci2*b[ijk-kk1] + ci3*b[ijk    ]))
//                           + cg1*(w[ijk    ] * (ci0*b[ijk-kk2] + ci1*b[ijk-kk1] + ci2*b[ijk    ] + ci3*b[ijk+kk1]))
//                           + cg2*(w[ijk+kk1] * (ci0*b[ijk-kk1] + ci1*b[ijk    ] + ci2*b[ijk+kk1] + ci3*b[ijk+kk2]))
//                           + cg3*(w[ijk+kk2] * (ti0*b[ijk-kk1] + ti1*b[ijk    ] + ti2*b[ijk+kk1] + ti3*b[ijk+kk2])) )
//                           * dzi4[kend-1];
//     }
//
//   master.sum(bpe_turb, grid.kcells);
//   master.sum(bpe_visc, grid.kcells);
//   master.sum(bpe_diss, grid.kcells);
//
//   int n = grid.itot*grid.jtot;
//   for (int k=grid.kstart; k<grid.kend; ++k)
//   {
//     bpe_turb[k] /= n;
//     bpe_visc[k] /= n;
//     bpe_diss[k] /= n;
//   }
// }
