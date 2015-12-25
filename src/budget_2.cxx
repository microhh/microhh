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
#include "diff.h"
#include "advec.h"
#include "stats.h"
#include <netcdfcpp.h>

#include "budget.h"
#include "budget_2.h"

using namespace Finite_difference::O2;

Budget_2::Budget_2(Input* inputin, Master* masterin, Grid* gridin, Fields* fieldsin, Thermo* thermoin, Diff* diffin, Advec* advecin, Stats* statsin) :
    Budget(inputin, masterin, gridin, fieldsin, thermoin, diffin, advecin, statsin)
{
    umodel = 0;
    vmodel = 0;

    // The LES flux budget requires one additional ghost cell in the horizontal
    if(diff.get_switch() == "smag2")
    {
        const int igc = 2;
        const int jgc = 2;
        const int kgc = 1;

        grid.set_minimum_ghost_cells(igc, jgc, kgc);
    }
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
    stats.add_prof("uw_shear" , "Shear production term in UW budget" , "m2 s-3", "zh");
    stats.add_prof("vw_shear" , "Shear production term in VW budget" , "m2 s-3", "zh");

    stats.add_prof("u2_turb" , "Turbulent transport term in U2 budget" , "m2 s-3", "z" );
    stats.add_prof("v2_turb" , "Turbulent transport term in V2 budget" , "m2 s-3", "z" );
    stats.add_prof("w2_turb" , "Turbulent transport term in W2 budget" , "m2 s-3", "zh");
    stats.add_prof("tke_turb", "Turbulent transport term in TKE budget", "m2 s-3", "z" );
    stats.add_prof("uw_turb" , "Turbulent transport term in UW budget" , "m2 s-3", "zh");
    stats.add_prof("vw_turb" , "Turbulent transport term in VW budget" , "m2 s-3", "zh");

    stats.add_prof("u2_visc" , "Viscous transport term in U2 budget" , "m2 s-3", "z" );
    stats.add_prof("v2_visc" , "Viscous transport term in V2 budget" , "m2 s-3", "z" );
    stats.add_prof("w2_visc" , "Viscous transport term in W2 budget" , "m2 s-3", "zh");
    stats.add_prof("tke_visc", "Viscous transport term in TKE budget", "m2 s-3", "z" );
    stats.add_prof("uw_visc" , "Viscous transport term in UW budget" , "m2 s-3", "zh");
    stats.add_prof("vw_visc" , "Viscous transport term in VW budget" , "m2 s-3", "zh");

    stats.add_prof("u2_diss" , "Dissipation term in U2 budget" , "m2 s-3", "z" );
    stats.add_prof("v2_diss" , "Dissipation term in V2 budget" , "m2 s-3", "z" );
    stats.add_prof("w2_diss" , "Dissipation term in W2 budget" , "m2 s-3", "zh");
    stats.add_prof("tke_diss", "Dissipation term in TKE budget", "m2 s-3", "z" );
    stats.add_prof("uw_diss" , "Dissipation term in UW budget" , "m2 s-3", "zh");
    stats.add_prof("vw_diss" , "Dissipation term in VW budget" , "m2 s-3", "zh");

    stats.add_prof("w2_pres" , "Pressure transport term in W2 budget" , "m2 s-3", "zh");
    stats.add_prof("tke_pres", "Pressure transport term in TKE budget", "m2 s-3", "z" );
    stats.add_prof("uw_pres" , "Pressure transport term in UW budget" , "m2 s-3", "zh");
    stats.add_prof("vw_pres" , "Pressure transport term in VW budget" , "m2 s-3", "zh");

    stats.add_prof("u2_rdstr", "Pressure redistribution term in U2 budget", "m2 s-3", "z" );
    stats.add_prof("v2_rdstr", "Pressure redistribution term in V2 budget", "m2 s-3", "z" );
    stats.add_prof("w2_rdstr", "Pressure redistribution term in W2 budget", "m2 s-3", "zh");
    stats.add_prof("uw_rdstr", "Pressure redistribution term in UW budget", "m2 s-3", "zh");
    stats.add_prof("vw_rdstr", "Pressure redistribution term in VW budget", "m2 s-3", "zh");

    if (thermo.get_switch() != "0")
    {
        stats.add_prof("w2_buoy" , "Buoyancy production/destruction term in W2 budget" , "m2 s-3", "zh");
        stats.add_prof("tke_buoy", "Buoyancy production/destruction term in TKE budget", "m2 s-3", "z" );
        stats.add_prof("uw_buoy" , "Buoyancy production/destruction term in UW budget" , "m2 s-3", "zh");
        stats.add_prof("vw_buoy" , "Buoyancy production/destruction term in VW budget" , "m2 s-3", "zh");
    }
}

void Budget_2::exec_stats(Mask* m)
{
    // calculate the mean of the fields
    grid.calc_mean(umodel, fields.u->data, grid.kcells);
    grid.calc_mean(vmodel, fields.v->data, grid.kcells);

    // Calculate kinetic and turbulent kinetic energy
    calc_kinetic_energy(m->profs["ke"].data, m->profs["tke"].data,
                        fields.u->data, fields.v->data, fields.w->data, umodel, vmodel, grid.utrans, grid.vtrans);

    if(advec.get_switch() != "0")
    {
        // Calculate the shear production and turbulent transport terms
        calc_advection_terms(m->profs["u2_shear"].data, m->profs["v2_shear"].data, m->profs["tke_shear"].data,
                             m->profs["uw_shear"].data, m->profs["vw_shear"].data,
                             m->profs["u2_turb"].data,  m->profs["v2_turb"].data,  m->profs["w2_turb"].data, m->profs["tke_turb"].data,
                             m->profs["uw_turb"].data, m->profs["vw_turb"].data,
                             fields.u->data, fields.v->data, fields.w->data, umodel, vmodel,
                             fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, grid.dzi, grid.dzhi);
    }

    // Calculate the pressure transport and redistribution terms
    calc_pressure_terms(m->profs["w2_pres"].data,  m->profs["tke_pres"].data,
                        m->profs["uw_pres"].data, m->profs["vw_pres"].data,
                        m->profs["u2_rdstr"].data, m->profs["v2_rdstr"].data, m->profs["w2_rdstr"].data,
                        m->profs["uw_rdstr"].data, m->profs["vw_rdstr"].data,
                        fields.u->data, fields.v->data, fields.w->data, fields.sd["p"]->data, umodel, vmodel,
                        grid.dzi, grid.dzhi, grid.dxi, grid.dyi);

    if(diff.get_switch() != "0")
    {
        // Calculate the diffusive transport and dissipation terms
        if(diff.get_switch() == "2" || diff.get_switch() == "4")
            calc_diffusion_terms_DNS(m->profs["u2_visc"].data, m->profs["v2_visc"].data, m->profs["w2_visc"].data, m->profs["tke_visc"].data, m->profs["uw_visc"].data,
                                     m->profs["u2_diss"].data, m->profs["v2_diss"].data, m->profs["w2_diss"].data, m->profs["tke_diss"].data, m->profs["uw_diss"].data,
                                     fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, fields.atmp["tmp3"]->data, fields.u->data, fields.v->data, fields.w->data, umodel, vmodel,
                                     grid.dzi, grid.dzhi, grid.dxi, grid.dyi, fields.visc);
        else if(diff.get_switch() == "smag2")
            calc_diffusion_terms_LES(m->profs["u2_visc"].data, m->profs["v2_visc"].data, m->profs["w2_visc"].data, m->profs["tke_visc"].data, m->profs["uw_visc"].data,
                                     m->profs["u2_diss"].data, m->profs["v2_diss"].data, m->profs["w2_diss"].data, m->profs["tke_diss"].data, m->profs["uw_diss"].data,
                                     m->profs["vw_diss"].data, fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, fields.atmp["tmp3"]->data,
                                     fields.u->data, fields.v->data, fields.w->data,
                                     fields.u->datafluxbot, fields.v->datafluxbot,
                                     fields.sd["evisc"]->data, umodel, vmodel,
                                     grid.dzi, grid.dzhi, grid.dxi, grid.dyi);
    }

    if(thermo.get_switch() != "0")
    {
        // Store the buoyancy in the tmp1 field
        thermo.get_thermo_field(fields.atmp["tmp1"], fields.atmp["tmp2"], "b");

        // Calculate mean fields
        grid.calc_mean(fields.atmp["tmp1"]->datamean, fields.atmp["tmp1"]->data, grid.kcells);
        grid.calc_mean(fields.sd["p"]->datamean, fields.sd["p"]->data, grid.kcells);

        // Calculate buoyancy terms
        calc_buoyancy_terms(m->profs["w2_buoy"].data, m->profs["tke_buoy"].data,
                            m->profs["uw_buoy"].data, m->profs["vw_buoy"].data,
                            fields.u->data, fields.v->data, fields.w->data, fields.atmp["tmp1"]->data,
                            umodel, vmodel, fields.atmp["tmp1"]->datamean);
    }
}

namespace
{
    // Double linear interpolation
    inline double interp2_4(const double a, const double b, const double c, const double d)
    {
        return 0.25 * (a + b + c + d);
    }
}

/**
 * Calculate the kinetic and turbulence kinetic energy
 * @param TO-DO
 */
void Budget_2::calc_kinetic_energy(double* const restrict ke, double* const restrict tke,
                                   const double* const restrict u, const double* const restrict v, const double* const restrict w,
                                   const double* const restrict umodel, const double* const restrict vmodel,
                                   const double utrans, const double vtrans)
{
    const int ii = 1;
    const int jj = grid.icells;
    const int kk = grid.ijcells;
    const int ijtot = grid.itot*grid.jtot;

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        ke[k]  = 0;
        tke[k] = 0;
    }

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                const double u2 = pow(interp2(u[ijk]+utrans, u[ijk+ii]+utrans), 2);
                const double v2 = pow(interp2(v[ijk]+vtrans, v[ijk+jj]+vtrans), 2);
                const double w2 = pow(interp2(w[ijk]       , w[ijk+ii]       ), 2);

                ke[k] += 0.5 * (u2 + v2 + w2);
            }

        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                const double u2 = pow(interp2(u[ijk]-umodel[k], u[ijk+ii]-umodel[k]), 2);
                const double v2 = pow(interp2(v[ijk]-vmodel[k], v[ijk+jj]-vmodel[k]), 2);
                const double w2 = pow(interp2(w[ijk]          , w[ijk+ii]          ), 2);

                tke[k] += 0.5 * (u2 + v2 + w2);
            }
    }

    // Calculate sum over all processes, and calc mean profiles
    master.sum(ke , grid.kcells);
    master.sum(tke, grid.kcells);

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
       ke[k]  /= ijtot;
       tke[k] /= ijtot;
    }
}

/**
 * Calculate the budget terms arrising from the advection term:
 * shear production (-2 u_i*u_j * d<u_i>/dx_j) and turbulent transport (-d(u_i^2*u_j)/dx_j)
 * @param TO-DO
 */
void Budget_2::calc_advection_terms(double* const restrict u2_shear, double* const restrict v2_shear,
                                    double* const restrict tke_shear,
                                    double* const restrict uw_shear, double* const restrict vw_shear,
                                    double* const restrict u2_turb,  double* const restrict v2_turb,
                                    double* const restrict w2_turb, double* const restrict tke_turb,
                                    double* const restrict uw_turb, double* const restrict vw_turb,
                                    const double* const restrict u, const double* const restrict v, const double* const restrict w,
                                    const double* const restrict umean, const double* const restrict vmean,
                                    double* const restrict wx, double* const restrict wy,
                                    const double* const restrict dzi, const double* const restrict dzhi)
{
    // Interpolate the vertical velocity to {xh,y,zh} (wx, below u) and {x,yh,zh} (wy, below v)
    const int wloc [3] = {0,0,1};
    const int wxloc[3] = {1,0,1};
    const int wyloc[3] = {0,1,1};

    grid.interpolate_2nd(wx, w, wloc, wxloc);
    grid.interpolate_2nd(wy, w, wloc, wyloc);

    const int ii = 1;
    const int jj = grid.icells;
    const int kk = grid.ijcells;
    const int ijtot = grid.itot * grid.jtot;

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        u2_shear [k] = 0;
        v2_shear [k] = 0;
        tke_shear[k] = 0;
        u2_turb  [k] = 0;
        v2_turb  [k] = 0;
        tke_turb [k] = 0;
    }

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        uw_shear [k] = 0;
        vw_shear [k] = 0;
        w2_turb  [k] = 0;
        uw_turb  [k] = 0;
        vw_turb  [k] = 0;
    }

    // Calculate shear terms (-2u_iw d<u_i>/dz)
    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        const double dudz = (interp2(umean[k], umean[k+1]) - interp2(umean[k-1], umean[k]) ) * dzi[k];
        const double dvdz = (interp2(vmean[k], vmean[k+1]) - interp2(vmean[k-1], vmean[k]) ) * dzi[k];

        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                u2_shear[k] -= 2 * (u[ijk]-umean[k]) * interp2(wx[ijk], wx[ijk+kk]) * dudz;

                v2_shear[k] -= 2 * (v[ijk]-vmean[k]) * interp2(wy[ijk], wy[ijk+kk]) * dvdz;

                uw_shear[k] -= pow(wx[ijk], 2) * (umean[k] - umean[k-1]) * dzhi[k];

                vw_shear[k] -= pow(wy[ijk], 2) * (vmean[k] - vmean[k-1]) * dzhi[k];
            }

        tke_shear[k] += 0.5*(u2_shear[k] + v2_shear[k]);
    }

    // Calculate turbulent transport terms (-d(u_i^2*w)/dz)
    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                u2_turb[k]  -= ( pow(interp2(u[ijk]-umean[k], u[ijk+kk]-umean[k+1]), 2) * wx[ijk+kk] -
                                 pow(interp2(u[ijk]-umean[k], u[ijk-kk]-umean[k-1]), 2) * wx[ijk   ] ) * dzi[k];

                v2_turb[k]  -= ( pow(interp2(v[ijk]-vmean[k], v[ijk+kk]-vmean[k+1]), 2) * wy[ijk+kk] -
                                 pow(interp2(v[ijk]-vmean[k], v[ijk-kk]-vmean[k-1]), 2) * wy[ijk   ] ) * dzi[k];

                tke_turb[k] -= 0.5 * ( pow(w[ijk+kk], 3) - pow(w[ijk], 3) ) * dzi[k];
            }
        tke_turb[k] += 0.5 * (u2_turb[k] + v2_turb[k]);
    }

    // Lower boundary kstart (z=0)
    int k = grid.kstart;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj + k*kk;

            // w^3 @ full level below sfc == -w^3 @ full level above sfc
            w2_turb[k] -= 2 * pow(interp2(w[ijk], w[ijk+kk]), 3) * dzhi[k];

            // w^2 @ full level below sfc == w^2 @ full level above sfc
            uw_turb[k] -= ( (u[ijk]   -umean[k  ]) * pow(interp2(wx[ijk], wx[ijk+kk]), 2) -
                            (u[ijk-kk]-umean[k-1]) * pow(interp2(wx[ijk], wx[ijk-kk]), 2) ) * dzhi[k];

            vw_turb[k] -= ( (v[ijk]   -vmean[k  ]) * pow(interp2(wy[ijk], wy[ijk+kk]), 2) -
                            (v[ijk-kk]-vmean[k-1]) * pow(interp2(wy[ijk], wy[ijk-kk]), 2) ) * dzhi[k];
        }

    // Top boundary kstart (z=zsize)
    k = grid.kend;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj + k*kk;

            // w^3 @ full level above top == -w^3 @ full level below top
            w2_turb[k] -= -2 * pow(interp2(w[ijk], w[ijk-kk]), 3) * dzhi[k];

            // w^2 @ full level above top == w^2 @ full level below top
            uw_turb[k] -= ( (u[ijk]   -umean[k  ]) * pow(interp2(wx[ijk], wx[ijk-kk]), 2) -
                            (u[ijk-kk]-umean[k-1]) * pow(interp2(wx[ijk], wx[ijk-kk]), 2) ) * dzhi[k];

            // w^2 @ full level above top == w^2 @ full level below top
            vw_turb[k] -= ( (v[ijk]   -vmean[k  ]) * pow(interp2(wy[ijk], wy[ijk-kk]), 2) -
                            (v[ijk-kk]-vmean[k-1]) * pow(interp2(wy[ijk], wy[ijk-kk]), 2) ) * dzhi[k];
        }

    // Inner domain
    for (int k=grid.kstart+1; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                w2_turb[k] -= ( pow(interp2(w[ijk], w[ijk+kk]), 3) -
                                pow(interp2(w[ijk], w[ijk-kk]), 3) ) * dzhi[k];

                uw_turb[k] -= ( (u[ijk]   -umean[k  ]) * pow(interp2(wx[ijk], wx[ijk+kk]), 2) -
                                (u[ijk-kk]-umean[k-1]) * pow(interp2(wx[ijk], wx[ijk-kk]), 2) ) * dzhi[k];

                vw_turb[k] -= ( (v[ijk]   -vmean[k  ]) * pow(interp2(wy[ijk], wy[ijk+kk]), 2) -
                                (v[ijk-kk]-vmean[k-1]) * pow(interp2(wy[ijk], wy[ijk-kk]), 2) ) * dzhi[k];
            }

    // Calculate sum over all processes, and calc mean profiles
    master.sum(u2_shear,  grid.kcells);
    master.sum(v2_shear,  grid.kcells);
    master.sum(tke_shear, grid.kcells);
    master.sum(uw_shear,  grid.kcells);
    master.sum(vw_shear,  grid.kcells);
    master.sum(u2_turb,   grid.kcells);
    master.sum(v2_turb,   grid.kcells);
    master.sum(w2_turb,   grid.kcells);
    master.sum(tke_turb,  grid.kcells);
    master.sum(uw_turb,   grid.kcells);
    master.sum(vw_turb,   grid.kcells);

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        u2_shear [k] /= ijtot;
        v2_shear [k] /= ijtot;
        tke_shear[k] /= ijtot;
        u2_turb  [k] /= ijtot;
        v2_turb  [k] /= ijtot;
        tke_turb [k] /= ijtot;
    }

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        uw_shear [k] /= ijtot;
        vw_shear [k] /= ijtot;
        w2_turb  [k] /= ijtot;
        uw_turb  [k] /= ijtot;
        vw_turb  [k] /= ijtot;
    }
}

/**
 * Calculate the budget terms arrising from pressure:
 * pressure transport (-2*dpu_i/dxi) and redistribution (2p*dui/dxi)
 * @param TO-DO
 */
void Budget_2::calc_pressure_terms(double* const restrict w2_pres,  double* const restrict tke_pres,
                                   double* const restrict uw_pres,  double* const restrict vw_pres,
                                   double* const restrict u2_rdstr, double* const restrict v2_rdstr, double* const restrict w2_rdstr,
                                   double* const restrict uw_rdstr, double* const restrict vw_rdstr,
                                   const double* const restrict u, const double* const restrict v,
                                   const double* const restrict w, const double* const restrict p,
                                   const double* const restrict umean, const double* const restrict vmean,
                                   const double* const restrict dzi, const double* const restrict dzhi, const double dxi, const double dyi)
{
    const int ii = 1;
    const int jj = grid.icells;
    const int kk = grid.ijcells;
    const int ijtot = grid.itot * grid.jtot;

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        tke_pres[k] = 0;
        u2_rdstr[k] = 0;
        v2_rdstr[k] = 0;
    }

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        w2_pres [k] = 0;
        uw_pres [k] = 0;
        vw_pres [k] = 0;
        w2_rdstr[k] = 0;
        uw_rdstr[k] = 0;
        vw_rdstr[k] = 0;
    }

    // Pressure transport term (-2*dpu_i/dxi)
    for (int k=grid.kstart; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                tke_pres[k] -= ( interp2(p[ijk], p[ijk+kk]) * w[ijk+kk] -
                                 interp2(p[ijk], p[ijk-kk]) * w[ijk   ] ) * dzi[k];

                uw_pres[k]  -= ( interp2(p[ijk   ], p[ijk-kk   ]) * w[ijk    ] -
                                 interp2(p[ijk-ii], p[ijk-ii-kk]) * w[ijk-ii ] ) * dxi +
                               ( interp2(p[ijk   ], p[ijk-ii   ]) * (u[ijk   ]-umean[k  ]) -
                                 interp2(p[ijk-kk], p[ijk-ii-kk]) * (u[ijk-kk]-umean[k-1]) ) * dzhi[k];

                vw_pres[k]  -= ( interp2(p[ijk-kk   ], p[ijk   ]) * w[ijk    ]  -
                                 interp2(p[ijk-jj-kk], p[ijk-jj]) * w[ijk-jj ] ) * dyi +
                               ( interp2(p[ijk-jj   ], p[ijk   ]) * (v[ijk   ]-vmean[k  ]) -
                                 interp2(p[ijk-jj-kk], p[ijk-kk]) * (v[ijk-kk]-vmean[k-1]) ) * dzhi[k];
            }

    // Lower boundary (z=0)
    int k = grid.kstart;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj + k*kk;

            // w @ full level below sfc == -w @ full level above sfc
            w2_pres[k] -= 2 * ( interp2(w[ijk], w[ijk+kk]) * p[ijk   ] -
                              - interp2(w[ijk], w[ijk+kk]) * p[ijk-kk] ) * dzhi[k];
        }

    // Top boundary (z=zsize)
    // TODO: what to do with w2_pres and uw_pres at the top boundary? Pressure at k=kend is undefined?

    // Inner domain
    for (int k=grid.kstart+1; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                w2_pres[k] -= 2 * ( interp2(w[ijk], w[ijk+kk]) * p[ijk   ] -
                                    interp2(w[ijk], w[ijk-kk]) * p[ijk-kk] ) * dzhi[k];
            }

    // Pressure redistribution term (2p*dui/dxi)
    for (int k=grid.kstart; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                u2_rdstr[k] += 2 * interp2(p[ijk], p[ijk-ii]) *
                                 ( interp2(u[ijk]-umean[k], u[ijk+ii]-umean[k]) -
                                   interp2(u[ijk]-umean[k], u[ijk-ii]-umean[k]) ) * dxi;

                v2_rdstr[k] += 2 * interp2(p[ijk], p[ijk-jj]) *
                                 ( interp2(v[ijk]-vmean[k], v[ijk+jj]-vmean[k]) -
                                   interp2(v[ijk]-vmean[k], v[ijk-jj]-vmean[k]) ) * dyi;

                uw_rdstr[k] += interp2_4(p[ijk], p[ijk-kk], p[ijk-ii-kk], p[ijk-ii]) *
                                 ( ((u[ijk]-umean[k]) - (u[ijk-kk]-umean[k-1])) * dzhi[k] + (w[ijk] - w[ijk-ii]) * dxi );

                vw_rdstr[k] += interp2_4(p[ijk], p[ijk-kk], p[ijk-jj-kk], p[ijk-jj]) *
                                 ( ((v[ijk]-vmean[k]) - (v[ijk-kk]-vmean[k-1])) * dzhi[k] + (w[ijk] - w[ijk-jj]) * dyi );
            }

    // Lower boundary (z=0)
    k = grid.kstart;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj + k*kk;

            // with w[kstart] == 0, dw/dz at surface equals (w[kstart+1] - w[kstart]) / dzi
            w2_rdstr[k] += 2 * interp2(p[ijk], p[ijk-kk]) * (w[ijk+kk] - w[ijk]) * dzi[k];
        }

    // Top boundary (z=zsize)
    // TODO: what to do with w2_rdstr at the top boundary? Pressure at k=kend is undefined?

    for (int k=grid.kstart+1; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                w2_rdstr[k] += 2 * interp2(p[ijk], p[ijk-kk]) *
                                 ( interp2(w[ijk], w[ijk+kk]) - interp2(w[ijk], w[ijk-kk]) ) * dzhi[k];
            }

    // Calculate sum over all processes, and calc mean profiles
    master.sum(w2_pres,  grid.kcells);
    master.sum(tke_pres, grid.kcells);
    master.sum(uw_pres,  grid.kcells);
    master.sum(vw_pres,  grid.kcells);
    master.sum(u2_rdstr, grid.kcells);
    master.sum(v2_rdstr, grid.kcells);
    master.sum(w2_rdstr, grid.kcells);
    master.sum(uw_rdstr, grid.kcells);
    master.sum(vw_rdstr, grid.kcells);

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        tke_pres[k] /= ijtot;
        u2_rdstr[k] /= ijtot;
        v2_rdstr[k] /= ijtot;
    }

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        w2_pres [k] /= ijtot;
        uw_pres [k] /= ijtot;
        vw_pres [k] /= ijtot;
        w2_rdstr[k] /= ijtot;
        uw_rdstr[k] /= ijtot;
        vw_rdstr[k] /= ijtot;
    }
}

/**
 * Calculate the budget terms arrising from diffusion, for a variable viscosity
 * molecular diffusion (nu*d/dxj(dui^2/dxj)) and dissipation (-2*nu*(dui/dxj)^2)
 * @param TO-DO
 */
void Budget_2::calc_diffusion_terms_LES(double* const restrict u2_visc, double* const restrict v2_visc,
                                        double* const restrict w2_visc, double* const restrict tke_visc, double* const restrict uw_visc,
                                        double* const restrict u2_diss, double* const restrict v2_diss,
                                        double* const restrict w2_diss, double* const restrict tke_diss,
                                        double* const restrict uw_diss, double* const restrict vw_diss,
                                        double* const restrict wz, double* const restrict evisch, double* const restrict wy,
                                        const double* const restrict u, const double* const restrict v,
                                        const double* const restrict w,
                                        const double* const restrict ufluxbot, const double* const restrict vfluxbot,
                                        const double* const restrict evisc,
                                        const double* const restrict umean, const double* const restrict vmean,
                                        const double* const restrict dzi, const double* const restrict dzhi,
                                        const double dxi, const double dyi)
{
    const int ii = 1;
    const int ii2 = 2;
    const int jj = grid.icells;
    const int jj2 = 2*grid.icells;
    const int kk = grid.ijcells;
    const int kk2 = 2*grid.ijcells;
    const int ijtot = grid.itot * grid.jtot;

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        u2_visc [k] = 0;
        v2_visc [k] = 0;
        tke_visc[k] = 0;
        u2_diss [k] = 0;
        v2_diss [k] = 0;
        tke_diss[k] = 0;
    }

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        w2_visc [k] = 0;
        uw_visc [k] = 0;
        w2_diss [k] = 0;
        uw_diss [k] = 0;
        vw_diss [k] = 0;
    }

    // Calculate w at full levels (grid center)
    for (int k=grid.kstart; k<grid.kend; ++k)
        for (int j=0; j<grid.jcells; ++j)
            #pragma ivdep
            for (int i=0; i<grid.icells; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                wz[ijk] = interp2(w[ijk], w[ijk+kk]);
            }

    // Set ghost cells such that the velocity interpolated to the boundaries is zero
    int ks = grid.kstart;
    int ke = grid.kend-1;
    for (int j=0; j<grid.jcells; ++j)
        #pragma ivdep
        for (int i=0; i<grid.icells; ++i)
        {
            const int ijks = i + j*jj + ks*kk;
            const int ijke = i + j*jj + ke*kk;
            wz[ijks-kk] = -wz[ijks];
            wz[ijke+kk] = -wz[ijke];
        }

    // Calculate evisc at half-half-half level
    for (int k=grid.kstart; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                evisch[ijk] = 0.125 * (evisc[ijk-ii-jj-kk] + evisc[ijk-ii-jj] + evisc[ijk-ii-kk] + evisc[ijk-ii] +
                                       evisc[ijk   -jj-kk] + evisc[ijk   -jj] + evisc[ijk   -kk] + evisc[ijk   ]);
            }
    grid.boundary_cyclic(evisch);

    if(true)
    {
        // -----------------------------
        // Test: directly calculate diffusion terms as 2 ui * d/dxj(visc * dui/dx + visc * duj/dxi)
        // Term is stored in xx_diss; xx_visc=0
        // -----------------------------
        for (int k=grid.kstart; k<grid.kend; ++k)
            for (int j=grid.jstart; j<grid.jend; ++j)
                #pragma ivdep
                for (int i=grid.istart; i<grid.iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    const double evisc_utop   = interp2_4(evisc[ijk], evisc[ijk+kk], evisc[ijk-ii+kk], evisc[ijk-ii]);
                    const double evisc_ubot   = interp2_4(evisc[ijk], evisc[ijk-kk], evisc[ijk-ii-kk], evisc[ijk-ii]);
                    const double evisc_unorth = interp2_4(evisc[ijk], evisc[ijk+jj], evisc[ijk+jj-ii], evisc[ijk-ii]);
                    const double evisc_usouth = interp2_4(evisc[ijk], evisc[ijk-jj], evisc[ijk-jj-ii], evisc[ijk-ii]);

                    const double evisc_vtop   = interp2_4(evisc[ijk], evisc[ijk+kk], evisc[ijk-jj+kk], evisc[ijk-jj]);
                    const double evisc_vbot   = interp2_4(evisc[ijk], evisc[ijk-kk], evisc[ijk-jj-kk], evisc[ijk-jj]);
                    const double evisc_veast  = interp2_4(evisc[ijk], evisc[ijk+ii], evisc[ijk+ii-jj], evisc[ijk-jj]);
                    const double evisc_vwest  = interp2_4(evisc[ijk], evisc[ijk-ii], evisc[ijk-ii-jj], evisc[ijk-jj]);

                    const double evisc_weast  = interp2_4(evisc[ijk], evisc[ijk+ii], evisc[ijk+ii-kk], evisc[ijk-kk]);
                    const double evisc_wwest  = interp2_4(evisc[ijk], evisc[ijk-ii], evisc[ijk-ii-kk], evisc[ijk-kk]);
                    const double evisc_wnorth = interp2_4(evisc[ijk], evisc[ijk+jj], evisc[ijk+jj-kk], evisc[ijk-kk]);
                    const double evisc_wsouth = interp2_4(evisc[ijk], evisc[ijk-jj], evisc[ijk-jj-kk], evisc[ijk-kk]);

                    // -----------------------------------------
                    // 2 * u * d/dx( visc * du/dx + visc * du/dx )
                    u2_diss[k] += 2 * (u[ijk]-umean[k]) * ( evisc[ijk   ] * (u[ijk+ii] - u[ijk   ]) * dxi -
                                                            evisc[ijk-ii] * (u[ijk   ] - u[ijk-ii]) * dxi ) * 2 * dxi;

                    // 2 * u * d/dy( visc * du/dy + visc * dv/dx)
                    u2_diss[k] += 2 * (u[ijk]-umean[k]) * ( evisc_unorth * (u[ijk+jj] - u[ijk      ]) * dyi -
                                                            evisc_usouth * (u[ijk   ] - u[ijk-jj   ]) * dyi +
                                                            evisc_unorth * (v[ijk+jj] - v[ijk+jj-ii]) * dxi -
                                                            evisc_usouth * (v[ijk   ] - v[ijk-ii   ]) * dxi ) * dyi;

                    // 2 * u * d/dz( visc * dw/dx )
                    u2_diss[k] += 2 * (u[ijk]-umean[k]) * ( evisc_utop * (w[ijk+kk] - w[ijk-ii+kk]) * dxi -
                                                            evisc_ubot * (w[ijk   ] - w[ijk-ii   ]) * dxi ) * dzi[k];

                    // -----------------------------------------
                    // 2 * v * d/dy( visc * dv/dy + visc * dv/dy )
                    v2_diss[k] += 2 * (v[ijk]-vmean[k]) * ( evisc[ijk   ] * (v[ijk+jj] - v[ijk   ]) * dyi -
                                                            evisc[ijk-jj] * (v[ijk   ] - v[ijk-jj]) * dyi ) * 2 * dyi;

                    // 2 * v * d/dx( visc * dv/dx + visc * du/dy )
                    v2_diss[k] += 2 * (v[ijk]-vmean[k]) * ( evisc_veast * (v[ijk+ii] - v[ijk      ]) * dxi -
                                                            evisc_vwest * (v[ijk   ] - v[ijk-ii   ]) * dxi +
                                                            evisc_veast * (u[ijk+ii] - u[ijk+ii-jj]) * dyi -
                                                            evisc_vwest * (u[ijk   ] - u[ijk-jj   ]) * dyi ) * dxi;

                    // 2 * v * d/dz( visc * dw/dy )
                    v2_diss[k] += 2 * (v[ijk]-vmean[k]) * ( evisc_vtop * (w[ijk+kk] - w[ijk-jj+kk]) * dyi -
                                                            evisc_vbot * (w[ijk   ] - w[ijk-jj   ]) * dyi ) * dzi[k];

                    // -----------------------------------------
                    // 2 * w * d/dx( visc * dw/dx )
                    w2_diss[k] += 2 * w[ijk] * ( evisc_weast * (w[ijk+ii] - w[ijk   ]) * dxi -
                                                 evisc_wwest * (w[ijk   ] - w[ijk-ii]) * dxi ) * dxi;

                    // 2 * w * d/dy( visc * dw/dy )
                    w2_diss[k] += 2 * w[ijk] * ( evisc_wnorth * (w[ijk+jj] - w[ijk   ]) * dyi -
                                                 evisc_wsouth * (w[ijk   ] - w[ijk-jj]) * dyi ) * dyi;

                    // -----------------------------------------
                    // 2 * w * d/dx( visc * dw/dx )
                    tke_diss[k] += wz[ijk] * ( interp2(evisc[ijk], evisc[ijk+ii]) * (wz[ijk+ii] - wz[ijk   ]) * dxi -
                                               interp2(evisc[ijk], evisc[ijk-ii]) * (wz[ijk   ] - wz[ijk-ii]) * dxi ) * dxi;

                    // 2 * w * d/dx( visc * du/dz )
                    tke_diss[k] += wz[ijk] * ( interp2(evisc[ijk], evisc[ijk+ii]) * (interp2(u[ijk+ii], u[ijk+ii+kk]) - interp2(u[ijk+ii], u[ijk+ii-kk])) * dzi[k] -
                                               interp2(evisc[ijk], evisc[ijk-ii]) * (interp2(u[ijk   ], u[ijk   +kk]) - interp2(u[ijk   ], u[ijk   -kk])) * dzi[k] ) * dxi;

                    // 2 * w * d/dy( visc * dw/dy )
                    tke_diss[k] += wz[ijk] * ( interp2(evisc[ijk], evisc[ijk+jj]) * (wz[ijk+jj] - wz[ijk   ]) * dyi -
                                               interp2(evisc[ijk], evisc[ijk-jj]) * (wz[ijk   ] - wz[ijk-jj]) * dyi ) * dyi;

                    // 2 * w * d/dy( visc * dv/dz )
                    tke_diss[k] += wz[ijk] * ( interp2(evisc[ijk], evisc[ijk+jj]) * (interp2(v[ijk+jj], v[ijk+jj+kk]) - interp2(v[ijk+jj], v[ijk+jj-kk])) * dzi[k] -
                                               interp2(evisc[ijk], evisc[ijk-jj]) * (interp2(v[ijk   ], v[ijk   +kk]) - interp2(v[ijk   ], v[ijk   -kk])) * dzi[k] ) * dyi;
                }

        for (int k=grid.kstart+1; k<grid.kend; ++k)
        {
            for (int j=grid.jstart; j<grid.jend; ++j)
                #pragma ivdep
                for (int i=grid.istart; i<grid.iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    const double evisc_utop   = interp2_4(evisc[ijk], evisc[ijk+kk], evisc[ijk-ii+kk], evisc[ijk-ii]);
                    const double evisc_ubot   = interp2_4(evisc[ijk], evisc[ijk-kk], evisc[ijk-ii-kk], evisc[ijk-ii]);

                    const double evisc_vtop   = interp2_4(evisc[ijk], evisc[ijk+kk], evisc[ijk-jj+kk], evisc[ijk-jj]);
                    const double evisc_vbot   = interp2_4(evisc[ijk], evisc[ijk-kk], evisc[ijk-jj-kk], evisc[ijk-jj]);

                    const double evisc_weast  = interp2_4(evisc[ijk], evisc[ijk+ii], evisc[ijk+ii-kk], evisc[ijk-kk]);
                    const double evisc_wwest  = interp2_4(evisc[ijk], evisc[ijk-ii], evisc[ijk-ii-kk], evisc[ijk-kk]);
                    const double evisc_wnorth = interp2_4(evisc[ijk], evisc[ijk+jj], evisc[ijk+jj-kk], evisc[ijk-kk]);
                    const double evisc_wsouth = interp2_4(evisc[ijk], evisc[ijk-jj], evisc[ijk-jj-kk], evisc[ijk-kk]);

                    // -----------------------------------------
                    // 2 * u * d/dz( visc * du/dz )
                    u2_diss[k] += 2 * (u[ijk]-umean[k]) * ( evisc_utop * (u[ijk+kk] - u[ijk   ]) * dzhi[k+1] -
                                                            evisc_ubot * (u[ijk   ] - u[ijk-kk]) * dzhi[k  ] ) * dzi[k];

                    // -----------------------------------------
                    // 2 * v * d/dz( visc * dv/dz )
                    v2_diss[k] += 2 * (v[ijk]-vmean[k]) * ( evisc_vtop * (v[ijk+kk] - v[ijk   ]) * dzhi[k+1] -
                                                            evisc_vbot * (v[ijk   ] - v[ijk-kk]) * dzhi[k  ] ) * dzi[k];

                    // -----------------------------------------
                    // 2 * w * d/dx( visc * du/dz )
                    w2_diss[k] += 2 * w[ijk] * ( evisc_weast * (u[ijk+ii] - u[ijk+ii-kk]) * dzhi[k] -
                                                 evisc_wwest * (u[ijk   ] - u[ijk   -kk]) * dzhi[k] ) * dxi;

                    // 2 * w * d/dy( visc * dv/dz )
                    w2_diss[k] += 2 * w[ijk] * ( evisc_wnorth * (v[ijk+jj] - v[ijk+jj-kk]) * dzhi[k] -
                                                 evisc_wsouth * (v[ijk   ] - v[ijk   -kk]) * dzhi[k] ) * dyi;

                    // 2 * w * d/dz( visc * dw/dz )
                    w2_diss[k] += 2 * w[ijk] * ( evisc[ijk   ] * (w[ijk+kk] - w[ijk   ]) * dzi[k  ] -
                                                 evisc[ijk-kk] * (w[ijk   ] - w[ijk-kk]) * dzi[k-1] ) * 2 * dzhi[k];

                    // -----------------------------------------
                    // 2 * w * d/dz( 2 * visc * dw/dz )
                    tke_diss[k] += wz[ijk] * ( interp2(evisc[ijk], evisc[ijk+kk]) * (wz[ijk+kk] - wz[ijk   ]) * dzhi[k+1] -
                                               interp2(evisc[ijk], evisc[ijk-kk]) * (wz[ijk   ] - wz[ijk-kk]) * dzhi[k  ] ) * 2 * dzi[k];

                    // -----------------------------------------
                    // w * d/dx(visc * du/dx + visc * du/dx)
                    uw_diss[k] += ( ( interp2(w[ijk-ii], w[ijk    ])
                                      * ( ( ( ( 2 * interp2(evisc[ijk    -kk], evisc[ijk        ]) )
                                          * ( interp2(u[ijk+ii-kk], u[ijk+ii    ]) - interp2(u[ijk    -kk], u[ijk        ]) ) )
                                        * dxi ) - ( ( ( 2 * interp2(evisc[ijk-ii-kk], evisc[ijk-ii    ]) )
                                          * ( interp2(u[ijk    -kk], u[ijk        ]) - interp2(u[ijk-ii-kk], u[ijk-ii    ]) ) )
                                        * dxi ) ) )
                                    * dxi );

                    // w * d/dy(visc * du/dy + visc * dv/dx)
                    uw_diss[k] += ( ( interp2(w[ijk-ii], w[ijk    ])
                                      * ( ( evisch[ijk+jj]
                                        * ( ( ( interp2(u[ijk+jj-kk], u[ijk+jj    ]) - interp2(u[ijk    -kk], u[ijk        ]) )
                                            * dyi )
                                          + ( ( interp2(v[ijk    +jj-kk], v[ijk    +jj    ]) - interp2(v[ijk-ii+jj-kk], v[ijk-ii+jj    ]) )
                                            * dxi ) ) ) - ( evisch[ijk    ]
                                        * ( ( ( interp2(u[ijk    -kk], u[ijk        ]) - interp2(u[ijk-jj-kk], u[ijk-jj    ]) )
                                            * dyi )
                                          + ( ( interp2(v[ijk        -kk], v[ijk            ]) - interp2(v[ijk-ii    -kk], v[ijk-ii        ]) )
                                            * dxi ) ) ) ) )
                                    * dyi );

                    // w * d/dz(visc * du/dz + visc * dw/dx)
                    uw_diss[k] += ( ( interp2(w[ijk-ii], w[ijk    ])
                                      * ( ( interp2(evisc[ijk-ii    ], evisc[ijk        ])
                                        * ( ( ( interp2(u[ijk    ], u[ijk+kk]) - interp2(u[ijk-kk], u[ijk    ]) )
                                            * dzi[k  ] )
                                          + ( ( interp2(w[ijk        ], w[ijk    +kk]) - interp2(w[ijk-ii    ], w[ijk-ii+kk]) )
                                            * dxi ) ) ) - ( interp2(evisc[ijk-ii-kk], evisc[ijk    -kk])
                                        * ( ( ( interp2(u[ijk-kk], u[ijk    ]) - interp2(u[ijk-kk2], u[ijk-kk]) )
                                            * dzi[k-1] )
                                          + ( ( interp2(w[ijk    -kk], w[ijk        ]) - interp2(w[ijk-ii-kk], w[ijk-ii    ]) )
                                            * dxi ) ) ) ) )
                                    * dzhi[k] );

                    // u * d/dx(visc * dw/dx + visc * du/dz)
                    uw_diss[k] += ( ( interp2(u[ijk-kk], u[ijk    ])
                                      * ( ( interp2(evisc[ijk    -kk], evisc[ijk        ])
                                        * ( ( ( interp2(w[ijk    ], w[ijk+ii]) - interp2(w[ijk-ii], w[ijk    ]) )
                                            * dxi )
                                          + ( ( interp2(u[ijk        ], u[ijk+ii    ]) - interp2(u[ijk    -kk], u[ijk+ii-kk]) )
                                            * dzhi[k] ) ) ) - ( interp2(evisc[ijk-ii-kk], evisc[ijk-ii    ])
                                        * ( ( ( interp2(w[ijk-ii], w[ijk    ]) - interp2(w[ijk-ii2], w[ijk-ii]) )
                                            * dxi )
                                          + ( ( interp2(u[ijk-ii    ], u[ijk        ]) - interp2(u[ijk-ii-kk], u[ijk    -kk]) )
                                            * dzhi[k] ) ) ) ) )
                                    * dxi );

                    // u * d/dy(visc * dw/dy + visc * dv/dz)
                    uw_diss[k] += ( ( interp2(u[ijk-kk], u[ijk    ])
                                      * ( ( evisch[ijk+jj]
                                        * ( ( ( interp2(w[ijk-ii+jj], w[ijk    +jj]) - interp2(w[ijk-ii    ], w[ijk        ]) )
                                            * dyi )
                                          + ( ( interp2(v[ijk-ii+jj    ], v[ijk    +jj    ]) - interp2(v[ijk-ii+jj-kk], v[ijk    +jj-kk]) )
                                            * dzhi[k] ) ) ) - ( evisch[ijk    ]
                                        * ( ( ( interp2(w[ijk-ii    ], w[ijk        ]) - interp2(w[ijk-ii-jj], w[ijk    -jj]) )
                                            * dyi )
                                          + ( ( interp2(v[ijk-ii        ], v[ijk            ]) - interp2(v[ijk-ii    -kk], v[ijk        -kk]) )
                                            * dzhi[k] ) ) ) ) )
                                    * dyi );

                    // u * d/dz(visc * dw/dz + visc * dw/dz)
                    uw_diss[k] += ( ( interp2(u[ijk-kk], u[ijk    ])
                                      * ( ( ( 2 * interp2(evisc[ijk-ii    ], evisc[ijk        ]) )
                                        * ( ( interp2(w[ijk-ii+kk], w[ijk    +kk]) - interp2(w[ijk-ii    ], w[ijk        ]) )
                                          * dzi[k  ] ) ) - ( ( 2 * interp2(evisc[ijk-ii-kk], evisc[ijk    -kk]) )
                                        * ( ( interp2(w[ijk-ii    ], w[ijk        ]) - interp2(w[ijk-ii-kk], w[ijk    -kk]) )
                                          * dzi[k-1] ) ) ) )
                                    * dzhi[k] );

                    // ------------------------------------------------
                    // w * d/dx(visc * dv/dx + visc * du/dy)
                    vw_diss[k] += ( ( interp2(w[ijk-jj], w[ijk    ])
                                    * ( ( evisch[ijk+ii]
                                      * ( ( ( interp2(v[ijk+ii-kk], v[ijk+ii    ]) - interp2(v[ijk    -kk], v[ijk        ]) )
                                          * dxi )
                                        + ( ( interp2(u[ijk+ii    -kk], u[ijk+ii        ]) - interp2(u[ijk+ii-jj-kk], u[ijk+ii-jj    ]) )
                                          * dyi ) ) ) - ( evisch[ijk    ]
                                      * ( ( ( interp2(v[ijk    -kk], v[ijk        ]) - interp2(v[ijk-ii-kk], v[ijk-ii    ]) )
                                          * dxi )
                                        + ( ( interp2(u[ijk        -kk], u[ijk            ]) - interp2(u[ijk    -jj-kk], u[ijk    -jj    ]) )
                                          * dyi ) ) ) ) )
                                  * dxi );

                    // w * d/dy(visc * dv/dy + visc * dv/dy)
                    vw_diss[k] += ( ( interp2(w[ijk-jj], w[ijk    ])
                                    * ( ( ( 2 * interp2(evisc[ijk    -kk], evisc[ijk        ]) )
                                      * ( interp2(v[ijk+jj-kk], v[ijk+jj    ]) - interp2(v[ijk    -kk], v[ijk        ]) ) ) - ( ( 2 * interp2(evisc[ijk-jj-kk], evisc[ijk-jj    ]) )
                                      * ( interp2(v[ijk    -kk], v[ijk        ]) - interp2(v[ijk-jj-kk], v[ijk-jj    ]) ) ) ) )
                                  * dyi );

                    // w * d/dz(visc * du/dz + visc * dw/dx)
                    vw_diss[k] += ( ( interp2(w[ijk-jj], w[ijk    ])
                                    * ( ( interp2(evisc[ijk-jj    ], evisc[ijk        ])
                                      * ( ( ( interp2(v[ijk    ], v[ijk+kk]) - interp2(v[ijk-kk], v[ijk    ]) )
                                          * dzi[k  ] )
                                        + ( ( interp2(w[ijk        ], w[ijk    +kk]) - interp2(w[ijk-jj    ], w[ijk-jj+kk]) )
                                          * dyi ) ) ) - ( interp2(evisc[ijk-jj-kk], evisc[ijk    -kk])
                                      * ( ( ( interp2(v[ijk-kk], v[ijk    ]) - interp2(v[ijk-kk2], v[ijk-kk]) )
                                          * dzi[k-1] )
                                        + ( ( interp2(w[ijk    -kk], w[ijk        ]) - interp2(w[ijk-jj-kk], w[ijk-jj    ]) )
                                          * dyi ) ) ) ) )
                                  * dzhi[k] );

                    // v * d/dx(visc * dw/dx + visc * du/dz)
                    vw_diss[k] += ( ( interp2(v[ijk-kk], v[ijk    ])
                                    * ( ( evisch[ijk+ii]
                                      * ( ( ( interp2(w[ijk+ii-jj], w[ijk+ii    ]) - interp2(w[ijk    -jj], w[ijk        ]) )
                                          * dxi )
                                        + ( ( interp2(u[ijk+ii-jj    ], u[ijk+ii        ]) - interp2(u[ijk+ii-jj-kk], u[ijk+ii    -kk]) )
                                          * dzhi[k] ) ) ) - ( evisch[ijk    ]
                                      * ( ( ( interp2(w[ijk    -jj], w[ijk        ]) - interp2(w[ijk-ii-jj], w[ijk-ii    ]) )
                                          * dxi )
                                        + ( ( interp2(u[ijk    -jj    ], u[ijk            ]) - interp2(u[ijk    -jj-kk], u[ijk        -kk]) )
                                          * dzhi[k] ) ) ) ) )
                                  * dxi );

                    // v * d/dy(visc * dw/dy + visc * dv/dz)
                    vw_diss[k] += ( ( interp2(v[ijk-kk], v[ijk    ])
                                    * ( ( interp2(evisc[ijk    -kk], evisc[ijk        ])
                                      * ( ( ( interp2(w[ijk    ], w[ijk+jj]) - interp2(w[ijk-jj], w[ijk    ]) )
                                          * dyi )
                                        + ( ( interp2(v[ijk        ], v[ijk+jj    ]) - interp2(v[ijk    -kk], v[ijk+jj-kk]) )
                                          * dzhi[k] ) ) ) - ( interp2(evisc[ijk-jj-kk], evisc[ijk-jj    ])
                                      * ( ( ( interp2(w[ijk-jj], w[ijk    ]) - interp2(w[ijk-jj2], w[ijk-jj]) )
                                          * dyi )
                                        + ( ( interp2(v[ijk-jj    ], v[ijk        ]) - interp2(v[ijk-jj-kk], v[ijk    -kk]) )
                                          * dzhi[k] ) ) ) ) )
                                  * dyi );

                    // v * d/dz(visc * dw/dz + visc * dw/dz)
                    vw_diss[k] += ( ( interp2(v[ijk-kk], v[ijk    ])
                                    * ( ( ( 2 * interp2(evisc[ijk-jj    ], evisc[ijk        ]) )
                                      * ( ( interp2(w[ijk-jj+kk], w[ijk    +kk]) - interp2(w[ijk-jj    ], w[ijk        ]) )
                                        * dzi[k  ] ) ) - ( ( 2 * interp2(evisc[ijk-jj-kk], evisc[ijk    -kk]) )
                                      * ( ( interp2(w[ijk-jj    ], w[ijk        ]) - interp2(w[ijk-jj-kk], w[ijk    -kk]) )
                                        * dzi[k-1] ) ) ) )
                                  * dzhi[k] );
                }
            tke_diss[k] += 0.5 * (u2_diss[k] + v2_diss[k]);
        }

        int k = grid.kstart;
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                const int ij  = i + j*jj;

                const double evisc_utop   = interp2_4(evisc[ijk], evisc[ijk+kk], evisc[ijk-ii+kk], evisc[ijk-ii]);
                const double evisc_vtop   = interp2_4(evisc[ijk], evisc[ijk+kk], evisc[ijk-jj+kk], evisc[ijk-jj]);

                // 2 u * d/dz( visc * du/dz )
                u2_diss[k] += 2 * (u[ijk]-umean[k]) * ( evisc_utop * (u[ijk+kk] - u[ijk   ]) * dzhi[k+1] + ufluxbot[ij]) * dzi[k];

                // 2 v * d/dz( visc * dv/dz )
                v2_diss[k] += 2 * (v[ijk]-vmean[k]) * ( evisc_vtop * (v[ijk+kk] - v[ijk   ]) * dzhi[k+1] + vfluxbot[ij]) * dzi[k];

                // 2 * w * d/dz( visc * dw/dz )
                // What to do with evisc at surface (term visc * dw/dz at surface)?
                tke_diss[k] += wz[ijk] * ( interp2(evisc[ijk], evisc[ijk+kk]) * (wz[ijk+kk] - wz[ijk   ]) * dzhi[k+1] ) * 2 * dzi[k];

                // uw_diss is zero at surface for no-slip case, unequal for free-slip...
            }
        tke_diss[k] += 0.5 * (u2_diss[k] + v2_diss[k]);
    }

    // Calculate diffusion terms by splitting them into transport and dissipation.
    if(false)
    {
        // Diffusion term (d/dxj(nu * dui^2/dxj))
        // Lower boundary
        int k = grid.kstart;
        for (int j=grid.jstart; j<grid.jend; ++j)
        {
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                const double evisc_utop = interp2_4(evisc[ijk], evisc[ijk+kk], evisc[ijk-ii+kk], evisc[ijk-ii]);
                const double evisc_vtop = interp2_4(evisc[ijk], evisc[ijk+kk], evisc[ijk-jj+kk], evisc[ijk-jj]);

                // d/dz(visc * du^2/dz)
                u2_visc[k] +=  ( evisc_utop * (pow(u[ijk+kk]-umean[k+1], 2) - pow(u[ijk   ]-umean[k  ], 2)) * dzhi[k+1] ) * dzi[k];

                // d/dz(visc * duw/dx)
                u2_visc[k] += 2* ( evisc_utop * (interp2_4(u[ijk]-umean[k], u[ijk+ii]-umean[k], u[ijk+ii+kk]-umean[k+1], u[ijk+kk]-umean[k+1]) * w[ijk+kk   ] -
                                                 interp2_4(u[ijk]-umean[k], u[ijk-ii]-umean[k], u[ijk-ii+kk]-umean[k+1], u[ijk+kk]-umean[k+1]) * w[ijk-ii+kk] ) * dxi ) * dzi[k];

                // d/dz(visc * dv^2/dz)
                v2_visc[k] +=  ( evisc_vtop * (pow(v[ijk+kk]-vmean[k+1], 2) - pow(v[ijk   ]-vmean[k  ], 2)) * dzhi[k+1] ) * dzi[k];

                // d/dz(visc * dvw/dy)
                v2_visc[k] += 2* ( evisc_vtop * (interp2_4(v[ijk]-vmean[k], v[ijk+jj]-vmean[k], v[ijk+jj+kk]-vmean[k+1], v[ijk+kk]-vmean[k+1]) * w[ijk+kk   ] -
                                                 interp2_4(v[ijk]-vmean[k], v[ijk-jj]-vmean[k], v[ijk-jj+kk]-vmean[k+1], v[ijk+kk]-vmean[k+1]) * w[ijk-jj+kk] ) * dyi ) * dzi[k];

                // d/dz(visc * dw^2/dz)
                tke_visc[k] += 1.5 * ( interp2(evisc[ijk], evisc[ijk+kk]) * (pow(wz[ijk+kk], 2) - pow(wz[ijk   ], 2)) * dzhi[k+1] ) * dzi[k];
            }
            tke_visc[k] += 0.5 * (u2_visc[k] + v2_visc[k]);
        }

        // TODO what to do with terms that require interpolated evisc at surface?
        // Diffusion term (d/dxj(nu * dui^2/dxj))
        for (int k=grid.kstart+1; k<grid.kend-1; ++k)
        {
            for (int j=grid.jstart; j<grid.jend; ++j)
                #pragma ivdep
                for (int i=grid.istart; i<grid.iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    const double evisc_utop = interp2_4(evisc[ijk], evisc[ijk+kk], evisc[ijk-ii+kk], evisc[ijk-ii]);
                    const double evisc_ubot = interp2_4(evisc[ijk], evisc[ijk-kk], evisc[ijk-ii-kk], evisc[ijk-ii]);

                    // d/dz(visc * du^2/dz)
                    u2_visc[k] +=  ( evisc_utop * (pow(u[ijk+kk]-umean[k+1], 2) - pow(u[ijk   ]-umean[k  ], 2)) * dzhi[k+1] -
                                     evisc_ubot * (pow(u[ijk   ]-umean[k  ], 2) - pow(u[ijk-kk]-umean[k-1], 2)) * dzhi[k  ] ) * dzi[k];

                    // d/dz(visc * duw/dx)
                    u2_visc[k] += 2* ( evisc_utop * (interp2_4(u[ijk]-umean[k], u[ijk+ii]-umean[k], u[ijk+ii+kk]-umean[k+1], u[ijk+kk]-umean[k+1]) * w[ijk+kk   ] -
                                                     interp2_4(u[ijk]-umean[k], u[ijk-ii]-umean[k], u[ijk-ii+kk]-umean[k+1], u[ijk+kk]-umean[k+1]) * w[ijk-ii+kk] ) * dxi -
                                       evisc_ubot * (interp2_4(u[ijk]-umean[k], u[ijk+ii]-umean[k], u[ijk+ii-kk]-umean[k-1], u[ijk-kk]-umean[k-1]) * w[ijk      ] -
                                                     interp2_4(u[ijk]-umean[k], u[ijk-ii]-umean[k], u[ijk-ii-kk]-umean[k-1], u[ijk-kk]-umean[k-1]) * w[ijk-ii   ] ) * dxi ) * dzi[k];

                    const double evisc_vtop = interp2_4(evisc[ijk], evisc[ijk+kk], evisc[ijk-jj+kk], evisc[ijk-jj]);
                    const double evisc_vbot = interp2_4(evisc[ijk], evisc[ijk-kk], evisc[ijk-jj-kk], evisc[ijk-jj]);

                    // d/dz(visc * dv^2/dz)
                    v2_visc[k] +=  ( evisc_vtop * (pow(v[ijk+kk]-vmean[k+1], 2) - pow(v[ijk   ]-vmean[k  ], 2)) * dzhi[k+1] -
                                     evisc_vbot * (pow(v[ijk   ]-vmean[k  ], 2) - pow(v[ijk-kk]-vmean[k-1], 2)) * dzhi[k  ] ) * dzi[k];

                    // d/dz(visc * dvw/dy)
                    v2_visc[k] += 2* ( evisc_vtop * (interp2_4(v[ijk]-vmean[k], v[ijk+jj]-vmean[k], v[ijk+jj+kk]-vmean[k+1], v[ijk+kk]-vmean[k+1]) * w[ijk+kk   ] -
                                                     interp2_4(v[ijk]-vmean[k], v[ijk-jj]-vmean[k], v[ijk-jj+kk]-vmean[k+1], v[ijk+kk]-vmean[k+1]) * w[ijk-jj+kk] ) * dyi -
                                       evisc_vbot * (interp2_4(v[ijk]-vmean[k], v[ijk+jj]-vmean[k], v[ijk+jj-kk]-vmean[k-1], v[ijk-kk]-vmean[k-1]) * w[ijk      ] -
                                                     interp2_4(v[ijk]-vmean[k], v[ijk-jj]-vmean[k], v[ijk-jj-kk]-vmean[k-1], v[ijk-kk]-vmean[k-1]) * w[ijk-jj   ] ) * dyi ) * dzi[k];

                    // d/dz(visc * dw^2/dz)
                    w2_visc[k] += 3 * (evisc[ijk   ] * (pow(w[ijk+kk], 2) - pow(w[ijk   ], 2)) * dzi[k  ] -
                                       evisc[ijk-kk] * (pow(w[ijk   ], 2) - pow(w[ijk-kk], 2)) * dzi[k-1]) * dzhi[k];

                    // d/dz(visc * dw^2/dz)
                    tke_visc[k] += 1.5 * ( interp2(evisc[ijk], evisc[ijk+kk]) * (pow(wz[ijk+kk], 2) - pow(wz[ijk   ], 2)) * dzhi[k+1] -
                                           interp2(evisc[ijk], evisc[ijk-kk]) * (pow(wz[ijk   ], 2) - pow(wz[ijk-kk], 2)) * dzhi[k  ] ) * dzi[k];
                }
            tke_visc[k] += 0.5 * (u2_visc[k] + v2_visc[k]);
        }

        for (int k=grid.kstart; k<grid.kend; ++k)
        {
            for (int j=grid.jstart; j<grid.jend; ++j)
                #pragma ivdep
                for (int i=grid.istart; i<grid.iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // -2 * visc * ((du/dx)^2 + (du/dy)^2 + (du/dz)^2)
                    u2_diss[k] -= 2 * interp2(evisc[ijk], evisc[ijk-ii]) *
                                             ( 2* pow( (interp2(u[ijk]-umean[k], u[ijk+ii]-umean[k  ]) - interp2(u[ijk   ]-umean[k], u[ijk-ii   ]-umean[k  ])) * dxi,    2) +
                                               pow( (interp2(u[ijk]-umean[k], u[ijk+jj]-umean[k  ]) - interp2(u[ijk   ]-umean[k], u[ijk-jj   ]-umean[k  ])) * dyi,    2) +
                                               pow( (interp2(u[ijk]-umean[k], u[ijk+kk]-umean[k+1]) - interp2(u[ijk   ]-umean[k], u[ijk-kk   ]-umean[k-1])) * dzi[k], 2) );

                    // -2 * visc * ((dv/dx)^2 + (dv/dy)^2 + (dv/dz)^2)
                    v2_diss[k] -= 2 * interp2(evisc[ijk], evisc[ijk-jj]) *
                                             ( pow( (interp2(v[ijk]-vmean[k], v[ijk+ii]-vmean[k  ]) - interp2(v[ijk   ]-vmean[k], v[ijk-ii   ]-vmean[k  ])) * dxi,    2) +
                                               pow( (interp2(v[ijk]-vmean[k], v[ijk+jj]-vmean[k  ]) - interp2(v[ijk   ]-vmean[k], v[ijk-jj   ]-vmean[k  ])) * dyi,    2) +
                                               pow( (interp2(v[ijk]-vmean[k], v[ijk+kk]-vmean[k+1]) - interp2(v[ijk   ]-vmean[k], v[ijk-kk   ]-vmean[k-1])) * dzi[k], 2) );

                    // -2 * visc * (du/dx*du/dx + du/dy*dv/dx + du/dz*dw/dx)
                    // TODO first term (du/dx*du/dx) could be moved to previous u2_diss term
                    u2_diss[k] -= 2 * interp2(evisc[ijk], evisc[ijk-ii]) *
                                             ( pow( (interp2(u[ijk]-umean[k], u[ijk+ii]-umean[k  ]) - interp2(u[ijk   ]-umean[k], u[ijk-ii   ]-umean[k  ])) * dxi, 2) +
                                                    (interp2(u[ijk]-umean[k], u[ijk+jj]-umean[k  ]) - interp2(u[ijk   ]-umean[k], u[ijk-jj   ]-umean[k  ])) * dyi     *
                                                    (interp2(v[ijk]-vmean[k], v[ijk+jj]-vmean[k  ]) - interp2(v[ijk-ii]-vmean[k], v[ijk-jj+jj]-vmean[k  ])) * dxi     +
                                                    (interp2(u[ijk]-umean[k], u[ijk+kk]-umean[k+1]) - interp2(u[ijk   ]-umean[k], u[ijk-kk   ]-umean[k-1])) * dzi[k]  *
                                                    (wz[ijk]                                        - wz[ijk-ii]                                          ) * dxi );

                    // -2 * visc * (dv/dx*du/dy + dv/dy*dv/dy + dv/dz*dw/dy)
                    // TODO second term (dv/dy*dv/dy) could be moved to previous v2_diss term
                    v2_diss[k] -= 2 * interp2(evisc[ijk], evisc[ijk-jj]) *
                                                  ( (interp2(v[ijk]-vmean[k], v[ijk+ii]-vmean[k  ]) - interp2(v[ijk   ]-vmean[k], v[ijk-ii   ]-vmean[k  ])) * dxi     *
                                                    (interp2(u[ijk]-umean[k], u[ijk+ii]-umean[k  ]) - interp2(u[ijk-jj]-umean[k], u[ijk+ii-jj]-umean[k  ])) * dyi     +
                                               pow( (interp2(v[ijk]-vmean[k], v[ijk+jj]-vmean[k  ]) - interp2(v[ijk   ]-vmean[k], v[ijk-jj   ]-vmean[k  ])) * dyi, 2) +
                                                    (interp2(v[ijk]-vmean[k], v[ijk+kk]-vmean[k+1]) - interp2(v[ijk   ]-vmean[k], v[ijk-kk   ]-vmean[k-1])) * dzi[k]  *
                                                    (wz[ijk]                                        - wz[ijk-jj]                                          ) * dyi );

                    // -2 * visc * ((dw/dx)^2 + (dw/dy)^2 + (dw/dz)^2)
                    tke_diss[k] -= evisc[ijk] * ( pow( (interp2(wz[ijk+ii], wz[ijk]) - interp2(wz[ijk], wz[ijk-ii])) * dxi, 2) +
                                                  pow( (interp2(wz[ijk+jj], wz[ijk]) - interp2(wz[ijk], wz[ijk-jj])) * dyi, 2) +
                                                  pow( (w[ijk+kk] - w[ijk]) * dzi[k], 2) );

                    // -2 * visc * (dw/dx*du/dz + dw/dy*dv/dz + dw/dz*dw/dz)
                    tke_diss[k] -= evisc[ijk] * ( (interp2  (wz[ijk], wz[ijk+ii]) - interp2(wz[ijk], wz[ijk-ii])) * dxi *
                                                  (interp2_4(u[ijk]-umean[k], u[ijk+ii]-umean[k], u[ijk+ii+kk]-umean[k+1], u[ijk+kk]-umean[k+1]) -
                                                   interp2_4(u[ijk]-umean[k], u[ijk+ii]-umean[k], u[ijk+ii-kk]-umean[k-1], u[ijk-kk]-umean[k-1]) ) * dzi[k] +
                                                  (interp2  (wz[ijk], wz[ijk+jj]) - interp2(wz[ijk], wz[ijk-jj])) * dyi *
                                                  (interp2_4(v[ijk]-vmean[k], v[ijk+jj]-vmean[k], v[ijk+jj+kk]-vmean[k+1], v[ijk+kk]-vmean[k+1]) -
                                                   interp2_4(v[ijk]-vmean[k], v[ijk+jj]-vmean[k], v[ijk+jj-kk]-vmean[k-1], v[ijk-kk]-vmean[k-1]) ) * dzi[k] +
                                                   pow((w[ijk+kk] - w[ijk]) * dzi[k], 2) );
                }
            tke_diss[k] += 0.5 * (u2_diss[k] + v2_diss[k]);
        }

        // TODO BvS: what to do with w2_diss at surface? Most terms are zero (because e.g. dw/dx, dw/dy are zero), leaving only the vsc*dw/dz terms..
        for (int k=grid.kstart+1; k<grid.kend-1; ++k)
            for (int j=grid.jstart; j<grid.jend; ++j)
                #pragma ivdep
                for (int i=grid.istart; i<grid.iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // -2 * visc * ((dw/dx)^2 + (dw/dy)^2 + (dw/dz)^2)
                    w2_diss[k] -= 2 * interp2(evisc[ijk], evisc[ijk-kk]) *
                                             ( pow( (interp2(w[ijk], w[ijk+ii]) - interp2(w[ijk], w[ijk-ii])) * dxi,     2) +
                                               pow( (interp2(w[ijk], w[ijk+jj]) - interp2(w[ijk], w[ijk-jj])) * dyi,     2) +
                                               pow( (interp2(w[ijk], w[ijk+kk]) - interp2(w[ijk], w[ijk-kk])) * dzhi[k], 2) );

                    // -2 * visc * (dw/dx*du/dz + dw/dy+dv/dz + dw/dz*dw/dz
                    w2_diss[k] -= 2 * interp2(evisc[ijk], evisc[ijk-kk]) *
                                      ( (interp2(w[ijk], w[ijk+ii]) - interp2(w[ijk], w[ijk-ii])) * dxi *
                                        (interp2(u[ijk]-umean[k], u[ijk+ii]-umean[k]) - interp2(u[ijk-kk]-umean[k-1], u[ijk-kk+ii]-umean[k-1])) * dzhi[k] +
                                        (interp2(w[ijk], w[ijk+jj]) - interp2(w[ijk], w[ijk-jj])) * dyi *
                                        (interp2(v[ijk]-vmean[k], v[ijk+jj]-vmean[k]) - interp2(v[ijk-kk]-vmean[k-1], v[ijk-kk+jj]-vmean[k-1])) * dzhi[k] +
                                        pow((wz[ijk] - wz[ijk-kk]) * dzhi[k], 2) );
                }
    }

    // Calculate sum over all processes, and calc mean profiles
    master.sum(u2_visc,  grid.kcells);
    master.sum(v2_visc,  grid.kcells);
    master.sum(w2_visc,  grid.kcells);
    master.sum(tke_visc, grid.kcells);
    master.sum(uw_visc,  grid.kcells);
    master.sum(u2_diss,  grid.kcells);
    master.sum(v2_diss,  grid.kcells);
    master.sum(w2_diss,  grid.kcells);
    master.sum(tke_diss, grid.kcells);
    master.sum(uw_diss,  grid.kcells);
    master.sum(vw_diss,  grid.kcells);

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        u2_visc[k]  /= ijtot;
        v2_visc[k]  /= ijtot;
        tke_visc[k] /= ijtot;
        u2_diss[k]  /= ijtot;
        v2_diss[k]  /= ijtot;
        tke_diss[k] /= ijtot;
    }

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        w2_visc[k]  /= ijtot;
        uw_visc[k]  /= ijtot;
        w2_diss[k]  /= ijtot;
        uw_diss[k]  /= ijtot;
        vw_diss[k]  /= ijtot;
    }
}

/**
 * Calculate the budget terms arrising from diffusion, for a fixed viscosity
 * molecular diffusion (nu*d/dxj(dui^2/dxj)) and dissipation (-2*nu*(dui/dxj)^2)
 * @param TO-DO
 */
void Budget_2::calc_diffusion_terms_DNS(double* const restrict u2_visc, double* const restrict v2_visc,
                                        double* const restrict w2_visc, double* const restrict tke_visc, double* const restrict uw_visc,
                                        double* const restrict u2_diss, double* const restrict v2_diss,
                                        double* const restrict w2_diss, double* const restrict tke_diss, double* const restrict uw_diss,
                                        double* const restrict wz, double* const restrict wx, double* const restrict wy,
                                        const double* const restrict u, const double* const restrict v,
                                        const double* const restrict w, const double* const restrict umean, const double* const restrict vmean,
                                        const double* const restrict dzi, const double* const restrict dzhi,
                                        const double dxi, const double dyi, const double visc)
{
    // Interpolate the vertical velocity to {xh,y,zh} (wx, below u) and {x,yh,zh} (wy, below v)
    const int wloc [3] = {0,0,1};
    const int wxloc[3] = {1,0,1};
    const int wyloc[3] = {0,1,1};

    grid.interpolate_2nd(wx, w, wloc, wxloc);
    grid.interpolate_2nd(wy, w, wloc, wyloc);

    const int ii = 1;
    const int jj = grid.icells;
    const int kk = grid.ijcells;
    const int ijtot = grid.itot * grid.jtot;

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        u2_visc [k] = 0;
        v2_visc [k] = 0;
        tke_visc[k] = 0;
        u2_diss [k] = 0;
        v2_diss [k] = 0;
        tke_diss[k] = 0;
    }

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        w2_visc [k] = 0;
        uw_visc [k] = 0;
        w2_diss [k] = 0;
        uw_diss [k] = 0;
    }

    // Calculate w at full levels (grid center)
    for (int k=grid.kstart; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                wz[ijk] = interp2(w[ijk], w[ijk+kk]);
            }

    // Set ghost cells such that the velocity interpolated to the boundaries is zero
    int ks = grid.kstart;
    int ke = grid.kend-1;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijks = i + j*jj + ks*kk;
            const int ijke = i + j*jj + ke*kk;
            wz[ijks-kk] = -wz[ijks];
            wz[ijke+kk] = -wz[ijke];
        }

    // Molecular diffusion term (nu*d/dxj(dui^2/dxj))
    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // visc * d/dz(du^2/dz)
                u2_visc[k] += visc * ( (pow(u[ijk+kk]-umean[k+1], 2) - pow(u[ijk   ]-umean[k  ], 2)) * dzhi[k+1] -
                                       (pow(u[ijk   ]-umean[k  ], 2) - pow(u[ijk-kk]-umean[k-1], 2)) * dzhi[k  ] ) * dzi[k];

                // visc * d/dz(dv^2/dz)
                v2_visc[k] += visc * ( (pow(v[ijk+kk]-vmean[k+1], 2) - pow(v[ijk   ]-vmean[k  ], 2)) * dzhi[k+1] -
                                       (pow(v[ijk   ]-vmean[k  ], 2) - pow(v[ijk-kk]-vmean[k-1], 2)) * dzhi[k  ] ) * dzi[k];

                // visc * d/dz(dw^2/dz)
                tke_visc[k] += 0.5 * visc * ( (pow(wz[ijk+kk], 2) - pow(wz[ijk   ], 2)) * dzhi[k+1] -
                                              (pow(wz[ijk   ], 2) - pow(wz[ijk-kk], 2)) * dzhi[k  ] ) * dzi[k];
            }
        tke_visc[k] += 0.5 * (u2_visc[k] + v2_visc[k]);
    }

    // Lower boundary (z=0)
    int k = grid.kstart;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj + k*kk;

            // visc * d/dz(dw^2/dz)
            // w[kstart-1] = -w[kstart+1]
            w2_visc[k] += visc * ( (pow(w[ijk+kk], 2) - pow( w[ijk   ], 2)) * dzi[k  ] -
                                   (pow(w[ijk   ], 2) - pow(-w[ijk+kk], 2)) * dzi[k-1] ) * dzhi[k];

            // visc * d/dz(duw/dz)
            // wx[kstart-1] = -wx[kstart+1]
            // Calculate u at dz below surface, extrapolating gradient between u[kstart] and u[kstart-1]
            const double utmp = 1.5*(u[ijk-kk]-umean[k-1]) - 0.5*(u[ijk]-umean[k]);
            uw_visc[k] += visc * ( ( interp2(u[ijk   ]-umean[k  ], u[ijk+kk   ]-umean[k+1]) *  wx[ijk+kk] -
                                     interp2(u[ijk   ]-umean[k  ], u[ijk-kk   ]-umean[k-1]) *  wx[ijk   ] ) * dzi[k  ] -
                                   ( interp2(u[ijk   ]-umean[k  ], u[ijk-kk   ]-umean[k-1]) *  wx[ijk   ] -
                                     utmp                                                   * -wx[ijk+kk] ) * dzi[k-1] ) * dzhi[k];
        }

    // Top boundary (z=zsize)
    k = grid.kend;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj + k*kk;

            // visc * d/dz(dw^2/dz)
            // w[kend+1] = -w[kend-1]
            w2_visc[k] += visc * ( (pow(-w[ijk-kk], 2) - pow(w[ijk   ], 2)) * dzi[k  ] -
                                   (pow( w[ijk   ], 2) - pow(w[ijk-kk], 2)) * dzi[k-1] ) * dzhi[k];

            // visc * d/dz(duw/dz)
            // wx[kend+1] = -wx[kend-1]
            // Calculate u at dz above top, extrapolating gradient between u[kend] and u[kend-1]
            const double utmp = 1.5*(u[ijk]-umean[k]) - 0.5*(u[ijk-kk]-umean[k-1]);
            uw_visc[k] += visc * ( ( utmp                                                   * -wx[ijk-kk] -
                                     interp2(u[ijk   ]-umean[k  ], u[ijk-kk   ]-umean[k-1]) *  wx[ijk   ] ) * dzi[k  ] -
                                   ( interp2(u[ijk   ]-umean[k  ], u[ijk-kk   ]-umean[k-1]) *  wx[ijk   ] -
                                     interp2(u[ijk-kk]-umean[k-1], u[ijk-kk-kk]-umean[k-2]) *  wx[ijk-kk] ) * dzi[k-1] ) * dzhi[k];
        }

    // Interior
    for (int k=grid.kstart+1; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // visc * d/dz(dw^2/dz)
                w2_visc[k] += visc * ( (pow(w[ijk+kk], 2) - pow(w[ijk   ], 2)) * dzi[k  ] -
                                       (pow(w[ijk   ], 2) - pow(w[ijk-kk], 2)) * dzi[k-1] ) * dzhi[k];

                // visc * d/dz(duw/dz)
                uw_visc[k] += visc * ( ( interp2(u[ijk   ]-umean[k  ], u[ijk+kk   ]-umean[k+1]) * wx[ijk+kk] -
                                         interp2(u[ijk   ]-umean[k  ], u[ijk-kk   ]-umean[k-1]) * wx[ijk   ] ) * dzi[k  ] -
                                       ( interp2(u[ijk   ]-umean[k  ], u[ijk-kk   ]-umean[k-1]) * wx[ijk   ] -
                                         interp2(u[ijk-kk]-umean[k-1], u[ijk-kk-kk]-umean[k-2]) * wx[ijk-kk] ) * dzi[k-1] ) * dzhi[k];
            }

    // Dissipation term (-2*nu*(dui/dxj)^2)
    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // -2 * visc * ((du/dx)^2 + (du/dy)^2 + (du/dz)^2)
                u2_diss[k] -= 2 * visc * ( pow( (interp2(u[ijk]-umean[k], u[ijk+ii]-umean[k  ]) - interp2(u[ijk]-umean[k], u[ijk-ii]-umean[k  ])) * dxi,    2) +
                                           pow( (interp2(u[ijk]-umean[k], u[ijk+jj]-umean[k  ]) - interp2(u[ijk]-umean[k], u[ijk-jj]-umean[k  ])) * dyi,    2) +
                                           pow( (interp2(u[ijk]-umean[k], u[ijk+kk]-umean[k+1]) - interp2(u[ijk]-umean[k], u[ijk-kk]-umean[k-1])) * dzi[k], 2) );

                // -2 * visc * ((dv/dx)^2 + (dv/dy)^2 + (dv/dz)^2)
                v2_diss[k] -= 2 * visc * ( pow( (interp2(v[ijk]-vmean[k], v[ijk+ii]-vmean[k  ]) - interp2(v[ijk]-vmean[k], v[ijk-ii]-vmean[k  ])) * dxi,    2) +
                                           pow( (interp2(v[ijk]-vmean[k], v[ijk+jj]-vmean[k  ]) - interp2(v[ijk]-vmean[k], v[ijk-jj]-vmean[k  ])) * dyi,    2) +
                                           pow( (interp2(v[ijk]-vmean[k], v[ijk+kk]-vmean[k+1]) - interp2(v[ijk]-vmean[k], v[ijk-kk]-vmean[k-1])) * dzi[k], 2) );

                // -2 * visc * ((dw/dx)^2 + (dw/dy)^2 + (dw/dz)^2)
                tke_diss[k] -=    visc * ( pow( (w[ijk+ii] - w[ijk]) * dxi,    2) +
                                           pow( (w[ijk+jj] - w[ijk]) * dyi,    2) +
                                           pow( (w[ijk+kk] - w[ijk]) * dzi[k], 2) );

                // -2 * visc * du/dx * dw/dx
                uw_diss[k] -= 2 * visc * ( interp2_4(u[ijk]-umean[k], u[ijk+ii]-umean[k], u[ijk+ii-kk]-umean[k-1], u[ijk-kk]-umean[k-1]) -
                                           interp2_4(u[ijk]-umean[k], u[ijk-ii]-umean[k], u[ijk-ii-kk]-umean[k-1], u[ijk-kk]-umean[k-1]) ) * dxi *
                                         ( w[ijk] - w[ijk-ii] ) * dxi;

                // -2 * visc * du/dy * dw/dy
                uw_diss[k] -= 2 * visc * ( interp2_4(u[ijk]-umean[k], u[ijk+jj]-umean[k], u[ijk+jj-kk]-umean[k-1], u[ijk-kk]-umean[k-1]) -
                                           interp2_4(u[ijk]-umean[k], u[ijk-jj]-umean[k], u[ijk-jj-kk]-umean[k-1], u[ijk-kk]-umean[k-1]) ) * dyi *
                                         ( interp2_4(w[ijk]         , w[ijk+jj]         , w[ijk+jj-ii]           , w[ijk-ii]           ) -
                                           interp2_4(w[ijk]         , w[ijk-jj]         , w[ijk-jj-ii]           , w[ijk-ii]           ) ) * dyi;
            }
            tke_diss[k] += 0.5 * (u2_diss[k] + v2_diss[k]);
    }

    // Bottom boundary (z=0)
    k = grid.kstart;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj + k*kk;

            // -2 * visc * ((dw/dx)^2 + (dw/dy)^2 + (dw/dz)^2)
            // w @ full level kstart-1 = -w @ full level kstart+1
            w2_diss[k] -= 2 * visc * ( pow( (interp2(w[ijk], w[ijk+ii]) - interp2(w[ijk], w[ijk-ii])) * dxi,     2) +
                                       pow( (interp2(w[ijk], w[ijk+jj]) - interp2(w[ijk], w[ijk-jj])) * dyi,     2) +
                                       pow( (2*interp2(w[ijk], w[ijk+kk])                           ) * dzhi[k], 2) );

            // -2 * visc * du/dz * dw/dz
            // w @ full level kstart-1 = -w @ full level kstart+1
            uw_diss[k] -= 2 * visc * ((u[ijk]-umean[k]) - (u[ijk-kk]-umean[k-1])) * dzhi[k] *
                                     2*interp2_4(w[ijk], w[ijk+kk], w[ijk+kk-ii], w[ijk-ii]) * dzhi[k];
        }

    // Top boundary (z=zsize)
    k = grid.kend;
    for (int j=grid.jstart; j<grid.jend; ++j)
        #pragma ivdep
        for (int i=grid.istart; i<grid.iend; ++i)
        {
            const int ijk = i + j*jj + k*kk;

            // -2 * visc * ((dw/dx)^2 + (dw/dy)^2 + (dw/dz)^2)
            // w @ full level kend = -w @ full level kend-1
            w2_diss[k] -= 2 * visc * ( pow( (interp2(w[ijk], w[ijk+ii]) - interp2(w[ijk], w[ijk-ii])) * dxi,     2) +
                                       pow( (interp2(w[ijk], w[ijk+jj]) - interp2(w[ijk], w[ijk-jj])) * dyi,     2) +
                                       pow( (                        -2 * interp2(w[ijk], w[ijk-kk])) * dzhi[k], 2) );

            // -2 * visc * du/dz * dw/dz
            // w @ full level kend = - w @ full level kend-1
            uw_diss[k] -= 2 * visc * ((u[ijk]-umean[k]) - (u[ijk-kk]-umean[k-1])) * dzhi[k] *
                                     -2*interp2_4(w[ijk], w[ijk-kk], w[ijk-kk-ii], w[ijk-ii]) * dzhi[k];
        }

    // Interior
    for (int k=grid.kstart+1; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // -2 * visc * ((dw/dx)^2 + (dw/dy)^2 + (dw/dz)^2)
                w2_diss[k] -= 2 * visc * ( pow( (interp2(w[ijk], w[ijk+ii]) - interp2(w[ijk], w[ijk-ii])) * dxi,     2) +
                                           pow( (interp2(w[ijk], w[ijk+jj]) - interp2(w[ijk], w[ijk-jj])) * dyi,     2) +
                                           pow( (interp2(w[ijk], w[ijk+kk]) - interp2(w[ijk], w[ijk-kk])) * dzhi[k], 2) );

                // -2 * visc * du/dz * dw/dz
                uw_diss[k] -= 2 * visc * ((u[ijk]-umean[k]) - (u[ijk-kk]-umean[k-1])) * dzhi[k] *
                                         ( interp2_4(w[ijk], w[ijk+kk], w[ijk+kk-ii], w[ijk-ii]) -
                                           interp2_4(w[ijk], w[ijk-kk], w[ijk-kk-ii], w[ijk-ii]) ) * dzhi[k];
            }

    // Calculate sum over all processes, and calc mean profiles
    master.sum(u2_visc,  grid.kcells);
    master.sum(v2_visc,  grid.kcells);
    master.sum(w2_visc,  grid.kcells);
    master.sum(tke_visc, grid.kcells);
    master.sum(uw_visc,  grid.kcells);
    master.sum(u2_diss,  grid.kcells);
    master.sum(v2_diss,  grid.kcells);
    master.sum(w2_diss,  grid.kcells);
    master.sum(tke_diss, grid.kcells);
    master.sum(uw_diss,  grid.kcells);

    for (int k=grid.kstart; k<grid.kend; ++k)
    {
        u2_visc[k]  /= ijtot;
        v2_visc[k]  /= ijtot;
        tke_visc[k] /= ijtot;
        u2_diss[k]  /= ijtot;
        v2_diss[k]  /= ijtot;
        tke_diss[k] /= ijtot;
    }

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        w2_visc[k]  /= ijtot;
        uw_visc[k]  /= ijtot;
        w2_diss[k]  /= ijtot;
        uw_diss[k]  /= ijtot;
    }
}

/**
 * Calculate the budget terms arrising from buoyancy
 * @param TO-DO
 */
void Budget_2::calc_buoyancy_terms(double* const restrict w2_buoy, double* const restrict tke_buoy,
                                   double* const restrict uw_buoy, double* const restrict vw_buoy,
                                   const double* const restrict u, const double* const restrict v,
                                   const double* const restrict w, const double* const restrict b,
                                   const double* const restrict umean, const double* const restrict vmean,
                                   const double* const restrict bmean)
{
    const int ii = 1;
    const int jj = grid.icells;
    const int kk = grid.ijcells;
    const int ijtot = grid.itot * grid.jtot;

    for (int k=grid.kstart; k<grid.kend; ++k)
        tke_buoy[k] = 0;

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        w2_buoy[k] = 0;
        uw_buoy[k] = 0;
        vw_buoy[k] = 0;
    }

    for (int k=grid.kstart; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // w'b'
                tke_buoy[k] += interp2(w[ijk], w[ijk+kk]) * (b[ijk] - bmean[k]);
            }

    for (int k=grid.kstart+1; k<grid.kend; ++k)
        for (int j=grid.jstart; j<grid.jend; ++j)
            #pragma ivdep
            for (int i=grid.istart; i<grid.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // w'b'
                w2_buoy[k] += 2 * interp2(b[ijk]-bmean[k], b[ijk-kk]-bmean[k-1]) * w[ijk];

                // u'b'
                uw_buoy[k] += interp2  (u[ijk]-umean[k], u[ijk-kk]-umean[k-1]) *
                              interp2_4(b[ijk]-bmean[k], b[ijk-ii]-bmean[k], b[ijk-ii-kk]-bmean[k-1], b[ijk-kk]-bmean[k-1]);

                // v'b'
                vw_buoy[k] += interp2  (v[ijk]-vmean[k], v[ijk-kk]-vmean[k-1]) *
                              interp2_4(b[ijk]-bmean[k], b[ijk-jj]-bmean[k], b[ijk-jj-kk]-bmean[k-1], b[ijk-kk]-bmean[k-1]);

            }

    // Calculate sum over all processes, and calc mean profiles
    master.sum(w2_buoy,  grid.kcells);
    master.sum(tke_buoy, grid.kcells);
    master.sum(uw_buoy,  grid.kcells);
    master.sum(vw_buoy,  grid.kcells);

    for (int k=grid.kstart; k<grid.kend; ++k)
        tke_buoy[k] /= ijtot;

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        w2_buoy[k] /= ijtot;
        uw_buoy[k] /= ijtot;
        vw_buoy[k] /= ijtot;
    }
}
