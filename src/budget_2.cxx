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
#include "stats.h"
#include <netcdfcpp.h>

#include "budget.h"
#include "budget_2.h"

using namespace Finite_difference::O2;

Budget_2::Budget_2(Input* inputin, Master* masterin, Grid* gridin, Fields* fieldsin, Thermo* thermoin, Diff* diffin, Stats* statsin) :
    Budget(inputin, masterin, gridin, fieldsin, thermoin, diffin, statsin)
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

    // Calculate the shear production and turbulent transport terms
    calc_advection_terms(m->profs["u2_shear"].data, m->profs["v2_shear"].data, m->profs["tke_shear"].data, m->profs["uw_shear"].data,
                         m->profs["u2_turb"].data,  m->profs["v2_turb"].data,  m->profs["w2_turb"].data, m->profs["tke_turb"].data, m->profs["uw_turb"].data,
                         fields.u->data, fields.v->data, fields.w->data, umodel, vmodel,
                         fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, grid.dzi, grid.dzhi);

    // Calculate the pressure transport and redistribution terms
    calc_pressure_terms(m->profs["w2_pres"].data,  m->profs["tke_pres"].data, m->profs["uw_pres"].data,
                        m->profs["u2_rdstr"].data, m->profs["v2_rdstr"].data, m->profs["w2_rdstr"].data, m->profs["uw_rdstr"].data,
                        fields.u->data, fields.v->data, fields.w->data, fields.sd["p"]->data, umodel, vmodel,
                        grid.dzi, grid.dzhi, grid.dxi, grid.dyi);

    // Calculate the diffusive transport and dissipation terms
    if(diff.get_name() == "2" || diff.get_name() == "4")
        calc_diffusion_terms_DNS(m->profs["u2_visc"].data, m->profs["v2_visc"].data, m->profs["w2_visc"].data, m->profs["tke_visc"].data, m->profs["uw_visc"].data,
                                 m->profs["u2_diss"].data, m->profs["v2_diss"].data, m->profs["w2_diss"].data, m->profs["tke_diss"].data, m->profs["uw_diss"].data,
                                 fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, fields.atmp["tmp3"]->data, fields.u->data, fields.v->data, fields.w->data, umodel, vmodel,
                                 grid.dzi, grid.dzhi, grid.dxi, grid.dyi, fields.visc);
    else if(diff.get_name() == "smag2")
        std::cout << "Ha! not implemented yet" << std::endl;

    if(thermo.get_switch() != "0")
    {
        // Store the buoyancy in the tmp1 field
        thermo.get_thermo_field(fields.atmp["tmp1"], fields.atmp["tmp2"], "b");

        // Calculate mean fields
        grid.calc_mean(fields.atmp["tmp1"]->datamean, fields.atmp["tmp1"]->data, grid.kcells);
        grid.calc_mean(fields.sd["p"]->datamean, fields.sd["p"]->data, grid.kcells);

        // Calculate buoyancy terms
        calc_buoyancy_terms(m->profs["w2_buoy"].data, m->profs["tke_buoy"].data, m->profs["uw_buoy"].data,
                            fields.u->data, fields.w->data, fields.atmp["tmp1"]->data,
                            umodel, fields.atmp["tmp1"]->datamean);
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
                                    double* const restrict tke_shear, double* const restrict uw_shear,
                                    double* const restrict u2_turb,  double* const restrict v2_turb,
                                    double* const restrict w2_turb, double* const restrict tke_turb, double* const restrict uw_turb,
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
        w2_turb  [k] = 0;
        uw_turb  [k] = 0;
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
                            (u[ijk-kk]-umean[k-1]) * pow(interp2(wx[ijk], wx[ijk+kk]), 2) ) * dzhi[k];
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
            }

    // Calculate sum over all processes, and calc mean profiles
    master.sum(u2_shear,  grid.kcells);
    master.sum(v2_shear,  grid.kcells);
    master.sum(tke_shear, grid.kcells);
    master.sum(uw_shear,  grid.kcells);
    master.sum(u2_turb,   grid.kcells);
    master.sum(v2_turb,   grid.kcells);
    master.sum(w2_turb,   grid.kcells);
    master.sum(tke_turb,  grid.kcells);
    master.sum(uw_turb,   grid.kcells);

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
        w2_turb  [k] /= ijtot;
        uw_turb  [k] /= ijtot;
    }
}

/**
 * Calculate the budget terms arrising from pressure:
 * pressure transport (-2*dpu_i/dxi) and redistribution (2p*dui/dxi)
 * @param TO-DO
 */
void Budget_2::calc_pressure_terms(double* const restrict w2_pres,  double* const restrict tke_pres, double* const restrict uw_pres,
                                   double* const restrict u2_rdstr, double* const restrict v2_rdstr, double* const restrict w2_rdstr, double* const restrict uw_rdstr,
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
        w2_rdstr[k] = 0;
        uw_rdstr[k] = 0;
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

                uw_pres[k]  -= ( interp2(p[ijk   ], p[ijk-kk   ]) * w[ijk   ] -
                                 interp2(p[ijk-ii], p[ijk-ii-kk]) * w[ijk-ii] ) * dxi +
                               ( interp2(p[ijk   ], p[ijk-ii   ]) * u[ijk   ] -
                                 interp2(p[ijk-kk], p[ijk-ii-kk]) * u[ijk-kk] ) * dzhi[k];
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
    master.sum(u2_rdstr, grid.kcells);
    master.sum(v2_rdstr, grid.kcells);
    master.sum(w2_rdstr, grid.kcells);
    master.sum(uw_rdstr, grid.kcells);

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
        w2_rdstr[k] /= ijtot;
        uw_rdstr[k] /= ijtot;
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
void Budget_2::calc_buoyancy_terms(double* const restrict w2_buoy, double* const restrict tke_buoy, double* const restrict uw_buoy,
                                   const double* const restrict u, const double* const restrict w, const double* const restrict b,
                                   const double* const restrict umean, const double* const restrict bmean)
{
    const int ii = 1;
    const int jj = grid.icells;
    const int kk = grid.ijcells;
    const int ijtot = grid.itot * grid.jtot;

    for (int k=grid.kstart; k<grid.kend; ++k)
        tke_buoy[k] = 0;

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        w2_buoy [k] = 0;
        uw_buoy [k] = 0;
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
                w2_buoy[k] += 2 * interp2(b[ijk], b[ijk-kk]) * w[ijk];

                // u'b'
                uw_buoy[k] += interp2  (u[ijk]-umean[k], u[ijk-kk]-umean[k]) *
                              interp2_4(b[ijk]-bmean[k], b[ijk-ii]-bmean[k], b[ijk-ii-kk]-bmean[k-1], b[ijk-kk]-bmean[k-1]);
            }

    // Calculate sum over all processes, and calc mean profiles
    master.sum(w2_buoy,  grid.kcells);
    master.sum(tke_buoy, grid.kcells);
    master.sum(uw_buoy,  grid.kcells);

    for (int k=grid.kstart; k<grid.kend; ++k)
        tke_buoy[k] /= ijtot;

    for (int k=grid.kstart; k<grid.kend+1; ++k)
    {
        w2_buoy[k] /= ijtot;
        uw_buoy[k] /= ijtot;
    }
}

