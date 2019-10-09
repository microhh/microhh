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

#include <iostream>
#include <cmath>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "finite_difference.h"
#include "fast_math.h"
#include "model.h"
#include "thermo.h"
#include "diff.h"
#include "advec.h"
#include "force.h"
#include "stats.h"

#include "budget.h"
#include "budget_2.h"

namespace
{
    using namespace Finite_difference::O2;

    /**
     * Calculate the kinetic and turbulence kinetic energy
     */
    template<typename TF>
    void calc_kinetic_energy(
            TF* const restrict ke, TF* const restrict tke,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict umodel, const TF* const restrict vmodel,
            const TF utrans, const TF vtrans,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        using Fast_math::pow2;

        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        for (int k=kstart; k<kend; ++k)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    const TF u2 = pow2(interp2(u[ijk]+utrans, u[ijk+ii]+utrans));
                    const TF v2 = pow2(interp2(v[ijk]+vtrans, v[ijk+jj]+vtrans));
                    const TF w2 = pow2(interp2(w[ijk]       , w[ijk+kk]       ));

                    ke[ijk] = 0.5 * (u2 + v2 + w2);
                }

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    const TF u2 = pow2(interp2(u[ijk]-umodel[k], u[ijk+ii]-umodel[k]));
                    const TF v2 = pow2(interp2(v[ijk]-vmodel[k], v[ijk+jj]-vmodel[k]));
                    const TF w2 = pow2(interp2(w[ijk]          , w[ijk+kk]          ));

                    tke[ijk] = TF(0.5) * (u2 + v2 + w2);
                }
        }
    }

    /**
     * Calculate the budget terms related to shear production (-2 u_i*u_j * d<u_i>/dx_j)
     */
    template<typename TF>
    void calc_shear_terms(
            TF* const restrict u2_shear, TF* const restrict v2_shear, TF* const restrict tke_shear,
            TF* const restrict uw_shear, TF* const restrict vw_shear,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict umean, const TF* const restrict vmean,
            const TF* const restrict wx, const TF* const restrict wy,
            const TF* const restrict dzi, const TF* const restrict dzhi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        // Calculate shear terms (-2u_iw d<u_i>/dz)
        for (int k=kstart; k<kend; ++k)
        {
            const TF dudz = (interp2(umean[k], umean[k+1]) - interp2(umean[k-1], umean[k]) ) * dzi[k];
            const TF dvdz = (interp2(vmean[k], vmean[k+1]) - interp2(vmean[k-1], vmean[k]) ) * dzi[k];

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    u2_shear[ijk] = -2 * (u[ijk]-umean[k]) * interp2(wx[ijk], wx[ijk+kk]) * dudz;
                    v2_shear[ijk] = -2 * (v[ijk]-vmean[k]) * interp2(wy[ijk], wy[ijk+kk]) * dvdz;

                    uw_shear[ijk] = -pow(wx[ijk], 2) * (umean[k] - umean[k-1]) * dzhi[k];
                    vw_shear[ijk] = -pow(wy[ijk], 2) * (vmean[k] - vmean[k-1]) * dzhi[k];

                    tke_shear[ijk] = 0.5*(u2_shear[ijk] + v2_shear[ijk]);
                }
        }
    }

    /**
     * Calculate the budget terms related to turbulent transport (-d(u_i^2*u_j)/dx_j)
     */
    template<typename TF>
    void calc_turb_terms(
            TF* const restrict u2_turb,  TF* const restrict v2_turb,
            TF* const restrict w2_turb, TF* const restrict tke_turb,
            TF* const restrict uw_turb, TF* const restrict vw_turb,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict umean, const TF* const restrict vmean,
            const TF* const restrict wx, const TF* const restrict wy,
            const TF* const restrict dzi, const TF* const restrict dzhi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        // Calculate turbulent transport terms (-d(u_i^2*w)/dz)
        for (int k=kstart; k<kend; ++k)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    u2_turb[ijk] = - ( pow(interp2(u[ijk]-umean[k], u[ijk+kk]-umean[k+1]), 2) * wx[ijk+kk]
                                     - pow(interp2(u[ijk]-umean[k], u[ijk-kk]-umean[k-1]), 2) * wx[ijk   ] ) * dzi[k];

                    v2_turb[ijk] = - ( pow(interp2(v[ijk]-vmean[k], v[ijk+kk]-vmean[k+1]), 2) * wy[ijk+kk]
                                     - pow(interp2(v[ijk]-vmean[k], v[ijk-kk]-vmean[k-1]), 2) * wy[ijk   ] ) * dzi[k];

                    tke_turb[ijk] = - TF(0.5) * ( pow(w[ijk+kk], 3) - pow(w[ijk], 3) ) * dzi[k]
                                    + TF(0.5) * (u2_turb[ijk] + v2_turb[ijk]);
                }
        }

        // Lower boundary kstart (z=0)
        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // w^3 @ full level below sfc == -w^3 @ full level above sfc
                w2_turb[ijk] = - TF(2.) * pow(interp2(w[ijk], w[ijk+kk]), 3) * dzhi[k];

                // w^2 @ full level below sfc == w^2 @ full level above sfc
                uw_turb[ijk] = - ( (u[ijk]   -umean[k  ]) * pow(interp2(wx[ijk], wx[ijk+kk]), 2)
                                 - (u[ijk-kk]-umean[k-1]) * pow(interp2(wx[ijk], wx[ijk-kk]), 2) ) * dzhi[k];

                vw_turb[ijk] = - ( (v[ijk]   -vmean[k  ]) * pow(interp2(wy[ijk], wy[ijk+kk]), 2)
                                 - (v[ijk-kk]-vmean[k-1]) * pow(interp2(wy[ijk], wy[ijk-kk]), 2) ) * dzhi[k];
            }

        // Top boundary kstart (z=zsize)
        k = kend;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // w^3 @ full level above top == -w^3 @ full level below top
                w2_turb[ijk] = - TF(2.) * pow(interp2(w[ijk], w[ijk-kk]), 3) * dzhi[k];

                // w^2 @ full level above top == w^2 @ full level below top
                uw_turb[ijk] = - ( (u[ijk]   -umean[k  ]) * pow(interp2(wx[ijk], wx[ijk-kk]), 2)
                                 - (u[ijk-kk]-umean[k-1]) * pow(interp2(wx[ijk], wx[ijk-kk]), 2) ) * dzhi[k];

                // w^2 @ full level above top == w^2 @ full level below top
                vw_turb[ijk] =  - ( (v[ijk]   -vmean[k  ]) * pow(interp2(wy[ijk], wy[ijk-kk]), 2)
                                  - (v[ijk-kk]-vmean[k-1]) * pow(interp2(wy[ijk], wy[ijk-kk]), 2) ) * dzhi[k];
            }

        // Inner domain
        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    w2_turb[ijk] = - ( pow(interp2(w[ijk], w[ijk+kk]), 3)
                                     - pow(interp2(w[ijk], w[ijk-kk]), 3) ) * dzhi[k];

                    uw_turb[ijk] = - ( (u[ijk]   -umean[k  ]) * pow(interp2(wx[ijk], wx[ijk+kk]), 2)
                                     - (u[ijk-kk]-umean[k-1]) * pow(interp2(wx[ijk], wx[ijk-kk]), 2) ) * dzhi[k];

                    vw_turb[ijk] = - ( (v[ijk]   -vmean[k  ]) * pow(interp2(wy[ijk], wy[ijk+kk]), 2)
                                     - (v[ijk-kk]-vmean[k-1]) * pow(interp2(wy[ijk], wy[ijk-kk]), 2) ) * dzhi[k];
                }
    }

    /**
     * Calculate the budget terms related to coriolis force.
     */
    template<typename TF>
    void calc_coriolis_terms(
            TF* const restrict u2_cor, TF* const restrict v2_cor,
            TF* const restrict uw_cor, TF* const restrict vw_cor,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict umean, const TF* const restrict vmean, const TF fc,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    u2_cor[ijk] = TF( 2.) * (u[ijk]-umean[k]) * (interp22(v[ijk-ii], v[ijk], v[ijk-ii+jj], v[ijk+jj])-vmean[k]) * fc;
                    v2_cor[ijk] = TF(-2.) * (v[ijk]-vmean[k]) * (interp22(u[ijk-jj], u[ijk], u[ijk+ii-jj], u[ijk+ii])-umean[k]) * fc;
                }

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    uw_cor[ijk] = interp2(w[ijk], w[ijk-ii]) *
                        interp2(interp22(v[ijk   ]-vmean[k], v[ijk-ii   ]-vmean[k], v[ijk-ii-kk   ]-vmean[k-1], v[ijk-kk   ]-vmean[k-1]),
                                interp22(v[ijk+jj]-vmean[k], v[ijk-ii+jj]-vmean[k], v[ijk-ii+jj-kk]-vmean[k-1], v[ijk+jj-kk]-vmean[k-1])) * fc;

                    vw_cor[ijk] = interp2(w[ijk], w[ijk-jj]) *
                        interp2(interp22(u[ijk   ]-umean[k], u[ijk-jj   ]-umean[k], u[ijk-jj-kk   ]-umean[k-1], u[ijk-kk   ]-umean[k-1]),
                                interp22(u[ijk+ii]-umean[k], u[ijk+ii-jj]-umean[k], u[ijk+ii-jj-kk]-umean[k-1], u[ijk+ii-kk]-umean[k-1])) * fc;
                }
    }

    /**
     * Calculate the budget terms arrising from pressure:
     * pressure transport (-2*dpu_i/dxi) and redistribution (2p*dui/dxi)
     */
    template<typename TF>
    void calc_pressure_terms(
            TF* const restrict w2_pres,  TF* const restrict tke_pres,
            TF* const restrict uw_pres,  TF* const restrict vw_pres,
            TF* const restrict u2_rdstr, TF* const restrict v2_rdstr, TF* const restrict w2_rdstr,
            TF* const restrict uw_rdstr, TF* const restrict vw_rdstr,
            const TF* const restrict u, const TF* const restrict v,
            const TF* const restrict w, const TF* const restrict p,
            const TF* const restrict umean, const TF* const restrict vmean,
            const TF* const restrict dzi, const TF* const restrict dzhi, const TF dxi, const TF dyi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        // Pressure transport term (-2*dpu_i/dxi)
        for (int k=kstart; k<kend; ++k)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    tke_pres[ijk] = - ( interp2(p[ijk], p[ijk+kk]) * w[ijk+kk] -
                                        interp2(p[ijk], p[ijk-kk]) * w[ijk   ] ) * dzi[k];

                    uw_pres[ijk] = - ( interp2(p[ijk   ], p[ijk-kk   ]) * w[ijk    ] -
                                       interp2(p[ijk-ii], p[ijk-ii-kk]) * w[ijk-ii ] ) * dxi +
                                     ( interp2(p[ijk   ], p[ijk-ii   ]) * (u[ijk   ]-umean[k  ]) -
                                       interp2(p[ijk-kk], p[ijk-ii-kk]) * (u[ijk-kk]-umean[k-1]) ) * dzhi[k];

                    vw_pres[ijk] = - ( interp2(p[ijk-kk   ], p[ijk   ]) * w[ijk    ]  -
                                       interp2(p[ijk-jj-kk], p[ijk-jj]) * w[ijk-jj ] ) * dyi +
                                     ( interp2(p[ijk-jj   ], p[ijk   ]) * (v[ijk   ]-vmean[k  ]) -
                                       interp2(p[ijk-jj-kk], p[ijk-kk]) * (v[ijk-kk]-vmean[k-1]) ) * dzhi[k];
                }
        }

        // Lower boundary (z=0)
        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // w @ full level below sfc == -w @ full level above sfc
                w2_pres[ijk] = TF(-2.) * ( interp2(w[ijk], w[ijk+kk]) * p[ijk   ] -
                                         - interp2(w[ijk], w[ijk+kk]) * p[ijk-kk] ) * dzhi[k];
            }

        // Top boundary (z=zsize)
        // TODO: what to do with w2_pres and uw_pres at the top boundary? Pressure at k=kend is undefined?

        // Inner domain
        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    w2_pres[ijk] = TF(-2.) * ( interp2(w[ijk], w[ijk+kk]) * p[ijk   ] -
                                               interp2(w[ijk], w[ijk-kk]) * p[ijk-kk] ) * dzhi[k];
                }

        // Pressure redistribution term (2p*dui/dxi)
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    u2_rdstr[ijk] = TF(2.) * interp2(p[ijk], p[ijk-ii]) *
                        ( interp2(u[ijk]-umean[k], u[ijk+ii]-umean[k]) -
                          interp2(u[ijk]-umean[k], u[ijk-ii]-umean[k]) ) * dxi;

                    v2_rdstr[ijk] = TF(2.) * interp2(p[ijk], p[ijk-jj]) *
                        ( interp2(v[ijk]-vmean[k], v[ijk+jj]-vmean[k]) -
                          interp2(v[ijk]-vmean[k], v[ijk-jj]-vmean[k]) ) * dyi;

                    uw_rdstr[ijk] = interp2_4(p[ijk], p[ijk-kk], p[ijk-ii-kk], p[ijk-ii]) *
                         ( ((u[ijk]-umean[k]) - (u[ijk-kk]-umean[k-1])) * dzhi[k] + (w[ijk] - w[ijk-ii]) * dxi );

                    vw_rdstr[ijk] = interp2_4(p[ijk], p[ijk-kk], p[ijk-jj-kk], p[ijk-jj]) *
                         ( ((v[ijk]-vmean[k]) - (v[ijk-kk]-vmean[k-1])) * dzhi[k] + (w[ijk] - w[ijk-jj]) * dyi );
                }

        // Lower boundary (z=0)
        k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // with w[kstart] == 0, dw/dz at surface equals (w[kstart+1] - w[kstart]) / dzi
                w2_rdstr[ijk] = TF(2.) * interp2(p[ijk], p[ijk-kk]) * (w[ijk+kk] - w[ijk]) * dzi[k];
            }


        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    w2_rdstr[ijk] = TF(2.) * interp2(p[ijk], p[ijk-kk]) *
                        ( interp2(w[ijk], w[ijk+kk]) - interp2(w[ijk], w[ijk-kk]) ) * dzhi[k];
                }
    }
}

template<typename TF>
Budget_2<TF>::Budget_2(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin,
        Thermo<TF>& thermoin, Diff<TF>& diffin, Advec<TF>& advecin, Force<TF>& forcein, Input& inputin) :
    Budget<TF>(masterin, gridin, fieldsin, thermoin, diffin, advecin, forcein, inputin),
    field3d_operators(masterin, gridin, fieldsin)
{
    // The LES flux budget requires one additional ghost cell in the horizontal.
    if (diff.get_switch() == Diffusion_type::Diff_smag2)
    {
        const int igc = 2;
        const int jgc = 2;
        const int kgc = 1;

        grid.set_minimum_ghost_cells(igc, jgc, kgc);
    }
}

template<typename TF>
Budget_2<TF>::~Budget_2()
{
}

template<typename TF>
void Budget_2<TF>::init()
{
    auto& gd = grid.get_grid_data();

    umodel.resize(gd.kcells);
    vmodel.resize(gd.kcells);
}

template<typename TF>
void Budget_2<TF>::create(Stats<TF>& stats)
{
    const std::string group_name = "budget";

    // Add the profiles for the kinetic energy to the statistics.
    stats.add_prof("ke" , "Kinetic energy" , "m2 s-2", "z", group_name);
    stats.add_prof("tke", "Turbulent kinetic energy" , "m2 s-2", "z", group_name);

    // Add the profiles for the kinetic energy budget to the statistics.
    stats.add_prof("u2_shear" , "Shear production term in U2 budget" , "m2 s-3", "z" , group_name);
    stats.add_prof("v2_shear" , "Shear production term in V2 budget" , "m2 s-3", "z" , group_name);
    stats.add_prof("tke_shear", "Shear production term in TKE budget", "m2 s-3", "z" , group_name);
    stats.add_prof("uw_shear" , "Shear production term in UW budget" , "m2 s-3", "zh", group_name);
    stats.add_prof("vw_shear" , "Shear production term in VW budget" , "m2 s-3", "zh", group_name);

    stats.add_prof("u2_turb" , "Turbulent transport term in U2 budget" , "m2 s-3", "z" , group_name);
    stats.add_prof("v2_turb" , "Turbulent transport term in V2 budget" , "m2 s-3", "z" , group_name);
    stats.add_prof("w2_turb" , "Turbulent transport term in W2 budget" , "m2 s-3", "zh", group_name);
    stats.add_prof("tke_turb", "Turbulent transport term in TKE budget", "m2 s-3", "z" , group_name);
    stats.add_prof("uw_turb" , "Turbulent transport term in UW budget" , "m2 s-3", "zh", group_name);
    stats.add_prof("vw_turb" , "Turbulent transport term in VW budget" , "m2 s-3", "zh", group_name);

    if (force.get_switch_lspres() == Large_scale_pressure_type::Geo_wind)
    {
        stats.add_prof("u2_cor", "Coriolis term in U2 budget", "m2 s-3", "z" , group_name);
        stats.add_prof("v2_cor", "Coriolis term in V2 budget", "m2 s-3", "z" , group_name);
        stats.add_prof("uw_cor", "Coriolis term in UW budget", "m2 s-3", "zh", group_name);
        stats.add_prof("vw_cor", "Coriolis term in VW budget", "m2 s-3", "zh", group_name);
    }

    /*
    if (diff.get_switch() != Diffusion_type::Disabled)
    {
        stats.add_prof("u2_diss" , "Dissipation term in U2 budget" , "m2 s-3", "z" , group_name);
        stats.add_prof("v2_diss" , "Dissipation term in V2 budget" , "m2 s-3", "z" , group_name);
        stats.add_prof("w2_diss" , "Dissipation term in W2 budget" , "m2 s-3", "zh", group_name);
        stats.add_prof("tke_diss", "Dissipation term in TKE budget", "m2 s-3", "z" , group_name);
        stats.add_prof("uw_diss" , "Dissipation term in UW budget" , "m2 s-3", "zh", group_name);
        stats.add_prof("vw_diss" , "Dissipation term in VW budget" , "m2 s-3", "zh", group_name);

        stats.add_prof("u2_visc" , "Viscous transport term in U2 budget" , "m2 s-3", "z" , group_name);
        stats.add_prof("v2_visc" , "Viscous transport term in V2 budget" , "m2 s-3", "z" , group_name);
        stats.add_prof("w2_visc" , "Viscous transport term in W2 budget" , "m2 s-3", "zh", group_name);
        stats.add_prof("tke_visc", "Viscous transport term in TKE budget", "m2 s-3", "z" , group_name);
        stats.add_prof("uw_visc" , "Viscous transport term in UW budget" , "m2 s-3", "zh", group_name);
        stats.add_prof("vw_visc" , "Viscous transport term in VW budget" , "m2 s-3", "zh", group_name);

        // For LES, add the total diffusive budget terms, which (unlike diss + visc) close
        if (diff.get_switch() == Diffusion_type::Diff_smag2)
        {
            stats.add_prof("u2_diff" , "Total diffusive term in U2 budget" , "m2 s-3", "z" , group_name);
            stats.add_prof("v2_diff" , "Total diffusive term in V2 budget" , "m2 s-3", "z" , group_name);
            stats.add_prof("w2_diff" , "Total diffusive term in W2 budget" , "m2 s-3", "zh", group_name);
            stats.add_prof("tke_diff", "Total diffusive term in TKE budget", "m2 s-3", "z" , group_name);
            stats.add_prof("uw_diff" , "Total diffusive term in UW budget" , "m2 s-3", "zh", group_name);
            stats.add_prof("vw_diff" , "Total diffusive term in VW budget" , "m2 s-3", "zh", group_name);
        }
    }


    if (thermo.get_switch() != "0")
    {
        stats.add_prof("w2_buoy" , "Buoyancy production/destruction term in W2 budget" , "m2 s-3", "zh", group_name);
        stats.add_prof("tke_buoy", "Buoyancy production/destruction term in TKE budget", "m2 s-3", "z" , group_name);
        stats.add_prof("uw_buoy" , "Buoyancy production/destruction term in UW budget" , "m2 s-3", "zh", group_name);
        stats.add_prof("vw_buoy" , "Buoyancy production/destruction term in VW budget" , "m2 s-3", "zh", group_name);

        stats.add_prof("b2_shear", "Shear production term in B2 budget"   , "m2 s-5", "z", group_name);
        stats.add_prof("b2_turb" , "Turbulent transport term in B2 budget", "m2 s-5", "z", group_name);

        stats.add_prof("bw_shear", "Shear production term in B2 budget"   , "m2 s-4", "zh", group_name);
        stats.add_prof("bw_turb" , "Turbulent transport term in B2 budget", "m2 s-4", "zh", group_name);

        if (diff.get_switch() != Diffusion_type::Disabled)
        {
            stats.add_prof("b2_visc" , "Viscous transport term in B2 budget", "m2 s-5", "z" , group_name);
            stats.add_prof("b2_diss" , "Dissipation term in B2 budget"      , "m2 s-5", "z" , group_name);
            stats.add_prof("bw_visc" , "Viscous transport term in BW budget", "m2 s-4", "zh", group_name);
            stats.add_prof("bw_diss" , "Dissipation term in BW budget"      , "m2 s-4", "zh", group_name);
        }

        stats.add_prof("bw_rdstr", "Redistribution term in BW budget"     , "m2 s-4", "zh", group_name);
        stats.add_prof("bw_buoy" , "Buoyancy term in BW budget"           , "m2 s-4", "zh", group_name);
        stats.add_prof("bw_pres" , "Pressure transport term in BW budget" , "m2 s-4", "zh", group_name);
    }

    stats.add_prof("w2_pres" , "Pressure transport term in W2 budget" , "m2 s-3", "zh", group_name);
    stats.add_prof("tke_pres", "Pressure transport term in TKE budget", "m2 s-3", "z" , group_name);
    stats.add_prof("uw_pres" , "Pressure transport term in UW budget" , "m2 s-3", "zh", group_name);
    stats.add_prof("vw_pres" , "Pressure transport term in VW budget" , "m2 s-3", "zh", group_name);

    stats.add_prof("u2_rdstr", "Pressure redistribution term in U2 budget", "m2 s-3", "z" , group_name);
    stats.add_prof("v2_rdstr", "Pressure redistribution term in V2 budget", "m2 s-3", "z" , group_name);
    stats.add_prof("w2_rdstr", "Pressure redistribution term in W2 budget", "m2 s-3", "zh", group_name);
    stats.add_prof("uw_rdstr", "Pressure redistribution term in UW budget", "m2 s-3", "zh", group_name);
    stats.add_prof("vw_rdstr", "Pressure redistribution term in VW budget", "m2 s-3", "zh", group_name);
    */
}

template<typename TF>
void Budget_2<TF>::exec_stats(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    // Calculate the mean of the fields.
    field3d_operators.calc_mean_profile(umodel.data(), fields.mp.at("u")->fld.data());
    field3d_operators.calc_mean_profile(vmodel.data(), fields.mp.at("v")->fld.data());

    // Calculate kinetic and turbulent kinetic energy
    auto ke  = fields.get_tmp();
    auto tke = fields.get_tmp();

    const TF no_offset = 0.;
    const TF no_threshold = 0.;

    calc_kinetic_energy(
            ke->fld.data(), tke->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            umodel.data(), vmodel.data(),
            grid.utrans, grid.vtrans,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    stats.calc_stats("ke" , *ke , no_offset, no_threshold);
    stats.calc_stats("tke", *tke, no_offset, no_threshold);

    auto wx = std::move(ke );
    auto wy = std::move(tke);

    // Interpolate w to the locations of u and v.
    const int wloc [3] = {0,0,1};
    const int wxloc[3] = {1,0,1};
    const int wyloc[3] = {0,1,1};

    grid.interpolate_2nd(wx->fld.data(), fields.mp.at("w")->fld.data(), wloc, wxloc);
    grid.interpolate_2nd(wy->fld.data(), fields.mp.at("w")->fld.data(), wloc, wyloc);

    auto u2_shear = fields.get_tmp();
    auto v2_shear = fields.get_tmp();
    auto tke_shear = fields.get_tmp();
    auto uw_shear = fields.get_tmp();
    auto vw_shear = fields.get_tmp();

    calc_shear_terms(
            u2_shear->fld.data(), v2_shear->fld.data(), tke_shear->fld.data(),
            uw_shear->fld.data(), vw_shear->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            umodel.data(), vmodel.data(),
            wx->fld.data(), wy->fld.data(),
            gd.dzi.data(), gd.dzhi.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    stats.calc_stats("u2_shear" , *u2_shear , no_offset, no_threshold);
    stats.calc_stats("v2_shear" , *v2_shear , no_offset, no_threshold);
    stats.calc_stats("tke_shear", *tke_shear, no_offset, no_threshold);
    stats.calc_stats("uw_shear" , *uw_shear , no_offset, no_threshold);
    stats.calc_stats("vw_shear" , *vw_shear , no_offset, no_threshold);

    auto u2_turb = std::move(u2_shear);
    auto v2_turb = std::move(v2_shear);
    auto w2_turb = fields.get_tmp();
    auto tke_turb = std::move(tke_shear);
    auto uw_turb = std::move(uw_shear);
    auto vw_turb = std::move(vw_shear);

    calc_turb_terms(
            u2_turb->fld.data(), v2_turb->fld.data(),
            w2_turb->fld.data(), tke_turb->fld.data(),
            uw_turb->fld.data(), vw_turb->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            umodel.data(), vmodel.data(),
            wx->fld.data(), wy->fld.data(),
            gd.dzi.data(), gd.dzhi.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    stats.calc_stats("u2_turb" , *u2_turb , no_offset, no_threshold);
    stats.calc_stats("v2_turb" , *v2_turb , no_offset, no_threshold);
    stats.calc_stats("w2_turb" , *w2_turb , no_offset, no_threshold);
    stats.calc_stats("tke_turb", *tke_turb, no_offset, no_threshold);
    stats.calc_stats("uw_turb" , *uw_turb , no_offset, no_threshold);
    stats.calc_stats("vw_turb" , *vw_turb , no_offset, no_threshold);

    fields.release_tmp(wx);
    fields.release_tmp(wy);

    fields.release_tmp(u2_turb);
    fields.release_tmp(v2_turb);
    fields.release_tmp(w2_turb);
    fields.release_tmp(tke_turb);
    fields.release_tmp(uw_turb);
    fields.release_tmp(vw_turb);

    if (force.get_switch_lspres() == Large_scale_pressure_type::Geo_wind)
    {
        auto u2_cor = fields.get_tmp();
        auto v2_cor = fields.get_tmp();
        auto uw_cor = fields.get_tmp();
        auto vw_cor = fields.get_tmp();

        const TF fc = force.get_coriolis_parameter();
        calc_coriolis_terms(
                u2_cor->fld.data(), v2_cor->fld.data(),
                uw_cor->fld.data(), vw_cor->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                umodel.data(), vmodel.data(), fc,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_stats("u2_cor", *u2_cor, no_offset, no_threshold);
        stats.calc_stats("v2_cor", *v2_cor, no_offset, no_threshold);
        stats.calc_stats("uw_cor", *uw_cor, no_offset, no_threshold);
        stats.calc_stats("vw_cor", *vw_cor, no_offset, no_threshold);

        fields.release_tmp(u2_cor);
        fields.release_tmp(v2_cor);
        fields.release_tmp(uw_cor);
        fields.release_tmp(vw_cor);
    }

    /*
    if(diff.get_switch() != "0")
    {
        // Calculate the diffusive transport and dissipation terms
        if(diff.get_switch() == "2" || diff.get_switch() == "4")
            calc_diffusion_terms_DNS(m->profs["u2_visc"].data, m->profs["v2_visc"].data, m->profs["w2_visc"].data, m->profs["tke_visc"].data, m->profs["uw_visc"].data,
                                     m->profs["u2_diss"].data, m->profs["v2_diss"].data, m->profs["w2_diss"].data, m->profs["tke_diss"].data, m->profs["uw_diss"].data,
                                     fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, fields.atmp["tmp3"]->data, fields.u->data, fields.v->data, fields.w->data, umodel, vmodel,
                                     grid.dzi, grid.dzhi, grid.dxi, grid.dyi, fields.visc);
        else if(diff.get_switch() == "smag2")
            calc_diffusion_terms_LES(m->profs["u2_diss"].data,  m->profs["v2_diss"].data, m->profs["w2_diss"].data,
                                     m->profs["tke_diss"].data, m->profs["uw_diss"].data, m->profs["vw_diss"].data,
                                     m->profs["u2_visc"].data,  m->profs["v2_visc"].data, m->profs["w2_visc"].data,
                                     m->profs["tke_visc"].data, m->profs["uw_visc"].data, m->profs["vw_visc"].data,
                                     m->profs["u2_diff"].data,  m->profs["v2_diff"].data, m->profs["w2_diff"].data,
                                     m->profs["tke_diff"].data, m->profs["uw_diff"].data, m->profs["vw_diff"].data,
                                     fields.atmp["tmp1"]->data, fields.atmp["tmp2"]->data, fields.atmp["tmp3"]->data,
                                     fields.u->data, fields.v->data, fields.w->data,
                                     fields.u->datafluxbot, fields.v->datafluxbot,
                                     fields.sd.at("evisc")->data, umodel, vmodel,
                                     grid.dzi, grid.dzhi, grid.dxi, grid.dyi);
    }

    if(thermo.get_switch() != "0")
    {
        // Get the buoyancy diffusivity from the thermo class
        const double diff_b = thermo.get_buoyancy_diffusivity();

        // Store the buoyancy in the tmp1 field
        thermo.get_thermo_field(fields.atmp["tmp1"], fields.atmp["tmp2"], "b", true);

        // Calculate mean fields
        grid.calc_mean(fields.atmp["tmp1"]->datamean, fields.atmp["tmp1"]->data, grid.kcells);
        grid.calc_mean(fields.sd.at("p")->datamean, fields.sd.at("p")->data, grid.kcells);

        // Calculate buoyancy terms
        calc_buoyancy_terms(m->profs["w2_buoy"].data, m->profs["tke_buoy"].data,
                            m->profs["uw_buoy"].data, m->profs["vw_buoy"].data,
                            fields.u->data, fields.v->data, fields.w->data, fields.atmp["tmp1"]->data,
                            umodel, vmodel, fields.atmp["tmp1"]->datamean);

        // Buoyancy variance and flux budgets
        calc_buoyancy_terms_scalar(m->profs["bw_buoy"].data,
                                   fields.atmp["tmp1"]->data, fields.atmp["tmp1"]->data,
                                   fields.atmp["tmp1"]->datamean, fields.atmp["tmp1"]->datamean);

        if (advec.get_switch() != "0")
            calc_advection_terms_scalar(m->profs["b2_shear"].data, m->profs["b2_turb"].data,
                                        m->profs["bw_shear"].data, m->profs["bw_turb"].data,
                                        fields.atmp["tmp1"]->data, fields.w->data, fields.atmp["tmp1"]->datamean,
                                        grid.dzi, grid.dzhi);

        if (diff.get_switch() == "2" || diff.get_switch() == "4")
            calc_diffusion_terms_scalar_DNS(m->profs["b2_visc"].data, m->profs["b2_diss"].data,
                                            m->profs["bw_visc"].data, m->profs["bw_diss"].data,
                                            fields.atmp["tmp1"]->data, fields.w->data,
                                            fields.atmp["tmp1"]->datamean,
                                            grid.dzi, grid.dzhi, grid.dxi, grid.dyi, fields.visc, diff_b);

        calc_pressure_terms_scalar(m->profs["bw_pres"].data,  m->profs["bw_rdstr"].data,
                                   fields.atmp["tmp1"]->data, fields.sd.at("p")->data,
                                   fields.atmp["tmp1"]->datamean, fields.sd.at("p")->datamean,
                                   grid.dzi, grid.dzhi);
    }

    // Calculate the pressure transport and redistribution terms
    calc_pressure_terms(m->profs["w2_pres"].data,     m->profs["tke_pres"].data,
                        m->profs["uw_pres"].data,     m->profs["vw_pres"].data,
                        m->profs["u2_rdstr"].data,    m->profs["v2_rdstr"].data, m->profs["w2_rdstr"].data,
                        m->profs["uw_rdstr"].data,    m->profs["vw_rdstr"].data,
                        fields.u->data, fields.v->data, fields.w->data, fields.sd.at("p")->data, umodel, vmodel,
                        grid.dzi, grid.dzhi, grid.dxi, grid.dyi);
    */
}

template class Budget_2<double>;
template class Budget_2<float>;
