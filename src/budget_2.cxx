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
#include "netcdf_interface.h"

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
            const TF* const restrict umodel, const TF* const restrict vmodel, const TF* const restrict wmodel,
            const TF utrans, const TF vtrans,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        using Fast_math::pow2;

        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        #pragma omp parallel for
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

                    ke[ijk] = TF(0.5) * (u2 + v2 + w2);
                }

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    const TF u2 = pow2(interp2(u[ijk]-umodel[k], u[ijk+ii]-umodel[k]));
                    const TF v2 = pow2(interp2(v[ijk]-vmodel[k], v[ijk+jj]-vmodel[k]));
                    const TF w2 = pow2(interp2(w[ijk]-wmodel[k], w[ijk+kk]-wmodel[k+1]));

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
            const TF* const restrict umean, const TF* const restrict vmean, const TF* const restrict wmean,
            const TF* const restrict wx, const TF* const restrict wy,
            const TF* const restrict dzi, const TF* const restrict dzhi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        // Calculate shear terms (-2u_iw d<u_i>/dz)
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
        {
            const TF dudz = (interp2(umean[k], umean[k+1]) - interp2(umean[k-1], umean[k]) ) * dzi[k];
            const TF dvdz = (interp2(vmean[k], vmean[k+1]) - interp2(vmean[k-1], vmean[k]) ) * dzi[k];

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    u2_shear[ijk] = TF(-2.) * (u[ijk]-umean[k]) * interp2(wx[ijk]-wmean[k], wx[ijk+kk]-wmean[k+1]) * dudz;
                    v2_shear[ijk] = TF(-2.) * (v[ijk]-vmean[k]) * interp2(wy[ijk]-wmean[k], wy[ijk+kk]-wmean[k+1]) * dvdz;

                    uw_shear[ijk] = -pow(wx[ijk], 2) * (umean[k] - umean[k-1]) * dzhi[k];
                    vw_shear[ijk] = -pow(wy[ijk], 2) * (vmean[k] - vmean[k-1]) * dzhi[k];

                    tke_shear[ijk] = TF(0.5)*(u2_shear[ijk] + v2_shear[ijk]);
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
            const TF* const restrict umean, const TF* const restrict vmean, const TF* const restrict wmean,
            const TF* const restrict wx, const TF* const restrict wy,
            const TF* const restrict dzi, const TF* const restrict dzhi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        // Calculate turbulent transport terms (-d(u_i^2*w)/dz)
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    u2_turb[ijk] = - ( pow(interp2(u[ijk]-umean[k], u[ijk+kk]-umean[k+1]), 2) * (wx[ijk+kk]-wmean[k+1])
                                     - pow(interp2(u[ijk]-umean[k], u[ijk-kk]-umean[k-1]), 2) * (wx[ijk   ]-wmean[k  ]) ) * dzi[k];

                    v2_turb[ijk] = - ( pow(interp2(v[ijk]-vmean[k], v[ijk+kk]-vmean[k+1]), 2) * (wy[ijk+kk]-wmean[k+1])
                                     - pow(interp2(v[ijk]-vmean[k], v[ijk-kk]-vmean[k-1]), 2) * (wy[ijk   ]-wmean[k  ]) ) * dzi[k];

                    tke_turb[ijk] = - TF(0.5) * ( pow(w[ijk+kk]-wmean[k+1], 3) - pow(w[ijk]-wmean[k], 3) ) * dzi[k]
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
                uw_turb[ijk] = - ( (u[ijk]   -umean[k  ]) * pow(interp2(wx[ijk]-wmean[k], wx[ijk+kk]-wmean[k+1]), 2)
                                 - (u[ijk-kk]-umean[k-1]) * pow(interp2(wx[ijk]-wmean[k], wx[ijk-kk]-wmean[k+1]), 2) ) * dzhi[k];

                vw_turb[ijk] = - ( (v[ijk]   -vmean[k  ]) * pow(interp2(wy[ijk]-wmean[k], wy[ijk+kk]-wmean[k+1]), 2)
                                 - (v[ijk-kk]-vmean[k-1]) * pow(interp2(wy[ijk]-wmean[k], wy[ijk-kk]-wmean[k+1]), 2) ) * dzhi[k];
            }

        // Top boundary kstart (z=zsize)
        k = kend;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // w^3 @ full level above top == -w^3 @ full level below top
                w2_turb[ijk] = - TF(2.) * pow(interp2(w[ijk]-wmean[k], w[ijk-kk]-wmean[k-1]), 3) * dzhi[k];

                // w^2 @ full level above top == w^2 @ full level below top
                uw_turb[ijk] = - ( (u[ijk]   -umean[k  ]) * pow(interp2(wx[ijk]-wmean[k], wx[ijk-kk]-wmean[k-1]), 2)
                                 - (u[ijk-kk]-umean[k-1]) * pow(interp2(wx[ijk]-wmean[k], wx[ijk-kk]-wmean[k-1]), 2) ) * dzhi[k];

                // w^2 @ full level above top == w^2 @ full level below top
                vw_turb[ijk] =  - ( (v[ijk]   -vmean[k  ]) * pow(interp2(wy[ijk]-wmean[k], wy[ijk-kk]-wmean[k-1]), 2)
                                  - (v[ijk-kk]-vmean[k-1]) * pow(interp2(wy[ijk]-wmean[k], wy[ijk-kk]-wmean[k-1]), 2) ) * dzhi[k];
            }

        // Inner domain
        #pragma omp parallel for
        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    w2_turb[ijk] = - ( pow(interp2(w[ijk]-wmean[k], w[ijk+kk]-wmean[k+1]), 3)
                                     - pow(interp2(w[ijk]-wmean[k], w[ijk-kk]), 3) ) * dzhi[k];

                    uw_turb[ijk] = - ( (u[ijk]   -umean[k  ]) * pow(interp2(wx[ijk]-wmean[k], wx[ijk+kk]-wmean[k+1]), 2)
                                     - (u[ijk-kk]-umean[k-1]) * pow(interp2(wx[ijk]-wmean[k], wx[ijk-kk]-wmean[k-1]), 2) ) * dzhi[k];

                    vw_turb[ijk] = - ( (v[ijk]   -vmean[k  ]) * pow(interp2(wy[ijk]-wmean[k], wy[ijk+kk]-wmean[k+1]), 2)
                                     - (v[ijk-kk]-vmean[k-1]) * pow(interp2(wy[ijk]-wmean[k], wy[ijk-kk]-wmean[k-1]), 2) ) * dzhi[k];
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
            const TF* const restrict umean, const TF* const restrict vmean, const TF* const restrict wmean, const TF fc,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    u2_cor[ijk] = TF( 2.) * (u[ijk]-umean[k]) * (interp22(v[ijk-ii], v[ijk], v[ijk-ii+jj], v[ijk+jj])-vmean[k]) * fc;
                    v2_cor[ijk] = TF(-2.) * (v[ijk]-vmean[k]) * (interp22(u[ijk-jj], u[ijk], u[ijk+ii-jj], u[ijk+ii])-umean[k]) * fc;
                }

        #pragma omp parallel for
        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    uw_cor[ijk] = interp2(w[ijk]-wmean[k], w[ijk-ii]-wmean[k]) *
                        interp2(interp22(v[ijk   ]-vmean[k], v[ijk-ii   ]-vmean[k], v[ijk-ii-kk   ]-vmean[k-1], v[ijk-kk   ]-vmean[k-1]),
                                interp22(v[ijk+jj]-vmean[k], v[ijk-ii+jj]-vmean[k], v[ijk-ii+jj-kk]-vmean[k-1], v[ijk+jj-kk]-vmean[k-1])) * fc;

                    vw_cor[ijk] = interp2(w[ijk]-wmean[k], w[ijk-jj]-wmean[k]) *
                        interp2(interp22(u[ijk   ]-umean[k], u[ijk-jj   ]-umean[k], u[ijk-jj-kk   ]-umean[k-1], u[ijk-kk   ]-umean[k-1]),
                                interp22(u[ijk+ii]-umean[k], u[ijk+ii-jj]-umean[k], u[ijk+ii-jj-kk]-umean[k-1], u[ijk+ii-kk]-umean[k-1])) * fc;
                }
    }

    /**
     * Calculate the budget terms related to pressure transport (-2*dpu_i/dxi)
     */
    template<typename TF>
    void calc_pressure_transport_terms(
            TF* const restrict w2_pres,  TF* const restrict tke_pres,
            TF* const restrict uw_pres,  TF* const restrict vw_pres,
            const TF* const restrict u, const TF* const restrict v,
            const TF* const restrict w, const TF* const restrict p,
            const TF* const restrict umean, const TF* const restrict vmean, const TF* const restrict wmean,
            const TF* const restrict dzi, const TF* const restrict dzhi, const TF dxi, const TF dyi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        // Pressure transport term (-2*dpu_i/dxi)
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    tke_pres[ijk] = - ( interp2(p[ijk], p[ijk+kk]) * (w[ijk+kk] -wmean[k+1])-
                                        interp2(p[ijk], p[ijk-kk]) * (w[ijk   ] -wmean[k  ])) * dzi[k];

                    uw_pres[ijk] = - ( interp2(p[ijk   ], p[ijk-kk   ]) * (w[ijk   ]-wmean[k  ])-
                                       interp2(p[ijk-ii], p[ijk-ii-kk]) * (w[ijk-ii]-wmean[k  ]) ) * dxi +
                                     ( interp2(p[ijk   ], p[ijk-ii   ]) * (u[ijk   ]-umean[k  ]) -
                                       interp2(p[ijk-kk], p[ijk-ii-kk]) * (u[ijk-kk]-umean[k-1]) ) * dzhi[k];

                    vw_pres[ijk] = - ( interp2(p[ijk-kk   ], p[ijk   ]) * (w[ijk   ]-wmean[k  ])  -
                                       interp2(p[ijk-jj-kk], p[ijk-jj]) * (w[ijk-jj]-wmean[k  ]) ) * dyi +
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
                w2_pres[ijk] = TF(-2.) * ( interp2(w[ijk]-wmean[k  ], w[ijk+kk]-wmean[k+1]) * p[ijk   ] -
                                         - interp2(w[ijk]-wmean[k  ], w[ijk+kk]-wmean[k+1]) * p[ijk-kk] ) * dzhi[k];
            }

        // Top boundary (z=zsize)
        // TODO: what to do with w2_pres and uw_pres at the top boundary? Pressure at k=kend is undefined?

        // Inner domain
        #pragma omp parallel for
        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    w2_pres[ijk] = TF(-2.) * ( interp2(w[ijk]-wmean[k  ], w[ijk+kk]-wmean[k+1]) * p[ijk   ] -
                                               interp2(w[ijk]-wmean[k  ], w[ijk-kk]-wmean[k-1]) * p[ijk-kk] ) * dzhi[k];
                }
    }

    /**
     * Calculate the budget terms related to redistribution (2p*dui/dxi)
     */
    template<typename TF>
    void calc_pressure_redistribution_terms(
            TF* const restrict u2_rdstr, TF* const restrict v2_rdstr, TF* const restrict w2_rdstr,
            TF* const restrict uw_rdstr, TF* const restrict vw_rdstr,
            const TF* const restrict u, const TF* const restrict v,
            const TF* const restrict w, const TF* const restrict p,
            const TF* const restrict umean, const TF* const restrict vmean, const TF* const restrict wmean,
            const TF* const restrict dzi, const TF* const restrict dzhi, const TF dxi, const TF dyi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        // Pressure redistribution term (2p*dui/dxi)
        #pragma omp parallel for
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

                    uw_rdstr[ijk] = interp22(p[ijk], p[ijk-kk], p[ijk-ii-kk], p[ijk-ii]) *
                         ( ((u[ijk]-umean[k]) - (u[ijk-kk]-umean[k-1])) * dzhi[k] + (w[ijk] - w[ijk-ii]) * dxi );

                    vw_rdstr[ijk] = interp22(p[ijk], p[ijk-kk], p[ijk-jj-kk], p[ijk-jj]) *
                         ( ((v[ijk]-vmean[k]) - (v[ijk-kk]-vmean[k-1])) * dzhi[k] + (w[ijk] - w[ijk-jj]) * dyi );
                }

        // Lower boundary (z=0)
        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // with w[kstart] == 0, dw/dz at surface equals (w[kstart+1] - w[kstart]) / dzi
                w2_rdstr[ijk] = TF(2.) * interp2(p[ijk], p[ijk-kk]) * (w[ijk+kk]-wmean[k+1] - (w[ijk]-wmean[k])) * dzi[k];
            }

        #pragma omp parallel for
        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    w2_rdstr[ijk] = TF(2.) * interp2(p[ijk], p[ijk-kk]) *
                        ( interp2(w[ijk]-wmean[k], w[ijk+kk]-wmean[k+1]) - interp2(w[ijk]-wmean[k], w[ijk-kk]-wmean[k-1]) ) * dzhi[k];
                }
    }

    /**
     * Calculate the budget terms arrising from diffusion, for a fixed viscosity
     * molecular diffusion (nu*d/dxj(dui^2/dxj)) and dissipation (-2*nu*(dui/dxj)^2)
     */
    template<typename TF>
    void calc_diffusion_transport_terms_dns(
            TF* const restrict u2_visc, TF* const restrict v2_visc,
            TF* const restrict w2_visc, TF* const restrict tke_visc, TF* const restrict uw_visc,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            TF* const restrict wz, TF* const restrict wx, TF* const restrict wy,
            const TF* const restrict umean, const TF* const restrict vmean, const TF* const restrict wmean,
            const TF* const restrict dzi, const TF* const restrict dzhi,
            const TF dxi, const TF dyi, const TF visc,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        // Calculate w at full levels (grid center)
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    wz[ijk] = interp2(w[ijk]-wmean[k], w[ijk+kk]-wmean[k+1]);
                }

        // Set ghost cells such that the velocity interpolated to the boundaries is zero
        int ks = kstart;
        int ke = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijks = i + j*jj + ks*kk;
                const int ijke = i + j*jj + ke*kk;
                wz[ijks-kk] = -wz[ijks];
                wz[ijke+kk] = -wz[ijke];
            }

        // Molecular diffusion term (nu*d/dxj(dui^2/dxj))
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // visc * d/dz(du^2/dz)
                    u2_visc[ijk] = visc * ( (pow(u[ijk+kk]-umean[k+1], 2) - pow(u[ijk   ]-umean[k  ], 2)) * dzhi[k+1] -
                                            (pow(u[ijk   ]-umean[k  ], 2) - pow(u[ijk-kk]-umean[k-1], 2)) * dzhi[k  ] ) * dzi[k];

                    // visc * d/dz(dv^2/dz)
                    v2_visc[ijk] = visc * ( (pow(v[ijk+kk]-vmean[k+1], 2) - pow(v[ijk   ]-vmean[k  ], 2)) * dzhi[k+1] -
                                            (pow(v[ijk   ]-vmean[k  ], 2) - pow(v[ijk-kk]-vmean[k-1], 2)) * dzhi[k  ] ) * dzi[k];

                    // visc * d/dz(dw^2/dz)
                    tke_visc[ijk] = TF(0.5) * visc * ( (pow(wz[ijk+kk], 2) - pow(wz[ijk   ], 2)) * dzhi[k+1] -
                                                       (pow(wz[ijk   ], 2) - pow(wz[ijk-kk], 2)) * dzhi[k  ] ) * dzi[k]
                                  + TF(0.5) * (u2_visc[ijk] + v2_visc[ijk]);
                }
        }

        // Lower boundary (z=0)
        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // visc * d/dz(dw^2/dz)
                // w[kstart-1] = -w[kstart+1]
                w2_visc[ijk] = visc * ( (pow(w[ijk+kk]-wmean[k+1], 2) - pow( w[ijk   ]-wmean[k  ], 2)) * dzi[k  ] -
                                        (pow(w[ijk   ]-wmean[k  ], 2) - pow( w[ijk+kk]-wmean[k+1], 2)) * dzi[k-1] ) * dzhi[k];

                // visc * d/dz(duw/dz)
                // wx[kstart-1] = -wx[kstart+1]
                // Calculate u at dz below surface, extrapolating gradient between u[kstart] and u[kstart-1]
                const TF utmp = TF(1.5)*(u[ijk-kk]-umean[k-1]) - TF(0.5)*(u[ijk]-umean[k]);
                uw_visc[ijk] = visc * ( ( interp2(u[ijk   ]-umean[k  ], u[ijk+kk   ]-umean[k+1]) *  (wx[ijk+kk] - wmean[k+1])-
                                          interp2(u[ijk   ]-umean[k  ], u[ijk-kk   ]-umean[k-1]) *  (wx[ijk   ] - wmean[k  ])) * dzi[k  ] -
                                        ( interp2(u[ijk   ]-umean[k  ], u[ijk-kk   ]-umean[k-1]) *  (wx[ijk   ] - wmean[k  ])-
                                        utmp *                                                     -(wx[ijk+kk] - wmean[k+1]) ) * dzi[k-1] ) * dzhi[k];
            }

        // Top boundary (z=zsize)
        k = kend;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // visc * d/dz(dw^2/dz)
                // w[kend+1] = -w[kend-1]
                w2_visc[ijk] = visc * ( (pow( w[ijk-kk]- wmean[k-1], 2) - pow(w[ijk   ]- wmean[k  ], 2)) * dzi[k  ] -
                                        (pow( w[ijk   ]- wmean[k  ], 2) - pow(w[ijk-kk]- wmean[k-1], 2)) * dzi[k-1] ) * dzhi[k];

                // visc * d/dz(duw/dz)
                // wx[kend+1] = -wx[kend-1]
                // Calculate u at dz above top, extrapolating gradient between u[kend] and u[kend-1]
                const TF utmp = TF(1.5)*(u[ijk]-umean[k]) - TF(0.5)*(u[ijk-kk]-umean[k-1]);
                uw_visc[ijk] = visc * ( ( utmp                                                   * -(wx[ijk-kk]- wmean[k-1]) -
                                          interp2(u[ijk   ]-umean[k  ], u[ijk-kk   ]-umean[k-1]) *  (wx[ijk   ]- wmean[k  ]) ) * dzi[k  ] -
                                        ( interp2(u[ijk   ]-umean[k  ], u[ijk-kk   ]-umean[k-1]) *  (wx[ijk   ]- wmean[k  ]) -
                                          interp2(u[ijk-kk]-umean[k-1], u[ijk-kk-kk]-umean[k-2]) *  (wx[ijk-kk]- wmean[k-1]) ) * dzi[k-1] ) * dzhi[k];
            }

        // Interior
        #pragma omp parallel for
        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // visc * d/dz(dw^2/dz)
                    w2_visc[ijk] = visc * ( (pow(w[ijk+kk]- wmean[k+1], 2) - pow(w[ijk   ]- wmean[k  ], 2)) * dzi[k  ] -
                                            (pow(w[ijk   ]- wmean[k  ], 2) - pow(w[ijk-kk]- wmean[k-1], 2)) * dzi[k-1] ) * dzhi[k];

                    // visc * d/dz(duw/dz)
                    uw_visc[ijk] = visc * ( ( interp2(u[ijk   ]-umean[k  ], u[ijk+kk   ]-umean[k+1]) * (wx[ijk+kk]- wmean[k+1]) -
                                              interp2(u[ijk   ]-umean[k  ], u[ijk-kk   ]-umean[k-1]) * (wx[ijk   ]- wmean[k  ]) ) * dzi[k  ] -
                                            ( interp2(u[ijk   ]-umean[k  ], u[ijk-kk   ]-umean[k-1]) * (wx[ijk   ]- wmean[k  ]) -
                                              interp2(u[ijk-kk]-umean[k-1], u[ijk-kk-kk]-umean[k-2]) * (wx[ijk-kk]- wmean[k-1]) ) * dzi[k-1] ) * dzhi[k];
                }
    }

    /**
     * Calculate the budget terms related to dissipation (-2*nu*(dui/dxj)^2)
     */
    template<typename TF>
    void calc_diffusion_dissipation_terms_dns(
            TF* const restrict u2_diss, TF* const restrict v2_diss,
            TF* const restrict w2_diss, TF* const restrict tke_diss, TF* const restrict uw_diss,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict umean, const TF* const restrict vmean, const TF* const restrict wmean,
            const TF* const restrict dzi, const TF* const restrict dzhi,
            const TF dxi, const TF dyi, const TF visc,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        // Dissipation term (-2*nu*(dui/dxj)^2)
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // -2 * visc * ((du/dx)^2 + (du/dy)^2 + (du/dz)^2)
                    u2_diss[ijk] = TF(-2.) * visc * ( pow( (interp2(u[ijk]-umean[k], u[ijk+ii]-umean[k  ]) - interp2(u[ijk]-umean[k], u[ijk-ii]-umean[k  ])) * dxi,    2) +
                                                      pow( (interp2(u[ijk]-umean[k], u[ijk+jj]-umean[k  ]) - interp2(u[ijk]-umean[k], u[ijk-jj]-umean[k  ])) * dyi,    2) +
                                                      pow( (interp2(u[ijk]-umean[k], u[ijk+kk]-umean[k+1]) - interp2(u[ijk]-umean[k], u[ijk-kk]-umean[k-1])) * dzi[k], 2) );

                    // -2 * visc * ((dv/dx)^2 + (dv/dy)^2 + (dv/dz)^2)
                    v2_diss[ijk] = TF(-2.) * visc * ( pow( (interp2(v[ijk]-vmean[k], v[ijk+ii]-vmean[k  ]) - interp2(v[ijk]-vmean[k], v[ijk-ii]-vmean[k  ])) * dxi,    2) +
                                                      pow( (interp2(v[ijk]-vmean[k], v[ijk+jj]-vmean[k  ]) - interp2(v[ijk]-vmean[k], v[ijk-jj]-vmean[k  ])) * dyi,    2) +
                                                      pow( (interp2(v[ijk]-vmean[k], v[ijk+kk]-vmean[k+1]) - interp2(v[ijk]-vmean[k], v[ijk-kk]-vmean[k-1])) * dzi[k], 2) );

                    // -2 * visc * ((dw/dx)^2 + (dw/dy)^2 + (dw/dz)^2)
                    tke_diss[ijk] = - visc * ( pow( (w[ijk+ii] - w[ijk]) * dxi,    2) +
                                               pow( (w[ijk+jj] - w[ijk]) * dyi,    2) +
                                               pow( (w[ijk+kk] - wmean[k+1] - ( w[ijk]-wmean[k])) * dzi[k], 2) )
                                  + TF(0.5) * (u2_diss[ijk] + v2_diss[ijk]);

                    // -2 * visc * du/dx * dw/dx
                    uw_diss[ijk] = TF(-2.) * visc * ( interp22(u[ijk]-umean[k], u[ijk+ii]-umean[k], u[ijk+ii-kk]-umean[k-1], u[ijk-kk]-umean[k-1]) -
                                                      interp22(u[ijk]-umean[k], u[ijk-ii]-umean[k], u[ijk-ii-kk]-umean[k-1], u[ijk-kk]-umean[k-1]) ) * dxi *
                                                    ( w[ijk] - w[ijk-ii] ) * dxi;

                    // -2 * visc * du/dy * dw/dy
                    uw_diss[ijk] += TF(-2.) * visc * ( interp22(u[ijk]-umean[k], u[ijk+jj]-umean[k], u[ijk+jj-kk]-umean[k-1], u[ijk-kk]-umean[k-1]) -
                                                       interp22(u[ijk]-umean[k], u[ijk-jj]-umean[k], u[ijk-jj-kk]-umean[k-1], u[ijk-kk]-umean[k-1]) ) * dyi *
                                                     ( interp22(w[ijk]         , w[ijk+jj]         , w[ijk+jj-ii]           , w[ijk-ii]           ) -
                                                       interp22(w[ijk]         , w[ijk-jj]         , w[ijk-jj-ii]           , w[ijk-ii]           ) ) * dyi;
                }
        }

        // Bottom boundary (z=0)
        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // -2 * visc * ((dw/dx)^2 + (dw/dy)^2 + (dw/dz)^2)
                // w @ full level kstart-1 = -w @ full level kstart+1
                w2_diss[ijk] = TF(-2.) * visc * ( pow( (interp2(w[ijk], w[ijk+ii]) - interp2(w[ijk], w[ijk-ii])) * dxi,     2) +
                                                  pow( (interp2(w[ijk], w[ijk+jj]) - interp2(w[ijk], w[ijk-jj])) * dyi,     2) +
                                                  pow( (TF(2.)*interp2(w[ijk], w[ijk+kk])                           ) * dzhi[k], 2) );

                // -2 * visc * du/dz * dw/dz
                // w @ full level kstart-1 = -w @ full level kstart+1
                uw_diss[ijk] = TF(-2.) * visc * ((u[ijk]-umean[k]) - (u[ijk-kk]-umean[k-1])) * dzhi[k] *
                                         TF(2.)*interp22(w[ijk]-wmean[k], w[ijk+kk]-wmean[k+1], w[ijk+kk-ii]-wmean[k+1], w[ijk-ii]-wmean[k]) * dzhi[k];
            }

        // Top boundary (z=zsize)
        k = kend;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                // -2 * visc * ((dw/dx)^2 + (dw/dy)^2 + (dw/dz)^2)
                // w @ full level kend = -w @ full level kend-1
                w2_diss[ijk] = TF(-2.) * visc * ( pow( (interp2(w[ijk], w[ijk+ii]) - interp2(w[ijk], w[ijk-ii])) * dxi,     2) +
                                                  pow( (interp2(w[ijk], w[ijk+jj]) - interp2(w[ijk], w[ijk-jj])) * dyi,     2) +
                                                  pow( (                   TF(-2.) * interp2(w[ijk]-wmean[k], w[ijk-kk]-wmean[k-1])) * dzhi[k], 2) );

                // -2 * visc * du/dz * dw/dz
                // w @ full level kend = - w @ full level kend-1
                uw_diss[ijk] = TF(-2.) * visc * ((u[ijk]-umean[k]) - (u[ijk-kk]-umean[k-1])) * dzhi[k] *
                                         - TF(2.)*interp22(w[ijk]-wmean[k], w[ijk-kk]-wmean[k-1], w[ijk-kk-ii]-wmean[k-1], w[ijk-ii]-wmean[k]) * dzhi[k];
            }

        // Interior
        #pragma omp parallel for
        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // -2 * visc * ((dw/dx)^2 + (dw/dy)^2 + (dw/dz)^2)
                    w2_diss[ijk] = TF(-2.) * visc * ( pow( (interp2(w[ijk], w[ijk+ii]) - interp2(w[ijk], w[ijk-ii])) * dxi,     2) +
                                                      pow( (interp2(w[ijk], w[ijk+jj]) - interp2(w[ijk], w[ijk-jj])) * dyi,     2) +
                                                      pow( (interp2(w[ijk]-wmean[k], w[ijk+kk]-wmean[k+1]) - interp2(w[ijk]-wmean[k], w[ijk-kk]-wmean[k-1])) * dzhi[k], 2) );

                    // -2 * visc * du/dz * dw/dz
                    uw_diss[ijk] = TF(-2.) * visc * ((u[ijk]-umean[k]) - (u[ijk-kk]-umean[k-1])) * dzhi[k] *
                                             ( interp22(w[ijk]-wmean[k], w[ijk+kk]-wmean[k+1], w[ijk+kk-ii]-wmean[k+1], w[ijk-ii]-wmean[k]) -
                                               interp22(w[ijk]-wmean[k], w[ijk-kk]-wmean[k-1], w[ijk-kk-ii]-wmean[k-1], w[ijk-ii]-wmean[k]) ) * dzhi[k];
                }
    }

    /**
     * Calculate the budget terms related to diffusion. In the approach here, we catch
     * molecular diffusion (nu*d/dxj(dui^2/dxj)) and dissipation (-2*nu*(dui/dxj)^2) in a
     * single term in order to ensure a closing budget.
     */
    template<typename TF>
    void calc_diffusion_terms_les(
            TF* const restrict u2_diff, TF* const restrict v2_diff,
            TF* const restrict w2_diff, TF* const restrict tke_diff,
            TF* const restrict uw_diff, TF* const restrict vw_diff,
            TF* const restrict wz, TF* const restrict evisch,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict ufluxbot, const TF* const restrict vfluxbot,
            const TF* const restrict evisc,
            const TF* const restrict umean, const TF* const restrict vmean, const TF* const restrict wmean,
            const TF* const restrict dzi, const TF* const restrict dzhi,
            const TF dxi, const TF dyi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int jcells, const int ijcells)
    {
        const int ii = 1;
        const int ii2 = 2;
        const int jj = icells;
        const int jj2 = 2*icells;
        const int kk = ijcells;
        const int kk2 = 2*ijcells;

        // Calculate w at full levels (center)
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    wz[ijk] = interp2(w[ijk]-wmean[k], w[ijk+kk]-wmean[k+1]);
                }

        // Set ghost cells such that the velocity interpolated to the boundaries is zero
        int ks = kstart;
        int ke = kend-1;
        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ijks = i + j*jj + ks*kk;
                const int ijke = i + j*jj + ke*kk;
                wz[ijks-kk] = -wz[ijks];
                wz[ijke+kk] = -wz[ijke];
            }

        // Calculate evisc at half-half-half level
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    evisch[ijk] = TF(0.125) * (evisc[ijk-ii-jj-kk] + evisc[ijk-ii-jj] + evisc[ijk-ii-kk] + evisc[ijk-ii] +
                                               evisc[ijk   -jj-kk] + evisc[ijk   -jj] + evisc[ijk   -kk] + evisc[ijk   ]);
                }

        // boundary_cyclic(evisch);

        // -----------------------------
        // Test: directly calculate diffusion terms as 2 ui * d/dxj(visc * dui/dx + visc * duj/dxi)
        // Term is stored in xx_diss; xx_visc=0
        // -----------------------------
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    const TF evisc_utop   = interp22(evisc[ijk], evisc[ijk+kk], evisc[ijk-ii+kk], evisc[ijk-ii]);
                    const TF evisc_ubot   = interp22(evisc[ijk], evisc[ijk-kk], evisc[ijk-ii-kk], evisc[ijk-ii]);
                    const TF evisc_unorth = interp22(evisc[ijk], evisc[ijk+jj], evisc[ijk+jj-ii], evisc[ijk-ii]);
                    const TF evisc_usouth = interp22(evisc[ijk], evisc[ijk-jj], evisc[ijk-jj-ii], evisc[ijk-ii]);

                    const TF evisc_vtop   = interp22(evisc[ijk], evisc[ijk+kk], evisc[ijk-jj+kk], evisc[ijk-jj]);
                    const TF evisc_vbot   = interp22(evisc[ijk], evisc[ijk-kk], evisc[ijk-jj-kk], evisc[ijk-jj]);
                    const TF evisc_veast  = interp22(evisc[ijk], evisc[ijk+ii], evisc[ijk+ii-jj], evisc[ijk-jj]);
                    const TF evisc_vwest  = interp22(evisc[ijk], evisc[ijk-ii], evisc[ijk-ii-jj], evisc[ijk-jj]);

                    const TF evisc_weast  = interp22(evisc[ijk], evisc[ijk+ii], evisc[ijk+ii-kk], evisc[ijk-kk]);
                    const TF evisc_wwest  = interp22(evisc[ijk], evisc[ijk-ii], evisc[ijk-ii-kk], evisc[ijk-kk]);
                    const TF evisc_wnorth = interp22(evisc[ijk], evisc[ijk+jj], evisc[ijk+jj-kk], evisc[ijk-kk]);
                    const TF evisc_wsouth = interp22(evisc[ijk], evisc[ijk-jj], evisc[ijk-jj-kk], evisc[ijk-kk]);

                    // -----------------------------------------
                    // 2 * u * d/dx( visc * du/dx + visc * du/dx )
                    u2_diff[ijk] = TF(2.) * (u[ijk]-umean[k]) * ( evisc[ijk   ] * (u[ijk+ii] - u[ijk   ]) * dxi -
                                                                  evisc[ijk-ii] * (u[ijk   ] - u[ijk-ii]) * dxi ) * TF(2.) * dxi;

                    // 2 * u * d/dy( visc * du/dy + visc * dv/dx)
                    u2_diff[ijk] += TF(2.) * (u[ijk]-umean[k]) * ( evisc_unorth * (u[ijk+jj] - u[ijk      ]) * dyi -
                                                                   evisc_usouth * (u[ijk   ] - u[ijk-jj   ]) * dyi +
                                                                   evisc_unorth * (v[ijk+jj] - v[ijk+jj-ii]) * dxi -
                                                                   evisc_usouth * (v[ijk   ] - v[ijk-ii   ]) * dxi ) * dyi;

                    // 2 * u * d/dz( visc * dw/dx )
                    u2_diff[ijk] += TF(2.) * (u[ijk]-umean[k]) * ( evisc_utop * (w[ijk+kk] - w[ijk-ii+kk]) * dxi -
                                                                   evisc_ubot * (w[ijk   ] - w[ijk-ii   ]) * dxi ) * dzi[k];

                    // -----------------------------------------
                    // 2 * v * d/dy( visc * dv/dy + visc * dv/dy )
                    v2_diff[ijk] = TF(2.) * (v[ijk]-vmean[k]) * ( evisc[ijk   ] * (v[ijk+jj] - v[ijk   ]) * dyi -
                                                                  evisc[ijk-jj] * (v[ijk   ] - v[ijk-jj]) * dyi ) * TF(2.) * dyi;

                    // 2 * v * d/dx( visc * dv/dx + visc * du/dy )
                    v2_diff[ijk] += TF(2.) * (v[ijk]-vmean[k]) * ( evisc_veast * (v[ijk+ii] - v[ijk      ]) * dxi -
                                                                   evisc_vwest * (v[ijk   ] - v[ijk-ii   ]) * dxi +
                                                                   evisc_veast * (u[ijk+ii] - u[ijk+ii-jj]) * dyi -
                                                                   evisc_vwest * (u[ijk   ] - u[ijk-jj   ]) * dyi ) * dxi;

                    // 2 * v * d/dz( visc * dw/dy )
                    v2_diff[ijk] += TF(2.) * (v[ijk]-vmean[k]) * ( evisc_vtop * (w[ijk+kk] - w[ijk-jj+kk]) * dyi -
                                                                   evisc_vbot * (w[ijk   ] - w[ijk-jj   ]) * dyi ) * dzi[k];

                    // -----------------------------------------
                    // 2 * w * d/dx( visc * dw/dx )
                    w2_diff[ijk] = TF(2.) * (w[ijk]-wmean[k]) * ( evisc_weast * (w[ijk+ii] - w[ijk   ]) * dxi -
                                                       evisc_wwest * (w[ijk   ] - w[ijk-ii]) * dxi ) * dxi;

                    // 2 * w * d/dy( visc * dw/dy )
                    w2_diff[ijk] += TF(2.) * (w[ijk]-wmean[k]) * ( evisc_wnorth * (w[ijk+jj] - w[ijk   ]) * dyi -
                                                        evisc_wsouth * (w[ijk   ] - w[ijk-jj]) * dyi ) * dyi;

                    // -----------------------------------------
                    // 2 * w * d/dx( visc * dw/dx )
                    tke_diff[ijk] = wz[ijk] * ( interp2(evisc[ijk], evisc[ijk+ii]) * (wz[ijk+ii] - wz[ijk   ]) * dxi -
                                                interp2(evisc[ijk], evisc[ijk-ii]) * (wz[ijk   ] - wz[ijk-ii]) * dxi ) * dxi;

                    // 2 * w * d/dx( visc * du/dz )
                    tke_diff[ijk] += (wz[ijk]-wmean[k]) * ( interp2(evisc[ijk], evisc[ijk+ii]) * (interp2(u[ijk+ii], u[ijk+ii+kk]) - interp2(u[ijk+ii], u[ijk+ii-kk])) * dzi[k] -
                                                 interp2(evisc[ijk], evisc[ijk-ii]) * (interp2(u[ijk   ], u[ijk   +kk]) - interp2(u[ijk   ], u[ijk   -kk])) * dzi[k] ) * dxi;

                    // 2 * w * d/dy( visc * dw/dy )
                    tke_diff[ijk] += (wz[ijk]-wmean[k]) * ( interp2(evisc[ijk], evisc[ijk+jj]) * (wz[ijk+jj] - wz[ijk   ]) * dyi -
                                                 interp2(evisc[ijk], evisc[ijk-jj]) * (wz[ijk   ] - wz[ijk-jj]) * dyi ) * dyi;

                    // 2 * w * d/dy( visc * dv/dz )
                    tke_diff[ijk] += (wz[ijk]-wmean[k]) * ( interp2(evisc[ijk], evisc[ijk+jj]) * (interp2(v[ijk+jj], v[ijk+jj+kk]) - interp2(v[ijk+jj], v[ijk+jj-kk])) * dzi[k] -
                                                 interp2(evisc[ijk], evisc[ijk-jj]) * (interp2(v[ijk   ], v[ijk   +kk]) - interp2(v[ijk   ], v[ijk   -kk])) * dzi[k] ) * dyi;
                }

        #pragma omp parallel for
        for (int k=kstart+1; k<kend; ++k)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    const TF evisc_utop   = interp22(evisc[ijk], evisc[ijk+kk], evisc[ijk-ii+kk], evisc[ijk-ii]);
                    const TF evisc_ubot   = interp22(evisc[ijk], evisc[ijk-kk], evisc[ijk-ii-kk], evisc[ijk-ii]);

                    const TF evisc_vtop   = interp22(evisc[ijk], evisc[ijk+kk], evisc[ijk-jj+kk], evisc[ijk-jj]);
                    const TF evisc_vbot   = interp22(evisc[ijk], evisc[ijk-kk], evisc[ijk-jj-kk], evisc[ijk-jj]);

                    const TF evisc_weast  = interp22(evisc[ijk], evisc[ijk+ii], evisc[ijk+ii-kk], evisc[ijk-kk]);
                    const TF evisc_wwest  = interp22(evisc[ijk], evisc[ijk-ii], evisc[ijk-ii-kk], evisc[ijk-kk]);
                    const TF evisc_wnorth = interp22(evisc[ijk], evisc[ijk+jj], evisc[ijk+jj-kk], evisc[ijk-kk]);
                    const TF evisc_wsouth = interp22(evisc[ijk], evisc[ijk-jj], evisc[ijk-jj-kk], evisc[ijk-kk]);

                    // -----------------------------------------
                    // 2 * u * d/dz( visc * du/dz )
                    u2_diff[ijk] += TF(2.) * (u[ijk]-umean[k]) * ( evisc_utop * (u[ijk+kk] - u[ijk   ]) * dzhi[k+1] -
                                                                  evisc_ubot * (u[ijk   ] - u[ijk-kk]) * dzhi[k  ] ) * dzi[k];

                    // -----------------------------------------
                    // 2 * v * d/dz( visc * dv/dz )
                    v2_diff[ijk] += TF(2.) * (v[ijk]-vmean[k]) * ( evisc_vtop * (v[ijk+kk] - v[ijk   ]) * dzhi[k+1] -
                                                                  evisc_vbot * (v[ijk   ] - v[ijk-kk]) * dzhi[k  ] ) * dzi[k];

                    // -----------------------------------------
                    // 2 * w * d/dx( visc * du/dz )
                    w2_diff[ijk] += TF(2.) * (w[ijk]-wmean[k]) * ( evisc_weast * (u[ijk+ii] - u[ijk+ii-kk]) * dzhi[k] -
                                                       evisc_wwest * (u[ijk   ] - u[ijk   -kk]) * dzhi[k] ) * dxi;

                    // 2 * w * d/dy( visc * dv/dz )
                    w2_diff[ijk] += TF(2.) * (w[ijk]-wmean[k]) * ( evisc_wnorth * (v[ijk+jj] - v[ijk+jj-kk]) * dzhi[k] -
                                                        evisc_wsouth * (v[ijk   ] - v[ijk   -kk]) * dzhi[k] ) * dyi;

                    // 2 * w * d/dz( visc * dw/dz )
                    w2_diff[ijk] += TF(2.) * (w[ijk]-wmean[k]) * ( evisc[ijk   ] * (w[ijk+kk] - w[ijk   ]) * dzi[k  ] -
                                                        evisc[ijk-kk] * (w[ijk   ] - w[ijk-kk]) * dzi[k-1] ) * TF(2.) * dzhi[k];

                    // -----------------------------------------
                    // 2 * w * d/dz( 2 * visc * dw/dz )
                    tke_diff[ijk] += wz[ijk] * ( interp2(evisc[ijk], evisc[ijk+kk]) * (wz[ijk+kk] - wz[ijk   ]) * dzhi[k+1] -
                                                 interp2(evisc[ijk], evisc[ijk-kk]) * (wz[ijk   ] - wz[ijk-kk]) * dzhi[k  ] ) * TF(2.) * dzi[k]
                                   + TF(0.5) * (u2_diff[ijk] + v2_diff[ijk]);

                    // -----------------------------------------
                    // w * d/dx(visc * du/dx + visc * du/dx)
                    uw_diff[ijk] = ( ( interp2(w[ijk-ii], w[ijk    ])
                                      * ( ( ( ( TF(2.) * interp2(evisc[ijk    -kk], evisc[ijk        ]) )
                                          * ( interp2(u[ijk+ii-kk], u[ijk+ii    ]) - interp2(u[ijk    -kk], u[ijk        ]) ) )
                                        * dxi ) - ( ( ( TF(2.) * interp2(evisc[ijk-ii-kk], evisc[ijk-ii    ]) )
                                          * ( interp2(u[ijk    -kk], u[ijk        ]) - interp2(u[ijk-ii-kk], u[ijk-ii    ]) ) )
                                        * dxi ) ) )
                                    * dxi );

                    // w * d/dy(visc * du/dy + visc * dv/dx)
                    uw_diff[ijk] += ( ( interp2(w[ijk-ii], w[ijk    ])
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
                    uw_diff[ijk] += ( ( interp2(w[ijk-ii], w[ijk    ])
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
                    uw_diff[ijk] += ( ( interp2(u[ijk-kk], u[ijk    ])
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
                    uw_diff[ijk] += ( ( interp2(u[ijk-kk], u[ijk    ])
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
                    uw_diff[ijk] += ( ( interp2(u[ijk-kk], u[ijk    ])
                                      * ( ( ( 2 * interp2(evisc[ijk-ii    ], evisc[ijk        ]) )
                                        * ( ( interp2(w[ijk-ii+kk], w[ijk    +kk]) - interp2(w[ijk-ii    ], w[ijk        ]) )
                                          * dzi[k  ] ) ) - ( ( TF(2.) * interp2(evisc[ijk-ii-kk], evisc[ijk    -kk]) )
                                        * ( ( interp2(w[ijk-ii    ], w[ijk        ]) - interp2(w[ijk-ii-kk], w[ijk    -kk]) )
                                          * dzi[k-1] ) ) ) )
                                    * dzhi[k] );

                    // ------------------------------------------------
                    // w * d/dx(visc * dv/dx + visc * du/dy)
                    vw_diff[ijk] = ( ( interp2(w[ijk-jj], w[ijk    ])
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
                    vw_diff[ijk] += ( ( interp2(w[ijk-jj], w[ijk    ])
                                    * ( ( ( TF(2.) * interp2(evisc[ijk    -kk], evisc[ijk        ]) )
                                      * ( interp2(v[ijk+jj-kk], v[ijk+jj    ]) - interp2(v[ijk    -kk], v[ijk        ]) ) ) - ( ( TF(2.) * interp2(evisc[ijk-jj-kk], evisc[ijk-jj    ]) )
                                      * ( interp2(v[ijk    -kk], v[ijk        ]) - interp2(v[ijk-jj-kk], v[ijk-jj    ]) ) ) ) )
                                  * dyi );

                    // w * d/dz(visc * du/dz + visc * dw/dx)
                    vw_diff[ijk] += ( ( interp2(w[ijk-jj], w[ijk    ])
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
                    vw_diff[ijk] += ( ( interp2(v[ijk-kk], v[ijk    ])
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
                    vw_diff[ijk] += ( ( interp2(v[ijk-kk], v[ijk    ])
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
                    vw_diff[ijk] += ( ( interp2(v[ijk-kk], v[ijk    ])
                                    * ( ( ( TF(2.) * interp2(evisc[ijk-jj    ], evisc[ijk        ]) )
                                      * ( ( interp2(w[ijk-jj+kk], w[ijk    +kk]) - interp2(w[ijk-jj    ], w[ijk        ]) )
                                        * dzi[k  ] ) ) - ( ( TF(2.) * interp2(evisc[ijk-jj-kk], evisc[ijk    -kk]) )
                                      * ( ( interp2(w[ijk-jj    ], w[ijk        ]) - interp2(w[ijk-jj-kk], w[ijk    -kk]) )
                                        * dzi[k-1] ) ) ) )
                                  * dzhi[k] );
                }
        }

        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                const int ij  = i + j*jj;

                const TF evisc_utop   = interp22(evisc[ijk], evisc[ijk+kk], evisc[ijk-ii+kk], evisc[ijk-ii]);
                const TF evisc_vtop   = interp22(evisc[ijk], evisc[ijk+kk], evisc[ijk-jj+kk], evisc[ijk-jj]);

                // 2 u * d/dz( visc * du/dz )
                u2_diff[ijk] = TF(2.) * (u[ijk]-umean[k]) * ( evisc_utop * (u[ijk+kk] - u[ijk   ]) * dzhi[k+1] + ufluxbot[ij]) * dzi[k];

                // 2 v * d/dz( visc * dv/dz )
                v2_diff[ijk] = TF(2.) * (v[ijk]-vmean[k]) * ( evisc_vtop * (v[ijk+kk] - v[ijk   ]) * dzhi[k+1] + vfluxbot[ij]) * dzi[k];

                // 2 * w * d/dz( visc * dw/dz )
                // What to do with evisc at surface (term visc * dw/dz at surface)?
                tke_diff[ijk] = (wz[ijk]-wmean[k]) * ( interp2(evisc[ijk], evisc[ijk+kk]) * (wz[ijk+kk] - wz[ijk   ]) * dzhi[k+1] ) * TF(2.) * dzi[k]
                              + TF(0.5) * (u2_diff[ijk] + v2_diff[ijk]);

                // uw_diff is zero at surface for no-slip case, unequal for free-slip...
            }
    }

    /**
     * Calculate the budget related to buoyancy.
     */
    template<typename TF>
    void calc_buoyancy_terms(
            TF* const restrict w2_buoy, TF* const restrict tke_buoy,
            TF* const restrict uw_buoy, TF* const restrict vw_buoy,
            const TF* const restrict u, const TF* const restrict v,
            const TF* const restrict w, const TF* const restrict b,
            const TF* const restrict umean, const TF* const restrict vmean, const TF* const restrict wmean,
            const TF* const restrict bmean,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // w'b'
                    tke_buoy[ijk] = interp2(w[ijk]-wmean[k], w[ijk+kk]-wmean[k+1]) * (b[ijk] - bmean[k]);
                }

        #pragma omp parallel for
        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // w'b'
                    w2_buoy[ijk] = TF(2.) * interp2(b[ijk]-bmean[k], b[ijk-kk]-bmean[k-1]) * (w[ijk]-wmean[k]);

                    // u'b'
                    uw_buoy[ijk] = interp2 (u[ijk]-umean[k], u[ijk-kk]-umean[k-1]) *
                                   interp22(b[ijk]-bmean[k], b[ijk-ii]-bmean[k], b[ijk-ii-kk]-bmean[k-1], b[ijk-kk]-bmean[k-1]);

                    // v'b'
                    vw_buoy[ijk] = interp2 (v[ijk]-vmean[k], v[ijk-kk]-vmean[k-1]) *
                                   interp22(b[ijk]-bmean[k], b[ijk-jj]-bmean[k], b[ijk-jj-kk]-bmean[k-1], b[ijk-kk]-bmean[k-1]);

                }
    }

    /**
     * Calculate the scalar budget terms arising from buoyancy
     */
    template<typename TF>
    void calc_buoyancy_terms_scalar(
            TF* const restrict sw_buoy,
            const TF* const restrict s, const TF* const restrict b,
            const TF* const restrict smean, const TF* const restrict bmean,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    sw_buoy[ijk] = interp2(s[ijk]-smean[k], s[ijk-kk]-smean[k-1]) * interp2(b[ijk]-bmean[k], b[ijk-kk]-bmean[k-1]);
                }
    }

    /**
     * Calculate the scalar budget terms arising from the advection term
     */
    template<typename TF>
    void calc_advection_terms_scalar(
            TF* const restrict s2_shear, TF* const restrict s2_turb,
            TF* const restrict sw_shear, TF* const restrict sw_turb,
            const TF* const restrict w, const TF* const restrict s,
            const TF* const restrict smean,
            const TF* const restrict dzi, const TF* const restrict dzhi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
        {
            const TF dsdz  = (interp2(smean[k], smean[k+1]) - interp2(smean[k], smean[k-1])) * dzi[k];
            const TF dsdzh = (smean[k] - smean[k-1]) * dzhi[k];

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    s2_shear[ijk] = - TF(2.) * (s[ijk] - smean[k]) * interp2(w[ijk], w[ijk+kk]) * dsdz;

                    s2_turb[ijk]  = - ((pow(interp2(s[ijk]-smean[k], s[ijk+kk]-smean[k+1]), 2) * w[ijk+kk]) -
                                       (pow(interp2(s[ijk]-smean[k], s[ijk-kk]-smean[k-1]), 2) * w[ijk   ])) * dzi[k];

                    sw_shear[ijk] = - pow(w[ijk], 2) * dsdzh;

                    sw_turb[ijk]  = - ((pow(interp2(w[ijk], w[ijk+kk]), 2) * (s[ijk   ]-smean[k  ]))-
                                       (pow(interp2(w[ijk], w[ijk-kk]), 2) * (s[ijk-kk]-smean[k-1]))) * dzhi[k];
                }
        }
    }

    /**
     * Calculate the budget terms arrising from diffusion, for a fixed viscosity
     * molecular diffusion (nu*d/dxj(dui^2/dxj)) and dissipation (-2*nu*(dui/dxj)^2)
     */
    template<typename TF>
    void calc_diffusion_terms_scalar_dns(
            TF* const restrict b2_visc, TF* const restrict b2_diss,
            TF* const restrict bw_visc, TF* const restrict bw_diss,
            const TF* const restrict w, const TF* const restrict b,
            const TF* const restrict bmean,
            const TF* const restrict dzi, const TF* const restrict dzhi,
            const TF dxi, const TF dyi, const TF visc, const TF diff,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    b2_visc[ijk] = diff * ( (std::pow(b[ijk+kk]-bmean[k+1], 2) - std::pow(b[ijk   ]-bmean[k  ], 2))*dzhi[k+1] -
                                            (std::pow(b[ijk   ]-bmean[k  ], 2) - std::pow(b[ijk-kk]-bmean[k-1], 2))*dzhi[k  ] ) * dzi[k];

                    b2_diss[ijk] = TF(-2.) * diff * (
                                               std::pow((interp2(b[ijk]-bmean[k], b[ijk+kk]-bmean[k+1]) - interp2(b[ijk]-bmean[k], b[ijk-kk]-bmean[k-1])) * dzi[k], 2) +
                                               std::pow((interp2(b[ijk]-bmean[k], b[ijk+ii]-bmean[k  ]) - interp2(b[ijk]-bmean[k], b[ijk-ii]-bmean[k  ])) * dxi,    2) +
                                               std::pow((interp2(b[ijk]-bmean[k], b[ijk+jj]-bmean[k  ]) - interp2(b[ijk]-bmean[k], b[ijk-jj]-bmean[k  ])) * dyi,    2)
                                             );
                }

        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                // with w[kstart-1] undefined, use gradient w over lowest point
                bw_diss[ijk] = TF(-2.) * visc * (w[ijk+kk]-w[ijk]) * dzi[k] * ((b[ijk]-bmean[k])-(b[ijk-kk]-bmean[k-1]))*dzhi[k];
            }

        k = kend;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                bw_diss[ijk] = TF(-2.) * visc * (w[ijk]-w[ijk-kk]) * dzi[k-1] * ((b[ijk]-bmean[k])-(b[ijk-kk]-bmean[k-1]))*dzhi[k];
            }

        #pragma omp parallel for
        for (int k=kstart+1; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    bw_visc[ijk] = visc * ( ( (w[ijk+kk] * interp2(b[ijk      ]-bmean[k  ], b[ijk+kk]-bmean[k+1])) -
                                              (w[ijk   ] * interp2(b[ijk-kk   ]-bmean[k-1], b[ijk   ]-bmean[k  ])) ) * dzi[k  ] -
                                            ( (w[ijk   ] * interp2(b[ijk-kk   ]-bmean[k-1], b[ijk   ]-bmean[k  ])) -
                                              (w[ijk-kk] * interp2(b[ijk-kk-kk]-bmean[k-2], b[ijk-kk]-bmean[k-1])) ) * dzi[k-1] ) * dzhi[k];

                    bw_diss[ijk] = TF(-2.) * visc * (
                                                (interp2(w[ijk+ii], w[ijk]) - interp2(w[ijk], w[ijk-ii])) * dxi *
                                                (interp22(b[ijk]-bmean[k], b[ijk+ii]-bmean[k], b[ijk+ii-kk]-bmean[k-1],b[ijk-kk]-bmean[k-1]) -
                                                 interp22(b[ijk]-bmean[k], b[ijk-ii]-bmean[k], b[ijk-ii-kk]-bmean[k-1],b[ijk-kk]-bmean[k-1])) * dxi +
                                                (interp2(w[ijk+jj], w[ijk]) - interp2(w[ijk], w[ijk-jj])) * dyi *
                                                (interp22(b[ijk]-bmean[k], b[ijk+jj]-bmean[k], b[ijk+jj-kk]-bmean[k-1],b[ijk-kk]-bmean[k-1]) -
                                                 interp22(b[ijk]-bmean[k], b[ijk-jj]-bmean[k], b[ijk-jj-kk]-bmean[k-1],b[ijk-kk]-bmean[k-1])) * dyi +
                                                (interp2(w[ijk+kk], w[ijk]) - interp2(w[ijk], w[ijk-kk])) * dzhi[k] *
                                                ((b[ijk]-bmean[k])-(b[ijk-kk]-bmean[k-1]))*dzhi[k]
                                             );
                }

        // The second derivative of the flux at the lower and top boundary can't be calculated; with a biased
        // second derivative the term at kstart and kend equals the term at kstart+1 and kend-1, respectively
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk_start = i + j*jj + kstart*kk;
                const int ijk_end = i + j*jj + kend*kk;

                bw_visc[ijk_start] = bw_visc[ijk_start+kk];
                bw_visc[ijk_end  ] = bw_visc[ijk_end  -kk];
            }
    }

    /**
     * Calculate the scalar budget terms arising from pressure
     */
    template<typename TF>
    void calc_pressure_terms_scalar(
            TF* const restrict sw_pres, TF* const restrict sw_rdstr,
            const TF* const restrict s, const TF* const restrict p,
            const TF* const restrict smean, const TF* const restrict pmean,
            const TF* const restrict dzi, const TF* const restrict dzhi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    sw_pres[ijk] = - ((p[ijk]-pmean[k]) * (s[ijk]-smean[k]) - (p[ijk-kk]-pmean[k-1]) * (s[ijk-kk]-smean[k-1])) * dzhi[k];
                    sw_rdstr[ijk] = interp2(p[ijk]-pmean[k], p[ijk-kk]-pmean[k-1]) * ((s[ijk]-smean[k])-(s[ijk-kk]-smean[k-1])) * dzhi[k];
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
    wmodel.resize(gd.kcells);
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

    stats.add_prof("w2_pres" , "Pressure transport term in W2 budget" , "m2 s-3", "zh", group_name);
    stats.add_prof("tke_pres", "Pressure transport term in TKE budget", "m2 s-3", "z" , group_name);
    stats.add_prof("uw_pres" , "Pressure transport term in UW budget" , "m2 s-3", "zh", group_name);
    stats.add_prof("vw_pres" , "Pressure transport term in VW budget" , "m2 s-3", "zh", group_name);

    stats.add_prof("u2_rdstr", "Pressure redistribution term in U2 budget", "m2 s-3", "z" , group_name);
    stats.add_prof("v2_rdstr", "Pressure redistribution term in V2 budget", "m2 s-3", "z" , group_name);
    stats.add_prof("w2_rdstr", "Pressure redistribution term in W2 budget", "m2 s-3", "zh", group_name);
    stats.add_prof("uw_rdstr", "Pressure redistribution term in UW budget", "m2 s-3", "zh", group_name);
    stats.add_prof("vw_rdstr", "Pressure redistribution term in VW budget", "m2 s-3", "zh", group_name);

    if (force.get_switch_lspres() == Large_scale_pressure_type::Geo_wind)
    {
        stats.add_prof("u2_cor", "Coriolis term in U2 budget", "m2 s-3", "z" , group_name);
        stats.add_prof("v2_cor", "Coriolis term in V2 budget", "m2 s-3", "z" , group_name);
        stats.add_prof("uw_cor", "Coriolis term in UW budget", "m2 s-3", "zh", group_name);
        stats.add_prof("vw_cor", "Coriolis term in VW budget", "m2 s-3", "zh", group_name);
    }

    if (diff.get_switch() != Diffusion_type::Disabled)
    {
        if (diff.get_switch() == Diffusion_type::Diff_2 || diff.get_switch() == Diffusion_type::Diff_4)
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
        }

        // For LES, add only the total diffusive budget terms, which (unlike diss + visc) close
        else if (diff.get_switch() == Diffusion_type::Diff_smag2)
        {
            stats.add_prof("u2_diff" , "Total diffusive term in U2 budget" , "m2 s-3", "z" , group_name);
            stats.add_prof("v2_diff" , "Total diffusive term in V2 budget" , "m2 s-3", "z" , group_name);
            stats.add_prof("w2_diff" , "Total diffusive term in W2 budget" , "m2 s-3", "zh", group_name);
            stats.add_prof("tke_diff", "Total diffusive term in TKE budget", "m2 s-3", "z" , group_name);
            stats.add_prof("uw_diff" , "Total diffusive term in UW budget" , "m2 s-3", "zh", group_name);
            stats.add_prof("vw_diff" , "Total diffusive term in VW budget" , "m2 s-3", "zh", group_name);
        }
    }

    if (thermo.get_switch() != Thermo_type::Disabled)
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
}

template<typename TF>
void Budget_2<TF>::exec_stats(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    auto& masks = stats.get_masks();

    // The loop over masks inside of budget is necessary, because the mask mean is 
    // required in order to compute the budget terms.
    for (auto& m : masks)
    {
        // Calculate the mean of the fields.
        stats.calc_mask_mean_profile(umodel, m, *fields.mp.at("u"));
        stats.calc_mask_mean_profile(vmodel, m, *fields.mp.at("v"));
        stats.calc_mask_mean_profile(wmodel, m, *fields.mp.at("w"));

        // field3d_operators.calc_mean_profile(umodel.data(), fields.mp.at("u")->fld.data());
        // field3d_operators.calc_mean_profile(vmodel.data(), fields.mp.at("v")->fld.data());

        // Calculate kinetic and turbulent kinetic energy
        auto ke  = fields.get_tmp();
        auto tke = fields.get_tmp();

        constexpr TF no_offset = 0.;
        constexpr TF no_threshold = 0.;

        calc_kinetic_energy(
                ke->fld.data(), tke->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                umodel.data(), vmodel.data(), wmodel.data(),
                gd.utrans, gd.vtrans,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_mask_stats(m, "ke" , *ke , no_offset, no_threshold);
        stats.calc_mask_stats(m, "tke", *tke, no_offset, no_threshold);

        auto wx = std::move(ke );
        auto wy = std::move(tke);

        // Interpolate w to the locations of u and v.
        constexpr int wloc [3] = {0,0,1};
        constexpr int wxloc[3] = {1,0,1};
        constexpr int wyloc[3] = {0,1,1};

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
                umodel.data(), vmodel.data(), wmodel.data(),
                wx->fld.data(), wy->fld.data(),
                gd.dzi.data(), gd.dzhi.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_mask_stats(m, "u2_shear" , *u2_shear , no_offset, no_threshold);
        stats.calc_mask_stats(m, "v2_shear" , *v2_shear , no_offset, no_threshold);
        stats.calc_mask_stats(m, "tke_shear", *tke_shear, no_offset, no_threshold);
        stats.calc_mask_stats(m, "uw_shear" , *uw_shear , no_offset, no_threshold);
        stats.calc_mask_stats(m, "vw_shear" , *vw_shear , no_offset, no_threshold);

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
                umodel.data(), vmodel.data(), wmodel.data(),
                wx->fld.data(), wy->fld.data(),
                gd.dzi.data(), gd.dzhi.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_mask_stats(m, "u2_turb" , *u2_turb , no_offset, no_threshold);
        stats.calc_mask_stats(m, "v2_turb" , *v2_turb , no_offset, no_threshold);
        stats.calc_mask_stats(m, "w2_turb" , *w2_turb , no_offset, no_threshold);
        stats.calc_mask_stats(m, "tke_turb", *tke_turb, no_offset, no_threshold);
        stats.calc_mask_stats(m, "uw_turb" , *uw_turb , no_offset, no_threshold);
        stats.calc_mask_stats(m, "vw_turb" , *vw_turb , no_offset, no_threshold);

        fields.release_tmp(u2_turb);
        fields.release_tmp(v2_turb);
        fields.release_tmp(w2_turb);
        fields.release_tmp(tke_turb);
        fields.release_tmp(uw_turb);
        fields.release_tmp(vw_turb);

        if (diff.get_switch() != Diffusion_type::Disabled)
        {
            // Calculate the diffusive transport and dissipation terms
            if (diff.get_switch() == Diffusion_type::Diff_2 || diff.get_switch() == Diffusion_type::Diff_4)
            {
                auto u2_visc = fields.get_tmp();
                auto v2_visc = fields.get_tmp();
                auto w2_visc = fields.get_tmp();
                auto tke_visc = fields.get_tmp();
                auto uw_visc = fields.get_tmp();

                auto wz = fields.get_tmp();

                calc_diffusion_transport_terms_dns(
                        u2_visc->fld.data(), v2_visc->fld.data(), w2_visc->fld.data(), tke_visc->fld.data(), uw_visc->fld.data(),
                        fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                        wx->fld.data(), wy->fld.data(), wz->fld.data(),
                        umodel.data(), vmodel.data(), wmodel.data(),
                        gd.dzi.data(), gd.dzhi.data(), gd.dxi, gd.dyi, fields.visc,
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                        gd.icells, gd.ijcells);

                stats.calc_mask_stats(m, "u2_visc" , *u2_visc , no_offset, no_threshold);
                stats.calc_mask_stats(m, "v2_visc" , *v2_visc , no_offset, no_threshold);
                stats.calc_mask_stats(m, "w2_visc" , *w2_visc , no_offset, no_threshold);
                stats.calc_mask_stats(m, "tke_visc", *tke_visc, no_offset, no_threshold);
                stats.calc_mask_stats(m, "uw_visc" , *uw_visc , no_offset, no_threshold);

                auto u2_diss = std::move(u2_visc);
                auto v2_diss = std::move(v2_visc);
                auto w2_diss = std::move(w2_visc);
                auto tke_diss = std::move(tke_visc);
                auto uw_diss = std::move(uw_visc);

                calc_diffusion_dissipation_terms_dns(
                        u2_diss->fld.data(), v2_diss->fld.data(), w2_diss->fld.data(), tke_diss->fld.data(), uw_diss->fld.data(),
                        fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                        umodel.data(), vmodel.data(), wmodel.data(),
                        gd.dzi.data(), gd.dzhi.data(), gd.dxi, gd.dyi, fields.visc,
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                        gd.icells, gd.ijcells);

                stats.calc_mask_stats(m, "u2_diss" , *u2_diss , no_offset, no_threshold);
                stats.calc_mask_stats(m, "v2_diss" , *v2_diss , no_offset, no_threshold);
                stats.calc_mask_stats(m, "w2_diss" , *w2_diss , no_offset, no_threshold);
                stats.calc_mask_stats(m, "tke_diss", *tke_diss, no_offset, no_threshold);
                stats.calc_mask_stats(m, "uw_diss" , *uw_diss , no_offset, no_threshold);

                fields.release_tmp(u2_diss);
                fields.release_tmp(v2_diss);
                fields.release_tmp(w2_diss);
                fields.release_tmp(tke_diss);
                fields.release_tmp(uw_diss);

                fields.release_tmp(wz);
            }

            else if (diff.get_switch() == Diffusion_type::Diff_smag2)
            {
                auto u2_diff = fields.get_tmp();
                auto v2_diff = fields.get_tmp();
                auto w2_diff = fields.get_tmp();
                auto tke_diff = fields.get_tmp();
                auto uw_diff = fields.get_tmp();
                auto vw_diff = fields.get_tmp();
                auto wz = fields.get_tmp();
                auto evisch = fields.get_tmp();
 
                calc_diffusion_terms_les(
                        u2_diff->fld.data(), v2_diff->fld.data(),
                        w2_diff->fld.data(), tke_diff->fld.data(),
                        uw_diff->fld.data(), vw_diff->fld.data(),
                        wz->fld.data(), evisch->fld.data(),
                        fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                        fields.mp.at("u")->flux_bot.data(), fields.mp.at("v")->flux_bot.data(),
                        fields.sd.at("evisc")->fld.data(),
                        umodel.data(), vmodel.data(), wmodel.data(),
                        gd.dzi.data(), gd.dzhi.data(),
                        gd.dxi, gd.dyi,
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                        gd.icells, gd.jcells, gd.ijcells);

                stats.calc_mask_stats(m, "u2_diff" , *u2_diff , no_offset, no_threshold);
                stats.calc_mask_stats(m, "v2_diff" , *v2_diff , no_offset, no_threshold);
                stats.calc_mask_stats(m, "w2_diff" , *w2_diff , no_offset, no_threshold);
                stats.calc_mask_stats(m, "tke_diff", *tke_diff, no_offset, no_threshold);
                stats.calc_mask_stats(m, "uw_diff" , *uw_diff , no_offset, no_threshold);
                stats.calc_mask_stats(m, "vw_diff" , *vw_diff , no_offset, no_threshold);

                fields.release_tmp(u2_diff);
                fields.release_tmp(v2_diff);
                fields.release_tmp(w2_diff);
                fields.release_tmp(tke_diff);
                fields.release_tmp(uw_diff);
                fields.release_tmp(vw_diff);
                fields.release_tmp(wz);
                fields.release_tmp(evisch);
            }
        }

        fields.release_tmp(wx);
        fields.release_tmp(wy);

        auto w2_pres = fields.get_tmp();
        auto tke_pres = fields.get_tmp();
        auto uw_pres = fields.get_tmp();
        auto vw_pres = fields.get_tmp();

        calc_pressure_transport_terms(
                w2_pres->fld.data(), tke_pres->fld.data(),
                uw_pres->fld.data(), vw_pres->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(),
                fields.mp.at("w")->fld.data(), fields.sd.at("p")->fld.data(),
                umodel.data(), vmodel.data(), wmodel.data(),
                gd.dzi.data(), gd.dzhi.data(), gd.dxi, gd.dyi,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_mask_stats(m, "w2_pres" , *w2_pres , no_offset, no_threshold);
        stats.calc_mask_stats(m, "tke_pres", *tke_pres, no_offset, no_threshold);
        stats.calc_mask_stats(m, "uw_pres" , *uw_pres , no_offset, no_threshold);
        stats.calc_mask_stats(m, "vw_pres" , *vw_pres , no_offset, no_threshold);

        auto u2_rdstr = fields.get_tmp();
        auto v2_rdstr = std::move(tke_pres);
        auto w2_rdstr = std::move(w2_pres);
        auto uw_rdstr = std::move(uw_pres);
        auto vw_rdstr = std::move(vw_pres);

        calc_pressure_redistribution_terms(
                u2_rdstr->fld.data(), v2_rdstr->fld.data(), w2_rdstr->fld.data(),
                uw_rdstr->fld.data(), vw_rdstr->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(),
                fields.mp.at("w")->fld.data(), fields.sd.at("p")->fld.data(),
                umodel.data(), vmodel.data(), wmodel.data(),
                gd.dzi.data(), gd.dzhi.data(), gd.dxi, gd.dyi,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_mask_stats(m, "u2_rdstr", *u2_rdstr , no_offset, no_threshold);
        stats.calc_mask_stats(m, "v2_rdstr", *v2_rdstr , no_offset, no_threshold);
        stats.calc_mask_stats(m, "w2_rdstr", *w2_rdstr , no_offset, no_threshold);
        stats.calc_mask_stats(m, "uw_rdstr", *uw_rdstr , no_offset, no_threshold);
        stats.calc_mask_stats(m, "vw_rdstr", *vw_rdstr , no_offset, no_threshold);

        fields.release_tmp(u2_rdstr);
        fields.release_tmp(v2_rdstr);
        fields.release_tmp(w2_rdstr);
        fields.release_tmp(uw_rdstr);
        fields.release_tmp(vw_rdstr);

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
                    umodel.data(), vmodel.data(), wmodel.data(), fc,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

            stats.calc_mask_stats(m, "u2_cor", *u2_cor, no_offset, no_threshold);
            stats.calc_mask_stats(m, "v2_cor", *v2_cor, no_offset, no_threshold);
            stats.calc_mask_stats(m, "uw_cor", *uw_cor, no_offset, no_threshold);
            stats.calc_mask_stats(m, "vw_cor", *vw_cor, no_offset, no_threshold);

            fields.release_tmp(u2_cor);
            fields.release_tmp(v2_cor);
            fields.release_tmp(uw_cor);
            fields.release_tmp(vw_cor);
        }

        if (thermo.get_switch() != Thermo_type::Disabled)
        {
            // Get the buoyancy diffusivity from the thermo class
            const TF diff_b = thermo.get_buoyancy_diffusivity();

            // Acquire the buoyancy, cyclic=true, is_stat=true.
            auto b = fields.get_tmp();
            thermo.get_thermo_field(*b, "b", true, true);

            // Calculate the mean of the fields.
            field3d_operators.calc_mean_profile(b->fld_mean.data(), b->fld.data());
            field3d_operators.calc_mean_profile(fields.sd.at("p")->fld_mean.data(), fields.sd.at("p")->fld.data());

            auto w2_buoy = fields.get_tmp();
            auto tke_buoy = fields.get_tmp();
            auto uw_buoy = fields.get_tmp();
            auto vw_buoy = fields.get_tmp();

            // Calculate buoyancy terms
            calc_buoyancy_terms(
                    w2_buoy->fld.data(), tke_buoy->fld.data(),
                    uw_buoy->fld.data(), vw_buoy->fld.data(),
                    fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(),
                    fields.mp.at("w")->fld.data(), b->fld.data(),
                    umodel.data(), vmodel.data(), wmodel.data(),
                    b->fld_mean.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

            stats.calc_mask_stats(m, "w2_buoy" , *w2_buoy , no_offset, no_threshold);
            stats.calc_mask_stats(m, "tke_buoy", *tke_buoy, no_offset, no_threshold);
            stats.calc_mask_stats(m, "uw_buoy" , *uw_buoy , no_offset, no_threshold);
            stats.calc_mask_stats(m, "vw_buoy" , *vw_buoy , no_offset, no_threshold);

            fields.release_tmp(w2_buoy);
            fields.release_tmp(tke_buoy);
            fields.release_tmp(uw_buoy);
            auto bw_buoy = std::move(vw_buoy);

            // Buoyancy variance and flux budgets
            calc_buoyancy_terms_scalar(
                    bw_buoy->fld.data(),
                    b->fld.data(), b->fld.data(),
                    b->fld_mean.data(), b->fld_mean.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

            stats.calc_mask_stats(m, "bw_buoy", *bw_buoy, no_offset, no_threshold);
            fields.release_tmp(bw_buoy);

            if (advec.get_switch() != Advection_type::Disabled)
            {
                auto b2_shear = fields.get_tmp();
                auto b2_turb = fields.get_tmp();
                auto bw_shear = fields.get_tmp();
                auto bw_turb = fields.get_tmp();

                calc_advection_terms_scalar(
                        b2_shear->fld.data(), b2_turb->fld.data(),
                        bw_shear->fld.data(), bw_turb->fld.data(),
                        fields.mp.at("w")->fld.data(), b->fld.data(),
                        b->fld_mean.data(),
                        gd.dzi.data(), gd.dzhi.data(),
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                        gd.icells, gd.ijcells);

                stats.calc_mask_stats(m, "b2_shear", *b2_shear, no_offset, no_threshold);
                stats.calc_mask_stats(m, "b2_turb" , *b2_turb , no_offset, no_threshold);
                stats.calc_mask_stats(m, "bw_shear", *bw_shear, no_offset, no_threshold);
                stats.calc_mask_stats(m, "bw_turb" , *bw_turb , no_offset, no_threshold);

                fields.release_tmp(b2_shear);
                fields.release_tmp(b2_turb );
                fields.release_tmp(bw_shear);
                fields.release_tmp(bw_turb );
            }

            if (diff.get_switch() == Diffusion_type::Diff_2 || diff.get_switch() == Diffusion_type::Diff_4)
            {
                auto b2_visc = fields.get_tmp();
                auto b2_diss = fields.get_tmp();
                auto bw_visc = fields.get_tmp();
                auto bw_diss = fields.get_tmp();

                calc_diffusion_terms_scalar_dns(
                        b2_visc->fld.data(), b2_diss->fld.data(),
                        bw_visc->fld.data(), bw_diss->fld.data(),
                        fields.mp.at("w")->fld.data(), b->fld.data(),
                        b->fld_mean.data(),
                        gd.dzi.data(), gd.dzhi.data(),
                        gd.dxi, gd.dyi, fields.visc, diff_b,
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                        gd.icells, gd.ijcells);

                stats.calc_mask_stats(m, "b2_visc", *b2_visc, no_offset, no_threshold);
                stats.calc_mask_stats(m, "b2_diss", *b2_diss, no_offset, no_threshold);
                stats.calc_mask_stats(m, "bw_visc", *bw_visc, no_offset, no_threshold);
                stats.calc_mask_stats(m, "bw_diss", *bw_diss, no_offset, no_threshold);

                fields.release_tmp(b2_visc);
                fields.release_tmp(b2_diss);
                fields.release_tmp(bw_visc);
                fields.release_tmp(bw_diss);
            }

            auto bw_pres = fields.get_tmp();
            auto bw_rdstr = fields.get_tmp();

            calc_pressure_terms_scalar(
                    bw_pres->fld.data(), bw_rdstr->fld.data(),
                    b->fld.data(), fields.sd.at("p")->fld.data(),
                    b->fld_mean.data(), fields.sd.at("p")->fld_mean.data(),
                    gd.dzi.data(), gd.dzhi.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

            stats.calc_mask_stats(m, "bw_pres", *bw_pres, no_offset, no_threshold);
            stats.calc_mask_stats(m, "bw_rdstr", *bw_rdstr, no_offset, no_threshold);

            fields.release_tmp(bw_pres);
            fields.release_tmp(bw_rdstr);
            fields.release_tmp(b);
        }
    }
}


#ifdef FLOAT_SINGLE
template class Budget_2<float>;
#else
template class Budget_2<double>;
#endif
