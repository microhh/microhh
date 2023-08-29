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

#include <cstdio>
#include <cmath>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "fast_math.h"
#include "finite_difference.h"
#include "model.h"
#include "thermo.h"
#include "diff.h"
#include "advec.h"
#include "force.h"
#include "stats.h"
#include "field3d_operators.h"
#include "constants.h"
#include "netcdf_interface.h"

#include "budget.h"
#include "budget_4.h"

namespace
{
    template<typename TF>
    void calc_ke(TF* restrict ke, TF* restrict tke,
                 const TF* restrict u, const TF* restrict v, const TF* restrict w,
                 const TF* restrict umodel, const TF* restrict vmodel, const TF* restrict wmodel,
                 const TF utrans, const TF vtrans,
                 const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
                 const int icells, const int ijcells)
    {
        using Fast_math::pow2;

        using namespace Finite_difference::O4;

        const int ii1 = 1;
        const int ii2 = 2;
        const int jj1 = 1*icells;
        const int jj2 = 2*icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        for (int k=kstart; k<kend; ++k)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    const TF u2 = ci0<TF>*pow2(u[ijk-ii1] + utrans) + ci1<TF>*pow2(u[ijk    ] + utrans)
                                + ci2<TF>*pow2(u[ijk+ii1] + utrans) + ci3<TF>*pow2(u[ijk+ii2] + utrans);
                    const TF v2 = ci0<TF>*pow2(v[ijk-jj1] + vtrans) + ci1<TF>*pow2(v[ijk    ] + vtrans)
                                + ci2<TF>*pow2(v[ijk+jj1] + vtrans) + ci3<TF>*pow2(v[ijk+jj2] + vtrans);
                    const TF w2 = ci0<TF>*pow2(w[ijk-kk1]) + ci1<TF>*pow2(w[ijk]) + ci2<TF>*pow2(w[ijk+kk1]) + ci3<TF>*pow2(w[ijk+kk2]);
                    ke[ijk] = TF(0.5)*(u2 + v2 + w2);
                }

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    const TF u2 = ci0<TF>*pow2(u[ijk-ii1] - umodel[k]) + ci1<TF>*pow2(u[ijk    ] - umodel[k])
                                + ci2<TF>*pow2(u[ijk+ii1] - umodel[k]) + ci3<TF>*pow2(u[ijk+ii2] - umodel[k]);
                    const TF v2 = ci0<TF>*pow2(v[ijk-jj1] - vmodel[k]) + ci1<TF>*pow2(v[ijk    ] - vmodel[k])
                                + ci2<TF>*pow2(v[ijk+jj1] - vmodel[k]) + ci3<TF>*pow2(v[ijk+jj2] - vmodel[k]);
                    const TF w2 = ci0<TF>*pow2(w[ijk-kk1] - wmodel[k-1]) + ci1<TF>*pow2(w[ijk] - wmodel[k]) + ci2<TF>*pow2(w[ijk+kk1] - wmodel[k+1]) + ci3<TF>*pow2(w[ijk+kk2] - wmodel[k+2]);
                    tke[ijk] = TF(0.5)*(u2 + v2 + w2);
                }
        }
    }

    template<typename TF>
    void calc_prime(
            TF* restrict a_prime, const TF* restrict a, const TF* restrict a_mean,
            const int icells, const int jcells, const int kcells,
            const int kk)
    {
        const int jj = icells;

        for (int k=0; k<kcells; ++k)
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    a_prime[ijk] = a[ijk] - a_mean[k];
                }
    }

    template<typename TF>
    void calc_tke_budget_shear(
            TF* restrict u2_shear, TF* restrict v2_shear, TF* restrict tke_shear, TF* restrict uw_shear,
            const TF* restrict u, const TF* restrict v, const TF* restrict w,
            const TF* restrict wx, const TF* restrict wy,
            const TF* restrict umean, const TF* restrict vmean,
            const TF* restrict dzi4, const TF* restrict dzhi4,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        using namespace Finite_difference::O4;

        const int jj1 = 1*icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        // CALCULATE THE SHEAR TERM u'w*dumean/dz

        // bottom boundary
        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                u2_shear[ijk] = -2.*(u[ijk]-umean[k])*(ci0<TF>*wx[ijk-kk1] + ci1<TF>*wx[ijk] + ci2<TF>*wx[ijk+kk1] + ci3<TF>*wx[ijk+kk2])
                              * ( cg0<TF>*(bi0<TF>*umean[k-2] + bi1<TF>*umean[k-1] + bi2<TF>*umean[k  ] + bi3<TF>*umean[k+1])
                                + cg1<TF>*(ci0<TF>*umean[k-2] + ci1<TF>*umean[k-1] + ci2<TF>*umean[k  ] + ci3<TF>*umean[k+1])
                                + cg2<TF>*(ci0<TF>*umean[k-1] + ci1<TF>*umean[k  ] + ci2<TF>*umean[k+1] + ci3<TF>*umean[k+2])
                                + cg3<TF>*(ci0<TF>*umean[k  ] + ci1<TF>*umean[k+1] + ci2<TF>*umean[k+2] + ci3<TF>*umean[k+3])) * dzi4[k];

                v2_shear[ijk] = -2.*(v[ijk]-vmean[k])*(ci0<TF>*wy[ijk-kk1] + ci1<TF>*wy[ijk] + ci2<TF>*wy[ijk+kk1] + ci3<TF>*wy[ijk+kk2])
                              * ( cg0<TF>*(bi0<TF>*vmean[k-2] + bi1<TF>*vmean[k-1] + bi2<TF>*vmean[k  ] + bi3<TF>*vmean[k+1])
                                + cg1<TF>*(ci0<TF>*vmean[k-2] + ci1<TF>*vmean[k-1] + ci2<TF>*vmean[k  ] + ci3<TF>*vmean[k+1])
                                + cg2<TF>*(ci0<TF>*vmean[k-1] + ci1<TF>*vmean[k  ] + ci2<TF>*vmean[k+1] + ci3<TF>*vmean[k+2])
                                + cg3<TF>*(ci0<TF>*vmean[k  ] + ci1<TF>*vmean[k+1] + ci2<TF>*vmean[k+2] + ci3<TF>*vmean[k+3])) * dzi4[k];

                tke_shear[ijk] = 0.5*(u2_shear[ijk] + v2_shear[ijk]);
            }

        // interior
        for (int k=kstart+1; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    u2_shear[ijk] = -2.*(u[ijk]-umean[k])*(ci0<TF>*wx[ijk-kk1] + ci1<TF>*wx[ijk] + ci2<TF>*wx[ijk+kk1] + ci3<TF>*wx[ijk+kk2])
                                  * ( cg0<TF>*(ci0<TF>*umean[k-3] + ci1<TF>*umean[k-2] + ci2<TF>*umean[k-1] + ci3<TF>*umean[k  ])
                                    + cg1<TF>*(ci0<TF>*umean[k-2] + ci1<TF>*umean[k-1] + ci2<TF>*umean[k  ] + ci3<TF>*umean[k+1])
                                    + cg2<TF>*(ci0<TF>*umean[k-1] + ci1<TF>*umean[k  ] + ci2<TF>*umean[k+1] + ci3<TF>*umean[k+2])
                                    + cg3<TF>*(ci0<TF>*umean[k  ] + ci1<TF>*umean[k+1] + ci2<TF>*umean[k+2] + ci3<TF>*umean[k+3])) * dzi4[k];

                    v2_shear[ijk] = -2.*(v[ijk]-vmean[k])*(ci0<TF>*wy[ijk-kk1] + ci1<TF>*wy[ijk] + ci2<TF>*wy[ijk+kk1] + ci3<TF>*wy[ijk+kk2])
                                  * ( cg0<TF>*(ci0<TF>*vmean[k-3] + ci1<TF>*vmean[k-2] + ci2<TF>*vmean[k-1] + ci3<TF>*vmean[k  ])
                                    + cg1<TF>*(ci0<TF>*vmean[k-2] + ci1<TF>*vmean[k-1] + ci2<TF>*vmean[k  ] + ci3<TF>*vmean[k+1])
                                    + cg2<TF>*(ci0<TF>*vmean[k-1] + ci1<TF>*vmean[k  ] + ci2<TF>*vmean[k+1] + ci3<TF>*vmean[k+2])
                                    + cg3<TF>*(ci0<TF>*vmean[k  ] + ci1<TF>*vmean[k+1] + ci2<TF>*vmean[k+2] + ci3<TF>*vmean[k+3])) * dzi4[k];

                    tke_shear[ijk] = 0.5*(u2_shear[ijk] + v2_shear[ijk]);
                }

        // top boundary
        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                u2_shear[ijk] = -2.*(u[ijk]-umean[k])*(ci0<TF>*wx[ijk-kk1] + ci1<TF>*wx[ijk] + ci2<TF>*wx[ijk+kk1] + ci3<TF>*wx[ijk+kk2])
                              * ( cg0<TF>*(ci0<TF>*umean[k-3] + ci1<TF>*umean[k-2] + ci2<TF>*umean[k-1] + ci3<TF>*umean[k  ])
                                + cg1<TF>*(ci0<TF>*umean[k-2] + ci1<TF>*umean[k-1] + ci2<TF>*umean[k  ] + ci3<TF>*umean[k+1])
                                + cg2<TF>*(ci0<TF>*umean[k-1] + ci1<TF>*umean[k  ] + ci2<TF>*umean[k+1] + ci3<TF>*umean[k+2])
                                + cg3<TF>*(ti0<TF>*umean[k  ] + ti1<TF>*umean[k+1] + ti2<TF>*umean[k+2] + ti3<TF>*umean[k+3])) * dzi4[k];

                v2_shear[ijk] = -2.*(v[ijk]-vmean[k])*(ci0<TF>*wy[ijk-kk1] + ci1<TF>*wy[ijk] + ci2<TF>*wy[ijk+kk1] + ci3<TF>*wy[ijk+kk2])
                              * ( cg0<TF>*(ci0<TF>*vmean[k-3] + ci1<TF>*vmean[k-2] + ci2<TF>*vmean[k-1] + ci3<TF>*vmean[k  ])
                                + cg1<TF>*(ci0<TF>*vmean[k-2] + ci1<TF>*vmean[k-1] + ci2<TF>*vmean[k  ] + ci3<TF>*vmean[k+1])
                                + cg2<TF>*(ci0<TF>*vmean[k-1] + ci1<TF>*vmean[k  ] + ci2<TF>*vmean[k+1] + ci3<TF>*vmean[k+2])
                                + cg3<TF>*(ti0<TF>*vmean[k-1] + ti1<TF>*vmean[k  ] + ti2<TF>*vmean[k+1] + ti3<TF>*vmean[k+2])) * dzi4[k];

                tke_shear[ijk] = 0.5*(u2_shear[ijk] + v2_shear[ijk]);
            }

        // Reynolds stresses
        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    uw_shear[ijk] = -( std::pow(wx[ijk],2)
                                  * ( cg0<TF>*umean[k-2] + cg1<TF>*umean[k-1] + cg2<TF>*umean[k] + cg3<TF>*umean[k+1] ) ) * dzhi4[k];
                }
    }

    template<typename TF>
    void calc_tke_budget_turb(
            TF* restrict u2_turb, TF* restrict v2_turb, TF* restrict w2_turb, TF* restrict tke_turb, TF* restrict uw_turb,
            const TF* restrict u, const TF* restrict v, const TF* restrict w,
            const TF* restrict wx, const TF* restrict wy,
            const TF* restrict umean, const TF* restrict vmean,
            const TF* restrict dzi4, const TF* restrict dzhi4,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        using namespace Finite_difference::O4;

        const int jj1 = 1*icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;
        const int kk3 = 3*ijcells;

        // CALCULATE TURBULENT FLUXES
        // bottom boundary
        int k = kstart;

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                u2_turb[ijk] = - ( cg0<TF>*((bi0<TF>*std::pow(u[ijk-kk2]-umean[k-2],2) + bi1<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + bi2<TF>*std::pow(u[ijk    ]-umean[k  ],2) + bi3<TF>*std::pow(u[ijk+kk1]-umean[k+1],2))*wx[ijk-kk1])
                                 + cg1<TF>*((ci0<TF>*std::pow(u[ijk-kk2]-umean[k-2],2) + ci1<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + ci2<TF>*std::pow(u[ijk    ]-umean[k  ],2) + ci3<TF>*std::pow(u[ijk+kk1]-umean[k+1],2))*wx[ijk    ])
                                 + cg2<TF>*((ci0<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + ci1<TF>*std::pow(u[ijk    ]-umean[k  ],2) + ci2<TF>*std::pow(u[ijk+kk1]-umean[k+1],2) + ci3<TF>*std::pow(u[ijk+kk2]-umean[k+2],2))*wx[ijk+kk1])
                                 + cg3<TF>*((ci0<TF>*std::pow(u[ijk    ]-umean[k  ],2) + ci1<TF>*std::pow(u[ijk+kk1]-umean[k+1],2) + ci2<TF>*std::pow(u[ijk+kk2]-umean[k+2],2) + ci3<TF>*std::pow(u[ijk+kk3]-umean[k+3],2))*wx[ijk+kk2]) ) * dzi4[k];

                v2_turb[ijk] = - ( cg0<TF>*((bi0<TF>*std::pow(v[ijk-kk2]-vmean[k-2],2) + bi1<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + bi2<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + bi3<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2))*wy[ijk-kk1])
                                 + cg1<TF>*((ci0<TF>*std::pow(v[ijk-kk2]-vmean[k-2],2) + ci1<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci2<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + ci3<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2))*wy[ijk    ])
                                 + cg2<TF>*((ci0<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci1<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + ci2<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2) + ci3<TF>*std::pow(v[ijk+kk2]-vmean[k+2],2))*wy[ijk+kk1])
                                 + cg3<TF>*((ci0<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + ci1<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2) + ci2<TF>*std::pow(v[ijk+kk2]-vmean[k+2],2) + ci3<TF>*std::pow(v[ijk+kk3]-vmean[k+3],2))*wy[ijk+kk2]) ) * dzi4[k];

                tke_turb[ijk] = -0.5*( cg0<TF>*std::pow(w[ijk-kk1], 3) + cg1<TF>*std::pow(w[ijk], 3) + cg2<TF>*std::pow(w[ijk+kk1], 3) + cg3<TF>*std::pow(w[ijk+kk2], 3)) * dzi4[k];

                tke_turb[ijk] += 0.5*(u2_turb[ijk] + v2_turb[ijk]);
            }

        // interior
        for (int k=kstart+1; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    u2_turb[ijk] = - ( cg0<TF>*((ci0<TF>*std::pow(u[ijk-kk3]-umean[k-3],2) + ci1<TF>*std::pow(u[ijk-kk2]-umean[k-2],2) + ci2<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + ci3<TF>*std::pow(u[ijk    ]-umean[k  ],2))*wx[ijk-kk1])
                                     + cg1<TF>*((ci0<TF>*std::pow(u[ijk-kk2]-umean[k-2],2) + ci1<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + ci2<TF>*std::pow(u[ijk    ]-umean[k  ],2) + ci3<TF>*std::pow(u[ijk+kk1]-umean[k+1],2))*wx[ijk    ])
                                     + cg2<TF>*((ci0<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + ci1<TF>*std::pow(u[ijk    ]-umean[k  ],2) + ci2<TF>*std::pow(u[ijk+kk1]-umean[k+1],2) + ci3<TF>*std::pow(u[ijk+kk2]-umean[k+2],2))*wx[ijk+kk1])
                                     + cg3<TF>*((ci0<TF>*std::pow(u[ijk    ]-umean[k  ],2) + ci1<TF>*std::pow(u[ijk+kk1]-umean[k+1],2) + ci2<TF>*std::pow(u[ijk+kk2]-umean[k+2],2) + ci3<TF>*std::pow(u[ijk+kk3]-umean[k+3],2))*wx[ijk+kk2]) ) * dzi4[k];

                    v2_turb[ijk] = - ( cg0<TF>*((ci0<TF>*std::pow(v[ijk-kk3]-vmean[k-3],2) + ci1<TF>*std::pow(v[ijk-kk2]-vmean[k-2],2) + ci2<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci3<TF>*std::pow(v[ijk    ]-vmean[k  ],2))*wy[ijk-kk1])
                                     + cg1<TF>*((ci0<TF>*std::pow(v[ijk-kk2]-vmean[k-2],2) + ci1<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci2<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + ci3<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2))*wy[ijk    ])
                                     + cg2<TF>*((ci0<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci1<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + ci2<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2) + ci3<TF>*std::pow(v[ijk+kk2]-vmean[k+2],2))*wy[ijk+kk1])
                                     + cg3<TF>*((ci0<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + ci1<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2) + ci2<TF>*std::pow(v[ijk+kk2]-vmean[k+2],2) + ci3<TF>*std::pow(v[ijk+kk3]-vmean[k+3],2))*wy[ijk+kk2]) ) * dzi4[k];

                    tke_turb[ijk] = -0.5*( cg0<TF>*std::pow(w[ijk-kk1], 3) + cg1<TF>*std::pow(w[ijk], 3) + cg2<TF>*std::pow(w[ijk+kk1], 3) + cg3<TF>*std::pow(w[ijk+kk2], 3)) * dzi4[k];

                    tke_turb[ijk] += 0.5*(u2_turb[ijk] + v2_turb[ijk]);
                }

        // top boundary
        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                u2_turb[ijk] = - ( cg0<TF>*((ci0<TF>*std::pow(u[ijk-kk3]-umean[k-3],2) + ci1<TF>*std::pow(u[ijk-kk2]-umean[k-2],2) + ci2<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + ci3<TF>*std::pow(u[ijk    ]-umean[k  ],2))*wx[ijk-kk1])
                                 + cg1<TF>*((ci0<TF>*std::pow(u[ijk-kk2]-umean[k-2],2) + ci1<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + ci2<TF>*std::pow(u[ijk    ]-umean[k  ],2) + ci3<TF>*std::pow(u[ijk+kk1]-umean[k+1],2))*wx[ijk    ])
                                 + cg2<TF>*((ci0<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + ci1<TF>*std::pow(u[ijk    ]-umean[k  ],2) + ci2<TF>*std::pow(u[ijk+kk1]-umean[k+1],2) + ci3<TF>*std::pow(u[ijk+kk2]-umean[k+2],2))*wx[ijk+kk1])
                                 + cg3<TF>*((ti0<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + ti1<TF>*std::pow(u[ijk    ]-umean[k  ],2) + ti2<TF>*std::pow(u[ijk+kk1]-umean[k+1],2) + ti3<TF>*std::pow(u[ijk+kk2]-umean[k+2],2))*wx[ijk+kk1]) ) * dzi4[k];

                v2_turb[ijk] = - ( cg0<TF>*((ci0<TF>*std::pow(v[ijk-kk3]-vmean[k-3],2) + ci1<TF>*std::pow(v[ijk-kk2]-vmean[k-2],2) + ci2<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci3<TF>*std::pow(v[ijk    ]-vmean[k  ],2))*wy[ijk-kk1])
                                 + cg1<TF>*((ci0<TF>*std::pow(v[ijk-kk2]-vmean[k-2],2) + ci1<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci2<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + ci3<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2))*wy[ijk    ])
                                 + cg2<TF>*((ci0<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + ci1<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + ci2<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2) + ci3<TF>*std::pow(v[ijk+kk2]-vmean[k+2],2))*wy[ijk+kk1])
                                 + cg3<TF>*((ti0<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + ti1<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + ti2<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2) + ti3<TF>*std::pow(v[ijk+kk2]-vmean[k+2],2))*wy[ijk+kk1]) ) * dzi4[k];

                tke_turb[ijk] = -0.5*( cg0<TF>*std::pow(w[ijk-kk1], 3) + cg1<TF>*std::pow(w[ijk], 3) + cg2<TF>*std::pow(w[ijk+kk1], 3) + cg3<TF>*std::pow(w[ijk+kk2], 3)) * dzi4[k];

                tke_turb[ijk] += 0.5*(u2_turb[ijk] + v2_turb[ijk]);
            }

        // calculate the vertical velocity term and the vertical reynold stresses
        k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                w2_turb[ijk] = - ( cg0<TF>*(bi0<TF>*std::pow(w[ijk-kk2],3) + bi1<TF>*std::pow(w[ijk-kk1],3) + bi2<TF>*std::pow(w[ijk    ],3) + bi3<TF>*std::pow(w[ijk+kk1],3))
                                 + cg1<TF>*(ci0<TF>*std::pow(w[ijk-kk2],3) + ci1<TF>*std::pow(w[ijk-kk1],3) + ci2<TF>*std::pow(w[ijk    ],3) + ci3<TF>*std::pow(w[ijk+kk1],3))
                                 + cg2<TF>*(ci0<TF>*std::pow(w[ijk-kk1],3) + ci1<TF>*std::pow(w[ijk    ],3) + ci2<TF>*std::pow(w[ijk+kk1],3) + ci3<TF>*std::pow(w[ijk+kk2],3))
                                 + cg3<TF>*(ci0<TF>*std::pow(w[ijk    ],3) + ci1<TF>*std::pow(w[ijk+kk1],3) + ci2<TF>*std::pow(w[ijk+kk2],3) + ci3<TF>*std::pow(w[ijk+kk3],3)) ) * dzhi4[k];

                uw_turb[ijk] = - ( ( cg0<TF>*( std::pow(bi0<TF>*wx[ijk-kk2] + bi1<TF>*wx[ijk-kk1] + bi2<TF>*wx[ijk    ] + bi3<TF>*wx[ijk+kk1], 2) * (u[ijk-kk2]-umean[k-2]) )
                                   + cg1<TF>*( std::pow(ci0<TF>*wx[ijk-kk2] + ci1<TF>*wx[ijk-kk1] + ci2<TF>*wx[ijk    ] + ci3<TF>*wx[ijk+kk1], 2) * (u[ijk-kk1]-umean[k-1]) )
                                   + cg2<TF>*( std::pow(ci0<TF>*wx[ijk-kk1] + ci1<TF>*wx[ijk    ] + ci2<TF>*wx[ijk+kk1] + ci3<TF>*wx[ijk+kk2], 2) * (u[ijk    ]-umean[k  ]) )
                                   + cg3<TF>*( std::pow(ci0<TF>*wx[ijk    ] + ci1<TF>*wx[ijk+kk1] + ci2<TF>*wx[ijk+kk2] + ci3<TF>*wx[ijk+kk3], 2) * (u[ijk+kk1]-umean[k+1]) ) )
                                 * dzhi4[k  ] );
            }

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    w2_turb[ijk] = - ( cg0<TF>*(ci0<TF>*std::pow(w[ijk-kk3],3) + ci1<TF>*std::pow(w[ijk-kk2],3) + ci2<TF>*std::pow(w[ijk-kk1],3) + ci3<TF>*std::pow(w[ijk    ],3))
                                     + cg1<TF>*(ci0<TF>*std::pow(w[ijk-kk2],3) + ci1<TF>*std::pow(w[ijk-kk1],3) + ci2<TF>*std::pow(w[ijk    ],3) + ci3<TF>*std::pow(w[ijk+kk1],3))
                                     + cg2<TF>*(ci0<TF>*std::pow(w[ijk-kk1],3) + ci1<TF>*std::pow(w[ijk    ],3) + ci2<TF>*std::pow(w[ijk+kk1],3) + ci3<TF>*std::pow(w[ijk+kk2],3))
                                     + cg3<TF>*(ci0<TF>*std::pow(w[ijk    ],3) + ci1<TF>*std::pow(w[ijk+kk1],3) + ci2<TF>*std::pow(w[ijk+kk2],3) + ci3<TF>*std::pow(w[ijk+kk3],3)) ) * dzhi4[k];

                    uw_turb[ijk] = - ( ( cg0<TF>*( std::pow(ci0<TF>*wx[ijk-kk3] + ci1<TF>*wx[ijk-kk2] + ci2<TF>*wx[ijk-kk1] + ci3<TF>*wx[ijk    ], 2) * (u[ijk-kk2]-umean[k-2]) )
                                       + cg1<TF>*( std::pow(ci0<TF>*wx[ijk-kk2] + ci1<TF>*wx[ijk-kk1] + ci2<TF>*wx[ijk    ] + ci3<TF>*wx[ijk+kk1], 2) * (u[ijk-kk1]-umean[k-1]) )
                                       + cg2<TF>*( std::pow(ci0<TF>*wx[ijk-kk1] + ci1<TF>*wx[ijk    ] + ci2<TF>*wx[ijk+kk1] + ci3<TF>*wx[ijk+kk2], 2) * (u[ijk    ]-umean[k  ]) )
                                       + cg3<TF>*( std::pow(ci0<TF>*wx[ijk    ] + ci1<TF>*wx[ijk+kk1] + ci2<TF>*wx[ijk+kk2] + ci3<TF>*wx[ijk+kk3], 2) * (u[ijk+kk1]-umean[k+1]) ) )
                                     * dzhi4[k  ] );
                }

        k = kend;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                w2_turb[ijk] = - ( cg0<TF>*(ci0<TF>*std::pow(w[ijk-kk3],3) + ci1<TF>*std::pow(w[ijk-kk2],3) + ci2<TF>*std::pow(w[ijk-kk1],3) + ci3<TF>*std::pow(w[ijk    ],3))
                                 + cg1<TF>*(ci0<TF>*std::pow(w[ijk-kk2],3) + ci1<TF>*std::pow(w[ijk-kk1],3) + ci2<TF>*std::pow(w[ijk    ],3) + ci3<TF>*std::pow(w[ijk+kk1],3))
                                 + cg2<TF>*(ci0<TF>*std::pow(w[ijk-kk1],3) + ci1<TF>*std::pow(w[ijk    ],3) + ci2<TF>*std::pow(w[ijk+kk1],3) + ci3<TF>*std::pow(w[ijk+kk2],3))
                                 + cg3<TF>*(ti0<TF>*std::pow(w[ijk-kk1],3) + ti1<TF>*std::pow(w[ijk    ],3) + ti2<TF>*std::pow(w[ijk+kk1],3) + ti3<TF>*std::pow(w[ijk+kk2],3)) ) * dzhi4[k];

                uw_turb[ijk] = - ( ( cg0<TF>*( ( ci0<TF>*wx[ijk-kk3] + ci1<TF>*wx[ijk-kk2] + ci2<TF>*wx[ijk-kk1] + ci3<TF>*wx[ijk    ] ) * ( u[ijk-kk2] - umean[k-2] ) )
                                   + cg1<TF>*( ( ci0<TF>*wx[ijk-kk2] + ci1<TF>*wx[ijk-kk1] + ci2<TF>*wx[ijk    ] + ci3<TF>*wx[ijk+kk1] ) * ( u[ijk-kk1] - umean[k-1] ) )
                                   + cg2<TF>*( ( ci0<TF>*wx[ijk-kk1] + ci1<TF>*wx[ijk    ] + ci2<TF>*wx[ijk+kk1] + ci3<TF>*wx[ijk+kk2] ) * ( u[ijk    ] - umean[k  ] ) )
                                   + cg3<TF>*( ( ti0<TF>*wx[ijk-kk1] + ti1<TF>*wx[ijk    ] + ti2<TF>*wx[ijk+kk1] + ti3<TF>*wx[ijk+kk2] ) * ( u[ijk+kk1] - umean[k+1] ) ) )
                                 * dzhi4[k  ] );
            }
    }

    template<typename TF>
    void calc_tke_budget_pres(
            TF* restrict w2_pres, TF* restrict tke_pres, TF* restrict uw_pres,
            const TF* restrict u, const TF* restrict v, const TF* restrict w, const TF* restrict p,
            const TF* restrict umean, const TF* restrict vmean,
            const TF* restrict dzi4, const TF* restrict dzhi4,
            const TF dxi, const TF dyi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        using namespace Finite_difference::O4;

        const int ii1 = 1;
        const int ii2 = 2;
        const int jj1 = 1*icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;
        const int kk3 = 3*ijcells;

        // CALCULATE THE PRESSURE TRANSPORT TERM
        // bottom boundary
        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                tke_pres[ijk] = - ( cg0<TF>*((bi0<TF>*p[ijk-kk2] + bi1<TF>*p[ijk-kk1] + bi2<TF>*p[ijk    ] + bi3<TF>*p[ijk+kk1])*w[ijk-kk1])
                                  + cg1<TF>*((ci0<TF>*p[ijk-kk2] + ci1<TF>*p[ijk-kk1] + ci2<TF>*p[ijk    ] + ci3<TF>*p[ijk+kk1])*w[ijk    ])
                                  + cg2<TF>*((ci0<TF>*p[ijk-kk1] + ci1<TF>*p[ijk    ] + ci2<TF>*p[ijk+kk1] + ci3<TF>*p[ijk+kk2])*w[ijk+kk1])
                                  + cg3<TF>*((ci0<TF>*p[ijk    ] + ci1<TF>*p[ijk+kk1] + ci2<TF>*p[ijk+kk2] + ci3<TF>*p[ijk+kk3])*w[ijk+kk2]) ) * dzi4[k];
            }

        // interior
        for (int k=kstart+1; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
    #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    tke_pres[ijk] = - ( cg0<TF>*((ci0<TF>*p[ijk-kk3] + ci1<TF>*p[ijk-kk2] + ci2<TF>*p[ijk-kk1] + ci3<TF>*p[ijk    ])*w[ijk-kk1])
                                      + cg1<TF>*((ci0<TF>*p[ijk-kk2] + ci1<TF>*p[ijk-kk1] + ci2<TF>*p[ijk    ] + ci3<TF>*p[ijk+kk1])*w[ijk    ])
                                      + cg2<TF>*((ci0<TF>*p[ijk-kk1] + ci1<TF>*p[ijk    ] + ci2<TF>*p[ijk+kk1] + ci3<TF>*p[ijk+kk2])*w[ijk+kk1])
                                      + cg3<TF>*((ci0<TF>*p[ijk    ] + ci1<TF>*p[ijk+kk1] + ci2<TF>*p[ijk+kk2] + ci3<TF>*p[ijk+kk3])*w[ijk+kk2]) ) * dzi4[k];
                }

        // top boundary
        k = kend-1;
        for (int j=jstart; j<jend; ++j)
    #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                tke_pres[ijk] = - ( cg0<TF>*((ci0<TF>*p[ijk-kk3] + ci1<TF>*p[ijk-kk2] + ci2<TF>*p[ijk-kk1] + ci3<TF>*p[ijk    ])*w[ijk-kk1])
                                  + cg1<TF>*((ci0<TF>*p[ijk-kk2] + ci1<TF>*p[ijk-kk1] + ci2<TF>*p[ijk    ] + ci3<TF>*p[ijk+kk1])*w[ijk    ])
                                  + cg2<TF>*((ci0<TF>*p[ijk-kk1] + ci1<TF>*p[ijk    ] + ci2<TF>*p[ijk+kk1] + ci3<TF>*p[ijk+kk2])*w[ijk+kk1])
                                  + cg3<TF>*((ti0<TF>*p[ijk-kk1] + ti1<TF>*p[ijk    ] + ti2<TF>*p[ijk+kk1] + ti3<TF>*p[ijk+kk2])*w[ijk+kk2]) ) * dzi4[k];
            }

        // calculate the vertical velocity pressure transport term
        // CvH: implement the proper BC as soon as the full BC's for pressure are added
        // bottom boundary
        k = kstart;
        for (int j=jstart; j<jend; ++j)
    #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                w2_pres[ijk] = - 0.* ( cg0<TF>*((bi0<TF>*w[ijk-kk2] + bi1<TF>*w[ijk-kk1] + bi2<TF>*w[ijk    ] + bi3<TF>*w[ijk+kk1])*p[ijk-kk2])
                                     + cg1<TF>*((ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1])*p[ijk-kk1])
                                     + cg2<TF>*((ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2])*p[ijk    ])
                                     + cg3<TF>*((ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3])*p[ijk+kk1]) ) * dzhi4[k];
            }

        // interior
        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
    #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    w2_pres[ijk] = - 2.*( cg0<TF>*((ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ])*p[ijk-kk2])
                                        + cg1<TF>*((ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1])*p[ijk-kk1])
                                        + cg2<TF>*((ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2])*p[ijk    ])
                                        + cg3<TF>*((ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3])*p[ijk+kk1]) ) * dzhi4[k];
                }

        // top boundary
        k = kend;
        for (int j=jstart; j<jend; ++j)
    #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                w2_pres[ijk] = - 0.*( cg0<TF>*((ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ])*p[ijk-kk2])
                                    + cg1<TF>*((ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1])*p[ijk-kk1])
                                    + cg2<TF>*((ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2])*p[ijk    ])
                                    + cg3<TF>*((ti0<TF>*w[ijk-kk1] + ti1<TF>*w[ijk    ] + ti2<TF>*w[ijk+kk1] + ti3<TF>*w[ijk+kk2])*p[ijk+kk1]) ) * dzhi4[k];
            }

        // UW term
        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
    #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    uw_pres[ijk] = - ( ( ( cg0<TF>*( ( u[ijk        -kk2] - umean[k-2] ) * ( ci0<TF>*p[ijk-ii2    -kk2] + ci1<TF>*p[ijk-ii1    -kk2] + ci2<TF>*p[ijk        -kk2] + ci3<TF>*p[ijk+ii1    -kk2] ) )
                                         + cg1<TF>*( ( u[ijk        -kk1] - umean[k-1] ) * ( ci0<TF>*p[ijk-ii2    -kk1] + ci1<TF>*p[ijk-ii1    -kk1] + ci2<TF>*p[ijk        -kk1] + ci3<TF>*p[ijk+ii1    -kk1] ) )
                                         + cg2<TF>*( ( u[ijk            ] - umean[k  ] ) * ( ci0<TF>*p[ijk-ii2        ] + ci1<TF>*p[ijk-ii1        ] + ci2<TF>*p[ijk            ] + ci3<TF>*p[ijk+ii1        ] ) )
                                         + cg3<TF>*( ( u[ijk        +kk1] - umean[k+1] ) * ( ci0<TF>*p[ijk-ii2    +kk1] + ci1<TF>*p[ijk-ii1    +kk1] + ci2<TF>*p[ijk        +kk1] + ci3<TF>*p[ijk+ii1    +kk1] ) ) )

                                       * dzhi4[k  ] )

                                     + ( ( cg0<TF>*( w[ijk-ii2        ] * ( ci0<TF>*p[ijk-ii2    -kk2] + ci1<TF>*p[ijk-ii2    -kk1] + ci2<TF>*p[ijk-ii2        ] + ci3<TF>*p[ijk-ii2    +kk1] ) )
                                         + cg1<TF>*( w[ijk-ii1        ] * ( ci0<TF>*p[ijk-ii1    -kk2] + ci1<TF>*p[ijk-ii1    -kk1] + ci2<TF>*p[ijk-ii1        ] + ci3<TF>*p[ijk-ii1    +kk1] ) )
                                         + cg2<TF>*( w[ijk            ] * ( ci0<TF>*p[ijk        -kk2] + ci1<TF>*p[ijk        -kk1] + ci2<TF>*p[ijk            ] + ci3<TF>*p[ijk        +kk1] ) )
                                         + cg3<TF>*( w[ijk+ii1        ] * ( ci0<TF>*p[ijk+ii1    -kk2] + ci1<TF>*p[ijk+ii1    -kk1] + ci2<TF>*p[ijk+ii1        ] + ci3<TF>*p[ijk+ii1    +kk1] ) ) )

                                       * dxi ) );
                }
    }

    template<typename TF>
    void calc_tke_budget_visc(
            TF* restrict u2_visc, TF* restrict v2_visc, TF* restrict w2_visc, TF* restrict tke_visc, TF* restrict uw_visc,
            TF* restrict wz, TF* restrict uz,
            const TF* restrict u, const TF* restrict v, const TF* restrict w,
            const TF* restrict umean, const TF* restrict vmean,
            const TF* restrict dzi4, const TF* restrict dzhi4,
            const TF dxi, const TF dyi, const TF dzhi4bot, const TF dzhi4top,
            const TF visc,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        using namespace Finite_difference::O4;

        const int ii1 = 1;
        const int ii2 = 2;
        const int jj1 = 1*icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;
        const int kk3 = 3*ijcells;
        const int kk4 = 4*ijcells;

        // CALCULATE THE VISCOUS TRANSPORT TERM
        // first, interpolate the vertical velocity to the scalar levels using temporary array wz
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    wz[ijk] = ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2];
                }

        // calculate the ghost cells at the bottom
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + kstart*kk1;
                wz[ijk-kk1] = - 2.*wz[ijk] + (1./3.)*wz[ijk+kk1];
                wz[ijk-kk2] = - 9.*wz[ijk] + 2.*wz[ijk+kk1];
            }

        // calculate the ghost cells at the top
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + (kend-1)*kk1;
                wz[ijk+kk1] = - 2.*wz[ijk] + (1./3.)*wz[ijk-kk1];
                wz[ijk+kk2] = - 9.*wz[ijk] + 2.*wz[ijk-kk1];
            }

        // first, interpolate the horizontal velocity to the flux levels using temporary array uz
        int k = kstart-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                uz[ijk] = bi0<TF>*u[ijk-kk1] + bi1<TF>*u[ijk] + bi2<TF>*u[ijk+kk1] + bi3<TF>*u[ijk+kk2];
            }

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    uz[ijk] = ci0<TF>*u[ijk-kk2] + ci1<TF>*u[ijk-kk1] + ci2<TF>*u[ijk] + ci3<TF>*u[ijk+kk1];
                }

        k = kend;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                uz[ijk] = ti0<TF>*u[ijk-kk2] + ti1<TF>*u[ijk-kk1] + ti2<TF>*u[ijk] + ti3<TF>*u[ijk+kk1];
            }


        // bottom boundary
        k = kstart;
        for (int j=jstart; j<jend; ++j)
    #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk  = i + j*jj1 + k*kk1;
                u2_visc[ijk] = visc * ( cg0<TF>*((bg0<TF>*std::pow(u[ijk-kk2]-umean[k-2],2) + bg1<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + bg2<TF>*std::pow(u[ijk    ]-umean[k  ],2) + bg3<TF>*std::pow(u[ijk+kk1]-umean[k+1],2)) * dzhi4[k-1])
                                      + cg1<TF>*((cg0<TF>*std::pow(u[ijk-kk2]-umean[k-2],2) + cg1<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + cg2<TF>*std::pow(u[ijk    ]-umean[k  ],2) + cg3<TF>*std::pow(u[ijk+kk1]-umean[k+1],2)) * dzhi4[k  ])
                                      + cg2<TF>*((cg0<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + cg1<TF>*std::pow(u[ijk    ]-umean[k  ],2) + cg2<TF>*std::pow(u[ijk+kk1]-umean[k+1],2) + cg3<TF>*std::pow(u[ijk+kk2]-umean[k+2],2)) * dzhi4[k+1])
                                      + cg3<TF>*((cg0<TF>*std::pow(u[ijk    ]-umean[k  ],2) + cg1<TF>*std::pow(u[ijk+kk1]-umean[k+1],2) + cg2<TF>*std::pow(u[ijk+kk2]-umean[k+2],2) + cg3<TF>*std::pow(u[ijk+kk3]-umean[k+3],2)) * dzhi4[k+2]) ) * dzi4[k];

                v2_visc[ijk] = visc * ( cg0<TF>*((bg0<TF>*std::pow(v[ijk-kk2]-vmean[k-2],2) + bg1<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + bg2<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + bg3<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2)) * dzhi4[k-1])
                                      + cg1<TF>*((cg0<TF>*std::pow(v[ijk-kk2]-vmean[k-2],2) + cg1<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg2<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + cg3<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2)) * dzhi4[k  ])
                                      + cg2<TF>*((cg0<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg1<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + cg2<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2) + cg3<TF>*std::pow(v[ijk+kk2]-vmean[k+2],2)) * dzhi4[k+1])
                                      + cg3<TF>*((cg0<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + cg1<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2) + cg2<TF>*std::pow(v[ijk+kk2]-vmean[k+2],2) + cg3<TF>*std::pow(v[ijk+kk3]-vmean[k+3],2)) * dzhi4[k+2]) ) * dzi4[k];

                tke_visc[ijk] = 0.5 * visc * ( cg0<TF>*((bg0<TF>*std::pow(wz[ijk-kk2],2) + bg1<TF>*std::pow(wz[ijk-kk1],2) + bg2<TF>*std::pow(wz[ijk    ],2) + bg3<TF>*std::pow(wz[ijk+kk1],2)) * dzhi4[k-1])
                                             + cg1<TF>*((cg0<TF>*std::pow(wz[ijk-kk2],2) + cg1<TF>*std::pow(wz[ijk-kk1],2) + cg2<TF>*std::pow(wz[ijk    ],2) + cg3<TF>*std::pow(wz[ijk+kk1],2)) * dzhi4[k  ])
                                             + cg2<TF>*((cg0<TF>*std::pow(wz[ijk-kk1],2) + cg1<TF>*std::pow(wz[ijk    ],2) + cg2<TF>*std::pow(wz[ijk+kk1],2) + cg3<TF>*std::pow(wz[ijk+kk2],2)) * dzhi4[k+1])
                                             + cg3<TF>*((cg0<TF>*std::pow(wz[ijk    ],2) + cg1<TF>*std::pow(wz[ijk+kk1],2) + cg2<TF>*std::pow(wz[ijk+kk2],2) + cg3<TF>*std::pow(wz[ijk+kk3],2)) * dzhi4[k+2]) ) * dzi4[k];

                tke_visc[ijk] += 0.5*(u2_visc[ijk] + v2_visc[ijk]);
            }

        // interior
        for (int k=kstart+1; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
    #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    u2_visc[ijk] = visc * ( cg0<TF>*((cg0<TF>*std::pow(u[ijk-kk3]-umean[k-3],2) + cg1<TF>*std::pow(u[ijk-kk2]-umean[k-2],2) + cg2<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + cg3<TF>*std::pow(u[ijk    ]-umean[k  ],2)) * dzhi4[k-1])
                                          + cg1<TF>*((cg0<TF>*std::pow(u[ijk-kk2]-umean[k-2],2) + cg1<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + cg2<TF>*std::pow(u[ijk    ]-umean[k  ],2) + cg3<TF>*std::pow(u[ijk+kk1]-umean[k+1],2)) * dzhi4[k  ])
                                          + cg2<TF>*((cg0<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + cg1<TF>*std::pow(u[ijk    ]-umean[k  ],2) + cg2<TF>*std::pow(u[ijk+kk1]-umean[k+1],2) + cg3<TF>*std::pow(u[ijk+kk2]-umean[k+2],2)) * dzhi4[k+1])
                                          + cg3<TF>*((cg0<TF>*std::pow(u[ijk    ]-umean[k  ],2) + cg1<TF>*std::pow(u[ijk+kk1]-umean[k+1],2) + cg2<TF>*std::pow(u[ijk+kk2]-umean[k+2],2) + cg3<TF>*std::pow(u[ijk+kk3]-umean[k+3],2)) * dzhi4[k+2]) ) * dzi4[k];

                    v2_visc[ijk] = visc * ( cg0<TF>*((cg0<TF>*std::pow(v[ijk-kk3]-vmean[k-3],2) + cg1<TF>*std::pow(v[ijk-kk2]-vmean[k-2],2) + cg2<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg3<TF>*std::pow(v[ijk    ]-vmean[k  ],2)) * dzhi4[k-1])
                                          + cg1<TF>*((cg0<TF>*std::pow(v[ijk-kk2]-vmean[k-2],2) + cg1<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg2<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + cg3<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2)) * dzhi4[k  ])
                                          + cg2<TF>*((cg0<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg1<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + cg2<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2) + cg3<TF>*std::pow(v[ijk+kk2]-vmean[k+2],2)) * dzhi4[k+1])
                                          + cg3<TF>*((cg0<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + cg1<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2) + cg2<TF>*std::pow(v[ijk+kk2]-vmean[k+2],2) + cg3<TF>*std::pow(v[ijk+kk3]-vmean[k+3],2)) * dzhi4[k+2]) ) * dzi4[k];

                    tke_visc[ijk] = 0.5 * visc * ( cg0<TF>*((cg0<TF>*std::pow(wz[ijk-kk3],2) + cg1<TF>*std::pow(wz[ijk-kk2],2) + cg2<TF>*std::pow(wz[ijk-kk1],2) + cg3<TF>*std::pow(wz[ijk    ],2)) * dzhi4[k-1])
                                                 + cg1<TF>*((cg0<TF>*std::pow(wz[ijk-kk2],2) + cg1<TF>*std::pow(wz[ijk-kk1],2) + cg2<TF>*std::pow(wz[ijk    ],2) + cg3<TF>*std::pow(wz[ijk+kk1],2)) * dzhi4[k  ])
                                                 + cg2<TF>*((cg0<TF>*std::pow(wz[ijk-kk1],2) + cg1<TF>*std::pow(wz[ijk    ],2) + cg2<TF>*std::pow(wz[ijk+kk1],2) + cg3<TF>*std::pow(wz[ijk+kk2],2)) * dzhi4[k+1])
                                                 + cg3<TF>*((cg0<TF>*std::pow(wz[ijk    ],2) + cg1<TF>*std::pow(wz[ijk+kk1],2) + cg2<TF>*std::pow(wz[ijk+kk2],2) + cg3<TF>*std::pow(wz[ijk+kk3],2)) * dzhi4[k+2]) ) * dzi4[k];

                    tke_visc[ijk] += 0.5*(u2_visc[ijk] + v2_visc[ijk]);
                }

        // top boundary
        k = kend-1;
        for (int j=jstart; j<jend; ++j)
    #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                u2_visc[ijk] = visc * ( cg0<TF>*((cg0<TF>*std::pow(u[ijk-kk3]-umean[k-3],2) + cg1<TF>*std::pow(u[ijk-kk2]-umean[k-2],2) + cg2<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + cg3<TF>*std::pow(u[ijk    ]-umean[k  ],2)) * dzhi4[k-1])
                                      + cg1<TF>*((cg0<TF>*std::pow(u[ijk-kk2]-umean[k-2],2) + cg1<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + cg2<TF>*std::pow(u[ijk    ]-umean[k  ],2) + cg3<TF>*std::pow(u[ijk+kk1]-umean[k+1],2)) * dzhi4[k  ])
                                      + cg2<TF>*((cg0<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + cg1<TF>*std::pow(u[ijk    ]-umean[k  ],2) + cg2<TF>*std::pow(u[ijk+kk1]-umean[k+1],2) + cg3<TF>*std::pow(u[ijk+kk2]-umean[k+2],2)) * dzhi4[k+1])
                                      + cg3<TF>*((tg0<TF>*std::pow(u[ijk-kk1]-umean[k-1],2) + tg1<TF>*std::pow(u[ijk    ]-umean[k  ],2) + tg2<TF>*std::pow(u[ijk+kk1]-umean[k+1],2) + tg3<TF>*std::pow(u[ijk+kk2]-umean[k+2],2)) * dzhi4[k+2]) ) * dzi4[k];

                v2_visc[ijk] = visc * ( cg0<TF>*((cg0<TF>*std::pow(v[ijk-kk3]-vmean[k-3],2) + cg1<TF>*std::pow(v[ijk-kk2]-vmean[k-2],2) + cg2<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg3<TF>*std::pow(v[ijk    ]-vmean[k  ],2)) * dzhi4[k-1])
                                      + cg1<TF>*((cg0<TF>*std::pow(v[ijk-kk2]-vmean[k-2],2) + cg1<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg2<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + cg3<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2)) * dzhi4[k  ])
                                      + cg2<TF>*((cg0<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + cg1<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + cg2<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2) + cg3<TF>*std::pow(v[ijk+kk2]-vmean[k+2],2)) * dzhi4[k+1])
                                      + cg3<TF>*((tg0<TF>*std::pow(v[ijk-kk1]-vmean[k-1],2) + tg1<TF>*std::pow(v[ijk    ]-vmean[k  ],2) + tg2<TF>*std::pow(v[ijk+kk1]-vmean[k+1],2) + tg3<TF>*std::pow(v[ijk+kk2]-vmean[k+2],2)) * dzhi4[k+2]) ) * dzi4[k];

                tke_visc[ijk] = 0.5 * visc * ( cg0<TF>*((cg0<TF>*std::pow(wz[ijk-kk3],2) + cg1<TF>*std::pow(wz[ijk-kk2],2) + cg2<TF>*std::pow(wz[ijk-kk1],2) + cg3<TF>*std::pow(wz[ijk    ],2)) * dzhi4[k-1])
                                             + cg1<TF>*((cg0<TF>*std::pow(wz[ijk-kk2],2) + cg1<TF>*std::pow(wz[ijk-kk1],2) + cg2<TF>*std::pow(wz[ijk    ],2) + cg3<TF>*std::pow(wz[ijk+kk1],2)) * dzhi4[k  ])
                                             + cg2<TF>*((cg0<TF>*std::pow(wz[ijk-kk1],2) + cg1<TF>*std::pow(wz[ijk    ],2) + cg2<TF>*std::pow(wz[ijk+kk1],2) + cg3<TF>*std::pow(wz[ijk+kk2],2)) * dzhi4[k+1])
                                             + cg3<TF>*((tg0<TF>*std::pow(wz[ijk-kk1],2) + tg1<TF>*std::pow(wz[ijk    ],2) + tg2<TF>*std::pow(wz[ijk+kk1],2) + tg3<TF>*std::pow(wz[ijk+kk2],2)) * dzhi4[k+2]) ) * dzi4[k];

                tke_visc[ijk] += 0.5*(u2_visc[ijk] + v2_visc[ijk]);
            }

        // Calculate the viscous transport of vertical velocity variance and the fluxes.
        // Interpolate

        // Bottom boundary.
        k = kstart;
        for (int j=jstart; j<jend; ++j)
    #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                w2_visc[ijk] = visc * ( bg0<TF>*((bg0<TF>*std::pow(w[ijk-kk1],2) + bg1<TF>*std::pow(w[ijk    ],2) + bg2<TF>*std::pow(w[ijk+kk1],2) + bg3<TF>*std::pow(w[ijk+kk2],2)) * dzi4[k-1])
                                      + bg1<TF>*((cg0<TF>*std::pow(w[ijk-kk1],2) + cg1<TF>*std::pow(w[ijk    ],2) + cg2<TF>*std::pow(w[ijk+kk1],2) + cg3<TF>*std::pow(w[ijk+kk2],2)) * dzi4[k  ])
                                      + bg2<TF>*((cg0<TF>*std::pow(w[ijk    ],2) + cg1<TF>*std::pow(w[ijk+kk1],2) + cg2<TF>*std::pow(w[ijk+kk2],2) + cg3<TF>*std::pow(w[ijk+kk3],2)) * dzi4[k+1])
                                      + bg3<TF>*((cg0<TF>*std::pow(w[ijk+kk1],2) + cg1<TF>*std::pow(w[ijk+kk2],2) + cg2<TF>*std::pow(w[ijk+kk3],2) + cg3<TF>*std::pow(w[ijk+kk4],2)) * dzi4[k+2]) ) * dzhi4bot;


                uw_visc[ijk] = ( ( visc
                                * ( bg0<TF>*( ( bg0<TF>*( uz[ijk        -kk1] * ( ci0<TF>*w[ijk-ii2    -kk1] + ci1<TF>*w[ijk-ii1    -kk1] + ci2<TF>*w[ijk        -kk1] + ci3<TF>*w[ijk+ii1    -kk1] ) )
                                          + bg1<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) )
                                          + bg2<TF>*( uz[ijk        +kk1] * ( ci0<TF>*w[ijk-ii2    +kk1] + ci1<TF>*w[ijk-ii1    +kk1] + ci2<TF>*w[ijk        +kk1] + ci3<TF>*w[ijk+ii1    +kk1] ) )
                                          + bg3<TF>*( uz[ijk        +kk2] * ( ci0<TF>*w[ijk-ii2    +kk2] + ci1<TF>*w[ijk-ii1    +kk2] + ci2<TF>*w[ijk        +kk2] + ci3<TF>*w[ijk+ii1    +kk2] ) ) )

                                        * dzi4[k-1] )

                                  + bg1<TF>*( ( cg0<TF>*( uz[ijk        -kk1] * ( ci0<TF>*w[ijk-ii2    -kk1] + ci1<TF>*w[ijk-ii1    -kk1] + ci2<TF>*w[ijk        -kk1] + ci3<TF>*w[ijk+ii1    -kk1] ) )
                                          + cg1<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) )
                                          + cg2<TF>*( uz[ijk        +kk1] * ( ci0<TF>*w[ijk-ii2    +kk1] + ci1<TF>*w[ijk-ii1    +kk1] + ci2<TF>*w[ijk        +kk1] + ci3<TF>*w[ijk+ii1    +kk1] ) )
                                          + cg3<TF>*( uz[ijk        +kk2] * ( ci0<TF>*w[ijk-ii2    +kk2] + ci1<TF>*w[ijk-ii1    +kk2] + ci2<TF>*w[ijk        +kk2] + ci3<TF>*w[ijk+ii1    +kk2] ) ) )

                                        * dzi4[k  ] )

                                  + bg2<TF>*( ( cg0<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) )
                                          + cg1<TF>*( uz[ijk        +kk1] * ( ci0<TF>*w[ijk-ii2    +kk1] + ci1<TF>*w[ijk-ii1    +kk1] + ci2<TF>*w[ijk        +kk1] + ci3<TF>*w[ijk+ii1    +kk1] ) )
                                          + cg2<TF>*( uz[ijk        +kk2] * ( ci0<TF>*w[ijk-ii2    +kk2] + ci1<TF>*w[ijk-ii1    +kk2] + ci2<TF>*w[ijk        +kk2] + ci3<TF>*w[ijk+ii1    +kk2] ) )
                                          + cg3<TF>*( uz[ijk        +kk3] * ( ci0<TF>*w[ijk-ii2    +kk3] + ci1<TF>*w[ijk-ii1    +kk3] + ci2<TF>*w[ijk        +kk3] + ci3<TF>*w[ijk+ii1    +kk3] ) ) )

                                        * dzi4[k+1] )

                                  + bg3<TF>*( ( cg0<TF>*( uz[ijk        +kk1] * ( ci0<TF>*w[ijk-ii2    +kk1] + ci1<TF>*w[ijk-ii1    +kk1] + ci2<TF>*w[ijk        +kk1] + ci3<TF>*w[ijk+ii1    +kk1] ) )
                                          + cg1<TF>*( uz[ijk        +kk2] * ( ci0<TF>*w[ijk-ii2    +kk2] + ci1<TF>*w[ijk-ii1    +kk2] + ci2<TF>*w[ijk        +kk2] + ci3<TF>*w[ijk+ii1    +kk2] ) )
                                          + cg2<TF>*( uz[ijk        +kk3] * ( ci0<TF>*w[ijk-ii2    +kk3] + ci1<TF>*w[ijk-ii1    +kk3] + ci2<TF>*w[ijk        +kk3] + ci3<TF>*w[ijk+ii1    +kk3] ) )
                                          + cg3<TF>*( uz[ijk        +kk4] * ( ci0<TF>*w[ijk-ii2    +kk4] + ci1<TF>*w[ijk-ii1    +kk4] + ci2<TF>*w[ijk        +kk4] + ci3<TF>*w[ijk+ii1    +kk4] ) ) )

                                        * dzi4[k+2] ) ) )


                              * dzhi4bot );

            }

        // Bottom boundary + 1.
        k = kstart+1;
        for (int j=jstart; j<jend; ++j)
    #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                w2_visc[ijk] = visc * ( cg0<TF>*((bg0<TF>*std::pow(w[ijk-kk2],2) + bg1<TF>*std::pow(w[ijk-kk1],2) + bg2<TF>*std::pow(w[ijk    ],2) + bg3<TF>*std::pow(w[ijk+kk1],2)) * dzi4[k-2])
                                      + cg1<TF>*((cg0<TF>*std::pow(w[ijk-kk2],2) + cg1<TF>*std::pow(w[ijk-kk1],2) + cg2<TF>*std::pow(w[ijk    ],2) + cg3<TF>*std::pow(w[ijk+kk1],2)) * dzi4[k-1])
                                      + cg2<TF>*((cg0<TF>*std::pow(w[ijk-kk1],2) + cg1<TF>*std::pow(w[ijk    ],2) + cg2<TF>*std::pow(w[ijk+kk1],2) + cg3<TF>*std::pow(w[ijk+kk2],2)) * dzi4[k  ])
                                      + cg3<TF>*((cg0<TF>*std::pow(w[ijk    ],2) + cg1<TF>*std::pow(w[ijk+kk1],2) + cg2<TF>*std::pow(w[ijk+kk2],2) + cg3<TF>*std::pow(w[ijk+kk3],2)) * dzi4[k+1]) ) * dzhi4[k];


                uw_visc[ijk] = ( ( visc


                                * ( cg0<TF>*( ( bg0<TF>*( uz[ijk        -kk2] * ( ci0<TF>*w[ijk-ii2    -kk2] + ci1<TF>*w[ijk-ii1    -kk2] + ci2<TF>*w[ijk        -kk2] + ci3<TF>*w[ijk+ii1    -kk2] ) )
                                          + bg1<TF>*( uz[ijk        -kk1] * ( ci0<TF>*w[ijk-ii2    -kk1] + ci1<TF>*w[ijk-ii1    -kk1] + ci2<TF>*w[ijk        -kk1] + ci3<TF>*w[ijk+ii1    -kk1] ) )
                                          + bg2<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) )
                                          + bg3<TF>*( uz[ijk        +kk1] * ( ci0<TF>*w[ijk-ii2    +kk1] + ci1<TF>*w[ijk-ii1    +kk1] + ci2<TF>*w[ijk        +kk1] + ci3<TF>*w[ijk+ii1    +kk1] ) ) )

                                        * dzi4[k-2] )

                                  + cg1<TF>*( ( cg0<TF>*( uz[ijk        -kk2] * ( ci0<TF>*w[ijk-ii2    -kk2] + ci1<TF>*w[ijk-ii1    -kk2] + ci2<TF>*w[ijk        -kk2] + ci3<TF>*w[ijk+ii1    -kk2] ) )
                                          + cg1<TF>*( uz[ijk        -kk1] * ( ci0<TF>*w[ijk-ii2    -kk1] + ci1<TF>*w[ijk-ii1    -kk1] + ci2<TF>*w[ijk        -kk1] + ci3<TF>*w[ijk+ii1    -kk1] ) )
                                          + cg2<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) )
                                          + cg3<TF>*( uz[ijk        +kk1] * ( ci0<TF>*w[ijk-ii2    +kk1] + ci1<TF>*w[ijk-ii1    +kk1] + ci2<TF>*w[ijk        +kk1] + ci3<TF>*w[ijk+ii1    +kk1] ) ) )

                                        * dzi4[k-1] )

                                  + cg2<TF>*( ( cg0<TF>*( uz[ijk        -kk1] * ( ci0<TF>*w[ijk-ii2    -kk1] + ci1<TF>*w[ijk-ii1    -kk1] + ci2<TF>*w[ijk        -kk1] + ci3<TF>*w[ijk+ii1    -kk1] ) )
                                          + cg1<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) )
                                          + cg2<TF>*( uz[ijk        +kk1] * ( ci0<TF>*w[ijk-ii2    +kk1] + ci1<TF>*w[ijk-ii1    +kk1] + ci2<TF>*w[ijk        +kk1] + ci3<TF>*w[ijk+ii1    +kk1] ) )
                                          + cg3<TF>*( uz[ijk        +kk2] * ( ci0<TF>*w[ijk-ii2    +kk2] + ci1<TF>*w[ijk-ii1    +kk2] + ci2<TF>*w[ijk        +kk2] + ci3<TF>*w[ijk+ii1    +kk2] ) ) )

                                        * dzi4[k  ] )

                                  + cg3<TF>*( ( cg0<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) )
                                          + cg1<TF>*( uz[ijk        +kk1] * ( ci0<TF>*w[ijk-ii2    +kk1] + ci1<TF>*w[ijk-ii1    +kk1] + ci2<TF>*w[ijk        +kk1] + ci3<TF>*w[ijk+ii1    +kk1] ) )
                                          + cg2<TF>*( uz[ijk        +kk2] * ( ci0<TF>*w[ijk-ii2    +kk2] + ci1<TF>*w[ijk-ii1    +kk2] + ci2<TF>*w[ijk        +kk2] + ci3<TF>*w[ijk+ii1    +kk2] ) )
                                          + cg3<TF>*( uz[ijk        +kk3] * ( ci0<TF>*w[ijk-ii2    +kk3] + ci1<TF>*w[ijk-ii1    +kk3] + ci2<TF>*w[ijk        +kk3] + ci3<TF>*w[ijk+ii1    +kk3] ) ) )

                                        * dzi4[k+1] ) ) )


                              * dzhi4[k  ] );

            }

        // Interior.
        for (int k=kstart+2; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
    #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    w2_visc[ijk] = visc * ( cg0<TF>*((cg0<TF>*std::pow(w[ijk-kk3],2) + cg1<TF>*std::pow(w[ijk-kk2],2) + cg2<TF>*std::pow(w[ijk-kk1],2) + cg3<TF>*std::pow(w[ijk    ],2)) * dzi4[k-2])
                                          + cg1<TF>*((cg0<TF>*std::pow(w[ijk-kk2],2) + cg1<TF>*std::pow(w[ijk-kk1],2) + cg2<TF>*std::pow(w[ijk    ],2) + cg3<TF>*std::pow(w[ijk+kk1],2)) * dzi4[k-1])
                                          + cg2<TF>*((cg0<TF>*std::pow(w[ijk-kk1],2) + cg1<TF>*std::pow(w[ijk    ],2) + cg2<TF>*std::pow(w[ijk+kk1],2) + cg3<TF>*std::pow(w[ijk+kk2],2)) * dzi4[k  ])
                                          + cg3<TF>*((cg0<TF>*std::pow(w[ijk    ],2) + cg1<TF>*std::pow(w[ijk+kk1],2) + cg2<TF>*std::pow(w[ijk+kk2],2) + cg3<TF>*std::pow(w[ijk+kk3],2)) * dzi4[k+1]) ) * dzhi4[k];


                    uw_visc[ijk] = ( ( visc


                                    * ( cg0<TF>*( ( cg0<TF>*( uz[ijk        -kk3] * ( ci0<TF>*w[ijk-ii2    -kk3] + ci1<TF>*w[ijk-ii1    -kk3] + ci2<TF>*w[ijk        -kk3] + ci3<TF>*w[ijk+ii1    -kk3] ) )
                                              + cg1<TF>*( uz[ijk        -kk2] * ( ci0<TF>*w[ijk-ii2    -kk2] + ci1<TF>*w[ijk-ii1    -kk2] + ci2<TF>*w[ijk        -kk2] + ci3<TF>*w[ijk+ii1    -kk2] ) )
                                              + cg2<TF>*( uz[ijk        -kk1] * ( ci0<TF>*w[ijk-ii2    -kk1] + ci1<TF>*w[ijk-ii1    -kk1] + ci2<TF>*w[ijk        -kk1] + ci3<TF>*w[ijk+ii1    -kk1] ) )
                                              + cg3<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) ) )

                                            * dzi4[k-2] )

                                      + cg1<TF>*( ( cg0<TF>*( uz[ijk        -kk2] * ( ci0<TF>*w[ijk-ii2    -kk2] + ci1<TF>*w[ijk-ii1    -kk2] + ci2<TF>*w[ijk        -kk2] + ci3<TF>*w[ijk+ii1    -kk2] ) )
                                              + cg1<TF>*( uz[ijk        -kk1] * ( ci0<TF>*w[ijk-ii2    -kk1] + ci1<TF>*w[ijk-ii1    -kk1] + ci2<TF>*w[ijk        -kk1] + ci3<TF>*w[ijk+ii1    -kk1] ) )
                                              + cg2<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) )
                                              + cg3<TF>*( uz[ijk        +kk1] * ( ci0<TF>*w[ijk-ii2    +kk1] + ci1<TF>*w[ijk-ii1    +kk1] + ci2<TF>*w[ijk        +kk1] + ci3<TF>*w[ijk+ii1    +kk1] ) ) )

                                            * dzi4[k-1] )

                                      + cg2<TF>*( ( cg0<TF>*( uz[ijk        -kk1] * ( ci0<TF>*w[ijk-ii2    -kk1] + ci1<TF>*w[ijk-ii1    -kk1] + ci2<TF>*w[ijk        -kk1] + ci3<TF>*w[ijk+ii1    -kk1] ) )
                                              + cg1<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) )
                                              + cg2<TF>*( uz[ijk        +kk1] * ( ci0<TF>*w[ijk-ii2    +kk1] + ci1<TF>*w[ijk-ii1    +kk1] + ci2<TF>*w[ijk        +kk1] + ci3<TF>*w[ijk+ii1    +kk1] ) )
                                              + cg3<TF>*( uz[ijk        +kk2] * ( ci0<TF>*w[ijk-ii2    +kk2] + ci1<TF>*w[ijk-ii1    +kk2] + ci2<TF>*w[ijk        +kk2] + ci3<TF>*w[ijk+ii1    +kk2] ) ) )

                                            * dzi4[k  ] )

                                      + cg3<TF>*( ( cg0<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) )
                                              + cg1<TF>*( uz[ijk        +kk1] * ( ci0<TF>*w[ijk-ii2    +kk1] + ci1<TF>*w[ijk-ii1    +kk1] + ci2<TF>*w[ijk        +kk1] + ci3<TF>*w[ijk+ii1    +kk1] ) )
                                              + cg2<TF>*( uz[ijk        +kk2] * ( ci0<TF>*w[ijk-ii2    +kk2] + ci1<TF>*w[ijk-ii1    +kk2] + ci2<TF>*w[ijk        +kk2] + ci3<TF>*w[ijk+ii1    +kk2] ) )
                                              + cg3<TF>*( uz[ijk        +kk3] * ( ci0<TF>*w[ijk-ii2    +kk3] + ci1<TF>*w[ijk-ii1    +kk3] + ci2<TF>*w[ijk        +kk3] + ci3<TF>*w[ijk+ii1    +kk3] ) ) )

                                            * dzi4[k+1] ) ) )


                                  * dzhi4[k  ] );
                }

        // Top boundary - 1.
        k = kend-1;
        for (int j=jstart; j<jend; ++j)
    #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                w2_visc[ijk] = visc * ( cg0<TF>*((cg0<TF>*std::pow(w[ijk-kk3],2) + cg1<TF>*std::pow(w[ijk-kk2],2) + cg2<TF>*std::pow(w[ijk-kk1],2) + cg3<TF>*std::pow(w[ijk    ],2)) * dzi4[k-2])
                                      + cg1<TF>*((cg0<TF>*std::pow(w[ijk-kk2],2) + cg1<TF>*std::pow(w[ijk-kk1],2) + cg2<TF>*std::pow(w[ijk    ],2) + cg3<TF>*std::pow(w[ijk+kk1],2)) * dzi4[k-1])
                                      + cg2<TF>*((cg0<TF>*std::pow(w[ijk-kk1],2) + cg1<TF>*std::pow(w[ijk    ],2) + cg2<TF>*std::pow(w[ijk+kk1],2) + cg3<TF>*std::pow(w[ijk+kk2],2)) * dzi4[k  ])
                                      + cg3<TF>*((tg0<TF>*std::pow(w[ijk-kk1],2) + tg1<TF>*std::pow(w[ijk    ],2) + tg2<TF>*std::pow(w[ijk+kk1],2) + tg3<TF>*std::pow(w[ijk+kk2],2)) * dzi4[k+1]) ) * dzhi4[k];


                uw_visc[ijk] = ( ( visc


                                * ( cg0<TF>*( ( cg0<TF>*( uz[ijk        -kk3] * ( ci0<TF>*w[ijk-ii2    -kk3] + ci1<TF>*w[ijk-ii1    -kk3] + ci2<TF>*w[ijk        -kk3] + ci3<TF>*w[ijk+ii1    -kk3] ) )
                                          + cg1<TF>*( uz[ijk        -kk2] * ( ci0<TF>*w[ijk-ii2    -kk2] + ci1<TF>*w[ijk-ii1    -kk2] + ci2<TF>*w[ijk        -kk2] + ci3<TF>*w[ijk+ii1    -kk2] ) )
                                          + cg2<TF>*( uz[ijk        -kk1] * ( ci0<TF>*w[ijk-ii2    -kk1] + ci1<TF>*w[ijk-ii1    -kk1] + ci2<TF>*w[ijk        -kk1] + ci3<TF>*w[ijk+ii1    -kk1] ) )
                                          + cg3<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) ) )

                                        * dzi4[k-2] )

                                  + cg1<TF>*( ( cg0<TF>*( uz[ijk        -kk2] * ( ci0<TF>*w[ijk-ii2    -kk2] + ci1<TF>*w[ijk-ii1    -kk2] + ci2<TF>*w[ijk        -kk2] + ci3<TF>*w[ijk+ii1    -kk2] ) )
                                          + cg1<TF>*( uz[ijk        -kk1] * ( ci0<TF>*w[ijk-ii2    -kk1] + ci1<TF>*w[ijk-ii1    -kk1] + ci2<TF>*w[ijk        -kk1] + ci3<TF>*w[ijk+ii1    -kk1] ) )
                                          + cg2<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) )
                                          + cg3<TF>*( uz[ijk        +kk1] * ( ci0<TF>*w[ijk-ii2    +kk1] + ci1<TF>*w[ijk-ii1    +kk1] + ci2<TF>*w[ijk        +kk1] + ci3<TF>*w[ijk+ii1    +kk1] ) ) )

                                        * dzi4[k-1] )

                                  + cg2<TF>*( ( cg0<TF>*( uz[ijk        -kk1] * ( ci0<TF>*w[ijk-ii2    -kk1] + ci1<TF>*w[ijk-ii1    -kk1] + ci2<TF>*w[ijk        -kk1] + ci3<TF>*w[ijk+ii1    -kk1] ) )
                                          + cg1<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) )
                                          + cg2<TF>*( uz[ijk        +kk1] * ( ci0<TF>*w[ijk-ii2    +kk1] + ci1<TF>*w[ijk-ii1    +kk1] + ci2<TF>*w[ijk        +kk1] + ci3<TF>*w[ijk+ii1    +kk1] ) )
                                          + cg3<TF>*( uz[ijk        +kk2] * ( ci0<TF>*w[ijk-ii2    +kk2] + ci1<TF>*w[ijk-ii1    +kk2] + ci2<TF>*w[ijk        +kk2] + ci3<TF>*w[ijk+ii1    +kk2] ) ) )

                                        * dzi4[k  ] )

                                  + cg3<TF>*( ( tg0<TF>*( uz[ijk        -kk1] * ( ci0<TF>*w[ijk-ii2    -kk1] + ci1<TF>*w[ijk-ii1    -kk1] + ci2<TF>*w[ijk        -kk1] + ci3<TF>*w[ijk+ii1    -kk1] ) )
                                          + tg1<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) )
                                          + tg2<TF>*( uz[ijk        +kk1] * ( ci0<TF>*w[ijk-ii2    +kk1] + ci1<TF>*w[ijk-ii1    +kk1] + ci2<TF>*w[ijk        +kk1] + ci3<TF>*w[ijk+ii1    +kk1] ) )
                                          + tg3<TF>*( uz[ijk        +kk2] * ( ci0<TF>*w[ijk-ii2    +kk2] + ci1<TF>*w[ijk-ii1    +kk2] + ci2<TF>*w[ijk        +kk2] + ci3<TF>*w[ijk+ii1    +kk2] ) ) )

                                        * dzi4[k+1] ) ) )


                              * dzhi4[k  ] );
            }

        // Top boundary.
        k = kend;
        for (int j=jstart; j<jend; ++j)
    #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                w2_visc[ijk] = visc * ( tg0<TF>*((cg0<TF>*std::pow(w[ijk-kk4],2) + cg1<TF>*std::pow(w[ijk-kk3],2) + cg2<TF>*std::pow(w[ijk-kk2],2) + cg3<TF>*std::pow(w[ijk-kk1],2)) * dzi4[k-3])
                                      + tg1<TF>*((cg0<TF>*std::pow(w[ijk-kk3],2) + cg1<TF>*std::pow(w[ijk-kk2],2) + cg2<TF>*std::pow(w[ijk-kk1],2) + cg3<TF>*std::pow(w[ijk    ],2)) * dzi4[k-2])
                                      + tg2<TF>*((cg0<TF>*std::pow(w[ijk-kk2],2) + cg1<TF>*std::pow(w[ijk-kk1],2) + cg2<TF>*std::pow(w[ijk    ],2) + cg3<TF>*std::pow(w[ijk+kk1],2)) * dzi4[k-1])
                                      + tg3<TF>*((tg0<TF>*std::pow(w[ijk-kk2],2) + tg1<TF>*std::pow(w[ijk-kk1],2) + tg2<TF>*std::pow(w[ijk    ],2) + tg3<TF>*std::pow(w[ijk+kk1],2)) * dzi4[k  ]) ) * dzhi4top;


                uw_visc[ijk] += ( ( visc


                                * ( tg0<TF>*( ( cg0<TF>*( uz[ijk        -kk4] * ( ci0<TF>*w[ijk-ii2    -kk4] + ci1<TF>*w[ijk-ii1    -kk4] + ci2<TF>*w[ijk        -kk4] + ci3<TF>*w[ijk+ii1    -kk4] ) )
                                          + cg1<TF>*( uz[ijk        -kk3] * ( ci0<TF>*w[ijk-ii2    -kk3] + ci1<TF>*w[ijk-ii1    -kk3] + ci2<TF>*w[ijk        -kk3] + ci3<TF>*w[ijk+ii1    -kk3] ) )
                                          + cg2<TF>*( uz[ijk        -kk2] * ( ci0<TF>*w[ijk-ii2    -kk2] + ci1<TF>*w[ijk-ii1    -kk2] + ci2<TF>*w[ijk        -kk2] + ci3<TF>*w[ijk+ii1    -kk2] ) )
                                          + cg3<TF>*( uz[ijk        -kk1] * ( ci0<TF>*w[ijk-ii2    -kk1] + ci1<TF>*w[ijk-ii1    -kk1] + ci2<TF>*w[ijk        -kk1] + ci3<TF>*w[ijk+ii1    -kk1] ) ) )

                                        * dzi4[k-3] )

                                  + tg1<TF>*( ( cg0<TF>*( uz[ijk        -kk3] * ( ci0<TF>*w[ijk-ii2    -kk3] + ci1<TF>*w[ijk-ii1    -kk3] + ci2<TF>*w[ijk        -kk3] + ci3<TF>*w[ijk+ii1    -kk3] ) )
                                          + cg1<TF>*( uz[ijk        -kk2] * ( ci0<TF>*w[ijk-ii2    -kk2] + ci1<TF>*w[ijk-ii1    -kk2] + ci2<TF>*w[ijk        -kk2] + ci3<TF>*w[ijk+ii1    -kk2] ) )
                                          + cg2<TF>*( uz[ijk        -kk1] * ( ci0<TF>*w[ijk-ii2    -kk1] + ci1<TF>*w[ijk-ii1    -kk1] + ci2<TF>*w[ijk        -kk1] + ci3<TF>*w[ijk+ii1    -kk1] ) )
                                          + cg3<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) ) )

                                        * dzi4[k-2] )

                                  + tg2<TF>*( ( cg0<TF>*( uz[ijk        -kk2] * ( ci0<TF>*w[ijk-ii2    -kk2] + ci1<TF>*w[ijk-ii1    -kk2] + ci2<TF>*w[ijk        -kk2] + ci3<TF>*w[ijk+ii1    -kk2] ) )
                                          + cg1<TF>*( uz[ijk        -kk1] * ( ci0<TF>*w[ijk-ii2    -kk1] + ci1<TF>*w[ijk-ii1    -kk1] + ci2<TF>*w[ijk        -kk1] + ci3<TF>*w[ijk+ii1    -kk1] ) )
                                          + cg2<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) )
                                          + cg3<TF>*( uz[ijk        +kk1] * ( ci0<TF>*w[ijk-ii2    +kk1] + ci1<TF>*w[ijk-ii1    +kk1] + ci2<TF>*w[ijk        +kk1] + ci3<TF>*w[ijk+ii1    +kk1] ) ) )

                                        * dzi4[k-1] )

                                  + tg3<TF>*( ( tg0<TF>*( uz[ijk        -kk2] * ( ci0<TF>*w[ijk-ii2    -kk2] + ci1<TF>*w[ijk-ii1    -kk2] + ci2<TF>*w[ijk        -kk2] + ci3<TF>*w[ijk+ii1    -kk2] ) )
                                          + tg1<TF>*( uz[ijk        -kk1] * ( ci0<TF>*w[ijk-ii2    -kk1] + ci1<TF>*w[ijk-ii1    -kk1] + ci2<TF>*w[ijk        -kk1] + ci3<TF>*w[ijk+ii1    -kk1] ) )
                                          + tg2<TF>*( uz[ijk            ] * ( ci0<TF>*w[ijk-ii2        ] + ci1<TF>*w[ijk-ii1        ] + ci2<TF>*w[ijk            ] + ci3<TF>*w[ijk+ii1        ] ) )
                                          + tg3<TF>*( uz[ijk        +kk1] * ( ci0<TF>*w[ijk-ii2    +kk1] + ci1<TF>*w[ijk-ii1    +kk1] + ci2<TF>*w[ijk        +kk1] + ci3<TF>*w[ijk+ii1    +kk1] ) ) )

                                        * dzi4[k  ] ) ) )


                              * dzhi4top );

            }
    }

    template<typename TF>
    void calc_tke_budget_diss(
            TF* restrict u2_diss, TF* restrict v2_diss, TF* restrict w2_diss, TF* restrict tke_diss, TF* restrict uw_diss,
            const TF* restrict u, const TF* restrict v, const TF* restrict w,
            const TF* restrict umean, const TF* restrict vmean,
            const TF* restrict dzi4, const TF* restrict dzhi4,
            const TF dxi, const TF dyi, const TF dzhi4bot, const TF dzhi4top,
            const TF visc,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        using namespace Finite_difference::O4;

        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj1 = 1*icells;
        const int jj2 = 2*icells;
        const int jj3 = 3*icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;
        const int kk3 = 3*ijcells;
        const int kk4 = 4*ijcells;

        // 6. CALCULATE THE DISSIPATION TERM

        // bottom boundary
        int k = kstart;
        for (int j=jstart; j<jend; ++j)
             #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                u2_diss[ijk] = -2.*visc * (
                                 std::pow( ( cg0<TF>*((ci0<TF>*(u[ijk-ii3]-umean[k]) + ci1<TF>*(u[ijk-ii2]-umean[k]) + ci2<TF>*(u[ijk-ii1]-umean[k]) + ci3<TF>*(u[ijk    ]-umean[k])))
                                           + cg1<TF>*((ci0<TF>*(u[ijk-ii2]-umean[k]) + ci1<TF>*(u[ijk-ii1]-umean[k]) + ci2<TF>*(u[ijk    ]-umean[k]) + ci3<TF>*(u[ijk+ii1]-umean[k])))
                                           + cg2<TF>*((ci0<TF>*(u[ijk-ii1]-umean[k]) + ci1<TF>*(u[ijk    ]-umean[k]) + ci2<TF>*(u[ijk+ii1]-umean[k]) + ci3<TF>*(u[ijk+ii2]-umean[k])))
                                           + cg3<TF>*((ci0<TF>*(u[ijk    ]-umean[k]) + ci1<TF>*(u[ijk+ii1]-umean[k]) + ci2<TF>*(u[ijk+ii2]-umean[k]) + ci3<TF>*(u[ijk+ii3]-umean[k]))) ) * dxi, 2)

                               + std::pow( ( cg0<TF>*((ci0<TF>*(u[ijk-jj3]-umean[k]) + ci1<TF>*(u[ijk-jj2]-umean[k]) + ci2<TF>*(u[ijk-jj1]-umean[k]) + ci3<TF>*(u[ijk    ]-umean[k])))
                                           + cg1<TF>*((ci0<TF>*(u[ijk-jj2]-umean[k]) + ci1<TF>*(u[ijk-jj1]-umean[k]) + ci2<TF>*(u[ijk    ]-umean[k]) + ci3<TF>*(u[ijk+jj1]-umean[k])))
                                           + cg2<TF>*((ci0<TF>*(u[ijk-jj1]-umean[k]) + ci1<TF>*(u[ijk    ]-umean[k]) + ci2<TF>*(u[ijk+jj1]-umean[k]) + ci3<TF>*(u[ijk+jj2]-umean[k])))
                                           + cg3<TF>*((ci0<TF>*(u[ijk    ]-umean[k]) + ci1<TF>*(u[ijk+jj1]-umean[k]) + ci2<TF>*(u[ijk+jj2]-umean[k]) + ci3<TF>*(u[ijk+jj3]-umean[k]))) ) * dyi, 2)

                               + std::pow( ( cg0<TF>*((bi0<TF>*(u[ijk-kk2]-umean[k-2]) + bi1<TF>*(u[ijk-kk1]-umean[k-1]) + bi2<TF>*(u[ijk    ]-umean[k  ]) + bi3<TF>*(u[ijk+kk1]-umean[k+1])))
                                           + cg1<TF>*((ci0<TF>*(u[ijk-kk2]-umean[k-2]) + ci1<TF>*(u[ijk-kk1]-umean[k-1]) + ci2<TF>*(u[ijk    ]-umean[k  ]) + ci3<TF>*(u[ijk+kk1]-umean[k+1])))
                                           + cg2<TF>*((ci0<TF>*(u[ijk-kk1]-umean[k-1]) + ci1<TF>*(u[ijk    ]-umean[k  ]) + ci2<TF>*(u[ijk+kk1]-umean[k+1]) + ci3<TF>*(u[ijk+kk2]-umean[k+2])))
                                           + cg3<TF>*((ci0<TF>*(u[ijk    ]-umean[k  ]) + ci1<TF>*(u[ijk+kk1]-umean[k+1]) + ci2<TF>*(u[ijk+kk2]-umean[k+2]) + ci3<TF>*(u[ijk+kk3]-umean[k+3]))) ) * dzi4[k], 2) );

                v2_diss[ijk] = -2.*visc * (
                                 std::pow( ( cg0<TF>*((ci0<TF>*(v[ijk-ii3]-vmean[k]) + ci1<TF>*(v[ijk-ii2]-vmean[k]) + ci2<TF>*(v[ijk-ii1]-vmean[k]) + ci3<TF>*(v[ijk    ]-vmean[k])))
                                           + cg1<TF>*((ci0<TF>*(v[ijk-ii2]-vmean[k]) + ci1<TF>*(v[ijk-ii1]-vmean[k]) + ci2<TF>*(v[ijk    ]-vmean[k]) + ci3<TF>*(v[ijk+ii1]-vmean[k])))
                                           + cg2<TF>*((ci0<TF>*(v[ijk-ii1]-vmean[k]) + ci1<TF>*(v[ijk    ]-vmean[k]) + ci2<TF>*(v[ijk+ii1]-vmean[k]) + ci3<TF>*(v[ijk+ii2]-vmean[k])))
                                           + cg3<TF>*((ci0<TF>*(v[ijk    ]-vmean[k]) + ci1<TF>*(v[ijk+ii1]-vmean[k]) + ci2<TF>*(v[ijk+ii2]-vmean[k]) + ci3<TF>*(v[ijk+ii3]-vmean[k]))) ) * dxi, 2)

                               + std::pow( ( cg0<TF>*((ci0<TF>*(v[ijk-jj3]-vmean[k]) + ci1<TF>*(v[ijk-jj2]-vmean[k]) + ci2<TF>*(v[ijk-jj1]-vmean[k]) + ci3<TF>*(v[ijk    ]-vmean[k])))
                                           + cg1<TF>*((ci0<TF>*(v[ijk-jj2]-vmean[k]) + ci1<TF>*(v[ijk-jj1]-vmean[k]) + ci2<TF>*(v[ijk    ]-vmean[k]) + ci3<TF>*(v[ijk+jj1]-vmean[k])))
                                           + cg2<TF>*((ci0<TF>*(v[ijk-jj1]-vmean[k]) + ci1<TF>*(v[ijk    ]-vmean[k]) + ci2<TF>*(v[ijk+jj1]-vmean[k]) + ci3<TF>*(v[ijk+jj2]-vmean[k])))
                                           + cg3<TF>*((ci0<TF>*(v[ijk    ]-vmean[k]) + ci1<TF>*(v[ijk+jj1]-vmean[k]) + ci2<TF>*(v[ijk+jj2]-vmean[k]) + ci3<TF>*(v[ijk+jj3]-vmean[k]))) ) * dyi, 2)

                               + std::pow( ( cg0<TF>*((bi0<TF>*(v[ijk-kk2]-vmean[k-2]) + bi1<TF>*(v[ijk-kk1]-vmean[k-1]) + bi2<TF>*(v[ijk    ]-vmean[k  ]) + bi3<TF>*(v[ijk+kk1]-vmean[k+1])))
                                           + cg1<TF>*((ci0<TF>*(v[ijk-kk2]-vmean[k-2]) + ci1<TF>*(v[ijk-kk1]-vmean[k-1]) + ci2<TF>*(v[ijk    ]-vmean[k  ]) + ci3<TF>*(v[ijk+kk1]-vmean[k+1])))
                                           + cg2<TF>*((ci0<TF>*(v[ijk-kk1]-vmean[k-1]) + ci1<TF>*(v[ijk    ]-vmean[k  ]) + ci2<TF>*(v[ijk+kk1]-vmean[k+1]) + ci3<TF>*(v[ijk+kk2]-vmean[k+2])))
                                           + cg3<TF>*((ci0<TF>*(v[ijk    ]-vmean[k  ]) + ci1<TF>*(v[ijk+kk1]-vmean[k+1]) + ci2<TF>*(v[ijk+kk2]-vmean[k+2]) + ci3<TF>*(v[ijk+kk3]-vmean[k+3]))) ) * dzi4[k], 2) );

                tke_diss[ijk] = -visc * (
                                        std::pow( (cg0<TF>*w[ijk-ii1] + cg1<TF>*w[ijk] + cg2<TF>*w[ijk+ii1] + cg3<TF>*w[ijk+ii2]) * dxi, 2)
                                      + std::pow( (cg0<TF>*w[ijk-jj1] + cg1<TF>*w[ijk] + cg2<TF>*w[ijk+jj1] + cg3<TF>*w[ijk+jj2]) * dyi, 2)
                                      + std::pow( (cg0<TF>*w[ijk-kk1] + cg1<TF>*w[ijk] + cg2<TF>*w[ijk+kk1] + cg3<TF>*w[ijk+kk2]) * dzi4[k], 2) );

                tke_diss[ijk] += 0.5*(u2_diss[ijk] + v2_diss[ijk]);
            }

        // interior
        for (int k=kstart+1; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    u2_diss[ijk] = -2.*visc * (
                                     std::pow( ( cg0<TF>*((ci0<TF>*(u[ijk-ii3]-umean[k]) + ci1<TF>*(u[ijk-ii2]-umean[k]) + ci2<TF>*(u[ijk-ii1]-umean[k]) + ci3<TF>*(u[ijk    ]-umean[k])))
                                               + cg1<TF>*((ci0<TF>*(u[ijk-ii2]-umean[k]) + ci1<TF>*(u[ijk-ii1]-umean[k]) + ci2<TF>*(u[ijk    ]-umean[k]) + ci3<TF>*(u[ijk+ii1]-umean[k])))
                                               + cg2<TF>*((ci0<TF>*(u[ijk-ii1]-umean[k]) + ci1<TF>*(u[ijk    ]-umean[k]) + ci2<TF>*(u[ijk+ii1]-umean[k]) + ci3<TF>*(u[ijk+ii2]-umean[k])))
                                               + cg3<TF>*((ci0<TF>*(u[ijk    ]-umean[k]) + ci1<TF>*(u[ijk+ii1]-umean[k]) + ci2<TF>*(u[ijk+ii2]-umean[k]) + ci3<TF>*(u[ijk+ii3]-umean[k]))) ) * dxi, 2)

                                   + std::pow( ( cg0<TF>*((ci0<TF>*(u[ijk-jj3]-umean[k]) + ci1<TF>*(u[ijk-jj2]-umean[k]) + ci2<TF>*(u[ijk-jj1]-umean[k]) + ci3<TF>*(u[ijk    ]-umean[k])))
                                               + cg1<TF>*((ci0<TF>*(u[ijk-jj2]-umean[k]) + ci1<TF>*(u[ijk-jj1]-umean[k]) + ci2<TF>*(u[ijk    ]-umean[k]) + ci3<TF>*(u[ijk+jj1]-umean[k])))
                                               + cg2<TF>*((ci0<TF>*(u[ijk-jj1]-umean[k]) + ci1<TF>*(u[ijk    ]-umean[k]) + ci2<TF>*(u[ijk+jj1]-umean[k]) + ci3<TF>*(u[ijk+jj2]-umean[k])))
                                               + cg3<TF>*((ci0<TF>*(u[ijk    ]-umean[k]) + ci1<TF>*(u[ijk+jj1]-umean[k]) + ci2<TF>*(u[ijk+jj2]-umean[k]) + ci3<TF>*(u[ijk+jj3]-umean[k]))) ) * dyi, 2)

                                   + std::pow( ( cg0<TF>*((ci0<TF>*(u[ijk-kk3]-umean[k-3]) + ci1<TF>*(u[ijk-kk2]-umean[k-2]) + ci2<TF>*(u[ijk-kk1]-umean[k-1]) + ci3<TF>*(u[ijk    ]-umean[k  ])))
                                               + cg1<TF>*((ci0<TF>*(u[ijk-kk2]-umean[k-2]) + ci1<TF>*(u[ijk-kk1]-umean[k-1]) + ci2<TF>*(u[ijk    ]-umean[k  ]) + ci3<TF>*(u[ijk+kk1]-umean[k+1])))
                                               + cg2<TF>*((ci0<TF>*(u[ijk-kk1]-umean[k-1]) + ci1<TF>*(u[ijk    ]-umean[k  ]) + ci2<TF>*(u[ijk+kk1]-umean[k+1]) + ci3<TF>*(u[ijk+kk2]-umean[k+2])))
                                               + cg3<TF>*((ci0<TF>*(u[ijk    ]-umean[k  ]) + ci1<TF>*(u[ijk+kk1]-umean[k+1]) + ci2<TF>*(u[ijk+kk2]-umean[k+2]) + ci3<TF>*(u[ijk+kk3]-umean[k+3]))) ) * dzi4[k], 2) );

                    v2_diss[ijk] = -2.*visc * (
                                     std::pow( ( cg0<TF>*((ci0<TF>*(v[ijk-ii3]-vmean[k]) + ci1<TF>*(v[ijk-ii2]-vmean[k]) + ci2<TF>*(v[ijk-ii1]-vmean[k]) + ci3<TF>*(v[ijk    ]-vmean[k])))
                                               + cg1<TF>*((ci0<TF>*(v[ijk-ii2]-vmean[k]) + ci1<TF>*(v[ijk-ii1]-vmean[k]) + ci2<TF>*(v[ijk    ]-vmean[k]) + ci3<TF>*(v[ijk+ii1]-vmean[k])))
                                               + cg2<TF>*((ci0<TF>*(v[ijk-ii1]-vmean[k]) + ci1<TF>*(v[ijk    ]-vmean[k]) + ci2<TF>*(v[ijk+ii1]-vmean[k]) + ci3<TF>*(v[ijk+ii2]-vmean[k])))
                                               + cg3<TF>*((ci0<TF>*(v[ijk    ]-vmean[k]) + ci1<TF>*(v[ijk+ii1]-vmean[k]) + ci2<TF>*(v[ijk+ii2]-vmean[k]) + ci3<TF>*(v[ijk+ii3]-vmean[k]))) ) * dxi, 2)

                                   + std::pow( ( cg0<TF>*((ci0<TF>*(v[ijk-jj3]-vmean[k]) + ci1<TF>*(v[ijk-jj2]-vmean[k]) + ci2<TF>*(v[ijk-jj1]-vmean[k]) + ci3<TF>*(v[ijk    ]-vmean[k])))
                                               + cg1<TF>*((ci0<TF>*(v[ijk-jj2]-vmean[k]) + ci1<TF>*(v[ijk-jj1]-vmean[k]) + ci2<TF>*(v[ijk    ]-vmean[k]) + ci3<TF>*(v[ijk+jj1]-vmean[k])))
                                               + cg2<TF>*((ci0<TF>*(v[ijk-jj1]-vmean[k]) + ci1<TF>*(v[ijk    ]-vmean[k]) + ci2<TF>*(v[ijk+jj1]-vmean[k]) + ci3<TF>*(v[ijk+jj2]-vmean[k])))
                                               + cg3<TF>*((ci0<TF>*(v[ijk    ]-vmean[k]) + ci1<TF>*(v[ijk+jj1]-vmean[k]) + ci2<TF>*(v[ijk+jj2]-vmean[k]) + ci3<TF>*(v[ijk+jj3]-vmean[k]))) ) * dyi, 2)

                                   + std::pow( ( cg0<TF>*((ci0<TF>*(v[ijk-kk3]-vmean[k-3]) + ci1<TF>*(v[ijk-kk2]-vmean[k-2]) + ci2<TF>*(v[ijk-kk1]-vmean[k-1]) + ci3<TF>*(v[ijk    ]-vmean[k  ])))
                                               + cg1<TF>*((ci0<TF>*(v[ijk-kk2]-vmean[k-2]) + ci1<TF>*(v[ijk-kk1]-vmean[k-1]) + ci2<TF>*(v[ijk    ]-vmean[k  ]) + ci3<TF>*(v[ijk+kk1]-vmean[k+1])))
                                               + cg2<TF>*((ci0<TF>*(v[ijk-kk1]-vmean[k-1]) + ci1<TF>*(v[ijk    ]-vmean[k  ]) + ci2<TF>*(v[ijk+kk1]-vmean[k+1]) + ci3<TF>*(v[ijk+kk2]-vmean[k+2])))
                                               + cg3<TF>*((ci0<TF>*(v[ijk    ]-vmean[k  ]) + ci1<TF>*(v[ijk+kk1]-vmean[k+1]) + ci2<TF>*(v[ijk+kk2]-vmean[k+2]) + ci3<TF>*(v[ijk+kk3]-vmean[k+3]))) ) * dzi4[k], 2) );

                    tke_diss[ijk] = -visc * (
                                       std::pow( (cg0<TF>*w[ijk-ii1] + cg1<TF>*w[ijk] + cg2<TF>*w[ijk+ii1] + cg3<TF>*w[ijk+ii2]) * dxi, 2)
                                     + std::pow( (cg0<TF>*w[ijk-jj1] + cg1<TF>*w[ijk] + cg2<TF>*w[ijk+jj1] + cg3<TF>*w[ijk+jj2]) * dyi, 2)
                                     + std::pow( (cg0<TF>*w[ijk-kk1] + cg1<TF>*w[ijk] + cg2<TF>*w[ijk+kk1] + cg3<TF>*w[ijk+kk2]) * dzi4[k], 2) );

                    tke_diss[ijk] += 0.5*(u2_diss[ijk] + v2_diss[ijk]);
                }

        // top boundary
        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                u2_diss[ijk] = - 2.*visc * (
                                 std::pow( ( cg0<TF>*((ci0<TF>*(u[ijk-ii3]-umean[k]) + ci1<TF>*(u[ijk-ii2]-umean[k]) + ci2<TF>*(u[ijk-ii1]-umean[k]) + ci3<TF>*(u[ijk    ]-umean[k])))
                                           + cg1<TF>*((ci0<TF>*(u[ijk-ii2]-umean[k]) + ci1<TF>*(u[ijk-ii1]-umean[k]) + ci2<TF>*(u[ijk    ]-umean[k]) + ci3<TF>*(u[ijk+ii1]-umean[k])))
                                           + cg2<TF>*((ci0<TF>*(u[ijk-ii1]-umean[k]) + ci1<TF>*(u[ijk    ]-umean[k]) + ci2<TF>*(u[ijk+ii1]-umean[k]) + ci3<TF>*(u[ijk+ii2]-umean[k])))
                                           + cg3<TF>*((ci0<TF>*(u[ijk    ]-umean[k]) + ci1<TF>*(u[ijk+ii1]-umean[k]) + ci2<TF>*(u[ijk+ii2]-umean[k]) + ci3<TF>*(u[ijk+ii3]-umean[k]))) ) * dxi, 2)

                               + std::pow( ( cg0<TF>*((ci0<TF>*(u[ijk-jj3]-umean[k]) + ci1<TF>*(u[ijk-jj2]-umean[k]) + ci2<TF>*(u[ijk-jj1]-umean[k]) + ci3<TF>*(u[ijk    ]-umean[k])))
                                           + cg1<TF>*((ci0<TF>*(u[ijk-jj2]-umean[k]) + ci1<TF>*(u[ijk-jj1]-umean[k]) + ci2<TF>*(u[ijk    ]-umean[k]) + ci3<TF>*(u[ijk+jj1]-umean[k])))
                                           + cg2<TF>*((ci0<TF>*(u[ijk-jj1]-umean[k]) + ci1<TF>*(u[ijk    ]-umean[k]) + ci2<TF>*(u[ijk+jj1]-umean[k]) + ci3<TF>*(u[ijk+jj2]-umean[k])))
                                           + cg3<TF>*((ci0<TF>*(u[ijk    ]-umean[k]) + ci1<TF>*(u[ijk+jj1]-umean[k]) + ci2<TF>*(u[ijk+jj2]-umean[k]) + ci3<TF>*(u[ijk+jj3]-umean[k]))) ) * dyi, 2)

                               + std::pow( ( cg0<TF>*((ci0<TF>*(u[ijk-kk3]-umean[k-3]) + ci1<TF>*(u[ijk-kk2]-umean[k-2]) + ci2<TF>*(u[ijk-kk1]-umean[k-1]) + ci3<TF>*(u[ijk    ]-umean[k  ])))
                                           + cg1<TF>*((ci0<TF>*(u[ijk-kk2]-umean[k-2]) + ci1<TF>*(u[ijk-kk1]-umean[k-1]) + ci2<TF>*(u[ijk    ]-umean[k  ]) + ci3<TF>*(u[ijk+kk1]-umean[k+1])))
                                           + cg2<TF>*((ci0<TF>*(u[ijk-kk1]-umean[k-1]) + ci1<TF>*(u[ijk    ]-umean[k  ]) + ci2<TF>*(u[ijk+kk1]-umean[k+1]) + ci3<TF>*(u[ijk+kk2]-umean[k+2])))
                                           + cg3<TF>*((ti0<TF>*(u[ijk-kk1]-umean[k-1]) + ti1<TF>*(u[ijk    ]-umean[k  ]) + ti2<TF>*(u[ijk+kk1]-umean[k+1]) + ti3<TF>*(u[ijk+kk2]-umean[k+2]))) ) * dzi4[k], 2) );

                v2_diss[ijk] = - 2.*visc * (
                                 std::pow( ( cg0<TF>*((ci0<TF>*(v[ijk-ii3]-vmean[k]) + ci1<TF>*(v[ijk-ii2]-vmean[k]) + ci2<TF>*(v[ijk-ii1]-vmean[k]) + ci3<TF>*(v[ijk    ]-vmean[k])))
                                           + cg1<TF>*((ci0<TF>*(v[ijk-ii2]-vmean[k]) + ci1<TF>*(v[ijk-ii1]-vmean[k]) + ci2<TF>*(v[ijk    ]-vmean[k]) + ci3<TF>*(v[ijk+ii1]-vmean[k])))
                                           + cg2<TF>*((ci0<TF>*(v[ijk-ii1]-vmean[k]) + ci1<TF>*(v[ijk    ]-vmean[k]) + ci2<TF>*(v[ijk+ii1]-vmean[k]) + ci3<TF>*(v[ijk+ii2]-vmean[k])))
                                           + cg3<TF>*((ci0<TF>*(v[ijk    ]-vmean[k]) + ci1<TF>*(v[ijk+ii1]-vmean[k]) + ci2<TF>*(v[ijk+ii2]-vmean[k]) + ci3<TF>*(v[ijk+ii3]-vmean[k]))) ) * dxi, 2)

                               + std::pow( ( cg0<TF>*((ci0<TF>*(v[ijk-jj3]-vmean[k]) + ci1<TF>*(v[ijk-jj2]-vmean[k]) + ci2<TF>*(v[ijk-jj1]-vmean[k]) + ci3<TF>*(v[ijk    ]-vmean[k])))
                                           + cg1<TF>*((ci0<TF>*(v[ijk-jj2]-vmean[k]) + ci1<TF>*(v[ijk-jj1]-vmean[k]) + ci2<TF>*(v[ijk    ]-vmean[k]) + ci3<TF>*(v[ijk+jj1]-vmean[k])))
                                           + cg2<TF>*((ci0<TF>*(v[ijk-jj1]-vmean[k]) + ci1<TF>*(v[ijk    ]-vmean[k]) + ci2<TF>*(v[ijk+jj1]-vmean[k]) + ci3<TF>*(v[ijk+jj2]-vmean[k])))
                                           + cg3<TF>*((ci0<TF>*(v[ijk    ]-vmean[k]) + ci1<TF>*(v[ijk+jj1]-vmean[k]) + ci2<TF>*(v[ijk+jj2]-vmean[k]) + ci3<TF>*(v[ijk+jj3]-vmean[k]))) ) * dyi, 2)

                               + std::pow( ( cg0<TF>*((ci0<TF>*(v[ijk-kk3]-vmean[k-3]) + ci1<TF>*(v[ijk-kk2]-vmean[k-2]) + ci2<TF>*(v[ijk-kk1]-vmean[k-1]) + ci3<TF>*(v[ijk    ]-vmean[k  ])))
                                           + cg1<TF>*((ci0<TF>*(v[ijk-kk2]-vmean[k-2]) + ci1<TF>*(v[ijk-kk1]-vmean[k-1]) + ci2<TF>*(v[ijk    ]-vmean[k  ]) + ci3<TF>*(v[ijk+kk1]-vmean[k+1])))
                                           + cg2<TF>*((ci0<TF>*(v[ijk-kk1]-vmean[k-1]) + ci1<TF>*(v[ijk    ]-vmean[k  ]) + ci2<TF>*(v[ijk+kk1]-vmean[k+1]) + ci3<TF>*(v[ijk+kk2]-vmean[k+2])))
                                           + cg3<TF>*((ti0<TF>*(v[ijk-kk1]-vmean[k-1]) + ti1<TF>*(v[ijk    ]-vmean[k  ]) + ti2<TF>*(v[ijk+kk1]-vmean[k+1]) + ti3<TF>*(v[ijk+kk2]-vmean[k+2]))) ) * dzi4[k], 2) );

                tke_diss[ijk] = - visc * (
                                   std::pow( (cg0<TF>*w[ijk-ii1] + cg1<TF>*w[ijk] + cg2<TF>*w[ijk+ii1] + cg3<TF>*w[ijk+ii2]) * dxi, 2)
                                 + std::pow( (cg0<TF>*w[ijk-jj1] + cg1<TF>*w[ijk] + cg2<TF>*w[ijk+jj1] + cg3<TF>*w[ijk+jj2]) * dyi, 2)
                                 + std::pow( (cg0<TF>*w[ijk-kk1] + cg1<TF>*w[ijk] + cg2<TF>*w[ijk+kk1] + cg3<TF>*w[ijk+kk2]) * dzi4[k], 2) );

                tke_diss[ijk] += 0.5*(u2_diss[ijk] + v2_diss[ijk]);
            }

        // calculate the w2 budget term
        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
    #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    w2_diss[ijk] = - 2.*visc * (
                                     std::pow( ( cg0<TF>*(ci0<TF>*w[ijk-ii3] + ci1<TF>*w[ijk-ii2] + ci2<TF>*w[ijk-ii1] + ci3<TF>*w[ijk    ])
                                               + cg1<TF>*(ci0<TF>*w[ijk-ii2] + ci1<TF>*w[ijk-ii1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1])
                                               + cg2<TF>*(ci0<TF>*w[ijk-ii1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+ii1] + ci3<TF>*w[ijk+ii2])
                                               + cg3<TF>*(ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+ii1] + ci2<TF>*w[ijk+ii2] + ci3<TF>*w[ijk+ii3]) ) * dxi, 2)

                                   + std::pow( ( cg0<TF>*(ci0<TF>*w[ijk-jj3] + ci1<TF>*w[ijk-jj2] + ci2<TF>*w[ijk-jj1] + ci3<TF>*w[ijk    ])
                                               + cg1<TF>*(ci0<TF>*w[ijk-jj2] + ci1<TF>*w[ijk-jj1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1])
                                               + cg2<TF>*(ci0<TF>*w[ijk-jj1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+jj1] + ci3<TF>*w[ijk+jj2])
                                               + cg3<TF>*(ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+jj1] + ci2<TF>*w[ijk+jj2] + ci3<TF>*w[ijk+jj3]) ) * dyi, 2)

                                   + std::pow( ( cg0<TF>*(ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ])
                                               + cg1<TF>*(ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1])
                                               + cg2<TF>*(ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2])
                                               + cg3<TF>*(ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3]) ) * dzhi4[k], 2) );
                }

        k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;

                uw_diss[ijk] = - ( ( 2 * visc )


                              * ( ( ( ( cg0<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-ii3-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-ii2-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk-ii1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk    -kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-ii3-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-ii2-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk-ii1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk    -kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-ii3    ] - umean[k  ] ) + ci1<TF>*( u[ijk-ii2    ] - umean[k  ] ) + ci2<TF>*( u[ijk-ii1    ] - umean[k  ] ) + ci3<TF>*( u[ijk        ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-ii3+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-ii2+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk-ii1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk    +kk1] - umean[k+1] ) ) )

                                      + cg1<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-ii2-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-ii1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+ii1-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-ii2-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-ii1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+ii1-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-ii2    ] - umean[k  ] ) + ci1<TF>*( u[ijk-ii1    ] - umean[k  ] ) + ci2<TF>*( u[ijk        ] - umean[k  ] ) + ci3<TF>*( u[ijk+ii1    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-ii2+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-ii1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+ii1+kk1] - umean[k+1] ) ) )

                                      + cg2<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-ii1-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+ii1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+ii2-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-ii1-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+ii1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+ii2-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-ii1    ] - umean[k  ] ) + ci1<TF>*( u[ijk        ] - umean[k  ] ) + ci2<TF>*( u[ijk+ii1    ] - umean[k  ] ) + ci3<TF>*( u[ijk+ii2    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-ii1+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+ii1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+ii2+kk1] - umean[k+1] ) ) )

                                      + cg3<TF>*( ci0<TF>*( ci0<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci1<TF>*( u[ijk+ii1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+ii2-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+ii3-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci1<TF>*( u[ijk+ii1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+ii2-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+ii3-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk        ] - umean[k  ] ) + ci1<TF>*( u[ijk+ii1    ] - umean[k  ] ) + ci2<TF>*( u[ijk+ii2    ] - umean[k  ] ) + ci3<TF>*( u[ijk+ii3    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci1<TF>*( u[ijk+ii1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+ii2+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+ii3+kk1] - umean[k+1] ) ) ) )


                                    * dxi )


                                  * ( cg0<TF>*w[ijk-ii2] + cg1<TF>*w[ijk-ii1] + cg2<TF>*w[ijk    ] + cg3<TF>*w[ijk+ii1] ) )


                                * dxi ) );

                uw_diss[ijk] = - ( ( 2 * visc )

                              * ( ( ( ( cg0<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-jj3-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-jj2-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk-jj1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk    -kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-jj3-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-jj2-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk-jj1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk    -kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-jj3    ] - umean[k  ] ) + ci1<TF>*( u[ijk-jj2    ] - umean[k  ] ) + ci2<TF>*( u[ijk-jj1    ] - umean[k  ] ) + ci3<TF>*( u[ijk        ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-jj3+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-jj2+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk-jj1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk    +kk1] - umean[k+1] ) ) )

                                      + cg1<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-jj2-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-jj1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+jj1-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-jj2-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-jj1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+jj1-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-jj2    ] - umean[k  ] ) + ci1<TF>*( u[ijk-jj1    ] - umean[k  ] ) + ci2<TF>*( u[ijk        ] - umean[k  ] ) + ci3<TF>*( u[ijk+jj1    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-jj2+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-jj1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+jj1+kk1] - umean[k+1] ) ) )

                                      + cg2<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-jj1-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+jj1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+jj2-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-jj1-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+jj1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+jj2-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-jj1    ] - umean[k  ] ) + ci1<TF>*( u[ijk        ] - umean[k  ] ) + ci2<TF>*( u[ijk+jj1    ] - umean[k  ] ) + ci3<TF>*( u[ijk+jj2    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-jj1+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+jj1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+jj2+kk1] - umean[k+1] ) ) )

                                      + cg3<TF>*( ci0<TF>*( ci0<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci1<TF>*( u[ijk+jj1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+jj2-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+jj3-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci1<TF>*( u[ijk+jj1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+jj2-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+jj3-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk        ] - umean[k  ] ) + ci1<TF>*( u[ijk+jj1    ] - umean[k  ] ) + ci2<TF>*( u[ijk+jj2    ] - umean[k  ] ) + ci3<TF>*( u[ijk+jj3    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci1<TF>*( u[ijk+jj1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+jj2+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+jj3+kk1] - umean[k+1] ) ) ) )

                                    * dyi )


                                  * ( cg0<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-jj3] + ci1<TF>*w[ijk-ii1-jj3] + ci2<TF>*w[ijk    -jj3] + ci3<TF>*w[ijk+ii1-jj3] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2-jj2] + ci1<TF>*w[ijk-ii1-jj2] + ci2<TF>*w[ijk    -jj2] + ci3<TF>*w[ijk+ii1-jj2] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2-jj1] + ci1<TF>*w[ijk-ii1-jj1] + ci2<TF>*w[ijk    -jj1] + ci3<TF>*w[ijk+ii1-jj1] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] ) )

                                    + cg1<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-jj2] + ci1<TF>*w[ijk-ii1-jj2] + ci2<TF>*w[ijk    -jj2] + ci3<TF>*w[ijk+ii1-jj2] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2-jj1] + ci1<TF>*w[ijk-ii1-jj1] + ci2<TF>*w[ijk    -jj1] + ci3<TF>*w[ijk+ii1-jj1] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+jj1] + ci1<TF>*w[ijk-ii1+jj1] + ci2<TF>*w[ijk    +jj1] + ci3<TF>*w[ijk+ii1+jj1] ) )

                                    + cg2<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-jj1] + ci1<TF>*w[ijk-ii1-jj1] + ci2<TF>*w[ijk    -jj1] + ci3<TF>*w[ijk+ii1-jj1] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2+jj1] + ci1<TF>*w[ijk-ii1+jj1] + ci2<TF>*w[ijk    +jj1] + ci3<TF>*w[ijk+ii1+jj1] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+jj2] + ci1<TF>*w[ijk-ii1+jj2] + ci2<TF>*w[ijk    +jj2] + ci3<TF>*w[ijk+ii1+jj2] ) )

                                    + cg3<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2+jj1] + ci1<TF>*w[ijk-ii1+jj1] + ci2<TF>*w[ijk    +jj1] + ci3<TF>*w[ijk+ii1+jj1] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2+jj2] + ci1<TF>*w[ijk-ii1+jj2] + ci2<TF>*w[ijk    +jj2] + ci3<TF>*w[ijk+ii1+jj2] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+jj3] + ci1<TF>*w[ijk-ii1+jj3] + ci2<TF>*w[ijk    +jj3] + ci3<TF>*w[ijk+ii1+jj3] ) ) ) )


                                * dyi ) );

                uw_diss[ijk] = - ( ( 2 * visc )


                              * ( ( ( ( cg0<TF>*( u[ijk-kk2] - umean[k-2] ) + cg1<TF>*( u[ijk-kk1] - umean[k-1] ) + cg2<TF>*( u[ijk    ] - umean[k  ] ) + cg3<TF>*( u[ijk+kk1] - umean[k+1] ) ) * dzhi4[k] )


                                  * ( bg0<TF>*( bi0<TF>*( ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk    -kk1] + ci3<TF>*w[ijk+ii1-kk1] )
                                          + bi1<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + bi2<TF>*( ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk    +kk1] + ci3<TF>*w[ijk+ii1+kk1] )
                                          + bi3<TF>*( ci0<TF>*w[ijk-ii2+kk2] + ci1<TF>*w[ijk-ii1+kk2] + ci2<TF>*w[ijk    +kk2] + ci3<TF>*w[ijk+ii1+kk2] ) )

                                    + bg1<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk    -kk1] + ci3<TF>*w[ijk+ii1-kk1] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk    +kk1] + ci3<TF>*w[ijk+ii1+kk1] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+kk2] + ci1<TF>*w[ijk-ii1+kk2] + ci2<TF>*w[ijk    +kk2] + ci3<TF>*w[ijk+ii1+kk2] ) )

                                    + bg2<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk    +kk1] + ci3<TF>*w[ijk+ii1+kk1] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2+kk2] + ci1<TF>*w[ijk-ii1+kk2] + ci2<TF>*w[ijk    +kk2] + ci3<TF>*w[ijk+ii1+kk2] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+kk3] + ci1<TF>*w[ijk-ii1+kk3] + ci2<TF>*w[ijk    +kk3] + ci3<TF>*w[ijk+ii1+kk3] ) )

                                    + bg3<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk    +kk1] + ci3<TF>*w[ijk+ii1+kk1] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2+kk2] + ci1<TF>*w[ijk-ii1+kk2] + ci2<TF>*w[ijk    +kk2] + ci3<TF>*w[ijk+ii1+kk2] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2+kk3] + ci1<TF>*w[ijk-ii1+kk3] + ci2<TF>*w[ijk    +kk3] + ci3<TF>*w[ijk+ii1+kk3] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+kk4] + ci1<TF>*w[ijk-ii1+kk4] + ci2<TF>*w[ijk    +kk4] + ci3<TF>*w[ijk+ii1+kk4] ) ) ) )


                                * dzhi4bot ) );
            }

        k = kstart+1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;

                uw_diss[ijk] = - ( ( 2 * visc )


                              * ( ( ( ( cg0<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-ii3-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-ii2-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk-ii1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk    -kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-ii3-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-ii2-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk-ii1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk    -kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-ii3    ] - umean[k  ] ) + ci1<TF>*( u[ijk-ii2    ] - umean[k  ] ) + ci2<TF>*( u[ijk-ii1    ] - umean[k  ] ) + ci3<TF>*( u[ijk        ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-ii3+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-ii2+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk-ii1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk    +kk1] - umean[k+1] ) ) )

                                      + cg1<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-ii2-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-ii1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+ii1-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-ii2-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-ii1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+ii1-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-ii2    ] - umean[k  ] ) + ci1<TF>*( u[ijk-ii1    ] - umean[k  ] ) + ci2<TF>*( u[ijk        ] - umean[k  ] ) + ci3<TF>*( u[ijk+ii1    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-ii2+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-ii1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+ii1+kk1] - umean[k+1] ) ) )

                                      + cg2<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-ii1-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+ii1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+ii2-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-ii1-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+ii1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+ii2-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-ii1    ] - umean[k  ] ) + ci1<TF>*( u[ijk        ] - umean[k  ] ) + ci2<TF>*( u[ijk+ii1    ] - umean[k  ] ) + ci3<TF>*( u[ijk+ii2    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-ii1+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+ii1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+ii2+kk1] - umean[k+1] ) ) )

                                      + cg3<TF>*( ci0<TF>*( ci0<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci1<TF>*( u[ijk+ii1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+ii2-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+ii3-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci1<TF>*( u[ijk+ii1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+ii2-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+ii3-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk        ] - umean[k  ] ) + ci1<TF>*( u[ijk+ii1    ] - umean[k  ] ) + ci2<TF>*( u[ijk+ii2    ] - umean[k  ] ) + ci3<TF>*( u[ijk+ii3    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci1<TF>*( u[ijk+ii1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+ii2+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+ii3+kk1] - umean[k+1] ) ) ) )


                                    * dxi )


                                  * ( cg0<TF>*w[ijk-ii2] + cg1<TF>*w[ijk-ii1] + cg2<TF>*w[ijk    ] + cg3<TF>*w[ijk+ii1] ) )


                                * dxi ) );

                uw_diss[ijk] = - ( ( 2 * visc )


                              * ( ( ( ( cg0<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-jj3-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-jj2-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk-jj1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk    -kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-jj3-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-jj2-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk-jj1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk    -kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-jj3    ] - umean[k  ] ) + ci1<TF>*( u[ijk-jj2    ] - umean[k  ] ) + ci2<TF>*( u[ijk-jj1    ] - umean[k  ] ) + ci3<TF>*( u[ijk        ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-jj3+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-jj2+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk-jj1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk    +kk1] - umean[k+1] ) ) )

                                      + cg1<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-jj2-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-jj1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+jj1-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-jj2-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-jj1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+jj1-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-jj2    ] - umean[k  ] ) + ci1<TF>*( u[ijk-jj1    ] - umean[k  ] ) + ci2<TF>*( u[ijk        ] - umean[k  ] ) + ci3<TF>*( u[ijk+jj1    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-jj2+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-jj1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+jj1+kk1] - umean[k+1] ) ) )

                                      + cg2<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-jj1-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+jj1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+jj2-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-jj1-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+jj1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+jj2-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-jj1    ] - umean[k  ] ) + ci1<TF>*( u[ijk        ] - umean[k  ] ) + ci2<TF>*( u[ijk+jj1    ] - umean[k  ] ) + ci3<TF>*( u[ijk+jj2    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-jj1+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+jj1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+jj2+kk1] - umean[k+1] ) ) )

                                      + cg3<TF>*( ci0<TF>*( ci0<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci1<TF>*( u[ijk+jj1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+jj2-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+jj3-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci1<TF>*( u[ijk+jj1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+jj2-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+jj3-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk        ] - umean[k  ] ) + ci1<TF>*( u[ijk+jj1    ] - umean[k  ] ) + ci2<TF>*( u[ijk+jj2    ] - umean[k  ] ) + ci3<TF>*( u[ijk+jj3    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci1<TF>*( u[ijk+jj1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+jj2+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+jj3+kk1] - umean[k+1] ) ) ) )


                                    * dyi )


                                  * ( cg0<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-jj3] + ci1<TF>*w[ijk-ii1-jj3] + ci2<TF>*w[ijk    -jj3] + ci3<TF>*w[ijk+ii1-jj3] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2-jj2] + ci1<TF>*w[ijk-ii1-jj2] + ci2<TF>*w[ijk    -jj2] + ci3<TF>*w[ijk+ii1-jj2] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2-jj1] + ci1<TF>*w[ijk-ii1-jj1] + ci2<TF>*w[ijk    -jj1] + ci3<TF>*w[ijk+ii1-jj1] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] ) )

                                    + cg1<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-jj2] + ci1<TF>*w[ijk-ii1-jj2] + ci2<TF>*w[ijk    -jj2] + ci3<TF>*w[ijk+ii1-jj2] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2-jj1] + ci1<TF>*w[ijk-ii1-jj1] + ci2<TF>*w[ijk    -jj1] + ci3<TF>*w[ijk+ii1-jj1] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+jj1] + ci1<TF>*w[ijk-ii1+jj1] + ci2<TF>*w[ijk    +jj1] + ci3<TF>*w[ijk+ii1+jj1] ) )

                                    + cg2<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-jj1] + ci1<TF>*w[ijk-ii1-jj1] + ci2<TF>*w[ijk    -jj1] + ci3<TF>*w[ijk+ii1-jj1] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2+jj1] + ci1<TF>*w[ijk-ii1+jj1] + ci2<TF>*w[ijk    +jj1] + ci3<TF>*w[ijk+ii1+jj1] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+jj2] + ci1<TF>*w[ijk-ii1+jj2] + ci2<TF>*w[ijk    +jj2] + ci3<TF>*w[ijk+ii1+jj2] ) )

                                    + cg3<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2+jj1] + ci1<TF>*w[ijk-ii1+jj1] + ci2<TF>*w[ijk    +jj1] + ci3<TF>*w[ijk+ii1+jj1] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2+jj2] + ci1<TF>*w[ijk-ii1+jj2] + ci2<TF>*w[ijk    +jj2] + ci3<TF>*w[ijk+ii1+jj2] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+jj3] + ci1<TF>*w[ijk-ii1+jj3] + ci2<TF>*w[ijk    +jj3] + ci3<TF>*w[ijk+ii1+jj3] ) ) ) )


                                * dyi ) );

                uw_diss[ijk] = - ( ( 2 * visc )


                              * ( ( ( ( cg0<TF>*( u[ijk-kk2] - umean[k-2] ) + cg1<TF>*( u[ijk-kk1] - umean[k-1] ) + cg2<TF>*( u[ijk    ] - umean[k  ] ) + cg3<TF>*( u[ijk+kk1] - umean[k+1] ) ) * dzhi4[k] )


                                  * ( cg0<TF>*( bi0<TF>*( ci0<TF>*w[ijk-ii2-kk2] + ci1<TF>*w[ijk-ii1-kk2] + ci2<TF>*w[ijk    -kk2] + ci3<TF>*w[ijk+ii1-kk2] )
                                          + bi1<TF>*( ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk    -kk1] + ci3<TF>*w[ijk+ii1-kk1] )
                                          + bi2<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + bi3<TF>*( ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk    +kk1] + ci3<TF>*w[ijk+ii1+kk1] ) )

                                    + cg1<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-kk2] + ci1<TF>*w[ijk-ii1-kk2] + ci2<TF>*w[ijk    -kk2] + ci3<TF>*w[ijk+ii1-kk2] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk    -kk1] + ci3<TF>*w[ijk+ii1-kk1] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk    +kk1] + ci3<TF>*w[ijk+ii1+kk1] ) )

                                    + cg2<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk    -kk1] + ci3<TF>*w[ijk+ii1-kk1] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk    +kk1] + ci3<TF>*w[ijk+ii1+kk1] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+kk2] + ci1<TF>*w[ijk-ii1+kk2] + ci2<TF>*w[ijk    +kk2] + ci3<TF>*w[ijk+ii1+kk2] ) )

                                    + cg3<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk    +kk1] + ci3<TF>*w[ijk+ii1+kk1] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2+kk2] + ci1<TF>*w[ijk-ii1+kk2] + ci2<TF>*w[ijk    +kk2] + ci3<TF>*w[ijk+ii1+kk2] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+kk3] + ci1<TF>*w[ijk-ii1+kk3] + ci2<TF>*w[ijk    +kk3] + ci3<TF>*w[ijk+ii1+kk3] ) ) ) )


                                * dzhi4[k] ) );
            }

        for (int k=kstart+2; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;

                    uw_diss[ijk] = - ( ( 2 * visc )


                                  * ( ( ( ( cg0<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-ii3-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-ii2-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk-ii1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk    -kk2] - umean[k-2] ) )
                                                + ci1<TF>*( ci0<TF>*( u[ijk-ii3-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-ii2-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk-ii1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk    -kk1] - umean[k-1] ) )
                                                + ci2<TF>*( ci0<TF>*( u[ijk-ii3    ] - umean[k  ] ) + ci1<TF>*( u[ijk-ii2    ] - umean[k  ] ) + ci2<TF>*( u[ijk-ii1    ] - umean[k  ] ) + ci3<TF>*( u[ijk        ] - umean[k  ] ) )
                                                + ci3<TF>*( ci0<TF>*( u[ijk-ii3+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-ii2+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk-ii1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk    +kk1] - umean[k+1] ) ) )

                                          + cg1<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-ii2-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-ii1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+ii1-kk2] - umean[k-2] ) )
                                                + ci1<TF>*( ci0<TF>*( u[ijk-ii2-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-ii1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+ii1-kk1] - umean[k-1] ) )
                                                + ci2<TF>*( ci0<TF>*( u[ijk-ii2    ] - umean[k  ] ) + ci1<TF>*( u[ijk-ii1    ] - umean[k  ] ) + ci2<TF>*( u[ijk        ] - umean[k  ] ) + ci3<TF>*( u[ijk+ii1    ] - umean[k  ] ) )
                                                + ci3<TF>*( ci0<TF>*( u[ijk-ii2+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-ii1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+ii1+kk1] - umean[k+1] ) ) )

                                          + cg2<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-ii1-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+ii1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+ii2-kk2] - umean[k-2] ) )
                                                + ci1<TF>*( ci0<TF>*( u[ijk-ii1-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+ii1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+ii2-kk1] - umean[k-1] ) )
                                                + ci2<TF>*( ci0<TF>*( u[ijk-ii1    ] - umean[k  ] ) + ci1<TF>*( u[ijk        ] - umean[k  ] ) + ci2<TF>*( u[ijk+ii1    ] - umean[k  ] ) + ci3<TF>*( u[ijk+ii2    ] - umean[k  ] ) )
                                                + ci3<TF>*( ci0<TF>*( u[ijk-ii1+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+ii1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+ii2+kk1] - umean[k+1] ) ) )

                                          + cg3<TF>*( ci0<TF>*( ci0<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci1<TF>*( u[ijk+ii1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+ii2-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+ii3-kk2] - umean[k-2] ) )
                                                + ci1<TF>*( ci0<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci1<TF>*( u[ijk+ii1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+ii2-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+ii3-kk1] - umean[k-1] ) )
                                                + ci2<TF>*( ci0<TF>*( u[ijk        ] - umean[k  ] ) + ci1<TF>*( u[ijk+ii1    ] - umean[k  ] ) + ci2<TF>*( u[ijk+ii2    ] - umean[k  ] ) + ci3<TF>*( u[ijk+ii3    ] - umean[k  ] ) )
                                                + ci3<TF>*( ci0<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci1<TF>*( u[ijk+ii1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+ii2+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+ii3+kk1] - umean[k+1] ) ) ) )


                                        * dxi )


                                      * ( cg0<TF>*w[ijk-ii2] + cg1<TF>*w[ijk-ii1] + cg2<TF>*w[ijk    ] + cg3<TF>*w[ijk+ii1] ) )


                                    * dxi ) );

                    uw_diss[ijk] = - ( ( 2 * visc )


                                  * ( ( ( ( cg0<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-jj3-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-jj2-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk-jj1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk    -kk2] - umean[k-2] ) )
                                                + ci1<TF>*( ci0<TF>*( u[ijk-jj3-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-jj2-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk-jj1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk    -kk1] - umean[k-1] ) )
                                                + ci2<TF>*( ci0<TF>*( u[ijk-jj3    ] - umean[k  ] ) + ci1<TF>*( u[ijk-jj2    ] - umean[k  ] ) + ci2<TF>*( u[ijk-jj1    ] - umean[k  ] ) + ci3<TF>*( u[ijk        ] - umean[k  ] ) )
                                                + ci3<TF>*( ci0<TF>*( u[ijk-jj3+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-jj2+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk-jj1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk    +kk1] - umean[k+1] ) ) )

                                          + cg1<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-jj2-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-jj1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+jj1-kk2] - umean[k-2] ) )
                                                + ci1<TF>*( ci0<TF>*( u[ijk-jj2-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-jj1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+jj1-kk1] - umean[k-1] ) )
                                                + ci2<TF>*( ci0<TF>*( u[ijk-jj2    ] - umean[k  ] ) + ci1<TF>*( u[ijk-jj1    ] - umean[k  ] ) + ci2<TF>*( u[ijk        ] - umean[k  ] ) + ci3<TF>*( u[ijk+jj1    ] - umean[k  ] ) )
                                                + ci3<TF>*( ci0<TF>*( u[ijk-jj2+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-jj1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+jj1+kk1] - umean[k+1] ) ) )

                                          + cg2<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-jj1-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+jj1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+jj2-kk2] - umean[k-2] ) )
                                                + ci1<TF>*( ci0<TF>*( u[ijk-jj1-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+jj1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+jj2-kk1] - umean[k-1] ) )
                                                + ci2<TF>*( ci0<TF>*( u[ijk-jj1    ] - umean[k  ] ) + ci1<TF>*( u[ijk        ] - umean[k  ] ) + ci2<TF>*( u[ijk+jj1    ] - umean[k  ] ) + ci3<TF>*( u[ijk+jj2    ] - umean[k  ] ) )
                                                + ci3<TF>*( ci0<TF>*( u[ijk-jj1+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+jj1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+jj2+kk1] - umean[k+1] ) ) )

                                          + cg3<TF>*( ci0<TF>*( ci0<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci1<TF>*( u[ijk+jj1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+jj2-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+jj3-kk2] - umean[k-2] ) )
                                                + ci1<TF>*( ci0<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci1<TF>*( u[ijk+jj1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+jj2-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+jj3-kk1] - umean[k-1] ) )
                                                + ci2<TF>*( ci0<TF>*( u[ijk        ] - umean[k  ] ) + ci1<TF>*( u[ijk+jj1    ] - umean[k  ] ) + ci2<TF>*( u[ijk+jj2    ] - umean[k  ] ) + ci3<TF>*( u[ijk+jj3    ] - umean[k  ] ) )
                                                + ci3<TF>*( ci0<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci1<TF>*( u[ijk+jj1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+jj2+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+jj3+kk1] - umean[k+1] ) ) ) )


                                        * dyi )


                                      * ( cg0<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-jj3] + ci1<TF>*w[ijk-ii1-jj3] + ci2<TF>*w[ijk    -jj3] + ci3<TF>*w[ijk+ii1-jj3] )
                                              + ci1<TF>*( ci0<TF>*w[ijk-ii2-jj2] + ci1<TF>*w[ijk-ii1-jj2] + ci2<TF>*w[ijk    -jj2] + ci3<TF>*w[ijk+ii1-jj2] )
                                              + ci2<TF>*( ci0<TF>*w[ijk-ii2-jj1] + ci1<TF>*w[ijk-ii1-jj1] + ci2<TF>*w[ijk    -jj1] + ci3<TF>*w[ijk+ii1-jj1] )
                                              + ci3<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] ) )

                                        + cg1<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-jj2] + ci1<TF>*w[ijk-ii1-jj2] + ci2<TF>*w[ijk    -jj2] + ci3<TF>*w[ijk+ii1-jj2] )
                                              + ci1<TF>*( ci0<TF>*w[ijk-ii2-jj1] + ci1<TF>*w[ijk-ii1-jj1] + ci2<TF>*w[ijk    -jj1] + ci3<TF>*w[ijk+ii1-jj1] )
                                              + ci2<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                              + ci3<TF>*( ci0<TF>*w[ijk-ii2+jj1] + ci1<TF>*w[ijk-ii1+jj1] + ci2<TF>*w[ijk    +jj1] + ci3<TF>*w[ijk+ii1+jj1] ) )

                                        + cg2<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-jj1] + ci1<TF>*w[ijk-ii1-jj1] + ci2<TF>*w[ijk    -jj1] + ci3<TF>*w[ijk+ii1-jj1] )
                                              + ci1<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                              + ci2<TF>*( ci0<TF>*w[ijk-ii2+jj1] + ci1<TF>*w[ijk-ii1+jj1] + ci2<TF>*w[ijk    +jj1] + ci3<TF>*w[ijk+ii1+jj1] )
                                              + ci3<TF>*( ci0<TF>*w[ijk-ii2+jj2] + ci1<TF>*w[ijk-ii1+jj2] + ci2<TF>*w[ijk    +jj2] + ci3<TF>*w[ijk+ii1+jj2] ) )

                                        + cg3<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                              + ci1<TF>*( ci0<TF>*w[ijk-ii2+jj1] + ci1<TF>*w[ijk-ii1+jj1] + ci2<TF>*w[ijk    +jj1] + ci3<TF>*w[ijk+ii1+jj1] )
                                              + ci2<TF>*( ci0<TF>*w[ijk-ii2+jj2] + ci1<TF>*w[ijk-ii1+jj2] + ci2<TF>*w[ijk    +jj2] + ci3<TF>*w[ijk+ii1+jj2] )
                                              + ci3<TF>*( ci0<TF>*w[ijk-ii2+jj3] + ci1<TF>*w[ijk-ii1+jj3] + ci2<TF>*w[ijk    +jj3] + ci3<TF>*w[ijk+ii1+jj3] ) ) ) )


                                    * dyi ) );

                    uw_diss[ijk] = - ( ( 2 * visc )


                                  * ( ( ( ( cg0<TF>*( u[ijk-kk2] - umean[k-2] ) + cg1<TF>*( u[ijk-kk1] - umean[k-1] ) + cg2<TF>*( u[ijk    ] - umean[k  ] ) + cg3<TF>*( u[ijk+kk1] - umean[k+1] ) ) * dzhi4[k] )


                                      * ( cg0<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-kk3] + ci1<TF>*w[ijk-ii1-kk3] + ci2<TF>*w[ijk    -kk3] + ci3<TF>*w[ijk+ii1-kk3] )
                                              + ci1<TF>*( ci0<TF>*w[ijk-ii2-kk2] + ci1<TF>*w[ijk-ii1-kk2] + ci2<TF>*w[ijk    -kk2] + ci3<TF>*w[ijk+ii1-kk2] )
                                              + ci2<TF>*( ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk    -kk1] + ci3<TF>*w[ijk+ii1-kk1] )
                                              + ci3<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] ) )

                                        + cg1<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-kk2] + ci1<TF>*w[ijk-ii1-kk2] + ci2<TF>*w[ijk    -kk2] + ci3<TF>*w[ijk+ii1-kk2] )
                                              + ci1<TF>*( ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk    -kk1] + ci3<TF>*w[ijk+ii1-kk1] )
                                              + ci2<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                              + ci3<TF>*( ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk    +kk1] + ci3<TF>*w[ijk+ii1+kk1] ) )

                                        + cg2<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk    -kk1] + ci3<TF>*w[ijk+ii1-kk1] )
                                              + ci1<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                              + ci2<TF>*( ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk    +kk1] + ci3<TF>*w[ijk+ii1+kk1] )
                                              + ci3<TF>*( ci0<TF>*w[ijk-ii2+kk2] + ci1<TF>*w[ijk-ii1+kk2] + ci2<TF>*w[ijk    +kk2] + ci3<TF>*w[ijk+ii1+kk2] ) )

                                        + cg3<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                              + ci1<TF>*( ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk    +kk1] + ci3<TF>*w[ijk+ii1+kk1] )
                                              + ci2<TF>*( ci0<TF>*w[ijk-ii2+kk2] + ci1<TF>*w[ijk-ii1+kk2] + ci2<TF>*w[ijk    +kk2] + ci3<TF>*w[ijk+ii1+kk2] )
                                              + ci3<TF>*( ci0<TF>*w[ijk-ii2+kk3] + ci1<TF>*w[ijk-ii1+kk3] + ci2<TF>*w[ijk    +kk3] + ci3<TF>*w[ijk+ii1+kk3] ) ) ) )


                                    * dzhi4[k] ) );

                }

        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;

                uw_diss[ijk] = - ( ( 2 * visc )


                              * ( ( ( ( cg0<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-ii3-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-ii2-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk-ii1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk    -kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-ii3-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-ii2-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk-ii1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk    -kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-ii3    ] - umean[k  ] ) + ci1<TF>*( u[ijk-ii2    ] - umean[k  ] ) + ci2<TF>*( u[ijk-ii1    ] - umean[k  ] ) + ci3<TF>*( u[ijk        ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-ii3+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-ii2+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk-ii1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk    +kk1] - umean[k+1] ) ) )

                                      + cg1<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-ii2-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-ii1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+ii1-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-ii2-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-ii1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+ii1-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-ii2    ] - umean[k  ] ) + ci1<TF>*( u[ijk-ii1    ] - umean[k  ] ) + ci2<TF>*( u[ijk        ] - umean[k  ] ) + ci3<TF>*( u[ijk+ii1    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-ii2+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-ii1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+ii1+kk1] - umean[k+1] ) ) )

                                      + cg2<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-ii1-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+ii1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+ii2-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-ii1-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+ii1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+ii2-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-ii1    ] - umean[k  ] ) + ci1<TF>*( u[ijk        ] - umean[k  ] ) + ci2<TF>*( u[ijk+ii1    ] - umean[k  ] ) + ci3<TF>*( u[ijk+ii2    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-ii1+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+ii1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+ii2+kk1] - umean[k+1] ) ) )

                                      + cg3<TF>*( ci0<TF>*( ci0<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci1<TF>*( u[ijk+ii1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+ii2-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+ii3-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci1<TF>*( u[ijk+ii1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+ii2-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+ii3-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk        ] - umean[k  ] ) + ci1<TF>*( u[ijk+ii1    ] - umean[k  ] ) + ci2<TF>*( u[ijk+ii2    ] - umean[k  ] ) + ci3<TF>*( u[ijk+ii3    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci1<TF>*( u[ijk+ii1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+ii2+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+ii3+kk1] - umean[k+1] ) ) ) )


                                    * dxi )


                                  * ( cg0<TF>*w[ijk-ii2] + cg1<TF>*w[ijk-ii1] + cg2<TF>*w[ijk    ] + cg3<TF>*w[ijk+ii1] ) )


                                * dxi ) );

                uw_diss[ijk] = - ( ( 2 * visc )


                              * ( ( ( ( cg0<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-jj3-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-jj2-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk-jj1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk    -kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-jj3-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-jj2-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk-jj1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk    -kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-jj3    ] - umean[k  ] ) + ci1<TF>*( u[ijk-jj2    ] - umean[k  ] ) + ci2<TF>*( u[ijk-jj1    ] - umean[k  ] ) + ci3<TF>*( u[ijk        ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-jj3+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-jj2+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk-jj1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk    +kk1] - umean[k+1] ) ) )

                                      + cg1<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-jj2-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-jj1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+jj1-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-jj2-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-jj1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+jj1-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-jj2    ] - umean[k  ] ) + ci1<TF>*( u[ijk-jj1    ] - umean[k  ] ) + ci2<TF>*( u[ijk        ] - umean[k  ] ) + ci3<TF>*( u[ijk+jj1    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-jj2+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-jj1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+jj1+kk1] - umean[k+1] ) ) )

                                      + cg2<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-jj1-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+jj1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+jj2-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-jj1-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+jj1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+jj2-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-jj1    ] - umean[k  ] ) + ci1<TF>*( u[ijk        ] - umean[k  ] ) + ci2<TF>*( u[ijk+jj1    ] - umean[k  ] ) + ci3<TF>*( u[ijk+jj2    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-jj1+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+jj1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+jj2+kk1] - umean[k+1] ) ) )

                                      + cg3<TF>*( ci0<TF>*( ci0<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci1<TF>*( u[ijk+jj1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+jj2-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+jj3-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci1<TF>*( u[ijk+jj1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+jj2-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+jj3-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk        ] - umean[k  ] ) + ci1<TF>*( u[ijk+jj1    ] - umean[k  ] ) + ci2<TF>*( u[ijk+jj2    ] - umean[k  ] ) + ci3<TF>*( u[ijk+jj3    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci1<TF>*( u[ijk+jj1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+jj2+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+jj3+kk1] - umean[k+1] ) ) ) )


                                    * dyi )


                                  * ( cg0<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-jj3] + ci1<TF>*w[ijk-ii1-jj3] + ci2<TF>*w[ijk    -jj3] + ci3<TF>*w[ijk+ii1-jj3] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2-jj2] + ci1<TF>*w[ijk-ii1-jj2] + ci2<TF>*w[ijk    -jj2] + ci3<TF>*w[ijk+ii1-jj2] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2-jj1] + ci1<TF>*w[ijk-ii1-jj1] + ci2<TF>*w[ijk    -jj1] + ci3<TF>*w[ijk+ii1-jj1] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] ) )

                                    + cg1<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-jj2] + ci1<TF>*w[ijk-ii1-jj2] + ci2<TF>*w[ijk    -jj2] + ci3<TF>*w[ijk+ii1-jj2] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2-jj1] + ci1<TF>*w[ijk-ii1-jj1] + ci2<TF>*w[ijk    -jj1] + ci3<TF>*w[ijk+ii1-jj1] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+jj1] + ci1<TF>*w[ijk-ii1+jj1] + ci2<TF>*w[ijk    +jj1] + ci3<TF>*w[ijk+ii1+jj1] ) )

                                    + cg2<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-jj1] + ci1<TF>*w[ijk-ii1-jj1] + ci2<TF>*w[ijk    -jj1] + ci3<TF>*w[ijk+ii1-jj1] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2+jj1] + ci1<TF>*w[ijk-ii1+jj1] + ci2<TF>*w[ijk    +jj1] + ci3<TF>*w[ijk+ii1+jj1] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+jj2] + ci1<TF>*w[ijk-ii1+jj2] + ci2<TF>*w[ijk    +jj2] + ci3<TF>*w[ijk+ii1+jj2] ) )

                                    + cg3<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2+jj1] + ci1<TF>*w[ijk-ii1+jj1] + ci2<TF>*w[ijk    +jj1] + ci3<TF>*w[ijk+ii1+jj1] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2+jj2] + ci1<TF>*w[ijk-ii1+jj2] + ci2<TF>*w[ijk    +jj2] + ci3<TF>*w[ijk+ii1+jj2] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+jj3] + ci1<TF>*w[ijk-ii1+jj3] + ci2<TF>*w[ijk    +jj3] + ci3<TF>*w[ijk+ii1+jj3] ) ) ) )


                                * dyi ) );

                uw_diss[ijk] = - ( ( 2 * visc )


                              * ( ( ( ( cg0<TF>*( u[ijk-kk2] - umean[k-2] ) + cg1<TF>*( u[ijk-kk1] - umean[k-1] ) + cg2<TF>*( u[ijk    ] - umean[k  ] ) + cg3<TF>*( u[ijk+kk1] - umean[k+1] ) ) * dzhi4[k] )


                                  * ( cg0<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-kk3] + ci1<TF>*w[ijk-ii1-kk3] + ci2<TF>*w[ijk    -kk3] + ci3<TF>*w[ijk+ii1-kk3] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2-kk2] + ci1<TF>*w[ijk-ii1-kk2] + ci2<TF>*w[ijk    -kk2] + ci3<TF>*w[ijk+ii1-kk2] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk    -kk1] + ci3<TF>*w[ijk+ii1-kk1] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] ) )

                                    + cg1<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-kk2] + ci1<TF>*w[ijk-ii1-kk2] + ci2<TF>*w[ijk    -kk2] + ci3<TF>*w[ijk+ii1-kk2] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk    -kk1] + ci3<TF>*w[ijk+ii1-kk1] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk    +kk1] + ci3<TF>*w[ijk+ii1+kk1] ) )

                                    + cg2<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk    -kk1] + ci3<TF>*w[ijk+ii1-kk1] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk    +kk1] + ci3<TF>*w[ijk+ii1+kk1] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+kk2] + ci1<TF>*w[ijk-ii1+kk2] + ci2<TF>*w[ijk    +kk2] + ci3<TF>*w[ijk+ii1+kk2] ) )

                                    + cg3<TF>*( ti0<TF>*( ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk    -kk1] + ci3<TF>*w[ijk+ii1-kk1] )
                                          + ti1<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ti2<TF>*( ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk    +kk1] + ci3<TF>*w[ijk+ii1+kk1] )
                                          + ti3<TF>*( ci0<TF>*w[ijk-ii2+kk2] + ci1<TF>*w[ijk-ii1+kk2] + ci2<TF>*w[ijk    +kk2] + ci3<TF>*w[ijk+ii1+kk2] ) ) ) )


                                * dzhi4[k] ) );
            }

        k = kend;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;

                uw_diss[ijk] = - ( ( 2 * visc )


                              * ( ( ( ( cg0<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-ii3-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-ii2-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk-ii1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk    -kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-ii3-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-ii2-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk-ii1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk    -kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-ii3    ] - umean[k  ] ) + ci1<TF>*( u[ijk-ii2    ] - umean[k  ] ) + ci2<TF>*( u[ijk-ii1    ] - umean[k  ] ) + ci3<TF>*( u[ijk        ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-ii3+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-ii2+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk-ii1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk    +kk1] - umean[k+1] ) ) )

                                      + cg1<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-ii2-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-ii1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+ii1-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-ii2-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-ii1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+ii1-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-ii2    ] - umean[k  ] ) + ci1<TF>*( u[ijk-ii1    ] - umean[k  ] ) + ci2<TF>*( u[ijk        ] - umean[k  ] ) + ci3<TF>*( u[ijk+ii1    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-ii2+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-ii1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+ii1+kk1] - umean[k+1] ) ) )

                                      + cg2<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-ii1-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+ii1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+ii2-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-ii1-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+ii1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+ii2-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-ii1    ] - umean[k  ] ) + ci1<TF>*( u[ijk        ] - umean[k  ] ) + ci2<TF>*( u[ijk+ii1    ] - umean[k  ] ) + ci3<TF>*( u[ijk+ii2    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-ii1+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+ii1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+ii2+kk1] - umean[k+1] ) ) )

                                      + cg3<TF>*( ci0<TF>*( ci0<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci1<TF>*( u[ijk+ii1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+ii2-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+ii3-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci1<TF>*( u[ijk+ii1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+ii2-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+ii3-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk        ] - umean[k  ] ) + ci1<TF>*( u[ijk+ii1    ] - umean[k  ] ) + ci2<TF>*( u[ijk+ii2    ] - umean[k  ] ) + ci3<TF>*( u[ijk+ii3    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci1<TF>*( u[ijk+ii1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+ii2+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+ii3+kk1] - umean[k+1] ) ) ) )


                                    * dxi )


                                  * ( cg0<TF>*w[ijk-ii2] + cg1<TF>*w[ijk-ii1] + cg2<TF>*w[ijk    ] + cg3<TF>*w[ijk+ii1] ) )


                                * dxi ) );

                uw_diss[ijk] = - ( ( 2 * visc )


                              * ( ( ( ( cg0<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-jj3-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-jj2-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk-jj1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk    -kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-jj3-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-jj2-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk-jj1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk    -kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-jj3    ] - umean[k  ] ) + ci1<TF>*( u[ijk-jj2    ] - umean[k  ] ) + ci2<TF>*( u[ijk-jj1    ] - umean[k  ] ) + ci3<TF>*( u[ijk        ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-jj3+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-jj2+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk-jj1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk    +kk1] - umean[k+1] ) ) )

                                      + cg1<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-jj2-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-jj1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+jj1-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-jj2-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk-jj1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+jj1-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-jj2    ] - umean[k  ] ) + ci1<TF>*( u[ijk-jj1    ] - umean[k  ] ) + ci2<TF>*( u[ijk        ] - umean[k  ] ) + ci3<TF>*( u[ijk+jj1    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-jj2+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk-jj1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+jj1+kk1] - umean[k+1] ) ) )

                                      + cg2<TF>*( ci0<TF>*( ci0<TF>*( u[ijk-jj1-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+jj1-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+jj2-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk-jj1-kk1] - umean[k-1] ) + ci1<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+jj1-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+jj2-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk-jj1    ] - umean[k  ] ) + ci1<TF>*( u[ijk        ] - umean[k  ] ) + ci2<TF>*( u[ijk+jj1    ] - umean[k  ] ) + ci3<TF>*( u[ijk+jj2    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk-jj1+kk1] - umean[k+1] ) + ci1<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+jj1+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+jj2+kk1] - umean[k+1] ) ) )

                                      + cg3<TF>*( ci0<TF>*( ci0<TF>*( u[ijk    -kk2] - umean[k-2] ) + ci1<TF>*( u[ijk+jj1-kk2] - umean[k-2] ) + ci2<TF>*( u[ijk+jj2-kk2] - umean[k-2] ) + ci3<TF>*( u[ijk+jj3-kk2] - umean[k-2] ) )
                                            + ci1<TF>*( ci0<TF>*( u[ijk    -kk1] - umean[k-1] ) + ci1<TF>*( u[ijk+jj1-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk+jj2-kk1] - umean[k-1] ) + ci3<TF>*( u[ijk+jj3-kk1] - umean[k-1] ) )
                                            + ci2<TF>*( ci0<TF>*( u[ijk        ] - umean[k  ] ) + ci1<TF>*( u[ijk+jj1    ] - umean[k  ] ) + ci2<TF>*( u[ijk+jj2    ] - umean[k  ] ) + ci3<TF>*( u[ijk+jj3    ] - umean[k  ] ) )
                                            + ci3<TF>*( ci0<TF>*( u[ijk    +kk1] - umean[k+1] ) + ci1<TF>*( u[ijk+jj1+kk1] - umean[k+1] ) + ci2<TF>*( u[ijk+jj2+kk1] - umean[k+1] ) + ci3<TF>*( u[ijk+jj3+kk1] - umean[k+1] ) ) ) )


                                    * dyi )


                                  * ( cg0<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-jj3] + ci1<TF>*w[ijk-ii1-jj3] + ci2<TF>*w[ijk    -jj3] + ci3<TF>*w[ijk+ii1-jj3] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2-jj2] + ci1<TF>*w[ijk-ii1-jj2] + ci2<TF>*w[ijk    -jj2] + ci3<TF>*w[ijk+ii1-jj2] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2-jj1] + ci1<TF>*w[ijk-ii1-jj1] + ci2<TF>*w[ijk    -jj1] + ci3<TF>*w[ijk+ii1-jj1] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] ) )

                                    + cg1<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-jj2] + ci1<TF>*w[ijk-ii1-jj2] + ci2<TF>*w[ijk    -jj2] + ci3<TF>*w[ijk+ii1-jj2] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2-jj1] + ci1<TF>*w[ijk-ii1-jj1] + ci2<TF>*w[ijk    -jj1] + ci3<TF>*w[ijk+ii1-jj1] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+jj1] + ci1<TF>*w[ijk-ii1+jj1] + ci2<TF>*w[ijk    +jj1] + ci3<TF>*w[ijk+ii1+jj1] ) )

                                    + cg2<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-jj1] + ci1<TF>*w[ijk-ii1-jj1] + ci2<TF>*w[ijk    -jj1] + ci3<TF>*w[ijk+ii1-jj1] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2+jj1] + ci1<TF>*w[ijk-ii1+jj1] + ci2<TF>*w[ijk    +jj1] + ci3<TF>*w[ijk+ii1+jj1] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+jj2] + ci1<TF>*w[ijk-ii1+jj2] + ci2<TF>*w[ijk    +jj2] + ci3<TF>*w[ijk+ii1+jj2] ) )

                                    + cg3<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2+jj1] + ci1<TF>*w[ijk-ii1+jj1] + ci2<TF>*w[ijk    +jj1] + ci3<TF>*w[ijk+ii1+jj1] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2+jj2] + ci1<TF>*w[ijk-ii1+jj2] + ci2<TF>*w[ijk    +jj2] + ci3<TF>*w[ijk+ii1+jj2] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+jj3] + ci1<TF>*w[ijk-ii1+jj3] + ci2<TF>*w[ijk    +jj3] + ci3<TF>*w[ijk+ii1+jj3] ) ) ) )


                                * dyi ) );

                uw_diss[ijk] = - ( ( 2 * visc )


                              * ( ( ( ( cg0<TF>*( u[ijk-kk2] - umean[k-2] ) + cg1<TF>*( u[ijk-kk1] - umean[k-1] ) + cg2<TF>*( u[ijk    ] - umean[k  ] ) + cg3<TF>*( u[ijk+kk1] - umean[k+1] ) ) * dzhi4[k] )


                                  * ( tg0<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-kk4] + ci1<TF>*w[ijk-ii1-kk4] + ci2<TF>*w[ijk    -kk4] + ci3<TF>*w[ijk+ii1-kk4] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2-kk3] + ci1<TF>*w[ijk-ii1-kk3] + ci2<TF>*w[ijk    -kk3] + ci3<TF>*w[ijk+ii1-kk3] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2-kk2] + ci1<TF>*w[ijk-ii1-kk2] + ci2<TF>*w[ijk    -kk2] + ci3<TF>*w[ijk+ii1-kk2] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk    -kk1] + ci3<TF>*w[ijk+ii1-kk1] ) )

                                    + tg1<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-kk3] + ci1<TF>*w[ijk-ii1-kk3] + ci2<TF>*w[ijk    -kk3] + ci3<TF>*w[ijk+ii1-kk3] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2-kk2] + ci1<TF>*w[ijk-ii1-kk2] + ci2<TF>*w[ijk    -kk2] + ci3<TF>*w[ijk+ii1-kk2] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk    -kk1] + ci3<TF>*w[ijk+ii1-kk1] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] ) )

                                    + tg2<TF>*( ci0<TF>*( ci0<TF>*w[ijk-ii2-kk2] + ci1<TF>*w[ijk-ii1-kk2] + ci2<TF>*w[ijk    -kk2] + ci3<TF>*w[ijk+ii1-kk2] )
                                          + ci1<TF>*( ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk    -kk1] + ci3<TF>*w[ijk+ii1-kk1] )
                                          + ci2<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ci3<TF>*( ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk    +kk1] + ci3<TF>*w[ijk+ii1+kk1] ) )

                                    + tg3<TF>*( ti0<TF>*( ci0<TF>*w[ijk-ii2-kk2] + ci1<TF>*w[ijk-ii1-kk2] + ci2<TF>*w[ijk    -kk2] + ci3<TF>*w[ijk+ii1-kk2] )
                                          + ti1<TF>*( ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk    -kk1] + ci3<TF>*w[ijk+ii1-kk1] )
                                          + ti2<TF>*( ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk        ] + ci3<TF>*w[ijk+ii1    ] )
                                          + ti3<TF>*( ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk    +kk1] + ci3<TF>*w[ijk+ii1+kk1] ) ) ) )


                                * dzhi4top ) );
            }
    }

    template<typename TF>
    void calc_tke_budget_rdstr(
            TF* restrict u2_rdstr, TF* restrict v2_rdstr, TF* restrict w2_rdstr, TF* restrict uw_rdstr,
            const TF* restrict u, const TF* restrict v, const TF* restrict w, const TF* restrict p,
            const TF* restrict umean, const TF* restrict vmean,
            const TF* restrict dzi4, const TF* restrict dzhi4,
            const TF dxi, const TF dyi,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        using namespace Finite_difference::O4;

        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj1 = 1*icells;
        const int jj2 = 2*icells;
        const int jj3 = 3*icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;
        const int kk3 = 3*ijcells;

        // CALCULATE THE PRESSURE REDISTRIBUTION TERM
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    u2_rdstr[ijk] = 2.*(ci0<TF>*p[ijk-ii2] + ci1<TF>*p[ijk-ii1] + ci2<TF>*p[ijk] + ci3<TF>*p[ijk+ii1])
                                  * ( cg0<TF>*((ci0<TF>*(u[ijk-ii3]-umean[k]) + ci1<TF>*(u[ijk-ii2]-umean[k]) + ci2<TF>*(u[ijk-ii1]-umean[k]) + ci3<TF>*(u[ijk    ]-umean[k])))
                                    + cg1<TF>*((ci0<TF>*(u[ijk-ii2]-umean[k]) + ci1<TF>*(u[ijk-ii1]-umean[k]) + ci2<TF>*(u[ijk    ]-umean[k]) + ci3<TF>*(u[ijk+ii1]-umean[k])))
                                    + cg2<TF>*((ci0<TF>*(u[ijk-ii1]-umean[k]) + ci1<TF>*(u[ijk    ]-umean[k]) + ci2<TF>*(u[ijk+ii1]-umean[k]) + ci3<TF>*(u[ijk+ii2]-umean[k])))
                                    + cg3<TF>*((ci0<TF>*(u[ijk    ]-umean[k]) + ci1<TF>*(u[ijk+ii1]-umean[k]) + ci2<TF>*(u[ijk+ii2]-umean[k]) + ci3<TF>*(u[ijk+ii3]-umean[k]))) ) * dxi;
                    v2_rdstr[ijk] = 2.*(ci0<TF>*p[ijk-jj2] + ci1<TF>*p[ijk-jj1] + ci2<TF>*p[ijk] + ci3<TF>*p[ijk+jj1])
                                  * ( cg0<TF>*((ci0<TF>*(v[ijk-jj3]-vmean[k]) + ci1<TF>*(v[ijk-jj2]-vmean[k]) + ci2<TF>*(v[ijk-jj1]-vmean[k]) + ci3<TF>*(v[ijk    ]-vmean[k])))
                                    + cg1<TF>*((ci0<TF>*(v[ijk-jj2]-vmean[k]) + ci1<TF>*(v[ijk-jj1]-vmean[k]) + ci2<TF>*(v[ijk    ]-vmean[k]) + ci3<TF>*(v[ijk+jj1]-vmean[k])))
                                    + cg2<TF>*((ci0<TF>*(v[ijk-jj1]-vmean[k]) + ci1<TF>*(v[ijk    ]-vmean[k]) + ci2<TF>*(v[ijk+jj1]-vmean[k]) + ci3<TF>*(v[ijk+jj2]-vmean[k])))
                                    + cg3<TF>*((ci0<TF>*(v[ijk    ]-vmean[k]) + ci1<TF>*(v[ijk+jj1]-vmean[k]) + ci2<TF>*(v[ijk+jj2]-vmean[k]) + ci3<TF>*(v[ijk+jj3]-vmean[k]))) ) * dyi;
                }

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    w2_rdstr[ijk] = 2.*(ci0<TF>*p[ijk-kk2] + ci1<TF>*p[ijk-kk1] + ci2<TF>*p[ijk] + ci3<TF>*p[ijk+kk1])
                                  * ( cg0<TF>*(ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ])
                                    + cg1<TF>*(ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1])
                                    + cg2<TF>*(ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2])
                                    + cg3<TF>*(ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3]) ) * dzhi4[k];
                }

        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
    #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    uw_rdstr[ijk] = ( ( ci0<TF>*( ci0<TF>*p[ijk-ii2-kk2] + ci1<TF>*p[ijk-ii1-kk2] + ci2<TF>*p[ijk-kk2] + ci3<TF>*p[ijk+ii1-kk2] )
                                      + ci1<TF>*( ci0<TF>*p[ijk-ii2-kk1] + ci1<TF>*p[ijk-ii1-kk1] + ci2<TF>*p[ijk-kk1] + ci3<TF>*p[ijk+ii1-kk1] )
                                      + ci2<TF>*( ci0<TF>*p[ijk-ii2    ] + ci1<TF>*p[ijk-ii1    ] + ci2<TF>*p[ijk    ] + ci3<TF>*p[ijk+ii1    ] )
                                      + ci3<TF>*( ci0<TF>*p[ijk-ii2+kk1] + ci1<TF>*p[ijk-ii1+kk1] + ci2<TF>*p[ijk+kk1] + ci3<TF>*p[ijk+ii1+kk1] ) )

                                    * ( ( ( cg0<TF>*( u[ijk-kk2] - umean[k-2] ) + cg1<TF>*( u[ijk-kk1] - umean[k-1] ) + cg2<TF>*( u[ijk] - umean[k] ) + cg3<TF>*( u[ijk+kk1] - umean[k+1] ) ) * dzhi4[k] ) + ( ( cg0<TF>*w[ijk-ii2] + cg1<TF>*w[ijk-ii1] + cg2<TF>*w[ijk] + cg3<TF>*w[ijk+ii1] ) * dxi ) ) );
                }
    }

    template<typename TF>
    void calc_tke_budget_buoy(
            TF* restrict w2_buoy, TF* restrict tke_buoy, TF* restrict uw_buoy,
            const TF* restrict u, const TF* restrict w, const TF* restrict b,
            const TF* restrict umean, const TF* restrict bmean,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        using namespace Finite_difference::O4;

        const int ii1 = 1;
        const int ii2 = 2;
        const int jj1 = 1*icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        // CALCULATE THE BUOYANCY TERM.
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    tke_buoy[ijk] = (ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2])*( b[ijk] - bmean[k] );
                }

        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    w2_buoy[ijk] = 2.*(ci0<TF>*b[ijk-kk2] + ci1<TF>*b[ijk-kk1] + ci2<TF>*b[ijk] + ci3<TF>*b[ijk+kk1])*w[ijk];

                    uw_buoy[ijk] = ( ( ci0<TF>*( u[ijk-kk2] - umean[k-2] ) + ci1<TF>*( u[ijk-kk1] - umean[k-1] ) + ci2<TF>*( u[ijk    ] - umean[k  ] ) + ci3<TF>*( u[ijk+kk1] - umean[k+1] ) )

                                   * ( ci0<TF>*( ci0<TF>*( b[ijk-ii2-kk2] - bmean[k-2] ) + ci1<TF>*( b[ijk-ii1-kk2] - bmean[k-2] ) + ci2<TF>*( b[ijk    -kk2] - bmean[k-2] ) + ci3<TF>*( b[ijk+ii1-kk2] - bmean[k-2] ) )
                                     + ci1<TF>*( ci0<TF>*( b[ijk-ii2-kk1] - bmean[k-1] ) + ci1<TF>*( b[ijk-ii1-kk1] - bmean[k-1] ) + ci2<TF>*( b[ijk    -kk1] - bmean[k-1] ) + ci3<TF>*( b[ijk+ii1-kk1] - bmean[k-1] ) )
                                     + ci2<TF>*( ci0<TF>*( b[ijk-ii2    ] - bmean[k  ] ) + ci1<TF>*( b[ijk-ii1    ] - bmean[k  ] ) + ci2<TF>*( b[ijk        ] - bmean[k  ] ) + ci3<TF>*( b[ijk+ii1    ] - bmean[k  ] ) )
                                     + ci3<TF>*( ci0<TF>*( b[ijk-ii2+kk1] - bmean[k+1] ) + ci1<TF>*( b[ijk-ii1+kk1] - bmean[k+1] ) + ci2<TF>*( b[ijk    +kk1] - bmean[k+1] ) + ci3<TF>*( b[ijk+ii1+kk1] - bmean[k+1] ) ) ) );
                }
    }

    template<typename TF>
    void calc_b2_budget(
            TF* restrict b2_shear, TF* restrict b2_turb, TF* restrict b2_visc, TF* restrict b2_diss,
            const TF* restrict w, const TF* restrict b,
            const TF* restrict bmean,
            const TF* restrict dzi4, const TF* restrict dzhi4,
            const TF dxi, const TF dyi,
            const TF visc,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        using namespace Finite_difference::O4;

        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj1 = 1*icells;
        const int jj2 = 2*icells;
        const int jj3 = 3*icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;
        const int kk3 = 3*ijcells;

        // 1. CALCULATE THE GRADIENT PRODUCTION TERM
        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                b2_shear[ijk] = - ( ( ( ( 2.0 * ( b[ijk] - bmean[k] ) ) * ( ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2] ) )

                                 * ( cg0<TF>*( bi0<TF>*bmean[k-2] + bi1<TF>*bmean[k-1] + bi2<TF>*bmean[k  ] + bi3<TF>*bmean[k+1] )
                                   + cg1<TF>*( ci0<TF>*bmean[k-2] + ci1<TF>*bmean[k-1] + ci2<TF>*bmean[k  ] + ci3<TF>*bmean[k+1] )
                                   + cg2<TF>*( ci0<TF>*bmean[k-1] + ci1<TF>*bmean[k  ] + ci2<TF>*bmean[k+1] + ci3<TF>*bmean[k+2] )
                                   + cg3<TF>*( ci0<TF>*bmean[k  ] + ci1<TF>*bmean[k+1] + ci2<TF>*bmean[k+2] + ci3<TF>*bmean[k+3] ) ) )

                               * dzi4[k] );
            }

        for (int k=kstart+1; k<kend-1; ++k)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    b2_shear[ijk] = - ( ( ( ( 2.0 * ( b[ijk] - bmean[k] ) ) * ( ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2] ) )

                                     * ( cg0<TF>*( ci0<TF>*bmean[k-3] + ci1<TF>*bmean[k-2] + ci2<TF>*bmean[k-1] + ci3<TF>*bmean[k  ] )
                                       + cg1<TF>*( ci0<TF>*bmean[k-2] + ci1<TF>*bmean[k-1] + ci2<TF>*bmean[k  ] + ci3<TF>*bmean[k+1] )
                                       + cg2<TF>*( ci0<TF>*bmean[k-1] + ci1<TF>*bmean[k  ] + ci2<TF>*bmean[k+1] + ci3<TF>*bmean[k+2] )
                                       + cg3<TF>*( ci0<TF>*bmean[k  ] + ci1<TF>*bmean[k+1] + ci2<TF>*bmean[k+2] + ci3<TF>*bmean[k+3] ) ) )

                                   * dzi4[k] );
                }
        }

        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                b2_shear[ijk] = - ( ( ( ( 2.0 * ( b[ijk] - bmean[k] ) ) * ( ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2] ) )

                                 * ( cg0<TF>*( ci0<TF>*bmean[k-3] + ci1<TF>*bmean[k-2] + ci2<TF>*bmean[k-1] + ci3<TF>*bmean[k  ] )
                                   + cg1<TF>*( ci0<TF>*bmean[k-2] + ci1<TF>*bmean[k-1] + ci2<TF>*bmean[k  ] + ci3<TF>*bmean[k+1] )
                                   + cg2<TF>*( ci0<TF>*bmean[k-1] + ci1<TF>*bmean[k  ] + ci2<TF>*bmean[k+1] + ci3<TF>*bmean[k+2] )
                                   + cg3<TF>*( ti0<TF>*bmean[k-1] + ti1<TF>*bmean[k  ] + ti2<TF>*bmean[k+1] + ti3<TF>*bmean[k+2] ) ) )

                               * dzi4[k] );
            }


        // 2. CALCULATE THE TURBULENT TRANSPORT TERM
        k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                b2_turb[ijk] = - ( ( cg0<TF>*( std::pow( ( bi0<TF>*( b[ijk-kk2] - bmean[k-2] ) + bi1<TF>*( b[ijk-kk1] - bmean[k-1] ) + bi2<TF>*( b[ijk    ] - bmean[k  ] ) + bi3<TF>*( b[ijk+kk1] - bmean[k+1] ) ), 2 ) * w[ijk-kk1] )
                                   + cg1<TF>*( std::pow( ( ci0<TF>*( b[ijk-kk2] - bmean[k-2] ) + ci1<TF>*( b[ijk-kk1] - bmean[k-1] ) + ci2<TF>*( b[ijk    ] - bmean[k  ] ) + ci3<TF>*( b[ijk+kk1] - bmean[k+1] ) ), 2 ) * w[ijk    ] )
                                   + cg2<TF>*( std::pow( ( ci0<TF>*( b[ijk-kk1] - bmean[k-1] ) + ci1<TF>*( b[ijk    ] - bmean[k  ] ) + ci2<TF>*( b[ijk+kk1] - bmean[k+1] ) + ci3<TF>*( b[ijk+kk2] - bmean[k+2] ) ), 2 ) * w[ijk+kk1] )
                                   + cg3<TF>*( std::pow( ( ci0<TF>*( b[ijk    ] - bmean[k  ] ) + ci1<TF>*( b[ijk+kk1] - bmean[k+1] ) + ci2<TF>*( b[ijk+kk2] - bmean[k+2] ) + ci3<TF>*( b[ijk+kk3] - bmean[k+3] ) ), 2 ) * w[ijk+kk2] ) )

                                 * dzi4[k] );
            }

        for (int k=kstart+1; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    b2_turb[ijk] = - ( ( cg0<TF>*( std::pow( ( ci0<TF>*( b[ijk-kk3] - bmean[k-3] ) + ci1<TF>*( b[ijk-kk2] - bmean[k-2] ) + ci2<TF>*( b[ijk-kk1] - bmean[k-1] ) + ci3<TF>*( b[ijk    ] - bmean[k  ] ) ), 2 ) * w[ijk-kk1] )
                                       + cg1<TF>*( std::pow( ( ci0<TF>*( b[ijk-kk2] - bmean[k-2] ) + ci1<TF>*( b[ijk-kk1] - bmean[k-1] ) + ci2<TF>*( b[ijk    ] - bmean[k  ] ) + ci3<TF>*( b[ijk+kk1] - bmean[k+1] ) ), 2 ) * w[ijk    ] )
                                       + cg2<TF>*( std::pow( ( ci0<TF>*( b[ijk-kk1] - bmean[k-1] ) + ci1<TF>*( b[ijk    ] - bmean[k  ] ) + ci2<TF>*( b[ijk+kk1] - bmean[k+1] ) + ci3<TF>*( b[ijk+kk2] - bmean[k+2] ) ), 2 ) * w[ijk+kk1] )
                                       + cg3<TF>*( std::pow( ( ci0<TF>*( b[ijk    ] - bmean[k  ] ) + ci1<TF>*( b[ijk+kk1] - bmean[k+1] ) + ci2<TF>*( b[ijk+kk2] - bmean[k+2] ) + ci3<TF>*( b[ijk+kk3] - bmean[k+3] ) ), 2 ) * w[ijk+kk2] ) )

                                     * dzi4[k] );
                }

        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                b2_turb[ijk] = - ( ( cg0<TF>*( std::pow( ( ci0<TF>*( b[ijk-kk3] - bmean[k-3] ) + ci1<TF>*( b[ijk-kk2] - bmean[k-2] ) + ci2<TF>*( b[ijk-kk1] - bmean[k-1] ) + ci3<TF>*( b[ijk    ] - bmean[k  ] ) ), 2 ) * w[ijk-kk1] )
                                   + cg1<TF>*( std::pow( ( ci0<TF>*( b[ijk-kk2] - bmean[k-2] ) + ci1<TF>*( b[ijk-kk1] - bmean[k-1] ) + ci2<TF>*( b[ijk    ] - bmean[k  ] ) + ci3<TF>*( b[ijk+kk1] - bmean[k+1] ) ), 2 ) * w[ijk    ] )
                                   + cg2<TF>*( std::pow( ( ci0<TF>*( b[ijk-kk1] - bmean[k-1] ) + ci1<TF>*( b[ijk    ] - bmean[k  ] ) + ci2<TF>*( b[ijk+kk1] - bmean[k+1] ) + ci3<TF>*( b[ijk+kk2] - bmean[k+2] ) ), 2 ) * w[ijk+kk1] )
                                   + cg3<TF>*( std::pow( ( ti0<TF>*( b[ijk-kk1] - bmean[k-1] ) + ti1<TF>*( b[ijk    ] - bmean[k  ] ) + ti2<TF>*( b[ijk+kk1] - bmean[k+1] ) + ti3<TF>*( b[ijk+kk2] - bmean[k+2] ) ), 2 ) * w[ijk+kk2] ) )

                                 * dzi4[k] );
            }


        // 3. CALCULATE THE VISCOUS TRANSPORT TERM
        k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                b2_visc[ijk] = ( ( visc

                                 * ( cg0<TF>*( ( bg0<TF>*std::pow( ( b[ijk-kk2] - bmean[k-2] ), 2 ) + bg1<TF>*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + bg2<TF>*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + bg3<TF>*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) ) * dzhi4[k-1] )
                                   + cg1<TF>*( ( cg0<TF>*std::pow( ( b[ijk-kk2] - bmean[k-2] ), 2 ) + cg1<TF>*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + cg2<TF>*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + cg3<TF>*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) ) * dzhi4[k  ] )
                                   + cg2<TF>*( ( cg0<TF>*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + cg1<TF>*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + cg2<TF>*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) + cg3<TF>*std::pow( ( b[ijk+kk2] - bmean[k+2] ), 2 ) ) * dzhi4[k+1] )
                                   + cg3<TF>*( ( cg0<TF>*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + cg1<TF>*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) + cg2<TF>*std::pow( ( b[ijk+kk2] - bmean[k+2] ), 2 ) + cg3<TF>*std::pow( ( b[ijk+kk3] - bmean[k+3] ), 2 ) ) * dzhi4[k+2] ) ) )

                               * dzi4[k] );
            }

        for (int k=kstart+1; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    b2_visc[ijk] = ( ( visc

                                     * ( cg0<TF>*( ( cg0<TF>*std::pow( ( b[ijk-kk3] - bmean[k-3] ), 2 ) + cg1<TF>*std::pow( ( b[ijk-kk2] - bmean[k-2] ), 2 ) + cg2<TF>*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + cg3<TF>*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) ) * dzhi4[k-1] )
                                       + cg1<TF>*( ( cg0<TF>*std::pow( ( b[ijk-kk2] - bmean[k-2] ), 2 ) + cg1<TF>*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + cg2<TF>*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + cg3<TF>*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) ) * dzhi4[k  ] )
                                       + cg2<TF>*( ( cg0<TF>*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + cg1<TF>*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + cg2<TF>*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) + cg3<TF>*std::pow( ( b[ijk+kk2] - bmean[k+2] ), 2 ) ) * dzhi4[k+1] )
                                       + cg3<TF>*( ( cg0<TF>*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + cg1<TF>*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) + cg2<TF>*std::pow( ( b[ijk+kk2] - bmean[k+2] ), 2 ) + cg3<TF>*std::pow( ( b[ijk+kk3] - bmean[k+3] ), 2 ) ) * dzhi4[k+2] ) ) )

                                   * dzi4[k] );
                }

        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                b2_visc[ijk] = ( ( visc

                                 * ( cg0<TF>*( ( cg0<TF>*std::pow( ( b[ijk-kk3] - bmean[k-3] ), 2 ) + cg1<TF>*std::pow( ( b[ijk-kk2] - bmean[k-2] ), 2 ) + cg2<TF>*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + cg3<TF>*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) ) * dzhi4[k-1] )
                                   + cg1<TF>*( ( cg0<TF>*std::pow( ( b[ijk-kk2] - bmean[k-2] ), 2 ) + cg1<TF>*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + cg2<TF>*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + cg3<TF>*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) ) * dzhi4[k  ] )
                                   + cg2<TF>*( ( cg0<TF>*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + cg1<TF>*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + cg2<TF>*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) + cg3<TF>*std::pow( ( b[ijk+kk2] - bmean[k+2] ), 2 ) ) * dzhi4[k+1] )
                                   + cg3<TF>*( ( tg0<TF>*std::pow( ( b[ijk-kk1] - bmean[k-1] ), 2 ) + tg1<TF>*std::pow( ( b[ijk    ] - bmean[k  ] ), 2 ) + tg2<TF>*std::pow( ( b[ijk+kk1] - bmean[k+1] ), 2 ) + tg3<TF>*std::pow( ( b[ijk+kk2] - bmean[k+2] ), 2 ) ) * dzhi4[k+2] ) ) )

                               * dzi4[k] );
            }


        // 4. CALCULATE THE DISSIPATION TERM
        k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                b2_diss[ijk] = - ( ( 2.0 * visc )

                              * ( ( std::pow( ( ( cg0<TF>*( ci0<TF>*( b[ijk-ii3] - bmean[k] ) + ci1<TF>*( b[ijk-ii2] - bmean[k] ) + ci2<TF>*( b[ijk-ii1] - bmean[k] ) + ci3<TF>*( b[ijk    ] - bmean[k] ) )
                                                + cg1<TF>*( ci0<TF>*( b[ijk-ii2] - bmean[k] ) + ci1<TF>*( b[ijk-ii1] - bmean[k] ) + ci2<TF>*( b[ijk    ] - bmean[k] ) + ci3<TF>*( b[ijk+ii1] - bmean[k] ) )
                                                + cg2<TF>*( ci0<TF>*( b[ijk-ii1] - bmean[k] ) + ci1<TF>*( b[ijk    ] - bmean[k] ) + ci2<TF>*( b[ijk+ii1] - bmean[k] ) + ci3<TF>*( b[ijk+ii2] - bmean[k] ) )
                                                + cg3<TF>*( ci0<TF>*( b[ijk    ] - bmean[k] ) + ci1<TF>*( b[ijk+ii1] - bmean[k] ) + ci2<TF>*( b[ijk+ii2] - bmean[k] ) + ci3<TF>*( b[ijk+ii3] - bmean[k] ) ) )

                                              * dxi )

                                    , 2 )

                                  + std::pow( ( ( cg0<TF>*( ci0<TF>*( b[ijk-jj3] - bmean[k] ) + ci1<TF>*( b[ijk-jj2] - bmean[k] ) + ci2<TF>*( b[ijk-jj1] - bmean[k] ) + ci3<TF>*( b[ijk    ] - bmean[k] ) )
                                                + cg1<TF>*( ci0<TF>*( b[ijk-jj2] - bmean[k] ) + ci1<TF>*( b[ijk-jj1] - bmean[k] ) + ci2<TF>*( b[ijk    ] - bmean[k] ) + ci3<TF>*( b[ijk+jj1] - bmean[k] ) )
                                                + cg2<TF>*( ci0<TF>*( b[ijk-jj1] - bmean[k] ) + ci1<TF>*( b[ijk    ] - bmean[k] ) + ci2<TF>*( b[ijk+jj1] - bmean[k] ) + ci3<TF>*( b[ijk+jj2] - bmean[k] ) )
                                                + cg3<TF>*( ci0<TF>*( b[ijk    ] - bmean[k] ) + ci1<TF>*( b[ijk+jj1] - bmean[k] ) + ci2<TF>*( b[ijk+jj2] - bmean[k] ) + ci3<TF>*( b[ijk+jj3] - bmean[k] ) ) )

                                              * dyi )

                                    , 2 ) )

                                + std::pow( ( ( cg0<TF>*( bi0<TF>*( b[ijk-kk2] - bmean[k-2] ) + bi1<TF>*( b[ijk-kk1] - bmean[k-1] ) + bi2<TF>*( b[ijk    ] - bmean[k  ] ) + bi3<TF>*( b[ijk+kk1] - bmean[k+1] ) )
                                              + cg1<TF>*( ci0<TF>*( b[ijk-kk2] - bmean[k-2] ) + ci1<TF>*( b[ijk-kk1] - bmean[k-1] ) + ci2<TF>*( b[ijk    ] - bmean[k  ] ) + ci3<TF>*( b[ijk+kk1] - bmean[k+1] ) )
                                              + cg2<TF>*( ci0<TF>*( b[ijk-kk1] - bmean[k-1] ) + ci1<TF>*( b[ijk    ] - bmean[k  ] ) + ci2<TF>*( b[ijk+kk1] - bmean[k+1] ) + ci3<TF>*( b[ijk+kk2] - bmean[k+2] ) )
                                              + cg3<TF>*( ci0<TF>*( b[ijk    ] - bmean[k  ] ) + ci1<TF>*( b[ijk+kk1] - bmean[k+1] ) + ci2<TF>*( b[ijk+kk2] - bmean[k+2] ) + ci3<TF>*( b[ijk+kk3] - bmean[k+3] ) ) )

                                            * dzi4[k] )

                                  , 2 ) ) );
            }

        for (int k=kstart+1; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    b2_diss[ijk] = - ( ( 2.0 * visc )

                                  * ( ( std::pow( ( ( cg0<TF>*( ci0<TF>*( b[ijk-ii3] - bmean[k] ) + ci1<TF>*( b[ijk-ii2] - bmean[k] ) + ci2<TF>*( b[ijk-ii1] - bmean[k] ) + ci3<TF>*( b[ijk    ] - bmean[k] ) )
                                                    + cg1<TF>*( ci0<TF>*( b[ijk-ii2] - bmean[k] ) + ci1<TF>*( b[ijk-ii1] - bmean[k] ) + ci2<TF>*( b[ijk    ] - bmean[k] ) + ci3<TF>*( b[ijk+ii1] - bmean[k] ) )
                                                    + cg2<TF>*( ci0<TF>*( b[ijk-ii1] - bmean[k] ) + ci1<TF>*( b[ijk    ] - bmean[k] ) + ci2<TF>*( b[ijk+ii1] - bmean[k] ) + ci3<TF>*( b[ijk+ii2] - bmean[k] ) )
                                                    + cg3<TF>*( ci0<TF>*( b[ijk    ] - bmean[k] ) + ci1<TF>*( b[ijk+ii1] - bmean[k] ) + ci2<TF>*( b[ijk+ii2] - bmean[k] ) + ci3<TF>*( b[ijk+ii3] - bmean[k] ) ) )

                                                  * dxi )

                                        , 2 )

                                      + std::pow( ( ( cg0<TF>*( ci0<TF>*( b[ijk-jj3] - bmean[k] ) + ci1<TF>*( b[ijk-jj2] - bmean[k] ) + ci2<TF>*( b[ijk-jj1] - bmean[k] ) + ci3<TF>*( b[ijk    ] - bmean[k] ) )
                                                    + cg1<TF>*( ci0<TF>*( b[ijk-jj2] - bmean[k] ) + ci1<TF>*( b[ijk-jj1] - bmean[k] ) + ci2<TF>*( b[ijk    ] - bmean[k] ) + ci3<TF>*( b[ijk+jj1] - bmean[k] ) )
                                                    + cg2<TF>*( ci0<TF>*( b[ijk-jj1] - bmean[k] ) + ci1<TF>*( b[ijk    ] - bmean[k] ) + ci2<TF>*( b[ijk+jj1] - bmean[k] ) + ci3<TF>*( b[ijk+jj2] - bmean[k] ) )
                                                    + cg3<TF>*( ci0<TF>*( b[ijk    ] - bmean[k] ) + ci1<TF>*( b[ijk+jj1] - bmean[k] ) + ci2<TF>*( b[ijk+jj2] - bmean[k] ) + ci3<TF>*( b[ijk+jj3] - bmean[k] ) ) )

                                                  * dyi )

                                        , 2 ) )

                                    + std::pow( ( ( cg0<TF>*( ci0<TF>*( b[ijk-kk3] - bmean[k-3] ) + ci1<TF>*( b[ijk-kk2] - bmean[k-2] ) + ci2<TF>*( b[ijk-kk1] - bmean[k-1] ) + ci3<TF>*( b[ijk    ] - bmean[k  ] ) )
                                                  + cg1<TF>*( ci0<TF>*( b[ijk-kk2] - bmean[k-2] ) + ci1<TF>*( b[ijk-kk1] - bmean[k-1] ) + ci2<TF>*( b[ijk    ] - bmean[k  ] ) + ci3<TF>*( b[ijk+kk1] - bmean[k+1] ) )
                                                  + cg2<TF>*( ci0<TF>*( b[ijk-kk1] - bmean[k-1] ) + ci1<TF>*( b[ijk    ] - bmean[k  ] ) + ci2<TF>*( b[ijk+kk1] - bmean[k+1] ) + ci3<TF>*( b[ijk+kk2] - bmean[k+2] ) )
                                                  + cg3<TF>*( ci0<TF>*( b[ijk    ] - bmean[k  ] ) + ci1<TF>*( b[ijk+kk1] - bmean[k+1] ) + ci2<TF>*( b[ijk+kk2] - bmean[k+2] ) + ci3<TF>*( b[ijk+kk3] - bmean[k+3] ) ) )

                                                * dzi4[k] )

                                      , 2 ) ) );
                }

        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                b2_diss[ijk] = - ( ( 2.0 * visc )

                              * ( ( std::pow( ( ( cg0<TF>*( ci0<TF>*( b[ijk-ii3] - bmean[k] ) + ci1<TF>*( b[ijk-ii2] - bmean[k] ) + ci2<TF>*( b[ijk-ii1] - bmean[k] ) + ci3<TF>*( b[ijk    ] - bmean[k] ) )
                                                + cg1<TF>*( ci0<TF>*( b[ijk-ii2] - bmean[k] ) + ci1<TF>*( b[ijk-ii1] - bmean[k] ) + ci2<TF>*( b[ijk    ] - bmean[k] ) + ci3<TF>*( b[ijk+ii1] - bmean[k] ) )
                                                + cg2<TF>*( ci0<TF>*( b[ijk-ii1] - bmean[k] ) + ci1<TF>*( b[ijk    ] - bmean[k] ) + ci2<TF>*( b[ijk+ii1] - bmean[k] ) + ci3<TF>*( b[ijk+ii2] - bmean[k] ) )
                                                + cg3<TF>*( ci0<TF>*( b[ijk    ] - bmean[k] ) + ci1<TF>*( b[ijk+ii1] - bmean[k] ) + ci2<TF>*( b[ijk+ii2] - bmean[k] ) + ci3<TF>*( b[ijk+ii3] - bmean[k] ) ) )

                                              * dxi )

                                    , 2 )

                                  + std::pow( ( ( cg0<TF>*( ci0<TF>*( b[ijk-jj3] - bmean[k] ) + ci1<TF>*( b[ijk-jj2] - bmean[k] ) + ci2<TF>*( b[ijk-jj1] - bmean[k] ) + ci3<TF>*( b[ijk    ] - bmean[k] ) )
                                                + cg1<TF>*( ci0<TF>*( b[ijk-jj2] - bmean[k] ) + ci1<TF>*( b[ijk-jj1] - bmean[k] ) + ci2<TF>*( b[ijk    ] - bmean[k] ) + ci3<TF>*( b[ijk+jj1] - bmean[k] ) )
                                                + cg2<TF>*( ci0<TF>*( b[ijk-jj1] - bmean[k] ) + ci1<TF>*( b[ijk    ] - bmean[k] ) + ci2<TF>*( b[ijk+jj1] - bmean[k] ) + ci3<TF>*( b[ijk+jj2] - bmean[k] ) )
                                                + cg3<TF>*( ci0<TF>*( b[ijk    ] - bmean[k] ) + ci1<TF>*( b[ijk+jj1] - bmean[k] ) + ci2<TF>*( b[ijk+jj2] - bmean[k] ) + ci3<TF>*( b[ijk+jj3] - bmean[k] ) ) )

                                              * dyi )

                                    , 2 ) )

                                + std::pow( ( ( cg0<TF>*( ci0<TF>*( b[ijk-kk3] - bmean[k-3] ) + ci1<TF>*( b[ijk-kk2] - bmean[k-2] ) + ci2<TF>*( b[ijk-kk1] - bmean[k-1] ) + ci3<TF>*( b[ijk    ] - bmean[k  ] ) )
                                              + cg1<TF>*( ci0<TF>*( b[ijk-kk2] - bmean[k-2] ) + ci1<TF>*( b[ijk-kk1] - bmean[k-1] ) + ci2<TF>*( b[ijk    ] - bmean[k  ] ) + ci3<TF>*( b[ijk+kk1] - bmean[k+1] ) )
                                              + cg2<TF>*( ci0<TF>*( b[ijk-kk1] - bmean[k-1] ) + ci1<TF>*( b[ijk    ] - bmean[k  ] ) + ci2<TF>*( b[ijk+kk1] - bmean[k+1] ) + ci3<TF>*( b[ijk+kk2] - bmean[k+2] ) )
                                              + cg3<TF>*( ti0<TF>*( b[ijk-kk1] - bmean[k-1] ) + ti1<TF>*( b[ijk    ] - bmean[k  ] ) + ti2<TF>*( b[ijk+kk1] - bmean[k+1] ) + ti3<TF>*( b[ijk+kk2] - bmean[k+2] ) ) )

                                            * dzi4[k] )

                                  , 2 ) ) );
            }
    }

    template<typename TF>
    void calc_bw_budget_shear_turb_visc(
            TF* restrict bw_shear, TF* restrict bw_turb, TF* restrict bw_visc,
            TF* restrict bz,
            const TF* restrict w, const TF* restrict p, const TF* restrict b,
            const TF* restrict pmean, const TF* restrict bmean,
            const TF* restrict dzi4, const TF* restrict dzhi4,
            const TF dxi, const TF dyi, const TF dzhi4bot, const TF dzhi4top,
            const TF visc,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int jcells, const int ijcells)
    {
        using namespace Finite_difference::O4;

        const int jj1 = 1*icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;
        const int kk3 = 3*ijcells;
        const int kk4 = 4*ijcells;

        // 0. Create an interpolated field for b on the cell face.
        int k = kstart-1;
        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bz[ijk] = ( bi0<TF>*( b[ijk-kk1] - bmean[k-1] ) + bi1<TF>*( b[ijk    ] - bmean[k  ] ) + bi2<TF>*( b[ijk+kk1] - bmean[k+1] ) + bi3<TF>*( b[ijk+kk2] - bmean[k+2] ) );
            }

        for (int k=kstart; k<kend+1; ++k)
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    bz[ijk] = ( ci0<TF>*( b[ijk-kk2] - bmean[k-2] ) + ci1<TF>*( b[ijk-kk1] - bmean[k-1] ) + ci2<TF>*( b[ijk    ] - bmean[k  ] ) + ci3<TF>*( b[ijk+kk1] - bmean[k+1] ) );
                }

        k = kend+1;
        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bz[ijk] = ( ti0<TF>*( b[ijk-kk3] - bmean[k-3] ) + ti1<TF>*( b[ijk-kk2] - bmean[k-2] ) + ti2<TF>*( b[ijk-kk1] - bmean[k-1] ) + ti3<TF>*( b[ijk    ] - bmean[k  ] ) );
            }


        // 1. CALCULATE THE GRADIENT PRODUCTION TERM
        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    bw_shear[ijk] = - ( ( std::pow( w[ijk], 2 ) * ( cg0<TF>*bmean[k-2] + cg1<TF>*bmean[k-1] + cg2<TF>*bmean[k  ] + cg3<TF>*bmean[k+1] ) ) * dzhi4[k] );
                }

        // 2. CALCULATE THE TURBULENT TRANSPORT TERM
        k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_turb[ijk] = -( ( bg0<TF>*( std::pow( ( bi0<TF>*w[ijk-kk1] + bi1<TF>*w[ijk    ] + bi2<TF>*w[ijk+kk1] + bi3<TF>*w[ijk+kk2] ), 2 ) * ( b[ijk-kk1] - bmean[k-1] ) )
                                  + bg1<TF>*( std::pow( ( ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2] ), 2 ) * ( b[ijk    ] - bmean[k  ] ) )
                                  + bg2<TF>*( std::pow( ( ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3] ), 2 ) * ( b[ijk+kk1] - bmean[k+1] ) )
                                  + bg3<TF>*( std::pow( ( ci0<TF>*w[ijk+kk1] + ci1<TF>*w[ijk+kk2] + ci2<TF>*w[ijk+kk3] + ci3<TF>*w[ijk+kk4] ), 2 ) * ( b[ijk+kk2] - bmean[k+2] ) ) )

                                * dzhi4bot );
            }

        k = kstart+1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_turb[ijk] = -( ( cg0<TF>*( std::pow( ( bi0<TF>*w[ijk-kk2] + bi1<TF>*w[ijk-kk1] + bi2<TF>*w[ijk    ] + bi3<TF>*w[ijk+kk1] ), 2 ) * ( b[ijk-kk2] - bmean[k-2] ) )
                                  + cg1<TF>*( std::pow( ( ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1] ), 2 ) * ( b[ijk-kk1] - bmean[k-1] ) )
                                  + cg2<TF>*( std::pow( ( ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2] ), 2 ) * ( b[ijk    ] - bmean[k  ] ) )
                                  + cg3<TF>*( std::pow( ( ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3] ), 2 ) * ( b[ijk+kk1] - bmean[k+1] ) ) )

                                * dzhi4[k] );
            }

        for (int k=kstart+2; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    bw_turb[ijk] = - ( ( cg0<TF>*( std::pow( ( ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ] ), 2 ) * ( b[ijk-kk2] - bmean[k-2] ) )
                                       + cg1<TF>*( std::pow( ( ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1] ), 2 ) * ( b[ijk-kk1] - bmean[k-1] ) )
                                       + cg2<TF>*( std::pow( ( ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2] ), 2 ) * ( b[ijk    ] - bmean[k  ] ) )
                                       + cg3<TF>*( std::pow( ( ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3] ), 2 ) * ( b[ijk+kk1] - bmean[k+1] ) ) )

                                     * dzhi4[k] );
                }

        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_turb[ijk] = - ( ( cg0<TF>*( std::pow( ( ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ] ), 2 ) * ( b[ijk-kk2] - bmean[k-2] ) )
                                   + cg1<TF>*( std::pow( ( ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1] ), 2 ) * ( b[ijk-kk1] - bmean[k-1] ) )
                                   + cg2<TF>*( std::pow( ( ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2] ), 2 ) * ( b[ijk    ] - bmean[k  ] ) )
                                   + cg3<TF>*( std::pow( ( ti0<TF>*w[ijk-kk1] + ti1<TF>*w[ijk    ] + ti2<TF>*w[ijk+kk1] + ti3<TF>*w[ijk+kk2] ), 2 ) * ( b[ijk+kk1] - bmean[k+1] ) ) )

                                 * dzhi4[k] );
            }

        k = kend;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_turb[ijk] = - ( ( tg0<TF>*( std::pow( ( ci0<TF>*w[ijk-kk4] + ci1<TF>*w[ijk-kk3] + ci2<TF>*w[ijk-kk2] + ci3<TF>*w[ijk-kk1] ), 2 ) * ( b[ijk-kk3] - bmean[k-3] ) )
                                   + tg1<TF>*( std::pow( ( ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ] ), 2 ) * ( b[ijk-kk2] - bmean[k-2] ) )
                                   + tg2<TF>*( std::pow( ( ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1] ), 2 ) * ( b[ijk-kk1] - bmean[k-1] ) )
                                   + tg3<TF>*( std::pow( ( ti0<TF>*w[ijk-kk2] + ti1<TF>*w[ijk-kk1] + ti2<TF>*w[ijk    ] + ti3<TF>*w[ijk+kk1] ), 2 ) * ( b[ijk    ] - bmean[k  ] ) ) )

                                 * dzhi4top );
            }


        // 3. CALCULATE THE VISCOUS TRANSPORT TERM
        k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_visc[ijk] = ( ( visc

                                * ( bg0<TF>*( ( bg0<TF>*( w[ijk-kk1] * bz[ijk-kk1] ) + bg1<TF>*( w[ijk    ] * bz[ijk    ] ) + bg2<TF>*( w[ijk+kk1] * bz[ijk+kk1] ) + bg3<TF>*( w[ijk+kk2] * bz[ijk+kk2] ) ) * dzi4[k-1] )
                                  + bg1<TF>*( ( cg0<TF>*( w[ijk-kk1] * bz[ijk-kk1] ) + cg1<TF>*( w[ijk    ] * bz[ijk    ] ) + cg2<TF>*( w[ijk+kk1] * bz[ijk+kk1] ) + cg3<TF>*( w[ijk+kk2] * bz[ijk+kk2] ) ) * dzi4[k  ] )
                                  + bg2<TF>*( ( cg0<TF>*( w[ijk    ] * bz[ijk    ] ) + cg1<TF>*( w[ijk+kk1] * bz[ijk+kk1] ) + cg2<TF>*( w[ijk+kk2] * bz[ijk+kk2] ) + cg3<TF>*( w[ijk+kk3] * bz[ijk+kk3] ) ) * dzi4[k+1] )
                                  + bg3<TF>*( ( cg0<TF>*( w[ijk+kk1] * bz[ijk+kk1] ) + cg1<TF>*( w[ijk+kk2] * bz[ijk+kk2] ) + cg2<TF>*( w[ijk+kk3] * bz[ijk+kk3] ) + cg3<TF>*( w[ijk+kk4] * bz[ijk+kk4] ) ) * dzi4[k+2] ) ) )

                              * dzhi4bot );
            }

        k = kstart+1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_visc[ijk] = ( ( visc

                                * ( cg0<TF>*( ( bg0<TF>*( w[ijk-kk2] * bz[ijk-kk2] ) + bg1<TF>*( w[ijk-kk1] * bz[ijk-kk1] ) + bg2<TF>*( w[ijk    ] * bz[ijk    ] ) + bg3<TF>*( w[ijk+kk1] * bz[ijk+kk1] ) ) * dzi4[k-2] )
                                  + cg1<TF>*( ( cg0<TF>*( w[ijk-kk2] * bz[ijk-kk2] ) + cg1<TF>*( w[ijk-kk1] * bz[ijk-kk1] ) + cg2<TF>*( w[ijk    ] * bz[ijk    ] ) + cg3<TF>*( w[ijk+kk1] * bz[ijk+kk1] ) ) * dzi4[k-1] )
                                  + cg2<TF>*( ( cg0<TF>*( w[ijk-kk1] * bz[ijk-kk1] ) + cg1<TF>*( w[ijk    ] * bz[ijk    ] ) + cg2<TF>*( w[ijk+kk1] * bz[ijk+kk1] ) + cg3<TF>*( w[ijk+kk2] * bz[ijk+kk2] ) ) * dzi4[k  ] )
                                  + cg3<TF>*( ( cg0<TF>*( w[ijk    ] * bz[ijk    ] ) + cg1<TF>*( w[ijk+kk1] * bz[ijk+kk1] ) + cg2<TF>*( w[ijk+kk2] * bz[ijk+kk2] ) + cg3<TF>*( w[ijk+kk3] * bz[ijk+kk3] ) ) * dzi4[k+1] ) ) )

                              * dzhi4[k] );
            }

        for (int k=kstart+2; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    bw_visc[ijk] = ( ( visc

                                    * ( cg0<TF>*( ( cg0<TF>*( w[ijk-kk3] * bz[ijk-kk3] ) + cg1<TF>*( w[ijk-kk2] * bz[ijk-kk2] ) + cg2<TF>*( w[ijk-kk1] * bz[ijk-kk1] ) + cg3<TF>*( w[ijk    ] * bz[ijk    ] ) ) * dzi4[k-2] )
                                      + cg1<TF>*( ( cg0<TF>*( w[ijk-kk2] * bz[ijk-kk2] ) + cg1<TF>*( w[ijk-kk1] * bz[ijk-kk1] ) + cg2<TF>*( w[ijk    ] * bz[ijk    ] ) + cg3<TF>*( w[ijk+kk1] * bz[ijk+kk1] ) ) * dzi4[k-1] )
                                      + cg2<TF>*( ( cg0<TF>*( w[ijk-kk1] * bz[ijk-kk1] ) + cg1<TF>*( w[ijk    ] * bz[ijk    ] ) + cg2<TF>*( w[ijk+kk1] * bz[ijk+kk1] ) + cg3<TF>*( w[ijk+kk2] * bz[ijk+kk2] ) ) * dzi4[k  ] )
                                      + cg3<TF>*( ( cg0<TF>*( w[ijk    ] * bz[ijk    ] ) + cg1<TF>*( w[ijk+kk1] * bz[ijk+kk1] ) + cg2<TF>*( w[ijk+kk2] * bz[ijk+kk2] ) + cg3<TF>*( w[ijk+kk3] * bz[ijk+kk3] ) ) * dzi4[k+1] ) ) )

                                  * dzhi4[k] );
                }

        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_visc[ijk] = ( ( visc

                                * ( cg0<TF>*( ( cg0<TF>*( w[ijk-kk3] * bz[ijk-kk3] ) + cg1<TF>*( w[ijk-kk2] * bz[ijk-kk2] ) + cg2<TF>*( w[ijk-kk1] * bz[ijk-kk1] ) + cg3<TF>*( w[ijk    ] * bz[ijk    ] ) ) * dzi4[k-2] )
                                  + cg1<TF>*( ( cg0<TF>*( w[ijk-kk2] * bz[ijk-kk2] ) + cg1<TF>*( w[ijk-kk1] * bz[ijk-kk1] ) + cg2<TF>*( w[ijk    ] * bz[ijk    ] ) + cg3<TF>*( w[ijk+kk1] * bz[ijk+kk1] ) ) * dzi4[k-1] )
                                  + cg2<TF>*( ( cg0<TF>*( w[ijk-kk1] * bz[ijk-kk1] ) + cg1<TF>*( w[ijk    ] * bz[ijk    ] ) + cg2<TF>*( w[ijk+kk1] * bz[ijk+kk1] ) + cg3<TF>*( w[ijk+kk2] * bz[ijk+kk2] ) ) * dzi4[k  ] )
                                  + cg3<TF>*( ( tg0<TF>*( w[ijk-kk1] * bz[ijk-kk1] ) + tg1<TF>*( w[ijk    ] * bz[ijk    ] ) + tg2<TF>*( w[ijk+kk1] * bz[ijk+kk1] ) + tg3<TF>*( w[ijk+kk2] * bz[ijk+kk2] ) ) * dzi4[k+1] ) ) )

                              * dzhi4[k] );
            }

        k = kend;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_visc[ijk] = ( ( visc

                                * ( tg0<TF>*( ( cg0<TF>*( w[ijk-kk4] * bz[ijk-kk4] ) + cg1<TF>*( w[ijk-kk3] * bz[ijk-kk3] ) + cg2<TF>*( w[ijk-kk2] * bz[ijk-kk2] ) + cg3<TF>*( w[ijk-kk1] * bz[ijk-kk1] ) ) * dzi4[k-3] )
                                  + tg1<TF>*( ( cg0<TF>*( w[ijk-kk3] * bz[ijk-kk3] ) + cg1<TF>*( w[ijk-kk2] * bz[ijk-kk2] ) + cg2<TF>*( w[ijk-kk1] * bz[ijk-kk1] ) + cg3<TF>*( w[ijk    ] * bz[ijk    ] ) ) * dzi4[k-2] )
                                  + tg2<TF>*( ( cg0<TF>*( w[ijk-kk2] * bz[ijk-kk2] ) + cg1<TF>*( w[ijk-kk1] * bz[ijk-kk1] ) + cg2<TF>*( w[ijk    ] * bz[ijk    ] ) + cg3<TF>*( w[ijk+kk1] * bz[ijk+kk1] ) ) * dzi4[k-1] )
                                  + tg3<TF>*( ( tg0<TF>*( w[ijk-kk2] * bz[ijk-kk2] ) + tg1<TF>*( w[ijk-kk1] * bz[ijk-kk1] ) + tg2<TF>*( w[ijk    ] * bz[ijk    ] ) + tg3<TF>*( w[ijk+kk1] * bz[ijk+kk1] ) ) * dzi4[k  ] ) ) )

                              * dzhi4top );
            }
    }

    template<typename TF>
    void calc_bw_budget_buoy_rdstr_diss_pres(
            TF* restrict bw_buoy, TF* restrict bw_rdstr, TF* restrict bw_diss, TF* restrict bw_pres,
            TF* restrict bz,
            const TF* restrict w, const TF* restrict p, const TF* restrict b,
            const TF* restrict pmean, const TF* restrict bmean,
            const TF* restrict dzi4, const TF* restrict dzhi4,
            const TF dxi, const TF dyi, const TF dzhi4bot, const TF dzhi4top,
            const TF visc,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int jcells, const int ijcells)
    {
        using namespace Finite_difference::O4;

        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj1 = 1*icells;
        const int jj2 = 2*icells;
        const int jj3 = 3*icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;
        const int kk3 = 3*ijcells;
        const int kk4 = 4*ijcells;

        // 4. CALCULATE THE BUOYANCY TERM
        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    bw_buoy[ijk] = std::pow( bz[ijk], 2 );
                }

        // 5. CALCULATE THE REDISTRIBUTION TERM
        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    bw_rdstr[ijk] = ( ( ( ci0<TF>*( p[ijk-kk2] - pmean[k-2] ) + ci1<TF>*( p[ijk-kk1] - pmean[k-1] ) + ci2<TF>*( p[ijk    ] - pmean[k  ] ) + ci3<TF>*( p[ijk+kk1] - pmean[k+1] ) ) * ( cg0<TF>*( b[ijk-kk2] - bmean[k-2] ) + cg1<TF>*( b[ijk-kk1] - bmean[k-1] ) + cg2<TF>*( b[ijk    ] - bmean[k  ] ) + cg3<TF>*( b[ijk+kk1] - bmean[k+1] ) ) ) * dzhi4[k] );
                }

        // 6. CALCULATE THE DISSIPATION TERM
        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_diss[ijk] = - ( ( 2.0 * visc )

                              * ( ( ( ( ( ( cg0<TF>*( ci0<TF>*w[ijk-ii3] + ci1<TF>*w[ijk-ii2] + ci2<TF>*w[ijk-ii1] + ci3<TF>*w[ijk    ] )
                                          + cg1<TF>*( ci0<TF>*w[ijk-ii2] + ci1<TF>*w[ijk-ii1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1] )
                                          + cg2<TF>*( ci0<TF>*w[ijk-ii1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+ii1] + ci3<TF>*w[ijk+ii2] )
                                          + cg3<TF>*( ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+ii1] + ci2<TF>*w[ijk+ii2] + ci3<TF>*w[ijk+ii3] ) )

                                        * dxi )

                                      * ( cg0<TF>*( ci0<TF>*bz[ijk-ii3] + ci1<TF>*bz[ijk-ii2] + ci2<TF>*bz[ijk-ii1] + ci3<TF>*bz[ijk    ] )
                                        + cg1<TF>*( ci0<TF>*bz[ijk-ii2] + ci1<TF>*bz[ijk-ii1] + ci2<TF>*bz[ijk    ] + ci3<TF>*bz[ijk+ii1] )
                                        + cg2<TF>*( ci0<TF>*bz[ijk-ii1] + ci1<TF>*bz[ijk    ] + ci2<TF>*bz[ijk+ii1] + ci3<TF>*bz[ijk+ii2] )
                                        + cg3<TF>*( ci0<TF>*bz[ijk    ] + ci1<TF>*bz[ijk+ii1] + ci2<TF>*bz[ijk+ii2] + ci3<TF>*bz[ijk+ii3] ) ) )

                                    * dxi )

                                  + ( ( ( ( cg0<TF>*( ci0<TF>*w[ijk-jj3] + ci1<TF>*w[ijk-jj2] + ci2<TF>*w[ijk-jj1] + ci3<TF>*w[ijk    ] )
                                          + cg1<TF>*( ci0<TF>*w[ijk-jj2] + ci1<TF>*w[ijk-jj1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1] )
                                          + cg2<TF>*( ci0<TF>*w[ijk-jj1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+jj1] + ci3<TF>*w[ijk+jj2] )
                                          + cg3<TF>*( ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+jj1] + ci2<TF>*w[ijk+jj2] + ci3<TF>*w[ijk+jj3] ) )

                                        * dyi )

                                      * ( cg0<TF>*( ci0<TF>*bz[ijk-jj3] + ci1<TF>*bz[ijk-jj2] + ci2<TF>*bz[ijk-jj1] + ci3<TF>*bz[ijk    ] )
                                        + cg1<TF>*( ci0<TF>*bz[ijk-jj2] + ci1<TF>*bz[ijk-jj1] + ci2<TF>*bz[ijk    ] + ci3<TF>*bz[ijk+jj1] )
                                        + cg2<TF>*( ci0<TF>*bz[ijk-jj1] + ci1<TF>*bz[ijk    ] + ci2<TF>*bz[ijk+jj1] + ci3<TF>*bz[ijk+jj2] )
                                        + cg3<TF>*( ci0<TF>*bz[ijk    ] + ci1<TF>*bz[ijk+jj1] + ci2<TF>*bz[ijk+jj2] + ci3<TF>*bz[ijk+jj3] ) ) )

                                    * dyi ) )

                                + ( ( ( ( bg0<TF>*( bi0<TF>*w[ijk-kk1] + bi1<TF>*w[ijk    ] + bi2<TF>*w[ijk+kk1] + bi3<TF>*w[ijk+kk2] )
                                        + bg1<TF>*( ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2] )
                                        + bg2<TF>*( ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3] )
                                        + bg3<TF>*( ci0<TF>*w[ijk+kk1] + ci1<TF>*w[ijk+kk2] + ci2<TF>*w[ijk+kk3] + ci3<TF>*w[ijk+kk4] ) )

                                      * dzhi4bot )

                                    * ( cg0<TF>*( b[ijk-kk2] - bmean[k-2] ) + cg1<TF>*( b[ijk-kk1] - bmean[k-1] ) + cg2<TF>*( b[ijk    ] - bmean[k  ] ) + cg3<TF>*( b[ijk+kk1] - bmean[k+1] ) ) )

                                  * dzhi4bot ) ) );
            }

        k = kstart+1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_diss[ijk] = - ( ( 2.0 * visc )

                              * ( ( ( ( ( ( cg0<TF>*( ci0<TF>*w[ijk-ii3] + ci1<TF>*w[ijk-ii2] + ci2<TF>*w[ijk-ii1] + ci3<TF>*w[ijk    ] )
                                          + cg1<TF>*( ci0<TF>*w[ijk-ii2] + ci1<TF>*w[ijk-ii1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1] )
                                          + cg2<TF>*( ci0<TF>*w[ijk-ii1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+ii1] + ci3<TF>*w[ijk+ii2] )
                                          + cg3<TF>*( ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+ii1] + ci2<TF>*w[ijk+ii2] + ci3<TF>*w[ijk+ii3] ) )

                                        * dxi )

                                      * ( cg0<TF>*( ci0<TF>*bz[ijk-ii3] + ci1<TF>*bz[ijk-ii2] + ci2<TF>*bz[ijk-ii1] + ci3<TF>*bz[ijk    ] )
                                        + cg1<TF>*( ci0<TF>*bz[ijk-ii2] + ci1<TF>*bz[ijk-ii1] + ci2<TF>*bz[ijk    ] + ci3<TF>*bz[ijk+ii1] )
                                        + cg2<TF>*( ci0<TF>*bz[ijk-ii1] + ci1<TF>*bz[ijk    ] + ci2<TF>*bz[ijk+ii1] + ci3<TF>*bz[ijk+ii2] )
                                        + cg3<TF>*( ci0<TF>*bz[ijk    ] + ci1<TF>*bz[ijk+ii1] + ci2<TF>*bz[ijk+ii2] + ci3<TF>*bz[ijk+ii3] ) ) )

                                    * dxi )

                                  + ( ( ( ( cg0<TF>*( ci0<TF>*w[ijk-jj3] + ci1<TF>*w[ijk-jj2] + ci2<TF>*w[ijk-jj1] + ci3<TF>*w[ijk    ] )
                                          + cg1<TF>*( ci0<TF>*w[ijk-jj2] + ci1<TF>*w[ijk-jj1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1] )
                                          + cg2<TF>*( ci0<TF>*w[ijk-jj1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+jj1] + ci3<TF>*w[ijk+jj2] )
                                          + cg3<TF>*( ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+jj1] + ci2<TF>*w[ijk+jj2] + ci3<TF>*w[ijk+jj3] ) )

                                        * dyi )

                                      * ( cg0<TF>*( ci0<TF>*bz[ijk-jj3] + ci1<TF>*bz[ijk-jj2] + ci2<TF>*bz[ijk-jj1] + ci3<TF>*bz[ijk    ] )
                                        + cg1<TF>*( ci0<TF>*bz[ijk-jj2] + ci1<TF>*bz[ijk-jj1] + ci2<TF>*bz[ijk    ] + ci3<TF>*bz[ijk+jj1] )
                                        + cg2<TF>*( ci0<TF>*bz[ijk-jj1] + ci1<TF>*bz[ijk    ] + ci2<TF>*bz[ijk+jj1] + ci3<TF>*bz[ijk+jj2] )
                                        + cg3<TF>*( ci0<TF>*bz[ijk    ] + ci1<TF>*bz[ijk+jj1] + ci2<TF>*bz[ijk+jj2] + ci3<TF>*bz[ijk+jj3] ) ) )

                                    * dyi ) )

                                + ( ( ( ( cg0<TF>*( bi0<TF>*w[ijk-kk2] + bi1<TF>*w[ijk-kk1] + bi2<TF>*w[ijk    ] + bi3<TF>*w[ijk+kk1] )
                                        + cg1<TF>*( ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1] )
                                        + cg2<TF>*( ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2] )
                                        + cg3<TF>*( ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3] ) )

                                      * dzhi4[k] )

                                    * ( cg0<TF>*( b[ijk-kk2] - bmean[k-2] ) + cg1<TF>*( b[ijk-kk1] - bmean[k-1] ) + cg2<TF>*( b[ijk    ] - bmean[k  ] ) + cg3<TF>*( b[ijk+kk1] - bmean[k+1] ) ) )

                                  * dzhi4[k] ) ) );
            }

        for (int k=kstart+2; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    bw_diss[ijk] = - ( ( 2.0 * visc )

                                  * ( ( ( ( ( ( cg0<TF>*( ci0<TF>*w[ijk-ii3] + ci1<TF>*w[ijk-ii2] + ci2<TF>*w[ijk-ii1] + ci3<TF>*w[ijk    ] )
                                              + cg1<TF>*( ci0<TF>*w[ijk-ii2] + ci1<TF>*w[ijk-ii1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1] )
                                              + cg2<TF>*( ci0<TF>*w[ijk-ii1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+ii1] + ci3<TF>*w[ijk+ii2] )
                                              + cg3<TF>*( ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+ii1] + ci2<TF>*w[ijk+ii2] + ci3<TF>*w[ijk+ii3] ) )

                                            * dxi )

                                          * ( cg0<TF>*( ci0<TF>*bz[ijk-ii3] + ci1<TF>*bz[ijk-ii2] + ci2<TF>*bz[ijk-ii1] + ci3<TF>*bz[ijk    ] )
                                            + cg1<TF>*( ci0<TF>*bz[ijk-ii2] + ci1<TF>*bz[ijk-ii1] + ci2<TF>*bz[ijk    ] + ci3<TF>*bz[ijk+ii1] )
                                            + cg2<TF>*( ci0<TF>*bz[ijk-ii1] + ci1<TF>*bz[ijk    ] + ci2<TF>*bz[ijk+ii1] + ci3<TF>*bz[ijk+ii2] )
                                            + cg3<TF>*( ci0<TF>*bz[ijk    ] + ci1<TF>*bz[ijk+ii1] + ci2<TF>*bz[ijk+ii2] + ci3<TF>*bz[ijk+ii3] ) ) )

                                        * dxi )

                                      + ( ( ( ( cg0<TF>*( ci0<TF>*w[ijk-jj3] + ci1<TF>*w[ijk-jj2] + ci2<TF>*w[ijk-jj1] + ci3<TF>*w[ijk    ] )
                                              + cg1<TF>*( ci0<TF>*w[ijk-jj2] + ci1<TF>*w[ijk-jj1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1] )
                                              + cg2<TF>*( ci0<TF>*w[ijk-jj1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+jj1] + ci3<TF>*w[ijk+jj2] )
                                              + cg3<TF>*( ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+jj1] + ci2<TF>*w[ijk+jj2] + ci3<TF>*w[ijk+jj3] ) )

                                            * dyi )

                                          * ( cg0<TF>*( ci0<TF>*bz[ijk-jj3] + ci1<TF>*bz[ijk-jj2] + ci2<TF>*bz[ijk-jj1] + ci3<TF>*bz[ijk    ] )
                                            + cg1<TF>*( ci0<TF>*bz[ijk-jj2] + ci1<TF>*bz[ijk-jj1] + ci2<TF>*bz[ijk    ] + ci3<TF>*bz[ijk+jj1] )
                                            + cg2<TF>*( ci0<TF>*bz[ijk-jj1] + ci1<TF>*bz[ijk    ] + ci2<TF>*bz[ijk+jj1] + ci3<TF>*bz[ijk+jj2] )
                                            + cg3<TF>*( ci0<TF>*bz[ijk    ] + ci1<TF>*bz[ijk+jj1] + ci2<TF>*bz[ijk+jj2] + ci3<TF>*bz[ijk+jj3] ) ) )

                                        * dyi ) )

                                    + ( ( ( ( cg0<TF>*( ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ] )
                                            + cg1<TF>*( ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1] )
                                            + cg2<TF>*( ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2] )
                                            + cg3<TF>*( ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3] ) )

                                          * dzhi4[k] )

                                        * ( cg0<TF>*( b[ijk-kk2] - bmean[k-2] ) + cg1<TF>*( b[ijk-kk1] - bmean[k-1] ) + cg2<TF>*( b[ijk    ] - bmean[k  ] ) + cg3<TF>*( b[ijk+kk1] - bmean[k+1] ) ) )

                                      * dzhi4[k] ) ) );
                }

        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_diss[ijk] = - ( ( 2.0 * visc )

                              * ( ( ( ( ( ( cg0<TF>*( ci0<TF>*w[ijk-ii3] + ci1<TF>*w[ijk-ii2] + ci2<TF>*w[ijk-ii1] + ci3<TF>*w[ijk    ] )
                                          + cg1<TF>*( ci0<TF>*w[ijk-ii2] + ci1<TF>*w[ijk-ii1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1] )
                                          + cg2<TF>*( ci0<TF>*w[ijk-ii1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+ii1] + ci3<TF>*w[ijk+ii2] )
                                          + cg3<TF>*( ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+ii1] + ci2<TF>*w[ijk+ii2] + ci3<TF>*w[ijk+ii3] ) )

                                        * dxi )

                                      * ( cg0<TF>*( ci0<TF>*bz[ijk-ii3] + ci1<TF>*bz[ijk-ii2] + ci2<TF>*bz[ijk-ii1] + ci3<TF>*bz[ijk    ] )
                                        + cg1<TF>*( ci0<TF>*bz[ijk-ii2] + ci1<TF>*bz[ijk-ii1] + ci2<TF>*bz[ijk    ] + ci3<TF>*bz[ijk+ii1] )
                                        + cg2<TF>*( ci0<TF>*bz[ijk-ii1] + ci1<TF>*bz[ijk    ] + ci2<TF>*bz[ijk+ii1] + ci3<TF>*bz[ijk+ii2] )
                                        + cg3<TF>*( ci0<TF>*bz[ijk    ] + ci1<TF>*bz[ijk+ii1] + ci2<TF>*bz[ijk+ii2] + ci3<TF>*bz[ijk+ii3] ) ) )

                                    * dxi )

                                  + ( ( ( ( cg0<TF>*( ci0<TF>*w[ijk-jj3] + ci1<TF>*w[ijk-jj2] + ci2<TF>*w[ijk-jj1] + ci3<TF>*w[ijk    ] )
                                          + cg1<TF>*( ci0<TF>*w[ijk-jj2] + ci1<TF>*w[ijk-jj1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1] )
                                          + cg2<TF>*( ci0<TF>*w[ijk-jj1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+jj1] + ci3<TF>*w[ijk+jj2] )
                                          + cg3<TF>*( ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+jj1] + ci2<TF>*w[ijk+jj2] + ci3<TF>*w[ijk+jj3] ) )

                                        * dyi )

                                      * ( cg0<TF>*( ci0<TF>*bz[ijk-jj3] + ci1<TF>*bz[ijk-jj2] + ci2<TF>*bz[ijk-jj1] + ci3<TF>*bz[ijk    ] )
                                        + cg1<TF>*( ci0<TF>*bz[ijk-jj2] + ci1<TF>*bz[ijk-jj1] + ci2<TF>*bz[ijk    ] + ci3<TF>*bz[ijk+jj1] )
                                        + cg2<TF>*( ci0<TF>*bz[ijk-jj1] + ci1<TF>*bz[ijk    ] + ci2<TF>*bz[ijk+jj1] + ci3<TF>*bz[ijk+jj2] )
                                        + cg3<TF>*( ci0<TF>*bz[ijk    ] + ci1<TF>*bz[ijk+jj1] + ci2<TF>*bz[ijk+jj2] + ci3<TF>*bz[ijk+jj3] ) ) )

                                    * dyi ) )

                                + ( ( ( ( cg0<TF>*( ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ] )
                                        + cg1<TF>*( ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1] )
                                        + cg2<TF>*( ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2] )
                                        + cg3<TF>*( ti0<TF>*w[ijk-kk1] + ti1<TF>*w[ijk    ] + ti2<TF>*w[ijk+kk1] + ti3<TF>*w[ijk+kk2] ) )

                                      * dzhi4[k] )

                                    * ( cg0<TF>*( b[ijk-kk2] - bmean[k-2] ) + cg1<TF>*( b[ijk-kk1] - bmean[k-1] ) + cg2<TF>*( b[ijk    ] - bmean[k  ] ) + cg3<TF>*( b[ijk+kk1] - bmean[k+1] ) ) )

                                  * dzhi4[k] ) ) );


            }

        k = kend;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                bw_diss[ijk] = - ( ( 2.0 * visc )

                              * ( ( ( ( ( ( cg0<TF>*( ci0<TF>*w[ijk-ii3] + ci1<TF>*w[ijk-ii2] + ci2<TF>*w[ijk-ii1] + ci3<TF>*w[ijk    ] )
                                          + cg1<TF>*( ci0<TF>*w[ijk-ii2] + ci1<TF>*w[ijk-ii1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1] )
                                          + cg2<TF>*( ci0<TF>*w[ijk-ii1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+ii1] + ci3<TF>*w[ijk+ii2] )
                                          + cg3<TF>*( ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+ii1] + ci2<TF>*w[ijk+ii2] + ci3<TF>*w[ijk+ii3] ) )

                                        * dxi )

                                      * ( cg0<TF>*( ci0<TF>*bz[ijk-ii3] + ci1<TF>*bz[ijk-ii2] + ci2<TF>*bz[ijk-ii1] + ci3<TF>*bz[ijk    ] )
                                        + cg1<TF>*( ci0<TF>*bz[ijk-ii2] + ci1<TF>*bz[ijk-ii1] + ci2<TF>*bz[ijk    ] + ci3<TF>*bz[ijk+ii1] )
                                        + cg2<TF>*( ci0<TF>*bz[ijk-ii1] + ci1<TF>*bz[ijk    ] + ci2<TF>*bz[ijk+ii1] + ci3<TF>*bz[ijk+ii2] )
                                        + cg3<TF>*( ci0<TF>*bz[ijk    ] + ci1<TF>*bz[ijk+ii1] + ci2<TF>*bz[ijk+ii2] + ci3<TF>*bz[ijk+ii3] ) ) )

                                    * dxi )

                                  + ( ( ( ( cg0<TF>*( ci0<TF>*w[ijk-jj3] + ci1<TF>*w[ijk-jj2] + ci2<TF>*w[ijk-jj1] + ci3<TF>*w[ijk    ] )
                                          + cg1<TF>*( ci0<TF>*w[ijk-jj2] + ci1<TF>*w[ijk-jj1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1] )
                                          + cg2<TF>*( ci0<TF>*w[ijk-jj1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+jj1] + ci3<TF>*w[ijk+jj2] )
                                          + cg3<TF>*( ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+jj1] + ci2<TF>*w[ijk+jj2] + ci3<TF>*w[ijk+jj3] ) )

                                        * dyi )

                                      * ( cg0<TF>*( ci0<TF>*bz[ijk-jj3] + ci1<TF>*bz[ijk-jj2] + ci2<TF>*bz[ijk-jj1] + ci3<TF>*bz[ijk    ] )
                                        + cg1<TF>*( ci0<TF>*bz[ijk-jj2] + ci1<TF>*bz[ijk-jj1] + ci2<TF>*bz[ijk    ] + ci3<TF>*bz[ijk+jj1] )
                                        + cg2<TF>*( ci0<TF>*bz[ijk-jj1] + ci1<TF>*bz[ijk    ] + ci2<TF>*bz[ijk+jj1] + ci3<TF>*bz[ijk+jj2] )
                                        + cg3<TF>*( ci0<TF>*bz[ijk    ] + ci1<TF>*bz[ijk+jj1] + ci2<TF>*bz[ijk+jj2] + ci3<TF>*bz[ijk+jj3] ) ) )

                                    * dyi ) )

                                + ( ( ( ( tg0<TF>*( ci0<TF>*w[ijk-kk4] + ci1<TF>*w[ijk-kk3] + ci2<TF>*w[ijk-kk2] + ci3<TF>*w[ijk-kk1] )
                                        + tg1<TF>*( ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ] )
                                        + tg2<TF>*( ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1] )
                                        + tg3<TF>*( ti0<TF>*w[ijk-kk2] + ti1<TF>*w[ijk-kk1] + ti2<TF>*w[ijk    ] + ti3<TF>*w[ijk+kk1] ) )

                                      * dzhi4top )

                                    * ( cg0<TF>*( b[ijk-kk2] - bmean[k-2] ) + cg1<TF>*( b[ijk-kk1] - bmean[k-1] ) + cg2<TF>*( b[ijk    ] - bmean[k  ] ) + cg3<TF>*( b[ijk+kk1] - bmean[k+1] ) ) )

                                  * dzhi4top ) ) );
            }

        // 7. CALCULATE THE PRESSURE TRANSPORT TERM
        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    bw_pres[ijk] = - ( ( cg0<TF>*( ( p[ijk-kk2] - pmean[k-2] ) * ( b[ijk-kk2] - bmean[k-2] ) ) + cg1<TF>*( ( p[ijk-kk1] - pmean[k-1] ) * ( b[ijk-kk1] - bmean[k-1] ) ) + cg2<TF>*( ( p[ijk    ] - pmean[k  ] ) * ( b[ijk    ] - bmean[k  ] ) ) + cg3<TF>*( ( p[ijk+kk1] - pmean[k+1] ) * ( b[ijk+kk1] - bmean[k+1] ) ) ) * dzhi4[k] );
                }
    }

    template<typename TF>
    void calc_sorted_prof(
            TF* restrict data, TF* restrict bin, TF* restrict prof,
            const TF* restrict z, const TF* restrict dz,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells,
            const int itot, const int jtot, const int nmax,
            Grid_order grid_order, Master& master)
    {
        const int jj = icells;
        const int kk = ijcells;

        TF minval =  Constants::dhuge;
        TF maxval = -Constants::dhuge;

        // First, get min and max.
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    if (data[ijk] < minval)
                        minval = data[ijk];
                    if (data[ijk] > maxval)
                        maxval = data[ijk];
                }

        master.min(&minval, 1);
        master.max(&maxval, 1);

        // make sure that the max ends up in the last bin (introduce 1E-9 error)
        maxval *= (1.+Constants::dsmall);

        const TF range = maxval-minval;

        // In case the field is entirely uniform, dbin becomes zero. In that case we set the profile to the minval.
        if (range < 1.e-16)
        {
            for (int k=kstart; k<kend; ++k)
                prof[k] = minval;
        }
        else
        {
            // create bins, equal to the number of grid cells per proc
            // make sure that bins is not larger than the memory of one 3d field
            const int bins = nmax;

            // calculate bin width, subtract one to make the minimum and maximum
            // are in the middle of the bin range and add half a bin size on both sides
            // |----x----|----x----|----x----|
            const TF dbin = range / (TF)(bins-1);

            minval -= 0.5*dbin;
            maxval += 0.5*dbin;

            // set the bin array to zero
            for (int n=0; n<bins; ++n)
                bin[n] = 0;

            // calculate the division factor of one equivalent height unit
            // (the total volume saved is itot*jtot*zsize)
            const TF nslice = (TF)(itot*jtot);

            // check in which bin each value falls and increment the bin count
            for (int k=kstart; k<kend; ++k)
            {
                const TF dzslice = dz[k] / nslice;
                for (int j=jstart; j<jend; ++j)
                    // do not add a ivdep pragma here, because multiple instances could write the same bin[index]
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        const int index = (int)((data[ijk] - minval) / dbin);
                        bin[index] += dzslice;
                    }
            }

            // get the bin count
            master.sum(bin, bins);

            // set the starting values of the loop
            int index = 0;
            TF zbin = 0.5*bin[index];
            TF profval = minval + 0.5*dbin;

            for (int k=kstart; k<kend; ++k)
            {
                // Integrate the profile up to the bin count.
                // Escape the while loop when the integrated profile
                // exceeds the next grid point.
                while (zbin < z[k])
                {
                    zbin += 0.5*(bin[index]+bin[index+1]);
                    profval += dbin;
                    ++index;
                }

                // In case the first bin is larger than the grid spacing, which can happen
                // in the inital phase of an MPI run, make sure that no out-of-bounds reads
                // happen.
                if (index == 0)
                    prof[k] = profval;
                else
                {
                    const TF dzfrac = (zbin-z[k]) / (0.5*(bin[index-1]+bin[index]));
                    prof[k] = profval - dzfrac*dbin;
                }
            }
        }

        // now calculate the ghost cells
        // \TODO this might not be accurate enough, extrapolate properly
        TF profbot = minval;
        TF proftop = maxval;

        if (grid_order == Grid_order::Second)
        {
            prof[kstart-1] = 2.*profbot - prof[kstart];
            prof[kend]     = 2.*proftop - prof[kend-1];
        }
        else if (grid_order == Grid_order::Fourth)
        {
            prof[kstart-1] = (8./3.)*profbot - 2.*prof[kstart] + (1./3.)*prof[kstart+1];
            prof[kstart-2] = 8.*profbot      - 9.*prof[kstart] + 2.*prof[kstart+1];
            prof[kend]     = (8./3.)*proftop - 2.*prof[kend-1] + (1./3.)*prof[kend-2];
            prof[kend+1]   = 8.*proftop      - 9.*prof[kend-1] + 2.*prof[kend-2];
        }
    }
}

template<typename TF>
Budget_4<TF>::Budget_4(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin,
        Thermo<TF>& thermoin, Diff<TF>& diffin, Advec<TF>& advecin, Force<TF>& forcein, Input& inputin) :
    Budget<TF>(masterin, gridin, fieldsin, thermoin, diffin, advecin, forcein, inputin),
    field3d_operators(masterin, gridin, fieldsin)
{}

template<typename TF>
Budget_4<TF>::~Budget_4()
{}

template<typename TF>
void Budget_4<TF>::init()
{
    auto& gd = grid.get_grid_data();

    umodel.resize(gd.kcells);
    vmodel.resize(gd.kcells);
    wmodel.resize(gd.kcells);
}

template<typename TF>
void Budget_4<TF>::create(Stats<TF>& stats)
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

    stats.add_prof("u2_turb" , "Turbulent transport term in U2 budget" , "m2 s-3", "z" , group_name);
    stats.add_prof("v2_turb" , "Turbulent transport term in V2 budget" , "m2 s-3", "z" , group_name);
    stats.add_prof("w2_turb" , "Turbulent transport term in W2 budget" , "m2 s-3", "zh", group_name);
    stats.add_prof("tke_turb", "Turbulent transport term in TKE budget", "m2 s-3", "z" , group_name);
    stats.add_prof("uw_turb" , "Turbulent transport term in UW budget" , "m2 s-3", "zh", group_name);

    stats.add_prof("u2_visc" , "Viscous transport term in U2 budget" , "m2 s-3", "z" , group_name);
    stats.add_prof("v2_visc" , "Viscous transport term in V2 budget" , "m2 s-3", "z" , group_name);
    stats.add_prof("w2_visc" , "Viscous transport term in W2 budget" , "m2 s-3", "zh", group_name);
    stats.add_prof("tke_visc", "Viscous transport term in TKE budget", "m2 s-3", "z" , group_name);
    stats.add_prof("uw_visc" , "Viscous transport term in UW budget" , "m2 s-3", "zh", group_name);

    stats.add_prof("u2_diss" , "Dissipation term in U2 budget" , "m2 s-3", "z" , group_name);
    stats.add_prof("v2_diss" , "Dissipation term in V2 budget" , "m2 s-3", "z" , group_name);
    stats.add_prof("w2_diss" , "Dissipation term in W2 budget" , "m2 s-3", "zh", group_name);
    stats.add_prof("tke_diss", "Dissipation term in TKE budget", "m2 s-3", "z" , group_name);
    stats.add_prof("uw_diss" , "Dissipation term in UW budget" , "m2 s-3", "zh", group_name);

    stats.add_prof("w2_pres" , "Pressure transport term in W2 budget" , "m2 s-3", "zh", group_name);
    stats.add_prof("tke_pres", "Pressure transport term in TKE budget", "m2 s-3", "z" , group_name);
    stats.add_prof("uw_pres" , "Pressure transport term in UW budget" , "m2 s-3", "zh", group_name);

    stats.add_prof("u2_rdstr", "Pressure redistribution term in U2 budget", "m2 s-3", "z" , group_name);
    stats.add_prof("v2_rdstr", "Pressure redistribution term in V2 budget", "m2 s-3", "z" , group_name);
    stats.add_prof("w2_rdstr", "Pressure redistribution term in W2 budget", "m2 s-3", "zh", group_name);
    stats.add_prof("uw_rdstr", "Pressure redistribution term in UW budget", "m2 s-3", "zh", group_name);

    if (thermo.get_switch() != "0")
    {
        stats.add_prof("w2_buoy" , "Buoyancy production/destruction term in W2 budget" , "m2 s-3", "zh", group_name);
        stats.add_prof("tke_buoy", "Buoyancy production/destruction term in TKE budget", "m2 s-3", "z" , group_name);
        stats.add_prof("uw_buoy" , "Buoyancy production/destruction term in UW budget" , "m2 s-3", "zh", group_name);

        stats.add_prof("b2_shear", "Shear production term in B2 budget"   , "m2 s-5", "z", group_name);
        stats.add_prof("b2_turb" , "Turbulent transport term in B2 budget", "m2 s-5", "z", group_name);
        stats.add_prof("b2_visc" , "Viscous transport term in B2 budget"  , "m2 s-5", "z", group_name);
        stats.add_prof("b2_diss" , "Dissipation term in B2 budget"        , "m2 s-5", "z", group_name);

        stats.add_prof("bw_shear", "Shear production term in BW budget"   , "m2 s-4", "zh", group_name);
        stats.add_prof("bw_turb" , "Turbulent transport term in BW budget", "m2 s-4", "zh", group_name);
        stats.add_prof("bw_visc" , "Viscous transport term in BW budget"  , "m2 s-4", "zh", group_name);
        stats.add_prof("bw_rdstr", "Redistribution term in BW budget"     , "m2 s-4", "zh", group_name);
        stats.add_prof("bw_buoy" , "Buoyancy term in BW budget"           , "m2 s-4", "zh", group_name);
        stats.add_prof("bw_diss" , "Dissipation term in BW budget"        , "m2 s-4", "zh", group_name);
        stats.add_prof("bw_pres" , "Pressure transport term in BW budget" , "m2 s-4", "zh", group_name);

        stats.add_prof("b_sort", "Sorted buoyancy", "m s-2", "z", group_name);
    }

    /*
    if (thermo.get_switch() != "0")
    {
        // Add the profiles for the potential energy budget to the statistics.
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
    */
}

template<typename TF>
void Budget_4<TF>::exec_stats(Stats<TF>& stats)
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

        // Calculate the TKE budget.
        auto ke  = fields.get_tmp();
        auto tke = fields.get_tmp();

        const TF no_offset = 0.;
        const TF no_threshold = 0.;

        calc_ke(ke->fld.data(), tke->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                umodel.data(), vmodel.data(), wmodel.data(),
                gd.utrans, gd.vtrans,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_mask_stats(m, "ke" , *ke , no_offset, no_threshold);
        stats.calc_mask_stats(m, "tke", *tke, no_offset, no_threshold);

        // Subtract mean
        auto w_prime = fields.get_tmp();
        calc_prime(
                w_prime->fld.data(), fields.mp.at("w")->fld.data(), wmodel.data(),
                gd.icells, gd.jcells, gd.kcells,
                gd.ijcells);

        auto wx = std::move(ke );
        auto wy = std::move(tke);

        // Interpolate w to the locations of u and v.
        const int wloc [3] = {0,0,1};
        const int wxloc[3] = {1,0,1};
        const int wyloc[3] = {0,1,1};

        grid.interpolate_4th(wx->fld.data(), w_prime->fld.data(), wloc, wxloc);
        grid.interpolate_4th(wy->fld.data(), w_prime->fld.data(), wloc, wyloc);

        auto u2_shear = fields.get_tmp();
        auto v2_shear = fields.get_tmp();
        auto tke_shear = fields.get_tmp();
        auto uw_shear = fields.get_tmp();

        calc_tke_budget_shear(
                u2_shear->fld.data(), v2_shear->fld.data(), tke_shear->fld.data(), uw_shear->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), w_prime->fld.data(),
                wx->fld.data(), wy->fld.data(),
                umodel.data(), vmodel.data(),
                gd.dzi4.data(), gd.dzhi4.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_mask_stats(m, "u2_shear" , *u2_shear , no_offset, no_threshold);
        stats.calc_mask_stats(m, "v2_shear" , *v2_shear , no_offset, no_threshold);
        stats.calc_mask_stats(m, "tke_shear", *tke_shear, no_offset, no_threshold);
        stats.calc_mask_stats(m, "uw_shear" , *uw_shear , no_offset, no_threshold);

        auto u2_turb = std::move(u2_shear);
        auto v2_turb = std::move(v2_shear);
        auto w2_turb = fields.get_tmp();
        auto tke_turb = std::move(tke_shear);
        auto uw_turb = std::move(uw_shear);

        calc_tke_budget_turb(
                u2_turb->fld.data(), v2_turb->fld.data(), w2_turb->fld.data(), tke_turb->fld.data(), uw_turb->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), w_prime->fld.data(),
                wx->fld.data(), wy->fld.data(),
                umodel.data(), vmodel.data(),
                gd.dzi4.data(), gd.dzhi4.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_mask_stats(m, "u2_turb" , *u2_turb , no_offset, no_threshold);
        stats.calc_mask_stats(m, "v2_turb" , *v2_turb , no_offset, no_threshold);
        stats.calc_mask_stats(m, "w2_turb" , *w2_turb , no_offset, no_threshold);
        stats.calc_mask_stats(m, "tke_turb", *tke_turb, no_offset, no_threshold);
        stats.calc_mask_stats(m, "uw_turb" , *uw_turb , no_offset, no_threshold);

        auto w2_pres  = std::move(w2_turb);
        auto tke_pres = std::move(tke_turb);
        auto uw_pres  = std::move(uw_turb);

        calc_tke_budget_pres(
                w2_pres->fld.data(), tke_pres->fld.data(), uw_pres->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(),
                w_prime->fld.data(), fields.sd.at("p")->fld.data(),
                umodel.data(), vmodel.data(),
                gd.dzi4.data(), gd.dzhi4.data(),
                gd.dxi, gd.dyi,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_mask_stats(m, "w2_pres" , *w2_pres , no_offset, no_threshold);
        stats.calc_mask_stats(m, "tke_pres", *tke_pres, no_offset, no_threshold);
        stats.calc_mask_stats(m, "uw_pres" , *uw_pres , no_offset, no_threshold);

        auto u2_visc  = std::move(u2_turb);
        auto v2_visc  = std::move(v2_turb);
        auto w2_visc  = std::move(w2_pres);
        auto tke_visc = std::move(tke_pres);
        auto uw_visc  = std::move(uw_pres);

        auto wz = std::move(wx);
        auto uz = std::move(wy);

        calc_tke_budget_visc(
                u2_visc->fld.data(), v2_visc->fld.data(), w2_visc->fld.data(), tke_visc->fld.data(), uw_visc->fld.data(),
                wz->fld.data(), uz->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), w_prime->fld.data(),
                umodel.data(), vmodel.data(),
                gd.dzi4.data(), gd.dzhi4.data(),
                gd.dxi, gd.dyi, gd.dzhi4bot, gd.dzhi4top,
                fields.visc,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_mask_stats(m, "u2_visc" , *u2_visc , no_offset, no_threshold);
        stats.calc_mask_stats(m, "v2_visc" , *v2_visc , no_offset, no_threshold);
        stats.calc_mask_stats(m, "w2_visc" , *w2_visc , no_offset, no_threshold);
        stats.calc_mask_stats(m, "tke_visc", *tke_visc, no_offset, no_threshold);
        stats.calc_mask_stats(m, "uw_visc" , *uw_visc , no_offset, no_threshold);

        auto u2_diss  = std::move(u2_visc);
        auto v2_diss  = std::move(v2_visc);
        auto w2_diss  = std::move(w2_visc);
        auto tke_diss = std::move(tke_visc);
        auto uw_diss  = std::move(uw_visc);

        calc_tke_budget_diss(
                u2_diss->fld.data(), v2_diss->fld.data(), w2_diss->fld.data(), tke_diss->fld.data(), uw_diss->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), w_prime->fld.data(),
                umodel.data(), vmodel.data(),
                gd.dzi4.data(), gd.dzhi4.data(),
                gd.dxi, gd.dyi, gd.dzhi4bot, gd.dzhi4top,
                fields.visc,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_mask_stats(m, "u2_diss" , *u2_diss , no_offset, no_threshold);
        stats.calc_mask_stats(m, "v2_diss" , *v2_diss , no_offset, no_threshold);
        stats.calc_mask_stats(m, "w2_diss" , *w2_diss , no_offset, no_threshold);
        stats.calc_mask_stats(m, "tke_diss", *tke_diss, no_offset, no_threshold);
        stats.calc_mask_stats(m, "uw_diss" , *uw_diss , no_offset, no_threshold);

        auto u2_rdstr = std::move(u2_diss);
        auto v2_rdstr = std::move(v2_diss);
        auto w2_rdstr = std::move(w2_diss);
        auto uw_rdstr = std::move(uw_diss);

        calc_tke_budget_rdstr(
                u2_rdstr->fld.data(), v2_rdstr->fld.data(), w2_rdstr->fld.data(), uw_rdstr->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(),
                w_prime->fld.data(), fields.sd.at("p")->fld.data(),
                umodel.data(), vmodel.data(),
                gd.dzi4.data(), gd.dzhi4.data(),
                gd.dxi, gd.dyi,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_mask_stats(m, "u2_rdstr", *u2_rdstr, no_offset, no_threshold);
        stats.calc_mask_stats(m, "v2_rdstr", *v2_rdstr, no_offset, no_threshold);
        stats.calc_mask_stats(m, "w2_rdstr", *w2_rdstr, no_offset, no_threshold);
        stats.calc_mask_stats(m, "uw_rdstr", *uw_rdstr, no_offset, no_threshold);

        // Release the tmp arrays that are still in use.
        fields.release_tmp(uz);
        fields.release_tmp(wz);
        fields.release_tmp(u2_rdstr);
        fields.release_tmp(v2_rdstr);
        fields.release_tmp(w2_rdstr);
        fields.release_tmp(tke_diss);
        fields.release_tmp(uw_rdstr);

        // Calculate the buoyancy term of the TKE budget.
        if (thermo.get_switch() != "0")
        {
            auto b = fields.get_tmp();

            // Compute the buoyancy, cyclic is true, and stat is true.
            thermo.get_thermo_field(*b, "b", true, true);

            field3d_operators.calc_mean_profile(b->fld_mean.data(), b->fld.data());
            field3d_operators.calc_mean_profile(fields.sd.at("p")->fld_mean.data(), b->fld.data());

            auto w2_buoy  = fields.get_tmp();
            auto tke_buoy = fields.get_tmp();
            auto uw_buoy  = fields.get_tmp();

            calc_tke_budget_buoy(
                    w2_buoy->fld.data(), tke_buoy->fld.data(), uw_buoy->fld.data(),
                    fields.mp.at("u")->fld.data(), w_prime->fld.data(), b->fld.data(),
                    umodel.data(), b->fld_mean.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

            stats.calc_mask_stats(m, "w2_buoy" , *w2_buoy , no_offset, no_threshold);
            stats.calc_mask_stats(m, "tke_buoy", *tke_buoy, no_offset, no_threshold);
            stats.calc_mask_stats(m, "uw_buoy" , *uw_buoy , no_offset, no_threshold);

            auto b2_shear = std::move(w2_buoy);
            auto b2_turb = std::move(tke_buoy);
            auto b2_visc = std::move(uw_buoy);
            auto b2_diss = fields.get_tmp();

            calc_b2_budget(
                    b2_shear->fld.data(), b2_turb->fld.data(), b2_visc->fld.data(), b2_diss->fld.data(),
                    w_prime->fld.data(), b->fld.data(),
                    b->fld_mean.data(),
                    gd.dzi4.data(), gd.dzhi4.data(),
                    gd.dxi, gd.dyi,
                    thermo.get_buoyancy_diffusivity(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

            stats.calc_mask_stats(m, "b2_shear", *b2_shear, no_offset, no_threshold);
            stats.calc_mask_stats(m, "b2_turb" , *b2_turb , no_offset, no_threshold);
            stats.calc_mask_stats(m, "b2_visc" , *b2_visc , no_offset, no_threshold);
            stats.calc_mask_stats(m, "b2_diss" , *b2_diss , no_offset, no_threshold);

            auto bw_shear = std::move(b2_shear);
            auto bw_turb  = std::move(b2_turb);
            auto bw_visc  = std::move(b2_visc);
            auto bz       = std::move(b2_diss);

            calc_bw_budget_shear_turb_visc(
                    bw_shear->fld.data(), bw_turb->fld.data(), bw_visc->fld.data(),
                    bz->fld.data(),
                    w_prime->fld.data(), fields.sd.at("p")->fld.data(), b->fld.data(),
                    fields.sd.at("p")->fld_mean.data(), b->fld_mean.data(),
                    gd.dzi4.data(), gd.dzhi4.data(),
                    gd.dxi, gd.dyi, gd.dzhi4bot, gd.dzhi4top,
                    thermo.get_buoyancy_diffusivity(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.jcells, gd.ijcells);

            stats.calc_mask_stats(m, "bw_shear", *bw_shear, no_offset, no_threshold);
            stats.calc_mask_stats(m, "bw_turb" , *bw_turb , no_offset, no_threshold);
            stats.calc_mask_stats(m, "bw_visc" , *bw_visc , no_offset, no_threshold);

            auto bw_buoy  = std::move(bw_shear);
            auto bw_rdstr = std::move(bw_turb);
            auto bw_diss  = std::move(bw_visc);
            auto bw_pres  = fields.get_tmp();

            calc_bw_budget_buoy_rdstr_diss_pres(
                    bw_buoy->fld.data(), bw_rdstr->fld.data(), bw_diss->fld.data(), bw_pres->fld.data(),
                    bz->fld.data(),
                    w_prime->fld.data(), fields.sd.at("p")->fld.data(), b->fld.data(),
                    fields.sd.at("p")->fld_mean.data(), b->fld_mean.data(),
                    gd.dzi4.data(), gd.dzhi4.data(),
                    gd.dxi, gd.dyi, gd.dzhi4bot, gd.dzhi4top,
                    thermo.get_buoyancy_diffusivity(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.jcells, gd.ijcells);

            stats.calc_mask_stats(m, "bw_buoy" , *bw_buoy , no_offset, no_threshold);
            stats.calc_mask_stats(m, "bw_rdstr", *bw_rdstr, no_offset, no_threshold);
            stats.calc_mask_stats(m, "bw_diss" , *bw_diss , no_offset, no_threshold);
            stats.calc_mask_stats(m, "bw_pres" , *bw_pres , no_offset, no_threshold);

            fields.release_tmp(bw_buoy);
            fields.release_tmp(bw_rdstr);
            fields.release_tmp(bw_diss);
            fields.release_tmp(bw_pres);

            auto b_sort = std::move(bz);

            // Calculate the sorted buoyancy profile.
            // CvH: This does not work out well with the masking.
            calc_sorted_prof(
                    b->fld.data(), b_sort->fld.data(), b_sort->fld_mean.data(),
                    gd.z.data(), gd.dz.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells,
                    gd.itot, gd.jtot, gd.nmax,
                    grid.get_spatial_order(), master);

            stats.set_prof("b_sort", b_sort->fld_mean);

            fields.release_tmp(b);
            fields.release_tmp(b_sort);
        }

        fields.release_tmp(w_prime);
    }
}


#ifdef FLOAT_SINGLE
template class Budget_4<float>;
#else
template class Budget_4<double>;
#endif
