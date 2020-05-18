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

#include <algorithm>
#include <cmath>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "advec_4.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"

namespace
{
    using namespace Finite_difference::O4;
}

template<typename TF>
Advec_4<TF>::Advec_4(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Advec<TF>(masterin, gridin, fieldsin, inputin)
{
}

template<typename TF>
Advec_4<TF>::~Advec_4() {}

namespace
{
    template<typename TF>
    TF calc_cfl(
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict dzi, const TF dx, const TF dy,
            const TF dt, Master& master,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        const TF dxi = TF(1.)/dx;
        const TF dyi = TF(1.)/dy;

        TF cfl = 0;

        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    cfl = std::max(cfl, std::abs(interp4c(u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2]))*dxi
                                      + std::abs(interp4c(v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2]))*dyi
                                      + std::abs(interp4c(w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2]))*dzi[k]);
                }

        master.max(&cfl, 1);

        cfl = cfl*dt;

        return cfl;
    }

    template<typename TF, bool dim3>
    void advec_u(
            TF* const restrict ut,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict dzi4, const TF dx, const TF dy,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int jj3 = 3*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;
        const int kk3 = 3*kk;

        const TF dxi = TF(1.)/dx;
        const TF dyi = TF(1.)/dy;

        // bottom boundary
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj1 + kstart*kk1;
                ut[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii3] + ci1<TF>*u[ijk-ii2] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk    ]) * (ci0<TF>*u[ijk-ii3] + ci1<TF>*u[ijk-ii2] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk    ]))
                           + cg1<TF>*((ci0<TF>*u[ijk-ii2] + ci1<TF>*u[ijk-ii1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+ii1]) * (ci0<TF>*u[ijk-ii2] + ci1<TF>*u[ijk-ii1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+ii1]))
                           + cg2<TF>*((ci0<TF>*u[ijk-ii1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii2]) * (ci0<TF>*u[ijk-ii1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii2]))
                           + cg3<TF>*((ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+ii1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii3]) * (ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+ii1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii3])) ) * dxi;

                if (dim3)
                {
                    ut[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-ii2-jj1] + ci1<TF>*v[ijk-ii1-jj1] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk+ii1-jj1]) * (ci0<TF>*u[ijk-jj3] + ci1<TF>*u[ijk-jj2] + ci2<TF>*u[ijk-jj1] + ci3<TF>*u[ijk    ]))
                               + cg1<TF>*((ci0<TF>*v[ijk-ii2    ] + ci1<TF>*v[ijk-ii1    ] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+ii1    ]) * (ci0<TF>*u[ijk-jj2] + ci1<TF>*u[ijk-jj1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+jj1]))
                               + cg2<TF>*((ci0<TF>*v[ijk-ii2+jj1] + ci1<TF>*v[ijk-ii1+jj1] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+ii1+jj1]) * (ci0<TF>*u[ijk-jj1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+jj1] + ci3<TF>*u[ijk+jj2]))
                               + cg3<TF>*((ci0<TF>*v[ijk-ii2+jj2] + ci1<TF>*v[ijk-ii1+jj2] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+ii1+jj2]) * (ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+jj1] + ci2<TF>*u[ijk+jj2] + ci3<TF>*u[ijk+jj3])) ) * dyi;
                }

                ut[ijk] -= ( cg0<TF>*((ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk+ii1-kk1]) * (bi0<TF>*u[ijk-kk2] + bi1<TF>*u[ijk-kk1] + bi2<TF>*u[ijk    ] + bi3<TF>*u[ijk+kk1]))
                           + cg1<TF>*((ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1    ]) * (ci0<TF>*u[ijk-kk2] + ci1<TF>*u[ijk-kk1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+kk1]))
                           + cg2<TF>*((ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+ii1+kk1]) * (ci0<TF>*u[ijk-kk1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+kk1] + ci3<TF>*u[ijk+kk2]))
                           + cg3<TF>*((ci0<TF>*w[ijk-ii2+kk2] + ci1<TF>*w[ijk-ii1+kk2] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+ii1+kk2]) * (ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+kk1] + ci2<TF>*u[ijk+kk2] + ci3<TF>*u[ijk+kk3])) )
                         * dzi4[kstart];
            }

        for (int k=kstart+1; k<kend-1; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    ut[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii3] + ci1<TF>*u[ijk-ii2] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk    ]) * (ci0<TF>*u[ijk-ii3] + ci1<TF>*u[ijk-ii2] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk    ]))
                               + cg1<TF>*((ci0<TF>*u[ijk-ii2] + ci1<TF>*u[ijk-ii1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+ii1]) * (ci0<TF>*u[ijk-ii2] + ci1<TF>*u[ijk-ii1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+ii1]))
                               + cg2<TF>*((ci0<TF>*u[ijk-ii1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii2]) * (ci0<TF>*u[ijk-ii1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii2]))
                               + cg3<TF>*((ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+ii1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii3]) * (ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+ii1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii3])) ) * dxi;

                    if (dim3)
                    {
                        ut[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-ii2-jj1] + ci1<TF>*v[ijk-ii1-jj1] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk+ii1-jj1]) * (ci0<TF>*u[ijk-jj3] + ci1<TF>*u[ijk-jj2] + ci2<TF>*u[ijk-jj1] + ci3<TF>*u[ijk    ]))
                                   + cg1<TF>*((ci0<TF>*v[ijk-ii2    ] + ci1<TF>*v[ijk-ii1    ] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+ii1    ]) * (ci0<TF>*u[ijk-jj2] + ci1<TF>*u[ijk-jj1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+jj1]))
                                   + cg2<TF>*((ci0<TF>*v[ijk-ii2+jj1] + ci1<TF>*v[ijk-ii1+jj1] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+ii1+jj1]) * (ci0<TF>*u[ijk-jj1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+jj1] + ci3<TF>*u[ijk+jj2]))
                                   + cg3<TF>*((ci0<TF>*v[ijk-ii2+jj2] + ci1<TF>*v[ijk-ii1+jj2] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+ii1+jj2]) * (ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+jj1] + ci2<TF>*u[ijk+jj2] + ci3<TF>*u[ijk+jj3])) ) * dyi;
                    }

                    ut[ijk] -= ( cg0<TF>*((ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk+ii1-kk1]) * (ci0<TF>*u[ijk-kk3] + ci1<TF>*u[ijk-kk2] + ci2<TF>*u[ijk-kk1] + ci3<TF>*u[ijk    ]))
                               + cg1<TF>*((ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1    ]) * (ci0<TF>*u[ijk-kk2] + ci1<TF>*u[ijk-kk1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+kk1]))
                               + cg2<TF>*((ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+ii1+kk1]) * (ci0<TF>*u[ijk-kk1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+kk1] + ci3<TF>*u[ijk+kk2]))
                               + cg3<TF>*((ci0<TF>*w[ijk-ii2+kk2] + ci1<TF>*w[ijk-ii1+kk2] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+ii1+kk2]) * (ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+kk1] + ci2<TF>*u[ijk+kk2] + ci3<TF>*u[ijk+kk3])) )
                             * dzi4[k];
                }

        // top boundary
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj1 + (kend-1)*kk1;
                ut[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii3] + ci1<TF>*u[ijk-ii2] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk    ]) * (ci0<TF>*u[ijk-ii3] + ci1<TF>*u[ijk-ii2] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk    ]))
                           + cg1<TF>*((ci0<TF>*u[ijk-ii2] + ci1<TF>*u[ijk-ii1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+ii1]) * (ci0<TF>*u[ijk-ii2] + ci1<TF>*u[ijk-ii1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+ii1]))
                           + cg2<TF>*((ci0<TF>*u[ijk-ii1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii2]) * (ci0<TF>*u[ijk-ii1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii2]))
                           + cg3<TF>*((ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+ii1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii3]) * (ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+ii1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii3])) ) * dxi;

                if (dim3)
                {
                    ut[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-ii2-jj1] + ci1<TF>*v[ijk-ii1-jj1] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk+ii1-jj1]) * (ci0<TF>*u[ijk-jj3] + ci1<TF>*u[ijk-jj2] + ci2<TF>*u[ijk-jj1] + ci3<TF>*u[ijk    ]))
                               + cg1<TF>*((ci0<TF>*v[ijk-ii2    ] + ci1<TF>*v[ijk-ii1    ] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+ii1    ]) * (ci0<TF>*u[ijk-jj2] + ci1<TF>*u[ijk-jj1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+jj1]))
                               + cg2<TF>*((ci0<TF>*v[ijk-ii2+jj1] + ci1<TF>*v[ijk-ii1+jj1] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+ii1+jj1]) * (ci0<TF>*u[ijk-jj1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+jj1] + ci3<TF>*u[ijk+jj2]))
                               + cg3<TF>*((ci0<TF>*v[ijk-ii2+jj2] + ci1<TF>*v[ijk-ii1+jj2] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+ii1+jj2]) * (ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+jj1] + ci2<TF>*u[ijk+jj2] + ci3<TF>*u[ijk+jj3])) ) * dyi;
                }

                ut[ijk] -= ( cg0<TF>*((ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk+ii1-kk1]) * (ci0<TF>*u[ijk-kk3] + ci1<TF>*u[ijk-kk2] + ci2<TF>*u[ijk-kk1] + ci3<TF>*u[ijk    ]))
                           + cg1<TF>*((ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1    ]) * (ci0<TF>*u[ijk-kk2] + ci1<TF>*u[ijk-kk1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+kk1]))
                           + cg2<TF>*((ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+ii1+kk1]) * (ci0<TF>*u[ijk-kk1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+kk1] + ci3<TF>*u[ijk+kk2]))
                           + cg3<TF>*((ci0<TF>*w[ijk-ii2+kk2] + ci1<TF>*w[ijk-ii1+kk2] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+ii1+kk2]) * (ti0<TF>*u[ijk-kk1] + ti1<TF>*u[ijk    ] + ti2<TF>*u[ijk+kk1] + ti3<TF>*u[ijk+kk2])) )
                         * dzi4[kend-1];
            }
    }

    template<typename TF, bool dim3>
    void advec_v(
            TF* const restrict vt,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict dzi4, const TF dx, const TF dy,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int jj3 = 3*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;
        const int kk3 = 3*kk;

        const TF dxi = TF(1.)/dx;
        const TF dyi = TF(1.)/dy;

        // bottom boundary
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj1 + kstart*kk1;
                vt[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii1-jj2] + ci1<TF>*u[ijk-ii1-jj1] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk-ii1+jj1]) * (ci0<TF>*v[ijk-ii3] + ci1<TF>*v[ijk-ii2] + ci2<TF>*v[ijk-ii1] + ci3<TF>*v[ijk    ]))
                           + cg1<TF>*((ci0<TF>*u[ijk    -jj2] + ci1<TF>*u[ijk    -jj1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk    +jj1]) * (ci0<TF>*v[ijk-ii2] + ci1<TF>*v[ijk-ii1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+ii1]))
                           + cg2<TF>*((ci0<TF>*u[ijk+ii1-jj2] + ci1<TF>*u[ijk+ii1-jj1] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii1+jj1]) * (ci0<TF>*v[ijk-ii1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+ii1] + ci3<TF>*v[ijk+ii2]))
                           + cg3<TF>*((ci0<TF>*u[ijk+ii2-jj2] + ci1<TF>*u[ijk+ii2-jj1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii2+jj1]) * (ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+ii1] + ci2<TF>*v[ijk+ii2] + ci3<TF>*v[ijk+ii3])) ) * dxi;

                if (dim3)
                {
                    vt[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-jj3] + ci1<TF>*v[ijk-jj2] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk    ]) * (ci0<TF>*v[ijk-jj3] + ci1<TF>*v[ijk-jj2] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk    ]))
                               + cg1<TF>*((ci0<TF>*v[ijk-jj2] + ci1<TF>*v[ijk-jj1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+jj1]) * (ci0<TF>*v[ijk-jj2] + ci1<TF>*v[ijk-jj1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+jj1]))
                               + cg2<TF>*((ci0<TF>*v[ijk-jj1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj2]) * (ci0<TF>*v[ijk-jj1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj2]))
                               + cg3<TF>*((ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+jj1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj3]) * (ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+jj1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj3])) ) * dyi;
                }

                vt[ijk] -= ( cg0<TF>*((ci0<TF>*w[ijk-jj2-kk1] + ci1<TF>*w[ijk-jj1-kk1] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk+jj1-kk1]) * (bi0<TF>*v[ijk-kk2] + bi1<TF>*v[ijk-kk1] + bi2<TF>*v[ijk    ] + bi3<TF>*v[ijk+kk1]))
                           + cg1<TF>*((ci0<TF>*w[ijk-jj2    ] + ci1<TF>*w[ijk-jj1    ] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1    ]) * (ci0<TF>*v[ijk-kk2] + ci1<TF>*v[ijk-kk1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+kk1]))
                           + cg2<TF>*((ci0<TF>*w[ijk-jj2+kk1] + ci1<TF>*w[ijk-jj1+kk1] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+jj1+kk1]) * (ci0<TF>*v[ijk-kk1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+kk1] + ci3<TF>*v[ijk+kk2]))
                           + cg3<TF>*((ci0<TF>*w[ijk-jj2+kk2] + ci1<TF>*w[ijk-jj1+kk2] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+jj1+kk2]) * (ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+kk1] + ci2<TF>*v[ijk+kk2] + ci3<TF>*v[ijk+kk3])) )
                         * dzi4[kstart];
            }

        for (int k=kstart+1; k<kend-1; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    vt[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii1-jj2] + ci1<TF>*u[ijk-ii1-jj1] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk-ii1+jj1]) * (ci0<TF>*v[ijk-ii3] + ci1<TF>*v[ijk-ii2] + ci2<TF>*v[ijk-ii1] + ci3<TF>*v[ijk    ]))
                               + cg1<TF>*((ci0<TF>*u[ijk    -jj2] + ci1<TF>*u[ijk    -jj1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk    +jj1]) * (ci0<TF>*v[ijk-ii2] + ci1<TF>*v[ijk-ii1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+ii1]))
                               + cg2<TF>*((ci0<TF>*u[ijk+ii1-jj2] + ci1<TF>*u[ijk+ii1-jj1] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii1+jj1]) * (ci0<TF>*v[ijk-ii1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+ii1] + ci3<TF>*v[ijk+ii2]))
                               + cg3<TF>*((ci0<TF>*u[ijk+ii2-jj2] + ci1<TF>*u[ijk+ii2-jj1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii2+jj1]) * (ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+ii1] + ci2<TF>*v[ijk+ii2] + ci3<TF>*v[ijk+ii3])) ) * dxi;

                    if (dim3)
                    {
                        vt[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-jj3] + ci1<TF>*v[ijk-jj2] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk    ]) * (ci0<TF>*v[ijk-jj3] + ci1<TF>*v[ijk-jj2] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk    ]))
                                   + cg1<TF>*((ci0<TF>*v[ijk-jj2] + ci1<TF>*v[ijk-jj1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+jj1]) * (ci0<TF>*v[ijk-jj2] + ci1<TF>*v[ijk-jj1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+jj1]))
                                   + cg2<TF>*((ci0<TF>*v[ijk-jj1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj2]) * (ci0<TF>*v[ijk-jj1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj2]))
                                   + cg3<TF>*((ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+jj1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj3]) * (ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+jj1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj3])) ) * dyi;
                    }

                    vt[ijk] -= ( cg0<TF>*((ci0<TF>*w[ijk-jj2-kk1] + ci1<TF>*w[ijk-jj1-kk1] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk+jj1-kk1]) * (ci0<TF>*v[ijk-kk3] + ci1<TF>*v[ijk-kk2] + ci2<TF>*v[ijk-kk1] + ci3<TF>*v[ijk    ]))
                               + cg1<TF>*((ci0<TF>*w[ijk-jj2    ] + ci1<TF>*w[ijk-jj1    ] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1    ]) * (ci0<TF>*v[ijk-kk2] + ci1<TF>*v[ijk-kk1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+kk1]))
                               + cg2<TF>*((ci0<TF>*w[ijk-jj2+kk1] + ci1<TF>*w[ijk-jj1+kk1] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+jj1+kk1]) * (ci0<TF>*v[ijk-kk1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+kk1] + ci3<TF>*v[ijk+kk2]))
                               + cg3<TF>*((ci0<TF>*w[ijk-jj2+kk2] + ci1<TF>*w[ijk-jj1+kk2] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+jj1+kk2]) * (ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+kk1] + ci2<TF>*v[ijk+kk2] + ci3<TF>*v[ijk+kk3])) )
                             * dzi4[k];
                }

        // top boundary
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj1 + (kend-1)*kk1;
                vt[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii1-jj2] + ci1<TF>*u[ijk-ii1-jj1] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk-ii1+jj1]) * (ci0<TF>*v[ijk-ii3] + ci1<TF>*v[ijk-ii2] + ci2<TF>*v[ijk-ii1] + ci3<TF>*v[ijk    ]))
                           + cg1<TF>*((ci0<TF>*u[ijk    -jj2] + ci1<TF>*u[ijk    -jj1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk    +jj1]) * (ci0<TF>*v[ijk-ii2] + ci1<TF>*v[ijk-ii1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+ii1]))
                           + cg2<TF>*((ci0<TF>*u[ijk+ii1-jj2] + ci1<TF>*u[ijk+ii1-jj1] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii1+jj1]) * (ci0<TF>*v[ijk-ii1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+ii1] + ci3<TF>*v[ijk+ii2]))
                           + cg3<TF>*((ci0<TF>*u[ijk+ii2-jj2] + ci1<TF>*u[ijk+ii2-jj1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii2+jj1]) * (ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+ii1] + ci2<TF>*v[ijk+ii2] + ci3<TF>*v[ijk+ii3])) ) * dxi;

                if (dim3)
                {
                    vt[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-jj3] + ci1<TF>*v[ijk-jj2] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk    ]) * (ci0<TF>*v[ijk-jj3] + ci1<TF>*v[ijk-jj2] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk    ]))
                               + cg1<TF>*((ci0<TF>*v[ijk-jj2] + ci1<TF>*v[ijk-jj1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+jj1]) * (ci0<TF>*v[ijk-jj2] + ci1<TF>*v[ijk-jj1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+jj1]))
                               + cg2<TF>*((ci0<TF>*v[ijk-jj1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj2]) * (ci0<TF>*v[ijk-jj1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj2]))
                               + cg3<TF>*((ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+jj1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj3]) * (ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+jj1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj3])) ) * dyi;
                }

                vt[ijk] -= ( cg0<TF>*((ci0<TF>*w[ijk-jj2-kk1] + ci1<TF>*w[ijk-jj1-kk1] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk+jj1-kk1]) * (ci0<TF>*v[ijk-kk3] + ci1<TF>*v[ijk-kk2] + ci2<TF>*v[ijk-kk1] + ci3<TF>*v[ijk    ]))
                           + cg1<TF>*((ci0<TF>*w[ijk-jj2    ] + ci1<TF>*w[ijk-jj1    ] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1    ]) * (ci0<TF>*v[ijk-kk2] + ci1<TF>*v[ijk-kk1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+kk1]))
                           + cg2<TF>*((ci0<TF>*w[ijk-jj2+kk1] + ci1<TF>*w[ijk-jj1+kk1] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+jj1+kk1]) * (ci0<TF>*v[ijk-kk1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+kk1] + ci3<TF>*v[ijk+kk2]))
                           + cg3<TF>*((ci0<TF>*w[ijk-jj2+kk2] + ci1<TF>*w[ijk-jj1+kk2] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+jj1+kk2]) * (ti0<TF>*v[ijk-kk1] + ti1<TF>*v[ijk    ] + ti2<TF>*v[ijk+kk1] + ti3<TF>*v[ijk+kk2])) )
                         * dzi4[kend-1];
            }
    }

    template<typename TF, bool dim3>
    void advec_w(
            TF* const restrict wt,
            const TF* const restrict u, const TF* const restrict v, TF* const restrict w,
            const TF* const restrict dzhi4, const TF dx, const TF dy,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int jj3 = 3*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;
        const int kk3 = 3*kk;

        const TF dxi = TF(1.)/dx;
        const TF dyi = TF(1.)/dy;

        // bottom boundary
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj1 + (kstart+1)*kk1;
                wt[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii1-kk2] + ci1<TF>*u[ijk-ii1-kk1] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk-ii1+kk1]) * (ci0<TF>*w[ijk-ii3] + ci1<TF>*w[ijk-ii2] + ci2<TF>*w[ijk-ii1] + ci3<TF>*w[ijk    ]))
                           + cg1<TF>*((ci0<TF>*u[ijk    -kk2] + ci1<TF>*u[ijk    -kk1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk    +kk1]) * (ci0<TF>*w[ijk-ii2] + ci1<TF>*w[ijk-ii1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1]))
                           + cg2<TF>*((ci0<TF>*u[ijk+ii1-kk2] + ci1<TF>*u[ijk+ii1-kk1] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii1+kk1]) * (ci0<TF>*w[ijk-ii1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+ii1] + ci3<TF>*w[ijk+ii2]))
                           + cg3<TF>*((ci0<TF>*u[ijk+ii2-kk2] + ci1<TF>*u[ijk+ii2-kk1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii2+kk1]) * (ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+ii1] + ci2<TF>*w[ijk+ii2] + ci3<TF>*w[ijk+ii3])) ) * dxi;

                if (dim3)
                {
                    wt[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-jj1-kk2] + ci1<TF>*v[ijk-jj1-kk1] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk-jj1+kk1]) * (ci0<TF>*w[ijk-jj3] + ci1<TF>*w[ijk-jj2] + ci2<TF>*w[ijk-jj1] + ci3<TF>*w[ijk    ]))
                               + cg1<TF>*((ci0<TF>*v[ijk    -kk2] + ci1<TF>*v[ijk    -kk1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk    +kk1]) * (ci0<TF>*w[ijk-jj2] + ci1<TF>*w[ijk-jj1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1]))
                               + cg2<TF>*((ci0<TF>*v[ijk+jj1-kk2] + ci1<TF>*v[ijk+jj1-kk1] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj1+kk1]) * (ci0<TF>*w[ijk-jj1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+jj1] + ci3<TF>*w[ijk+jj2]))
                               + cg3<TF>*((ci0<TF>*v[ijk+jj2-kk2] + ci1<TF>*v[ijk+jj2-kk1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj2+kk1]) * (ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+jj1] + ci2<TF>*w[ijk+jj2] + ci3<TF>*w[ijk+jj3])) ) * dyi;
                }

                wt[ijk] -= ( cg0<TF>*((bi0<TF>*w[ijk-kk2] + bi1<TF>*w[ijk-kk1] + bi2<TF>*w[ijk    ] + bi3<TF>*w[ijk+kk1]) * (bi0<TF>*w[ijk-kk2] + bi1<TF>*w[ijk-kk1] + bi2<TF>*w[ijk    ] + bi3<TF>*w[ijk+kk1]))
                           + cg1<TF>*((ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1]) * (ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1]))
                           + cg2<TF>*((ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2]) * (ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2]))
                           + cg3<TF>*((ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3]) * (ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3])) )
                         * dzhi4[kstart+1];
            }

        for (int k=kstart+2; k<kend-1; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    wt[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii1-kk2] + ci1<TF>*u[ijk-ii1-kk1] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk-ii1+kk1]) * (ci0<TF>*w[ijk-ii3] + ci1<TF>*w[ijk-ii2] + ci2<TF>*w[ijk-ii1] + ci3<TF>*w[ijk    ]))
                               + cg1<TF>*((ci0<TF>*u[ijk    -kk2] + ci1<TF>*u[ijk    -kk1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk    +kk1]) * (ci0<TF>*w[ijk-ii2] + ci1<TF>*w[ijk-ii1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1]))
                               + cg2<TF>*((ci0<TF>*u[ijk+ii1-kk2] + ci1<TF>*u[ijk+ii1-kk1] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii1+kk1]) * (ci0<TF>*w[ijk-ii1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+ii1] + ci3<TF>*w[ijk+ii2]))
                               + cg3<TF>*((ci0<TF>*u[ijk+ii2-kk2] + ci1<TF>*u[ijk+ii2-kk1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii2+kk1]) * (ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+ii1] + ci2<TF>*w[ijk+ii2] + ci3<TF>*w[ijk+ii3])) ) * dxi;

                    if (dim3)
                    {
                        wt[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-jj1-kk2] + ci1<TF>*v[ijk-jj1-kk1] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk-jj1+kk1]) * (ci0<TF>*w[ijk-jj3] + ci1<TF>*w[ijk-jj2] + ci2<TF>*w[ijk-jj1] + ci3<TF>*w[ijk    ]))
                                   + cg1<TF>*((ci0<TF>*v[ijk    -kk2] + ci1<TF>*v[ijk    -kk1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk    +kk1]) * (ci0<TF>*w[ijk-jj2] + ci1<TF>*w[ijk-jj1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1]))
                                   + cg2<TF>*((ci0<TF>*v[ijk+jj1-kk2] + ci1<TF>*v[ijk+jj1-kk1] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj1+kk1]) * (ci0<TF>*w[ijk-jj1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+jj1] + ci3<TF>*w[ijk+jj2]))
                                   + cg3<TF>*((ci0<TF>*v[ijk+jj2-kk2] + ci1<TF>*v[ijk+jj2-kk1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj2+kk1]) * (ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+jj1] + ci2<TF>*w[ijk+jj2] + ci3<TF>*w[ijk+jj3])) ) * dyi;
                    }

                    wt[ijk] -= ( cg0<TF>*((ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ]) * (ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ]))
                               + cg1<TF>*((ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1]) * (ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1]))
                               + cg2<TF>*((ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2]) * (ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2]))
                               + cg3<TF>*((ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3]) * (ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3])) )
                             * dzhi4[k];
                }

        // top boundary
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj1 + (kend-1)*kk1;
                wt[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii1-kk2] + ci1<TF>*u[ijk-ii1-kk1] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk-ii1+kk1]) * (ci0<TF>*w[ijk-ii3] + ci1<TF>*w[ijk-ii2] + ci2<TF>*w[ijk-ii1] + ci3<TF>*w[ijk    ]))
                           + cg1<TF>*((ci0<TF>*u[ijk    -kk2] + ci1<TF>*u[ijk    -kk1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk    +kk1]) * (ci0<TF>*w[ijk-ii2] + ci1<TF>*w[ijk-ii1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1]))
                           + cg2<TF>*((ci0<TF>*u[ijk+ii1-kk2] + ci1<TF>*u[ijk+ii1-kk1] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii1+kk1]) * (ci0<TF>*w[ijk-ii1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+ii1] + ci3<TF>*w[ijk+ii2]))
                           + cg3<TF>*((ci0<TF>*u[ijk+ii2-kk2] + ci1<TF>*u[ijk+ii2-kk1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii2+kk1]) * (ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+ii1] + ci2<TF>*w[ijk+ii2] + ci3<TF>*w[ijk+ii3])) ) * dxi;

                if (dim3)
                {
                    wt[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-jj1-kk2] + ci1<TF>*v[ijk-jj1-kk1] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk-jj1+kk1]) * (ci0<TF>*w[ijk-jj3] + ci1<TF>*w[ijk-jj2] + ci2<TF>*w[ijk-jj1] + ci3<TF>*w[ijk    ]))
                               + cg1<TF>*((ci0<TF>*v[ijk    -kk2] + ci1<TF>*v[ijk    -kk1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk    +kk1]) * (ci0<TF>*w[ijk-jj2] + ci1<TF>*w[ijk-jj1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1]))
                               + cg2<TF>*((ci0<TF>*v[ijk+jj1-kk2] + ci1<TF>*v[ijk+jj1-kk1] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj1+kk1]) * (ci0<TF>*w[ijk-jj1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+jj1] + ci3<TF>*w[ijk+jj2]))
                               + cg3<TF>*((ci0<TF>*v[ijk+jj2-kk2] + ci1<TF>*v[ijk+jj2-kk1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj2+kk1]) * (ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+jj1] + ci2<TF>*w[ijk+jj2] + ci3<TF>*w[ijk+jj3])) ) * dyi;
                }

                wt[ijk] -= ( cg0<TF>*((ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ]) * (ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ]))
                           + cg1<TF>*((ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1]) * (ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1]))
                           + cg2<TF>*((ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2]) * (ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2]))
                           + cg3<TF>*((ti0<TF>*w[ijk-kk1] + ti1<TF>*w[ijk    ] + ti2<TF>*w[ijk+kk1] + ti3<TF>*w[ijk+kk2]) * (ti0<TF>*w[ijk-kk1] + ti1<TF>*w[ijk    ] + ti2<TF>*w[ijk+kk1] + ti3<TF>*w[ijk+kk2])) )
                         * dzhi4[kend-1];
            }
    }

    template<typename TF, bool dim3>
    void advec_s(
            TF* const restrict st, const TF* const restrict s,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict dzi4, const TF dx, const TF dy,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int jj3 = 3*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;
        const int kk3 = 3*kk;

        const TF dxi = TF(1.)/dx;
        const TF dyi = TF(1.)/dy;

        // bottom boundary
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj1 + kstart*kk1;
                st[ijk] -= ( cg0<TF>*(u[ijk-ii1] * (ci0<TF>*s[ijk-ii3] + ci1<TF>*s[ijk-ii2] + ci2<TF>*s[ijk-ii1] + ci3<TF>*s[ijk    ]))
                           + cg1<TF>*(u[ijk    ] * (ci0<TF>*s[ijk-ii2] + ci1<TF>*s[ijk-ii1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+ii1]))
                           + cg2<TF>*(u[ijk+ii1] * (ci0<TF>*s[ijk-ii1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+ii1] + ci3<TF>*s[ijk+ii2]))
                           + cg3<TF>*(u[ijk+ii2] * (ci0<TF>*s[ijk    ] + ci1<TF>*s[ijk+ii1] + ci2<TF>*s[ijk+ii2] + ci3<TF>*s[ijk+ii3])) ) * dxi;

                if (dim3)
                {
                    st[ijk] -= ( cg0<TF>*(v[ijk-jj1] * (ci0<TF>*s[ijk-jj3] + ci1<TF>*s[ijk-jj2] + ci2<TF>*s[ijk-jj1] + ci3<TF>*s[ijk    ]))
                               + cg1<TF>*(v[ijk    ] * (ci0<TF>*s[ijk-jj2] + ci1<TF>*s[ijk-jj1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+jj1]))
                               + cg2<TF>*(v[ijk+jj1] * (ci0<TF>*s[ijk-jj1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+jj1] + ci3<TF>*s[ijk+jj2]))
                               + cg3<TF>*(v[ijk+jj2] * (ci0<TF>*s[ijk    ] + ci1<TF>*s[ijk+jj1] + ci2<TF>*s[ijk+jj2] + ci3<TF>*s[ijk+jj3])) ) * dyi;
                }

                st[ijk] -= ( cg0<TF>*(w[ijk-kk1] * (bi0<TF>*s[ijk-kk2] + bi1<TF>*s[ijk-kk1] + bi2<TF>*s[ijk    ] + bi3<TF>*s[ijk+kk1]))
                           + cg1<TF>*(w[ijk    ] * (ci0<TF>*s[ijk-kk2] + ci1<TF>*s[ijk-kk1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+kk1]))
                           + cg2<TF>*(w[ijk+kk1] * (ci0<TF>*s[ijk-kk1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+kk1] + ci3<TF>*s[ijk+kk2]))
                           + cg3<TF>*(w[ijk+kk2] * (ci0<TF>*s[ijk    ] + ci1<TF>*s[ijk+kk1] + ci2<TF>*s[ijk+kk2] + ci3<TF>*s[ijk+kk3])) )
                         * dzi4[kstart];
            }

        for (int k=kstart+1; k<kend-1; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    st[ijk] -= ( cg0<TF>*(u[ijk-ii1] * (ci0<TF>*s[ijk-ii3] + ci1<TF>*s[ijk-ii2] + ci2<TF>*s[ijk-ii1] + ci3<TF>*s[ijk    ]))
                               + cg1<TF>*(u[ijk    ] * (ci0<TF>*s[ijk-ii2] + ci1<TF>*s[ijk-ii1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+ii1]))
                               + cg2<TF>*(u[ijk+ii1] * (ci0<TF>*s[ijk-ii1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+ii1] + ci3<TF>*s[ijk+ii2]))
                               + cg3<TF>*(u[ijk+ii2] * (ci0<TF>*s[ijk    ] + ci1<TF>*s[ijk+ii1] + ci2<TF>*s[ijk+ii2] + ci3<TF>*s[ijk+ii3])) ) * dxi;

                    if (dim3)
                    {
                        st[ijk] -= ( cg0<TF>*(v[ijk-jj1] * (ci0<TF>*s[ijk-jj3] + ci1<TF>*s[ijk-jj2] + ci2<TF>*s[ijk-jj1] + ci3<TF>*s[ijk    ]))
                                   + cg1<TF>*(v[ijk    ] * (ci0<TF>*s[ijk-jj2] + ci1<TF>*s[ijk-jj1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+jj1]))
                                   + cg2<TF>*(v[ijk+jj1] * (ci0<TF>*s[ijk-jj1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+jj1] + ci3<TF>*s[ijk+jj2]))
                                   + cg3<TF>*(v[ijk+jj2] * (ci0<TF>*s[ijk    ] + ci1<TF>*s[ijk+jj1] + ci2<TF>*s[ijk+jj2] + ci3<TF>*s[ijk+jj3])) ) * dyi;
                    }

                    st[ijk] -= ( cg0<TF>*(w[ijk-kk1] * (ci0<TF>*s[ijk-kk3] + ci1<TF>*s[ijk-kk2] + ci2<TF>*s[ijk-kk1] + ci3<TF>*s[ijk    ]))
                               + cg1<TF>*(w[ijk    ] * (ci0<TF>*s[ijk-kk2] + ci1<TF>*s[ijk-kk1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+kk1]))
                               + cg2<TF>*(w[ijk+kk1] * (ci0<TF>*s[ijk-kk1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+kk1] + ci3<TF>*s[ijk+kk2]))
                               + cg3<TF>*(w[ijk+kk2] * (ci0<TF>*s[ijk    ] + ci1<TF>*s[ijk+kk1] + ci2<TF>*s[ijk+kk2] + ci3<TF>*s[ijk+kk3])) )
                             * dzi4[k];
                }

        // top boundary
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ijk = i + j*jj1 + (kend-1)*kk1;
                st[ijk] -= ( cg0<TF>*(u[ijk-ii1] * (ci0<TF>*s[ijk-ii3] + ci1<TF>*s[ijk-ii2] + ci2<TF>*s[ijk-ii1] + ci3<TF>*s[ijk    ]))
                           + cg1<TF>*(u[ijk    ] * (ci0<TF>*s[ijk-ii2] + ci1<TF>*s[ijk-ii1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+ii1]))
                           + cg2<TF>*(u[ijk+ii1] * (ci0<TF>*s[ijk-ii1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+ii1] + ci3<TF>*s[ijk+ii2]))
                           + cg3<TF>*(u[ijk+ii2] * (ci0<TF>*s[ijk    ] + ci1<TF>*s[ijk+ii1] + ci2<TF>*s[ijk+ii2] + ci3<TF>*s[ijk+ii3])) ) * dxi;

                if (dim3)
                {
                    st[ijk] -= ( cg0<TF>*(v[ijk-jj1] * (ci0<TF>*s[ijk-jj3] + ci1<TF>*s[ijk-jj2] + ci2<TF>*s[ijk-jj1] + ci3<TF>*s[ijk    ]))
                               + cg1<TF>*(v[ijk    ] * (ci0<TF>*s[ijk-jj2] + ci1<TF>*s[ijk-jj1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+jj1]))
                               + cg2<TF>*(v[ijk+jj1] * (ci0<TF>*s[ijk-jj1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+jj1] + ci3<TF>*s[ijk+jj2]))
                               + cg3<TF>*(v[ijk+jj2] * (ci0<TF>*s[ijk    ] + ci1<TF>*s[ijk+jj1] + ci2<TF>*s[ijk+jj2] + ci3<TF>*s[ijk+jj3])) ) * dyi;
                }

                st[ijk] -= ( cg0<TF>*(w[ijk-kk1] * (ci0<TF>*s[ijk-kk3] + ci1<TF>*s[ijk-kk2] + ci2<TF>*s[ijk-kk1] + ci3<TF>*s[ijk    ]))
                           + cg1<TF>*(w[ijk    ] * (ci0<TF>*s[ijk-kk2] + ci1<TF>*s[ijk-kk1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+kk1]))
                           + cg2<TF>*(w[ijk+kk1] * (ci0<TF>*s[ijk-kk1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+kk1] + ci3<TF>*s[ijk+kk2]))
                           + cg3<TF>*(w[ijk+kk2] * (ti0<TF>*s[ijk-kk1] + ti1<TF>*s[ijk    ] + ti2<TF>*s[ijk+kk1] + ti3<TF>*s[ijk+kk2])) )
                         * dzi4[kend-1];
            }
        }

    template<typename TF>
    void advec_flux_u(
            TF* const restrict st, const TF* const restrict s, const TF* const restrict w,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] = (ci0<TF>*w[ijk-ii2] + ci1<TF>*w[ijk-ii1] + ci2<TF>*w[ijk] + ci3<TF>*w[ijk+ii1]) * (ci0<TF>*s[ijk-kk2] + ci1<TF>*s[ijk-kk1] + ci2<TF>*s[ijk] + ci3<TF>*s[ijk+kk1]);
                }
    }

    template<typename TF>
    void advec_flux_v(
            TF* const restrict st, const TF* const restrict s, const TF* const restrict w,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] = (ci0<TF>*w[ijk-jj2] + ci1<TF>*w[ijk-jj1] + ci2<TF>*w[ijk] + ci3<TF>*w[ijk+jj1]) * (ci0<TF>*s[ijk-kk2] + ci1<TF>*s[ijk-kk1] + ci2<TF>*s[ijk] + ci3<TF>*s[ijk+kk1]);
                }
    }

    template<typename TF>
    void advec_flux_s(
            TF* const restrict st, const TF* const restrict s, const TF* const restrict w,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        for (int k=kstart; k<kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] = w[ijk] * (ci0<TF>*s[ijk-kk2] + ci1<TF>*s[ijk-kk1] + ci2<TF>*s[ijk] + ci3<TF>*s[ijk+kk1]);
                }
    }
}

template<typename TF>
void Advec_4<TF>::create(Stats<TF>& stats)
{
    stats.add_tendency(*fields.mt.at("u"), "z", tend_name, tend_longname);
    stats.add_tendency(*fields.mt.at("v"), "z", tend_name, tend_longname);
    stats.add_tendency(*fields.mt.at("w"), "zh", tend_name, tend_longname);
    for (auto it : fields.st)
        stats.add_tendency(*it.second, "z", tend_name, tend_longname);
}

#ifndef USECUDA
template<typename TF>
double Advec_4<TF>::get_cfl(double dt)
{
    auto& gd = grid.get_grid_data();
    TF cfl = calc_cfl<TF>(fields.mp.at("u")->fld.data(),fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                          gd.dzi.data(), gd.dx, gd.dy,
                          dt, master,
                          gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                          gd.icells, gd.ijcells);

    // CFL is kept in double precision for time stepping accuracy
    return static_cast<double>(cfl);
}

template<typename TF>
unsigned long Advec_4<TF>::get_time_limit(unsigned long idt, double dt)
{
    // Calculate cfl and prevent zero divisons.
    auto& gd = grid.get_grid_data();

    double cfl = calc_cfl<TF>(fields.mp.at("u")->fld.data(),fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dx, gd.dy,
            dt, master,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    cfl = std::max(cflmin, cfl);
    return idt * cflmax / cfl;
}

template<typename TF>
void Advec_4<TF>::exec(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    if (gd.jtot == 1)
    {
        advec_u<TF,0>(
                fields.mt.at("u")->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                gd.dzi4.data(), gd.dx, gd.dy,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        advec_v<TF,0>(
                fields.mt.at("v")->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                gd.dzi4.data(), gd.dx, gd.dy,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        advec_w<TF,0>(
                fields.mt.at("w")->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                gd.dzhi4.data(), gd.dx, gd.dy,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        for (auto& it : fields.st)
            advec_s<TF,0>(
                    it.second->fld.data(), fields.sp.at(it.first)->fld.data(),
                    fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                    gd.dzi4.data(), gd.dx, gd.dy,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
    }
    else
    {
        advec_u<TF,1>(
                fields.mt.at("u")->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                gd.dzi4.data(), gd.dx, gd.dy,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        advec_v<TF,1>(
                fields.mt.at("v")->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                gd.dzi4.data(), gd.dx, gd.dy,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        advec_w<TF,1>(
                fields.mt.at("w")->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                gd.dzhi4.data(), gd.dx, gd.dy,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        for (auto& it : fields.st)
            advec_s<TF,1>(
                    it.second->fld.data(), fields.sp.at(it.first)->fld.data(),
                    fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                    gd.dzi4.data(), gd.dx, gd.dy,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
    }

    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);
    for (auto it : fields.st)
        stats.calc_tend(*it.second, tend_name);
}
#endif

template<typename TF>
void Advec_4<TF>::get_advec_flux(Field3d<TF>& advec_flux, const Field3d<TF>& fld)
{
    auto& gd = grid.get_grid_data();

    if (fld.loc == gd.uloc)
    {
        advec_flux_u(
                advec_flux.fld.data(), fld.fld.data(), fields.mp.at("w")->fld.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }
    else if (fld.loc == gd.vloc)
    {
        advec_flux_v(
                advec_flux.fld.data(), fld.fld.data(), fields.mp.at("w")->fld.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }
    else if (fld.loc == gd.sloc)
    {
        advec_flux_s(
                advec_flux.fld.data(), fld.fld.data(), fields.mp.at("w")->fld.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
    }
    else
        throw std::runtime_error("Advec_2 cannot deliver flux field at that location");
}

template class Advec_4<double>;
template class Advec_4<float>;
