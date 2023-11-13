/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

#include "advec_4.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "tools.h"
#include "constants.h"
#include "finite_difference.h"
#include "field3d_operators.h"

using namespace Finite_difference::O4;

namespace
{
    template<typename TF> __global__
    void advec_u_g(TF* __restrict__ ut, const TF* __restrict__ u,
                   const TF* __restrict__ v,  const TF* __restrict__ w,
                   const TF* __restrict__ dzi4, const TF dxi, const TF dyi,
                   const int jj, const int kk,
                   const int istart, const int jstart, const int kstart,
                   const int iend,   const int jend,   const int kend)
        {
            const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
            const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
            const int k = blockIdx.z + kstart;

            const int ii1 = 1;
            const int ii2 = 2;
            const int ii3 = 3;
            const int jj1 = 1*jj;
            const int jj2 = 2*jj;
            const int jj3 = 3*jj;
            const int kk1 = 1*kk;
            const int kk2 = 2*kk;
            const int kk3 = 3*kk;

            if (i < iend && j < jend && k > kstart && k < kend-1)
            {
                const int ijk = i + j*jj + k*kk;
                ut[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii3] + ci1<TF>*u[ijk-ii2] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk    ]) * (ci0<TF>*u[ijk-ii3] + ci1<TF>*u[ijk-ii2] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk    ]))
                           + cg1<TF>*((ci0<TF>*u[ijk-ii2] + ci1<TF>*u[ijk-ii1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+ii1]) * (ci0<TF>*u[ijk-ii2] + ci1<TF>*u[ijk-ii1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+ii1]))
                           + cg2<TF>*((ci0<TF>*u[ijk-ii1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii2]) * (ci0<TF>*u[ijk-ii1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii2]))
                           + cg3<TF>*((ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+ii1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii3]) * (ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+ii1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii3])) ) * dxi;

                ut[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-ii2-jj1] + ci1<TF>*v[ijk-ii1-jj1] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk+ii1-jj1]) * (ci0<TF>*u[ijk-jj3] + ci1<TF>*u[ijk-jj2] + ci2<TF>*u[ijk-jj1] + ci3<TF>*u[ijk    ]))
                           + cg1<TF>*((ci0<TF>*v[ijk-ii2    ] + ci1<TF>*v[ijk-ii1    ] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+ii1    ]) * (ci0<TF>*u[ijk-jj2] + ci1<TF>*u[ijk-jj1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+jj1]))
                           + cg2<TF>*((ci0<TF>*v[ijk-ii2+jj1] + ci1<TF>*v[ijk-ii1+jj1] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+ii1+jj1]) * (ci0<TF>*u[ijk-jj1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+jj1] + ci3<TF>*u[ijk+jj2]))
                           + cg3<TF>*((ci0<TF>*v[ijk-ii2+jj2] + ci1<TF>*v[ijk-ii1+jj2] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+ii1+jj2]) * (ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+jj1] + ci2<TF>*u[ijk+jj2] + ci3<TF>*u[ijk+jj3])) ) * dyi;

                ut[ijk] -= ( cg0<TF>*((ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk+ii1-kk1]) * (ci0<TF>*u[ijk-kk3] + ci1<TF>*u[ijk-kk2] + ci2<TF>*u[ijk-kk1] + ci3<TF>*u[ijk    ]))
                           + cg1<TF>*((ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1    ]) * (ci0<TF>*u[ijk-kk2] + ci1<TF>*u[ijk-kk1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+kk1]))
                           + cg2<TF>*((ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+ii1+kk1]) * (ci0<TF>*u[ijk-kk1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+kk1] + ci3<TF>*u[ijk+kk2]))
                           + cg3<TF>*((ci0<TF>*w[ijk-ii2+kk2] + ci1<TF>*w[ijk-ii1+kk2] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+ii1+kk2]) * (ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+kk1] + ci2<TF>*u[ijk+kk2] + ci3<TF>*u[ijk+kk3])) ) * dzi4[k];
            }
        }

    template<typename TF, int loc> __global__
    void advec_u_boundary_g(TF* __restrict__ ut, const TF* __restrict__ u,
                            const TF* __restrict__ v,  const TF* __restrict__ w,
                            const TF* __restrict__ dzi4, const TF dxi, const TF dyi,
                            const int jj, const int kk,
                            const int istart, const int jstart, const int kstart,
                            const int iend,   const int jend,   const int kend)
        {
            const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
            const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

            const int ii1 = 1;
            const int ii2 = 2;
            const int ii3 = 3;
            const int jj1 = 1*jj;
            const int jj2 = 2*jj;
            const int jj3 = 3*jj;
            const int kk1 = 1*kk;
            const int kk2 = 2*kk;
            const int kk3 = 3*kk;

            if (i < iend && j < jend)
            {
                if (loc == 0)
                {
                    const int k = kstart;
                    const int ijk = i + j*jj + k*kk;

                    ut[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii3] + ci1<TF>*u[ijk-ii2] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk    ]) * (ci0<TF>*u[ijk-ii3] + ci1<TF>*u[ijk-ii2] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk    ]))
                               + cg1<TF>*((ci0<TF>*u[ijk-ii2] + ci1<TF>*u[ijk-ii1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+ii1]) * (ci0<TF>*u[ijk-ii2] + ci1<TF>*u[ijk-ii1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+ii1]))
                               + cg2<TF>*((ci0<TF>*u[ijk-ii1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii2]) * (ci0<TF>*u[ijk-ii1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii2]))
                               + cg3<TF>*((ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+ii1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii3]) * (ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+ii1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii3])) ) * dxi;

                    ut[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-ii2-jj1] + ci1<TF>*v[ijk-ii1-jj1] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk+ii1-jj1]) * (ci0<TF>*u[ijk-jj3] + ci1<TF>*u[ijk-jj2] + ci2<TF>*u[ijk-jj1] + ci3<TF>*u[ijk    ]))
                               + cg1<TF>*((ci0<TF>*v[ijk-ii2    ] + ci1<TF>*v[ijk-ii1    ] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+ii1    ]) * (ci0<TF>*u[ijk-jj2] + ci1<TF>*u[ijk-jj1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+jj1]))
                               + cg2<TF>*((ci0<TF>*v[ijk-ii2+jj1] + ci1<TF>*v[ijk-ii1+jj1] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+ii1+jj1]) * (ci0<TF>*u[ijk-jj1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+jj1] + ci3<TF>*u[ijk+jj2]))
                               + cg3<TF>*((ci0<TF>*v[ijk-ii2+jj2] + ci1<TF>*v[ijk-ii1+jj2] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+ii1+jj2]) * (ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+jj1] + ci2<TF>*u[ijk+jj2] + ci3<TF>*u[ijk+jj3])) ) * dyi;

                    ut[ijk] -= ( cg0<TF>*((ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk+ii1-kk1]) * (bi0<TF>*u[ijk-kk2] + bi1<TF>*u[ijk-kk1] + bi2<TF>*u[ijk    ] + bi3<TF>*u[ijk+kk1]))
                               + cg1<TF>*((ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1    ]) * (ci0<TF>*u[ijk-kk2] + ci1<TF>*u[ijk-kk1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+kk1]))
                               + cg2<TF>*((ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+ii1+kk1]) * (ci0<TF>*u[ijk-kk1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+kk1] + ci3<TF>*u[ijk+kk2]))
                               + cg3<TF>*((ci0<TF>*w[ijk-ii2+kk2] + ci1<TF>*w[ijk-ii1+kk2] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+ii1+kk2]) * (ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+kk1] + ci2<TF>*u[ijk+kk2] + ci3<TF>*u[ijk+kk3])) ) * dzi4[k];
                }
                else if (loc == 1)
                {
                    const int k = kend-1;
                    const int ijk = i + j*jj + k*kk;

                    ut[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii3] + ci1<TF>*u[ijk-ii2] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk    ]) * (ci0<TF>*u[ijk-ii3] + ci1<TF>*u[ijk-ii2] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk    ]))
                               + cg1<TF>*((ci0<TF>*u[ijk-ii2] + ci1<TF>*u[ijk-ii1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+ii1]) * (ci0<TF>*u[ijk-ii2] + ci1<TF>*u[ijk-ii1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+ii1]))
                               + cg2<TF>*((ci0<TF>*u[ijk-ii1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii2]) * (ci0<TF>*u[ijk-ii1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii2]))
                               + cg3<TF>*((ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+ii1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii3]) * (ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+ii1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii3])) ) * dxi;

                    ut[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-ii2-jj1] + ci1<TF>*v[ijk-ii1-jj1] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk+ii1-jj1]) * (ci0<TF>*u[ijk-jj3] + ci1<TF>*u[ijk-jj2] + ci2<TF>*u[ijk-jj1] + ci3<TF>*u[ijk    ]))
                               + cg1<TF>*((ci0<TF>*v[ijk-ii2    ] + ci1<TF>*v[ijk-ii1    ] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+ii1    ]) * (ci0<TF>*u[ijk-jj2] + ci1<TF>*u[ijk-jj1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+jj1]))
                               + cg2<TF>*((ci0<TF>*v[ijk-ii2+jj1] + ci1<TF>*v[ijk-ii1+jj1] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+ii1+jj1]) * (ci0<TF>*u[ijk-jj1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+jj1] + ci3<TF>*u[ijk+jj2]))
                               + cg3<TF>*((ci0<TF>*v[ijk-ii2+jj2] + ci1<TF>*v[ijk-ii1+jj2] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+ii1+jj2]) * (ci0<TF>*u[ijk    ] + ci1<TF>*u[ijk+jj1] + ci2<TF>*u[ijk+jj2] + ci3<TF>*u[ijk+jj3])) ) * dyi;

                    ut[ijk] -= ( cg0<TF>*((ci0<TF>*w[ijk-ii2-kk1] + ci1<TF>*w[ijk-ii1-kk1] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk+ii1-kk1]) * (ci0<TF>*u[ijk-kk3] + ci1<TF>*u[ijk-kk2] + ci2<TF>*u[ijk-kk1] + ci3<TF>*u[ijk    ]))
                               + cg1<TF>*((ci0<TF>*w[ijk-ii2    ] + ci1<TF>*w[ijk-ii1    ] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1    ]) * (ci0<TF>*u[ijk-kk2] + ci1<TF>*u[ijk-kk1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk+kk1]))
                               + cg2<TF>*((ci0<TF>*w[ijk-ii2+kk1] + ci1<TF>*w[ijk-ii1+kk1] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+ii1+kk1]) * (ci0<TF>*u[ijk-kk1] + ci1<TF>*u[ijk    ] + ci2<TF>*u[ijk+kk1] + ci3<TF>*u[ijk+kk2]))
                               + cg3<TF>*((ci0<TF>*w[ijk-ii2+kk2] + ci1<TF>*w[ijk-ii1+kk2] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+ii1+kk2]) * (ti0<TF>*u[ijk-kk1] + ti1<TF>*u[ijk    ] + ti2<TF>*u[ijk+kk1] + ti3<TF>*u[ijk+kk2])) ) * dzi4[k];
                }
            }
        }


    template<typename TF> __global__
    void advec_v_g(TF* __restrict__ vt, const TF* __restrict__ u,
                   const TF* __restrict__ v,  const TF* __restrict__ w,
                   const TF* __restrict__ dzi4, const TF dxi, const TF dyi,
                   const int jj, const int kk,
                   const int istart, const int jstart, const int kstart,
                   const int iend,   const int jend,   const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int jj3 = 3*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;
        const int kk3 = 3*kk;

        if (i < iend && j < jend && k > kstart && k < kend-1)
        {
            const int ijk = i + j*jj + k*kk;

            vt[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii1-jj2] + ci1<TF>*u[ijk-ii1-jj1] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk-ii1+jj1]) * (ci0<TF>*v[ijk-ii3] + ci1<TF>*v[ijk-ii2] + ci2<TF>*v[ijk-ii1] + ci3<TF>*v[ijk    ]))
                       + cg1<TF>*((ci0<TF>*u[ijk    -jj2] + ci1<TF>*u[ijk    -jj1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk    +jj1]) * (ci0<TF>*v[ijk-ii2] + ci1<TF>*v[ijk-ii1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+ii1]))
                       + cg2<TF>*((ci0<TF>*u[ijk+ii1-jj2] + ci1<TF>*u[ijk+ii1-jj1] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii1+jj1]) * (ci0<TF>*v[ijk-ii1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+ii1] + ci3<TF>*v[ijk+ii2]))
                       + cg3<TF>*((ci0<TF>*u[ijk+ii2-jj2] + ci1<TF>*u[ijk+ii2-jj1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii2+jj1]) * (ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+ii1] + ci2<TF>*v[ijk+ii2] + ci3<TF>*v[ijk+ii3])) ) * dxi;

            vt[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-jj3] + ci1<TF>*v[ijk-jj2] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk    ]) * (ci0<TF>*v[ijk-jj3] + ci1<TF>*v[ijk-jj2] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk    ]))
                       + cg1<TF>*((ci0<TF>*v[ijk-jj2] + ci1<TF>*v[ijk-jj1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+jj1]) * (ci0<TF>*v[ijk-jj2] + ci1<TF>*v[ijk-jj1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+jj1]))
                       + cg2<TF>*((ci0<TF>*v[ijk-jj1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj2]) * (ci0<TF>*v[ijk-jj1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj2]))
                       + cg3<TF>*((ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+jj1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj3]) * (ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+jj1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj3])) ) * dyi;

            vt[ijk] -= ( cg0<TF>*((ci0<TF>*w[ijk-jj2-kk1] + ci1<TF>*w[ijk-jj1-kk1] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk+jj1-kk1]) * (ci0<TF>*v[ijk-kk3] + ci1<TF>*v[ijk-kk2] + ci2<TF>*v[ijk-kk1] + ci3<TF>*v[ijk    ]))
                       + cg1<TF>*((ci0<TF>*w[ijk-jj2    ] + ci1<TF>*w[ijk-jj1    ] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1    ]) * (ci0<TF>*v[ijk-kk2] + ci1<TF>*v[ijk-kk1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+kk1]))
                       + cg2<TF>*((ci0<TF>*w[ijk-jj2+kk1] + ci1<TF>*w[ijk-jj1+kk1] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+jj1+kk1]) * (ci0<TF>*v[ijk-kk1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+kk1] + ci3<TF>*v[ijk+kk2]))
                       + cg3<TF>*((ci0<TF>*w[ijk-jj2+kk2] + ci1<TF>*w[ijk-jj1+kk2] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+jj1+kk2]) * (ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+kk1] + ci2<TF>*v[ijk+kk2] + ci3<TF>*v[ijk+kk3])) ) * dzi4[k];
        }
    }

    template<typename TF, int loc> __global__
    void advec_v_boundary_g(TF* __restrict__ vt, const TF* __restrict__ u,
                            const TF* __restrict__ v,  const TF* __restrict__ w,
                            const TF* __restrict__ dzi4, const TF dxi, const TF dyi,
                            const int jj, const int kk,
                            const int istart, const int jstart, const int kstart,
                            const int iend,   const int jend,   const int kend)
        {
            const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
            const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

            const int ii1 = 1;
            const int ii2 = 2;
            const int ii3 = 3;
            const int jj1 = 1*jj;
            const int jj2 = 2*jj;
            const int jj3 = 3*jj;
            const int kk1 = 1*kk;
            const int kk2 = 2*kk;
            const int kk3 = 3*kk;

            if (i < iend && j < jend)
            {
                if (loc == 0)
                {
                    const int k = kstart;
                    const int ijk = i + j*jj + k*kk;

                    vt[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii1-jj2] + ci1<TF>*u[ijk-ii1-jj1] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk-ii1+jj1]) * (ci0<TF>*v[ijk-ii3] + ci1<TF>*v[ijk-ii2] + ci2<TF>*v[ijk-ii1] + ci3<TF>*v[ijk    ]))
                               + cg1<TF>*((ci0<TF>*u[ijk    -jj2] + ci1<TF>*u[ijk    -jj1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk    +jj1]) * (ci0<TF>*v[ijk-ii2] + ci1<TF>*v[ijk-ii1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+ii1]))
                               + cg2<TF>*((ci0<TF>*u[ijk+ii1-jj2] + ci1<TF>*u[ijk+ii1-jj1] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii1+jj1]) * (ci0<TF>*v[ijk-ii1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+ii1] + ci3<TF>*v[ijk+ii2]))
                               + cg3<TF>*((ci0<TF>*u[ijk+ii2-jj2] + ci1<TF>*u[ijk+ii2-jj1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii2+jj1]) * (ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+ii1] + ci2<TF>*v[ijk+ii2] + ci3<TF>*v[ijk+ii3])) ) * dxi;

                    vt[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-jj3] + ci1<TF>*v[ijk-jj2] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk    ]) * (ci0<TF>*v[ijk-jj3] + ci1<TF>*v[ijk-jj2] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk    ]))
                               + cg1<TF>*((ci0<TF>*v[ijk-jj2] + ci1<TF>*v[ijk-jj1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+jj1]) * (ci0<TF>*v[ijk-jj2] + ci1<TF>*v[ijk-jj1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+jj1]))
                               + cg2<TF>*((ci0<TF>*v[ijk-jj1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj2]) * (ci0<TF>*v[ijk-jj1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj2]))
                               + cg3<TF>*((ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+jj1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj3]) * (ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+jj1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj3])) ) * dyi;

                    vt[ijk] -= ( cg0<TF>*((ci0<TF>*w[ijk-jj2-kk1] + ci1<TF>*w[ijk-jj1-kk1] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk+jj1-kk1]) * (bi0<TF>*v[ijk-kk2] + bi1<TF>*v[ijk-kk1] + bi2<TF>*v[ijk    ] + bi3<TF>*v[ijk+kk1]))
                               + cg1<TF>*((ci0<TF>*w[ijk-jj2    ] + ci1<TF>*w[ijk-jj1    ] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1    ]) * (ci0<TF>*v[ijk-kk2] + ci1<TF>*v[ijk-kk1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+kk1]))
                               + cg2<TF>*((ci0<TF>*w[ijk-jj2+kk1] + ci1<TF>*w[ijk-jj1+kk1] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+jj1+kk1]) * (ci0<TF>*v[ijk-kk1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+kk1] + ci3<TF>*v[ijk+kk2]))
                               + cg3<TF>*((ci0<TF>*w[ijk-jj2+kk2] + ci1<TF>*w[ijk-jj1+kk2] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+jj1+kk2]) * (ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+kk1] + ci2<TF>*v[ijk+kk2] + ci3<TF>*v[ijk+kk3])) ) * dzi4[k];
                }
                else if (loc == 1)
                {
                    const int k = kend-1;
                    const int ijk = i + j*jj + k*kk;

                    vt[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii1-jj2] + ci1<TF>*u[ijk-ii1-jj1] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk-ii1+jj1]) * (ci0<TF>*v[ijk-ii3] + ci1<TF>*v[ijk-ii2] + ci2<TF>*v[ijk-ii1] + ci3<TF>*v[ijk    ]))
                               + cg1<TF>*((ci0<TF>*u[ijk    -jj2] + ci1<TF>*u[ijk    -jj1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk    +jj1]) * (ci0<TF>*v[ijk-ii2] + ci1<TF>*v[ijk-ii1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+ii1]))
                               + cg2<TF>*((ci0<TF>*u[ijk+ii1-jj2] + ci1<TF>*u[ijk+ii1-jj1] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii1+jj1]) * (ci0<TF>*v[ijk-ii1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+ii1] + ci3<TF>*v[ijk+ii2]))
                               + cg3<TF>*((ci0<TF>*u[ijk+ii2-jj2] + ci1<TF>*u[ijk+ii2-jj1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii2+jj1]) * (ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+ii1] + ci2<TF>*v[ijk+ii2] + ci3<TF>*v[ijk+ii3])) ) * dxi;

                    vt[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-jj3] + ci1<TF>*v[ijk-jj2] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk    ]) * (ci0<TF>*v[ijk-jj3] + ci1<TF>*v[ijk-jj2] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk    ]))
                               + cg1<TF>*((ci0<TF>*v[ijk-jj2] + ci1<TF>*v[ijk-jj1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+jj1]) * (ci0<TF>*v[ijk-jj2] + ci1<TF>*v[ijk-jj1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+jj1]))
                               + cg2<TF>*((ci0<TF>*v[ijk-jj1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj2]) * (ci0<TF>*v[ijk-jj1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj2]))
                               + cg3<TF>*((ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+jj1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj3]) * (ci0<TF>*v[ijk    ] + ci1<TF>*v[ijk+jj1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj3])) ) * dyi;

                    vt[ijk] -= ( cg0<TF>*((ci0<TF>*w[ijk-jj2-kk1] + ci1<TF>*w[ijk-jj1-kk1] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk+jj1-kk1]) * (ci0<TF>*v[ijk-kk3] + ci1<TF>*v[ijk-kk2] + ci2<TF>*v[ijk-kk1] + ci3<TF>*v[ijk    ]))
                               + cg1<TF>*((ci0<TF>*w[ijk-jj2    ] + ci1<TF>*w[ijk-jj1    ] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1    ]) * (ci0<TF>*v[ijk-kk2] + ci1<TF>*v[ijk-kk1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk+kk1]))
                               + cg2<TF>*((ci0<TF>*w[ijk-jj2+kk1] + ci1<TF>*w[ijk-jj1+kk1] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+jj1+kk1]) * (ci0<TF>*v[ijk-kk1] + ci1<TF>*v[ijk    ] + ci2<TF>*v[ijk+kk1] + ci3<TF>*v[ijk+kk2]))
                               + cg3<TF>*((ci0<TF>*w[ijk-jj2+kk2] + ci1<TF>*w[ijk-jj1+kk2] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+jj1+kk2]) * (ti0<TF>*v[ijk-kk1] + ti1<TF>*v[ijk    ] + ti2<TF>*v[ijk+kk1] + ti3<TF>*v[ijk+kk2])) ) * dzi4[k];
                }
            }
        }

    template<typename TF> __global__
    void advec_w_g(TF* __restrict__ wt, const TF* __restrict__ u,
                   const TF* __restrict__ v,  const TF* __restrict__ w,
                   const TF* __restrict__ dzhi4, const TF dxi, const TF dyi,
                   const int jj, const int kk,
                   const int istart, const int jstart, const int kstart,
                   const int iend,   const int jend,   const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart + 1;

        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int jj3 = 3*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;
        const int kk3 = 3*kk;

        if(i < iend && j < jend && k > kstart+1 && k < kend-1)
        {
            const int ijk = i + j*jj + k*kk;

            wt[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii1-kk2] + ci1<TF>*u[ijk-ii1-kk1] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk-ii1+kk1]) * (ci0<TF>*w[ijk-ii3] + ci1<TF>*w[ijk-ii2] + ci2<TF>*w[ijk-ii1] + ci3<TF>*w[ijk    ]))
                       + cg1<TF>*((ci0<TF>*u[ijk    -kk2] + ci1<TF>*u[ijk    -kk1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk    +kk1]) * (ci0<TF>*w[ijk-ii2] + ci1<TF>*w[ijk-ii1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1]))
                       + cg2<TF>*((ci0<TF>*u[ijk+ii1-kk2] + ci1<TF>*u[ijk+ii1-kk1] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii1+kk1]) * (ci0<TF>*w[ijk-ii1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+ii1] + ci3<TF>*w[ijk+ii2]))
                       + cg3<TF>*((ci0<TF>*u[ijk+ii2-kk2] + ci1<TF>*u[ijk+ii2-kk1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii2+kk1]) * (ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+ii1] + ci2<TF>*w[ijk+ii2] + ci3<TF>*w[ijk+ii3])) ) * dxi;

            wt[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-jj1-kk2] + ci1<TF>*v[ijk-jj1-kk1] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk-jj1+kk1]) * (ci0<TF>*w[ijk-jj3] + ci1<TF>*w[ijk-jj2] + ci2<TF>*w[ijk-jj1] + ci3<TF>*w[ijk    ]))
                       + cg1<TF>*((ci0<TF>*v[ijk    -kk2] + ci1<TF>*v[ijk    -kk1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk    +kk1]) * (ci0<TF>*w[ijk-jj2] + ci1<TF>*w[ijk-jj1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1]))
                       + cg2<TF>*((ci0<TF>*v[ijk+jj1-kk2] + ci1<TF>*v[ijk+jj1-kk1] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj1+kk1]) * (ci0<TF>*w[ijk-jj1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+jj1] + ci3<TF>*w[ijk+jj2]))
                       + cg3<TF>*((ci0<TF>*v[ijk+jj2-kk2] + ci1<TF>*v[ijk+jj2-kk1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj2+kk1]) * (ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+jj1] + ci2<TF>*w[ijk+jj2] + ci3<TF>*w[ijk+jj3])) ) * dyi;

            wt[ijk] -= ( cg0<TF>*((ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ]) * (ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ]))
                       + cg1<TF>*((ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1]) * (ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1]))
                       + cg2<TF>*((ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2]) * (ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2]))
                       + cg3<TF>*((ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3]) * (ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3])) ) * dzhi4[k];
        }
    }

    template<typename TF, int loc> __global__
    void advec_w_boundary_g(TF* __restrict__ wt, const TF* __restrict__ u,
                            const TF* __restrict__ v,  const TF* __restrict__ w,
                            const TF* __restrict__ dzhi4, const TF dxi, const TF dyi,
                            const int jj, const int kk,
                            const int istart, const int jstart, const int kstart,
                            const int iend,   const int jend,   const int kend)
        {
            const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
            const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

            const int ii1 = 1;
            const int ii2 = 2;
            const int ii3 = 3;
            const int jj1 = 1*jj;
            const int jj2 = 2*jj;
            const int jj3 = 3*jj;
            const int kk1 = 1*kk;
            const int kk2 = 2*kk;
            const int kk3 = 3*kk;

            if (i < iend && j < jend)
            {
                if (loc == 0) //== kstart+1
                {
                    const int k = kstart+1;
                    const int ijk = i + j*jj + k*kk;

                    wt[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii1-kk2] + ci1<TF>*u[ijk-ii1-kk1] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk-ii1+kk1]) * (ci0<TF>*w[ijk-ii3] + ci1<TF>*w[ijk-ii2] + ci2<TF>*w[ijk-ii1] + ci3<TF>*w[ijk    ]))
                               + cg1<TF>*((ci0<TF>*u[ijk    -kk2] + ci1<TF>*u[ijk    -kk1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk    +kk1]) * (ci0<TF>*w[ijk-ii2] + ci1<TF>*w[ijk-ii1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1]))
                               + cg2<TF>*((ci0<TF>*u[ijk+ii1-kk2] + ci1<TF>*u[ijk+ii1-kk1] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii1+kk1]) * (ci0<TF>*w[ijk-ii1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+ii1] + ci3<TF>*w[ijk+ii2]))
                               + cg3<TF>*((ci0<TF>*u[ijk+ii2-kk2] + ci1<TF>*u[ijk+ii2-kk1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii2+kk1]) * (ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+ii1] + ci2<TF>*w[ijk+ii2] + ci3<TF>*w[ijk+ii3])) ) * dxi;

                    wt[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-jj1-kk2] + ci1<TF>*v[ijk-jj1-kk1] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk-jj1+kk1]) * (ci0<TF>*w[ijk-jj3] + ci1<TF>*w[ijk-jj2] + ci2<TF>*w[ijk-jj1] + ci3<TF>*w[ijk    ]))
                               + cg1<TF>*((ci0<TF>*v[ijk    -kk2] + ci1<TF>*v[ijk    -kk1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk    +kk1]) * (ci0<TF>*w[ijk-jj2] + ci1<TF>*w[ijk-jj1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1]))
                               + cg2<TF>*((ci0<TF>*v[ijk+jj1-kk2] + ci1<TF>*v[ijk+jj1-kk1] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj1+kk1]) * (ci0<TF>*w[ijk-jj1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+jj1] + ci3<TF>*w[ijk+jj2]))
                               + cg3<TF>*((ci0<TF>*v[ijk+jj2-kk2] + ci1<TF>*v[ijk+jj2-kk1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj2+kk1]) * (ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+jj1] + ci2<TF>*w[ijk+jj2] + ci3<TF>*w[ijk+jj3])) ) * dyi;

                    wt[ijk] -= ( cg0<TF>*((bi0<TF>*w[ijk-kk2] + bi1<TF>*w[ijk-kk1] + bi2<TF>*w[ijk    ] + bi3<TF>*w[ijk+kk1]) * (bi0<TF>*w[ijk-kk2] + bi1<TF>*w[ijk-kk1] + bi2<TF>*w[ijk    ] + bi3<TF>*w[ijk+kk1]))
                               + cg1<TF>*((ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1]) * (ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1]))
                               + cg2<TF>*((ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2]) * (ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2]))
                               + cg3<TF>*((ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3]) * (ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+kk1] + ci2<TF>*w[ijk+kk2] + ci3<TF>*w[ijk+kk3])) ) * dzhi4[k];
                }
                else if (loc == 1) //==kend-1
                {
                    const int k = kend-1;
                    const int ijk = i + j*jj + k*kk;

                    wt[ijk] -= ( cg0<TF>*((ci0<TF>*u[ijk-ii1-kk2] + ci1<TF>*u[ijk-ii1-kk1] + ci2<TF>*u[ijk-ii1] + ci3<TF>*u[ijk-ii1+kk1]) * (ci0<TF>*w[ijk-ii3] + ci1<TF>*w[ijk-ii2] + ci2<TF>*w[ijk-ii1] + ci3<TF>*w[ijk    ]))
                               + cg1<TF>*((ci0<TF>*u[ijk    -kk2] + ci1<TF>*u[ijk    -kk1] + ci2<TF>*u[ijk    ] + ci3<TF>*u[ijk    +kk1]) * (ci0<TF>*w[ijk-ii2] + ci1<TF>*w[ijk-ii1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+ii1]))
                               + cg2<TF>*((ci0<TF>*u[ijk+ii1-kk2] + ci1<TF>*u[ijk+ii1-kk1] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii1+kk1]) * (ci0<TF>*w[ijk-ii1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+ii1] + ci3<TF>*w[ijk+ii2]))
                               + cg3<TF>*((ci0<TF>*u[ijk+ii2-kk2] + ci1<TF>*u[ijk+ii2-kk1] + ci2<TF>*u[ijk+ii2] + ci3<TF>*u[ijk+ii2+kk1]) * (ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+ii1] + ci2<TF>*w[ijk+ii2] + ci3<TF>*w[ijk+ii3])) ) * dxi;

                    wt[ijk] -= ( cg0<TF>*((ci0<TF>*v[ijk-jj1-kk2] + ci1<TF>*v[ijk-jj1-kk1] + ci2<TF>*v[ijk-jj1] + ci3<TF>*v[ijk-jj1+kk1]) * (ci0<TF>*w[ijk-jj3] + ci1<TF>*w[ijk-jj2] + ci2<TF>*w[ijk-jj1] + ci3<TF>*w[ijk    ]))
                               + cg1<TF>*((ci0<TF>*v[ijk    -kk2] + ci1<TF>*v[ijk    -kk1] + ci2<TF>*v[ijk    ] + ci3<TF>*v[ijk    +kk1]) * (ci0<TF>*w[ijk-jj2] + ci1<TF>*w[ijk-jj1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+jj1]))
                               + cg2<TF>*((ci0<TF>*v[ijk+jj1-kk2] + ci1<TF>*v[ijk+jj1-kk1] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj1+kk1]) * (ci0<TF>*w[ijk-jj1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+jj1] + ci3<TF>*w[ijk+jj2]))
                               + cg3<TF>*((ci0<TF>*v[ijk+jj2-kk2] + ci1<TF>*v[ijk+jj2-kk1] + ci2<TF>*v[ijk+jj2] + ci3<TF>*v[ijk+jj2+kk1]) * (ci0<TF>*w[ijk    ] + ci1<TF>*w[ijk+jj1] + ci2<TF>*w[ijk+jj2] + ci3<TF>*w[ijk+jj3])) ) * dyi;

                    wt[ijk] -= ( cg0<TF>*((ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ]) * (ci0<TF>*w[ijk-kk3] + ci1<TF>*w[ijk-kk2] + ci2<TF>*w[ijk-kk1] + ci3<TF>*w[ijk    ]))
                               + cg1<TF>*((ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1]) * (ci0<TF>*w[ijk-kk2] + ci1<TF>*w[ijk-kk1] + ci2<TF>*w[ijk    ] + ci3<TF>*w[ijk+kk1]))
                               + cg2<TF>*((ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2]) * (ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk    ] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2]))
                               + cg3<TF>*((ti0<TF>*w[ijk-kk1] + ti1<TF>*w[ijk    ] + ti2<TF>*w[ijk+kk1] + ti3<TF>*w[ijk+kk2]) * (ti0<TF>*w[ijk-kk1] + ti1<TF>*w[ijk    ] + ti2<TF>*w[ijk+kk1] + ti3<TF>*w[ijk+kk2])) ) * dzhi4[k];
                }
            }
        }

    template<typename TF> __global__
    void advec_s_g(TF* __restrict__ st, const TF* __restrict__ s,
                   const TF* __restrict__ u,  const TF* __restrict__ v, const TF* __restrict__ w,
                   const TF* __restrict__ dzi4, const TF dxi, const TF dyi,
                   const int jj, const int kk,
                   const int istart, const int jstart, const int kstart,
                   const int iend,   const int jend,   const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int ii1 = 1;
        const int ii2 = 2;
        const int ii3 = 3;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int jj3 = 3*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;
        const int kk3 = 3*kk;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;

            if (k == kstart)
            {
                st[ijk] -= ( cg0<TF>*(u[ijk-ii1] * (ci0<TF>*s[ijk-ii3] + ci1<TF>*s[ijk-ii2] + ci2<TF>*s[ijk-ii1] + ci3<TF>*s[ijk    ]))
                           + cg1<TF>*(u[ijk    ] * (ci0<TF>*s[ijk-ii2] + ci1<TF>*s[ijk-ii1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+ii1]))
                           + cg2<TF>*(u[ijk+ii1] * (ci0<TF>*s[ijk-ii1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+ii1] + ci3<TF>*s[ijk+ii2]))
                           + cg3<TF>*(u[ijk+ii2] * (ci0<TF>*s[ijk    ] + ci1<TF>*s[ijk+ii1] + ci2<TF>*s[ijk+ii2] + ci3<TF>*s[ijk+ii3])) ) * dxi;

                st[ijk] -= ( cg0<TF>*(v[ijk-jj1] * (ci0<TF>*s[ijk-jj3] + ci1<TF>*s[ijk-jj2] + ci2<TF>*s[ijk-jj1] + ci3<TF>*s[ijk    ]))
                           + cg1<TF>*(v[ijk    ] * (ci0<TF>*s[ijk-jj2] + ci1<TF>*s[ijk-jj1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+jj1]))
                           + cg2<TF>*(v[ijk+jj1] * (ci0<TF>*s[ijk-jj1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+jj1] + ci3<TF>*s[ijk+jj2]))
                           + cg3<TF>*(v[ijk+jj2] * (ci0<TF>*s[ijk    ] + ci1<TF>*s[ijk+jj1] + ci2<TF>*s[ijk+jj2] + ci3<TF>*s[ijk+jj3])) ) * dyi;

                st[ijk] -= ( cg0<TF>*(w[ijk-kk1] * (bi0<TF>*s[ijk-kk2] + bi1<TF>*s[ijk-kk1] + bi2<TF>*s[ijk    ] + bi3<TF>*s[ijk+kk1]))
                           + cg1<TF>*(w[ijk    ] * (ci0<TF>*s[ijk-kk2] + ci1<TF>*s[ijk-kk1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+kk1]))
                           + cg2<TF>*(w[ijk+kk1] * (ci0<TF>*s[ijk-kk1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+kk1] + ci3<TF>*s[ijk+kk2]))
                           + cg3<TF>*(w[ijk+kk2] * (ci0<TF>*s[ijk    ] + ci1<TF>*s[ijk+kk1] + ci2<TF>*s[ijk+kk2] + ci3<TF>*s[ijk+kk3])) ) * dzi4[k];
            }
            else if (k == kend-1)
            {
                st[ijk] -= ( cg0<TF>*(u[ijk-ii1] * (ci0<TF>*s[ijk-ii3] + ci1<TF>*s[ijk-ii2] + ci2<TF>*s[ijk-ii1] + ci3<TF>*s[ijk    ]))
                           + cg1<TF>*(u[ijk    ] * (ci0<TF>*s[ijk-ii2] + ci1<TF>*s[ijk-ii1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+ii1]))
                           + cg2<TF>*(u[ijk+ii1] * (ci0<TF>*s[ijk-ii1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+ii1] + ci3<TF>*s[ijk+ii2]))
                           + cg3<TF>*(u[ijk+ii2] * (ci0<TF>*s[ijk    ] + ci1<TF>*s[ijk+ii1] + ci2<TF>*s[ijk+ii2] + ci3<TF>*s[ijk+ii3])) ) * dxi;

                st[ijk] -= ( cg0<TF>*(v[ijk-jj1] * (ci0<TF>*s[ijk-jj3] + ci1<TF>*s[ijk-jj2] + ci2<TF>*s[ijk-jj1] + ci3<TF>*s[ijk    ]))
                           + cg1<TF>*(v[ijk    ] * (ci0<TF>*s[ijk-jj2] + ci1<TF>*s[ijk-jj1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+jj1]))
                           + cg2<TF>*(v[ijk+jj1] * (ci0<TF>*s[ijk-jj1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+jj1] + ci3<TF>*s[ijk+jj2]))
                           + cg3<TF>*(v[ijk+jj2] * (ci0<TF>*s[ijk    ] + ci1<TF>*s[ijk+jj1] + ci2<TF>*s[ijk+jj2] + ci3<TF>*s[ijk+jj3])) ) * dyi;

                st[ijk] -= ( cg0<TF>*(w[ijk-kk1] * (ci0<TF>*s[ijk-kk3] + ci1<TF>*s[ijk-kk2] + ci2<TF>*s[ijk-kk1] + ci3<TF>*s[ijk    ]))
                           + cg1<TF>*(w[ijk    ] * (ci0<TF>*s[ijk-kk2] + ci1<TF>*s[ijk-kk1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+kk1]))
                           + cg2<TF>*(w[ijk+kk1] * (ci0<TF>*s[ijk-kk1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+kk1] + ci3<TF>*s[ijk+kk2]))
                           + cg3<TF>*(w[ijk+kk2] * (ti0<TF>*s[ijk-kk1] + ti1<TF>*s[ijk    ] + ti2<TF>*s[ijk+kk1] + ti3<TF>*s[ijk+kk2])) ) * dzi4[k];
            }
            else
            {
                st[ijk] -= ( cg0<TF>*(u[ijk-ii1] * (ci0<TF>*s[ijk-ii3] + ci1<TF>*s[ijk-ii2] + ci2<TF>*s[ijk-ii1] + ci3<TF>*s[ijk    ]))
                           + cg1<TF>*(u[ijk    ] * (ci0<TF>*s[ijk-ii2] + ci1<TF>*s[ijk-ii1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+ii1]))
                           + cg2<TF>*(u[ijk+ii1] * (ci0<TF>*s[ijk-ii1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+ii1] + ci3<TF>*s[ijk+ii2]))
                           + cg3<TF>*(u[ijk+ii2] * (ci0<TF>*s[ijk    ] + ci1<TF>*s[ijk+ii1] + ci2<TF>*s[ijk+ii2] + ci3<TF>*s[ijk+ii3])) ) * dxi;

                st[ijk] -= ( cg0<TF>*(v[ijk-jj1] * (ci0<TF>*s[ijk-jj3] + ci1<TF>*s[ijk-jj2] + ci2<TF>*s[ijk-jj1] + ci3<TF>*s[ijk    ]))
                           + cg1<TF>*(v[ijk    ] * (ci0<TF>*s[ijk-jj2] + ci1<TF>*s[ijk-jj1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+jj1]))
                           + cg2<TF>*(v[ijk+jj1] * (ci0<TF>*s[ijk-jj1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+jj1] + ci3<TF>*s[ijk+jj2]))
                           + cg3<TF>*(v[ijk+jj2] * (ci0<TF>*s[ijk    ] + ci1<TF>*s[ijk+jj1] + ci2<TF>*s[ijk+jj2] + ci3<TF>*s[ijk+jj3])) ) * dyi;

                st[ijk] -= ( cg0<TF>*(w[ijk-kk1] * (ci0<TF>*s[ijk-kk3] + ci1<TF>*s[ijk-kk2] + ci2<TF>*s[ijk-kk1] + ci3<TF>*s[ijk    ]))
                           + cg1<TF>*(w[ijk    ] * (ci0<TF>*s[ijk-kk2] + ci1<TF>*s[ijk-kk1] + ci2<TF>*s[ijk    ] + ci3<TF>*s[ijk+kk1]))
                           + cg2<TF>*(w[ijk+kk1] * (ci0<TF>*s[ijk-kk1] + ci1<TF>*s[ijk    ] + ci2<TF>*s[ijk+kk1] + ci3<TF>*s[ijk+kk2]))
                           + cg3<TF>*(w[ijk+kk2] * (ci0<TF>*s[ijk    ] + ci1<TF>*s[ijk+kk1] + ci2<TF>*s[ijk+kk2] + ci3<TF>*s[ijk+kk3])) ) * dzi4[k];
            }
        }
    }

    template<typename TF> __global__
    void calc_cfl_g(TF* const __restrict__ tmp1,
                    const TF* __restrict__ u, const TF* __restrict__ v, const TF* __restrict__ w,
                    const TF* __restrict__ dzi, const TF dxi, const TF dyi,
                    const int jj, const int kk,
                    const int istart, const int jstart, const int kstart,
                    const int iend,   const int jend,   const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        const int ii1 = 1;
        const int ii2 = 2;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        const int ijk = i + j*jj + k*kk;

        if (i < iend && j < jend && k < kend)
            tmp1[ijk] = std::abs(ci0<TF>*u[ijk-ii1] + ci1<TF>*u[ijk] + ci2<TF>*u[ijk+ii1] + ci3<TF>*u[ijk+ii2])*dxi +
                        std::abs(ci0<TF>*v[ijk-jj1] + ci1<TF>*v[ijk] + ci2<TF>*v[ijk+jj1] + ci3<TF>*v[ijk+jj2])*dyi +
                        std::abs(ci0<TF>*w[ijk-kk1] + ci1<TF>*w[ijk] + ci2<TF>*w[ijk+kk1] + ci3<TF>*w[ijk+kk2])*dzi[k];
    }
}
#ifdef USECUDA
template<typename TF>
unsigned long Advec_4<TF>::get_time_limit(unsigned long idt, double dt)
{
    // Calculate cfl and prevent zero divisons.
    double cfl = get_cfl(dt);
    cfl = std::max(cflmin, cfl);
    const unsigned long idtlim = idt * cflmax / cfl;
    return idtlim;
}

template<typename TF>
double Advec_4<TF>::get_cfl(const double dt)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    auto cfl_3d = fields.get_tmp_g();

    calc_cfl_g<TF><<<gridGPU, blockGPU>>>(
        cfl_3d->fld_g,
        fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
        gd.dzi_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    TF cfl = field3d_operators.calc_max_g(cfl_3d->fld_g);
    // TO DO communicate.

    cfl = cfl*dt;

    fields.release_tmp_g(cfl_3d);

    return cfl;
}

template<typename TF>
void Advec_4<TF>::exec(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    dim3 gridGPU2D (gridi, gridj, 1);
    dim3 blockGPU2D(blocki, blockj, 1);

    // Top and bottom boundary:
    advec_u_boundary_g<TF,0><<<gridGPU2D, blockGPU2D>>>(
        fields.mt.at("u")->fld_g, fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g,
        fields.mp.at("w")->fld_g, gd.dzi4_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    advec_u_boundary_g<TF,1><<<gridGPU2D, blockGPU2D>>>(
        fields.mt.at("u")->fld_g, fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g,
        fields.mp.at("w")->fld_g, gd.dzi4_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    // Interior:
    advec_u_g<TF><<<gridGPU, blockGPU>>>(
        fields.mt.at("u")->fld_g, fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g,
        fields.mp.at("w")->fld_g, gd.dzi4_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    // Top and bottom boundary:
    advec_v_boundary_g<TF,0><<<gridGPU2D, blockGPU2D>>>(
        fields.mt.at("v")->fld_g, fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g,
        fields.mp.at("w")->fld_g, gd.dzi4_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    advec_v_boundary_g<TF,1><<<gridGPU2D, blockGPU2D>>>(
        fields.mt.at("v")->fld_g, fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g,
        fields.mp.at("w")->fld_g, gd.dzi4_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    // Interior
    advec_v_g<TF><<<gridGPU, blockGPU>>>(
        fields.mt.at("v")->fld_g, fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g,
        fields.mp.at("w")->fld_g, gd.dzi4_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    // Top and bottom boundary:
    advec_w_boundary_g<TF,0><<<gridGPU2D, blockGPU2D>>>(
        fields.mt.at("w")->fld_g, fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g,
        fields.mp.at("w")->fld_g, gd.dzhi4_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    advec_w_boundary_g<TF,1><<<gridGPU2D, blockGPU2D>>>(
        fields.mt.at("w")->fld_g, fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g,
        fields.mp.at("w")->fld_g, gd.dzhi4_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    // Interior:
    advec_w_g<TF><<<gridGPU, blockGPU>>>(
        fields.mt.at("w")->fld_g, fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g,
        fields.mp.at("w")->fld_g, gd.dzhi4_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    for (auto& it : fields.st)
        advec_s_g<TF><<<gridGPU, blockGPU>>>(
            it.second->fld_g, fields.sp.at(it.first)->fld_g,
            fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
            gd.dzi4_g, gd.dxi, gd.dyi,
            gd.icells, gd.ijcells,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    cudaDeviceSynchronize();

    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);

    for (auto it : fields.st)
        stats.calc_tend(*it.second, tend_name);
}
#endif


#ifdef FLOAT_SINGLE
template class Advec_4<float>;
#else
template class Advec_4<double>;
#endif
