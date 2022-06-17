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

#include "advec_4m.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "tools.h"
#include "constants.h"
#include "finite_difference.h"
#include "field3d_operators.h"

using namespace Finite_difference::O2;
using namespace Finite_difference::O4;

namespace
{
    template <typename TF> __global__
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

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;

            if (k == kstart)
            {
                ut[ijk] +=
                    - grad4(  interp4c(u[ijk-ii3],     u[ijk-ii2],     u[ijk-ii1], u[ijk    ])     * interp2(u[ijk-ii3], u[ijk    ]),
                              interp4c(u[ijk-ii2],     u[ijk-ii1],     u[ijk    ], u[ijk+ii1])     * interp2(u[ijk-ii1], u[ijk    ]),
                              interp4c(u[ijk-ii1],     u[ijk    ],     u[ijk+ii1], u[ijk+ii2])     * interp2(u[ijk    ], u[ijk+ii1]),
                              interp4c(u[ijk    ],     u[ijk+ii1],     u[ijk+ii2], u[ijk+ii3])     * interp2(u[ijk    ], u[ijk+ii3])) * dxi

                    - grad4(  interp4c(v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1]) * interp2(u[ijk-jj3], u[ijk    ]),
                              interp4c(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp2(u[ijk-jj1], u[ijk    ]),
                              interp4c(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp2(u[ijk    ], u[ijk+jj1]),
                              interp4c(v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2]) * interp2(u[ijk    ], u[ijk+jj3])) * dyi

                    // boundary condition
                    - grad4(-interp4c(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp2(u[ijk-kk1], u[ijk+kk2]),
                             interp4c(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp2(u[ijk-kk1], u[ijk    ]),
                             interp4c(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp2(u[ijk    ], u[ijk+kk1]),
                             interp4c(w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2]) * interp2(u[ijk    ], u[ijk+kk3])) * dzi4[kstart];
            }
            else if (k == kend-1)
            {
                ut[ijk] +=
                    - grad4(  interp4c(u[ijk-ii3],     u[ijk-ii2],     u[ijk-ii1], u[ijk    ])     * interp2(u[ijk-ii3], u[ijk    ]),
                              interp4c(u[ijk-ii2],     u[ijk-ii1],     u[ijk    ], u[ijk+ii1])     * interp2(u[ijk-ii1], u[ijk    ]),
                              interp4c(u[ijk-ii1],     u[ijk    ],     u[ijk+ii1], u[ijk+ii2])     * interp2(u[ijk    ], u[ijk+ii1]),
                              interp4c(u[ijk    ],     u[ijk+ii1],     u[ijk+ii2], u[ijk+ii3])     * interp2(u[ijk    ], u[ijk+ii3])) * dxi

                    - grad4(  interp4c(v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1]) * interp2(u[ijk-jj3], u[ijk    ]),
                              interp4c(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp2(u[ijk-jj1], u[ijk    ]),
                              interp4c(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp2(u[ijk    ], u[ijk+jj1]),
                              interp4c(v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2]) * interp2(u[ijk    ], u[ijk+jj3])) * dyi

                    - grad4(  interp4c(w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1]) * interp2(u[ijk-kk3], u[ijk    ]),
                              interp4c(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp2(u[ijk-kk1], u[ijk    ]),
                              interp4c(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp2(u[ijk    ], u[ijk+kk1]),
                             -interp4c(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp2(u[ijk-kk2], u[ijk+kk1])) * dzi4[kend-1];
            }
            else
            {
                ut[ijk] +=
                    - grad4(interp4c(u[ijk-ii3],     u[ijk-ii2],     u[ijk-ii1], u[ijk    ])     * interp2(u[ijk-ii3], u[ijk    ]),
                            interp4c(u[ijk-ii2],     u[ijk-ii1],     u[ijk    ], u[ijk+ii1])     * interp2(u[ijk-ii1], u[ijk    ]),
                            interp4c(u[ijk-ii1],     u[ijk    ],     u[ijk+ii1], u[ijk+ii2])     * interp2(u[ijk    ], u[ijk+ii1]),
                            interp4c(u[ijk    ],     u[ijk+ii1],     u[ijk+ii2], u[ijk+ii3])     * interp2(u[ijk    ], u[ijk+ii3])) * dxi

                    - grad4(interp4c(v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1]) * interp2(u[ijk-jj3], u[ijk    ]),
                            interp4c(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp2(u[ijk-jj1], u[ijk    ]),
                            interp4c(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp2(u[ijk    ], u[ijk+jj1]),
                            interp4c(v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2]) * interp2(u[ijk    ], u[ijk+jj3])) * dyi

                    - grad4(interp4c(w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1]) * interp2(u[ijk-kk3], u[ijk    ]),
                            interp4c(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp2(u[ijk-kk1], u[ijk    ]),
                            interp4c(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp2(u[ijk    ], u[ijk+kk1]),
                            interp4c(w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2]) * interp2(u[ijk    ], u[ijk+kk3])) * dzi4[k];
            }
        }
    }

    template <typename TF> __global__
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

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;

            if (k == kstart)
            {
                vt[ijk] +=
                    - grad4(  interp4c(u[ijk-ii1-jj2], u[ijk-ii1-jj1], u[ijk-ii1], u[ijk-ii1+jj1]) * interp2(v[ijk-ii3], v[ijk    ]),
                              interp4c(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp2(v[ijk-ii1], v[ijk    ]),
                              interp4c(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp2(v[ijk    ], v[ijk+ii1]),
                              interp4c(u[ijk+ii2-jj2], u[ijk+ii2-jj1], u[ijk+ii2], u[ijk+ii2+jj1]) * interp2(v[ijk    ], v[ijk+ii3])) * dxi

                    - grad4(  interp4c(v[ijk-jj3],     v[ijk-jj2],     v[ijk-jj1], v[ijk    ])     * interp2(v[ijk-jj3], v[ijk    ]),
                              interp4c(v[ijk-jj2],     v[ijk-jj1],     v[ijk    ], v[ijk+jj1])     * interp2(v[ijk-jj1], v[ijk    ]),
                              interp4c(v[ijk-jj1],     v[ijk    ],     v[ijk+jj1], v[ijk+jj2])     * interp2(v[ijk    ], v[ijk+jj1]),
                              interp4c(v[ijk    ],     v[ijk+jj1],     v[ijk+jj2], v[ijk+jj3])     * interp2(v[ijk    ], v[ijk+jj3])) * dyi

                    - grad4( -interp4c(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp2(v[ijk-kk1], v[ijk+kk2]),
                              interp4c(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp2(v[ijk-kk1], v[ijk    ]),
                              interp4c(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp2(v[ijk    ], v[ijk+kk1]),
                              interp4c(w[ijk-jj2+kk2], w[ijk-jj1+kk2], w[ijk+kk2], w[ijk+jj1+kk2]) * interp2(v[ijk    ], v[ijk+kk3])) * dzi4[kstart];
            }
            else if (k == kend-1)
            {
                vt[ijk] +=
                    - grad4(  interp4c(u[ijk-ii1-jj2], u[ijk-ii1-jj1], u[ijk-ii1], u[ijk-ii1+jj1]) * interp2(v[ijk-ii3], v[ijk    ]),
                              interp4c(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp2(v[ijk-ii1], v[ijk    ]),
                              interp4c(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp2(v[ijk    ], v[ijk+ii1]),
                              interp4c(u[ijk+ii2-jj2], u[ijk+ii2-jj1], u[ijk+ii2], u[ijk+ii2+jj1]) * interp2(v[ijk    ], v[ijk+ii3])) * dxi

                    - grad4(  interp4c(v[ijk-jj3],     v[ijk-jj2],     v[ijk-jj1], v[ijk    ])     * interp2(v[ijk-jj3], v[ijk    ]),
                              interp4c(v[ijk-jj2],     v[ijk-jj1],     v[ijk    ], v[ijk+jj1])     * interp2(v[ijk-jj1], v[ijk    ]),
                              interp4c(v[ijk-jj1],     v[ijk    ],     v[ijk+jj1], v[ijk+jj2])     * interp2(v[ijk    ], v[ijk+jj1]),
                              interp4c(v[ijk    ],     v[ijk+jj1],     v[ijk+jj2], v[ijk+jj3])     * interp2(v[ijk    ], v[ijk+jj3])) * dyi

                    - grad4(  interp4c(w[ijk-jj2-kk1], w[ijk-jj1-kk1], w[ijk-kk1], w[ijk+jj1-kk1]) * interp2(v[ijk-kk3], v[ijk    ]),
                              interp4c(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp2(v[ijk-kk1], v[ijk    ]),
                              interp4c(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp2(v[ijk    ], v[ijk+kk1]),
                             -interp4c(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp2(v[ijk-kk2], v[ijk+kk1])) * dzi4[kend-1];
            }
            else
            {
                vt[ijk] +=
                    - grad4( interp4c(u[ijk-ii1-jj2], u[ijk-ii1-jj1], u[ijk-ii1], u[ijk-ii1+jj1]) * interp2(v[ijk-ii3], v[ijk    ]),
                             interp4c(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp2(v[ijk-ii1], v[ijk    ]),
                             interp4c(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp2(v[ijk    ], v[ijk+ii1]),
                             interp4c(u[ijk+ii2-jj2], u[ijk+ii2-jj1], u[ijk+ii2], u[ijk+ii2+jj1]) * interp2(v[ijk    ], v[ijk+ii3])) * dxi

                    - grad4( interp4c(v[ijk-jj3],     v[ijk-jj2],     v[ijk-jj1], v[ijk    ])     * interp2(v[ijk-jj3], v[ijk    ]),
                             interp4c(v[ijk-jj2],     v[ijk-jj1],     v[ijk    ], v[ijk+jj1])     * interp2(v[ijk-jj1], v[ijk    ]),
                             interp4c(v[ijk-jj1],     v[ijk    ],     v[ijk+jj1], v[ijk+jj2])     * interp2(v[ijk    ], v[ijk+jj1]),
                             interp4c(v[ijk    ],     v[ijk+jj1],     v[ijk+jj2], v[ijk+jj3])     * interp2(v[ijk    ], v[ijk+jj3])) * dyi

                    - grad4( interp4c(w[ijk-jj2-kk1], w[ijk-jj1-kk1], w[ijk-kk1], w[ijk+jj1-kk1]) * interp2(v[ijk-kk3], v[ijk    ]),
                             interp4c(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp2(v[ijk-kk1], v[ijk    ]),
                             interp4c(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp2(v[ijk    ], v[ijk+kk1]),
                             interp4c(w[ijk-jj2+kk2], w[ijk-jj1+kk2], w[ijk+kk2], w[ijk+jj1+kk2]) * interp2(v[ijk    ], v[ijk+kk3])) * dzi4[k];
            }
        }
    }

    template <typename TF> __global__
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

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;

            wt[ijk] +=
                - grad4( interp4c(u[ijk-ii1-kk2], u[ijk-ii1-kk1], u[ijk-ii1], u[ijk-ii1+kk1]) * interp2(w[ijk-ii3], w[ijk    ]),
                         interp4c(u[ijk    -kk2], u[ijk    -kk1], u[ijk    ], u[ijk    +kk1]) * interp2(w[ijk-ii1], w[ijk    ]),
                         interp4c(u[ijk+ii1-kk2], u[ijk+ii1-kk1], u[ijk+ii1], u[ijk+ii1+kk1]) * interp2(w[ijk    ], w[ijk+ii1]),
                         interp4c(u[ijk+ii2-kk2], u[ijk+ii2-kk1], u[ijk+ii2], u[ijk+ii2+kk1]) * interp2(w[ijk    ], w[ijk+ii3])) * dxi

                - grad4( interp4c(v[ijk-jj1-kk2], v[ijk-jj1-kk1], v[ijk-jj1], v[ijk-jj1+kk1]) * interp2(w[ijk-jj3], w[ijk    ]),
                         interp4c(v[ijk    -kk2], v[ijk    -kk1], v[ijk    ], v[ijk    +kk1]) * interp2(w[ijk-jj1], w[ijk    ]),
                         interp4c(v[ijk+jj1-kk2], v[ijk+jj1-kk1], v[ijk+jj1], v[ijk+jj1+kk1]) * interp2(w[ijk    ], w[ijk+jj1]),
                         interp4c(v[ijk+jj2-kk2], v[ijk+jj2-kk1], v[ijk+jj2], v[ijk+jj2+kk1]) * interp2(w[ijk    ], w[ijk+jj3])) * dyi

                - grad4( interp4c(w[ijk-kk3],     w[ijk-kk2],     w[ijk-kk1], w[ijk    ])     * interp2(w[ijk-kk3], w[ijk    ]),
                         interp4c(w[ijk-kk2],     w[ijk-kk1],     w[ijk    ], w[ijk+kk1])     * interp2(w[ijk-kk1], w[ijk    ]),
                         interp4c(w[ijk-kk1],     w[ijk    ],     w[ijk+kk1], w[ijk+kk2])     * interp2(w[ijk    ], w[ijk+kk1]),
                         interp4c(w[ijk    ],     w[ijk+kk1],     w[ijk+kk2], w[ijk+kk3])     * interp2(w[ijk    ], w[ijk+kk3])) * dzhi4[k];
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
                st[ijk] +=
                    - grad4( u[ijk-ii1] * interp2(s[ijk-ii3], s[ijk    ]),
                             u[ijk    ] * interp2(s[ijk-ii1], s[ijk    ]),
                             u[ijk+ii1] * interp2(s[ijk    ], s[ijk+ii1]),
                             u[ijk+ii2] * interp2(s[ijk    ], s[ijk+ii3])) * dxi

                    - grad4( v[ijk-jj1] * interp2(s[ijk-jj3], s[ijk    ]),
                             v[ijk    ] * interp2(s[ijk-jj1], s[ijk    ]),
                             v[ijk+jj1] * interp2(s[ijk    ], s[ijk+jj1]),
                             v[ijk+jj2] * interp2(s[ijk    ], s[ijk+jj3])) * dyi

                     -grad4(-w[ijk+kk1] * interp2(s[ijk-kk1], s[ijk+kk2]),
                             w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]),
                             w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1]),
                             w[ijk+kk2] * interp2(s[ijk    ], s[ijk+kk3])) * dzi4[kstart];
            }
            else if (k == kend-1)
            {
                st[ijk] +=
                    - grad4( u[ijk-ii1] * interp2(s[ijk-ii3], s[ijk    ]),
                             u[ijk    ] * interp2(s[ijk-ii1], s[ijk    ]),
                             u[ijk+ii1] * interp2(s[ijk    ], s[ijk+ii1]),
                             u[ijk+ii2] * interp2(s[ijk    ], s[ijk+ii3])) * dxi

                    - grad4( v[ijk-jj1] * interp2(s[ijk-jj3], s[ijk    ]),
                             v[ijk    ] * interp2(s[ijk-jj1], s[ijk    ]),
                             v[ijk+jj1] * interp2(s[ijk    ], s[ijk+jj1]),
                             v[ijk+jj2] * interp2(s[ijk    ], s[ijk+jj3])) * dyi

                    - grad4( w[ijk-kk1] * interp2(s[ijk-kk3], s[ijk    ]),
                             w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]),
                             w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1]),
                            -w[ijk    ] * interp2(s[ijk-kk2], s[ijk+kk1])) * dzi4[kend-1];
            }
            else
            {
                st[ijk] +=
                    - grad4(u[ijk-ii1] * interp2(s[ijk-ii3], s[ijk    ]),
                            u[ijk    ] * interp2(s[ijk-ii1], s[ijk    ]),
                            u[ijk+ii1] * interp2(s[ijk    ], s[ijk+ii1]),
                            u[ijk+ii2] * interp2(s[ijk    ], s[ijk+ii3])) * dxi

                    - grad4(v[ijk-jj1] * interp2(s[ijk-jj3], s[ijk    ]),
                            v[ijk    ] * interp2(s[ijk-jj1], s[ijk    ]),
                            v[ijk+jj1] * interp2(s[ijk    ], s[ijk+jj1]),
                            v[ijk+jj2] * interp2(s[ijk    ], s[ijk+jj3])) * dyi

                    - grad4(w[ijk-kk1] * interp2(s[ijk-kk3], s[ijk    ]),
                            w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]),
                            w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1]),
                            w[ijk+kk2] * interp2(s[ijk    ], s[ijk+kk3])) * dzi4[k];
            }
        }
    }

    template <typename TF> __global__
    void calc_cfl_g(TF* const __restrict__ tmp1,
                    const TF* __restrict__ u,   const TF* __restrict__ v, const TF* __restrict__ w,
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
unsigned long Advec_4m<TF>::get_time_limit(unsigned long idt, double dt)
{
    // Calculate cfl and prevent zero divisons.
    double cfl = get_cfl(dt);
    cfl = std::max(cflmin, cfl);
    const unsigned long idtlim = idt * cflmax / cfl;
    return idtlim;
}

template<typename TF>
double Advec_4m<TF>::get_cfl(const double dt)
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
void Advec_4m<TF>::exec(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    advec_u_g<TF><<<gridGPU, blockGPU>>>(
        fields.mt.at("u")->fld_g, fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g,
        fields.mp.at("w")->fld_g, gd.dzi4_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    advec_v_g<TF><<<gridGPU, blockGPU>>>(
        fields.mt.at("v")->fld_g, fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g,
        fields.mp.at("w")->fld_g, gd.dzi4_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

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

template class Advec_4m<double>;
template class Advec_4m<float>;
