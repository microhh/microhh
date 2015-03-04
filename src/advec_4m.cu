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

#include <cmath>
#include "advec_4m.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "fd.h"
#include "tools.h"
#include "constants.h"

using namespace fd::o2;
using namespace fd::o4;

namespace Advec_4m_g
{
    __global__ 
    void advec_u(double* __restrict__ ut, double* __restrict__ u, 
                 double* __restrict__ v,  double* __restrict__ w,
                 double* __restrict__ dzi4, double dxi, double dyi, 
                 int jj, int kk,
                 int istart, int jstart, int kstart,
                 int iend,   int jend,   int kend)
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
                    - grad4(  interp4(u[ijk-ii3],     u[ijk-ii2],     u[ijk-ii1], u[ijk    ])     * interp2(u[ijk-ii3], u[ijk    ]),
                              interp4(u[ijk-ii2],     u[ijk-ii1],     u[ijk    ], u[ijk+ii1])     * interp2(u[ijk-ii1], u[ijk    ]),
                              interp4(u[ijk-ii1],     u[ijk    ],     u[ijk+ii1], u[ijk+ii2])     * interp2(u[ijk    ], u[ijk+ii1]),
                              interp4(u[ijk    ],     u[ijk+ii1],     u[ijk+ii2], u[ijk+ii3])     * interp2(u[ijk    ], u[ijk+ii3]), dxi)

                    - grad4(  interp4(v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1]) * interp2(u[ijk-jj3], u[ijk    ]),
                              interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp2(u[ijk-jj1], u[ijk    ]),
                              interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp2(u[ijk    ], u[ijk+jj1]),
                              interp4(v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2]) * interp2(u[ijk    ], u[ijk+jj3]), dyi)

                    // boundary condition
                    - grad4x(-interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp2(u[ijk-kk1], u[ijk+kk2]),
                              interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp2(u[ijk-kk1], u[ijk    ]),
                              interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp2(u[ijk    ], u[ijk+kk1]),
                              interp4(w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2]) * interp2(u[ijk    ], u[ijk+kk3])) * dzi4[kstart];
            }
            else if (k == kend-1)
            {
                ut[ijk] +=
                    - grad4(  interp4(u[ijk-ii3],     u[ijk-ii2],     u[ijk-ii1], u[ijk    ])     * interp2(u[ijk-ii3], u[ijk    ]),
                              interp4(u[ijk-ii2],     u[ijk-ii1],     u[ijk    ], u[ijk+ii1])     * interp2(u[ijk-ii1], u[ijk    ]),
                              interp4(u[ijk-ii1],     u[ijk    ],     u[ijk+ii1], u[ijk+ii2])     * interp2(u[ijk    ], u[ijk+ii1]),
                              interp4(u[ijk    ],     u[ijk+ii1],     u[ijk+ii2], u[ijk+ii3])     * interp2(u[ijk    ], u[ijk+ii3]), dxi)

                    - grad4(  interp4(v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1]) * interp2(u[ijk-jj3], u[ijk    ]),
                              interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp2(u[ijk-jj1], u[ijk    ]),
                              interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp2(u[ijk    ], u[ijk+jj1]),
                              interp4(v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2]) * interp2(u[ijk    ], u[ijk+jj3]), dyi)

                    - grad4x( interp4(w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1]) * interp2(u[ijk-kk3], u[ijk    ]),
                              interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp2(u[ijk-kk1], u[ijk    ]),
                              interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp2(u[ijk    ], u[ijk+kk1]),
                             -interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp2(u[ijk-kk2], u[ijk+kk1])) * dzi4[kend-1];
            }
            else
            {
                ut[ijk] +=
                    - grad4(interp4(u[ijk-ii3],     u[ijk-ii2],     u[ijk-ii1], u[ijk    ])     * interp2(u[ijk-ii3], u[ijk    ]),
                            interp4(u[ijk-ii2],     u[ijk-ii1],     u[ijk    ], u[ijk+ii1])     * interp2(u[ijk-ii1], u[ijk    ]),
                            interp4(u[ijk-ii1],     u[ijk    ],     u[ijk+ii1], u[ijk+ii2])     * interp2(u[ijk    ], u[ijk+ii1]),
                            interp4(u[ijk    ],     u[ijk+ii1],     u[ijk+ii2], u[ijk+ii3])     * interp2(u[ijk    ], u[ijk+ii3]), dxi)

                    - grad4(interp4(v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1]) * interp2(u[ijk-jj3], u[ijk    ]),
                            interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp2(u[ijk-jj1], u[ijk    ]),
                            interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp2(u[ijk    ], u[ijk+jj1]),
                            interp4(v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2]) * interp2(u[ijk    ], u[ijk+jj3]), dyi)

                    - grad4x(interp4(w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1]) * interp2(u[ijk-kk3], u[ijk    ]),
                            interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp2(u[ijk-kk1], u[ijk    ]),
                            interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp2(u[ijk    ], u[ijk+kk1]),
                            interp4(w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2]) * interp2(u[ijk    ], u[ijk+kk3])) * dzi4[k];
            }
        }
    }

    __global__ 
    void advec_v(double* __restrict__ vt, double* __restrict__ u, 
                 double* __restrict__ v,  double* __restrict__ w,
                 double* __restrict__ dzi4, double dxi, double dyi, 
                 int jj, int kk,
                 int istart, int jstart, int kstart,
                 int iend,   int jend,   int kend)
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
                    - grad4(  interp4(u[ijk-ii1-jj2], u[ijk-ii1-jj1], u[ijk-ii1], u[ijk-ii1+jj1]) * interp2(v[ijk-ii3], v[ijk    ]),
                              interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp2(v[ijk-ii1], v[ijk    ]),
                              interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp2(v[ijk    ], v[ijk+ii1]),
                              interp4(u[ijk+ii2-jj2], u[ijk+ii2-jj1], u[ijk+ii2], u[ijk+ii2+jj1]) * interp2(v[ijk    ], v[ijk+ii3]), dxi)

                    - grad4(  interp4(v[ijk-jj3],     v[ijk-jj2],     v[ijk-jj1], v[ijk    ])     * interp2(v[ijk-jj3], v[ijk    ]),
                              interp4(v[ijk-jj2],     v[ijk-jj1],     v[ijk    ], v[ijk+jj1])     * interp2(v[ijk-jj1], v[ijk    ]),
                              interp4(v[ijk-jj1],     v[ijk    ],     v[ijk+jj1], v[ijk+jj2])     * interp2(v[ijk    ], v[ijk+jj1]),
                              interp4(v[ijk    ],     v[ijk+jj1],     v[ijk+jj2], v[ijk+jj3])     * interp2(v[ijk    ], v[ijk+jj3]), dyi)

                    - grad4x(-interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp2(v[ijk-kk1], v[ijk+kk2]),
                              interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp2(v[ijk-kk1], v[ijk    ]),
                              interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp2(v[ijk    ], v[ijk+kk1]),
                              interp4(w[ijk-jj2+kk2], w[ijk-jj1+kk2], w[ijk+kk2], w[ijk+jj1+kk2]) * interp2(v[ijk    ], v[ijk+kk3])) * dzi4[kstart];
            }
            else if (k == kend-1)
            {
                vt[ijk] +=
                    - grad4(  interp4(u[ijk-ii1-jj2], u[ijk-ii1-jj1], u[ijk-ii1], u[ijk-ii1+jj1]) * interp2(v[ijk-ii3], v[ijk    ]),
                              interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp2(v[ijk-ii1], v[ijk    ]),
                              interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp2(v[ijk    ], v[ijk+ii1]),
                              interp4(u[ijk+ii2-jj2], u[ijk+ii2-jj1], u[ijk+ii2], u[ijk+ii2+jj1]) * interp2(v[ijk    ], v[ijk+ii3]), dxi)

                    - grad4(  interp4(v[ijk-jj3],     v[ijk-jj2],     v[ijk-jj1], v[ijk    ])     * interp2(v[ijk-jj3], v[ijk    ]),
                              interp4(v[ijk-jj2],     v[ijk-jj1],     v[ijk    ], v[ijk+jj1])     * interp2(v[ijk-jj1], v[ijk    ]),
                              interp4(v[ijk-jj1],     v[ijk    ],     v[ijk+jj1], v[ijk+jj2])     * interp2(v[ijk    ], v[ijk+jj1]),
                              interp4(v[ijk    ],     v[ijk+jj1],     v[ijk+jj2], v[ijk+jj3])     * interp2(v[ijk    ], v[ijk+jj3]), dyi)

                    - grad4x( interp4(w[ijk-jj2-kk1], w[ijk-jj1-kk1], w[ijk-kk1], w[ijk+jj1-kk1]) * interp2(v[ijk-kk3], v[ijk    ]),
                              interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp2(v[ijk-kk1], v[ijk    ]),
                              interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp2(v[ijk    ], v[ijk+kk1]),
                             -interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp2(v[ijk-kk2], v[ijk+kk1])) * dzi4[kend-1];
            }
            else
            {
                vt[ijk] +=
                    - grad4( interp4(u[ijk-ii1-jj2], u[ijk-ii1-jj1], u[ijk-ii1], u[ijk-ii1+jj1]) * interp2(v[ijk-ii3], v[ijk    ]),
                             interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp2(v[ijk-ii1], v[ijk    ]),
                             interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp2(v[ijk    ], v[ijk+ii1]),
                             interp4(u[ijk+ii2-jj2], u[ijk+ii2-jj1], u[ijk+ii2], u[ijk+ii2+jj1]) * interp2(v[ijk    ], v[ijk+ii3]), dxi)

                    - grad4( interp4(v[ijk-jj3],     v[ijk-jj2],     v[ijk-jj1], v[ijk    ])     * interp2(v[ijk-jj3], v[ijk    ]),
                             interp4(v[ijk-jj2],     v[ijk-jj1],     v[ijk    ], v[ijk+jj1])     * interp2(v[ijk-jj1], v[ijk    ]),
                             interp4(v[ijk-jj1],     v[ijk    ],     v[ijk+jj1], v[ijk+jj2])     * interp2(v[ijk    ], v[ijk+jj1]),
                             interp4(v[ijk    ],     v[ijk+jj1],     v[ijk+jj2], v[ijk+jj3])     * interp2(v[ijk    ], v[ijk+jj3]), dyi)

                    - grad4x(interp4(w[ijk-jj2-kk1], w[ijk-jj1-kk1], w[ijk-kk1], w[ijk+jj1-kk1]) * interp2(v[ijk-kk3], v[ijk    ]),
                             interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp2(v[ijk-kk1], v[ijk    ]),
                             interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp2(v[ijk    ], v[ijk+kk1]),
                             interp4(w[ijk-jj2+kk2], w[ijk-jj1+kk2], w[ijk+kk2], w[ijk+jj1+kk2]) * interp2(v[ijk    ], v[ijk+kk3])) * dzi4[k];
            }
        }
    }

    __global__ 
    void advec_w(double* __restrict__ wt, double* __restrict__ u, 
                 double* __restrict__ v,  double* __restrict__ w,
                 double* __restrict__ dzhi4, double dxi, double dyi, 
                 int jj, int kk, 
                 int istart, int jstart, int kstart,
                 int iend,   int jend,   int kend)
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
                - grad4( interp4(u[ijk-ii1-kk2], u[ijk-ii1-kk1], u[ijk-ii1], u[ijk-ii1+kk1]) * interp2(w[ijk-ii3], w[ijk    ]),
                         interp4(u[ijk    -kk2], u[ijk    -kk1], u[ijk    ], u[ijk    +kk1]) * interp2(w[ijk-ii1], w[ijk    ]),
                         interp4(u[ijk+ii1-kk2], u[ijk+ii1-kk1], u[ijk+ii1], u[ijk+ii1+kk1]) * interp2(w[ijk    ], w[ijk+ii1]),
                         interp4(u[ijk+ii2-kk2], u[ijk+ii2-kk1], u[ijk+ii2], u[ijk+ii2+kk1]) * interp2(w[ijk    ], w[ijk+ii3]), dxi)

                - grad4( interp4(v[ijk-jj1-kk2], v[ijk-jj1-kk1], v[ijk-jj1], v[ijk-jj1+kk1]) * interp2(w[ijk-jj3], w[ijk    ]),
                         interp4(v[ijk    -kk2], v[ijk    -kk1], v[ijk    ], v[ijk    +kk1]) * interp2(w[ijk-jj1], w[ijk    ]),
                         interp4(v[ijk+jj1-kk2], v[ijk+jj1-kk1], v[ijk+jj1], v[ijk+jj1+kk1]) * interp2(w[ijk    ], w[ijk+jj1]),
                         interp4(v[ijk+jj2-kk2], v[ijk+jj2-kk1], v[ijk+jj2], v[ijk+jj2+kk1]) * interp2(w[ijk    ], w[ijk+jj3]), dyi)

                - grad4x(interp4(w[ijk-kk3],     w[ijk-kk2],     w[ijk-kk1], w[ijk    ])     * interp2(w[ijk-kk3], w[ijk    ]),
                         interp4(w[ijk-kk2],     w[ijk-kk1],     w[ijk    ], w[ijk+kk1])     * interp2(w[ijk-kk1], w[ijk    ]),
                         interp4(w[ijk-kk1],     w[ijk    ],     w[ijk+kk1], w[ijk+kk2])     * interp2(w[ijk    ], w[ijk+kk1]),
                         interp4(w[ijk    ],     w[ijk+kk1],     w[ijk+kk2], w[ijk+kk3])     * interp2(w[ijk    ], w[ijk+kk3])) * dzhi4[k];
        }
    }

    __global__ 
    void advec_s(double* __restrict__ st, double* __restrict__ s, 
                 double* __restrict__ u,  double* __restrict__ v, double* __restrict__ w,
                 double* __restrict__ dzi4, double dxi, double dyi, 
                 int jj, int kk,
                 int istart, int jstart, int kstart,
                 int iend,   int jend,   int kend)
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
                             u[ijk+ii2] * interp2(s[ijk    ], s[ijk+ii3]), dxi)

                    - grad4( v[ijk-jj1] * interp2(s[ijk-jj3], s[ijk    ]),
                             v[ijk    ] * interp2(s[ijk-jj1], s[ijk    ]),
                             v[ijk+jj1] * interp2(s[ijk    ], s[ijk+jj1]),
                             v[ijk+jj2] * interp2(s[ijk    ], s[ijk+jj3]), dyi)

                    -grad4x(-w[ijk+kk1] * interp2(s[ijk-kk1], s[ijk+kk2]),
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
                             u[ijk+ii2] * interp2(s[ijk    ], s[ijk+ii3]), dxi)

                    - grad4( v[ijk-jj1] * interp2(s[ijk-jj3], s[ijk    ]),
                             v[ijk    ] * interp2(s[ijk-jj1], s[ijk    ]),
                             v[ijk+jj1] * interp2(s[ijk    ], s[ijk+jj1]),
                             v[ijk+jj2] * interp2(s[ijk    ], s[ijk+jj3]), dyi)

                   - grad4x( w[ijk-kk1] * interp2(s[ijk-kk3], s[ijk    ]),
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
                            u[ijk+ii2] * interp2(s[ijk    ], s[ijk+ii3]), dxi)

                    - grad4(v[ijk-jj1] * interp2(s[ijk-jj3], s[ijk    ]),
                            v[ijk    ] * interp2(s[ijk-jj1], s[ijk    ]),
                            v[ijk+jj1] * interp2(s[ijk    ], s[ijk+jj1]),
                            v[ijk+jj2] * interp2(s[ijk    ], s[ijk+jj3]), dyi)

                   - grad4x(w[ijk-kk1] * interp2(s[ijk-kk3], s[ijk    ]),
                            w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]),
                            w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1]),
                            w[ijk+kk2] * interp2(s[ijk    ], s[ijk+kk3])) * dzi4[k];
            }
        }
    }

    __global__ 
    void calc_cfl(double* const __restrict__ tmp1,
                  const double* const __restrict__ u,   const double* const __restrict__ v, const double* const __restrict__ w, 
                  const double* const __restrict__ dzi, const double dxi, const double dyi,
                  const int jj, const int kk,
                  const int istart, const int jstart, const int kstart,
                  const int iend,   const int jend,   const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y; 
        const int k = blockIdx.z; 

        const int ii1 = 1;
        const int ii2 = 2;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        const int ijk = i + j*jj + k*kk;

        if (i < iend && j < jend && k < kend)
            tmp1[ijk] = std::abs(ci0*u[ijk-ii1] + ci1*u[ijk] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2])*dxi + 
                        std::abs(ci0*v[ijk-jj1] + ci1*v[ijk] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2])*dyi + 
                        std::abs(ci0*w[ijk-kk1] + ci1*w[ijk] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])*dzi[k];
    }
}

#ifdef USECUDA
unsigned long Advec_4m::get_time_limit(unsigned long idt, double dt)
{
    // Calculate cfl and prevent zero divisons.
    double cfl = get_cfl(dt);
    cfl = std::max(cflmin, cfl);
    const unsigned long idtlim = idt * cflmax / cfl;
    return idtlim;
}

double Advec_4m::get_cfl(const double dt)
{
    const int blocki = grid->iThreadBlock;
    const int blockj = grid->jThreadBlock;
    const int gridi  = grid->icells/blocki + (grid->icells%blocki > 0);
    const int gridj  = grid->jcells/blockj + (grid->jcells%blockj > 0);

    dim3 gridGPU (gridi, gridj, grid->kcells);
    dim3 blockGPU(blocki, blockj, 1);

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    const int offs = grid->memoffset;

    Advec_4m_g::calc_cfl<<<gridGPU, blockGPU>>>(
        &fields->atmp["tmp1"]->data_g[offs],
        &fields->u->data_g[offs], &fields->v->data_g[offs], &fields->w->data_g[offs],
        grid->dzi_g, dxi, dyi,
        grid->icellsp, grid->ijcellsp,
        grid->istart,  grid->jstart, grid->kstart,
        grid->iend,    grid->jend,   grid->kend);
    cudaCheckError(); 

    double cfl = grid->getMax_g(&fields->atmp["tmp1"]->data_g[offs], fields->atmp["tmp2"]->data_g); 
    grid->getMax(&cfl); 
    cfl = cfl*dt;

    return cfl;
}

void Advec_4m::exec()
{
    const int blocki = grid->iThreadBlock;
    const int blockj = grid->jThreadBlock;
    const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
    const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, grid->kmax);
    dim3 blockGPU(blocki, blockj, 1);

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    const int offs = grid->memoffset;

    Advec_4m_g::advec_u<<<gridGPU, blockGPU>>>(
        &fields->ut->data_g[offs], &fields->u->data_g[offs], &fields->v->data_g[offs], 
        &fields->w->data_g[offs], grid->dzi4_g, dxi, dyi,
        grid->icellsp, grid->ijcellsp,
        grid->istart,  grid->jstart, grid->kstart,
        grid->iend,    grid->jend,   grid->kend);
    cudaCheckError(); 

    Advec_4m_g::advec_v<<<gridGPU, blockGPU>>>(
        &fields->vt->data_g[offs], &fields->u->data_g[offs], &fields->v->data_g[offs], 
        &fields->w->data_g[offs], grid->dzi4_g, dxi, dyi,
        grid->icellsp, grid->ijcellsp,
        grid->istart,  grid->jstart, grid->kstart,
        grid->iend,    grid->jend,   grid->kend);
    cudaCheckError(); 

    Advec_4m_g::advec_w<<<gridGPU, blockGPU>>>(
        &fields->wt->data_g[offs], &fields->u->data_g[offs], &fields->v->data_g[offs], 
        &fields->w->data_g[offs], grid->dzhi4_g, dxi, dyi,
        grid->icellsp, grid->ijcellsp,
        grid->istart,  grid->jstart, grid->kstart,
        grid->iend,    grid->jend,   grid->kend);
    cudaCheckError(); 

    for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
        Advec_4m_g::advec_s<<<gridGPU, blockGPU>>>(
            &it->second->data_g[offs], &fields->sp[it->first]->data_g[offs], 
            &fields->u->data_g[offs], &fields->v->data_g[offs], &fields->w->data_g[offs], 
            grid->dzi4_g, dxi, dyi,
            grid->icellsp, grid->ijcellsp,
            grid->istart,  grid->jstart, grid->kstart,
            grid->iend,    grid->jend,   grid->kend);
    cudaCheckError(); 
}
#endif
