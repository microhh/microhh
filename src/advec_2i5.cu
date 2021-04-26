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

#include "advec_2i5.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "tools.h"
#include "constants.h"
#include "finite_difference.h"
#include "field3d_operators.h"


namespace
{
    using namespace Finite_difference::O2;
    using namespace Finite_difference::O4;
    using namespace Finite_difference::O6;


    template<typename TF> __global__
    void advec_u_g(TF* __restrict__ ut, const TF* __restrict__ u,
                   const TF* __restrict__ v,  const TF* __restrict__ w,
                   const TF* __restrict__ rhoref, const TF* __restrict__ rhorefh,
                   const TF* __restrict__ dzi, const TF dxi, const TF dyi,
                   const int jj, int kk,
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

            ut[ijk] +=
                // u*du/dx
                - ( interp2(u[ijk        ], u[ijk+ii1]) * interp6_ws(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3])
                  - interp2(u[ijk-ii1    ], u[ijk    ]) * interp6_ws(u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]) ) * dxi

                + ( fabs(interp2(u[ijk        ], u[ijk+ii1])) * interp5_ws(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3])
                  - fabs(interp2(u[ijk-ii1    ], u[ijk    ])) * interp5_ws(u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]) ) * dxi

                // v*du/dy
                - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp6_ws(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2], u[ijk+jj3])
                  - interp2(v[ijk-ii1    ], v[ijk    ]) * interp6_ws(u[ijk-jj3], u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2]) ) * dyi

                + ( fabs(interp2(v[ijk-ii1+jj1], v[ijk+jj1])) * interp5_ws(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2], u[ijk+jj3])
                  - fabs(interp2(v[ijk-ii1    ], v[ijk    ])) * interp5_ws(u[ijk-jj3], u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2]) ) * dyi;

            if (k == kstart)
            {
                // w*du/dz -> second order interpolation for fluxtop, fluxbot = 0. as w=0
                ut[ijk] +=
                    - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp2(u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }
            else if (k == kstart+1)
            {
                ut[ijk] +=
                    // w*du/dz -> second order interpolation for fluxbot, fourth order for fluxtop
                    - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp4_ws(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2])
                      - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp2(   u[ijk-kk1], u[ijk    ]) ) / rhoref[k] * dzi[k]

                    + ( rhorefh[k+1] * fabs(interp2(w[ijk-ii1+kk1], w[ijk+kk1])) * interp3_ws(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]) ) / rhoref[k] * dzi[k];
            }
            else if (k == kstart+2)
            {
                ut[ijk] +=
                    // w*du/dz -> fourth order interpolation for fluxbot
                    - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp6_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2], u[ijk+kk3])
                      - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp4_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k]

                    + ( rhorefh[k+1] * fabs(interp2(w[ijk-ii1+kk1], w[ijk+kk1])) * interp5_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2], u[ijk+kk3])
                      - rhorefh[k  ] * fabs(interp2(w[ijk-ii1    ], w[ijk    ])) * interp3_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }
            else if (k == kend-3)
            {
                ut[ijk] +=
                    // w*du/dz -> fourth order interpolation for fluxtop
                    - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp4_ws(u[ijk-kk1   ], u[ijk    ], u[ijk+kk1], u[ijk+kk2])
                      - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp6_ws(u[ijk-kk3], u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]) ) / rhoref[k] * dzi[k]

                    + ( rhorefh[k+1] * fabs(interp2(w[ijk-ii1+kk1], w[ijk+kk1])) * interp3_ws(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2])
                      - rhorefh[k  ] * fabs(interp2(w[ijk-ii1    ], w[ijk    ])) * interp5_ws(u[ijk-kk3], u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]) ) / rhoref[k] * dzi[k];
            }
            else if (k == kend-2)
            {
                ut[ijk] +=
                    // w*du/dz -> second order interpolation for fluxtop, fourth order for fluxbot
                    - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp2(u[ijk    ], u[ijk+kk1])
                      - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp4_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k]

                    - ( rhorefh[k  ] * fabs(interp2(w[ijk-ii1    ], w[ijk    ])) * interp3_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }
            else if (k == kend-1)
            {
                ut[ijk] +=
                    // w*du/dz -> second order interpolation for fluxbot, fluxtop=0 as w=0
                    - ( -rhorefh[k] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp2(u[ijk-kk1], u[ijk    ]) ) / rhoref[k] * dzi[k];
            }
            else
            {
                ut[ijk] +=
                    // w*du/dz
                    - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp6_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2], u[ijk+kk3])
                      - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp6_ws(u[ijk-kk3], u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]) ) / rhoref[k] * dzi[k]

                    + ( rhorefh[k+1] * fabs(interp2(w[ijk-ii1+kk1], w[ijk+kk1])) * interp5_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2], u[ijk+kk3])
                      - rhorefh[k  ] * fabs(interp2(w[ijk-ii1    ], w[ijk    ])) * interp5_ws(u[ijk-kk3], u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]) ) / rhoref[k] * dzi[k];
            }
        }
    }


    template<typename TF> __global__
    void advec_v_g(TF* __restrict__ vt, const TF* __restrict__ u,
                   const TF* __restrict__ v,  const TF* __restrict__ w,
                   const TF* __restrict__ rhoref, const TF* __restrict__ rhorefh,
                   const TF* __restrict__ dzi, const TF dxi, const TF dyi,
                   const int jj, int kk,
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

            vt[ijk] +=
                // u*dv/dx
                - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp6_ws(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2], v[ijk+ii3])
                  - interp2(u[ijk    -jj1], u[ijk    ]) * interp6_ws(v[ijk-ii3], v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2]) ) * dxi

                + ( fabs(interp2(u[ijk+ii1-jj1], u[ijk+ii1])) * interp5_ws(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2], v[ijk+ii3])
                  - fabs(interp2(u[ijk    -jj1], u[ijk    ])) * interp5_ws(v[ijk-ii3], v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2]) ) * dxi

                // v*dv/dy
                - ( interp2(v[ijk        ], v[ijk+jj1]) * interp6_ws(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3])
                  - interp2(v[ijk-jj1    ], v[ijk    ]) * interp6_ws(v[ijk-jj3], v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]) ) * dyi

                + ( fabs(interp2(v[ijk        ], v[ijk+jj1])) * interp5_ws(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3])
                  - fabs(interp2(v[ijk-jj1    ], v[ijk    ])) * interp5_ws(v[ijk-jj3], v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]) ) * dyi;

            if (k == kstart)
            {
                vt[ijk] +=
                    // w*dv/dz -> second order interpolation for fluxtop, fluxbot=0 as w=0
                    - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp2(v[ijk    ], v[ijk+kk1])) / rhoref[k] * dzi[k];
            }
            else if (k == kstart+1)
            {
                vt[ijk] +=
                    // w*dv/dz -> second order interpolation for fluxbot, fourth order for fluxtop
                    - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp4_ws(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])
                      - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp2(v[ijk-kk1], v[ijk    ]) ) / rhoref[k] * dzi[k]

                    + ( rhorefh[k+1] * fabs(interp2(w[ijk-jj1+kk1], w[ijk+kk1])) * interp3_ws(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])) / rhoref[k] * dzi[k];
            }
            else if (k == kstart+2)
            {
                vt[ijk] +=
                    // w*dv/dz -> fourth order interpolation for fluxbot, sixth for fluxtop
                    - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp6_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2], v[ijk+kk3])
                      - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp4_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) / rhoref[k] * dzi[k]

                    + ( rhorefh[k+1] * fabs(interp2(w[ijk-jj1+kk1], w[ijk+kk1])) * interp5_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2], v[ijk+kk3])
                      - rhorefh[k  ] * fabs(interp2(w[ijk-jj1    ], w[ijk    ])) * interp3_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }
            else if (k == kend-3)
            {
                vt[ijk] +=
                    // w*dv/dz -> fourth order interpolation for fluxtop
                    - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp4_ws(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])
                      - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp6_ws(v[ijk-kk3], v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2]) ) / rhoref[k] * dzi[k]

                    + ( rhorefh[k+1] * fabs(interp2(w[ijk-jj1+kk1], w[ijk+kk1])) * interp3_ws(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])
                      - rhorefh[k  ] * fabs(interp2(w[ijk-jj1    ], w[ijk    ])) * interp5_ws(v[ijk-kk3], v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2]) ) / rhoref[k] * dzi[k];
            }
            else if (k == kend-2)
            {
                vt[ijk] +=
                    // w*dv/dz -> second order interpolation for fluxtop, fourth order for fluxbot
                    - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp2(v[ijk    ], v[ijk+kk1])
                      - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp4_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) / rhoref[k] * dzi[k]

                    - ( rhorefh[k  ] * fabs(interp2(w[ijk-jj1    ], w[ijk    ])) * interp3_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }
            else if (k == kend-1)
            {
                vt[ijk] +=
                    // w*dv/dz -> second order interpolation for fluxbot, fluxtop=0 as w=0
                    - ( -rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp2(v[ijk-kk1], v[ijk    ]) ) / rhoref[k] * dzi[k];
            }
            else
            {
                vt[ijk] +=
                    - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp6_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2], v[ijk+kk3])
                      - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp6_ws(v[ijk-kk3], v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2]) ) / rhoref[k] * dzi[k]

                    + ( rhorefh[k+1] * fabs(interp2(w[ijk-jj1+kk1], w[ijk+kk1])) * interp5_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2], v[ijk+kk3])
                      - rhorefh[k  ] * fabs(interp2(w[ijk-jj1    ], w[ijk    ])) * interp5_ws(v[ijk-kk3], v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2]) ) / rhoref[k] * dzi[k];
            }
        }
    }


    template<typename TF> __global__
    void advec_w_g(TF* __restrict__ wt, const TF* __restrict__ u,
                   const TF* __restrict__ v,  const TF* __restrict__ w,
                   const TF* __restrict__ rhoref, const TF* __restrict__ rhorefh,
                   const TF* __restrict__ dzhi, const TF dxi, const TF dyi,
                   const int jj, int kk,
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
                // u*dw/dx
                - ( interp2(u[ijk+ii1-kk], u[ijk+ii1]) * interp6_ws(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2], w[ijk+ii3])
                  - interp2(u[ijk    -kk], u[ijk    ]) * interp6_ws(w[ijk-ii3], w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2]) ) * dxi

                + ( fabs(interp2(u[ijk+ii1-kk], u[ijk+ii1])) * interp5_ws(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2], w[ijk+ii3])
                  - fabs(interp2(u[ijk    -kk], u[ijk    ])) * interp5_ws(w[ijk-ii3], w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2]) ) * dxi

                // v*dw/dy
                - ( interp2(v[ijk+jj1-kk], v[ijk+jj1]) * interp6_ws(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2], w[ijk+jj3])
                  - interp2(v[ijk    -kk], v[ijk    ]) * interp6_ws(w[ijk-jj3], w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2]) ) * dyi

                + ( fabs(interp2(v[ijk+jj1-kk], v[ijk+jj1])) * interp5_ws(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2], w[ijk+jj3])
                  - fabs(interp2(v[ijk    -kk], v[ijk    ])) * interp5_ws(w[ijk-jj3], w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2]) ) * dyi;

            if (k == kstart+1)
            {
                wt[ijk] +=
                    // w*dv/dz -> second order interpolation for fluxbot, fourth order for fluxtop
                    - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp4_ws(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2])
                      - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp2(w[ijk-kk1], w[ijk    ]) ) / rhorefh[k] * dzhi[k]

                    + ( rhoref[k  ] * fabs(interp2(w[ijk        ], w[ijk+kk1])) * interp3_ws(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) ) / rhorefh[k] * dzhi[k];
            }
            if (k == kstart+2)
            {
                wt[ijk] +=
                    // w*dv/dz -> fourth order interpolation for fluxbot, sixth order for fluxtop
                    - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp6_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3])
                      - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp4_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) / rhorefh[k] * dzhi[k]

                    + ( rhoref[k  ] * fabs(interp2(w[ijk        ], w[ijk+kk1])) * interp5_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3])
                      - rhoref[k-1] * fabs(interp2(w[ijk-kk1    ], w[ijk    ])) * interp3_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) / rhorefh[k] * dzhi[k];
            }
            else if (k == kend-2)
            {
                wt[ijk] +=
                    // w*dv/dz -> sixth order interpolation for fluxbot, fourth order for fluxtop
                    - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp4_ws(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2])
                      - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp6_ws(w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) ) / rhorefh[k] * dzhi[k]

                    + ( rhoref[k  ] * fabs(interp2(w[ijk        ], w[ijk+kk1])) * interp3_ws(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2])
                      - rhoref[k-1] * fabs(interp2(w[ijk-kk1    ], w[ijk    ])) * interp5_ws(w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) ) / rhorefh[k] * dzhi[k];
            }
            else if (k == kend-1)
            {
                wt[ijk] +=
                    // w*dv/dz -> fourth order interpolation for fluxbot, second order for fluxtop
                    - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp2(w[ijk    ], w[ijk+kk1])
                      - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp4_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) / rhorefh[k] * dzhi[k]

                    - ( rhoref[k-1] * fabs(interp2(w[ijk-kk1    ], w[ijk    ])) * interp3_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) / rhorefh[k] * dzhi[k];
            }
            else if ( (k >= kstart+2) && (k < kend-1) )
            {
                wt[ijk] +=
                    // w*dw/dz
                    - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp6_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3])
                      - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp6_ws(w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) ) / rhorefh[k] * dzhi[k]

                    + ( rhoref[k  ] * fabs(interp2(w[ijk        ], w[ijk+kk1])) * interp5_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3])
                      - rhoref[k-1] * fabs(interp2(w[ijk-kk1    ], w[ijk    ])) * interp5_ws(w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) ) / rhorefh[k] * dzhi[k];
            }
        }
    }


    template<typename TF> __global__
    void advec_s_g(TF* __restrict__ st, const TF* __restrict__ s, const TF* __restrict__ u,
                   const TF* __restrict__ v,  const TF* __restrict__ w,
                   const TF* __restrict__ rhoref, const TF* __restrict__ rhorefh,
                   const TF* __restrict__ dzi, const TF dxi, const TF dyi,
                   const int jj, int kk,
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

            st[ijk] +=
                - ( u[ijk+ii1] * interp4_ws(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                  - u[ijk    ] * interp4_ws(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                + ( fabs(u[ijk+ii1]) * interp3_ws(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                  - fabs(u[ijk    ]) * interp3_ws(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                - ( v[ijk+jj1] * interp4_ws(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                  - v[ijk    ] * interp4_ws(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi

                + ( fabs(v[ijk+jj1]) * interp3_ws(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                  - fabs(v[ijk    ]) * interp3_ws(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi;

            if (k == kstart)
            {
                st[ijk] +=
                    - ( rhorefh[k+1] * w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }
            else if (k == kstart+1)
            {
                st[ijk] +=
                    - ( rhorefh[k+1] * w[ijk+kk1] * interp4_ws(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                      - rhorefh[k  ] * w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]) ) / rhoref[k] * dzi[k];

                    + ( rhorefh[k+1] * fabs(w[ijk+kk1]) * interp3_ws(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2]) ) / rhoref[k] * dzi[k];
            }
            else if (k == kend-2)
            {
                st[ijk] +=
                    - ( rhorefh[k+1] * w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1])
                      - rhorefh[k  ] * w[ijk    ] * interp4_ws(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];

                    - ( rhorefh[k  ] * fabs(w[ijk    ]) * interp3_ws(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }
            else if (k == kend-1)
            {
                st[ijk] +=
                    - (-rhorefh[k  ] * w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]) ) / rhoref[k] * dzi[k];
            }
            else
            {
                st[ijk] +=
                    - ( rhorefh[k+1] * w[ijk+kk1] * interp4_ws(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                      - rhorefh[k  ] * w[ijk    ] * interp4_ws(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];

                    + ( rhorefh[k+1] * fabs(w[ijk+kk1]) * interp3_ws(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                      - rhorefh[k  ] * fabs(w[ijk    ]) * interp3_ws(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
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
        const int ii3 = 3;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int jj3 = 3*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;
        const int kk3 = 3*kk;

        const int ijk = i + j*jj + k*kk;

        if (i < iend && j < jend && k < kend)
        {
            if (k == kstart || k == kend-1)
                tmp1[ijk] = fabs(interp6_ws(u[ijk-ii2], u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]))*dxi
                          + fabs(interp6_ws(v[ijk-jj2], v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]))*dyi
                          + fabs(interp2(w[ijk], w[ijk+kk1]))*dzi[k];
            else if (k == kstart+1 || k == kend-2)
                tmp1[ijk] = fabs(interp6_ws(u[ijk-ii2], u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]))*dxi
                          + fabs(interp6_ws(v[ijk-jj2], v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]))*dyi
                          + fabs(interp4_ws(w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2]))*dzi[k];
            else
                tmp1[ijk] = fabs(interp6_ws(u[ijk-ii2], u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]))*dxi
                          + fabs(interp6_ws(v[ijk-jj2], v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]))*dyi
                          + fabs(interp6_ws(w[ijk-kk2], w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3]))*dzi[k];
        }
    }
}


#ifdef USECUDA
template<typename TF>
unsigned long Advec_2i5<TF>::get_time_limit(unsigned long idt, double dt)
{
    // Calculate cfl and prevent zero divisons.
    double cfl = get_cfl(dt);
    cfl = std::max(cflmin, cfl);
    const unsigned long idtlim = idt * cflmax / cfl;

    return idtlim;
}


template<typename TF>
double Advec_2i5<TF>::get_cfl(const double dt)
{
    const Grid_data<TF>& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    auto tmp1 = fields.get_tmp_g();

    calc_cfl_g<<<gridGPU, blockGPU>>>(
        tmp1->fld_g,
        fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
        gd.dzi_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend, gd.jend, gd.kend);
    cuda_check_error();

    TF cfl = field3d_operators.calc_max_g(tmp1->fld_g);
    fields.release_tmp_g(tmp1);

    cfl = cfl*dt;

    return static_cast<double>(cfl);
}


template<typename TF>
void Advec_2i5<TF>::exec(Stats<TF>& stats)
{
    const Grid_data<TF>& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    advec_u_g<TF><<<gridGPU, blockGPU>>>(
        fields.mt.at("u")->fld_g,
        fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
        fields.rhoref_g, fields.rhorefh_g, gd.dzi_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend, gd.jend, gd.kend);
    cuda_check_error();

    advec_v_g<TF><<<gridGPU, blockGPU>>>(
        fields.mt.at("v")->fld_g,
        fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
        fields.rhoref_g, fields.rhorefh_g, gd.dzi_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend, gd.jend, gd.kend);
    cuda_check_error();

    advec_w_g<TF><<<gridGPU, blockGPU>>>(
        fields.mt.at("w")->fld_g,
        fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
        fields.rhoref_g, fields.rhorefh_g, gd.dzhi_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend, gd.jend, gd.kend);
    cuda_check_error();

    for (auto& it : fields.st)
        advec_s_g<TF><<<gridGPU, blockGPU>>>(
            it.second->fld_g, fields.sp.at(it.first)->fld_g,
            fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
            fields.rhoref_g, fields.rhorefh_g, gd.dzi_g, gd.dxi, gd.dyi,
            gd.icells, gd.ijcells,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend, gd.jend, gd.kend);
    cuda_check_error();

    cudaDeviceSynchronize();
    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);
    for (auto it : fields.st)
        stats.calc_tend(*it.second, tend_name);
}
#endif


template class Advec_2i5<double>;
template class Advec_2i5<float>;
