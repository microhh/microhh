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
#include "advec_2i4.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"

template<typename TF>
Advec_2i4<TF>::Advec_2i4(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Advec<TF>(masterin, gridin, fieldsin, inputin)
{
    const int igc = 2;
    const int jgc = 2;
    const int kgc = 2;
    grid.set_minimum_ghost_cells(igc, jgc, kgc);
}

template<typename TF>
Advec_2i4<TF>::~Advec_2i4() {}

namespace
{
    using namespace Finite_difference::O2;
    using Finite_difference::O4::interp4c;

    template<typename TF>
    TF calc_cfl(
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict dzi, const TF dxi, const TF dyi,
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

        TF cfl = 0;

        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                cfl = std::max(cfl, std::abs(interp4c(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]))*dxi
                                  + std::abs(interp4c(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]))*dyi
                                  + std::abs(interp2(w[ijk    ], w[ijk+kk1]))*dzi[k]);
            }

        for (k=kstart+1; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    cfl = std::max(cfl, std::abs(interp4c(u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2]))*dxi
                                      + std::abs(interp4c(v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2]))*dyi
                                      + std::abs(interp4c(w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2]))*dzi[k]);
                }

        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk  = i + j*jj1 + k*kk1;
                cfl = std::max(cfl, std::abs(interp4c(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]))*dxi
                                  + std::abs(interp4c(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]))*dyi
                                  + std::abs(interp2(w[ijk    ], w[ijk+kk1]))*dzi[k]);
            }

        master.max(&cfl, 1);

        cfl = cfl*dt;

        return cfl;
    }

    template<typename TF>
    void advec_u(
            TF* const restrict ut,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict dzi, const TF dxi, const TF dyi,
            const TF* const restrict rhoref, const TF* const restrict rhorefh,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        int k = kstart;

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                ut[ijk] +=
                         // u*du/dx
                         - ( interp2(u[ijk        ], u[ijk+ii1]) * interp4c(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                           - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4c(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

                         // v*du/dy
                         - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4c(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                           - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4c(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi

                         // w*du/dz -> second order interpolation for fluxtop, fluxbot = 0. as w=0
                         - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp2(u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }

        k = kstart+1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                ut[ijk] +=
                         // u*du/dx
                         - ( interp2(u[ijk        ], u[ijk+ii1]) * interp4c(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                           - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4c(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

                         // v*du/dy
                         - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4c(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                           - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4c(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi

                         // w*du/dz -> second order interpolation for fluxbot
                         - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp4c(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2])
                           - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp2(u[ijk-kk1], u[ijk    ]) ) / rhoref[k] * dzi[k];
            }

        for (k=kstart+2; k<kend-2; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    ut[ijk] +=
                             // u*du/dx
                             - ( interp2(u[ijk        ], u[ijk+ii1]) * interp4c(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                               - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4c(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

                             // v*du/dy
                             - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4c(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                               - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4c(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi

                             // w*du/dz
                             - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp4c(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2])
                               - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp4c(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k];
                }

        k = kend-2;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                ut[ijk] +=
                         // u*du/dx
                         - ( interp2(u[ijk        ], u[ijk+ii1]) * interp4c(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                           - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4c(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

                         // v*du/dy
                         - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4c(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                           - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4c(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi

                         // w*du/dz -> second order interpolation for fluxtop
                         - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp2(u[ijk    ], u[ijk+kk1])
                           - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp4c(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }

        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                ut[ijk] +=
                         // u*du/dx
                         - ( interp2(u[ijk        ], u[ijk+ii1]) * interp4c(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                           - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4c(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

                         // v*du/dy
                         - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4c(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                           - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4c(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi

                         // w*du/dz -> second order interpolation for fluxbot, fluxtop=0 as w=0
                         - ( -rhorefh[k] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp2(u[ijk-kk1], u[ijk    ]) ) / rhoref[k] * dzi[k];
            }
    }

    template<typename TF>
    void advec_v(
            TF* const restrict vt,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict dzi, const TF dxi, const TF dyi,
            const TF* const restrict rhoref, const TF* const restrict rhorefh,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                vt[ijk] +=
                         // u*dv/dx
                         - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4c(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                           - interp2(u[ijk    -jj1], u[ijk    ]) * interp4c(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

                         // v*dv/dy
                         - ( interp2(v[ijk        ], v[ijk+jj1]) * interp4c(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                           - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4c(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

                         // w*dv/dz
                         - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp2(v[ijk    ], v[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }

        k = kstart+1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                vt[ijk] +=
                         // u*dv/dx
                         - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4c(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                           - interp2(u[ijk    -jj1], u[ijk    ]) * interp4c(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

                         // v*dv/dy
                         - ( interp2(v[ijk        ], v[ijk+jj1]) * interp4c(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                           - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4c(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

                         // w*dv/dz
                         - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp4c(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])
                           - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp2(v[ijk-kk1], v[ijk    ]) ) / rhoref[k] * dzi[k];
            }

        for (k=kstart+2; k<kend-2; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    vt[ijk] +=
                             // u*dv/dx
                             - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4c(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                               - interp2(u[ijk    -jj1], u[ijk    ]) * interp4c(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

                             // v*dv/dy
                             - ( interp2(v[ijk        ], v[ijk+jj1]) * interp4c(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                               - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4c(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

                             // w*dv/dz
                             - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp4c(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])
                               - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp4c(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) / rhoref[k] * dzi[k];
                }

        k = kend-2;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                vt[ijk] +=
                         // u*dv/dx
                         - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4c(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                           - interp2(u[ijk    -jj1], u[ijk    ]) * interp4c(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

                         // v*dv/dy
                         - ( interp2(v[ijk        ], v[ijk+jj1]) * interp4c(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                           - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4c(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

                         // w*dv/dz
                         - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp2(v[ijk    ], v[ijk+kk1])
                           - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp4c(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }

        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                vt[ijk] +=
                         // u*dv/dx
                         - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4c(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                           - interp2(u[ijk    -jj1], u[ijk    ]) * interp4c(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

                         // v*dv/dy
                         - ( interp2(v[ijk        ], v[ijk+jj1]) * interp4c(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                           - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4c(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

                         // w*dv/dz
                         - (- rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp2(v[ijk-kk1], v[ijk    ]) ) / rhoref[k] * dzi[k];
            }
    }

    template<typename TF>
    void advec_w(
            TF* const restrict wt,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict dzhi, const TF dxi, const TF dyi,
            const TF* const restrict rhoref, const TF* const restrict rhorefh,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        int k = kstart+1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                wt[ijk] +=
                         // u*dw/dx
                         - ( interp2(u[ijk+ii1-kk1], u[ijk+ii1]) * interp4c(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
                           - interp2(u[ijk    -kk1], u[ijk    ]) * interp4c(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi

                         // v*dw/dy
                         - ( interp2(v[ijk+jj1-kk1], v[ijk+jj1]) * interp4c(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
                           - interp2(v[ijk    -kk1], v[ijk    ]) * interp4c(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi

                         // w*dw/dz
                         - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp4c(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2])
                           - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp2(w[ijk-kk1], w[ijk    ]) ) / rhorefh[k] * dzhi[k];
            }

        for (k=kstart+2; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    wt[ijk] +=
                             // u*dw/dx
                             - ( interp2(u[ijk+ii1-kk1], u[ijk+ii1]) * interp4c(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
                               - interp2(u[ijk    -kk1], u[ijk    ]) * interp4c(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi

                             // v*dw/dy
                             - ( interp2(v[ijk+jj1-kk1], v[ijk+jj1]) * interp4c(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
                               - interp2(v[ijk    -kk1], v[ijk    ]) * interp4c(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi

                             // w*dw/dz
                             - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp4c(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2])
                               - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp4c(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) / rhorefh[k] * dzhi[k];
                }

        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                wt[ijk] +=
                         // u*dw/dx
                         - ( interp2(u[ijk+ii1-kk1], u[ijk+ii1]) * interp4c(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
                           - interp2(u[ijk    -kk1], u[ijk    ]) * interp4c(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi

                         // v*dw/dy
                         - ( interp2(v[ijk+jj1-kk1], v[ijk+jj1]) * interp4c(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
                           - interp2(v[ijk    -kk1], v[ijk    ]) * interp4c(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi

                         // w*dw/dz
                         - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp2(w[ijk    ], w[ijk+kk1])
                           - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp4c(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) / rhorefh[k] * dzhi[k];
            }
    }

    template<typename TF>
    void advec_s(
            TF* const restrict st, const TF* const restrict s,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict dzi, const TF dxi, const TF dyi,
            const TF* const restrict rhoref, const TF* const restrict rhorefh,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii1 = 1;
        const int ii2 = 2;
        const int jj1 = 1*jj;
        const int jj2 = 2*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        // assume that w at the boundary equals zero...
        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                st[ijk] +=
                         - ( u[ijk+ii1] * interp4c(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                           - u[ijk    ] * interp4c(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                         - ( v[ijk+jj1] * interp4c(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                           - v[ijk    ] * interp4c(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi

                         - ( rhorefh[k+1] * w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }

        k = kstart+1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                st[ijk] +=
                         - ( u[ijk+ii1] * interp4c(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                           - u[ijk    ] * interp4c(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                         - ( v[ijk+jj1] * interp4c(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                           - v[ijk    ] * interp4c(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi

                         - ( rhorefh[k+1] * w[ijk+kk1] * interp4c(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                           - rhorefh[k  ] * w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]) ) / rhoref[k] * dzi[k];
            }

        for (k=kstart+2; k<kend-2; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    st[ijk] +=
                             - ( u[ijk+ii1] * interp4c(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                               - u[ijk    ] * interp4c(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                             - ( v[ijk+jj1] * interp4c(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                               - v[ijk    ] * interp4c(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi

                             - ( rhorefh[k+1] * w[ijk+kk1] * interp4c(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                               - rhorefh[k  ] * w[ijk    ] * interp4c(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
                }

        k = kend-2;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                st[ijk] +=
                         - ( u[ijk+ii1] * interp4c(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                           - u[ijk    ] * interp4c(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                         - ( v[ijk+jj1] * interp4c(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                           - v[ijk    ] * interp4c(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi

                         - ( rhorefh[k+1] * w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1])
                           - rhorefh[k  ] * w[ijk    ] * interp4c(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }

        // assume that w at the boundary equals zero...
        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                st[ijk] +=
                         - ( u[ijk+ii1] * interp4c(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                           - u[ijk    ] * interp4c(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                         - ( v[ijk+jj1] * interp4c(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                           - v[ijk    ] * interp4c(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi

                         - (- rhorefh[k  ] * w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]) ) / rhoref[k] * dzi[k];
            }
        }

    template<typename TF>
    void advec_flux_u(
            TF* const restrict st, const TF* const restrict s, const TF* const restrict w,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii1 = 1;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        int k = kstart+1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                st[ijk-kk1] = 0; // Impose no flux through bottom wall.
                st[ijk] = interp2(w[ijk-ii1], w[ijk]) * interp2(s[ijk-kk1], s[ijk]);
            }

        for (int k=kstart+2; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] = interp2(w[ijk-ii1], w[ijk]) * interp4c(s[ijk-kk2], s[ijk-kk1], s[ijk], s[ijk+kk1]);
                }

        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                st[ijk] = interp2(w[ijk-ii1], w[ijk]) * interp2(s[ijk-kk1], s[ijk]);
                st[ijk+kk1] = 0; // Impose no flux through top wall.
            }
    }

    template<typename TF>
    void advec_flux_v(
            TF* const restrict st, const TF* const restrict s, const TF* const restrict w,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int jj1 = 1*jj;
        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        int k = kstart+1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                st[ijk-kk1] = 0; // Impose no flux through bottom wall.
                st[ijk] = interp2(w[ijk-jj1], w[ijk]) * interp2(s[ijk-kk1], s[ijk]);
            }

        for (int k=kstart+2; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] = interp2(w[ijk-jj1], w[ijk]) * interp4c(s[ijk-kk2], s[ijk-kk1], s[ijk], s[ijk+kk1]);
                }

        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                st[ijk] = interp2(w[ijk-jj1], w[ijk]) * interp2(s[ijk-kk1], s[ijk]);
                st[ijk+kk1] = 0; // Impose no flux through top wall.
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

        int k = kstart+1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                st[ijk-kk1] = 0; // Impose no flux through bottom wall.
                st[ijk] = w[ijk] * interp2(s[ijk-kk1], s[ijk]);
            }

        for (int k=kstart+2; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    st[ijk] = w[ijk] * interp4c(s[ijk-kk2], s[ijk-kk1], s[ijk], s[ijk+kk1]);
                }

        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                st[ijk] = w[ijk] * interp2(s[ijk-kk1], s[ijk]);
                st[ijk+kk1] = 0; // Impose no flux through top wall.
            }
    }
}

template<typename TF>
void Advec_2i4<TF>::create(Stats<TF>& stats)
{
    stats.add_tendency(*fields.mt.at("u"), "z", tend_name, tend_longname);
    stats.add_tendency(*fields.mt.at("v"), "z", tend_name, tend_longname);
    stats.add_tendency(*fields.mt.at("w"), "zh", tend_name, tend_longname);
    for (auto it : fields.st)
        stats.add_tendency(*it.second, "z", tend_name, tend_longname);
}

#ifndef USECUDA
template<typename TF>
double Advec_2i4<TF>::get_cfl(double dt)
{
    auto& gd = grid.get_grid_data();
    TF cfl = calc_cfl<TF>(
            fields.mp.at("u")->fld.data(),fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dxi, gd.dyi,
            dt, master,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    // CFL is kept in double precision for time stepping accuracy
    return static_cast<double>(cfl);
}

template<typename TF>
unsigned long Advec_2i4<TF>::get_time_limit(unsigned long idt, double dt)
{
    // Calculate cfl and prevent zero divisons.
    auto& gd = grid.get_grid_data();

    double cfl = calc_cfl<TF>(
            fields.mp.at("u")->fld.data(),fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dxi, gd.dyi,
            dt, master,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    cfl = std::max(cflmin, cfl);
    return idt * cflmax / cfl;
}

template<typename TF>
void Advec_2i4<TF>::exec(Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();
    advec_u(fields.mt.at("u")->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dxi, gd.dyi,
            fields.rhoref.data(), fields.rhorefh.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    advec_v(fields.mt.at("v")->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dxi, gd.dyi,
            fields.rhoref.data(), fields.rhorefh.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    advec_w(fields.mt.at("w")->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            gd.dzhi.data(), gd.dxi, gd.dyi,
            fields.rhoref.data(), fields.rhorefh.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    for (auto& it : fields.st)
        advec_s(it.second->fld.data(), fields.sp.at(it.first)->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                gd.dzi.data(), gd.dxi, gd.dyi,
                fields.rhoref.data(), fields.rhorefh.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);
    for (auto it : fields.st)
        stats.calc_tend(*it.second, tend_name);
}
#endif

template<typename TF>
void Advec_2i4<TF>::get_advec_flux(Field3d<TF>& advec_flux, const Field3d<TF>& fld)
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

template class Advec_2i4<double>;
template class Advec_2i4<float>;
