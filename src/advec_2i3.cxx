/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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
#include "advec_2i3.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"

namespace
{
    using namespace Finite_difference::O2;
}

template<typename TF>
Advec_2i3<TF>::Advec_2i3(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Advec<TF>(masterin, gridin, fieldsin, inputin)
{
    swadvec = "2i3";
}

template<typename TF>
Advec_2i3<TF>::~Advec_2i3() {}

#ifndef USECUDA
namespace
{
    template<typename TF>
    inline TF interp4_ws(const TF a, const TF b, const TF c, const TF d) 
    {
        const TF c0 = 7./12.;
        const TF c1 = 1./12.;
        return c0*(b + c) - c1*(a + d);
    }

    template<typename TF>
    inline TF interp3_ws(const TF a, const TF b, const TF c, const TF d) 
    {
        const TF c0 = 3./12.;
        const TF c1 = 1./12.;
        return c0*(c - b) - c1*(d - a);
    }

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

        const TF dxi = 1./dx;
        const TF dyi = 1./dy;

        TF cfl = 0;
    
        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                cfl = std::max(cfl, std::abs(interp4_ws(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]))*dxi 
                                  + std::abs(interp4_ws(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]))*dyi 
                                  + std::abs(interp2(w[ijk    ], w[ijk+kk1]))*dzi[k]);
            }
    
        for (k=kstart+1; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    cfl = std::max(cfl, std::abs(interp4_ws(u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2]))*dxi 
                                      + std::abs(interp4_ws(v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2]))*dyi 
                                      + std::abs(interp4_ws(w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2]))*dzi[k]);
                }
    
        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk  = i + j*jj1 + k*kk1;
                cfl = std::max(cfl, std::abs(interp4_ws(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]))*dxi 
                                  + std::abs(interp4_ws(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]))*dyi 
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
            const TF* const restrict dzi, const TF dx, const TF dy,
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

        const TF dxi = 1./dx;
        const TF dyi = 1./dy;

        int k = kstart; 
    
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                ut[ijk] += 
                         // u*du/dx
                         - ( interp2(u[ijk        ], u[ijk+ii1]) * interp4_ws(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                           - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4_ws(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi
    
                         + ( std::abs(interp2(u[ijk        ], u[ijk+ii1])) * interp3_ws(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                           - std::abs(interp2(u[ijk-ii1    ], u[ijk    ])) * interp3_ws(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi
    
                         // v*du/dy
                         - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4_ws(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                           - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4_ws(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 
    
                         + ( std::abs(interp2(v[ijk-ii1+jj1], v[ijk+jj1])) * interp3_ws(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                           - std::abs(interp2(v[ijk-ii1    ], v[ijk    ])) * interp3_ws(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 
    
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
                         - ( interp2(u[ijk        ], u[ijk+ii1]) * interp4_ws(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                           - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4_ws(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi
    
                         + ( std::abs(interp2(u[ijk        ], u[ijk+ii1])) * interp3_ws(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                           - std::abs(interp2(u[ijk-ii1    ], u[ijk    ])) * interp3_ws(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi
    
                         // v*du/dy
                         - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4_ws(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                           - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4_ws(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 
    
                         + ( std::abs(interp2(v[ijk-ii1+jj1], v[ijk+jj1])) * interp3_ws(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                           - std::abs(interp2(v[ijk-ii1    ], v[ijk    ])) * interp3_ws(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 
    
                         // w*du/dz -> second order interpolation for fluxbot
                         - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp4_ws(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2])
                           - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp2(u[ijk-kk1], u[ijk    ]) ) / rhoref[k] * dzi[k]
    
                         + ( rhorefh[k+1] * std::abs(interp2(w[ijk-ii1+kk1], w[ijk+kk1])) * interp3_ws(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]) ) / rhoref[k] * dzi[k];
            }
    
        for (k=kstart+2; k<kend-2; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    ut[ijk] += 
                             // u*du/dx
                             - ( interp2(u[ijk        ], u[ijk+ii1]) * interp4_ws(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                               - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4_ws(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi
    
                             + ( std::abs(interp2(u[ijk        ], u[ijk+ii1])) * interp3_ws(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                               - std::abs(interp2(u[ijk-ii1    ], u[ijk    ])) * interp3_ws(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi
    
                             // v*du/dy
                             - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4_ws(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                               - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4_ws(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 
    
                             + ( std::abs(interp2(v[ijk-ii1+jj1], v[ijk+jj1])) * interp3_ws(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                               - std::abs(interp2(v[ijk-ii1    ], v[ijk    ])) * interp3_ws(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 
    
                             // w*du/dz
                             - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp4_ws(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2])
                               - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp4_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k]
    
                             + ( rhorefh[k+1] * std::abs(interp2(w[ijk-ii1+kk1], w[ijk+kk1])) * interp3_ws(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2])
                               - rhorefh[k  ] * std::abs(interp2(w[ijk-ii1    ], w[ijk    ])) * interp3_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k];
                }
    
        k = kend-2; 
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                ut[ijk] += 
                         // u*du/dx
                         - ( interp2(u[ijk        ], u[ijk+ii1]) * interp4_ws(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                           - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4_ws(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi
    
                         + ( std::abs(interp2(u[ijk        ], u[ijk+ii1])) * interp3_ws(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                           - std::abs(interp2(u[ijk-ii1    ], u[ijk    ])) * interp3_ws(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi
    
                         // v*du/dy
                         - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4_ws(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                           - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4_ws(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 
    
                         + ( std::abs(interp2(v[ijk-ii1+jj1], v[ijk+jj1])) * interp3_ws(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                           - std::abs(interp2(v[ijk-ii1    ], v[ijk    ])) * interp3_ws(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 
    
                         // w*du/dz -> second order interpolation for fluxtop
                         - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp2(u[ijk    ], u[ijk+kk1])
                           - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp4_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k]
    
                         - ( rhorefh[k  ] * std::abs(interp2(w[ijk-ii1    ], w[ijk    ])) * interp3_ws(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }
    
        k = kend-1; 
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                ut[ijk] += 
                         // u*du/dx
                         - ( interp2(u[ijk        ], u[ijk+ii1]) * interp4_ws(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                           - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4_ws(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi
    
                         + ( std::abs(interp2(u[ijk        ], u[ijk+ii1])) * interp3_ws(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                           - std::abs(interp2(u[ijk-ii1    ], u[ijk    ])) * interp3_ws(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi
    
                         // v*du/dy
                         - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4_ws(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                           - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4_ws(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 
    
                         + ( std::abs(interp2(v[ijk-ii1+jj1], v[ijk+jj1])) * interp3_ws(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                           - std::abs(interp2(v[ijk-ii1    ], v[ijk    ])) * interp3_ws(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 
    
                         // w*du/dz -> second order interpolation for fluxbot, fluxtop=0 as w=0
                         - ( -rhorefh[k] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp2(u[ijk-kk1], u[ijk    ]) ) / rhoref[k] * dzi[k];
            }
    }

    template<typename TF>
    void advec_v(
            TF* const restrict vt,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict dzi, const TF dx, const TF dy,
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

        const TF dxi = 1./dx;
        const TF dyi = 1./dy;

        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                vt[ijk] += 
                         // u*dv/dx
                         - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4_ws(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                           - interp2(u[ijk    -jj1], u[ijk    ]) * interp4_ws(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi
    
                         + ( std::abs(interp2(u[ijk+ii1-jj1], u[ijk+ii1])) * interp3_ws(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                           - std::abs(interp2(u[ijk    -jj1], u[ijk    ])) * interp3_ws(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi
    
                         // v*dv/dy
                         - ( interp2(v[ijk        ], v[ijk+jj1]) * interp4_ws(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                           - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4_ws(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi
    
                         + ( std::abs(interp2(v[ijk        ], v[ijk+jj1])) * interp3_ws(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                           - std::abs(interp2(v[ijk-jj1    ], v[ijk    ])) * interp3_ws(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi
    
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
                         - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4_ws(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                           - interp2(u[ijk    -jj1], u[ijk    ]) * interp4_ws(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi
    
                         + ( std::abs(interp2(u[ijk+ii1-jj1], u[ijk+ii1])) * interp3_ws(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                           - std::abs(interp2(u[ijk    -jj1], u[ijk    ])) * interp3_ws(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi
    
                         // v*dv/dy
                         - ( interp2(v[ijk        ], v[ijk+jj1]) * interp4_ws(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                           - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4_ws(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi
    
                         + ( std::abs(interp2(v[ijk        ], v[ijk+jj1])) * interp3_ws(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                           - std::abs(interp2(v[ijk-jj1    ], v[ijk    ])) * interp3_ws(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi
    
                         // w*dv/dz
                         - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp4_ws(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])
                           - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp2(v[ijk-kk1], v[ijk    ]) ) / rhoref[k] * dzi[k]
    
                         + ( rhorefh[k+1] * std::abs(interp2(w[ijk-jj1+kk1], w[ijk+kk1])) * interp3_ws(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2]) ) / rhoref[k] * dzi[k];
            }
    
        for (k=kstart+2; k<kend-2; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    vt[ijk] += 
                             // u*dv/dx
                             - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4_ws(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                               - interp2(u[ijk    -jj1], u[ijk    ]) * interp4_ws(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi
    
                             + ( std::abs(interp2(u[ijk+ii1-jj1], u[ijk+ii1])) * interp3_ws(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                               - std::abs(interp2(u[ijk    -jj1], u[ijk    ])) * interp3_ws(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi
    
                             // v*dv/dy
                             - ( interp2(v[ijk        ], v[ijk+jj1]) * interp4_ws(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                               - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4_ws(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi
    
                             + ( std::abs(interp2(v[ijk        ], v[ijk+jj1])) * interp3_ws(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                               - std::abs(interp2(v[ijk-jj1    ], v[ijk    ])) * interp3_ws(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi
    
                             // w*dv/dz
                             - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp4_ws(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])
                               - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp4_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) / rhoref[k] * dzi[k]
    
                             + ( rhorefh[k+1] * std::abs(interp2(w[ijk-jj1+kk1], w[ijk+kk1])) * interp3_ws(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])
                               - rhorefh[k  ] * std::abs(interp2(w[ijk-jj1    ], w[ijk    ])) * interp3_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) / rhoref[k] * dzi[k];
                }
    
        k = kend-2;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                vt[ijk] += 
                         // u*dv/dx
                         - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4_ws(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                           - interp2(u[ijk    -jj1], u[ijk    ]) * interp4_ws(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi
    
                         + ( std::abs(interp2(u[ijk+ii1-jj1], u[ijk+ii1])) * interp3_ws(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                           - std::abs(interp2(u[ijk    -jj1], u[ijk    ])) * interp3_ws(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi
    
                         // v*dv/dy
                         - ( interp2(v[ijk        ], v[ijk+jj1]) * interp4_ws(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                           - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4_ws(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi
    
                         + ( std::abs(interp2(v[ijk        ], v[ijk+jj1])) * interp3_ws(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                           - std::abs(interp2(v[ijk-jj1    ], v[ijk    ])) * interp3_ws(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi
    
                         // w*dv/dz
                         - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp2(v[ijk    ], v[ijk+kk1])
                           - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp4_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) / rhoref[k] * dzi[k]
    
                         - ( rhorefh[k  ] * std::abs(interp2(w[ijk-jj1    ], w[ijk    ])) * interp3_ws(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }
    
        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                vt[ijk] +=
                         // u*dv/dx
                         - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4_ws(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                           - interp2(u[ijk    -jj1], u[ijk    ]) * interp4_ws(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi
    
                         + ( std::abs(interp2(u[ijk+ii1-jj1], u[ijk+ii1])) * interp3_ws(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                           - std::abs(interp2(u[ijk    -jj1], u[ijk    ])) * interp3_ws(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi
    
                         // v*dv/dy
                         - ( interp2(v[ijk        ], v[ijk+jj1]) * interp4_ws(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                           - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4_ws(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi
    
                         + ( std::abs(interp2(v[ijk        ], v[ijk+jj1])) * interp3_ws(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                           - std::abs(interp2(v[ijk-jj1    ], v[ijk    ])) * interp3_ws(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi
    
                         // w*dv/dz
                         - (- rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp2(v[ijk-kk1], v[ijk    ]) ) / rhoref[k] * dzi[k];
            }
    }

    template<typename TF>
    void advec_w(
            TF* const restrict wt,
            const TF* const restrict u, const TF* const restrict v, TF* const restrict w,
            const TF* const restrict dzhi, const TF dx, const TF dy,
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

        const TF dxi = 1./dx;
        const TF dyi = 1./dy;

        int k = kstart+1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                wt[ijk] += 
                         // u*dw/dx 
                         - ( interp2(u[ijk+ii1-kk1], u[ijk+ii1]) * interp4_ws(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
                           - interp2(u[ijk    -kk1], u[ijk    ]) * interp4_ws(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi
    
                         + ( std::abs(interp2(u[ijk+ii1-kk1], u[ijk+ii1])) * interp3_ws(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
                           - std::abs(interp2(u[ijk    -kk1], u[ijk    ])) * interp3_ws(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi
    
                         // v*dw/dy 
                         - ( interp2(v[ijk+jj1-kk1], v[ijk+jj1]) * interp4_ws(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
                           - interp2(v[ijk    -kk1], v[ijk    ]) * interp4_ws(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi
    
                         + ( std::abs(interp2(v[ijk+jj1-kk1], v[ijk+jj1])) * interp3_ws(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
                           - std::abs(interp2(v[ijk    -kk1], v[ijk    ])) * interp3_ws(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi
    
                         // w*dw/dz 
                         - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp4_ws(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2])
                           - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp2(w[ijk-kk1], w[ijk    ]) ) / rhorefh[k] * dzhi[k]
    
                         + ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp4_ws(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) ) / rhorefh[k] * dzhi[k];
            }
    
        for (k=kstart+2; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    wt[ijk] +=
                             // u*dw/dx 
                             - ( interp2(u[ijk+ii1-kk1], u[ijk+ii1]) * interp4_ws(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
                               - interp2(u[ijk    -kk1], u[ijk    ]) * interp4_ws(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi
    
                             + ( std::abs(interp2(u[ijk+ii1-kk1], u[ijk+ii1])) * interp3_ws(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
                               - std::abs(interp2(u[ijk    -kk1], u[ijk    ])) * interp3_ws(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi
    
                             // v*dw/dy 
                             - ( interp2(v[ijk+jj1-kk1], v[ijk+jj1]) * interp4_ws(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
                               - interp2(v[ijk    -kk1], v[ijk    ]) * interp4_ws(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi
    
                             + ( std::abs(interp2(v[ijk+jj1-kk1], v[ijk+jj1])) * interp3_ws(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
                               - std::abs(interp2(v[ijk    -kk1], v[ijk    ])) * interp3_ws(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi
    
                             // w*dw/dz 
                             - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp4_ws(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2])
                               - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp4_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) / rhorefh[k] * dzhi[k]
    
                             + ( rhoref[k  ] * std::abs(interp2(w[ijk        ], w[ijk+kk1])) * interp3_ws(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2])
                               - rhoref[k-1] * std::abs(interp2(w[ijk-kk1    ], w[ijk    ])) * interp3_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) / rhorefh[k] * dzhi[k];
                }
    
        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                wt[ijk] += 
                         // u*dw/dx 
                         - ( interp2(u[ijk+ii1-kk1], u[ijk+ii1]) * interp4_ws(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
                           - interp2(u[ijk    -kk1], u[ijk    ]) * interp4_ws(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi
    
                         + ( std::abs(interp2(u[ijk+ii1-kk1], u[ijk+ii1])) * interp3_ws(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
                           - std::abs(interp2(u[ijk    -kk1], u[ijk    ])) * interp3_ws(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi
    
                         // v*dw/dy 
                         - ( interp2(v[ijk+jj1-kk1], v[ijk+jj1]) * interp4_ws(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
                           - interp2(v[ijk    -kk1], v[ijk    ]) * interp4_ws(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi
    
                         + ( std::abs(interp2(v[ijk+jj1-kk1], v[ijk+jj1])) * interp3_ws(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
                           - std::abs(interp2(v[ijk    -kk1], v[ijk    ])) * interp3_ws(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi
    
                         // w*dw/dz 
                         - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp2(w[ijk    ], w[ijk+kk1])
                           - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp4_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) / rhorefh[k] * dzhi[k]
    
                         - ( rhoref[k-1] * std::abs(interp2(w[ijk-kk1    ], w[ijk    ])) * interp3_ws(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) / rhorefh[k] * dzhi[k];
            }
    }

    template<typename TF>
    void advec_s(
            TF* const restrict st, const TF* const restrict s,
            const TF* const restrict u, const TF* const restrict v, const TF* const restrict w,
            const TF* const restrict dzi, const TF dx, const TF dy,
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

        const TF dxi = 1./dx;
        const TF dyi = 1./dy;

        // assume that w at the boundary equals zero...
        int k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                st[ijk] += 
                         - ( u[ijk+ii1] * interp4_ws(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                           - u[ijk    ] * interp4_ws(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi
    
                         + ( std::abs(u[ijk+ii1]) * interp3_ws(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                           - std::abs(u[ijk    ]) * interp3_ws(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi
    
                         - ( v[ijk+jj1] * interp4_ws(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                           - v[ijk    ] * interp4_ws(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 
    
                         + ( std::abs(v[ijk+jj1]) * interp3_ws(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                           - std::abs(v[ijk    ]) * interp3_ws(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 
    
                         - ( rhorefh[k+1] * w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }
    
        k = kstart+1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                st[ijk] += 
                         - ( u[ijk+ii1] * interp4_ws(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                           - u[ijk    ] * interp4_ws(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi
    
                         + ( std::abs(u[ijk+ii1]) * interp3_ws(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                           - std::abs(u[ijk    ]) * interp3_ws(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi
    
                         - ( v[ijk+jj1] * interp4_ws(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                           - v[ijk    ] * interp4_ws(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 
    
                         + ( std::abs(v[ijk+jj1]) * interp3_ws(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                           - std::abs(v[ijk    ]) * interp3_ws(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 
    
                         - ( rhorefh[k+1] * w[ijk+kk1] * interp4_ws(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                           - rhorefh[k  ] * w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]) ) / rhoref[k] * dzi[k]
    
                         + ( rhorefh[k+1] * std::abs(w[ijk+kk1]) * interp3_ws(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2]) ) / rhoref[k] * dzi[k];
            }
    
        for (k=kstart+2; k<kend-2; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj1 + k*kk1;
                    st[ijk] += 
                             - ( u[ijk+ii1] * interp4_ws(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                               - u[ijk    ] * interp4_ws(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi
    
                             + ( std::abs(u[ijk+ii1]) * interp3_ws(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                               - std::abs(u[ijk    ]) * interp3_ws(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi
    
                             - ( v[ijk+jj1] * interp4_ws(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                               - v[ijk    ] * interp4_ws(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 
    
                             + ( std::abs(v[ijk+jj1]) * interp3_ws(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                               - std::abs(v[ijk    ]) * interp3_ws(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 
    
                             - ( rhorefh[k+1] * w[ijk+kk1] * interp4_ws(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                               - rhorefh[k  ] * w[ijk    ] * interp4_ws(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k]
    
                             + ( rhorefh[k+1] * std::abs(w[ijk+kk1]) * interp3_ws(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                               - rhorefh[k  ] * std::abs(w[ijk    ]) * interp3_ws(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
                }
    
        k = kend-2;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                st[ijk] += 
                         - ( u[ijk+ii1] * interp4_ws(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                           - u[ijk    ] * interp4_ws(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi
    
                         + ( std::abs(u[ijk+ii1]) * interp3_ws(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                           - std::abs(u[ijk    ]) * interp3_ws(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi
    
                         - ( v[ijk+jj1] * interp4_ws(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                           - v[ijk    ] * interp4_ws(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 
    
                         + ( std::abs(v[ijk+jj1]) * interp3_ws(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                           - std::abs(v[ijk    ]) * interp3_ws(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 
    
                         - ( rhorefh[k+1] * w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1])
                           - rhorefh[k  ] * w[ijk    ] * interp4_ws(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }
    
        // assume that w at the boundary equals zero...
        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                st[ijk] += 
                         - ( u[ijk+ii1] * interp4_ws(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                           - u[ijk    ] * interp4_ws(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi
    
                         + ( std::abs(u[ijk+ii1]) * interp3_ws(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                           - std::abs(u[ijk    ]) * interp3_ws(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi
    
                         - ( v[ijk+jj1] * interp4_ws(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                           - v[ijk    ] * interp4_ws(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 
    
                         + ( std::abs(v[ijk+jj1]) * interp3_ws(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                           - std::abs(v[ijk    ]) * interp3_ws(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 
    
                         - (- rhorefh[k  ] * w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]) ) / rhoref[k] * dzi[k];
            }
        }
}

template<typename TF>
double Advec_2i3<TF>::get_cfl(double dt)
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
unsigned long Advec_2i3<TF>::get_time_limit(unsigned long idt, double dt)
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
void Advec_2i3<TF>::exec()
{
    auto& gd = grid.get_grid_data();
    advec_u(fields.mt.at("u")->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dx, gd.dy,
            fields.rhoref.data(), fields.rhorefh.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    advec_v(fields.mt.at("v")->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            gd.dzi.data(), gd.dx, gd.dy,
            fields.rhoref.data(), fields.rhorefh.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    advec_w(fields.mt.at("w")->fld.data(),
            fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
            gd.dzhi.data(), gd.dx, gd.dy,
            fields.rhoref.data(), fields.rhorefh.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    for (auto& it : fields.st)
        advec_s(it.second->fld.data(), fields.sp.at(it.first)->fld.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), fields.mp.at("w")->fld.data(),
                gd.dzi.data(), gd.dx, gd.dy,
                fields.rhoref.data(), fields.rhorefh.data(),
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
}
#endif

template class Advec_2i3<double>;
template class Advec_2i3<float>;
