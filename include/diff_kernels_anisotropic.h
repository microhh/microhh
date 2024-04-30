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

#ifndef DIFF_KERNELS_ANISO_H
#define DIFF_KERNELS_ANISO_H

#include "fast_math.h"
#include "boundary.h"
#include "constants.h"

namespace Diff_kernels_anisotropic
{
    namespace fm = Fast_math;

    template <typename TF>
    void diff_u(
            TF* const restrict ut,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict dzi,
            const TF* const restrict dzhi,
            const TF dxi, const TF dyi,
            const TF* const restrict evisc_h,
            const TF* const restrict evisc_v,
            const TF* const restrict fluxbot,
            const TF* const restrict fluxtop,
            const TF* const restrict rhoref,
            const TF* const restrict rhorefh,
            const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii = 1;

        // bottom boundary
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;

                // Horizontal evisc.
                const TF evisce_h = evisc_h[ijk   ] + visc;
                const TF eviscw_h = evisc_h[ijk-ii] + visc;
                const TF eviscn_h = TF(0.25)*(evisc_h[ijk-ii   ] + evisc_h[ijk   ] + evisc_h[ijk-ii+jj] + evisc_h[ijk+jj]) + visc;
                const TF eviscs_h = TF(0.25)*(evisc_h[ijk-ii-jj] + evisc_h[ijk-jj] + evisc_h[ijk-ii   ] + evisc_h[ijk   ]) + visc;
                const TF evisct_h = TF(0.25)*(evisc_h[ijk-ii   ] + evisc_h[ijk   ] + evisc_h[ijk-ii+kk] + evisc_h[ijk+kk]) + visc;

                // Vertical evisc.
                const TF evisct_v = TF(0.25)*(evisc_v[ijk-ii   ] + evisc_v[ijk   ] + evisc_v[ijk-ii+kk] + evisc_v[ijk+kk]) + visc;

                ut[ijk] +=
                         // du/dx + du/dx
                         + ( evisce_h*(u[ijk+ii]-u[ijk   ])*dxi
                           - eviscw_h*(u[ijk   ]-u[ijk-ii])*dxi ) * TF(2.)*dxi
                         // du/dy + dv/dx
                         + ( eviscn_h*((u[ijk+jj]-u[ijk   ])*dyi + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                           - eviscs_h*((u[ijk   ]-u[ijk-jj])*dyi + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                         // du/dz + dw/dx
                         + ( rhorefh[kstart+1] * (evisct_v*(u[ijk+kk]-u[ijk])*dzhi[kstart+1] + evisct_h * (w[ijk+kk]-w[ijk-ii+kk])*dxi  )
                           + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
            }

        // top boundary
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + (kend-1)*kk;

                // Horizontal evisc.
                const TF evisce_h = evisc_h[ijk   ] + visc;
                const TF eviscw_h = evisc_h[ijk-ii] + visc;
                const TF eviscn_h = TF(0.25)*(evisc_h[ijk-ii   ] + evisc_h[ijk   ] + evisc_h[ijk-ii+jj] + evisc_h[ijk+jj]) + visc;
                const TF eviscs_h = TF(0.25)*(evisc_h[ijk-ii-jj] + evisc_h[ijk-jj] + evisc_h[ijk-ii   ] + evisc_h[ijk   ]) + visc;
                const TF eviscb_h = TF(0.25)*(evisc_h[ijk-ii-kk] + evisc_h[ijk-kk] + evisc_h[ijk-ii   ] + evisc_h[ijk   ]) + visc;

                // Vertical evisc.
                const TF eviscb_v = TF(0.25)*(evisc_v[ijk-ii-kk] + evisc_v[ijk-kk] + evisc_v[ijk-ii   ] + evisc_v[ijk   ]) + visc;

                ut[ijk] +=
                         // du/dx + du/dx
                         + ( evisce_h*(u[ijk+ii]-u[ijk   ])*dxi
                           - eviscw_h*(u[ijk   ]-u[ijk-ii])*dxi ) * TF(2.)*dxi
                         // du/dy + dv/dx
                         + ( eviscn_h*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                           - eviscs_h*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                         // du/dz + dw/dx
                         + (- rhorefh[kend  ] * fluxtop[ij]
                            - rhorefh[kend-1] * (eviscb_v*(u[ijk   ]-u[ijk-kk])*dzhi[kend-1] + eviscb_h*(w[ijk   ]-w[ijk-ii   ])*dxi))
                         / rhoref[kend-1] * dzi[kend-1];
            }

        for (int k=kstart+1; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // Horizontal evisc.
                    const TF evisce_h = evisc_h[ijk   ] + visc;
                    const TF eviscw_h = evisc_h[ijk-ii] + visc;
                    const TF eviscn_h = TF(0.25)*(evisc_h[ijk-ii   ] + evisc_h[ijk   ] + evisc_h[ijk-ii+jj] + evisc_h[ijk+jj]) + visc;
                    const TF eviscs_h = TF(0.25)*(evisc_h[ijk-ii-jj] + evisc_h[ijk-jj] + evisc_h[ijk-ii   ] + evisc_h[ijk   ]) + visc;
                    const TF evisct_h = TF(0.25)*(evisc_h[ijk-ii   ] + evisc_h[ijk   ] + evisc_h[ijk-ii+kk] + evisc_h[ijk+kk]) + visc;
                    const TF eviscb_h = TF(0.25)*(evisc_h[ijk-ii-kk] + evisc_h[ijk-kk] + evisc_h[ijk-ii   ] + evisc_h[ijk   ]) + visc;

                    // Vertical evisc.
                    const TF evisct_v = TF(0.25)*(evisc_v[ijk-ii   ] + evisc_v[ijk   ] + evisc_v[ijk-ii+kk] + evisc_v[ijk+kk]) + visc;
                    const TF eviscb_v = TF(0.25)*(evisc_v[ijk-ii-kk] + evisc_v[ijk-kk] + evisc_v[ijk-ii   ] + evisc_v[ijk   ]) + visc;

                    ut[ijk] +=
                             // du/dx + du/dx
                             + ( evisce_h*(u[ijk+ii]-u[ijk   ])*dxi
                               - eviscw_h*(u[ijk   ]-u[ijk-ii])*dxi ) * TF(2.)*dxi
                             // du/dy + dv/dx
                             + ( eviscn_h*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                               - eviscs_h*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                             // du/dz + dw/dx
                             + ( rhorefh[k+1] * (evisct_v*(u[ijk+kk]-u[ijk   ])*dzhi[k+1] + evisct_h*(w[ijk+kk]-w[ijk-ii+kk])*dxi)
                               - rhorefh[k  ] * (eviscb_v*(u[ijk   ]-u[ijk-kk])*dzhi[k  ] + eviscb_h*(w[ijk   ]-w[ijk-ii   ])*dxi) ) / rhoref[k] * dzi[k];
                }
    }


    template <typename TF>
    void diff_v(
            TF* const restrict vt,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict dzi,
            const TF* const restrict dzhi,
            const TF dxi, const TF dyi,
            const TF* const restrict evisc_h,
            const TF* const restrict evisc_v,
            const TF* const restrict fluxbot,
            const TF* const restrict fluxtop,
            const TF* const restrict rhoref,
            const TF* const restrict rhorefh,
            const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)

    {
        const int ii = 1;

        // bottom boundary
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;

                // Horizontal evisc.
                const TF evisce_h = TF(0.25)*(evisc_h[ijk   -jj] + evisc_h[ijk   ] + evisc_h[ijk+ii-jj] + evisc_h[ijk+ii]) + visc;
                const TF eviscw_h = TF(0.25)*(evisc_h[ijk-ii-jj] + evisc_h[ijk-ii] + evisc_h[ijk   -jj] + evisc_h[ijk   ]) + visc;
                const TF eviscn_h = evisc_h[ijk   ] + visc;
                const TF eviscs_h = evisc_h[ijk-jj] + visc;
                const TF evisct_h = TF(0.25)*(evisc_h[ijk   -jj] + evisc_h[ijk   ] + evisc_h[ijk+kk-jj] + evisc_h[ijk+kk]) + visc;

                // Vertical evisc.
                const TF evisct_v = TF(0.25)*(evisc_v[ijk   -jj] + evisc_v[ijk   ] + evisc_v[ijk+kk-jj] + evisc_v[ijk+kk]) + visc;

                vt[ijk] +=
                         // dv/dx + du/dy
                         + ( evisce_h*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                           - eviscw_h*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                         // dv/dy + dv/dy
                         + ( eviscn_h*(v[ijk+jj]-v[ijk   ])*dyi
                           - eviscs_h*(v[ijk   ]-v[ijk-jj])*dyi ) * TF(2.)*dyi
                         // dv/dz + dw/dy
                         + ( rhorefh[kstart+1] * (evisct_v*(v[ijk+kk]-v[ijk   ])*dzhi[kstart+1] + evisct_h*(w[ijk+kk]-w[ijk-jj+kk])*dyi)
                           + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
            }

        // top boundary
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + (kend-1)*kk;

                // Horizontal evisc.
                const TF evisce_h = TF(0.25)*(evisc_h[ijk   -jj] + evisc_h[ijk   ] + evisc_h[ijk+ii-jj] + evisc_h[ijk+ii]) + visc;
                const TF eviscw_h = TF(0.25)*(evisc_h[ijk-ii-jj] + evisc_h[ijk-ii] + evisc_h[ijk   -jj] + evisc_h[ijk   ]) + visc;
                const TF eviscn_h = evisc_h[ijk   ] + visc;
                const TF eviscs_h = evisc_h[ijk-jj] + visc;
                const TF eviscb_h = TF(0.25)*(evisc_h[ijk-kk-jj] + evisc_h[ijk-kk] + evisc_h[ijk   -jj] + evisc_h[ijk   ]) + visc;

                // Vertical evisc.
                const TF eviscb_v = TF(0.25)*(evisc_v[ijk-kk-jj] + evisc_v[ijk-kk] + evisc_v[ijk   -jj] + evisc_v[ijk   ]) + visc;

                vt[ijk] +=
                         // dv/dx + du/dy
                         + ( evisce_h*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                           - eviscw_h*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                         // dv/dy + dv/dy
                         + ( eviscn_h*(v[ijk+jj]-v[ijk   ])*dyi
                           - eviscs_h*(v[ijk   ]-v[ijk-jj])*dyi ) * TF(2.)*dyi
                         // dv/dz + dw/dy
                         + (- rhorefh[kend  ] * fluxtop[ij]
                            - rhorefh[kend-1] * (eviscb_v*(v[ijk   ]-v[ijk-kk])*dzhi[kend-1] + eviscb_h*(w[ijk   ]-w[ijk-jj   ])*dyi)
                            ) / rhoref[kend-1] * dzi[kend-1];
            }
        

        for (int k=kstart+1; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // Horizontal evisc.
                    const TF evisce_h = TF(0.25)*(evisc_h[ijk   -jj] + evisc_h[ijk   ] + evisc_h[ijk+ii-jj] + evisc_h[ijk+ii]) + visc;
                    const TF eviscw_h = TF(0.25)*(evisc_h[ijk-ii-jj] + evisc_h[ijk-ii] + evisc_h[ijk   -jj] + evisc_h[ijk   ]) + visc;
                    const TF eviscn_h = evisc_h[ijk   ] + visc;
                    const TF eviscs_h = evisc_h[ijk-jj] + visc;
                    const TF evisct_h = TF(0.25)*(evisc_h[ijk   -jj] + evisc_h[ijk   ] + evisc_h[ijk+kk-jj] + evisc_h[ijk+kk]) + visc;
                    const TF eviscb_h = TF(0.25)*(evisc_h[ijk-kk-jj] + evisc_h[ijk-kk] + evisc_h[ijk   -jj] + evisc_h[ijk   ]) + visc;

                    // Vertical evisc.
                    const TF evisct_v = TF(0.25)*(evisc_v[ijk   -jj] + evisc_v[ijk   ] + evisc_v[ijk+kk-jj] + evisc_v[ijk+kk]) + visc;
                    const TF eviscb_v = TF(0.25)*(evisc_v[ijk-kk-jj] + evisc_v[ijk-kk] + evisc_v[ijk   -jj] + evisc_v[ijk   ]) + visc;

                    vt[ijk] +=
                             // dv/dx + du/dy
                             + ( evisce_h*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                               - eviscw_h*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                             // dv/dy + dv/dy
                             + ( eviscn_h*(v[ijk+jj]-v[ijk   ])*dyi
                               - eviscs_h*(v[ijk   ]-v[ijk-jj])*dyi ) * TF(2.)*dyi
                             // dv/dz + dw/dy
                             + ( rhorefh[k+1] * (evisct_v*(v[ijk+kk]-v[ijk   ])*dzhi[k+1] + evisct_h*(w[ijk+kk]-w[ijk-jj+kk])*dyi)
                               - rhorefh[k  ] * (eviscb_v*(v[ijk   ]-v[ijk-kk])*dzhi[k  ] + eviscb_h*(w[ijk   ]-w[ijk-jj   ])*dyi) ) / rhoref[k] * dzi[k];
                }
    }


    template <typename TF>
    void diff_w(
            TF* const restrict wt,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict dzi,
            const TF* const restrict dzhi,
            const TF dxi, const TF dyi,
            const TF* const restrict evisc_h,
            const TF* const restrict evisc_v,
            const TF* const restrict rhoref,
            const TF* const restrict rhorefh,
            const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii = 1;

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // Horizontal evisc.
                    const TF evisce_h = TF(0.25)*(evisc_h[ijk   -kk] + evisc_h[ijk   ] + evisc_h[ijk+ii-kk] + evisc_h[ijk+ii]) + visc;
                    const TF eviscw_h = TF(0.25)*(evisc_h[ijk-ii-kk] + evisc_h[ijk-ii] + evisc_h[ijk   -kk] + evisc_h[ijk   ]) + visc;
                    const TF eviscn_h = TF(0.25)*(evisc_h[ijk   -kk] + evisc_h[ijk   ] + evisc_h[ijk+jj-kk] + evisc_h[ijk+jj]) + visc;
                    const TF eviscs_h = TF(0.25)*(evisc_h[ijk-jj-kk] + evisc_h[ijk-jj] + evisc_h[ijk   -kk] + evisc_h[ijk   ]) + visc;

                    // Vertical evisc.
                    const TF evisce_v = TF(0.25)*(evisc_v[ijk   -kk] + evisc_v[ijk   ] + evisc_v[ijk+ii-kk] + evisc_v[ijk+ii]) + visc;
                    const TF eviscw_v = TF(0.25)*(evisc_v[ijk-ii-kk] + evisc_v[ijk-ii] + evisc_v[ijk   -kk] + evisc_v[ijk   ]) + visc;
                    const TF eviscn_v = TF(0.25)*(evisc_v[ijk   -kk] + evisc_v[ijk   ] + evisc_v[ijk+jj-kk] + evisc_v[ijk+jj]) + visc;
                    const TF eviscs_v = TF(0.25)*(evisc_v[ijk-jj-kk] + evisc_v[ijk-jj] + evisc_v[ijk   -kk] + evisc_v[ijk   ]) + visc;
                    const TF evisct_v = evisc_v[ijk   ] + visc;
                    const TF eviscb_v = evisc_v[ijk-kk] + visc;

                    wt[ijk] +=
                             // dw/dx + du/dz
                             + ( (evisce_h*(w[ijk+ii]-w[ijk   ])*dxi + evisce_v*(u[ijk+ii]-u[ijk+ii-kk])*dzhi[k])
                               - (eviscw_h*(w[ijk   ]-w[ijk-ii])*dxi + eviscw_v*(u[ijk   ]-u[ijk+  -kk])*dzhi[k]) ) * dxi
                             // dw/dy + dv/dz
                             + ( (eviscn_h*(w[ijk+jj]-w[ijk   ])*dyi + eviscn_v*(v[ijk+jj]-v[ijk+jj-kk])*dzhi[k])
                               - (eviscs_h*(w[ijk   ]-w[ijk-jj])*dyi + eviscs_v*(v[ijk   ]-v[ijk+  -kk])*dzhi[k]) ) * dyi
                             // dw/dz + dw/dz
                             + ( rhoref[k  ] * evisct_v*(w[ijk+kk]-w[ijk   ])*dzi[k  ]
                               - rhoref[k-1] * eviscb_v*(w[ijk   ]-w[ijk-kk])*dzi[k-1] ) / rhorefh[k] * TF(2.)*dzhi[k];
                }
    }


    template <typename TF>
    void diff_c(
            TF* const restrict at,
            const TF* const restrict a,
            const TF* const restrict dzi,
            const TF* const restrict dzhi,
            const TF dxidxi, const TF dyidyi,
            const TF* const restrict evisc_h,
            const TF* const restrict evisc_v,
            const TF* const restrict fluxbot,
            const TF* const restrict fluxtop,
            const TF* const restrict rhoref,
            const TF* const restrict rhorefh,
            const TF tPr, const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii = 1;
        const TF tPr_i = TF(1)/tPr;

        // bottom boundary
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;

                const TF evisce_h = TF(0.5)*(evisc_h[ijk   ]+evisc_h[ijk+ii]) * tPr_i + visc;
                const TF eviscw_h = TF(0.5)*(evisc_h[ijk-ii]+evisc_h[ijk   ]) * tPr_i + visc;
                const TF eviscn_h = TF(0.5)*(evisc_h[ijk   ]+evisc_h[ijk+jj]) * tPr_i + visc;
                const TF eviscs_h = TF(0.5)*(evisc_h[ijk-jj]+evisc_h[ijk   ]) * tPr_i + visc;

                const TF evisct_v = TF(0.5)*(evisc_v[ijk   ]+evisc_v[ijk+kk]) * tPr_i + visc;

                at[ijk] +=
                         + ( evisce_h*(a[ijk+ii]-a[ijk   ])
                           - eviscw_h*(a[ijk   ]-a[ijk-ii]) ) * dxidxi
                         + ( eviscn_h*(a[ijk+jj]-a[ijk   ])
                           - eviscs_h*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                         + ( rhorefh[kstart+1] * evisct_v*(a[ijk+kk]-a[ijk   ])*dzhi[kstart+1]
                           + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
            }

        // top boundary
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + (kend-1)*kk;

                const TF evisce_h = TF(0.5)*(evisc_h[ijk   ]+evisc_h[ijk+ii]) * tPr_i + visc;
                const TF eviscw_h = TF(0.5)*(evisc_h[ijk-ii]+evisc_h[ijk   ]) * tPr_i + visc;
                const TF eviscn_h = TF(0.5)*(evisc_h[ijk   ]+evisc_h[ijk+jj]) * tPr_i + visc;
                const TF eviscs_h = TF(0.5)*(evisc_h[ijk-jj]+evisc_h[ijk   ]) * tPr_i + visc;

                const TF eviscb_v = TF(0.5)*(evisc_v[ijk-kk]+evisc_v[ijk   ]) * tPr_i + visc;

                at[ijk] +=
                         + ( evisce_h*(a[ijk+ii]-a[ijk   ])
                           - eviscw_h*(a[ijk   ]-a[ijk-ii]) ) * dxidxi
                         + ( eviscn_h*(a[ijk+jj]-a[ijk   ])
                           - eviscs_h*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                         + (-rhorefh[kend  ] * fluxtop[ij]
                           - rhorefh[kend-1] * eviscb_v*(a[ijk   ]-a[ijk-kk])*dzhi[kend-1] ) / rhoref[kend-1] * dzi[kend-1];
            }
        

        for (int k=kstart+1; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    const TF evisce_h = TF(0.5)*(evisc_h[ijk   ]+evisc_h[ijk+ii]) * tPr_i + visc;
                    const TF eviscw_h = TF(0.5)*(evisc_h[ijk-ii]+evisc_h[ijk   ]) * tPr_i + visc;
                    const TF eviscn_h = TF(0.5)*(evisc_h[ijk   ]+evisc_h[ijk+jj]) * tPr_i + visc;
                    const TF eviscs_h = TF(0.5)*(evisc_h[ijk-jj]+evisc_h[ijk   ]) * tPr_i + visc;

                    const TF evisct_v = TF(0.5)*(evisc_v[ijk   ]+evisc_v[ijk+kk]) * tPr_i + visc;
                    const TF eviscb_v = TF(0.5)*(evisc_v[ijk-kk]+evisc_v[ijk   ]) * tPr_i + visc;

                    at[ijk] +=
                             + ( evisce_h*(a[ijk+ii]-a[ijk   ])
                               - eviscw_h*(a[ijk   ]-a[ijk-ii]) ) * dxidxi
                             + ( eviscn_h*(a[ijk+jj]-a[ijk   ])
                               - eviscs_h*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                             + ( rhorefh[k+1] * evisct_v*(a[ijk+kk]-a[ijk   ])*dzhi[k+1]
                               - rhorefh[k  ] * eviscb_v*(a[ijk   ]-a[ijk-kk])*dzhi[k]  ) / rhoref[k] * dzi[k];
                }
    }


    template<typename TF>
    TF calc_dnmul(
            const TF* const restrict evisc_h,
            const TF* const restrict evisc_v,
            const TF* const restrict dzi,
            const TF dxidxi, const TF dyidyi,
            const TF tPr,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const TF tPrfac_i = TF(1)/std::min(TF(1.), tPr);

        // get the maximum time step for diffusion
        TF dnmul = 0;
        for (int k=kstart; k<kend; ++k)
        {
            const TF dzidzi = dzi[k] * dzi[k];

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    dnmul = std::max(dnmul, std::abs(
                                evisc_h[ijk] * tPrfac_i * (dxidxi + dyidyi) +
                                evisc_v[ijk] * tPrfac_i * dzidzi));
                }
        }

        return dnmul;
    }


    template <typename TF>
    void calc_diff_flux_c(
            TF* const restrict out,
            const TF* const restrict data,
            const TF* const restrict evisc_v,
            const TF* const restrict dzhi,
            const TF tPr, const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const TF tPr_i = TF(1)/tPr;

        #pragma omp parallel for
        for (int k=kstart+1; k<kend; ++k)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const TF eviscc = 0.5*(evisc_v[ijk-kk]+evisc_v[ijk]) * tPr_i + visc;

                    out[ijk] = - eviscc*(data[ijk] - data[ijk-kk])*dzhi[k];
                }
        }
    }


    template <typename TF>
    void calc_diff_flux_u(
            TF* const restrict out,
            const TF* const restrict data,
            const TF* const restrict w,
            const TF* const evisc_h,
            const TF* const evisc_v,
            const TF dxi, const TF* const dzhi,
            const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int ii = 1;

        #pragma omp parallel for
        for (int k=kstart+1; k<kend; ++k)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;

                    const TF eviscu_h = 0.25*(evisc_h[ijk-ii-ijcells]+evisc_h[ijk-ii]+evisc_h[ijk-ijcells]+evisc_h[ijk]) + visc;
                    const TF eviscu_v = 0.25*(evisc_v[ijk-ii-ijcells]+evisc_v[ijk-ii]+evisc_v[ijk-ijcells]+evisc_v[ijk]) + visc;

                    out[ijk] = - (eviscu_v*(data[ijk]-data[ijk-ijcells])*dzhi[k] + eviscu_h*(w[ijk]-w[ijk-ii])*dxi );
                }
        }
    }


    template <typename TF>
    void calc_diff_flux_v(
            TF* const restrict out,
            const TF* const restrict data,
            const TF* const restrict w,
            const TF* const evisc_h,
            const TF* const evisc_v,
            const TF dyi, const TF* const dzhi,
            const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        #pragma omp parallel for
        for (int k=kstart+1; k<kend; ++k)
        {
                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*icells + k*ijcells;

                        const TF eviscv_h = 0.25*(evisc_h[ijk-icells-ijcells]+evisc_h[ijk-icells]+evisc_h[ijk-ijcells]+evisc_h[ijk]) + visc;
                        const TF eviscv_v = 0.25*(evisc_v[ijk-icells-ijcells]+evisc_v[ijk-icells]+evisc_v[ijk-ijcells]+evisc_v[ijk]) + visc;

                        out[ijk] = - (eviscv_v*(data[ijk]-data[ijk-ijcells])*dzhi[k] + eviscv_h*(w[ijk]-w[ijk-icells])*dyi );
                    }
        }
    }
}
#endif
