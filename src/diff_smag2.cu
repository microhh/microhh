/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
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
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "master.h"
#include "diff_smag2.h"
#include "boundary_surface.h"
#include "defines.h"
#include "constants.h"
#include "thermo.h"
#include "model.h"
#include "tools.h"
#include "monin_obukhov.h"

namespace
{
    namespace most = Monin_obukhov;

    template<typename TF> __global__ 
    void strain2_g(TF* __restrict__ strain2,
                   TF* __restrict__ u,  TF* __restrict__ v,  TF* __restrict__ w,
                   TF* __restrict__ ufluxbot, TF* __restrict__ vfluxbot,
                   TF* __restrict__ ustar, TF* __restrict__ obuk, 
                   TF* __restrict__ z, TF* __restrict__ dzi, TF* __restrict__ dzhi, const TF dxi, const TF dyi, 
                   const int istart, const int jstart, const int kstart, 
                   const int iend,   const int jend,   const int kend, 
                   const int jj,     const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;
        const int ii = 1;

        if (i < iend && j < jend && k < kend)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + k*kk;

            if (k == kstart)
            {
                strain2[ijk] = 2.*(
                   // du/dz
                   + 0.5*pow(-0.5*(ufluxbot[ij]+ufluxbot[ij+ii])/(Constants::kappa*z[k]*ustar[ij])*most::phim(z[k]/obuk[ij]), 2)
                   // dv/dz
                   + 0.5*pow(-0.5*(vfluxbot[ij]+vfluxbot[ij+jj])/(Constants::kappa*z[k]*ustar[ij])*most::phim(z[k]/obuk[ij]), 2) );
                // add a small number to avoid zero divisions
                strain2[ijk] += Constants::dsmall;  
            }
            else
            {
                strain2[ijk] = 2.*(
                    // du/dx + du/dx
                    + pow((u[ijk+ii]-u[ijk])*dxi, 2)
                    // dv/dy + dv/dy
                    + pow((v[ijk+jj]-v[ijk])*dyi, 2)
                    // dw/dz + dw/dz
                    + pow((w[ijk+kk]-w[ijk])*dzi[k], 2)
                    // du/dy + dv/dx
                    + 0.125*pow((u[ijk      ]-u[ijk   -jj])*dyi  + (v[ijk      ]-v[ijk-ii   ])*dxi, 2)
                    + 0.125*pow((u[ijk+ii   ]-u[ijk+ii-jj])*dyi  + (v[ijk+ii   ]-v[ijk      ])*dxi, 2)
                    + 0.125*pow((u[ijk   +jj]-u[ijk      ])*dyi  + (v[ijk   +jj]-v[ijk-ii+jj])*dxi, 2)
                    + 0.125*pow((u[ijk+ii+jj]-u[ijk+ii   ])*dyi  + (v[ijk+ii+jj]-v[ijk   +jj])*dxi, 2)
                    // du/dz + dw/dx
                    + 0.125*pow((u[ijk      ]-u[ijk   -kk])*dzhi[k  ] + (w[ijk      ]-w[ijk-ii   ])*dxi, 2)
                    + 0.125*pow((u[ijk+ii   ]-u[ijk+ii-kk])*dzhi[k  ] + (w[ijk+ii   ]-w[ijk      ])*dxi, 2)
                    + 0.125*pow((u[ijk   +kk]-u[ijk      ])*dzhi[k+1] + (w[ijk   +kk]-w[ijk-ii+kk])*dxi, 2)
                    + 0.125*pow((u[ijk+ii+kk]-u[ijk+ii   ])*dzhi[k+1] + (w[ijk+ii+kk]-w[ijk   +kk])*dxi, 2)
                    // dv/dz + dw/dy
                    + 0.125*pow((v[ijk      ]-v[ijk   -kk])*dzhi[k  ] + (w[ijk      ]-w[ijk-jj   ])*dyi, 2)
                    + 0.125*pow((v[ijk+jj   ]-v[ijk+jj-kk])*dzhi[k  ] + (w[ijk+jj   ]-w[ijk      ])*dyi, 2)
                    + 0.125*pow((v[ijk   +kk]-v[ijk      ])*dzhi[k+1] + (w[ijk   +kk]-w[ijk-jj+kk])*dyi, 2)
                    + 0.125*pow((v[ijk+jj+kk]-v[ijk+jj   ])*dzhi[k+1] + (w[ijk+jj+kk]-w[ijk   +kk])*dyi, 2) );
                // add a small number to avoid zero divisions
                strain2[ijk] += Constants::dsmall;
            }
        }
    }

    template<typename TF> __global__ 
    void evisc_g(TF* __restrict__ evisc, TF* __restrict__ N2,
                 TF* __restrict__ bfluxbot, TF* __restrict__ ustar, TF* __restrict__ obuk,
                 TF* __restrict__ mlen,
                 const TF tPri, const TF z0m, const TF zsl,
                 const int istart, const int jstart, const int kstart,
                 const int iend,   const int jend,   const int kend, 
                 const int jj,     const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + k*kk;

            if (k == kstart)
            {
                // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
                TF RitPrratio = -bfluxbot[ij]/(Constants::kappa*zsl*ustar[ij])*most::phih(zsl/obuk[ij]) / evisc[ijk] * tPri;
                RitPrratio = fmin(RitPrratio, 1.-Constants::dsmall);
                evisc[ijk] = mlen[k] * sqrt(evisc[ijk] * (1.-RitPrratio));
            }
            else
            {
                // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
                TF RitPrratio = N2[ijk] / evisc[ijk] * tPri;
                RitPrratio = fmin(RitPrratio, 1.-Constants::dsmall);
                evisc[ijk] = mlen[k] * sqrt(evisc[ijk] * (1.-RitPrratio));
            }
        }
    }

    template<typename TF> __global__ 
    void evisc_neutral_g(TF* __restrict__ evisc, TF* __restrict__ mlen,
                         const int istart, const int jstart, const int kstart, 
                         const int iend,   const int jend,   const int kend, 
                         const int jj,     const int kk)

    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            evisc[ijk]    = mlen[k] * sqrt(evisc[ijk]);
        }
    }

    template<typename TF> __global__ 
    void diff_uvw_g(TF* __restrict__ ut, TF* __restrict__ vt, TF* __restrict__ wt, 
                    TF* __restrict__ evisc,
                    TF* __restrict__ u, TF* __restrict__ v, TF* __restrict__ w,
                    TF* __restrict__ fluxbotu, TF* __restrict__ fluxtopu, 
                    TF* __restrict__ fluxbotv, TF* __restrict__ fluxtopv, 
                    TF* __restrict__ dzi, TF* __restrict__ dzhi, const TF dxi, const TF dyi,
                    TF* __restrict__ rhoref, TF* __restrict__ rhorefh, 
                    const int istart, const int jstart, const int kstart, 
                    const int iend,   const int jend,   const int kend, 
                    const int jj,     const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ii  = 1;
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + k*kk;

            // U
            const TF eviscnu = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
            const TF eviscsu = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
            const TF evisctu = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
            const TF eviscbu = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);

            // V
            const TF eviscev = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
            const TF eviscwv = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
            const TF evisctv = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
            const TF eviscbv = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);

            // W
            const TF eviscew = 0.25*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+ii-kk] + evisc[ijk+ii]);
            const TF eviscww = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-ii] + evisc[ijk   -kk] + evisc[ijk   ]);
            const TF eviscnw = 0.25*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+jj-kk] + evisc[ijk+jj]);
            const TF eviscsw = 0.25*(evisc[ijk-jj-kk] + evisc[ijk-jj] + evisc[ijk   -kk] + evisc[ijk   ]);

            if (k == kstart)
            {
                ut[ijk] +=
                    // du/dx + du/dx
                    + (  evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
                       - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
                    // du/dy + dv/dx
                    + (  eviscnu*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                       - eviscsu*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                    // du/dz + dw/dx
                    + (  rhorefh[kstart+1] * evisctu*((u[ijk+kk]-u[ijk   ])* dzhi[kstart+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
                       + rhorefh[kstart  ] * fluxbotu[ij] ) / rhoref[kstart] * dzi[kstart];

                vt[ijk] +=
                    // dv/dx + du/dy
                    + (  eviscev*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                       - eviscwv*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                    // dv/dy + dv/dy
                    + (  evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
                       - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
                    // dv/dz + dw/dy
                    + (  rhorefh[k+1] * evisctv*((v[ijk+kk]-v[ijk   ])*dzhi[k+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
                       + rhorefh[k  ] * fluxbotv[ij] ) / rhoref[k] * dzi[k];
            }
            else if (k == kend-1)
            {
                ut[ijk] +=
                    // du/dx + du/dx
                    + (  evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
                       - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
                    // du/dy + dv/dx
                    + (  eviscnu*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                       - eviscsu*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                    // du/dz + dw/dx
                    + (- rhorefh[kend  ] * fluxtopu[ij]
                       - rhorefh[kend-1] * eviscbu*((u[ijk   ]-u[ijk-kk])* dzhi[kend-1] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) / rhoref[kend-1] * dzi[kend-1];

                vt[ijk] +=
                    // dv/dx + du/dy
                    + (  eviscev*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                       - eviscwv*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                    // dv/dy + dv/dy
                    + (  evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
                       - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
                    // dv/dz + dw/dy
                    + (- rhorefh[k  ] * fluxtopv[ij]
                       - rhorefh[k-1] * eviscbv*((v[ijk   ]-v[ijk-kk])*dzhi[k-1] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) / rhoref[k-1] * dzi[k-1];

                wt[ijk] +=
                    // dw/dx + du/dz
                    + (  eviscew*((w[ijk+ii]-w[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-kk])*dzhi[k])
                       - eviscww*((w[ijk   ]-w[ijk-ii])*dxi + (u[ijk   ]-u[ijk+  -kk])*dzhi[k]) ) * dxi
                    // dw/dy + dv/dz
                    + (  eviscnw*((w[ijk+jj]-w[ijk   ])*dyi + (v[ijk+jj]-v[ijk+jj-kk])*dzhi[k])
                       - eviscsw*((w[ijk   ]-w[ijk-jj])*dyi + (v[ijk   ]-v[ijk+  -kk])*dzhi[k]) ) * dyi
                    // dw/dz + dw/dz
                    + (  rhoref[k  ] * evisc[ijk   ]*(w[ijk+kk]-w[ijk   ])*dzi[k  ]
                       - rhoref[k-1] * evisc[ijk-kk]*(w[ijk   ]-w[ijk-kk])*dzi[k-1] ) / rhorefh[k] * 2.* dzhi[k];
            }
            else
            {
                ut[ijk] +=
                    // du/dx + du/dx
                    + (  evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
                       - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
                    // du/dy + dv/dx
                    + (  eviscnu*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                       - eviscsu*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                    // du/dz + dw/dx
                    + (  rhorefh[k+1] * evisctu*((u[ijk+kk]-u[ijk   ])* dzhi[k+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
                       - rhorefh[k  ] * eviscbu*((u[ijk   ]-u[ijk-kk])* dzhi[k  ] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) / rhoref[k] * dzi[k];

                vt[ijk] +=
                    // dv/dx + du/dy
                    + (  eviscev*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                       - eviscwv*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                    // dv/dy + dv/dy
                    + (  evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
                       - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
                    // dv/dz + dw/dy
                    + (  rhorefh[k+1] * evisctv*((v[ijk+kk]-v[ijk   ])*dzhi[k+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
                       - rhorefh[k  ] * eviscbv*((v[ijk   ]-v[ijk-kk])*dzhi[k  ] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) / rhoref[k] * dzi[k];

                wt[ijk] +=
                    // dw/dx + du/dz
                    + (  eviscew*((w[ijk+ii]-w[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-kk])*dzhi[k])
                       - eviscww*((w[ijk   ]-w[ijk-ii])*dxi + (u[ijk   ]-u[ijk+  -kk])*dzhi[k]) ) * dxi
                    // dw/dy + dv/dz
                    + (  eviscnw*((w[ijk+jj]-w[ijk   ])*dyi + (v[ijk+jj]-v[ijk+jj-kk])*dzhi[k])
                       - eviscsw*((w[ijk   ]-w[ijk-jj])*dyi + (v[ijk   ]-v[ijk+  -kk])*dzhi[k]) ) * dyi
                    // dw/dz + dw/dz
                    + (  rhoref[k  ] * evisc[ijk   ]*(w[ijk+kk]-w[ijk   ])*dzi[k  ]
                       - rhoref[k-1] * evisc[ijk-kk]*(w[ijk   ]-w[ijk-kk])*dzi[k-1] ) / rhorefh[k] * 2.* dzhi[k];
            }
        }
    }

    template<typename TF> __global__ 
    void diff_c_g(TF* __restrict__ at, TF* __restrict__ a, TF* __restrict__ evisc,
                  TF* __restrict__ fluxbot, TF* __restrict__ fluxtop, 
                  TF* __restrict__ dzi, TF* __restrict__ dzhi, const TF dxidxi, const TF dyidyi,
                  TF* __restrict__ rhoref, TF* __restrict__ rhorefh, const TF tPri, 
                  const int istart, const int jstart, const int kstart, 
                  const int iend,   const int jend,   const int kend, 
                  const int jj,     const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ii  = 1;
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + k*kk;

            if (k == kstart)
            {
                const TF evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])*tPri;
                const TF eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])*tPri;
                const TF eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])*tPri;
                const TF eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])*tPri;
                const TF evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])*tPri;

                at[ijk] +=
                    + (  evisce*(a[ijk+ii]-a[ijk   ]) 
                       - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi 
                    + (  eviscn*(a[ijk+jj]-a[ijk   ]) 
                       - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                    + (  rhorefh[k+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[k+1]
                       + rhorefh[k  ] * fluxbot[ij] ) / rhoref[k] * dzi[k];
            }
            else if (k == kend-1)
            {
                const TF evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])*tPri;
                const TF eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])*tPri;
                const TF eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])*tPri;
                const TF eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])*tPri;
                const TF eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])*tPri;

                at[ijk] +=
                    + (  evisce*(a[ijk+ii]-a[ijk   ]) 
                       - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi 
                    + (  eviscn*(a[ijk+jj]-a[ijk   ]) 
                       - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                    + (- rhorefh[k  ] * fluxtop[ij]
                       - rhorefh[k-1] * eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[k-1] ) / rhoref[k-1] * dzi[k-1];
            }
            else
            {
                const TF evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])*tPri;
                const TF eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])*tPri;
                const TF eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])*tPri;
                const TF eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])*tPri;
                const TF evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])*tPri;
                const TF eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])*tPri;

                at[ijk] +=
                    + (  evisce*(a[ijk+ii]-a[ijk   ]) 
                       - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi 
                    + (  eviscn*(a[ijk+jj]-a[ijk   ]) 
                       - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                    + (  rhorefh[k+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[k+1]
                       - rhorefh[k  ] * eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[k]  ) / rhoref[k] * dzi[k];
            }
        }
    }

    template<typename TF> __global__ 
    void calc_dnmul_g(TF* __restrict__ dnmul, TF* __restrict__ evisc, 
                      TF* __restrict__ dzi, TF tPrfac, const TF dxidxi, const TF dyidyi,
                      const int istart, const int jstart, const int kstart, 
                      const int iend,   const int jend,   const int kend, 
                      const int jj,     const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            dnmul[ijk] = fabs(tPrfac*evisc[ijk]*(dxidxi + dyidyi + dzi[k]*dzi[k]));
        }
    }
}

/* Calculate the mixing length (mlen) offline, and put on GPU */
#ifdef USECUDA
template<typename TF>
void Diff_smag2<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();

    const TF n=2.;
    TF mlen0;
    TF *mlen = new TF[gd.kcells];
    for (int k=0; k<gd.kcells; ++k) 
    {
        mlen0   = cs * pow(gd.dx*gd.dy*gd.dz[k], 1./3.);
        mlen[k] = pow(pow(1./(1./pow(mlen0, n) + 1./(pow(Constants::kappa*(gd.z[k]+boundaryptr->z0m), n))), 1./n), 2);
    }

    const int nmemsize = gd.kcells*sizeof(TF);
    cuda_safe_call(cudaMalloc(&mlen_g, nmemsize));
    cuda_safe_call(cudaMemcpy(mlen_g, mlen, nmemsize, cudaMemcpyHostToDevice));

    delete[] mlen;
}
#endif

template<typename TF>
void Diff_smag2<TF>::clear_device()
{
    cuda_safe_call(cudaFree(mlen_g));
}

#ifdef USECUDA
template<typename TF>
void Diff_smag2<TF>::exec_viscosity(Boundary<TF>& boundary, Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    // Calculate total strain rate
    strain2_g<<<gridGPU, blockGPU>>>(
        fields.sd["evisc"]->fld_g, 
        fields.mp["u"]->fld_g, fields.mp["v"]->fld_g, fields.mp["w"]->fld_g,
        fields.mp["u"]->flux_bot_g, fields.mp["v"]->flux_bot_g,
        boundary.ustar_g, boundary.obuk_g,
        gd.z_g, gd.dzi_g, gd.dzhi_g, gd.dxi, gd.dyi,
        gd.istart, gd.jstart, gd.kstart, 
        gd.iend,   gd.jend,   gd.kend,
        gd.icells, gd.ijcells);  
    cuda_check_error();

    // start with retrieving the stability information
    if (thermo.get_switch() == "0")
    {
        evisc_neutral_g<<<gridGPU, blockGPU>>>(
            fields.sd["evisc"]->fld_g, mlen_g,
            gd.istart, gd.jstart, gd.kstart, 
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.ijcells);  
        cuda_check_error();

        boundary_cyclic.exec_g(fields.sd["evisc"]->fld_g);
    }
    // assume buoyancy calculation is needed
    else
    {
        // store the buoyancyflux in datafluxbot of tmp1
        thermo.get_buoyancy_fluxbot(fields.atmp["tmp1"]);
        // store the Brunt-vaisala frequency in data of tmp1 
        thermo.get_thermo_field(fields.atmp["tmp1"], fields.atmp["tmp2"], "N2", false);

        // Calculate eddy viscosity
        TF tPri = 1./tPr;
        evisc_g<<<gridGPU, blockGPU>>>(
            fields.sd["evisc"]->fld_g, &fields.atmp["tmp1"]->fld_g, 
            fields.atmp["tmp1"]->flux_bot_g, boundaryptr->ustar_g, &boundaryptr->obuk_g,
            mlen_g, tPri, boundaryptr->z0m, gd.z[gd.kstart],
            gd.istart,  gd.jstart, gd.kstart, 
            gd.iend,    gd.jend,   gd.kend,
            gd.icells,  gd.ijcells);  
        cuda_check_error();

        gd.boundary_cyclic_g(fields.sd["evisc"]->fld_g);
    }
}
#endif

#ifdef USECUDA
template<typename TF>
void Diff_smag2<TF>::exec()
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    const TF dxidxi = 1./(gd.dx * gd.dx);
    const TF dyidyi = 1./(gd.dy * gd.dy);
    const TF tPri = 1./tPr;

    diff_uvw_g<<<gridGPU, blockGPU>>>(
            fields.mt["u"]->fld_g, fields.mt["v"]->fld_g, fields.mt["w"]->fld_g,
            fields.sd["evisc"]->fld_g, 
            fields.mp["u"]->fld_g, fields.mp["v"]->fld_g, fields.mp["w"]->fld_g,
            fields.mp["u"]->flux_bot_g, fields.mp["u"]->flux_top_g,
            fields.mp["v"]->flux_bot_g, fields.mp["v"]->flux_top_g,
            gd.dzi_g, gd.dzhi_g, gd.dxi, gd.dyi,
            fields.rhoref_g, fields.rhorefh_g,
            gd.istart,  gd.jstart, gd.kstart, 
            gd.iend,    gd.jend,   gd.kend,
            gd.icells,  gd.ijcells);  
    cuda_check_error();

    for (auto it : fields.st)
        diff_c_g<<<gridGPU, blockGPU>>>(it.second->fld_g, fields.sp[it.first]->fld_g, fields.sd["evisc"]->fld_g, 
                fields.sp[it.first]->flux_bot_g, fields.sp[it.first]->flux_top_g,
                gd.dzi_g, gd.dzhi_g, dxidxi, dyidyi,
                fields.rhoref_g, fields.rhorefh_g, tPri,
                gd.istart, gd.jstart, gd.kstart, 
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells);  
    cuda_check_error();
}
#endif

#ifdef USECUDA
template<typename TF>
unsigned long Diff_smag2<TF>::get_time_limit(unsigned long idt, double dt)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    const TF dxidxi = 1./(gd.dx * gd.dx);
    const TF dyidyi = 1./(gd.dy * gd.dy);
    const TF tPrfac = std::min(1., tPr);

    auto tmp1 = fields.get_tmp_g();
    auto tmp2 = fields.get_tmp_g();

    // Calculate dnmul in tmp1 field
    calc_dnmul_g<<<gridGPU, blockGPU>>>(
        tmp1->fld_g, fields.sd["evisc"]->fld_g,
        gd.dzi_g, tPrfac, dxidxi, dyidyi,  
        gd.istart,  gd.jstart, gd.kstart, 
        gd.iend,    gd.jend,   gd.kend,
        gd.icells, gd.ijcells);  
    cuda_check_error();

    // Get maximum from tmp1 field
    double dnmul = field3d_operators.get_max_g(tmp1->fld_g, tmp2->fld_g); 
    dnmul = std::max(Constants::dsmall, dnmul);
    const unsigned long idtlim = idt * dnmax/(dnmul*dt);

    return idtlim;
}
#endif

#ifdef USECUDA
template<typename TF>
double Diff_smag2<TF>::get_dn(double dt)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    const TF dxidxi = 1./(gd.dx * gd.dx);
    const TF dyidyi = 1./(gd.dy * gd.dy);
    const TF tPrfac = std::min(1., tPr);

    // Calculate dnmul in tmp1 field
    auto dnmul_tmp = fields.get_tmp_g();

    calc_dnmul_g<<<gridGPU, blockGPU>>>(
        dnmul_tmp->fld_g, fields.sd["evisc"]->fld_g,
        gd.dzi_g, tPrfac, dxidxi, dyidyi,  
        gd.istart, gd.jstart, gd.kstart, 
        gd.iend,   gd.jend,   gd.kend,
        gd.icells, gd.ijcells);  
    cuda_check_error();

    // Get maximum from tmp1 field
    auto tmp = fields.get_tmp_g();
    double dnmul = field3d_operators.get_max_g(dnmul_tmp->fld_g, tmp->fld_g); 

    return dnmul*dt;
}
#endif

template class Diff_smag2<double>;
template class Diff_smag2<float>;
