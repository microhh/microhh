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

#ifndef DIFF_KERNELS_CUH
#define DIFF_KERNELS_CUH

#include "fast_math.h"
#include "boundary.h"
#include "constants.h"

namespace Diff_kernels_g
{
    namespace fm = Fast_math;

    template<typename TF, Surface_model surface_model> __global__
    void calc_strain2_g(
            TF* const __restrict__ strain2,
            const TF* const __restrict__ u,
            const TF* const __restrict__ v,
            const TF* const __restrict__ w,
            const TF* const __restrict__ dudz,
            const TF* const __restrict__ dvdz,
            const TF* const __restrict__ dzi,
            const TF* const __restrict__ dzhi,
            const TF dxi, const TF dyi,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
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

            if (k == kstart && surface_model == Surface_model::Enabled)
            {
                strain2[ijk] = TF(2.)*(
                    // du/dx + du/dx
                    + fm::pow2((u[ijk+ii]-u[ijk])*dxi)

                    // dv/dy + dv/dy
                    + fm::pow2((v[ijk+jj]-v[ijk])*dyi)

                    // dw/dz + dw/dz
                    + fm::pow2((w[ijk+kk]-w[ijk])*dzi[k])

                    // du/dy + dv/dx
                    + TF(0.125)*fm::pow2((u[ijk      ]-u[ijk   -jj])*dyi  + (v[ijk      ]-v[ijk-ii   ])*dxi)
                    + TF(0.125)*fm::pow2((u[ijk+ii   ]-u[ijk+ii-jj])*dyi  + (v[ijk+ii   ]-v[ijk      ])*dxi)
                    + TF(0.125)*fm::pow2((u[ijk   +jj]-u[ijk      ])*dyi  + (v[ijk   +jj]-v[ijk-ii+jj])*dxi)
                    + TF(0.125)*fm::pow2((u[ijk+ii+jj]-u[ijk+ii   ])*dyi  + (v[ijk+ii+jj]-v[ijk   +jj])*dxi)

                    // du/dz + dw/dx
                    + TF(0.5)*fm::pow2(dudz[ij])
                    + TF(0.125)*fm::pow2((w[ijk      ]-w[ijk-ii   ])*dxi)
                    + TF(0.125)*fm::pow2((w[ijk+ii   ]-w[ijk      ])*dxi)
                    + TF(0.125)*fm::pow2((w[ijk   +kk]-w[ijk-ii+kk])*dxi)
                    + TF(0.125)*fm::pow2((w[ijk+ii+kk]-w[ijk   +kk])*dxi)

                    // dv/dz + dw/dy
                    + TF(0.5)*fm::pow2(dvdz[ij])
                    + TF(0.125)*fm::pow2((w[ijk      ]-w[ijk-jj   ])*dyi)
                    + TF(0.125)*fm::pow2((w[ijk+jj   ]-w[ijk      ])*dyi)
                    + TF(0.125)*fm::pow2((w[ijk   +kk]-w[ijk-jj+kk])*dyi)
                    + TF(0.125)*fm::pow2((w[ijk+jj+kk]-w[ijk   +kk])*dyi) );

                // add a small number to avoid zero divisions
                strain2[ijk] += Constants::dsmall;
            }
            else
            {
                strain2[ijk] = TF(2.)*(
                    // du/dx + du/dx
                    + fm::pow2((u[ijk+ii]-u[ijk])*dxi)
                    // dv/dy + dv/dy
                    + fm::pow2((v[ijk+jj]-v[ijk])*dyi)
                    // dw/dz + dw/dz
                    + fm::pow2((w[ijk+kk]-w[ijk])*dzi[k])
                    // du/dy + dv/dx
                    + TF(0.125)*fm::pow2((u[ijk      ]-u[ijk   -jj])*dyi  + (v[ijk      ]-v[ijk-ii   ])*dxi)
                    + TF(0.125)*fm::pow2((u[ijk+ii   ]-u[ijk+ii-jj])*dyi  + (v[ijk+ii   ]-v[ijk      ])*dxi)
                    + TF(0.125)*fm::pow2((u[ijk   +jj]-u[ijk      ])*dyi  + (v[ijk   +jj]-v[ijk-ii+jj])*dxi)
                    + TF(0.125)*fm::pow2((u[ijk+ii+jj]-u[ijk+ii   ])*dyi  + (v[ijk+ii+jj]-v[ijk   +jj])*dxi)
                    // du/dz + dw/dx
                    + TF(0.125)*fm::pow2((u[ijk      ]-u[ijk   -kk])*dzhi[k  ] + (w[ijk      ]-w[ijk-ii   ])*dxi)
                    + TF(0.125)*fm::pow2((u[ijk+ii   ]-u[ijk+ii-kk])*dzhi[k  ] + (w[ijk+ii   ]-w[ijk      ])*dxi)
                    + TF(0.125)*fm::pow2((u[ijk   +kk]-u[ijk      ])*dzhi[k+1] + (w[ijk   +kk]-w[ijk-ii+kk])*dxi)
                    + TF(0.125)*fm::pow2((u[ijk+ii+kk]-u[ijk+ii   ])*dzhi[k+1] + (w[ijk+ii+kk]-w[ijk   +kk])*dxi)
                    // dv/dz + dw/dy
                    + TF(0.125)*fm::pow2((v[ijk      ]-v[ijk   -kk])*dzhi[k  ] + (w[ijk      ]-w[ijk-jj   ])*dyi)
                    + TF(0.125)*fm::pow2((v[ijk+jj   ]-v[ijk+jj-kk])*dzhi[k  ] + (w[ijk+jj   ]-w[ijk      ])*dyi)
                    + TF(0.125)*fm::pow2((v[ijk   +kk]-v[ijk      ])*dzhi[k+1] + (w[ijk   +kk]-w[ijk-jj+kk])*dyi)
                    + TF(0.125)*fm::pow2((v[ijk+jj+kk]-v[ijk+jj   ])*dzhi[k+1] + (w[ijk+jj+kk]-w[ijk   +kk])*dyi) );

                // add a small number to avoid zero divisions
                strain2[ijk] += Constants::dsmall;
            }
        }
    }

    template<typename TF, Surface_model surface_model> __global__
    void diff_uvw_g(
            TF* const __restrict__ ut,
            TF* const __restrict__ vt,
            TF* const __restrict__ wt,
            const TF* const __restrict__ evisc,
            const TF* const __restrict__ u,
            const TF* const __restrict__ v,
            const TF* const __restrict__ w,
            const TF* const __restrict__ fluxbotu,
            const TF* const __restrict__ fluxtopu,
            const TF* const __restrict__ fluxbotv,
            const TF* const __restrict__ fluxtopv,
            const TF* const __restrict__ dzi,
            const TF* const __restrict__ dzhi,
            const TF* const __restrict__ rhoref,
            const TF* const __restrict__ rhorefh,
            const TF dxi,
            const TF dyi,
            const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
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
            const TF evisceu = evisc[ijk   ] + visc;
            const TF eviscwu = evisc[ijk-ii] + visc;
            const TF eviscnu = TF(0.25)*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]) + visc;
            const TF eviscsu = TF(0.25)*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]) + visc;
            const TF evisctu = TF(0.25)*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]) + visc;
            const TF eviscbu = TF(0.25)*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]) + visc;

            // V
            const TF eviscev = TF(0.25)*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]) + visc;
            const TF eviscwv = TF(0.25)*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]) + visc;
            const TF eviscnv = evisc[ijk   ] + visc;
            const TF eviscsv = evisc[ijk-jj] + visc;
            const TF evisctv = TF(0.25)*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]) + visc;
            const TF eviscbv = TF(0.25)*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]) + visc;

            // W
            const TF eviscew = TF(0.25)*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+ii-kk] + evisc[ijk+ii]) + visc;
            const TF eviscww = TF(0.25)*(evisc[ijk-ii-kk] + evisc[ijk-ii] + evisc[ijk   -kk] + evisc[ijk   ]) + visc;
            const TF eviscnw = TF(0.25)*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+jj-kk] + evisc[ijk+jj]) + visc;
            const TF eviscsw = TF(0.25)*(evisc[ijk-jj-kk] + evisc[ijk-jj] + evisc[ijk   -kk] + evisc[ijk   ]) + visc;
            const TF evisctw = evisc[ijk   ] + visc;
            const TF eviscbw = evisc[ijk-kk] + visc;

            if (k == kstart && surface_model == Surface_model::Enabled)
            {
                ut[ijk] +=
                    // du/dx + du/dx
                    + (  evisceu*(u[ijk+ii]-u[ijk   ])*dxi
                       - eviscwu*(u[ijk   ]-u[ijk-ii])*dxi ) * TF(2.)*dxi
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
                    + (  eviscnv*(v[ijk+jj]-v[ijk   ])*dyi
                       - eviscsv*(v[ijk   ]-v[ijk-jj])*dyi ) * TF(2.)*dyi
                    // dv/dz + dw/dy
                    + (  rhorefh[k+1] * evisctv*((v[ijk+kk]-v[ijk   ])*dzhi[k+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
                       + rhorefh[k  ] * fluxbotv[ij] ) / rhoref[k] * dzi[k];
            }
            else if (k == kend-1 && surface_model == Surface_model::Enabled)
            {
                ut[ijk] +=
                    // du/dx + du/dx
                    + (  evisceu*(u[ijk+ii]-u[ijk   ])*dxi
                       - eviscwu*(u[ijk   ]-u[ijk-ii])*dxi ) * TF(2.)*dxi
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
                    + (  eviscnv*(v[ijk+jj]-v[ijk   ])*dyi
                       - eviscsv*(v[ijk   ]-v[ijk-jj])*dyi ) * TF(2.)*dyi
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
                    + (  rhoref[k  ] * evisctw*(w[ijk+kk]-w[ijk   ])*dzi[k  ]
                       - rhoref[k-1] * eviscbw*(w[ijk   ]-w[ijk-kk])*dzi[k-1] ) / rhorefh[k] * TF(2.) * dzhi[k];
            }
            else
            {
                ut[ijk] +=
                    // du/dx + du/dx
                    + (  evisceu*(u[ijk+ii]-u[ijk   ])*dxi
                       - eviscwu*(u[ijk   ]-u[ijk-ii])*dxi ) * TF(2.)*dxi
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
                    + (  eviscnv*(v[ijk+jj]-v[ijk   ])*dyi
                       - eviscsv*(v[ijk   ]-v[ijk-jj])*dyi ) * TF(2.)*dyi
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
                    + (  rhoref[k  ] * evisctw*(w[ijk+kk]-w[ijk   ])*dzi[k  ]
                       - rhoref[k-1] * eviscbw*(w[ijk   ]-w[ijk-kk])*dzi[k-1] ) / rhorefh[k] * TF(2.) * dzhi[k];
            }
        }
    }

    template<typename TF, Surface_model surface_model> __global__
    void diff_c_g(
            TF* const __restrict__ at,
            const TF* const __restrict__ a,
            const TF* const __restrict__ evisc,
            const TF* const __restrict__ fluxbot,
            const TF* const __restrict__ fluxtop,
            const TF* const __restrict__ dzi,
            const TF* const __restrict__ dzhi,
            const TF* const __restrict__ rhoref,
            const TF* const __restrict__ rhorefh,
            const TF dxidxi,
            const TF dyidyi,
            const TF tPri,
            const TF visc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
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

            if (k == kstart && surface_model == Surface_model::Enabled)
            {
                const TF evisce = TF(0.5)*(evisc[ijk   ]+evisc[ijk+ii])*tPri + visc;
                const TF eviscw = TF(0.5)*(evisc[ijk-ii]+evisc[ijk   ])*tPri + visc;
                const TF eviscn = TF(0.5)*(evisc[ijk   ]+evisc[ijk+jj])*tPri + visc;
                const TF eviscs = TF(0.5)*(evisc[ijk-jj]+evisc[ijk   ])*tPri + visc;
                const TF evisct = TF(0.5)*(evisc[ijk   ]+evisc[ijk+kk])*tPri + visc;

                at[ijk] +=
                    + (  evisce*(a[ijk+ii]-a[ijk   ])
                       - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi
                    + (  eviscn*(a[ijk+jj]-a[ijk   ])
                       - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                    + (  rhorefh[k+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[k+1]
                       + rhorefh[k  ] * fluxbot[ij] ) / rhoref[k] * dzi[k];
            }
            else if (k == kend-1 && surface_model == Surface_model::Enabled)
            {
                const TF evisce = TF(0.5)*(evisc[ijk   ]+evisc[ijk+ii])*tPri + visc;
                const TF eviscw = TF(0.5)*(evisc[ijk-ii]+evisc[ijk   ])*tPri + visc;
                const TF eviscn = TF(0.5)*(evisc[ijk   ]+evisc[ijk+jj])*tPri + visc;
                const TF eviscs = TF(0.5)*(evisc[ijk-jj]+evisc[ijk   ])*tPri + visc;
                const TF eviscb = TF(0.5)*(evisc[ijk-kk]+evisc[ijk   ])*tPri + visc;

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
                const TF evisce = TF(0.5)*(evisc[ijk   ]+evisc[ijk+ii])*tPri + visc;
                const TF eviscw = TF(0.5)*(evisc[ijk-ii]+evisc[ijk   ])*tPri + visc;
                const TF eviscn = TF(0.5)*(evisc[ijk   ]+evisc[ijk+jj])*tPri + visc;
                const TF eviscs = TF(0.5)*(evisc[ijk-jj]+evisc[ijk   ])*tPri + visc;
                const TF evisct = TF(0.5)*(evisc[ijk   ]+evisc[ijk+kk])*tPri + visc;
                const TF eviscb = TF(0.5)*(evisc[ijk-kk]+evisc[ijk   ])*tPri + visc;

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
    void calc_dnmul_g(
            TF* const __restrict__ dnmul,
            const TF* const __restrict__ evisc,
            const TF* const __restrict__ dzi,
            const TF tPrfac_i,
            const TF dxidxi,
            const TF dyidyi,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj,     const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            dnmul[ijk] = fabs(evisc[ijk]*tPrfac_i*(dxidxi + dyidyi + dzi[k]*dzi[k]));
        }
    }
}
#endif
