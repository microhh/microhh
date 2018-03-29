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

    template<TF> __global__ 
    void strain2_g(double* __restrict__ strain2,
                   double* __restrict__ u,  double* __restrict__ v,  double* __restrict__ w,
                   double* __restrict__ ufluxbot, double* __restrict__ vfluxbot,
                   double* __restrict__ ustar, double* __restrict__ obuk, 
                   double* __restrict__ z, double* __restrict__ dzi, double* __restrict__ dzhi, const double dxi, const double dyi, 
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
    void evisc_g(double* __restrict__ evisc, double* __restrict__ N2,
                 double* __restrict__ bfluxbot, double* __restrict__ ustar, double* __restrict__ obuk,
                 double* __restrict__ mlen,
                 const double tPri, const double z0m, const double zsl,
                 const int istart,  const int jstart, const int kstart,
                 const int iend,    const int jend,   const int kend, 
                 const int jj,      const int kk)
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
                double RitPrratio = -bfluxbot[ij]/(Constants::kappa*zsl*ustar[ij])*most::phih(zsl/obuk[ij]) / evisc[ijk] * tPri;
                RitPrratio        = fmin(RitPrratio, 1.-Constants::dsmall);
                evisc[ijk]        = mlen[k] * sqrt(evisc[ijk] * (1.-RitPrratio));
            }
            else
            {
                // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
                double RitPrratio = N2[ijk] / evisc[ijk] * tPri;
                RitPrratio        = fmin(RitPrratio, 1.-Constants::dsmall);
                evisc[ijk]        = mlen[k] * sqrt(evisc[ijk] * (1.-RitPrratio));
            }
        }
    }

    template<typename TF> __global__ 
    void evisc_neutral_g(double* __restrict__ evisc, double* __restrict__ mlen,
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
    void diff_uvw_g(double* __restrict__ ut, double* __restrict__ vt, double* __restrict__ wt, 
                    double* __restrict__ evisc,
                    double* __restrict__ u, double* __restrict__ v, double* __restrict__ w,
                    double* __restrict__ fluxbotu, double* __restrict__ fluxtopu, 
                    double* __restrict__ fluxbotv, double* __restrict__ fluxtopv, 
                    double* __restrict__ dzi, double* __restrict__ dzhi, const double dxi, const double dyi,
                    double* __restrict__ rhoref, double* __restrict__ rhorefh, 
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
            const double eviscnu = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
            const double eviscsu = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
            const double evisctu = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
            const double eviscbu = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);

            // V
            const double eviscev = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
            const double eviscwv = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
            const double evisctv = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
            const double eviscbv = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);

            // W
            const double eviscew = 0.25*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+ii-kk] + evisc[ijk+ii]);
            const double eviscww = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-ii] + evisc[ijk   -kk] + evisc[ijk   ]);
            const double eviscnw = 0.25*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+jj-kk] + evisc[ijk+jj]);
            const double eviscsw = 0.25*(evisc[ijk-jj-kk] + evisc[ijk-jj] + evisc[ijk   -kk] + evisc[ijk   ]);

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
    void diff_c_g(double* __restrict__ at, double* __restrict__ a, double* __restrict__ evisc,
                  double* __restrict__ fluxbot, double* __restrict__ fluxtop, 
                  double* __restrict__ dzi, double* __restrict__ dzhi, const double dxidxi, const double dyidyi,
                  double* __restrict__ rhoref, double* __restrict__ rhorefh, const double tPri, 
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
                const double evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])*tPri;
                const double eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])*tPri;
                const double eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])*tPri;
                const double eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])*tPri;
                const double evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])*tPri;

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
                const double evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])*tPri;
                const double eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])*tPri;
                const double eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])*tPri;
                const double eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])*tPri;
                const double eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])*tPri;

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
                const double evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])*tPri;
                const double eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])*tPri;
                const double eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])*tPri;
                const double eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])*tPri;
                const double evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])*tPri;
                const double eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])*tPri;

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
    void calc_dnmul_g(double* __restrict__ dnmul, double* __restrict__ evisc, 
                      double* __restrict__ dzi, double tPrfac, const double dxidxi, const double dyidyi,
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
void Diff_smag_2<TF>::prepare_device()
{
    Boundary_surface *boundaryptr = static_cast<Boundary_surface *>(model->boundary);

    const double n=2.;
    double mlen0;
    double *mlen = new double[grid->kcells];
    for (int k=0; k<grid->kcells; ++k) 
    {
        mlen0   = cs * pow(grid->dx*grid->dy*grid->dz[k], 1./3.);
        mlen[k] = pow(pow(1./(1./pow(mlen0, n) + 1./(pow(Constants::kappa*(grid->z[k]+boundaryptr->z0m), n))), 1./n), 2);
    }

    const int nmemsize = grid->kcells*sizeof(double);
    cuda_safe_call(cudaMalloc(&mlen_g, nmemsize));
    cuda_safe_call(cudaMemcpy(mlen_g, mlen, nmemsize, cudaMemcpyHostToDevice));

    delete[] mlen;
}
#endif

template<typename TF>
void Diff_smag_2<TF>::clear_device()
{
    cuda_safe_call(cudaFree(mlen_g));
}

#ifdef USECUDA
template<typename TF>
void Diff_smag_2<TF>::exec_viscosity()
{
    // do a cast because the base boundary class does not have the MOST related variables
    Boundary_surface *boundaryptr = static_cast<Boundary_surface *>(model->boundary);

    const int blocki = grid->ithread_block;
    const int blockj = grid->jthread_block;
    const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
    const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, grid->kcells);
    dim3 blockGPU(blocki, blockj, 1);

    const int offs = grid->memoffset;

    // Calculate total strain rate
    strain2_g<<<gridGPU, blockGPU>>>(
        &fields->sd["evisc"]->data_g[offs], 
        &fields->u->data_g[offs],  &fields->v->data_g[offs],  &fields->w->data_g[offs],
        &fields->u->datafluxbot_g[offs],  &fields->v->datafluxbot_g[offs],
        &boundaryptr->ustar_g[offs], &boundaryptr->obuk_g[offs],
        grid->z_g, grid->dzi_g, grid->dzhi_g, grid->dxi, grid->dyi,
        grid->istart,  grid->jstart, grid->kstart, 
        grid->iend,    grid->jend,   grid->kend,
        grid->icellsp, grid->ijcellsp);  
    cuda_check_error();

    // start with retrieving the stability information
    if (model->thermo->get_switch() == "0")
    {
        evisc_neutral_g<<<gridGPU, blockGPU>>>(
            &fields->sd["evisc"]->data_g[offs], mlen_g,
            grid->istart,  grid->jstart, grid->kstart, 
            grid->iend,    grid->jend,   grid->kend,
            grid->icellsp, grid->ijcellsp);  
        cuda_check_error();

        grid->boundary_cyclic_g(&fields->sd["evisc"]->data_g[offs]);
    }
    // assume buoyancy calculation is needed
    else
    {
        // store the buoyancyflux in datafluxbot of tmp1
        model->thermo->get_buoyancy_fluxbot(fields->atmp["tmp1"]);
        // store the Brunt-vaisala frequency in data of tmp1 
        model->thermo->get_thermo_field(fields->atmp["tmp1"], fields->atmp["tmp2"], "N2", false);

        // Calculate eddy viscosity
        double tPri = 1./tPr;
        evisc_g<<<gridGPU, blockGPU>>>(
            &fields->sd["evisc"]->data_g[offs], &fields->atmp["tmp1"]->data_g[offs], 
            &fields->atmp["tmp1"]->datafluxbot_g[offs], &boundaryptr->ustar_g[offs], &boundaryptr->obuk_g[offs],
            mlen_g, tPri, boundaryptr->z0m, grid->z[grid->kstart],
            grid->istart,  grid->jstart, grid->kstart, 
            grid->iend,    grid->jend,   grid->kend,
            grid->icellsp, grid->ijcellsp);  
        cuda_check_error();

        grid->boundary_cyclic_g(&fields->sd["evisc"]->data_g[offs]);
    }
}
#endif

#ifdef USECUDA
template<typename TF>
void Diff_smag_2<TF>::exec()
{
    const int blocki = grid->ithread_block;
    const int blockj = grid->jthread_block;
    const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
    const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, grid->kmax);
    dim3 blockGPU(blocki, blockj, 1);

    const int offs = grid->memoffset;
    const double dxidxi = 1./(grid->dx * grid->dx);
    const double dyidyi = 1./(grid->dy * grid->dy);
    const double tPri = 1./tPr;

    diff_uvw_g<<<gridGPU, blockGPU>>>(&fields->ut->data_g[offs], &fields->vt->data_g[offs], &fields->wt->data_g[offs],
            &fields->sd["evisc"]->data_g[offs], 
            &fields->u->data_g[offs],  &fields->v->data_g[offs],  &fields->w->data_g[offs],
            &fields->u->datafluxbot_g[offs], &fields->u->datafluxtop_g[offs],
            &fields->v->datafluxbot_g[offs], &fields->v->datafluxtop_g[offs],
            grid->dzi_g, grid->dzhi_g, grid->dxi, grid->dyi,
            fields->rhoref_g, fields->rhorefh_g,
            grid->istart,  grid->jstart, grid->kstart, 
            grid->iend,    grid->jend,   grid->kend,
            grid->icellsp, grid->ijcellsp);  
    cuda_check_error();

    for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); ++it)
        diff_c_g<<<gridGPU, blockGPU>>>(&it->second->data_g[offs], &fields->sp[it->first]->data_g[offs], &fields->sd["evisc"]->data_g[offs], 
                &fields->sp[it->first]->datafluxbot_g[offs], &fields->sp[it->first]->datafluxtop_g[offs],
                grid->dzi_g, grid->dzhi_g, dxidxi, dyidyi,
                fields->rhoref_g, fields->rhorefh_g, tPri,
                grid->istart,  grid->jstart, grid->kstart, 
                grid->iend,    grid->jend,   grid->kend,
                grid->icellsp, grid->ijcellsp);  
    cuda_check_error();
}
#endif

#ifdef USECUDA
template<typename TF>
unsigned long Diff_smag_2<TF>::get_time_limit(unsigned long idt, double dt)
{
    const int blocki = grid->ithread_block;
    const int blockj = grid->jthread_block;
    const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
    const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, grid->kmax);
    dim3 blockGPU(blocki, blockj, 1);

    const int offs = grid->memoffset;

    const double dxidxi = 1./(grid->dx * grid->dx);
    const double dyidyi = 1./(grid->dy * grid->dy);
    const double tPrfac = std::min(1., tPr);

    // Calculate dnmul in tmp1 field
    calc_dnmul_g<<<gridGPU, blockGPU>>>(
        &fields->atmp["tmp1"]->data_g[offs], &fields->sd["evisc"]->data_g[offs],
        grid->dzi_g, tPrfac, dxidxi, dyidyi,  
        grid->istart,  grid->jstart, grid->kstart, 
        grid->iend,    grid->jend,   grid->kend,
        grid->icellsp, grid->ijcellsp);  
    cuda_check_error();

    // Get maximum from tmp1 field
    double dnmul = grid->get_max_g(&fields->atmp["tmp1"]->data_g[offs], fields->atmp["tmp2"]->data_g); 
    dnmul = std::max(Constants::dsmall, dnmul);
    const unsigned long idtlim = idt * dnmax/(dnmul*dt);

    return idtlim;
}
#endif

#ifdef USECUDA
template<typename TF>
double Diff_smag_2<TF>::get_dn(double dt)
{
    const int blocki = grid->ithread_block;
    const int blockj = grid->jthread_block;
    const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
    const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, grid->kmax);
    dim3 blockGPU(blocki, blockj, 1);

    const int offs = grid->memoffset;

    const double dxidxi = 1./(grid->dx * grid->dx);
    const double dyidyi = 1./(grid->dy * grid->dy);
    const double tPrfac = std::min(1., tPr);

    // Calculate dnmul in tmp1 field
    calc_dnmul_g<<<gridGPU, blockGPU>>>(
        &fields->atmp["tmp1"]->data_g[offs], &fields->sd["evisc"]->data_g[offs],
        grid->dzi_g, tPrfac, dxidxi, dyidyi,  
        grid->istart,  grid->jstart, grid->kstart, 
        grid->iend,    grid->jend,   grid->kend,
        grid->icellsp, grid->ijcellsp);  
    cuda_check_error();

    // Get maximum from tmp1 field
    double dnmul = grid->get_max_g(&fields->atmp["tmp1"]->data_g[offs], fields->atmp["tmp2"]->data_g); 

    return dnmul*dt;
}
#endif
