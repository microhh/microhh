/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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
#include "diff_les2s.h"
#include "boundary_surface.h"
#include "defines.h"
#include "constants.h"
#include "thermo.h"
#include "model.h"

__device__ double diff_les2s_phim(double zeta)
{
  double phim;
  if(zeta <= 0.)
    phim = pow(1. + 3.6*pow(fabs(zeta), 2./3.), -1./2.);
  else
    phim = 1. + 5.*zeta;
  return phim;
}

__device__ double diff_les2s_phih(double zeta)
{
  double phih;
  if(zeta <= 0.)
    phih = pow(1. + 7.9*pow(fabs(zeta), 2./3.), -1./2.);
  else
    phih = 1. + 5.*zeta;
  return phih;
}

__global__ void diff_les2s_strain2(double * __restrict__ strain2,
                                   double * __restrict__ u,  double * __restrict__ v,  double * __restrict__ w,
                                   double * __restrict__ ufluxbot, double * __restrict__ vfluxbot,
                                   double * __restrict__ ustar, double * __restrict__ obuk, 
                                   double * __restrict__ z, double * __restrict__ dzi, double * __restrict__ dzhi, double dxi, double dyi, 
                                   int istart, int jstart, int kstart, int iend, int jend, int kend, 
                                   int jj, int kk)

{
  const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  const int k = blockIdx.z + kstart;
  const int ii = 1;

  if(i < iend && j < jend && k < kend)
  {
    const int ij  = i + j*jj;
    const int ijk = i + j*jj + k*kk;

    if(k == kstart)
    {
      strain2[ijk] = 2.*(
        // du/dz
        + 0.5*pow(-0.5*(ufluxbot[ij]+ufluxbot[ij+ii])/(constants::kappa*z[k]*ustar[ij])*diff_les2s_phim(z[k]/obuk[ij]), 2)
        // dv/dz
        + 0.5*pow(-0.5*(vfluxbot[ij]+vfluxbot[ij+jj])/(constants::kappa*z[k]*ustar[ij])*diff_les2s_phim(z[k]/obuk[ij]), 2) );
       // add a small number to avoid zero divisions
       strain2[ijk] += constants::dsmall;  
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
      strain2[ijk] += constants::dsmall;
    }
  }
}

__global__ void diff_les2s_evisc(double * __restrict__ evisc, double * __restrict__ N2,
                                 double * __restrict__ bfluxbot, double * __restrict__ ustar, double * __restrict__ obuk,
                                 double * __restrict__ mlen,
                                 double tPri, double z0m, double zsl,
                                 int istart, int jstart, int kstart, int iend, int jend, int kend, 
                                 int jj, int kk)

{
  //__shared__ double fac;
  
  const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  const int k = blockIdx.z + kstart;

  if(i < iend && j < jend && k < kend)
  {
    const int ij  = i + j*jj;
    const int ijk = i + j*jj + k*kk;

    if(k == kstart)
    {
      // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
      double RitPrratio = -bfluxbot[ij]/(constants::kappa*zsl*ustar[ij])*diff_les2s_phih(zsl/obuk[ij]) / evisc[ijk] * tPri;
      RitPrratio        = fmin(RitPrratio, 1.-constants::dsmall);
      evisc[ijk]        = mlen[k] * sqrt(evisc[ijk] * (1.-RitPrratio));
    }
    else
    {
      // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
      double RitPrratio = N2[ijk] / evisc[ijk] * tPri;
      RitPrratio        = fmin(RitPrratio, 1.-constants::dsmall);
      evisc[ijk]        = mlen[k] * sqrt(evisc[ijk] * (1.-RitPrratio));
    }
  }
}

__global__ void diff_les2s_diffu(double * __restrict__ ut, double * __restrict__ evisc,
                                 double * __restrict__ u, double * __restrict__ v, double * __restrict__ w,
                                 double * __restrict__ fluxbot, double * __restrict__ fluxtop, 
                                 double * __restrict__ dzi, double * __restrict__ dzhi, double dxi, double dyi,
                                 double * __restrict__ rhoref, double * __restrict__ rhorefh, 
                                 int istart, int jstart, int kstart, int iend, int jend, int kend, 
                                 int jj, int kk)

{
  const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  const int k = blockIdx.z + kstart;
  double eviscn, eviscs, eviscb, evisct;

  if(i < iend && j < jend && k < kend)
  {
    const int ii  = 1;
    const int ij  = i + j*jj;
    const int ijk = i + j*jj + k*kk;

    if(k == kstart)
    {
      eviscn = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
      eviscs = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
      evisct = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
      eviscb = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);
      ut[ijk] +=
            // du/dx + du/dx
            + (  evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
               - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
            // du/dy + dv/dx
            + (  eviscn*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
               - eviscs*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
            // du/dz + dw/dx
            + (  rhorefh[kstart+1] * evisct*((u[ijk+kk]-u[ijk   ])* dzhi[kstart+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
               + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
    }
    else if(k == kend-1)
    {
      eviscn = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
      eviscs = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
      evisct = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
      eviscb = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);
      ut[ijk] +=
            // du/dx + du/dx
            + (  evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
               - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
            // du/dy + dv/dx
            + (  eviscn*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
               - eviscs*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
            // du/dz + dw/dx
            + (- rhorefh[kend  ] * fluxtop[ij]
               - rhorefh[kend-1] * eviscb*((u[ijk   ]-u[ijk-kk])* dzhi[kend-1] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) / rhoref[kend-1] * dzi[kend-1];
    }
    else
    {
      eviscn = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
      eviscs = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
      evisct = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
      eviscb = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);
      ut[ijk] +=
            // du/dx + du/dx
            + (  evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
               - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
            // du/dy + dv/dx
            + (  eviscn*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
               - eviscs*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
            // du/dz + dw/dx
            + (  rhorefh[k+1] * evisct*((u[ijk+kk]-u[ijk   ])* dzhi[k+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
               - rhorefh[k  ] * eviscb*((u[ijk   ]-u[ijk-kk])* dzhi[k  ] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) / rhoref[k] * dzi[k];
    }
  }
}

__global__ void diff_les2s_diffv(double * __restrict__ vt, double * __restrict__ evisc,
                                 double * __restrict__ u, double * __restrict__ v, double * __restrict__ w,
                                 double * __restrict__ fluxbot, double * __restrict__ fluxtop, 
                                 double * __restrict__ dzi, double * __restrict__ dzhi, double dxi, double dyi,
                                 double * __restrict__ rhoref, double * __restrict__ rhorefh, 
                                 int istart, int jstart, int kstart, int iend, int jend, int kend, 
                                 int jj, int kk)

{
  const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  const int k = blockIdx.z + kstart;
  double evisce,eviscw,eviscb,evisct;

  if(i < iend && j < jend && k < kend)
  {
    const int ii  = 1;
    const int ij  = i + j*jj;
    const int ijk = i + j*jj + k*kk;

    if(k == kstart)
    {
      evisce = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
      eviscw = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
      evisct = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
      eviscb = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);
      vt[ijk] +=
            // dv/dx + du/dy
            + (  evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
               - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
            // dv/dy + dv/dy
            + (  evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
               - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
            // dv/dz + dw/dy
            + (  rhorefh[k+1] * evisct*((v[ijk+kk]-v[ijk   ])*dzhi[k+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
               + rhorefh[k  ] * fluxbot[ij] ) / rhoref[k] * dzi[k];
    }
    else if(k == kend-1)
    {
      evisce = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
      eviscw = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
      evisct = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
      eviscb = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);
      vt[ijk] +=
            // dv/dx + du/dy
            + (  evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
               - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
            // dv/dy + dv/dy
            + (  evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
               - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
            // dv/dz + dw/dy
            + (- rhorefh[k  ] * fluxtop[ij]
               - rhorefh[k-1] * eviscb*((v[ijk   ]-v[ijk-kk])*dzhi[k-1] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) / rhoref[k-1] * dzi[k-1];
    }
    else
    {
      evisce = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
      eviscw = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
      evisct = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
      eviscb = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);
      vt[ijk] +=
            // dv/dx + du/dy
            + (  evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
               - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
            // dv/dy + dv/dy
            + (  evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
               - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
            // dv/dz + dw/dy
            + (  rhorefh[k+1] * evisct*((v[ijk+kk]-v[ijk   ])*dzhi[k+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
               - rhorefh[k  ] * eviscb*((v[ijk   ]-v[ijk-kk])*dzhi[k  ] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) / rhoref[k] * dzi[k];
    }
  }
}

__global__ void diff_les2s_diffw(double * __restrict__ wt, double * __restrict__ evisc,
                                 double * __restrict__ u, double * __restrict__ v, double * __restrict__ w,
                                 double * __restrict__ dzi, double * __restrict__ dzhi, double dxi, double dyi,
                                 double * __restrict__ rhoref, double * __restrict__ rhorefh, 
                                 int istart, int jstart, int kstart, int iend, int jend, int kend, 
                                 int jj, int kk)

{
  const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  const int k = blockIdx.z + kstart+1;
  double evisce, eviscw, eviscn, eviscs;

  if(i < iend && j < jend && k < kend)
  {
    const int ii  = 1;
    const int ijk = i + j*jj + k*kk;

    evisce = 0.25*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+ii-kk] + evisc[ijk+ii]);
    eviscw = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-ii] + evisc[ijk   -kk] + evisc[ijk   ]);
    eviscn = 0.25*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+jj-kk] + evisc[ijk+jj]);
    eviscs = 0.25*(evisc[ijk-jj-kk] + evisc[ijk-jj] + evisc[ijk   -kk] + evisc[ijk   ]);
    wt[ijk] +=
          // dw/dx + du/dz
          + (  evisce*((w[ijk+ii]-w[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-kk])*dzhi[k])
             - eviscw*((w[ijk   ]-w[ijk-ii])*dxi + (u[ijk   ]-u[ijk+  -kk])*dzhi[k]) ) * dxi
          // dw/dy + dv/dz
          + (  eviscn*((w[ijk+jj]-w[ijk   ])*dyi + (v[ijk+jj]-v[ijk+jj-kk])*dzhi[k])
             - eviscs*((w[ijk   ]-w[ijk-jj])*dyi + (v[ijk   ]-v[ijk+  -kk])*dzhi[k]) ) * dyi
          // dw/dz + dw/dz
          + (  rhoref[k  ] * evisc[ijk   ]*(w[ijk+kk]-w[ijk   ])*dzi[k  ]
             - rhoref[k-1] * evisc[ijk-kk]*(w[ijk   ]-w[ijk-kk])*dzi[k-1] ) / rhorefh[k] * 2.* dzhi[k];
  }
}

__global__ void diff_les2s_diffuvw(double * __restrict__ ut, double * __restrict__ vt, double * __restrict__ wt, 
                                   double * __restrict__ evisc,
                                   double * __restrict__ u, double * __restrict__ v, double * __restrict__ w,
                                   double * __restrict__ fluxbotu, double * __restrict__ fluxtopu, 
                                   double * __restrict__ fluxbotv, double * __restrict__ fluxtopv, 
                                   double * __restrict__ dzi, double * __restrict__ dzhi, double dxi, double dyi,
                                   double * __restrict__ rhoref, double * __restrict__ rhorefh, 
                                   int istart, int jstart, int kstart, int iend, int jend, int kend, 
                                   int jj, int kk)

{
  __shared__ double s[12]; // Contains rhoref, rhorefh, dzi, dzhi at k-1, k, k+1
  double * rhorefs  = &s[0];
  double * rhorefhs = &s[3];
  double * dzis     = &s[6];
  double * dzhis    = &s[9];

  const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  const int k = blockIdx.z + kstart;
  double eviscnu, eviscsu, eviscbu, evisctu;
  double eviscev, eviscwv, eviscbv, evisctv;
  double eviscew, eviscww, eviscnw, eviscsw;

  const int kms = 0;
  const int ks  = 1;
  const int kps = 2;  

  if(threadIdx.x == 0 and threadIdx.y == 0)
  {
    rhorefs[kms]  = rhoref[k-1];
    rhorefs[ks]   = rhoref[k];
    rhorefs[kps]  = rhoref[k+1];
    rhorefhs[kms] = rhorefh[k-1];
    rhorefhs[ks]  = rhorefh[k];
    rhorefhs[kps] = rhorefh[k+1];
    dzis[kms]     = dzi[k-1];
    dzis[ks]      = dzi[k];
    dzis[kps]     = dzi[k+1];
    dzhis[kms]    = dzhi[k-1];
    dzhis[ks]     = dzhi[k];
    dzhis[kps]    = dzhi[k+1];
  }
  __syncthreads();

  if(i < iend && j < jend && k < kend)
  {
    const int ii  = 1;
    const int ij  = i + j*jj;
    const int ijk = i + j*jj + k*kk;

    // U
    eviscnu = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
    eviscsu = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
    evisctu = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
    eviscbu = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);

    // V
    eviscev = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
    eviscwv = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
    evisctv = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
    eviscbv = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);

    // W
    eviscew = 0.25*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+ii-kk] + evisc[ijk+ii]);
    eviscww = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-ii] + evisc[ijk   -kk] + evisc[ijk   ]);
    eviscnw = 0.25*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+jj-kk] + evisc[ijk+jj]);
    eviscsw = 0.25*(evisc[ijk-jj-kk] + evisc[ijk-jj] + evisc[ijk   -kk] + evisc[ijk   ]);


    if(k == kstart)
    {
      ut[ijk] +=
            // du/dx + du/dx
            + (  evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
               - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
            // du/dy + dv/dx
            + (  eviscnu*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
               - eviscsu*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
            // du/dz + dw/dx
            + (  rhorefhs[kps] * evisctu*((u[ijk+kk]-u[ijk   ])* dzhis[kps] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
               + rhorefhs[ks ] * fluxbotu[ij] ) / rhorefs[ks] * dzis[ks];

      vt[ijk] +=
            // dv/dx + du/dy
            + (  eviscev*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
               - eviscwv*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
            // dv/dy + dv/dy
            + (  evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
               - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
            // dv/dz + dw/dy
            + (  rhorefhs[kps] * evisctv*((v[ijk+kk]-v[ijk   ])*dzhis[kps] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
               + rhorefhs[ks ] * fluxbotv[ij] ) / rhorefs[ks] * dzis[ks];
    }
    else if(k == kend-1)
    {
      ut[ijk] +=
            // du/dx + du/dx
            + (  evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
               - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
            // du/dy + dv/dx
            + (  eviscnu*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
               - eviscsu*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
            // du/dz + dw/dx
            + (- rhorefhs[kps] * fluxtopu[ij]
               - rhorefhs[ks ] * eviscbu*((u[ijk   ]-u[ijk-kk])* dzhis[ks] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) / rhorefs[ks] * dzis[ks];

      vt[ijk] +=
            // dv/dx + du/dy
            + (  eviscev*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
               - eviscwv*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
            // dv/dy + dv/dy
            + (  evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
               - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
            // dv/dz + dw/dy
            + (- rhorefhs[kps] * fluxtopv[ij]
               - rhorefhs[ks] * eviscbv*((v[ijk   ]-v[ijk-kk])* dzhis[ks] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) / rhorefs[ks] * dzis[ks];

      wt[ijk] +=
            // dw/dx + du/dz
            + (  eviscew*((w[ijk+ii]-w[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-kk])*dzhi[k])
               - eviscww*((w[ijk   ]-w[ijk-ii])*dxi + (u[ijk   ]-u[ijk+  -kk])*dzhi[k]) ) * dxi
            // dw/dy + dv/dz
            + (  eviscnw*((w[ijk+jj]-w[ijk   ])*dyi + (v[ijk+jj]-v[ijk+jj-kk])*dzhi[k])
               - eviscsw*((w[ijk   ]-w[ijk-jj])*dyi + (v[ijk   ]-v[ijk+  -kk])*dzhi[k]) ) * dyi
            // dw/dz + dw/dz
            + (  rhorefs[ks ] * evisc[ijk   ]*(w[ijk+kk]-w[ijk   ])*dzis[ks ]
               - rhorefs[kms] * evisc[ijk-kk]*(w[ijk   ]-w[ijk-kk])*dzis[kms] ) / rhorefhs[ks] * 2.* dzhis[ks];
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
            + (  rhorefhs[kps] * evisctu*((u[ijk+kk]-u[ijk   ])* dzhis[kps] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
               - rhorefhs[ks ] * eviscbu*((u[ijk   ]-u[ijk-kk])* dzhis[ks ] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) / rhorefs[ks] * dzis[ks];

      vt[ijk] +=
            // dv/dx + du/dy
            + (  eviscev*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
               - eviscwv*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
            // dv/dy + dv/dy
            + (  evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
               - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
            // dv/dz + dw/dy
            + (  rhorefhs[kps] * evisctv*((v[ijk+kk]-v[ijk   ])*dzhis[kps] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
               - rhorefhs[ks ] * eviscbv*((v[ijk   ]-v[ijk-kk])*dzhis[ks ] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) / rhorefs[ks] * dzis[ks];

      wt[ijk] +=
            // dw/dx + du/dz
            + (  eviscew*((w[ijk+ii]-w[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-kk])*dzhi[k])
               - eviscww*((w[ijk   ]-w[ijk-ii])*dxi + (u[ijk   ]-u[ijk+  -kk])*dzhi[k]) ) * dxi
            // dw/dy + dv/dz
            + (  eviscnw*((w[ijk+jj]-w[ijk   ])*dyi + (v[ijk+jj]-v[ijk+jj-kk])*dzhi[k])
               - eviscsw*((w[ijk   ]-w[ijk-jj])*dyi + (v[ijk   ]-v[ijk+  -kk])*dzhi[k]) ) * dyi
            // dw/dz + dw/dz
            + (  rhorefs[ks ] * evisc[ijk   ]*(w[ijk+kk]-w[ijk   ])*dzis[ks ]
               - rhorefs[kms] * evisc[ijk-kk]*(w[ijk   ]-w[ijk-kk])*dzis[kms] ) / rhorefhs[ks] * 2.* dzhis[ks];
    }
  }
}


__global__ void diff_les2s_diffc(double * __restrict__ at, double * __restrict__ a, double * __restrict__ evisc,
                                 double * __restrict__ fluxbot, double * __restrict__ fluxtop, 
                                 double * __restrict__ dzi, double * __restrict__ dzhi, double dxidxi, double dyidyi,
                                 double * __restrict__ rhoref, double * __restrict__ rhorefh, double tPri, 
                                 int istart, int jstart, int kstart, int iend, int jend, int kend, 
                                 int jj, int kk)

{
  const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  const int k = blockIdx.z + kstart;
  double evisce,eviscw,eviscn,eviscs,evisct,eviscb;

  if(i < iend && j < jend && k < kend)
  {
    const int ii  = 1;
    const int ij  = i + j*jj;
    const int ijk = i + j*jj + k*kk;

    if(k == kstart)
    {
      evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])*tPri;
      eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])*tPri;
      eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])*tPri;
      eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])*tPri;
      evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])*tPri;
      eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])*tPri;

      at[ijk] +=
            + (  evisce*(a[ijk+ii]-a[ijk   ]) 
               - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi 
            + (  eviscn*(a[ijk+jj]-a[ijk   ]) 
               - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
            + (  rhorefh[k+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[k+1]
               + rhorefh[k  ] * fluxbot[ij] ) / rhoref[k] * dzi[k];
    }
    else if(k == kend-1)
    {
      evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])*tPri;
      eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])*tPri;
      eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])*tPri;
      eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])*tPri;
      evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])*tPri;
      eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])*tPri;

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
      evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])*tPri;
      eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])*tPri;
      eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])*tPri;
      eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])*tPri;
      evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])*tPri;
      eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])*tPri;

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

__global__ void diff_les2s_calcdnmul(double * __restrict__ dnmul, double * __restrict__ evisc, 
                                     double * __restrict__ dzi, double tPrfac, double dxidxi, double dyidyi,
                                     int istart, int jstart, int kstart, int iend, int jend, int kend, int jj, int kk)

{
  const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  const int k = blockIdx.z + kstart;

  if(i < iend && j < jend && k < kend)
  {
    const int ijk = i + j*jj + k*kk;
    dnmul[ijk] = fabs(tPrfac*evisc[ijk]*(dxidxi + dyidyi + dzi[k]*dzi[k]));
  }
}

/* Calculate the mixing length (mlen) offline, and put on GPU */
#ifdef USECUDA
int cdiff_les2s::prepareDevice()
{
  cboundary_surface *boundaryptr = static_cast<cboundary_surface *>(model->boundary);

  const double n=2.;
  double mlen0;
  double *mlen = new double[grid->kcells];
  for(int k=0; k<grid->kcells; ++k) 
  {
    mlen0   = cs * pow(grid->dx*grid->dy*grid->dz[k], 1./3.);
    mlen[k] = pow(pow(1./(1./pow(mlen0, n) + 1./(pow(constants::kappa*(grid->z[k]+boundaryptr->z0m), n))), 1./n), 2);
  }

  const int nmemsize = grid->kcells*sizeof(double);
  cudaMalloc(&mlen_g, nmemsize);
  cudaMemcpy(mlen_g, mlen, nmemsize, cudaMemcpyHostToDevice);

  delete[] mlen;

  return 0;
}
#endif

#ifdef USECUDA
int cdiff_les2s::execvisc()
{
  // do a cast because the base boundary class does not have the MOST related variables
  cboundary_surface *boundaryptr = static_cast<cboundary_surface *>(model->boundary);

  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kcells);
  dim3 blockGPU(blocki, blockj, 1);

  const int offs = grid->memoffset;

  // Calculate total strain rate
  diff_les2s_strain2<<<gridGPU, blockGPU>>>(&fields->s["evisc"]->data_g[offs], 
                                            &fields->u->data_g[offs],  &fields->v->data_g[offs],  &fields->w->data_g[offs],
                                            &fields->u->datafluxbot_g[offs],  &fields->v->datafluxbot_g[offs],
                                            &boundaryptr->ustar_g[offs], &boundaryptr->obuk_g[offs],
                                            grid->z_g, grid->dzi_g, grid->dzhi_g, grid->dxi, grid->dyi,
                                            grid->istart, grid->jstart, grid->kstart, grid->iend, grid->jend, grid->kend,
                                            grid->icellsp, grid->ijcellsp);  

  // start with retrieving the stability information
  if(model->thermo->getsw() == "0")
  {
    master->printMessage("diff_les2s without thermo not yet supported on GPU\n");
    //evisc_neutral(fields->s["evisc"]->data,
    //              fields->u->data, fields->v->data, fields->w->data,
    //              fields->u->datafluxbot, fields->v->datafluxbot,
    //              grid->z, grid->dz, boundaryptr->z0m);
  }
  // assume buoyancy calculation is needed
  else
  {
    // store the buoyancyflux in datafluxbot of tmp1
    model->thermo->getbuoyancyfluxbot(fields->sd["tmp1"]);
    // store the Brunt-vaisala frequency in data of tmp1 
    model->thermo->getthermofield(fields->sd["tmp1"], fields->sd["tmp2"], "N2");

    // Calculate eddy viscosity
    double tPri = 1./tPr;
    diff_les2s_evisc<<<gridGPU, blockGPU>>>(&fields->s["evisc"]->data_g[offs], &fields->s["tmp1"]->data_g[offs], 
                                            &fields->sd["tmp1"]->datafluxbot_g[offs], &boundaryptr->ustar_g[offs], &boundaryptr->obuk_g[offs],
                                            mlen_g, tPri, boundaryptr->z0m, grid->z[grid->kstart],
                                            grid->istart, grid->jstart, grid->kstart, grid->iend, grid->jend, grid->kend,
                                            grid->icellsp, grid->ijcellsp);  
    grid->boundary_cyclic_g(&fields->sd["evisc"]->data_g[offs]);

  }

  return 0;
}
#endif

#ifdef USECUDA
int cdiff_les2s::exec()
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const int offs = grid->memoffset;
  const double dxidxi = 1./(grid->dx * grid->dx);
  const double dyidyi = 1./(grid->dy * grid->dy);
  const double tPri = 1./tPr;

  //diff_les2s_diffu<<<gridGPU, blockGPU>>>(&fields->ut->data_g[offs], &fields->s["evisc"]->data_g[offs], 
  //                                        &fields->u->data_g[offs],  &fields->v->data_g[offs],  &fields->w->data_g[offs],
  //                                        &fields->u->datafluxbot_g[offs], &fields->u->datafluxtop_g[offs],
  //                                        grid->dzi_g, grid->dzhi_g, grid->dxi, grid->dyi,
  //                                        fields->rhoref_g, fields->rhorefh_g,
  //                                        grid->istart, grid->jstart, grid->kstart, grid->iend, grid->jend, grid->kend,
  //                                        grid->icellsp, grid->ijcellsp);  

  //diff_les2s_diffv<<<gridGPU, blockGPU>>>(&fields->vt->data_g[offs], &fields->s["evisc"]->data_g[offs], 
  //                                        &fields->u->data_g[offs],  &fields->v->data_g[offs],  &fields->w->data_g[offs],
  //                                        &fields->v->datafluxbot_g[offs], &fields->v->datafluxtop_g[offs],
  //                                        grid->dzi_g, grid->dzhi_g, grid->dxi, grid->dyi,
  //                                        fields->rhoref_g, fields->rhorefh_g,
  //                                        grid->istart, grid->jstart, grid->kstart, grid->iend, grid->jend, grid->kend,
  //                                        grid->icellsp, grid->ijcellsp);  

  //diff_les2s_diffw<<<gridGPU, blockGPU>>>(&fields->wt->data_g[offs], &fields->s["evisc"]->data_g[offs], 
  //                                        &fields->u->data_g[offs],  &fields->v->data_g[offs],  &fields->w->data_g[offs],
  //                                        grid->dzi_g, grid->dzhi_g, grid->dxi, grid->dyi,
  //                                        fields->rhoref_g, fields->rhorefh_g,
  //                                        grid->istart, grid->jstart, grid->kstart, grid->iend, grid->jend, grid->kend,
  //                                        grid->icellsp, grid->ijcellsp);  

  diff_les2s_diffuvw<<<gridGPU, blockGPU>>>(&fields->ut->data_g[offs], &fields->vt->data_g[offs], &fields->wt->data_g[offs],
                                            &fields->s["evisc"]->data_g[offs], 
                                            &fields->u->data_g[offs],  &fields->v->data_g[offs],  &fields->w->data_g[offs],
                                            &fields->u->datafluxbot_g[offs], &fields->u->datafluxtop_g[offs],
                                            &fields->v->datafluxbot_g[offs], &fields->v->datafluxtop_g[offs],
                                            grid->dzi_g, grid->dzhi_g, grid->dxi, grid->dyi,
                                            fields->rhoref_g, fields->rhorefh_g,
                                            grid->istart, grid->jstart, grid->kstart, grid->iend, grid->jend, grid->kend,
                                            grid->icellsp, grid->ijcellsp);  

  for(fieldmap::const_iterator it = fields->st.begin(); it!=fields->st.end(); ++it)
    diff_les2s_diffc<<<gridGPU, blockGPU>>>(&it->second->data_g[offs], &fields->s[it->first]->data_g[offs], &fields->s["evisc"]->data_g[offs], 
                                            &fields->s[it->first]->datafluxbot_g[offs], &fields->s[it->first]->datafluxtop_g[offs],
                                            grid->dzi_g, grid->dzhi_g, dxidxi, dyidyi,
                                            fields->rhoref_g, fields->rhorefh_g, tPri,
                                            grid->istart, grid->jstart, grid->kstart, grid->iend, grid->jend, grid->kend,
                                            grid->icellsp, grid->ijcellsp);  

  return 0;
}
#endif

#ifdef USECUDA
unsigned long cdiff_les2s::gettimelim(unsigned long idt, double dt)
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  double dnmul;
  unsigned long idtlim;
  const int offs = grid->memoffset;
  const double dxidxi = 1./(grid->dx * grid->dx);
  const double dyidyi = 1./(grid->dy * grid->dy);
  const double tPrfac = std::min(1., tPr);

  // Calculate dnmul in tmp1 field
  diff_les2s_calcdnmul<<<gridGPU, blockGPU>>>(&fields->s["tmp1"]->data_g[offs], &fields->s["evisc"]->data_g[offs],
                                              grid->dzi_g, tPrfac, dxidxi, dyidyi,  
                                              grid->istart, grid->jstart, grid->kstart, grid->iend, grid->jend, grid->kend,
                                              grid->icellsp, grid->ijcellsp);  

  // Get maximum from tmp1 field
  dnmul = grid->getmax_g(&fields->a["tmp1"]->data_g[offs], fields->a["tmp2"]->data_g); 
  dnmul = std::max(constants::dsmall, dnmul);
  idtlim = idt * dnmax/(dnmul*dt);

  return idtlim;
}
#endif

#ifdef USECUDA
double cdiff_les2s::getdn(double dt)
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  double dnmul;
  const int offs = grid->memoffset;
  const double dxidxi = 1./(grid->dx * grid->dx);
  const double dyidyi = 1./(grid->dy * grid->dy);
  const double tPrfac = std::min(1., tPr);

  // Calculate dnmul in tmp1 field
  diff_les2s_calcdnmul<<<gridGPU, blockGPU>>>(&fields->s["tmp1"]->data_g[offs], &fields->s["evisc"]->data_g[offs],
                                              grid->dzi_g, tPrfac, dxidxi, dyidyi,  
                                              grid->istart, grid->jstart, grid->kstart, grid->iend, grid->jend, grid->kend,
                                              grid->icellsp, grid->ijcellsp);  

  // Get maximum from tmp1 field
  dnmul = grid->getmax_g(&fields->a["tmp1"]->data_g[offs], fields->a["tmp2"]->data_g); 

  return dnmul*dt;
}
#endif

