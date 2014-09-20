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
        + pow((u[ijk+ii]-u[ijk])*dxi, 2.)
        // dv/dy + dv/dy
        + pow((v[ijk+jj]-v[ijk])*dyi, 2.)
        // dw/dz + dw/dz
        + pow((w[ijk+kk]-w[ijk])*dzi[k], 2.)
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
                                 double * __restrict__ z, double * __restrict__ dz, double dx, double dy,
                                 double cs, double tPr, double z0m,
                                 int istart, int jstart, int kstart, int iend, int jend, int kend, 
                                 int jj, int kk)

{
  const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  const int k = blockIdx.z + kstart;

  if(i < iend && j < jend && k < kend)
  {
    const int ij  = i + j*jj;
    const int ijk = i + j*jj + k*kk;
    double n = 2;

    /* BvS: TODO: pre-calculate fac offline. Now every thread within one vertical slice 
       has to calculate it.... */
    if(k == kstart)
    {
      // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
      double mlen0      = cs * pow(dx*dy*dz[k], 1./3.);
      double mlen       = pow(1./(1./pow(mlen0, n) + 1./(pow(constants::kappa*(z[k]+z0m), n))), 1./n);
      double fac        = pow(mlen, 2);
      double RitPrratio = -bfluxbot[ij]/(constants::kappa*z[k]*ustar[ij])*diff_les2s_phih(z[k]/obuk[ij]) / evisc[ijk] / tPr;
      RitPrratio        = fmin(RitPrratio, 1.-constants::dsmall);
      evisc[ijk]        = fac * sqrt(evisc[ijk]) * sqrt(1.-RitPrratio);
    }
    else
    {
      // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
      double mlen0      = cs * pow(dx*dy*dz[k], 1./3.);
      double mlen       = pow(1./(1./pow(mlen0, n) + 1./(pow(constants::kappa*(z[k]+z0m), n))), 1./n);
      double fac        = std::pow(mlen, 2.);
      double RitPrratio = N2[ijk] / evisc[ijk] / tPr;
      RitPrratio        = fmin(RitPrratio, 1.-constants::dsmall);
      evisc[ijk]        = fac * sqrt(evisc[ijk]) * sqrt(1.-RitPrratio);
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

__global__ void diff_les2s_diffc(double * __restrict__ at, double * __restrict__ a, double * __restrict__ evisc,
                                 double * __restrict__ fluxbot, double * __restrict__ fluxtop, 
                                 double * __restrict__ dzi, double * __restrict__ dzhi, double dxidxi, double dyidyi,
                                 double * __restrict__ rhoref, double * __restrict__ rhorefh, double tPr, 
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
      evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])/tPr;
      eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])/tPr;
      eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])/tPr;
      eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])/tPr;
      evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])/tPr;
      eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])/tPr;

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
      evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])/tPr;
      eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])/tPr;
      eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])/tPr;
      eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])/tPr;
      evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])/tPr;
      eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])/tPr;

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
      evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])/tPr;
      eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])/tPr;
      eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])/tPr;
      eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])/tPr;
      evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])/tPr;
      eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])/tPr;

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



/*
#ifdef USECUDA
int cdiff_les2s::execvisc()
{
  // do a cast because the base boundary class does not have the MOST related variables
  cboundary_surface *boundaryptr = static_cast<cboundary_surface *>(model->boundary);

  fields->forwardDevice();
  boundaryptr->forwardDevice();

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

  fields->backwardDevice();

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

    fields->forwardDevice();

    // Calculate eddy viscosity
    diff_les2s_evisc<<<gridGPU, blockGPU>>>(&fields->s["evisc"]->data_g[offs], &fields->s["tmp1"]->data_g[offs], 
                                            &fields->sd["tmp1"]->datafluxbot_g[offs], &boundaryptr->ustar_g[offs], &boundaryptr->obuk_g[offs],
                                            grid->z_g, grid->dz_g, grid->dx, grid->dy, cs, tPr, boundaryptr->z0m,
                                            grid->istart, grid->jstart, grid->kstart, grid->iend, grid->jend, grid->kend,
                                            grid->icellsp, grid->ijcellsp);  
    grid->boundary_cyclic_g(&fields->sd["evisc"]->data_g[offs]);

  }

  fields->backwardDevice();

  return 0;
}
#endif
*/

#ifdef USECUDA
int cdiff_les2s::exec()
{
  fields->forwardDevice();

  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const int offs = grid->memoffset;
  const double dxidxi = 1./(grid->dx * grid->dx);
  const double dyidyi = 1./(grid->dy * grid->dy);

  diff_les2s_diffu<<<gridGPU, blockGPU>>>(&fields->ut->data_g[offs], &fields->s["evisc"]->data_g[offs], 
                                          &fields->u->data_g[offs],  &fields->v->data_g[offs],  &fields->w->data_g[offs],
                                          &fields->u->datafluxbot_g[offs], &fields->u->datafluxtop_g[offs],
                                          grid->dzi_g, grid->dzhi_g, grid->dxi, grid->dyi,
                                          fields->rhoref_g, fields->rhorefh_g,
                                          grid->istart, grid->jstart, grid->kstart, grid->iend, grid->jend, grid->kend,
                                          grid->icellsp, grid->ijcellsp);  

  diff_les2s_diffv<<<gridGPU, blockGPU>>>(&fields->vt->data_g[offs], &fields->s["evisc"]->data_g[offs], 
                                          &fields->u->data_g[offs],  &fields->v->data_g[offs],  &fields->w->data_g[offs],
                                          &fields->v->datafluxbot_g[offs], &fields->v->datafluxtop_g[offs],
                                          grid->dzi_g, grid->dzhi_g, grid->dxi, grid->dyi,
                                          fields->rhoref_g, fields->rhorefh_g,
                                          grid->istart, grid->jstart, grid->kstart, grid->iend, grid->jend, grid->kend,
                                          grid->icellsp, grid->ijcellsp);  

  diff_les2s_diffw<<<gridGPU, blockGPU>>>(&fields->wt->data_g[offs], &fields->s["evisc"]->data_g[offs], 
                                          &fields->u->data_g[offs],  &fields->v->data_g[offs],  &fields->w->data_g[offs],
                                          grid->dzi_g, grid->dzhi_g, grid->dxi, grid->dyi,
                                          fields->rhoref_g, fields->rhorefh_g,
                                          grid->istart, grid->jstart, grid->kstart, grid->iend, grid->jend, grid->kend,
                                          grid->icellsp, grid->ijcellsp);  

  for(fieldmap::const_iterator it = fields->st.begin(); it!=fields->st.end(); ++it)
    diff_les2s_diffc<<<gridGPU, blockGPU>>>(&it->second->data_g[offs], &fields->s[it->first]->data_g[offs], &fields->s["evisc"]->data_g[offs], 
                                            &fields->s[it->first]->datafluxbot_g[offs], &fields->s[it->first]->datafluxtop_g[offs],
                                            grid->dzi_g, grid->dzhi_g, dxidxi, dyidyi,
                                            fields->rhoref_g, fields->rhorefh_g, tPr,
                                            grid->istart, grid->jstart, grid->kstart, grid->iend, grid->jend, grid->kend,
                                            grid->icellsp, grid->ijcellsp);  

  fields->backwardDevice();

  return 0;
}
#endif

