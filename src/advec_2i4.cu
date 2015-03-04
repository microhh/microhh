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
#include "advec_2i4.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "fd.h"
#include "tools.h"
#include "constants.h"

using namespace fd::o2;
using namespace fd::o4;

namespace Advec2i4_g
{
  __global__ void advecu(double * __restrict__ ut, double * __restrict__ u, 
                         double * __restrict__ v, double * __restrict__ w,
                         double * __restrict__ rhoref, double * __restrict__ rhorefh, 
                         double * __restrict__ dzi, double dxi, double dyi, 
                         int jj, int kk,
                         int istart, int jstart, int kstart,
                         int iend,   int jend,   int kend)
  {
    int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
    int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
    int k = blockIdx.z + kstart;
  
    int ii1 = 1;
    int ii2 = 2;
    int jj1 = 1*jj;
    int jj2 = 2*jj;
    int kk1 = 1*kk;
    int kk2 = 2*kk;
  
    if(i < iend && j < jend && k < kend)
    {
      int ijk = i + j*jj + k*kk;

      if(k == kstart)
      {
        ut[ijk] += 
              - (  interp2(u[ijk        ], u[ijk+ii1]) * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                 - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

              - (  interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                 - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 

              - (  rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp2(u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k];
      }
      else if(k == kstart+1)
      {
        ut[ijk] += 
                - (  interp2(u[ijk        ], u[ijk+ii1]) * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                   - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

                - (  interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                   - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 

                - (  rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp4(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2])
                   - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp2(u[ijk-kk1], u[ijk    ]) ) / rhoref[k] * dzi[k];
      }
      else if(k == kend-2)
      {
        ut[ijk] += 
                - (  interp2(u[ijk        ], u[ijk+ii1]) * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                   - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

                - (  interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                   - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 

                - (  rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp2(u[ijk    ], u[ijk+kk1])
                   - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp4(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k];
      }
      else if(k == kend-1)
      {
        ut[ijk] += 
                - (  interp2(u[ijk        ], u[ijk+ii1]) * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                   - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

                - (  interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                   - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 

                - ( -rhorefh[k] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp2(u[ijk-kk1], u[ijk    ]) ) / rhoref[k] * dzi[k];
      }
      else
      {
        ut[ijk] += 
                - (  interp2(u[ijk        ], u[ijk+ii1]) * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                   - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

                - (  interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                   - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 

                - (  rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp4(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2])
                   - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp4(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k];
      }
    }
  }
  
  __global__ void advecv(double * __restrict__ vt, double * __restrict__ u, 
                         double * __restrict__ v, double * __restrict__ w,
                         double * __restrict__ rhoref, double * __restrict__ rhorefh, 
                         double * __restrict__ dzi, double dxi, double dyi, 
                         int jj, int kk,
                         int istart, int jstart, int kstart,
                         int iend,   int jend,   int kend)
  {
    int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
    int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
    int k = blockIdx.z + kstart;
  
    int ii1 = 1;
    int ii2 = 2;
    int jj1 = 1*jj;
    int jj2 = 2*jj;
    int kk1 = 1*kk;
    int kk2 = 2*kk;
  
    if(i < iend && j < jend && k < kend)
    {
      int ijk = i + j*jj + k*kk;

      if(k == kstart)
      {
        vt[ijk] += 
                - (  interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                   - interp2(u[ijk    -jj1], u[ijk    ]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

                - (  interp2(v[ijk        ], v[ijk+jj1]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                   - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

                - (  rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp2(v[ijk    ], v[ijk+kk1]) ) / rhoref[k] * dzi[k];
      }
      else if(k == kstart+1)
      {
        vt[ijk] += 
                - (  interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                   - interp2(u[ijk    -jj1], u[ijk    ]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

                - (  interp2(v[ijk        ], v[ijk+jj1]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                   - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

                - (  rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp4(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])
                   - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp2(v[ijk-kk1], v[ijk    ]) ) / rhoref[k] * dzi[k];
      }
      else if(k == kend-2)
      {
        vt[ijk] += 
                - (  interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                   - interp2(u[ijk    -jj1], u[ijk    ]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

                - (  interp2(v[ijk        ], v[ijk+jj1]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                   - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

                - (  rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp2(v[ijk    ], v[ijk+kk1])
                   - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp4(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) / rhoref[k] * dzi[k];
      }
      else if(k == kend-1)
      {
        vt[ijk] +=
                - (  interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                   - interp2(u[ijk    -jj1], u[ijk    ]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

                - (  interp2(v[ijk        ], v[ijk+jj1]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                   - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

                - (- rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp2(v[ijk-kk1], v[ijk    ]) ) / rhoref[k] * dzi[k];
      }
      else
      {
        vt[ijk] += 
                - (  interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                   - interp2(u[ijk    -jj1], u[ijk    ]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

                - (  interp2(v[ijk        ], v[ijk+jj1]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                   - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

                - (  rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp4(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])
                   - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp4(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) / rhoref[k] * dzi[k];
      }
    }
  }
  
  __global__ void advecw(double * __restrict__ wt, double * __restrict__ u, 
                         double * __restrict__ v, double * __restrict__ w,
                         double * __restrict__ rhoref, double * __restrict__ rhorefh, 
                         double * __restrict__ dzhi, double dxi, double dyi, 
                         int jj, int kk, int istart,
                         int jstart, int kstart,
                         int iend,   int jend,   int kend)
  {
    int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
    int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
    int k = blockIdx.z + kstart + 1;
  
    int ii1 = 1;
    int ii2 = 2;
    int jj1 = 1*jj;
    int jj2 = 2*jj;
    int kk1 = 1*kk;
    int kk2 = 2*kk;
  
    if(i < iend && j < jend && k < kend)
    {
      int ijk = i + j*jj + k*kk;

      if(k == kstart+1)
      {
        wt[ijk] += 
                - (  interp2(u[ijk+ii1-kk1], u[ijk+ii1]) * interp4(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
                   - interp2(u[ijk    -kk1], u[ijk    ]) * interp4(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi

                - (  interp2(v[ijk+jj1-kk1], v[ijk+jj1]) * interp4(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
                   - interp2(v[ijk    -kk1], v[ijk    ]) * interp4(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi

                - (  rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp4(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2])
                   - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp2(w[ijk-kk1], w[ijk    ]) ) / rhorefh[k] * dzhi[k];
      }
      else if(k == kend-1)
      {
        wt[ijk] += 
                - (  interp2(u[ijk+ii1-kk1], u[ijk+ii1]) * interp4(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
                   - interp2(u[ijk    -kk1], u[ijk    ]) * interp4(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi

                - (  interp2(v[ijk+jj1-kk1], v[ijk+jj1]) * interp4(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
                   - interp2(v[ijk    -kk1], v[ijk    ]) * interp4(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi

                - (  rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp2(w[ijk    ], w[ijk+kk1])
                   - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp4(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) / rhorefh[k] * dzhi[k];
      }
      else
      {
        wt[ijk] +=
                - (  interp2(u[ijk+ii1-kk1], u[ijk+ii1]) * interp4(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
                   - interp2(u[ijk    -kk1], u[ijk    ]) * interp4(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi

                - (  interp2(v[ijk+jj1-kk1], v[ijk+jj1]) * interp4(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
                   - interp2(v[ijk    -kk1], v[ijk    ]) * interp4(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi

                - (  rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp4(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2])
                   - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp4(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) / rhorefh[k] * dzhi[k];
      }
    }
  }
  
  __global__ void advecs(double * __restrict__ st, double * __restrict__ s, 
                         double * __restrict__ u, double * __restrict__ v, double * __restrict__ w,
                         double * __restrict__ rhoref, double * __restrict__ rhorefh, 
                         double * __restrict__ dzi, double dxi, double dyi, 
                         int jj, int kk,
                         int istart, int jstart, int kstart,
                         int iend,   int jend,   int kend)
  {
    int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
    int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
    int k = blockIdx.z + kstart;
  
    int ii1 = 1;
    int ii2 = 2;
    int jj1 = 1*jj;
    int jj2 = 2*jj;
    int kk1 = 1*kk;
    int kk2 = 2*kk;
  
    if(i < iend && j < jend && k < kend)
    {
      int ijk = i + j*jj + k*kk;

      if(k == kstart)
      {
        st[ijk] += 
              - (  u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                 - u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

              - (  v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                 - v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 

              - (  rhorefh[k+1] * w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
      }
      else if(k == kstart+1)
      {
        st[ijk] += 
              - (  u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                 - u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

              - (  v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                 - v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 

              - (  rhorefh[k+1] * w[ijk+kk1] * interp4(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                 - rhorefh[k  ] * w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]) ) / rhoref[k] * dzi[k];
      }
      else if(k == kend-2)
      {
        st[ijk] += 
              - (  u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                 - u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

              - (  v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                 - v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 

              - (  rhorefh[k+1] * w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1])
                 - rhorefh[k  ] * w[ijk    ] * interp4(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
      }
      else if(k == kend-1)
      {
        st[ijk] += 
              - (  u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                 - u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

              - (  v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                 - v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 

              - (- rhorefh[k  ] * w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]) ) / rhoref[k] * dzi[k];
      }
      else
      {
        st[ijk] += 
              - (  u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                 - u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

              - (  v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                 - v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 

              - (  rhorefh[k+1] * w[ijk+kk1] * interp4(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                 - rhorefh[k  ] * w[ijk    ] * interp4(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
      }
    }
  }
  
  __global__ void calc_cfl(double * const __restrict__ tmp1,
                           const double * const __restrict__ u, const double * const __restrict__ v, const double * const __restrict__ w, 
                           const double * const __restrict__ dzi, const double dxi, const double dyi,
                           const int jj, const int kk,
                           const int istart, const int jstart, const int kstart,
                           const int iend, const int jend, const int kend)
  {
    const int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
    const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
    const int k = blockIdx.z + kstart; 
  
    const int ii1 = 1;
    const int ii2 = 2;
    const int jj1 = 1*jj;
    const int jj2 = 2*jj;
    const int kk1 = 1*kk;
    const int kk2 = 2*kk;
  
    const int ijk = i + j*jj + k*kk;
  
    if(i < iend && j < jend && k < kend)
    {
      if(k == kstart || k == kend-1)
        tmp1[ijk] = fabs(interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]))*dxi 
                  + fabs(interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]))*dyi 
                  + fabs(interp2(w[ijk    ], w[ijk+kk1]))*dzi[k];
      else
        tmp1[ijk] = fabs(interp4(u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2]))*dxi 
                  + fabs(interp4(v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2]))*dyi 
                  + fabs(interp4(w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2]))*dzi[k];
    }
  }
}

#ifdef USECUDA
unsigned long Advec_2i4::get_time_limit(unsigned long idt, double dt)
{
  unsigned long idtlim;
  double cfl;

  // Calculate cfl and prevent zero divisons.
  cfl = get_cfl(dt);
  cfl = std::max(cflmin, cfl);

  idtlim = idt * cflmax / cfl;
  return idtlim;
}
#endif

#ifdef USECUDA
double Advec_2i4::get_cfl(const double dt)
{
  const int blocki = grid->iThreadBlock;
  const int blockj = grid->jThreadBlock;
  const int gridi  = grid->icells/blocki + (grid->icells%blocki > 0);
  const int gridj  = grid->jcells/blockj + (grid->jcells%blockj > 0);

  double cfl = 0;

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxi = 1./grid->dx;
  const double dyi = 1./grid->dy;

  const int offs = grid->memoffset;

  Advec2i4_g::calc_cfl<<<gridGPU, blockGPU>>>(&fields->atmp["tmp1"]->data_g[offs],
                                              &fields->u->data_g[offs], &fields->v->data_g[offs], &fields->w->data_g[offs],
                                              grid->dzi_g, dxi, dyi,
                                              grid->icellsp, grid->ijcellsp,
                                              grid->istart, grid->jstart, grid->kstart,
                                              grid->iend,   grid->jend,   grid->kend);
  cudaCheckError(); 

  cfl = grid->getMax_g(&fields->atmp["tmp1"]->data_g[offs], fields->atmp["tmp2"]->data_g); 
  grid->getMax(&cfl); 
  cfl = cfl*dt;

  return cfl;
}
#endif

#ifdef USECUDA
void Advec_2i4::exec()
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

  Advec2i4_g::advecu<<<gridGPU, blockGPU>>>(&fields->ut->data_g[offs], &fields->u->data_g[offs], 
                                            &fields->v->data_g[offs], &fields->w->data_g[offs],
                                            fields->rhoref_g, fields->rhorefh_g,
                                            grid->dzi_g, dxi, dyi,
                                            grid->icellsp, grid->ijcellsp,
                                            grid->istart, grid->jstart, grid->kstart,
                                            grid->iend,   grid->jend, grid->kend);
  cudaCheckError(); 

  Advec2i4_g::advecv<<<gridGPU, blockGPU>>>(&fields->vt->data_g[offs], &fields->u->data_g[offs],
                                            &fields->v->data_g[offs],  &fields->w->data_g[offs],
                                            fields->rhoref_g, fields->rhorefh_g,
                                            grid->dzi_g, dxi, dyi,
                                            grid->icellsp, grid->ijcellsp,
                                            grid->istart, grid->jstart, grid->kstart,
                                            grid->iend,   grid->jend, grid->kend);
  cudaCheckError(); 

  Advec2i4_g::advecw<<<gridGPU, blockGPU>>>(&fields->wt->data_g[offs], &fields->u->data_g[offs], 
                                            &fields->v->data_g[offs], &fields->w->data_g[offs],
                                            fields->rhoref_g, fields->rhorefh_g,
                                            grid->dzhi_g, dxi, dyi,
                                            grid->icellsp, grid->ijcellsp,
                                            grid->istart, grid->jstart, grid->kstart,
                                            grid->iend,   grid->jend, grid->kend);
  cudaCheckError(); 

  for(FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    Advec2i4_g::advecs<<<gridGPU, blockGPU>>>(&it->second->data_g[offs], &fields->sp[it->first]->data_g[offs], 
                                              &fields->u->data_g[offs], &fields->v->data_g[offs], &fields->w->data_g[offs], 
                                              fields->rhoref_g, fields->rhorefh_g,
                                              grid->dzi_g, dxi, dyi,
                                              grid->icellsp, grid->ijcellsp,
                                              grid->istart, grid->jstart, grid->kstart,
                                              grid->iend,   grid->jend, grid->kend);
  cudaCheckError(); 
}
#endif
