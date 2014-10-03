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

#include "advec_4m.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "fd.h"
#include "tools.h"

using namespace fd::o4;

__device__ double Advec_4m_grad4(double a, double b, double c, double d, double dxi)
{
  return ( -(1./24.)*(d-a) + (27./24.)*(c-b) ) * dxi;
}

__device__ double Advec_4m_grad4x(double a, double b, double c, double d)
{
  return (-(d-a) + 27.*(c-b)); 
}

__device__ double Advec_4m_interp4(double a, double b, double c, double d) 
{
  return ci0*a + ci1*b + ci2*c + ci3*d;
}

__device__ double Advec_4m_interp2(double a, double b)
{
  return 0.5*(a + b);
}

__global__ void advec_4m_advecu(double * __restrict__ ut, double * __restrict__ u, 
                                double * __restrict__ v, double * __restrict__ w,
                                double * __restrict__ dzi4, double dxi, double dyi, 
                                int jj, int kk,
                                int istart, int jstart, int kstart,
                                int iend,   int jend,   int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart;

  int ii1 = 1;
  int ii2 = 2;
  int ii3 = 3;
  int jj1 = 1*jj;
  int jj2 = 2*jj;
  int jj3 = 3*jj;
  int kk1 = 1*kk;
  int kk2 = 2*kk;
  int kk3 = 3*kk;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;

    if(k == kstart)
    {
      ut[ijk] +=
       - Advec_4m_grad4(Advec_4m_interp4(u[ijk-ii3],     u[ijk-ii2],     u[ijk-ii1], u[ijk    ])     * Advec_4m_interp2(u[ijk-ii3], u[ijk    ]),
                         Advec_4m_interp4(u[ijk-ii2],     u[ijk-ii1],     u[ijk    ], u[ijk+ii1])     * Advec_4m_interp2(u[ijk-ii1], u[ijk    ]),
                         Advec_4m_interp4(u[ijk-ii1],     u[ijk    ],     u[ijk+ii1], u[ijk+ii2])     * Advec_4m_interp2(u[ijk    ], u[ijk+ii1]),
                         Advec_4m_interp4(u[ijk    ],     u[ijk+ii1],     u[ijk+ii2], u[ijk+ii3])     * Advec_4m_interp2(u[ijk    ], u[ijk+ii3]), dxi)

       - Advec_4m_grad4(Advec_4m_interp4(v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1]) * Advec_4m_interp2(u[ijk-jj3], u[ijk    ]),
                         Advec_4m_interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * Advec_4m_interp2(u[ijk-jj1], u[ijk    ]),
                         Advec_4m_interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * Advec_4m_interp2(u[ijk    ], u[ijk+jj1]),
                         Advec_4m_interp4(v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2]) * Advec_4m_interp2(u[ijk    ], u[ijk+jj3]), dyi)

       // boundary condition
     - Advec_4m_grad4x(-Advec_4m_interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * Advec_4m_interp2(u[ijk-kk1], u[ijk+kk2]),
                         Advec_4m_interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * Advec_4m_interp2(u[ijk-kk1], u[ijk    ]),
                         Advec_4m_interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * Advec_4m_interp2(u[ijk    ], u[ijk+kk1]),
                         Advec_4m_interp4(w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2]) * Advec_4m_interp2(u[ijk    ], u[ijk+kk3])) * dzi4[kstart];
    }
    else if(k == kend-1)
    {
      ut[ijk] +=
       - Advec_4m_grad4(Advec_4m_interp4(u[ijk-ii3],     u[ijk-ii2],     u[ijk-ii1], u[ijk    ])     * Advec_4m_interp2(u[ijk-ii3], u[ijk    ]),
                         Advec_4m_interp4(u[ijk-ii2],     u[ijk-ii1],     u[ijk    ], u[ijk+ii1])     * Advec_4m_interp2(u[ijk-ii1], u[ijk    ]),
                         Advec_4m_interp4(u[ijk-ii1],     u[ijk    ],     u[ijk+ii1], u[ijk+ii2])     * Advec_4m_interp2(u[ijk    ], u[ijk+ii1]),
                         Advec_4m_interp4(u[ijk    ],     u[ijk+ii1],     u[ijk+ii2], u[ijk+ii3])     * Advec_4m_interp2(u[ijk    ], u[ijk+ii3]), dxi)

       - Advec_4m_grad4(Advec_4m_interp4(v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1]) * Advec_4m_interp2(u[ijk-jj3], u[ijk    ]),
                         Advec_4m_interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * Advec_4m_interp2(u[ijk-jj1], u[ijk    ]),
                         Advec_4m_interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * Advec_4m_interp2(u[ijk    ], u[ijk+jj1]),
                         Advec_4m_interp4(v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2]) * Advec_4m_interp2(u[ijk    ], u[ijk+jj3]), dyi)

     - Advec_4m_grad4x( Advec_4m_interp4(w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1]) * Advec_4m_interp2(u[ijk-kk3], u[ijk    ]),
                         Advec_4m_interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * Advec_4m_interp2(u[ijk-kk1], u[ijk    ]),
                         Advec_4m_interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * Advec_4m_interp2(u[ijk    ], u[ijk+kk1]),
                        -Advec_4m_interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * Advec_4m_interp2(u[ijk-kk2], u[ijk+kk1])) * dzi4[kend-1];
    }
    else
    {
      ut[ijk] +=
       - Advec_4m_grad4(Advec_4m_interp4(u[ijk-ii3],     u[ijk-ii2],     u[ijk-ii1], u[ijk    ])     * Advec_4m_interp2(u[ijk-ii3], u[ijk    ]),
                         Advec_4m_interp4(u[ijk-ii2],     u[ijk-ii1],     u[ijk    ], u[ijk+ii1])     * Advec_4m_interp2(u[ijk-ii1], u[ijk    ]),
                         Advec_4m_interp4(u[ijk-ii1],     u[ijk    ],     u[ijk+ii1], u[ijk+ii2])     * Advec_4m_interp2(u[ijk    ], u[ijk+ii1]),
                         Advec_4m_interp4(u[ijk    ],     u[ijk+ii1],     u[ijk+ii2], u[ijk+ii3])     * Advec_4m_interp2(u[ijk    ], u[ijk+ii3]), dxi)

       - Advec_4m_grad4(Advec_4m_interp4(v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1]) * Advec_4m_interp2(u[ijk-jj3], u[ijk    ]),
                         Advec_4m_interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * Advec_4m_interp2(u[ijk-jj1], u[ijk    ]),
                         Advec_4m_interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * Advec_4m_interp2(u[ijk    ], u[ijk+jj1]),
                         Advec_4m_interp4(v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2]) * Advec_4m_interp2(u[ijk    ], u[ijk+jj3]), dyi)

      - Advec_4m_grad4x(Advec_4m_interp4(w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1]) * Advec_4m_interp2(u[ijk-kk3], u[ijk    ]),
                         Advec_4m_interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * Advec_4m_interp2(u[ijk-kk1], u[ijk    ]),
                         Advec_4m_interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * Advec_4m_interp2(u[ijk    ], u[ijk+kk1]),
                         Advec_4m_interp4(w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2]) * Advec_4m_interp2(u[ijk    ], u[ijk+kk3])) * dzi4[k];
    }
  }
}

__global__ void advec_4m_advecv(double * __restrict__ vt, double * __restrict__ u, 
                                double * __restrict__ v, double * __restrict__ w,
                                double * __restrict__ dzi4, double dxi, double dyi, 
                                int jj, int kk,
                                int istart, int jstart, int kstart,
                                int iend,   int jend,   int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart;

  int ii1 = 1;
  int ii2 = 2;
  int ii3 = 3;
  int jj1 = 1*jj;
  int jj2 = 2*jj;
  int jj3 = 3*jj;
  int kk1 = 1*kk;
  int kk2 = 2*kk;
  int kk3 = 3*kk;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    if(k == kstart)
    {
      vt[ijk] +=
       - Advec_4m_grad4(Advec_4m_interp4(u[ijk-ii1-jj2], u[ijk-ii1-jj1], u[ijk-ii1], u[ijk-ii1+jj1]) * Advec_4m_interp2(v[ijk-ii3], v[ijk    ]),
                         Advec_4m_interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * Advec_4m_interp2(v[ijk-ii1], v[ijk    ]),
                         Advec_4m_interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * Advec_4m_interp2(v[ijk    ], v[ijk+ii1]),
                         Advec_4m_interp4(u[ijk+ii2-jj2], u[ijk+ii2-jj1], u[ijk+ii2], u[ijk+ii2+jj1]) * Advec_4m_interp2(v[ijk    ], v[ijk+ii3]), dxi)

       - Advec_4m_grad4(Advec_4m_interp4(v[ijk-jj3],     v[ijk-jj2],     v[ijk-jj1], v[ijk    ])     * Advec_4m_interp2(v[ijk-jj3], v[ijk    ]),
                         Advec_4m_interp4(v[ijk-jj2],     v[ijk-jj1],     v[ijk    ], v[ijk+jj1])     * Advec_4m_interp2(v[ijk-jj1], v[ijk    ]),
                         Advec_4m_interp4(v[ijk-jj1],     v[ijk    ],     v[ijk+jj1], v[ijk+jj2])     * Advec_4m_interp2(v[ijk    ], v[ijk+jj1]),
                         Advec_4m_interp4(v[ijk    ],     v[ijk+jj1],     v[ijk+jj2], v[ijk+jj3])     * Advec_4m_interp2(v[ijk    ], v[ijk+jj3]), dyi)

     - Advec_4m_grad4x(-Advec_4m_interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * Advec_4m_interp2(v[ijk-kk1], v[ijk+kk2]),
                         Advec_4m_interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * Advec_4m_interp2(v[ijk-kk1], v[ijk    ]),
                         Advec_4m_interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * Advec_4m_interp2(v[ijk    ], v[ijk+kk1]),
                         Advec_4m_interp4(w[ijk-jj2+kk2], w[ijk-jj1+kk2], w[ijk+kk2], w[ijk+jj1+kk2]) * Advec_4m_interp2(v[ijk    ], v[ijk+kk3])) * dzi4[kstart];
    }
    else if(k == kend-1)
    {
      vt[ijk] +=
       - Advec_4m_grad4(Advec_4m_interp4(u[ijk-ii1-jj2], u[ijk-ii1-jj1], u[ijk-ii1], u[ijk-ii1+jj1]) * Advec_4m_interp2(v[ijk-ii3], v[ijk    ]),
                         Advec_4m_interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * Advec_4m_interp2(v[ijk-ii1], v[ijk    ]),
                         Advec_4m_interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * Advec_4m_interp2(v[ijk    ], v[ijk+ii1]),
                         Advec_4m_interp4(u[ijk+ii2-jj2], u[ijk+ii2-jj1], u[ijk+ii2], u[ijk+ii2+jj1]) * Advec_4m_interp2(v[ijk    ], v[ijk+ii3]), dxi)

       - Advec_4m_grad4(Advec_4m_interp4(v[ijk-jj3],     v[ijk-jj2],     v[ijk-jj1], v[ijk    ])     * Advec_4m_interp2(v[ijk-jj3], v[ijk    ]),
                         Advec_4m_interp4(v[ijk-jj2],     v[ijk-jj1],     v[ijk    ], v[ijk+jj1])     * Advec_4m_interp2(v[ijk-jj1], v[ijk    ]),
                         Advec_4m_interp4(v[ijk-jj1],     v[ijk    ],     v[ijk+jj1], v[ijk+jj2])     * Advec_4m_interp2(v[ijk    ], v[ijk+jj1]),
                         Advec_4m_interp4(v[ijk    ],     v[ijk+jj1],     v[ijk+jj2], v[ijk+jj3])     * Advec_4m_interp2(v[ijk    ], v[ijk+jj3]), dyi)

     - Advec_4m_grad4x( Advec_4m_interp4(w[ijk-jj2-kk1], w[ijk-jj1-kk1], w[ijk-kk1], w[ijk+jj1-kk1]) * Advec_4m_interp2(v[ijk-kk3], v[ijk    ]),
                         Advec_4m_interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * Advec_4m_interp2(v[ijk-kk1], v[ijk    ]),
                         Advec_4m_interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * Advec_4m_interp2(v[ijk    ], v[ijk+kk1]),
                        -Advec_4m_interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * Advec_4m_interp2(v[ijk-kk2], v[ijk+kk1])) * dzi4[kend-1];
    }
    else
    {
      vt[ijk] +=
       - Advec_4m_grad4(Advec_4m_interp4(u[ijk-ii1-jj2], u[ijk-ii1-jj1], u[ijk-ii1], u[ijk-ii1+jj1]) * Advec_4m_interp2(v[ijk-ii3], v[ijk    ]),
                         Advec_4m_interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * Advec_4m_interp2(v[ijk-ii1], v[ijk    ]),
                         Advec_4m_interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * Advec_4m_interp2(v[ijk    ], v[ijk+ii1]),
                         Advec_4m_interp4(u[ijk+ii2-jj2], u[ijk+ii2-jj1], u[ijk+ii2], u[ijk+ii2+jj1]) * Advec_4m_interp2(v[ijk    ], v[ijk+ii3]), dxi)

       - Advec_4m_grad4(Advec_4m_interp4(v[ijk-jj3],     v[ijk-jj2],     v[ijk-jj1], v[ijk    ])     * Advec_4m_interp2(v[ijk-jj3], v[ijk    ]),
                         Advec_4m_interp4(v[ijk-jj2],     v[ijk-jj1],     v[ijk    ], v[ijk+jj1])     * Advec_4m_interp2(v[ijk-jj1], v[ijk    ]),
                         Advec_4m_interp4(v[ijk-jj1],     v[ijk    ],     v[ijk+jj1], v[ijk+jj2])     * Advec_4m_interp2(v[ijk    ], v[ijk+jj1]),
                         Advec_4m_interp4(v[ijk    ],     v[ijk+jj1],     v[ijk+jj2], v[ijk+jj3])     * Advec_4m_interp2(v[ijk    ], v[ijk+jj3]), dyi)

      - Advec_4m_grad4x(Advec_4m_interp4(w[ijk-jj2-kk1], w[ijk-jj1-kk1], w[ijk-kk1], w[ijk+jj1-kk1]) * Advec_4m_interp2(v[ijk-kk3], v[ijk    ]),
                         Advec_4m_interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * Advec_4m_interp2(v[ijk-kk1], v[ijk    ]),
                         Advec_4m_interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * Advec_4m_interp2(v[ijk    ], v[ijk+kk1]),
                         Advec_4m_interp4(w[ijk-jj2+kk2], w[ijk-jj1+kk2], w[ijk+kk2], w[ijk+jj1+kk2]) * Advec_4m_interp2(v[ijk    ], v[ijk+kk3])) * dzi4[k];
    }
  }
}

__global__ void advec_4m_advecw(double * __restrict__ wt, double * __restrict__ u, 
                                double * __restrict__ v, double * __restrict__ w,
                                double * __restrict__ dzhi4, double dxi, double dyi, 
                                int jj, int kk, int istart,
                                int jstart, int kstart,
                                int iend,   int jend,   int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart + 1;

  int ii1 = 1;
  int ii2 = 2;
  int ii3 = 3;
  int jj1 = 1*jj;
  int jj2 = 2*jj;
  int jj3 = 3*jj;
  int kk1 = 1*kk;
  int kk2 = 2*kk;
  int kk3 = 3*kk;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;

    wt[ijk] +=
     - Advec_4m_grad4(Advec_4m_interp4(u[ijk-ii1-kk2], u[ijk-ii1-kk1], u[ijk-ii1], u[ijk-ii1+kk1]) * Advec_4m_interp2(w[ijk-ii3], w[ijk    ]),
                       Advec_4m_interp4(u[ijk    -kk2], u[ijk    -kk1], u[ijk    ], u[ijk    +kk1]) * Advec_4m_interp2(w[ijk-ii1], w[ijk    ]),
                       Advec_4m_interp4(u[ijk+ii1-kk2], u[ijk+ii1-kk1], u[ijk+ii1], u[ijk+ii1+kk1]) * Advec_4m_interp2(w[ijk    ], w[ijk+ii1]),
                       Advec_4m_interp4(u[ijk+ii2-kk2], u[ijk+ii2-kk1], u[ijk+ii2], u[ijk+ii2+kk1]) * Advec_4m_interp2(w[ijk    ], w[ijk+ii3]), dxi)

     - Advec_4m_grad4(Advec_4m_interp4(v[ijk-jj1-kk2], v[ijk-jj1-kk1], v[ijk-jj1], v[ijk-jj1+kk1]) * Advec_4m_interp2(w[ijk-jj3], w[ijk    ]),
                       Advec_4m_interp4(v[ijk    -kk2], v[ijk    -kk1], v[ijk    ], v[ijk    +kk1]) * Advec_4m_interp2(w[ijk-jj1], w[ijk    ]),
                       Advec_4m_interp4(v[ijk+jj1-kk2], v[ijk+jj1-kk1], v[ijk+jj1], v[ijk+jj1+kk1]) * Advec_4m_interp2(w[ijk    ], w[ijk+jj1]),
                       Advec_4m_interp4(v[ijk+jj2-kk2], v[ijk+jj2-kk1], v[ijk+jj2], v[ijk+jj2+kk1]) * Advec_4m_interp2(w[ijk    ], w[ijk+jj3]), dyi)

    - Advec_4m_grad4x(Advec_4m_interp4(w[ijk-kk3],     w[ijk-kk2],     w[ijk-kk1], w[ijk    ])     * Advec_4m_interp2(w[ijk-kk3], w[ijk    ]),
                       Advec_4m_interp4(w[ijk-kk2],     w[ijk-kk1],     w[ijk    ], w[ijk+kk1])     * Advec_4m_interp2(w[ijk-kk1], w[ijk    ]),
                       Advec_4m_interp4(w[ijk-kk1],     w[ijk    ],     w[ijk+kk1], w[ijk+kk2])     * Advec_4m_interp2(w[ijk    ], w[ijk+kk1]),
                       Advec_4m_interp4(w[ijk    ],     w[ijk+kk1],     w[ijk+kk2], w[ijk+kk3])     * Advec_4m_interp2(w[ijk    ], w[ijk+kk3])) * dzhi4[k];
  }
}

__global__ void advec_4m_advecs(double * __restrict__ st, double * __restrict__ s, 
                                double * __restrict__ u, double * __restrict__ v, double * __restrict__ w,
                                double * __restrict__ dzi4, double dxi, double dyi, 
                                int jj, int kk,
                                int istart, int jstart, int kstart,
                                int iend,   int jend,   int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  int k = blockIdx.z + kstart;

  int ii1 = 1;
  int ii2 = 2;
  int ii3 = 3;
  int jj1 = 1*jj;
  int jj2 = 2*jj;
  int jj3 = 3*jj;
  int kk1 = 1*kk;
  int kk2 = 2*kk;
  int kk3 = 3*kk;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    if(k == kstart)
    {
      st[ijk] +=
       - Advec_4m_grad4(u[ijk-ii1] * Advec_4m_interp2(s[ijk-ii3], s[ijk    ]),
                         u[ijk    ] * Advec_4m_interp2(s[ijk-ii1], s[ijk    ]),
                         u[ijk+ii1] * Advec_4m_interp2(s[ijk    ], s[ijk+ii1]),
                         u[ijk+ii2] * Advec_4m_interp2(s[ijk    ], s[ijk+ii3]), dxi)

       - Advec_4m_grad4(v[ijk-jj1] * Advec_4m_interp2(s[ijk-jj3], s[ijk    ]),
                         v[ijk    ] * Advec_4m_interp2(s[ijk-jj1], s[ijk    ]),
                         v[ijk+jj1] * Advec_4m_interp2(s[ijk    ], s[ijk+jj1]),
                         v[ijk+jj2] * Advec_4m_interp2(s[ijk    ], s[ijk+jj3]), dyi)

      - Advec_4m_grad4x(-w[ijk+kk1] * Advec_4m_interp2(s[ijk-kk1], s[ijk+kk2]),
                          w[ijk    ] * Advec_4m_interp2(s[ijk-kk1], s[ijk    ]),
                          w[ijk+kk1] * Advec_4m_interp2(s[ijk    ], s[ijk+kk1]),
                          w[ijk+kk2] * Advec_4m_interp2(s[ijk    ], s[ijk+kk3])) * dzi4[kstart];
    }
    else if(k == kend-1)
    {
      st[ijk] +=
       - Advec_4m_grad4(u[ijk-ii1] * Advec_4m_interp2(s[ijk-ii3], s[ijk    ]),
                         u[ijk    ] * Advec_4m_interp2(s[ijk-ii1], s[ijk    ]),
                         u[ijk+ii1] * Advec_4m_interp2(s[ijk    ], s[ijk+ii1]),
                         u[ijk+ii2] * Advec_4m_interp2(s[ijk    ], s[ijk+ii3]), dxi)

       - Advec_4m_grad4(v[ijk-jj1] * Advec_4m_interp2(s[ijk-jj3], s[ijk    ]),
                         v[ijk    ] * Advec_4m_interp2(s[ijk-jj1], s[ijk    ]),
                         v[ijk+jj1] * Advec_4m_interp2(s[ijk    ], s[ijk+jj1]),
                         v[ijk+jj2] * Advec_4m_interp2(s[ijk    ], s[ijk+jj3]), dyi)

     - Advec_4m_grad4x( w[ijk-kk1] * Advec_4m_interp2(s[ijk-kk3], s[ijk    ]),
                         w[ijk    ] * Advec_4m_interp2(s[ijk-kk1], s[ijk    ]),
                         w[ijk+kk1] * Advec_4m_interp2(s[ijk    ], s[ijk+kk1]),
                        -w[ijk    ] * Advec_4m_interp2(s[ijk-kk2], s[ijk+kk1])) * dzi4[kend-1];
    }
    else
    {
      st[ijk] +=
       - Advec_4m_grad4(u[ijk-ii1] * Advec_4m_interp2(s[ijk-ii3], s[ijk    ]),
                         u[ijk    ] * Advec_4m_interp2(s[ijk-ii1], s[ijk    ]),
                         u[ijk+ii1] * Advec_4m_interp2(s[ijk    ], s[ijk+ii1]),
                         u[ijk+ii2] * Advec_4m_interp2(s[ijk    ], s[ijk+ii3]), dxi)

       - Advec_4m_grad4(v[ijk-jj1] * Advec_4m_interp2(s[ijk-jj3], s[ijk    ]),
                         v[ijk    ] * Advec_4m_interp2(s[ijk-jj1], s[ijk    ]),
                         v[ijk+jj1] * Advec_4m_interp2(s[ijk    ], s[ijk+jj1]),
                         v[ijk+jj2] * Advec_4m_interp2(s[ijk    ], s[ijk+jj3]), dyi)

      - Advec_4m_grad4x(w[ijk-kk1] * Advec_4m_interp2(s[ijk-kk3], s[ijk    ]),
                         w[ijk    ] * Advec_4m_interp2(s[ijk-kk1], s[ijk    ]),
                         w[ijk+kk1] * Advec_4m_interp2(s[ijk    ], s[ijk+kk1]),
                         w[ijk+kk2] * Advec_4m_interp2(s[ijk    ], s[ijk+kk3])) * dzi4[k];
    }
  }
}

__global__ void advec_4m_calccfl(double * const __restrict__ tmp1,
                                 const double * const __restrict__ u, const double * const __restrict__ v, const double * const __restrict__ w, 
                                 const double * const __restrict__ dzi, const double dxi, const double dyi,
                                 const int jj, const int kk,
                                 const int istart, const int jstart, const int kstart,
                                 const int iend, const int jend, const int kend)
{
  const int i = blockIdx.x*blockDim.x + threadIdx.x; 
  const int j = blockIdx.y*blockDim.y + threadIdx.y; 
  const int k = blockIdx.z; 

  const int ii1 = 1;
  const int ii2 = 2;
  const int jj1 = 1*jj;
  const int jj2 = 2*jj;
  const int kk1 = 1*kk;
  const int kk2 = 2*kk;

  const int ijk = i + j*jj + k*kk;

  if(i < iend && j < jend && k < kend)
    tmp1[ijk] = std::abs(ci0*u[ijk-ii1] + ci1*u[ijk] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2])*dxi + 
                std::abs(ci0*v[ijk-jj1] + ci1*v[ijk] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2])*dyi + 
                std::abs(ci0*w[ijk-kk1] + ci1*w[ijk] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])*dzi[k];
}

#ifdef USECUDA
void Advec_4m::exec()
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxi = 1./grid->dx;
  const double dyi = 1./grid->dy;

  const int offs = grid->memoffset;

  advec_4m_advecu<<<gridGPU, blockGPU>>>(&fields->ut->data_g[offs], &fields->u->data_g[offs], &fields->v->data_g[offs], 
                                         &fields->w->data_g[offs], grid->dzi4_g, dxi, dyi,
                                         grid->icellsp, grid->ijcellsp,
                                         grid->istart, grid->jstart, grid->kstart,
                                         grid->iend,   grid->jend, grid->kend);
  cudaCheckError(); 

  advec_4m_advecv<<<gridGPU, blockGPU>>>(&fields->vt->data_g[offs], &fields->u->data_g[offs], &fields->v->data_g[offs], 
                                         &fields->w->data_g[offs], grid->dzi4_g, dxi, dyi,
                                         grid->icellsp, grid->ijcellsp,
                                         grid->istart, grid->jstart, grid->kstart,
                                         grid->iend,   grid->jend, grid->kend);
  cudaCheckError(); 

  advec_4m_advecw<<<gridGPU, blockGPU>>>(&fields->wt->data_g[offs], &fields->u->data_g[offs], &fields->v->data_g[offs], 
                                         &fields->w->data_g[offs], grid->dzhi4_g, dxi, dyi,
                                         grid->icellsp, grid->ijcellsp,
                                         grid->istart, grid->jstart, grid->kstart,
                                         grid->iend,   grid->jend, grid->kend);
  cudaCheckError(); 

  for(fieldmap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    advec_4m_advecs<<<gridGPU, blockGPU>>>(&it->second->data_g[offs], &fields->sp[it->first]->data_g[offs], 
                                           &fields->u->data_g[offs], &fields->v->data_g[offs], &fields->w->data_g[offs], 
                                           grid->dzi4_g, dxi, dyi,
                                           grid->icellsp, grid->ijcellsp,
                                           grid->istart, grid->jstart, grid->kstart,
                                           grid->iend,   grid->jend, grid->kend);
  cudaCheckError(); 
}
#endif

#ifdef USECUDA
double Advec_4m::calccfl(double * u, double * v, double * w, double * dzi, double dt)
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->icells/blocki + (grid->icells%blocki > 0);
  const int gridj  = grid->jcells/blockj + (grid->jcells%blockj > 0);
  double cfl = 0;

  dim3 gridGPU (gridi, gridj, grid->kcells);
  dim3 blockGPU(blocki, blockj, 1);

  const double dxi = 1./grid->dx;
  const double dyi = 1./grid->dy;

  const int offs = grid->memoffset;

  advec_4m_calccfl<<<gridGPU, blockGPU>>>(&fields->atmp["tmp1"]->data_g[offs],
                                          &fields->u->data_g[offs], &fields->v->data_g[offs], &fields->w->data_g[offs],
                                          grid->dzi_g, dxi, dyi,
                                          grid->icellsp, grid->ijcellsp,
                                          grid->istart, grid->jstart, grid->kstart,
                                          grid->iend,   grid->jend,   grid->kend);
  cudaCheckError(); 

  cfl = grid->getmax_g(&fields->atmp["tmp1"]->data_g[offs], fields->atmp["tmp2"]->data_g); 
  grid->getmax(&cfl); 
  cfl = cfl*dt;

  return cfl;
}
#endif
