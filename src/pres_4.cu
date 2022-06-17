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

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fftw3.h>
#include <cufft.h>
#include <iostream>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "finite_difference.h"
#include "pres.h"
#include "pres_4.h"
#include "defines.h"
#include "tools.h"
#include "constants.h"
#include "field3d_operators.h"
#include "stats.h"

namespace
{
    using namespace Finite_difference::O4;

    template<typename TF> __global__
    void ghost_cells_wt_g(TF* const __restrict__ wt,
                const int jj, const int kk,
                const int istart, const int jstart, const int kstart,
                const int iend,   const int jend,   const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            int ijk = i + j*jj + kstart*kk;
            wt[ijk-kk] = -wt[ijk+kk];

            ijk = i + j*jj + kend*kk;
            wt[ijk+kk] = -wt[ijk-kk];
        }
    }

    template<typename TF> __global__
    void pres_in_g(TF* const __restrict__ p,
                   const TF* const __restrict__ u , const TF* const __restrict__ v , const TF* const __restrict__ w ,
                   const TF* const __restrict__ ut, const TF* const __restrict__ vt, const TF* const __restrict__ wt,
                   const TF* const __restrict__ dzi4,
                   const TF dxi, const TF dyi, const TF dti,
                   const int jj,   const int kk,
                   const int jjp,  const int kkp,
                   const int imax, const int jmax, const int kmax,
                   const int igc,  const int jgc,  const int kgc)
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

        if (i < imax && j < jmax && k < kmax)
        {
            const int ijkp = i + j*jjp + k*kkp;
            const int ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;

            p[ijkp] = (cg0<TF>*(ut[ijk-ii1] + u[ijk-ii1]*dti) + cg1<TF>*(ut[ijk] + u[ijk]*dti) + cg2<TF>*(ut[ijk+ii1] + u[ijk+ii1]*dti) + cg3<TF>*(ut[ijk+ii2] + u[ijk+ii2]*dti)) * dxi
                    + (cg0<TF>*(vt[ijk-jj1] + v[ijk-jj1]*dti) + cg1<TF>*(vt[ijk] + v[ijk]*dti) + cg2<TF>*(vt[ijk+jj1] + v[ijk+jj1]*dti) + cg3<TF>*(vt[ijk+jj2] + v[ijk+jj2]*dti)) * dyi
                    + (cg0<TF>*(wt[ijk-kk1] + w[ijk-kk1]*dti) + cg1<TF>*(wt[ijk] + w[ijk]*dti) + cg2<TF>*(wt[ijk+kk1] + w[ijk+kk1]*dti) + cg3<TF>*(wt[ijk+kk2] + w[ijk+kk2]*dti)) * dzi4[k+kgc];
        }
    }

    template<typename TF> __global__
    void solve_in_g(const TF* const __restrict__ p,
                    const TF* const __restrict__ m1, const TF* const __restrict__ m2, const TF* const __restrict__ m3, const TF* const __restrict__ m4,
                    const TF* const __restrict__ m5, const TF* const __restrict__ m6, const TF* const __restrict__ m7,
                    TF* const __restrict__ m1temp, TF* const __restrict__ m2temp, TF* __restrict__ const m3temp, TF* const __restrict__ m4temp,
                    TF* const __restrict__ m5temp, TF* const __restrict__ m6temp, TF* __restrict__ const m7temp, TF* const __restrict__ ptemp,
                    const TF* const __restrict__ bmati, const TF* const __restrict__ bmatj,
                    const int mpicoordx, const int mpicoordy,
                    const int iblock, const int jblock,
                    const int kmax,
                    const int n, const int jslice)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int jj = iblock;
        const int kk = iblock*jblock;

        const int kki1 = 1*iblock*jslice;
        const int kki2 = 2*iblock*jslice;
        const int kki3 = 3*iblock*jslice;

        int ik,ijk,iindex,jindex;

        if (i < iblock && j < jslice)
        {
            // Swap the mpicoords, because domain is turned 90 degrees to avoid two mpi transposes.
            iindex = mpicoordy*iblock + i;
            jindex = mpicoordx*jblock + n*jslice + j;

            // Set a zero gradient bc at the bottom.
            ik = i + j*jj;
            m1temp[ik] = TF( 0.);
            m2temp[ik] = TF( 0.);
            m3temp[ik] = TF( 0.);
            m4temp[ik] = TF( 1.);
            m5temp[ik] = TF( 0.);
            m6temp[ik] = TF( 0.);
            m7temp[ik] = TF(-1.);
            ptemp [ik] = TF( 0.);

            m1temp[ik+kki1] = TF( 0.);
            m2temp[ik+kki1] = TF( 0.);
            m3temp[ik+kki1] = TF( 0.);
            m4temp[ik+kki1] = TF( 1.);
            m5temp[ik+kki1] = TF(-1.);
            m6temp[ik+kki1] = TF( 0.);
            m7temp[ik+kki1] = TF( 0.);
            ptemp [ik+kki1] = TF( 0.);

            for (int k=0; k<kmax; ++k)
            {
                // Swap the mpicoords, because domain is turned 90 degrees to avoid two mpi transposes.
                ijk = i + (j + n*jslice)*jj + k*kk;
                ik  = i + j*jj + k*kki1;
                m1temp[ik+kki2] = m1[k];
                m2temp[ik+kki2] = m2[k];
                m3temp[ik+kki2] = m3[k];
                m4temp[ik+kki2] = m4[k] + bmati[iindex] + bmatj[jindex];
                m5temp[ik+kki2] = m5[k];
                m6temp[ik+kki2] = m6[k];
                m7temp[ik+kki2] = m7[k];
                ptemp [ik+kki2] = p[ijk];
            }

            // Set the top boundary.
            ik = i + j*jj + kmax*kki1;
            if (iindex == 0 && jindex == 0)
            {
                m1temp[ik+kki2] = TF(   0.);
                m2temp[ik+kki2] = TF(-1/3.);
                m3temp[ik+kki2] = TF(   2.);
                m4temp[ik+kki2] = TF(   1.);

                m1temp[ik+kki3] = TF(  -2.);
                m2temp[ik+kki3] = TF(   9.);
                m3temp[ik+kki3] = TF(   0.);
                m4temp[ik+kki3] = TF(   1.);
            }

            // Set dp/dz at top to zero.
            else
            {
                m1temp[ik+kki2] = TF( 0.);
                m2temp[ik+kki2] = TF( 0.);
                m3temp[ik+kki2] = TF(-1.);
                m4temp[ik+kki2] = TF( 1.);

                m1temp[ik+kki3] = TF(-1.);
                m2temp[ik+kki3] = TF( 0.);
                m3temp[ik+kki3] = TF( 0.);
                m4temp[ik+kki3] = TF( 1.);
            }

            // Set the top boundary.
            m5temp[ik+kki2] = TF(0.);
            m6temp[ik+kki2] = TF(0.);
            m7temp[ik+kki2] = TF(0.);
            ptemp [ik+kki2] = TF(0.);

            m5temp[ik+kki3] = TF(0.);
            m6temp[ik+kki3] = TF(0.);
            m7temp[ik+kki3] = TF(0.);
            ptemp [ik+kki3] = TF(0.);
        }
    }

    template<typename TF> __global__
    void solve_put_back_g(TF* const __restrict__ p,
                          const TF* const __restrict__ ptemp,
                          const int iblock, const int jblock,
                          const int kmax,
                          const int n, const int jslice)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int jj = iblock;
        const int kk = iblock*jblock;

        const int kki1 = 1*iblock*jslice;
        const int kki2 = 2*iblock*jslice;

        if (i < iblock && j < jslice)
        {
            // Put back the solution.
            for (int k=0; k<kmax; ++k)
            {
                const int ik  = i + j*jj + k*kki1;
                const int ijk = i + (j + n*jslice)*jj + k*kk;
                p[ijk] = ptemp[ik+kki2];
            }
        }
    }

    template<typename TF> __global__
    void hdma_g(TF* const __restrict__ m1, TF* const __restrict__ m2, TF* const __restrict__ m3, TF* const __restrict__ m4,
                TF* const __restrict__ m5, TF* const __restrict__ m6, TF* const __restrict__ m7, TF* const __restrict__ p,
                const int iblock, const int kmax, const int jslice)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        const int jj = iblock;

        const int kk1 = 1*iblock*jslice;
        const int kk2 = 2*iblock*jslice;
        const int kk3 = 3*iblock*jslice;

        int k,ik;

        if (i < iblock && j < jslice)
        {
            // Use LU factorization.
            k = 0;
            ik = i + j*jj;
            m1[ik] = TF(1.);
            m2[ik] = TF(1.);
            m3[ik] = TF(1.)            / m4[ik];
            m4[ik] = TF(1.);
            m5[ik] = m5[ik]*m3[ik];
            m6[ik] = m6[ik]*m3[ik];
            m7[ik] = m7[ik]*m3[ik];

            k = 1;
            ik = i + j*jj + k*kk1;
            m1[ik] = TF(1.);
            m2[ik] = TF(1.);
            m3[ik] = m3[ik]                     / m4[ik-kk1];
            m4[ik] = m4[ik] - m3[ik]*m5[ik-kk1];
            m5[ik] = m5[ik] - m3[ik]*m6[ik-kk1];
            m6[ik] = m6[ik] - m3[ik]*m7[ik-kk1];

            k = 2;
            ik = i + j*jj + k*kk1;
            m1[ik] = TF(1.);
            m2[ik] =   m2[ik]                                           / m4[ik-kk2];
            m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] ) / m4[ik-kk1];
            m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2];
            m5[ik] =   m5[ik] - m3[ik]*m6[ik-kk1] - m2[ik]*m7[ik-kk2];
            m6[ik] =   m6[ik] - m3[ik]*m7[ik-kk1];

            for (k=3; k<kmax+2; ++k)
            {
                ik = i + j*jj + k*kk1;
                m1[ik] = ( m1[ik]                                                            ) / m4[ik-kk3];
                m2[ik] = ( m2[ik]                                         - m1[ik]*m5[ik-kk3]) / m4[ik-kk2];
                m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] - m1[ik]*m6[ik-kk3]) / m4[ik-kk1];
                m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2] - m1[ik]*m7[ik-kk3];
                m5[ik] =   m5[ik] - m3[ik]*m6[ik-kk1] - m2[ik]*m7[ik-kk2];
                m6[ik] =   m6[ik] - m3[ik]*m7[ik-kk1];
            }

            k = kmax+1;
            ik = i + j*jj + k*kk1;
            m7[ik] = TF(1.);

            k = kmax+2;
            ik = i + j*jj + k*kk1;
            m1[ik] = ( m1[ik]                                                            ) / m4[ik-kk3];
            m2[ik] = ( m2[ik]                                         - m1[ik]*m5[ik-kk3]) / m4[ik-kk2];
            m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] - m1[ik]*m6[ik-kk3]) / m4[ik-kk1];
            m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2] - m1[ik]*m7[ik-kk3];
            m5[ik] =   m5[ik] - m3[ik]*m6[ik-kk1] - m2[ik]*m7[ik-kk2];
            m6[ik] = TF(1.);
            m7[ik] = TF(1.);

            k = kmax+3;
            ik = i + j*jj + k*kk1;
            m1[ik] = ( m1[ik]                                                            ) / m4[ik-kk3];
            m2[ik] = ( m2[ik]                                         - m1[ik]*m5[ik-kk3]) / m4[ik-kk2];
            m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] - m1[ik]*m6[ik-kk3]) / m4[ik-kk1];
            m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2] - m1[ik]*m7[ik-kk3];
            m5[ik] = TF(1.);
            m6[ik] = TF(1.);
            m7[ik] = TF(1.);

            // Do the backward substitution.
            // First, solve Ly = p, forward.
            ik = i + j*jj;
            p[ik    ] =             p[ik    ]*m3[ik    ];
            p[ik+kk1] = p[ik+kk1] - p[ik    ]*m3[ik+kk1];
            p[ik+kk2] = p[ik+kk2] - p[ik+kk1]*m3[ik+kk2] - p[ik]*m2[ik+kk2];

            for (k=3; k<kmax+4; ++k)
            {
                ik = i + j*jj + k*kk1;
                p[ik] = p[ik] - p[ik-kk1]*m3[ik] - p[ik-kk2]*m2[ik] - p[ik-kk3]*m1[ik];
            }

            // Second, solve Ux=y, backward.
            k = kmax+3;
            ik = i + j*jj + k*kk1;
            p[ik    ] =   p[ik    ]                                             / m4[ik    ];
            p[ik-kk1] = ( p[ik-kk1] - p[ik    ]*m5[ik-kk1] )                    / m4[ik-kk1];
            p[ik-kk2] = ( p[ik-kk2] - p[ik-kk1]*m5[ik-kk2] - p[ik]*m6[ik-kk2] ) / m4[ik-kk2];

            for (k=kmax; k>=0; --k)
            {
                ik = i + j*jj + k*kk1;
                p[ik] = ( p[ik] - p[ik+kk1]*m5[ik] - p[ik+kk2]*m6[ik] - p[ik+kk3]*m7[ik] ) / m4[ik];
            }
        }
    }

    template<typename TF> __global__
    void solve_out_g(TF* __restrict__ p, TF* __restrict__ work3d,
                     const int jj, const int kk,
                     const int jjp, const int kkp,
                     const int istart, const int jstart, const int kstart,
                     const int imax,   const int jmax,   const int kmax)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;
        const int k = blockIdx.z;

        const int kkp1 = 1*kkp;
        const int kkp2 = 2*kkp;

        if (i < imax && j < jmax && k < kmax)
        {
            const int ijk  = i + j*jj + k*kk;
            const int ijkp = i+istart + (j+jstart)*jjp + (k+kstart)*kkp;

            p[ijkp] = work3d[ijk];

            // set the BC
            if (k == 0)
            {
                p[ijkp-kkp1] = p[ijkp     ];
                p[ijkp-kkp2] = p[ijkp+kkp1];
            }
            else if (k == kmax-1)
            {
                p[ijkp+kkp1] = p[ijkp     ];
                p[ijkp+kkp2] = p[ijkp-kkp1];
            }
        }
    }

    template<typename TF> __global__
    void pres_out_g(TF* const __restrict__ ut, TF* const __restrict__ vt, TF* const __restrict__ wt,
                    const TF* const __restrict__ p,
                    const TF* const __restrict__ dzhi4,
                    const TF dxi, const TF dyi,
                    const int jj,     const int kk,
                    const int istart, const int jstart, const int kstart,
                    const int iend,   const int jend,   const int kend)
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

        if (i < iend && j < jend && k == kstart)
        {
            const int ijk = i + j*jj + k*kk;
            ut[ijk] -= (cg0<TF>*p[ijk-ii2] + cg1<TF>*p[ijk-ii1] + cg2<TF>*p[ijk] + cg3<TF>*p[ijk+ii1]) * dxi;
            vt[ijk] -= (cg0<TF>*p[ijk-jj2] + cg1<TF>*p[ijk-jj1] + cg2<TF>*p[ijk] + cg3<TF>*p[ijk+jj1]) * dyi;
        }
        else if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj1 + k*kk1;
            ut[ijk] -= (cg0<TF>*p[ijk-ii2] + cg1<TF>*p[ijk-ii1] + cg2<TF>*p[ijk] + cg3<TF>*p[ijk+ii1]) * dxi;
            vt[ijk] -= (cg0<TF>*p[ijk-jj2] + cg1<TF>*p[ijk-jj1] + cg2<TF>*p[ijk] + cg3<TF>*p[ijk+jj1]) * dyi;
            wt[ijk] -= (cg0<TF>*p[ijk-kk2] + cg1<TF>*p[ijk-kk1] + cg2<TF>*p[ijk] + cg3<TF>*p[ijk+kk1]) * dzhi4[k];
        }
    }

    template<typename TF> __global__
    void calc_divergence_g(TF* __restrict__ div,
                           TF* __restrict__ u, TF* __restrict__ v, TF* __restrict__ w,
                           TF* __restrict__ dzi4,
                           TF dxi, TF dyi,
                           int jj,     int kk,
                           int istart, int jstart, int kstart,
                           int iend,   int jend,   int kend)
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

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            div[ijk] = (cg0<TF>*u[ijk-ii1] + cg1<TF>*u[ijk] + cg2<TF>*u[ijk+ii1] + cg3<TF>*u[ijk+ii2]) * dxi
                     + (cg0<TF>*v[ijk-jj1] + cg1<TF>*v[ijk] + cg2<TF>*v[ijk+jj1] + cg3<TF>*v[ijk+jj2]) * dyi
                     + (cg0<TF>*w[ijk-kk1] + cg1<TF>*w[ijk] + cg2<TF>*w[ijk+kk1] + cg3<TF>*w[ijk+kk2]) * dzi4[k];
        }
    }
} // End namespace.


#ifdef USECUDA
template<typename TF>
void Pres_4<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();

    const int kmemsize = gd.kmax*sizeof(TF);
    const int imemsize = gd.itot*sizeof(TF);
    const int jmemsize = gd.jtot*sizeof(TF);

    cuda_safe_call(cudaMalloc((void**)&bmati_g, imemsize));
    cuda_safe_call(cudaMalloc((void**)&bmatj_g, jmemsize));

    cuda_safe_call(cudaMalloc((void**)&m1_g, kmemsize));
    cuda_safe_call(cudaMalloc((void**)&m2_g, kmemsize));
    cuda_safe_call(cudaMalloc((void**)&m3_g, kmemsize));
    cuda_safe_call(cudaMalloc((void**)&m4_g, kmemsize));
    cuda_safe_call(cudaMalloc((void**)&m5_g, kmemsize));
    cuda_safe_call(cudaMalloc((void**)&m6_g, kmemsize));
    cuda_safe_call(cudaMalloc((void**)&m7_g, kmemsize));

    cuda_safe_call(cudaMemcpy(bmati_g, bmati.data(), imemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bmatj_g, bmatj.data(), jmemsize, cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy(m1_g, m1.data(), kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(m2_g, m2.data(), kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(m3_g, m3.data(), kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(m4_g, m4.data(), kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(m5_g, m5.data(), kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(m6_g, m6.data(), kmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(m7_g, m7.data(), kmemsize, cudaMemcpyHostToDevice));

    make_cufft_plan();
}

template<typename TF>
void Pres_4<TF>::clear_device()
{
    cuda_safe_call(cudaFree(bmati_g));
    cuda_safe_call(cudaFree(bmatj_g));

    cuda_safe_call(cudaFree(m1_g));
    cuda_safe_call(cudaFree(m2_g));
    cuda_safe_call(cudaFree(m3_g));
    cuda_safe_call(cudaFree(m4_g));
    cuda_safe_call(cudaFree(m5_g));
    cuda_safe_call(cudaFree(m6_g));
    cuda_safe_call(cudaFree(m7_g));
}

template<typename TF>
void Pres_4<TF>::exec(double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    // 1. Create the input for the pressure solver.
    const int blocki = 128;
    const int blockj = 2;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    dim3 grid2dGPU (gridi, gridj);
    dim3 block2dGPU(blocki, blockj);

    const TF dti = 1./dt;

    // Calculate the cyclic BCs first.
    boundary_cyclic.exec_g(fields.mt.at("u")->fld_g);
    boundary_cyclic.exec_g(fields.mt.at("v")->fld_g);
    boundary_cyclic.exec_g(fields.mt.at("w")->fld_g);

    ghost_cells_wt_g<TF><<<grid2dGPU, block2dGPU>>>(
        fields.mt.at("w")->fld_g,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    pres_in_g<TF><<<gridGPU, blockGPU>>>(
        fields.sd.at("p")->fld_g,
        fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
        fields.mt.at("u")->fld_g, fields.mt.at("v")->fld_g, fields.mt.at("w")->fld_g,
        gd.dzi4_g, gd.dxi, gd.dyi, dti,
        gd.icells, gd.ijcells,
        gd.imax, gd.imax*gd.jmax,
        gd.imax, gd.jmax, gd.kmax,
        gd.igc,  gd.jgc,  gd.kgc);
    cuda_check_error();

    // Get two free tmp fields on GPU.
    auto tmp1 = fields.get_tmp_g();
    auto tmp2 = fields.get_tmp_g();

    fft_forward(fields.sd.at("p")->fld_g, tmp1->fld_g, tmp2->fld_g);

    // Set jslice to a higher value
    const int jslice = std::max(gd.jblock/4, 1);

    const int blockis = 128;
    const int blockjs = 1;
    const int gridis  = gd.iblock/blockis + (gd.iblock%blockis > 0);
    const int gridjs  =    jslice/blockjs + (   jslice%blockjs > 0);

    dim3 grid2dsGPU (gridis , gridjs );
    dim3 block2dsGPU(blockis, blockjs);

    const int ns = gd.iblock*jslice*(gd.kmax+4);
    const int nj = gd.jblock/jslice;

    // Shortcuts for short calls.
    TF* tmp1_g = tmp1->fld_g;
    TF* tmp2_g = tmp2->fld_g;

    for (int n=0; n<nj; ++n)
    {
        // Prepare the fields that go into the matrix solver
        solve_in_g<TF><<<grid2dsGPU,block2dsGPU>>>(
            fields.sd.at("p")->fld_g,
            m1_g, m2_g, m3_g, m4_g,
            m5_g, m6_g, m7_g,
            &tmp1_g[0*ns], &tmp1_g[1*ns], &tmp1_g[2*ns], &tmp1_g[3*ns],
            &tmp2_g[0*ns], &tmp2_g[1*ns], &tmp2_g[2*ns], &tmp2_g[3*ns],
            bmati_g, bmatj_g,
            md.mpicoordx, md.mpicoordy,
            gd.iblock, gd.jblock,
            gd.kmax,
            n, jslice);
        cuda_check_error();

        // Solve the sevenbanded matrix
        hdma_g<<<grid2dsGPU,block2dsGPU>>>(
            &tmp1_g[0*ns], &tmp1_g[1*ns], &tmp1_g[2*ns], &tmp1_g[3*ns],
            &tmp2_g[0*ns], &tmp2_g[1*ns], &tmp2_g[2*ns], &tmp2_g[3*ns],
            gd.iblock, gd.kmax, jslice);
        cuda_check_error();

        // Put the solution back into the pressure field
        solve_put_back_g<TF><<<grid2dsGPU,block2dsGPU>>>(
            fields.sd.at("p")->fld_g,
            &tmp2_g[3*ns],
            gd.iblock, gd.jblock,
            gd.kmax,
            n, jslice);
        cuda_check_error();
    }

    fft_backward(fields.sd.at("p")->fld_g, tmp1->fld_g, tmp2->fld_g);

    cuda_safe_call(cudaMemcpy(tmp1->fld_g, fields.sd.at("p")->fld_g, gd.ncells*sizeof(TF), cudaMemcpyDeviceToDevice));

    solve_out_g<TF><<<gridGPU, blockGPU>>>(
        fields.sd.at("p")->fld_g, tmp1->fld_g,
        gd.imax, gd.imax*gd.jmax,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.imax,   gd.jmax,   gd.kmax);
    cuda_check_error();

    boundary_cyclic.exec_g(fields.sd.at("p")->fld_g);

    // 3. Get the pressure tendencies from the pressure field.
    pres_out_g<TF><<<gridGPU, blockGPU>>>(
        fields.mt.at("u")->fld_g, fields.mt.at("v")->fld_g, fields.mt.at("w")->fld_g,
        fields.sd.at("p")->fld_g,
        gd.dzhi4_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    fields.release_tmp_g(tmp1);
    fields.release_tmp_g(tmp2);

    cudaDeviceSynchronize();
    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);
}

template<typename TF>
TF Pres_4<TF>::check_divergence()
{
    auto& gd = grid.get_grid_data();

    const int blocki = 128;
    const int blockj = 2;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    auto div = fields.get_tmp_g();

    calc_divergence_g<TF><<<gridGPU, blockGPU>>>(
            div->fld_g,
            fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
            gd.dzi4_g,
            gd.dxi, gd.dyi,
            gd.icells, gd.ijcells,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend);
    cuda_check_error();

    TF divmax = field3d_operators.calc_max_g(div->fld_g);
    // TO-DO: add grid.get_max() or similar for future parallel versions

    fields.release_tmp_g(div);

    return divmax;
}
#endif

template class Pres_4<double>;
template class Pres_4<float>;
