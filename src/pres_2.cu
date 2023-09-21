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
#include "pres.h"
#include "pres_2.h"
#include "defines.h"
#include "tools.h"
#include "constants.h"
#include "field3d_operators.h"
#include "stats.h"

namespace
{
    template<typename TF> __global__
    void pres_in_g(TF* __restrict__ p,
                   TF* __restrict__ u ,  TF* __restrict__ v ,     TF* __restrict__ w ,
                   TF* __restrict__ ut,  TF* __restrict__ vt,     TF* __restrict__ wt,
                   TF* __restrict__ dzi, TF* __restrict__ rhoref, TF* __restrict__ rhorefh,
                   TF dxi, TF dyi, TF dti,
                   const int jj, const int kk,
                   const int jjp, const int kkp,
                   const int imax, const int jmax, const int kmax,
                   const int igc, const int jgc, const int kgc)
    {
        const int ii = 1;
        const int i  = blockIdx.x*blockDim.x + threadIdx.x;
        const int j  = blockIdx.y*blockDim.y + threadIdx.y;
        const int k  = blockIdx.z;

        if (i < imax && j < jmax && k < kmax)
        {
            const int ijkp = i + j*jjp + k*kkp;
            const int ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;

            p[ijkp] = rhoref [k+kgc]   * ( (ut[ijk+ii] + u[ijk+ii] * dti) - (ut[ijk] + u[ijk] * dti) ) * dxi
                    + rhoref [k+kgc]   * ( (vt[ijk+jj] + v[ijk+jj] * dti) - (vt[ijk] + v[ijk] * dti) ) * dyi
                  + ( rhorefh[k+kgc+1] * (  wt[ijk+kk] + w[ijk+kk] * dti)
                    - rhorefh[k+kgc  ] * (  wt[ijk   ] + w[ijk   ] * dti) ) * dzi[k+kgc];
        }
    }

    template<typename TF> __global__
    void pres_out_g(TF* __restrict__ ut, TF* __restrict__ vt, TF* __restrict__ wt,
                    TF* __restrict__ p,
                    TF* __restrict__ dzhi, const TF dxi, const TF dyi,
                    const int jj, const int kk,
                    const int istart, const int jstart, const int kstart,
                    const int iend, const int jend, const int kend)
    {
        const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k  = blockIdx.z + kstart;
        const int ii = 1;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            ut[ijk] -= (p[ijk] - p[ijk-ii]) * dxi;
            vt[ijk] -= (p[ijk] - p[ijk-jj]) * dyi;
            wt[ijk] -= (p[ijk] - p[ijk-kk]) * dzhi[k];
        }
    }

    template<typename TF> __global__
    void solve_out_g(TF* __restrict__ p, TF* __restrict__ work3d,
                     const int jj, const int kk,
                     const int jjp, const int kkp,
                     const int istart, const int jstart, const int kstart,
                     const int imax, const int jmax, const int kmax)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;
        const int k = blockIdx.z;

        if (i < imax && j < jmax && k < kmax)
        {
            const int ijk  = i + j*jj + k*kk;
            const int ijkp = i+istart + (j+jstart)*jjp + (k+kstart)*kkp;

            p[ijkp] = work3d[ijk];

            if (k == 0)
                p[ijkp-kkp] = p[ijkp];
        }
    }

    template<typename TF> __global__
    void solve_in_g(TF* __restrict__ p,
                    TF* __restrict__ work3d, TF* __restrict__ b,
                    TF* __restrict__ a, TF* __restrict__ c,
                    TF* __restrict__ dz, TF* __restrict__ rhoref,
                    TF* __restrict__ bmati, TF* __restrict__ bmatj,
                    const int jj, const int kk,
                    const int imax, const int jmax, const int kmax,
                    const int kstart)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;
        const int k = blockIdx.z;

        if (i < imax && j < jmax && k < kmax)
        {
            const int ijk = i + j*jj + k*kk;

            // CvH this needs to be taken into account in case of an MPI run
            // iindex = mpi->mpicoordy * iblock + i;
            // jindex = mpi->mpicoordx * jblock + j;
            // b[ijk] = dz[k+kgc]*dz[k+kgc] * (bmati[iindex]+bmatj[jindex]) - (a[k]+c[k]);
            //  if(iindex == 0 && jindex == 0)

            b[ijk] = dz[k+kstart]*dz[k+kstart] * rhoref[k+kstart]*(bmati[i]+bmatj[j]) - (a[k]+c[k]);
            p[ijk] = dz[k+kstart]*dz[k+kstart] * p[ijk];

            if (k == 0)
            {
                // substitute BC's
                // ijk = i + j*jj;
                b[ijk] += a[0];
            }
            else if (k == kmax-1)
            {
                // for wave number 0, which contains average, set pressure at top to zero
                if (i == 0 && j == 0)
                    b[ijk] -= c[k];
                // set dp/dz at top to zero
                else
                    b[ijk] += c[k];
            }
        }
    }

    template<typename TF> __global__
    void tdma_g(TF* __restrict__ a, TF* __restrict__ b, TF* __restrict__ c,
                TF* __restrict__ p, TF* __restrict__ work3d,
                const int jj, const int kk,
                const int imax, const int jmax, const int kmax)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < imax && j < jmax)
        {
            const int ij = i + j*jj;

            TF work2d = b[ij];
            p[ij] /= work2d;

            for (int k=1; k<kmax; k++)
            {
                const int ijk = ij + k*kk;
                work3d[ijk] = c[k-1] / work2d;
                work2d = b[ijk] - a[k]*work3d[ijk];
                p[ijk] -= a[k]*p[ijk-kk];
                p[ijk] /= work2d;
            }

            for (int k=kmax-2; k>=0; k--)
            {
                const int ijk = ij + k*kk;
                p[ijk] -= work3d[ijk+kk]*p[ijk+kk];
            }
        }
    }

    template<typename TF> __global__
    void calc_divergence_g(TF* __restrict__ u, TF* __restrict__ v, TF* __restrict__ w,
                           TF* __restrict__ div, TF* __restrict__ dzi,
                           TF* __restrict__ rhoref, TF* __restrict__ rhorefh,
                           TF dxi, TF dyi,
                           int jj, int kk, int istart, int jstart, int kstart,
                           int iend, int jend, int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;
        const int ii = 1;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            div[ijk] = rhoref[k]*((u[ijk+ii]-u[ijk])*dxi + (v[ijk+jj]-v[ijk])*dyi)
                    + (rhorefh[k+1]*w[ijk+kk]-rhorefh[k]*w[ijk])*dzi[k];
        }
    }
} // End namespace.

#ifdef USECUDA
template<typename TF>
void Pres_2<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();

    const int kmemsize = gd.kmax*sizeof(TF);
    const int imemsize = gd.itot*sizeof(TF);
    const int jmemsize = gd.jtot*sizeof(TF);

    const int ijmemsize = gd.imax*gd.jmax*sizeof(TF);

    cuda_safe_call(cudaMalloc((void**)&bmati_g,  imemsize));
    cuda_safe_call(cudaMalloc((void**)&bmatj_g,  jmemsize));
    cuda_safe_call(cudaMalloc((void**)&a_g,      kmemsize));
    cuda_safe_call(cudaMalloc((void**)&c_g,      kmemsize));
    cuda_safe_call(cudaMalloc((void**)&work2d_g, ijmemsize));

    cuda_safe_call(cudaMemcpy(bmati_g,  bmati.data(),  imemsize,  cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(bmatj_g,  bmatj.data(),  jmemsize,  cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(a_g,      a.data(),      kmemsize,  cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(c_g,      c.data(),      kmemsize,  cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(work2d_g, work2d.data(), ijmemsize, cudaMemcpyHostToDevice));

    make_cufft_plan();
}

template<typename TF>
void Pres_2<TF>::clear_device()
{
    cuda_safe_call(cudaFree(bmati_g ));
    cuda_safe_call(cudaFree(bmatj_g ));
    cuda_safe_call(cudaFree(a_g     ));
    cuda_safe_call(cudaFree(c_g     ));
    cuda_safe_call(cudaFree(work2d_g));
}

template<typename TF>
void Pres_2<TF>::exec(double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);
    const TF dti = TF(1.)/dt;

    // 3D grid
    dim3 gridGPU (gridi,  gridj,  gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    // 2D grid
    dim3 grid2dGPU (gridi,  gridj);
    dim3 block2dGPU(blocki, blockj);

    // Get two free tmp fields on gpu
    auto tmp1 = fields.get_tmp_g();
    auto tmp2 = fields.get_tmp_g();

    // calculate the cyclic BCs first
    boundary_cyclic.exec_g(fields.mt.at("u")->fld_g);
    boundary_cyclic.exec_g(fields.mt.at("v")->fld_g);
    boundary_cyclic.exec_g(fields.mt.at("w")->fld_g);

    pres_in_g<TF><<<gridGPU, blockGPU>>>(
        fields.sd.at("p")->fld_g,
        fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g,
        fields.mt.at("u")->fld_g, fields.mt.at("v")->fld_g, fields.mt.at("w")->fld_g,
        gd.dzi_g, fields.rhoref_g, fields.rhorefh_g, gd.dxi, gd.dyi, static_cast<TF>(dti),
        gd.icells, gd.ijcells,
        gd.imax, gd.imax*gd.jmax,
        gd.imax, gd.jmax, gd.kmax,
        gd.igc, gd.jgc, gd.kgc);
    cuda_check_error();

    fft_forward(fields.sd.at("p")->fld_g, tmp1->fld_g, tmp2->fld_g);

    solve_in_g<TF><<<gridGPU, blockGPU>>>(
        fields.sd.at("p")->fld_g,
        tmp1->fld_g, tmp2->fld_g,
        a_g, c_g, gd.dz_g, fields.rhoref_g, bmati_g, bmatj_g,
        gd.imax, gd.imax*gd.jmax,
        gd.imax, gd.jmax, gd.kmax,
        gd.kstart);
    cuda_check_error();

    tdma_g<TF><<<grid2dGPU, block2dGPU>>>(
        a_g, tmp2->fld_g, c_g,
        fields.sd.at("p")->fld_g, tmp1->fld_g,
        gd.imax, gd.imax*gd.jmax,
        gd.imax, gd.jmax, gd.kmax);
    cuda_check_error();

    fft_backward(fields.sd.at("p")->fld_g, tmp1->fld_g, tmp2->fld_g);

    cuda_safe_call(cudaMemcpy(tmp1->fld_g, fields.sd.at("p")->fld_g, gd.ncells*sizeof(TF), cudaMemcpyDeviceToDevice));

    solve_out_g<TF><<<gridGPU, blockGPU>>>(
        fields.sd.at("p")->fld_g, tmp1->fld_g,
        gd.imax, gd.imax*gd.jmax,
        gd.icells, gd.ijcells,
        gd.istart, gd.jstart, gd.kstart,
        gd.imax, gd.jmax, gd.kmax);
    cuda_check_error();

    boundary_cyclic.exec_g(fields.sd.at("p")->fld_g);

    pres_out_g<TF><<<gridGPU, blockGPU>>>(
        fields.mt.at("u")->fld_g, fields.mt.at("v")->fld_g, fields.mt.at("w")->fld_g,
        fields.sd.at("p")->fld_g,
        gd.dzhi_g, TF(1.)/gd.dx, TF(1.)/gd.dy,
        gd.icells, gd.ijcells,
        gd.istart,  gd.jstart, gd.kstart,
        gd.iend,    gd.jend,   gd.kend);
    cuda_check_error();

    fields.release_tmp_g(tmp1);
    fields.release_tmp_g(tmp2);
    
    cudaDeviceSynchronize();
    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);
}

template<typename TF>
TF Pres_2<TF>::check_divergence()
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    auto divergence = fields.get_tmp_g();

    calc_divergence_g<<<gridGPU, blockGPU>>>(
        fields.mp.at("u")->fld_g, fields.mp.at("v")->fld_g, fields.mp.at("w")->fld_g, divergence->fld_g,
        gd.dzi_g, fields.rhoref_g, fields.rhorefh_g, gd.dxi, gd.dyi,
        gd.icells, gd.ijcells,
        gd.istart,  gd.jstart, gd.kstart,
        gd.iend, gd.jend, gd.kend);
    cuda_check_error();

    TF divmax = field3d_operators.calc_max_g(divergence->fld_g);
    // TO-DO: add grid.get_max() or similar for future parallel versions

    fields.release_tmp_g(divergence);

    return divmax;
}
#endif


#ifdef FLOAT_SINGLE
template class Pres_2<float>;
#else
template class Pres_2<double>;
#endif
