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

#ifndef PRES_2_KERNELS_CUH
#define PRES_2_KERNELS_CUH

#include "cuda_tiling.h"

namespace Pres_2_kernels
{
    template<typename TF>
    struct pres_in_g
    {
        DEFINE_GRID_KERNEL("pres_2::pres_in", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g,
                const int i, const int j, const int k,
                const Level level,
                TF* const __restrict__ p,
                const TF* const __restrict__ u,
                const TF* const __restrict__ v,
                const TF* const __restrict__ w ,
                const TF* const __restrict__ ut,
                const TF* const __restrict__ vt,
                const TF* const __restrict__ wt,
                const TF* const __restrict__ dzi,
                const TF* const __restrict__ rhoref,
                const TF* const __restrict__ rhorefh,
                const TF dxi, const TF dyi, const TF dti,
                const int icells, const int ijcells,
                const int igc, const int jgc, const int kgc)
        {
            const int ii = 1;
            const int jj = icells;
            const int kk = ijcells;

            const int ijkp = g(i, j, k);
            const int ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;

            p[ijkp] = rhoref [k+kgc]   * ((ut[ijk+ii] + u[ijk+ii] * dti) - (ut[ijk] + u[ijk] * dti)) * dxi
                    + rhoref [k+kgc]   * ((vt[ijk+jj] + v[ijk+jj] * dti) - (vt[ijk] + v[ijk] * dti)) * dyi
                  + ( rhorefh[k+kgc+1] * ( wt[ijk+kk] + w[ijk+kk] * dti)
                    - rhorefh[k+kgc  ] * ( wt[ijk   ] + w[ijk   ] * dti)) * dzi[k+kgc];
        }
    };


    template<typename TF>
    struct solve_in_g
    {
        DEFINE_GRID_KERNEL("pres_2::solve_in", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g,
                const int i, const int j, const int k,
                const Level level,
                TF* const __restrict__ p,
                const TF* const __restrict__ work3d,
                TF* const __restrict__ b,
                const TF* const __restrict__ a,
                const TF* const __restrict__ c,
                const TF* const __restrict__ dz,
                const TF* const __restrict__ rhoref,
                const TF* const __restrict__ bmati,
                const TF* const __restrict__ bmatj,
                const int kstart, const int kmax)
        {
            const int ijk = g(i, j, k);

            // CvH this needs to be taken into account in case of an MPI run
            // iindex = mpi->mpicoordy * iblock + i;
            // jindex = mpi->mpicoordx * jblock + j;
            // b[ijk] = dz[k+kgc]*dz[k+kgc] * (bmati[iindex]+bmatj[jindex]) - (a[k]+c[k]);
            //  if(iindex == 0 && jindex == 0)

            b[ijk] = dz[k+kstart]*dz[k+kstart] * rhoref[k+kstart]*(bmati[i]+bmatj[j]) - (a[k]+c[k]);
            p[ijk] = dz[k+kstart]*dz[k+kstart] * p[ijk];

            if (level.distance_to_start() == 0)
            {
                // Substitute BC's
                // ijk = i + j*jj;
                b[ijk] += a[0];
            }
            else if (level.distance_to_end() == 0)
            {
                // For wave number 0, which contains average, set pressure at top to zero
                if (i == 0 && j == 0)
                    b[ijk] -= c[k];
                else
                    b[ijk] += c[k];
            }
        }
    };


    template<typename TF>
    struct tdma_g
    {
        DEFINE_GRID_KERNEL("pres_2::tdma_g", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g,
                const int i, const int j, const int k,
                const Level level,
                const TF* const __restrict__ a,
                const TF* const __restrict__ b,
                const TF* const __restrict__ c,
                TF* const __restrict__ p,
                TF* const __restrict__ work3d,
                const int kmax)
        {
            const int ij = g(i, j, k);  // k=0
            const int kk = g.kstride;

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
    };


    template<typename TF>
    struct solve_out_g
    {
        DEFINE_GRID_KERNEL("pres_2::solve_out", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g,
                const int i, const int j, const int k,
                const Level level,
                TF* const __restrict__ p,
                const TF* const __restrict__ work3d,
                const int istart, const int jstart, const int kstart,
                const int jj_gc, const int kk_gc)
        {
            // Strides without ghost cells;
            const int jj = g.jstride;
            const int kk = g.kstride;

            // Strides with ghost cells:
            const int jjp = jj_gc;
            const int kkp = kk_gc;

            const int ijk  = i + j*jj + k*kk;
            const int ijkp = i+istart + (j+jstart)*jjp + (k+kstart)*kkp;

            p[ijkp] = work3d[ijk];

            if (level.distance_to_start() == 0)
                p[ijkp-kkp] = p[ijkp];
        }
    };


    template<typename TF>
    struct pres_out_g
    {
        DEFINE_GRID_KERNEL("pres_2::pres_out", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g,
                const int i, const int j, const int k,
                const Level level,
                TF* const __restrict__ ut,
                TF* const __restrict__ vt,
                TF* const __restrict__ wt,
                const TF* const __restrict__ p,
                const TF* const __restrict__ dzhi,
                const TF dxi, const TF dyi)
        {
            const int ii = g.istride;
            const int jj = g.jstride;
            const int kk = g.kstride;

            const int ijk = g(i, j, k);

            ut[ijk] -= (p[ijk] - p[ijk-ii]) * dxi;
            vt[ijk] -= (p[ijk] - p[ijk-jj]) * dyi;
            wt[ijk] -= (p[ijk] - p[ijk-kk]) * dzhi[k];
        }
    };
}
#endif //MICROHHC_PRES_2_KERNELS_CUH
