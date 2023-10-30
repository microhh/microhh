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

namespace Pres_2_kernel
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
                const TF dxi, const TF dyi, const TF dti)
        {
            // Strides over field including ghost cells.
            const int ii = 1;
            const int jj = g.jstride;
            const int kk = g.kstride;

            // Strides over field excluding ghost cells.
            const int imax = g.iend-g.istart;
            const int jmax = g.jend-g.jstart;
            const int jjp = imax;
            const int kkp = imax*jmax;

            // Calculate ghost cells from start index.
            const int igc = g.istart;
            const int jgc = g.jstart;
            const int kgc = g.kstart;

            const int ijk  = i + j*jj + k*kk;
            const int ijkp = (i-igc) + (j-jgc)*jjp + (k-kgc)*kkp;

            p[ijkp] = rhoref [k+kgc]   * ((ut[ijk+ii] + u[ijk+ii] * dti) - (ut[ijk] + u[ijk] * dti)) * dxi
                    + rhoref [k+kgc]   * ((vt[ijk+jj] + v[ijk+jj] * dti) - (vt[ijk] + v[ijk] * dti)) * dyi
                  + ( rhorefh[k+kgc+1] * ( wt[ijk+kk] + w[ijk+kk] * dti)
                    - rhorefh[k+kgc  ] * ( wt[ijk   ] + w[ijk   ] * dti)) * dzi[k+kgc];
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
