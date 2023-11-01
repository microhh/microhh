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

#ifndef FORCE_KERNELS_CUH
#define FORCE_KERNELS_CUH

#include "finite_difference.h"
#include "cuda_tiling.h"

namespace Force_kernels
{
    using namespace Finite_difference::O2;

    template<typename TF>
    struct advec_wls_2nd_mean_g
    {
        DEFINE_GRID_KERNEL("force::advec_wls_2nd_mean", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g, const int i, const int j, const int k, const Level level,
                TF* const __restrict__ st,
                const TF* const __restrict__ s,
                const TF* const __restrict__ wls,
                const TF* const __restrict__ dzhi)
        {
            const int ijk = g(i, j, k);

            if (wls[k] > TF(0))
                st[ijk] -= wls[k] * (s[k]-s[k-1])*dzhi[k];
            else
                st[ijk] -= wls[k] * (s[k+1]-s[k])*dzhi[k+1];
        }
    };


    template<typename TF>
    struct advec_wls_2nd_local_g
    {
        DEFINE_GRID_KERNEL("force::advec_wls_2nd_local", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g, const int i, const int j, const int k, const Level level,
                TF* const __restrict__ st,
                const TF* const __restrict__ s,
                const TF* const __restrict__ wls,
                const TF* const __restrict__ dzhi)
        {
            const int kk = g.kstride;
            const int ijk = g(i, j, k);

            if (wls[k] > TF(0))
                st[ijk] -= wls[k] * (s[ijk]-s[ijk-kk])*dzhi[k];
            else
                st[ijk] -= wls[k] * (s[ijk+kk]-s[ijk])*dzhi[k+1];
        }
    };


    template<typename TF>
    struct advec_wls_2nd_local_w_g
    {
        DEFINE_GRID_KERNEL("force::advec_wls_2nd_local_w", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g, const int i, const int j, const int k, const Level level,
                TF* const __restrict__ st,
                const TF* const __restrict__ s,
                const TF* const __restrict__ wls,
                const TF* const __restrict__ dzi)
        {
            const int kk = g.kstride;
            const int ijk = g(i, j, k);

            if (interp2( wls[k-1], wls[k] ) > 0.)
                st[ijk] -= interp2( wls[k-1], wls[k] ) * (s[ijk]-s[ijk-kk])*dzi[k-1];
            else
                st[ijk] -= interp2( wls[k-1], wls[k] ) * (s[ijk+kk]-s[ijk])*dzi[k];
        }
    };


    template<typename TF>
    struct coriolis_2nd_g
    {
        DEFINE_GRID_KERNEL("force::coriolis_2nd", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g, const int i, const int j, const int k, const Level level,
                TF* const __restrict__ ut,
                TF* const __restrict__ vt,
                const TF* const __restrict__ u,
                const TF* const __restrict__ v,
                const TF* const __restrict__ ug,
                const TF* const __restrict__ vg,
                const TF fc,
                const TF ugrid,
                const TF vgrid)
        {
            const int ii = 1;
            const int jj = g.jstride;
            const int ijk = g(i, j, k);

            ut[ijk] += fc * (TF(0.25)*(v[ijk-ii] + v[ijk] + v[ijk-ii+jj] + v[ijk+jj]) + vgrid - vg[k]);
            vt[ijk] -= fc * (TF(0.25)*(u[ijk-jj] + u[ijk] + u[ijk+ii-jj] + u[ijk+ii]) + ugrid - ug[k]);
        }
    };


    template<typename TF>
    struct add_profile_g
    {
        DEFINE_GRID_KERNEL("force::add_profile", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g, const int i, const int j, const int k, const Level level,
                TF* const __restrict__ st,
                const TF* const __restrict__ prof)
        {
            const int ijk = g(i, j, k);
            st[ijk] += prof[k];
        }
    };
}

#endif //FORCE_KERNELS_CUH
