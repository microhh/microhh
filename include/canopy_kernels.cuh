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

#ifndef CANOPY_KERNELS_CUH
#define CANOPY_KERNELS_CUH

#include "cuda_tiling.h"
#include "fast_math.h"

namespace Canopy_kernels
{
    namespace fm = Fast_math;

    template<typename TF, bool sw_pad_3d>
    struct canopy_drag_u_g
    {
        DEFINE_GRID_KERNEL("canopy::canopy_drag_u", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g, const int i, const int j, const int k, const Level level,
                TF* const __restrict__ ut,
                const TF* const __restrict__ u,
                const TF* const __restrict__ v,
                const TF* const __restrict__ w,
                const TF* const __restrict__ pad,
                const TF utrans,
                const TF vtrans,
                const TF cd)
        {
            const int ii = 1;
            const int jj = g.jstride;
            const int kk = g.kstride;
            const int ijk = g(i, j, k);

            // Interpolate `pad` to `u` locations.
            const TF pad_u = sw_pad_3d ? TF(0.5) * (pad[ijk-ii] + pad[ijk]) : pad[k];

            // Interpolate `v` and `w` to `u` locations.
            const TF u_on_u = u[ijk] + utrans;
            const TF v_on_u = TF(0.25) * (v[ijk] + v[ijk-ii] + v[ijk-ii+jj] + v[ijk+jj]) + vtrans;
            const TF w_on_u = TF(0.25) * (w[ijk] + w[ijk-ii] + w[ijk-ii+kk] + w[ijk+kk]);

            const TF ftau = -cd * pad_u *
                sqrt( fm::pow2(u_on_u) +
                      fm::pow2(v_on_u) +
                      fm::pow2(w_on_u) );

            ut[ijk] += ftau * u_on_u;
        }
    };


    template<typename TF, bool sw_pad_3d>
    struct canopy_drag_v_g
    {
        DEFINE_GRID_KERNEL("canopy::canopy_drag_v", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g, const int i, const int j, const int k, const Level level,
                TF* const __restrict__ vt,
                const TF* const __restrict__ u,
                const TF* const __restrict__ v,
                const TF* const __restrict__ w,
                const TF* const __restrict__ pad,
                const TF utrans,
                const TF vtrans,
                const TF cd)
        {
            const int ii = 1;
            const int jj = g.jstride;
            const int kk = g.kstride;
            const int ijk = g(i, j, k);

            // Interpolate `pad` to `v` locations.
            const TF pad_v = sw_pad_3d ? TF(0.5) * (pad[ijk-jj] + pad[ijk]) : pad[k];

            // Interpolate `u` and `w` to `v` locations.
            const TF u_on_v = TF(0.25) * (u[ijk] + u[ijk+ii] + u[ijk+ii-jj] + u[ijk-jj]) + utrans;
            const TF v_on_v = v[ijk] + vtrans;
            const TF w_on_v = TF(0.25) * (w[ijk] + w[ijk+kk] + w[ijk-jj] + w[ijk-jj+kk]);

            const TF ftau = -cd * pad_v *
                sqrt( fm::pow2(u_on_v) +
                      fm::pow2(v_on_v) +
                      fm::pow2(w_on_v) );

            vt[ijk] += ftau * v_on_v;
        }
    };


    template<typename TF, bool sw_pad_3d>
    struct canopy_drag_w_g
    {
        DEFINE_GRID_KERNEL("canopy::canopy_drag_w", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g, const int i, const int j, const int k, const Level level,
                TF* const __restrict__ wt,
                const TF* const __restrict__ u,
                const TF* const __restrict__ v,
                const TF* const __restrict__ w,
                const TF* const __restrict__ pad,
                const TF utrans,
                const TF vtrans,
                const TF cd)
        {
            const int ii = 1;
            const int jj = g.jstride;
            const int kk = g.kstride;
            const int ijk = g(i, j, k);

            // Interpolate `pad` to `w` locations.
            const TF pad_w = sw_pad_3d ? TF(0.5) * (pad[ijk-kk] + pad[ijk]) : TF(0.5) * (pad[k-1] + pad[k]);

            // Interpolate `u` and `v` to `w` locations.
            const TF u_on_w = TF(0.25) * (u[ijk] + u[ijk+ii] + u[ijk+ii-kk] + u[ijk-kk]) + utrans;
            const TF v_on_w = TF(0.25) * (v[ijk] + v[ijk+jj] + v[ijk+jj-kk] + v[ijk-kk]) + vtrans;
            const TF w_on_w = w[ijk];

            const TF ftau = -cd * pad_w *
                sqrt( fm::pow2(u_on_w) +
                      fm::pow2(v_on_w) +
                      fm::pow2(w_on_w) );

            wt[ijk] += ftau * w_on_w;
        }
    };
}
#endif // CANOPY_KERNELS_CUH
