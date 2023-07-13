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

#ifndef TIMELOOP_KERNELS_CUH
#define TIMELOOP_KERNELS_CUH

#include "cuda_tiling.h"

namespace timeloop
{
    template<typename TF, int substep>
    struct rk3_g
    {
        DEFINE_GRID_KERNEL("timeloop::rk3", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g, const int i, const int j, const int k, const Level level,
                TF* __restrict__ a, TF* __restrict__ at, const double dt)
        {
            constexpr TF cA1 = -5./9.;
            constexpr TF cA2 = -153./128.;

            constexpr TF cB0 =  1./ 3.;
            constexpr TF cB1 = 15./16.;
            constexpr TF cB2 =  8./15.;

            const int ijk = g(i, j, k);

            switch (substep)
            {
                case 0:
                    a [ijk] = a[ijk] + cB0*dt*at[ijk];
                    at[ijk] = cA1*at[ijk];
                    break;
                case 1:
                    a [ijk] = a[ijk] + cB1*dt*at[ijk];
                    at[ijk] = cA2*at[ijk];
                    break;
                case 2:
                    a [ijk] = a[ijk] + cB2*dt*at[ijk];
                    at[ijk] = TF(0.);
                    break;
            }
        }
    };

    template<typename TF, int substep>
    struct rk4_g
    {
        DEFINE_GRID_KERNEL("timeloop::rk4", 0)

        template <typename Level>
        CUDA_DEVICE
        void operator()(
                Grid_layout g, const int i, const int j, const int k, const Level level,
                TF* __restrict__ a, TF* __restrict__ at, const double dt)
        {

            constexpr TF cA1 = - 567301805773./1357537059087.;
            constexpr TF cA2 = -2404267990393./2016746695238.;
            constexpr TF cA3 = -3550918686646./2091501179385.;
            constexpr TF cA4 = -1275806237668./ 842570457699.;

            constexpr TF cB0 = 1432997174477./ 9575080441755.;
            constexpr TF cB1 = 5161836677717./13612068292357.;
            constexpr TF cB2 = 1720146321549./ 2090206949498.;
            constexpr TF cB3 = 3134564353537./ 4481467310338.;
            constexpr TF cB4 = 2277821191437./14882151754819.;

            const int ijk = g(i, j, k);

            switch (substep)
            {
                case 0:
                    a [ijk] = a[ijk] + cB0*dt*at[ijk];
                    at[ijk] = cA1*at[ijk];
                    break;
                case 1:
                    a [ijk] = a[ijk] + cB1*dt*at[ijk];
                    at[ijk] = cA2*at[ijk];
                    break;
                case 2:
                    a [ijk] = a[ijk] + cB2*dt*at[ijk];
                    at[ijk] = cA3*at[ijk];
                    break;
                case 3:
                    a [ijk] = a[ijk] + cB3*dt*at[ijk];
                    at[ijk] = cA4*at[ijk];
                    break;
                case 4:
                    a [ijk] = a[ijk] + cB4*dt*at[ijk];
                    at[ijk] = 0;
                    break;
            }
        }
    };
}

#endif //MICROHHC_TIMELOOP_KERNELS_CUH
