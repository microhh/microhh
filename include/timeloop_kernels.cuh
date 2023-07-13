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
}

#endif //MICROHHC_TIMELOOP_KERNELS_CUH
