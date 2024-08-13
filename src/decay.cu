/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

//#include <cstdio>
#include <algorithm>
#include <iostream>
#include <math.h>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "tools.h"
#include "stats.h"
#include "decay.h"

namespace
{
    template<typename TF> __global__
    void enforce_exponential_decay_g(TF* const __restrict__ tend, TF* const __restrict__ var, const TF rate,
                            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int jj, const int kk)
    {
        const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k  = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            tend[ijk] -= rate * var[ijk];
        }
    }
}

#ifdef USECUDA
template <typename TF>
void Decay<TF>::exec(double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    for (auto& it : dmap)
    {
        if (it.second.type == Decay_type::exponential)
        {
            const TF rate = 1./(std::max(it.second.timescale, dt));
            enforce_exponential_decay_g<TF><<<gridGPU, blockGPU>>>(
                fields.st.at(it.first)->fld_g, fields.sp.at(it.first)->fld_g, rate,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
            cuda_check_error();

            cudaDeviceSynchronize();
            stats.calc_tend(*fields.st.at(it.first), tend_name);
        }
    }

}
#endif


#ifdef FLOAT_SINGLE
template class Decay<float>;
#else
template class Decay<double>;
#endif
