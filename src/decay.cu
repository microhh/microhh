/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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


template<typename TF>
void Decay<TF>::get_mask(Stats<TF>& stats, std::string mask_name)
{
    using namespace Tools_g;
    auto& gd = grid.get_grid_data();

    if (mask_name == "couvreux")
    {
        if (fields.sp.find("couvreux") == fields.sp.end())
        {
            std::string message = "Couvreux mask not available without couvreux scalar";
            throw std::runtime_error(message);
        }
        auto couvreux = fields.get_tmp_g();
        auto couvreuxh = fields.get_tmp_g();

        // Calculate mean and standard deviation of couvreux scalar

        // Subtract mean and 1 std dev from couvreux scalar

        // Interpolate to half levels
        auto threads = grid.get_dim_gpu(gd.imax, gd.jmax, gd.kmax);
        interpolate_2nd_g<<<threads.first, threads.second>>>(couvreuxh->fld_g.data(), couvreux->fld_g.data(), gd.sloc[0] - gd.wloc[0], gd.sloc[1] - gd.wloc[1], gd.sloc[2] - gd.wloc[2], gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);

        //Calculate mask
        stats.set_mask_thres(mask_name, *couvreux, *couvreuxh, 0., Stats_mask_type::Plus);

        fields.release_tmp_g(couvreux);
        fields.release_tmp_g(couvreuxh);
    }

    else
    {
        std::string message = "Decay can not provide mask: \"" + mask_name +"\"";
        throw std::runtime_error(message);
    }
}

#endif


#ifdef FLOAT_SINGLE
template class Decay<float>;
#else
template class Decay<double>;
#endif
