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

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "tools.h"
#include "budget3d.h"

namespace
{
    template<typename TF> __global__
    void accumulate_tend_g(
            TF* const __restrict__ acc, const TF* const __restrict__ tend, const TF sub_dt,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            acc[ijk] += sub_dt * tend[ijk];
        }
    }
}

#ifdef USECUDA
template<typename TF>
void Budget3d<TF>::exec(const double sub_dt, const double dt, const bool in_substep)
{
    if (!swbudget3d)
        return;

    auto& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    const TF sub_dt_tf = static_cast<TF>(sub_dt);

    for (auto& varname : tendencylist)
    {
        auto& acc  = *fields.sd.at(varname + "_tend");
        auto& tend = *fields.at.at(varname);

        accumulate_tend_g<TF><<<gridGPU, blockGPU>>>(
                acc.fld_g, tend.fld_g, sub_dt_tf,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
        cuda_check_error();
    }

    // Sum of cB over the RK substeps is not 1; count the full time step once.
    if (!in_substep)
        elapsed_time += dt;
}
#endif


#ifdef FLOAT_SINGLE
template class Budget3d<float>;
#else
template class Budget3d<double>;
#endif
