/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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
#include "thermo.h"
#include "boundary_cyclic.h"
#include "thermo_moist_functions.h"
#include "constants.h"
#include "tools.h"

#include "microphys.h"
#include "microphys_2mom_warm.h"

using namespace Constants;
using namespace Thermo_moist_functions;
using namespace Micro_2mom_warm_constants;
using namespace Micro_2mom_warm_functions;

namespace
{
    template<typename TF> __global__
    void remove_negative_values_g(TF* __restrict__ field,
                                  int istart, int jstart, int kstart,
                                  int iend,   int jend,   int kend,
                                  int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            field[ijk] = fmax(field[ijk], TF(0));
        }
    }
}

#ifdef USECUDA
template<typename TF>
void Microphys_2mom_warm<TF>::exec(Thermo<TF>& thermo, const double dt)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, gd.kmax);
    dim3 blockGPU(blocki, blockj, 1);

    // Remove negative values from the qr and nr fields
    remove_negative_values_g<<<gridGPU, blockGPU>>>(fields.ap.at("qr")->fld_g,
        gd.istart, gd.jstart, gd.kstart, gd.iend, gd.jend, gd.kend, gd.icells, gd.ijcells);
    cuda_check_error();

    remove_negative_values_g<<<gridGPU, blockGPU>>>(fields.ap.at("nr")->fld_g,
        gd.istart, gd.jstart, gd.kstart, gd.iend, gd.jend, gd.kend, gd.icells, gd.ijcells);
    cuda_check_error();

}
#endif

#ifdef USECUDA
template<typename TF>
unsigned long Microphys_2mom_warm<TF>::get_time_limit(unsigned long idt, const double dt)
{
    return Constants::ulhuge;
}
#endif

template class Microphys_2mom_warm<double>;
template class Microphys_2mom_warm<float>;
