/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
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

#include "tools.h"
#include "grid.h"
#include "timedep.h"
#include "timeloop.h"
namespace
{

    template<typename TF> __global__
    void calc_time_dependent_prof_g(TF* const __restrict__ prof, const TF* const __restrict__ data,
                                    const double fac0, const double fac1,
                                    const int index0,  const int index1,
                                    const int kmax,    const int kgc)
    {
        const int k = blockIdx.x*blockDim.x + threadIdx.x;
        const int kk = kmax;

        if (k < kmax)
            prof[k+kgc] = fac0*data[index0*kk+k] + fac1*data[index1*kk+k];
    }
} // end namespace


#ifdef USECUDA
template<typename TF>
Timedep<TF>::~Timedep()
{
    #ifdef USECUDA
    cuda_safe_call(cudaFree(data_g));
    #endif
}
#endif

#ifdef USECUDA
template <typename TF>
void Timedep<TF>::update_time_dependent_prof_g(TF* prof, Timeloop<TF>& timeloop)
{
    if (sw == Timedep_switch::disabled)
        return;

    auto& gd = grid.get_grid_data();
    const int blockk = 128;
    const int gridk  = gd.kmax/blockk + (gd.kmax%blockk > 0);

    // Get/calculate the interpolation indexes/factors
    int index0 = 0, index1 = 0;
    TF fac0 = 0., fac1 = 0.;
    timeloop.get_interpolation_factors(index0, index1, fac0, fac1, time);

    // Calculate the new vertical profile
    calc_time_dependent_prof_g<<<gridk, blockk>>>(
        prof, data_g, fac0, fac1, index0, index1, gd.kmax, gd.kgc);
    cuda_check_error();
}
#endif

template class Timedep<double>;
template class Timedep<float>;
