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

#include <algorithm>
#include <iostream>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "limiter.h"

template<typename TF>
Limiter<TF>::Limiter(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    limit_list = inputin.get_list<std::string>("limiter", "limitlist", "", std::vector<std::string>());
}

template <typename TF>
Limiter<TF>::~Limiter()
{
}

template <typename TF>
void Limiter<TF>::create(Stats<TF>& stats)
{
    for (const std::string& s : limit_list)
        stats.add_tendency(*fields.at.at(s), "z", tend_name, tend_longname);
}

namespace
{
    // This function produces a tendency that represents a source that avoids sub zero values.
    template<typename TF>
    void tendency_limiter(
            TF* restrict at, const TF* restrict a, const TF dt,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const TF dti = TF(1.)/dt;

        // Add epsilon, to make sure the final result ends just above zero.
        // NOTE: don't use `eps<TF>` here; `eps<float>` is too large
        //       as a lower limit for e.g. hydrometeors or chemical species.
        constexpr TF eps = std::numeric_limits<double>::epsilon();

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const TF a_new = a[ijk] + dt*at[ijk];
                    at[ijk] += (a_new < TF(0.)) ? (-a_new + eps) * dti : TF(0.);
                }
    }
}

#ifndef USECUDA
template <typename TF>
void Limiter<TF>::exec(double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    for (auto& name : limit_list)
    {
         tendency_limiter<TF>(
                fields.at.at(name)->fld.data(), fields.ap.at(name)->fld.data(), dt,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_tend(*fields.at.at(name), tend_name);
    }
}
#endif


#ifdef FLOAT_SINGLE
template class Limiter<float>;
#else
template class Limiter<double>;
#endif
