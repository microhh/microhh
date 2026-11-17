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

#include <algorithm>
#include <iostream>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "limiter.h"
#include "constants.h"
#include "diff.h"

namespace
{
    // This function produces a tendency that represents a source that avoids sub zero values.
    template<typename TF>
    void limit_tendency(
            TF* const restrict at,
            const TF* const restrict a,
            const TF min_value, const TF dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const TF dti = TF(1.)/dt;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    const TF a_new = a[ijk] + dt*at[ijk];
                    at[ijk] += (a_new < min_value) ? (-a_new + min_value) * dti : TF(0.);
                }
    }

    // This function hard clips a field to guarantee it is at/above the minimum value.
    template<typename TF>
    void clip_field(
            TF* const restrict a,
            const TF min_value,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    a[ijk] = std::max(a[ijk], min_value);
                }
    }
}

template<typename TF>
Limiter<TF>::Limiter(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Diff<TF>& diffin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), diff(diffin)
{
    limit_list = inputin.get_list<std::string>("limiter", "limitlist", "", std::vector<std::string>());
    clip_list = inputin.get_list<std::string>("limiter", "cliplist", "", std::vector<std::string>());
    limit_sgstke = (diffin.get_switch() == Diffusion_type::Diff_tke2);
}

template <typename TF>
Limiter<TF>::~Limiter()
{
}

template <typename TF>
void Limiter<TF>::create(Stats<TF>& stats)
{
    // Check if all limit fields are valid.
    for (auto& name : limit_list)
        if (fields.at.find(name) == fields.at.end())
        {
            std::string error = "Non-existing prognostic field \"" + name + "\" in limiter!";
            throw std::runtime_error(error);
        }

    // Check if all clip fields are valid.
    for (auto& name : clip_list)
        if (fields.ap.find(name) == fields.ap.end())
        {
            std::string error = "Non-existing prognostic field \"" + name + "\" in limiter cliplist!";
            throw std::runtime_error(error);
        }

    for (const std::string& s : limit_list)
        stats.add_tendency(*fields.at.at(s), "z", tend_name, tend_longname);

    if (limit_sgstke)
        stats.add_tendency(*fields.at.at("sgstke"), "z", tend_name, tend_longname);
}

#ifndef USECUDA
template <typename TF>
void Limiter<TF>::limit_tendencies(double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    // Add epsilon, to make sure the final result ends just above zero.
    // NOTE: don't use `eps<TF>` here; `eps<float>` is too large
    //       as a lower limit for e.g. hydrometeors or chemical species.
    constexpr TF min_value = std::numeric_limits<double>::epsilon();

    for (auto& name : limit_list)
    {
        limit_tendency<TF>(
                fields.at.at(name)->fld.data(),
                fields.ap.at(name)->fld.data(),
                min_value, dt,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        stats.calc_tend(*fields.at.at(name), tend_name);
    }

    if (limit_sgstke)
    {
        limit_tendency<TF>(
               fields.at.at("sgstke")->fld.data(),
               fields.ap.at("sgstke")->fld.data(),
               Constants::sgstke_min<TF>, dt,
               gd.istart, gd.iend,
               gd.jstart, gd.jend,
               gd.kstart, gd.kend,
               gd.icells, gd.ijcells);

        stats.calc_tend(*fields.at.at("sgstke"), tend_name);
    }
}

template <typename TF>
void Limiter<TF>::clip_fields()
{
    auto& gd = grid.get_grid_data();

    for (auto& name : clip_list)
        clip_field<TF>(
                fields.ap.at(name)->fld.data(),
                TF(0.),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
}
#endif


#ifdef FLOAT_SINGLE
template class Limiter<float>;
#else
template class Limiter<double>;
#endif
