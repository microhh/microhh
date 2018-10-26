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

#include <algorithm>

#include "data_block.h"
#include "timedep.h"

template<typename TF>
Timedep<TF>::Timedep(Master& masterin, Grid<TF>& gridin, const std::string varnamein, const bool is_timedep) :
    master(masterin), grid(gridin), varname(varnamein)
{
    if (is_timedep)
        sw = Timedep_switch::Enabled;
    else
        sw = Timedep_switch::Disabled;
}

template<typename TF>
Timedep<TF>::~Timedep()
{
}

template <typename TF>
void Timedep<TF>::create_timedep_prof(Netcdf_handle& input_nc)
{
    if (sw == Timedep_switch::Disabled)
        return;

    auto& gd = grid.get_grid_data();
    Data_block data_block(master, varname+".time");
    std::vector<std::string> headers = data_block.get_headers();

    // Sort the times
    headers.erase(headers.begin());
    std::sort(headers.begin(), headers.end());
    std::vector<TF> tmp(gd.kmax);
    for (auto& it : headers)
    {
        time.push_back(std::stod(it));
        data_block.get_vector(tmp, it, gd.kmax, 0, 0);
        data.insert(data.end(), tmp.begin(), tmp.end());
    }
    #ifdef USECUDA
    prepare_device();
    #endif
}

template <typename TF>
void Timedep<TF>::create_timedep(Netcdf_handle& input_nc)
{
    if (sw == Timedep_switch::Disabled)
        return;

    Data_block data_block(master, varname+".time");

    int length = data_block.get_vector_length("time");
    time.resize(length);
    data.resize(length);
    data_block.get_vector(time, "time", length, 0, 0);
    data_block.get_vector(data, varname, length, 0, 0);
    #ifdef USECUDA
    prepare_device();
    #endif
}

template <typename TF>
void Timedep<TF>::update_time_dependent_prof(std::vector<TF>& prof, Timeloop<TF>& timeloop)
{
    if (sw == Timedep_switch::Disabled)
        return;

    auto& gd = grid.get_grid_data();
    const int kk = gd.kmax;
    const int kgc = gd.kgc;

    // Get/calculate the interpolation indexes/factors
    Interpolation_factors<TF> ifac = timeloop.get_interpolation_factors(time);

    // Calculate the new vertical profile
    for (int k=0; k<gd.kmax; ++k)
        prof[k+kgc] = ifac.fac0 * data[ifac.index0*kk+k] + ifac.fac1 * data[ifac.index1*kk+k];
}

template <typename TF>
void Timedep<TF>::update_time_dependent(TF& val, Timeloop<TF>& timeloop)
{
    if (sw == Timedep_switch::Disabled)
        return;

    // Get/calculate the interpolation indexes/factors
    Interpolation_factors<TF> ifac = timeloop.get_interpolation_factors(time);
    val = ifac.fac0 * data[ifac.index0] + ifac.fac1 * data[ifac.index1];
    return;
}

template class Timedep<double>;
template class Timedep<float>;
