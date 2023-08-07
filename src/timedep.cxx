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

#include "netcdf_interface.h"
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
void Timedep<TF>::create_timedep_prof(
        Netcdf_handle& input_nc, const TF offset, const std::string time_dim)
{
    if (sw == Timedep_switch::Disabled)
        return;

    Netcdf_group& group_nc = input_nc.get_group("timedep");

    std::map<std::string, int> dims = group_nc.get_variable_dimensions(varname);
    int time_dim_length = group_nc.get_dimension_size(time_dim);

    time.resize(time_dim_length);
    group_nc.get_variable(time, time_dim, {0}, {time_dim_length});

    auto& gd = grid.get_grid_data();
    data.resize(time_dim_length*gd.ktot);

    group_nc.get_variable(data, varname, {0, 0}, {time_dim_length, gd.ktot});

    // Add offset
    for (int i=0; i<data.size(); ++i)
        data[i] += offset;

    #ifdef USECUDA
    prepare_device();
    #endif
}

template <typename TF>
void Timedep<TF>::create_timedep_background_prof(
        Netcdf_handle& input_nc, const TF offset, const std::string time_dim, const TF z_dim_length)
{
    if (sw == Timedep_switch::Disabled)
        return;

    Netcdf_group& group_nc = input_nc.get_group("timedep");

    std::map<std::string, int> dims = group_nc.get_variable_dimensions(varname);
    int time_dim_length = group_nc.get_dimension_size(time_dim);

    time.resize(time_dim_length);
    group_nc.get_variable(time, time_dim, {0}, {time_dim_length});

//    auto& gd = grid.get_grid_data();
    data.resize(time_dim_length*z_dim_length);

    group_nc.get_variable(data, varname, {0, 0}, {time_dim_length, int(z_dim_length)});

    // Add offset
    for (int i=0; i<data.size(); ++i)
        data[i] += offset;

#ifdef USECUDA
    prepare_device();
#endif
}

template <typename TF>
void Timedep<TF>::create_timedep(Netcdf_handle& input_nc, const std::string time_dim)
{
    if (sw == Timedep_switch::Disabled)
        return;

    Netcdf_group& group_nc = input_nc.get_group("timedep");

    std::map<std::string, int> dims = group_nc.get_variable_dimensions(varname);
    int time_dim_length = group_nc.get_dimension_size(time_dim);

    time.resize(time_dim_length);
    data.resize(time_dim_length);

    group_nc.get_variable(time, time_dim, {0}, {time_dim_length});
    group_nc.get_variable(data, varname,  {0}, {time_dim_length});


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
void Timedep<TF>::update_time_dependent_prof(std::vector<TF>& prof, Timeloop<TF>& timeloop, const TF z_dim_length)
{
    if (sw == Timedep_switch::Disabled)
        return;

    const int kk = int(z_dim_length);

    // Get/calculate the interpolation indexes/factors
    Interpolation_factors<TF> ifac = timeloop.get_interpolation_factors(time);

    // Calculate the new vertical profile
    for (int k=0; k<kk; ++k)
        prof[k] = ifac.fac0 * data[ifac.index0*kk+k] + ifac.fac1 * data[ifac.index1*kk+k];
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
