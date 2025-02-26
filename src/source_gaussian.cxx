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

#include <cmath>
#include <algorithm>

#include "grid.h"
#include "fields.h"
#include "master.h"
#include "netcdf_interface.h"

#include "source.h"
#include "source_gaussian.h"
#include "source_gaussian_kernels.h"

namespace sgk = Source_gaussian_kernels;


template<typename TF>
Source_gaussian<TF>::Source_gaussian(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Source<TF>(masterin, gridin, fieldsin, inputin)
{
    sourcelist = inputin.get_list<std::string>("source", "sourcelist", "");
    source_x0  = inputin.get_list<TF>("source", "source_x0", "");
    source_y0  = inputin.get_list<TF>("source", "source_y0", "");
    source_z0  = inputin.get_list<TF>("source", "source_z0", "");
    sigma_x    = inputin.get_list<TF>("source", "sigma_x",   "");
    sigma_y    = inputin.get_list<TF>("source", "sigma_y",   "");
    sigma_z    = inputin.get_list<TF>("source", "sigma_z",   "");
    strength   = inputin.get_list<TF>("source", "strength",  "");
    line_x     = inputin.get_list<TF>("source", "line_x",    "", std::vector<TF>());
    line_y     = inputin.get_list<TF>("source", "line_y",    "", std::vector<TF>());
    line_z     = inputin.get_list<TF>("source", "line_z",    "", std::vector<TF>());
    
    auto check_and_default = [&](std::vector<TF>& vec)
    {
        if (vec.size() > 0 && vec.size() != source_x0.size())
            throw std::runtime_error("Number of line input values doesn't match other source input values.");
        else if (vec.size() == 0 && source_x0.size() > 0)
        {
            vec.resize(source_x0.size());
            std::fill(vec.begin(), vec.end(), TF(0));
        }
    };
    
    // The `line_` input options are allowed to be empty. Set to zero if size differs from other input options.
    check_and_default(line_x);
    check_and_default(line_y);
    check_and_default(line_z);
    
    // Timedep source location
    swtimedep_location = inputin.get_item<bool>("source", "swtimedep_location", "", false);
    swtimedep_strength = inputin.get_item<bool>("source", "swtimedep_strength", "", false);
    
    // Switch between input in mass or volume ratio.
    // swvmr=true = kmol tracer s-1, swvmr=false = kg tracer s-1
    sw_vmr = inputin.get_list<bool>("source", "swvmr", "");
    
    // Option for (non-time dependent) profiles for vertical distribution emissions.
    sw_emission_profile = inputin.get_item<bool>("source", "sw_profile", "", false);
    
    if (sw_emission_profile)
    {
        // Find number of unique input profiles.
        profile_index = inputin.get_list<int>("source", "profile_index", "");
    
        for (int& i : profile_index)
            unique_profile_indexes.insert(i);
    }
    
    // Limits on options (for now...)
    if (sw_emission_profile && (swtimedep_location || swtimedep_strength))
        throw std::runtime_error("Emission profiles with time dependent location/strength are not (yet) supported!");
    
    if (sw_emission_profile)
    {
        for (int i=0; i<source_x0.size(); ++i)
            if (line_x[i] > 0 || line_y[i] > 0 || line_z[i] > 0)
                throw std::runtime_error("Emission profiles with line emissions are not (yet) supported!");
    }
}


template<typename TF>
Source_gaussian<TF>::~Source_gaussian()
{
}


template<typename TF>
void Source_gaussian<TF>::init()
{
    shape.resize(source_x0.size());
    norm.resize(source_x0.size());
}


template<typename TF>
void Source_gaussian<TF>::create(Input& input, Timeloop<TF>& timeloop, Netcdf_handle& input_nc)
{
    auto& gd = grid.get_grid_data();

    if (sw_emission_profile)
    {
        // Get emission profiles from NetCDF input.
        Netcdf_group& init_group = input_nc.get_group("init");

        for (int i : unique_profile_indexes)
        {
            std::vector<TF> prof = std::vector<TF>(gd.kcells);
            init_group.get_variable(prof, "emission_profile_" + std::to_string(i), {0}, {gd.ktot});
            std::rotate(prof.rbegin(), prof.rbegin() + gd.kstart, prof.rend());
            profile_z.emplace(i, prof);
        }
    }

    for (int n=0; n<source_x0.size(); ++n)
    {
        // Shape of the source in each direction
        shape[n].range_x = sgk::calc_shape(gd.x.data(), source_x0[n], sigma_x[n], line_x[n], gd.istart, gd.iend);
        shape[n].range_y = sgk::calc_shape(gd.y.data(), source_y0[n], sigma_y[n], line_y[n], gd.jstart, gd.jend);

        if (sw_emission_profile)
            shape[n].range_z = sgk::calc_shape_profile(profile_z.at(profile_index[n]).data(), gd.z.data(), gd.kstart, gd.kend);
        else
            shape[n].range_z = sgk::calc_shape(gd.z.data(), source_z0[n], sigma_z[n], line_z[n], gd.kstart, gd.kend);

        if (sw_emission_profile)
            norm[n] = calc_norm(
                    gd.x.data(), source_x0[n], sigma_x[n], line_x[n],
                    gd.y.data(), source_y0[n], sigma_y[n], line_y[n],
                    gd.z.data(), source_z0[n], sigma_z[n], line_z[n],
                    profile_z.at(profile_index[n]).data(),
                    shape[n].range_x, shape[n].range_y, shape[n].range_z,
                    fields.rhoref.data(), sw_vmr[n], true);
        else
            norm[n] = calc_norm(
                    gd.x.data(), source_x0[n], sigma_x[n], line_x[n],
                    gd.y.data(), source_y0[n], sigma_y[n], line_y[n],
                    gd.z.data(), source_z0[n], sigma_z[n], line_z[n],
                    nullptr,
                    shape[n].range_x, shape[n].range_y, shape[n].range_z,
                    fields.rhoref.data(), sw_vmr[n], false);
    }

    // Create timedep
    if (swtimedep_location)
    {
        std::string timedep_dim = "time_source";

        for (int n=0; n<source_x0.size(); ++n)
        {
            std::string name_x = "source_x0_" + std::to_string(n);
            std::string name_y = "source_y0_" + std::to_string(n);
            std::string name_z = "source_z0_" + std::to_string(n);

            tdep_source_x0.emplace(name_x, new Timedep<TF>(master, grid, name_x, true));
            tdep_source_y0.emplace(name_y, new Timedep<TF>(master, grid, name_y, true));
            tdep_source_z0.emplace(name_z, new Timedep<TF>(master, grid, name_z, true));

            tdep_source_x0.at(name_x)->create_timedep(input_nc, timedep_dim);
            tdep_source_y0.at(name_y)->create_timedep(input_nc, timedep_dim);
            tdep_source_z0.at(name_z)->create_timedep(input_nc, timedep_dim);
        }
    }

    if (swtimedep_strength)
    {
        std::string timedep_dim = "time_source";

        for (int n=0; n<source_x0.size(); ++n)
        {
            std::string name_strength = "source_strength_" + std::to_string(n);
            tdep_source_strength.emplace(name_strength, new Timedep<TF>(master, grid, name_strength, true));
            tdep_source_strength.at(name_strength)->create_timedep(input_nc, timedep_dim);
        }
    }
}


#ifndef USECUDA
template<typename TF>
void Source_gaussian<TF>::exec()
{
    auto& gd = grid.get_grid_data();

    for (int n=0; n<sourcelist.size(); ++n)
    {
        if (sw_emission_profile)
            sgk::calc_source<TF, true>(
                    fields.st[sourcelist[n]]->fld.data(),
                    gd.x.data(), gd.y.data(), gd.z.data(),
                    profile_z.at(profile_index[n]).data(),
                    source_x0[n], sigma_x[n], line_x[n],
                    source_y0[n], sigma_y[n], line_y[n],
                    source_z0[n], sigma_z[n], line_z[n],
                    strength[n], norm[n],
                    shape[n].range_x[0], shape[n].range_x[1],
                    shape[n].range_y[0], shape[n].range_y[1],
                    shape[n].range_z[0], shape[n].range_z[1],
                    gd.icells, gd.ijcells);
        else
            sgk::calc_source<TF, false>(
                    fields.st[sourcelist[n]]->fld.data(),
                    gd.x.data(), gd.y.data(), gd.z.data(),
                    nullptr,
                    source_x0[n], sigma_x[n], line_x[n],
                    source_y0[n], sigma_y[n], line_y[n],
                    source_z0[n], sigma_z[n], line_z[n],
                    strength[n], norm[n],
                    shape[n].range_x[0], shape[n].range_x[1],
                    shape[n].range_y[0], shape[n].range_y[1],
                    shape[n].range_z[0], shape[n].range_z[1],
                    gd.icells, gd.ijcells);
    }
}


template<typename TF>
void Source_gaussian<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    auto& gd = grid.get_grid_data();

    if (swtimedep_location)
    {
        // Update source locations, and calculate new norm's
        for (int n=0; n<sourcelist.size(); ++n)
        {
            std::string name_x = "source_x0_" + std::to_string(n);
            std::string name_y = "source_y0_" + std::to_string(n);
            std::string name_z = "source_z0_" + std::to_string(n);

            tdep_source_x0.at(name_x)->update_time_dependent(source_x0[n], timeloop);
            tdep_source_y0.at(name_y)->update_time_dependent(source_y0[n], timeloop);
            tdep_source_z0.at(name_z)->update_time_dependent(source_z0[n], timeloop);

            // Shape of the source in each direction
            shape[n].range_x = sgk::calc_shape(gd.x.data(), source_x0[n], sigma_x[n], line_x[n], gd.istart, gd.iend);
            shape[n].range_y = sgk::calc_shape(gd.y.data(), source_y0[n], sigma_y[n], line_y[n], gd.jstart, gd.jend);
            shape[n].range_z = sgk::calc_shape(gd.z.data(), source_z0[n], sigma_z[n], line_z[n], gd.kstart, gd.kend);

            if (sw_emission_profile)
                throw std::runtime_error("Emission profiles with time dependent location/strength are not (yet) supported!");
            else
                norm[n] = calc_norm(
                        gd.x.data(), source_x0[n], sigma_x[n], line_x[n],
                        gd.y.data(), source_y0[n], sigma_y[n], line_y[n],
                        gd.z.data(), source_z0[n], sigma_z[n], line_z[n],
                        nullptr,
                        shape[n].range_x, shape[n].range_y, shape[n].range_z,
                        fields.rhoref.data(), sw_vmr[n], false);
        }
    }

    if (swtimedep_strength)
    {
        // Update source locations, and calculate new norm's
        for (int n=0; n<sourcelist.size(); ++n)
        {
            std::string name_strength = "source_strength_" + std::to_string(n);
            tdep_source_strength.at(name_strength)->update_time_dependent(strength[n], timeloop);
        }
    }
}
#endif


template<typename TF>
TF Source_gaussian<TF>::calc_norm(
        const TF* const restrict x, const TF x0, const TF sigma_x, const TF line_x,
        const TF* const restrict y, const TF y0, const TF sigma_y, const TF line_y,
        const TF* const restrict z, const TF z0, const TF sigma_z, const TF line_z,
        const TF* const restrict emission_profile,
        std::vector<int> range_x, std::vector<int>range_y, std::vector<int> range_z,
        const TF* const restrict rhoref, const bool sw_vmr, const bool use_profile)
{
    namespace fm = Fast_math;
    auto& gd = grid.get_grid_data();

    TF sum = 0.;
    TF blob_norm = 0.;
    TF scaling;
    const int jj = gd.icells;
    const int kk = gd.ijcells;

    for (int k=gd.kstart; k<gd.kend; ++k)
    {
        if (sw_vmr)
            // Emissions come in [kmol tracers s-1] and are added to grid boxes in [VMR s-1] unit.
            // rhoref [kg m-3] divided by xmair [kg kmol-1] transfers to units [kmol(tracer) / kmol(air) / s].
            scaling = rhoref[k]/Constants::xmair<TF>;
        else
            // Emissions come in [kg tracer s-1]. [kg tracer s-1 / (m3 * kg m-3)] results in
            // emissions in units [kg tracer / kg air / s].
            scaling = rhoref[k];

        for (int j=gd.jstart; j<gd.jend; ++j)
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;

                if (i>=range_x[0] && i<=range_x[1] && j>=range_y[0] && j<=range_y[1] && k>=range_z[0] && k<=range_z[1])
                {
                    if (line_x != 0)
                    {
                        if (x[i] >= x0+line_x)
                            blob_norm = exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                    - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                        else if (x[i]<=x0)
                            blob_norm = exp(
                                    - fm::pow2(x[i]-x0)       /fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                    - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                        else
                            blob_norm = exp(
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                    - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                    }
                    else if (line_y != 0)
                    {
                        if (y[j] >= y0+line_y)
                            blob_norm = exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                    - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                        else if (y[j]<=y0)
                            blob_norm = exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0)       /fm::pow2(sigma_y)
                                    - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                        else
                            blob_norm = exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));

                    }
                    else if (line_z != 0)
                    {
                        if (z[k] >= z0+line_z)
                            blob_norm = exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                    - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                        else if (z[k]<=z0)
                            blob_norm = exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                    - fm::pow2(z[k]-z0)       /fm::pow2(sigma_z));
                        else
                            blob_norm = exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y));
                    }
                    else
                    {
                        if (use_profile)
                            blob_norm = exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y))
                                    * emission_profile[k];
                        else
                            blob_norm = exp(
                                    - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                    - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                    - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                    }
                }
                else
                    blob_norm = TF(0);

                sum += blob_norm*gd.dx*gd.dy*gd.dz[k]*scaling;
            }
    }

    master.sum(&sum, 1);

    return sum;
}


#ifdef FLOAT_SINGLE
template class Source_gaussian<float>;
#else
template class Source_gaussian<double>;
#endif
