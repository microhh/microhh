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

#include <iostream>
#include <algorithm>

#include "aerosol.h"
#include "timeloop.h"
#include "input.h"
#include "grid.h"
#include "netcdf_interface.h"
#include "stats.h"

namespace
{
    // Kernels...
}

template<typename TF>
Aerosol<TF>::Aerosol(
    Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    // Read `.ini` settings.
    sw_aerosol = inputin.get_item<bool>("aerosol", "swaerosol", "", 0);
    sw_timedep = inputin.get_item<bool>("aerosol", "swtimedep", "", 0);
}

template <typename TF>
Aerosol<TF>::~Aerosol()
{
}

template <typename TF>
void Aerosol<TF>::init()
{
    // Allocate (`.resize`) arrays.
    if (!sw_aerosol)
        return;

    auto& gd = grid.get_grid_data();

    std::cout << "init() aerosols" << std::endl;

    bla.resize(gd.kcells);
}

template <typename TF>
void Aerosol<TF>::create(Input& inputin, Netcdf_handle& input_nc, Stats<TF>& stats)
{
    // Read input from NetCDF and prepare statistics output.
    if (!sw_aerosol)
        return;

    auto& gd = grid.get_grid_data();

    std::cout << "create() aerosols" << std::endl;

    // Read NetCDF input.
    Netcdf_group& group_nc = input_nc.get_group("init");
    group_nc.get_variable(bla, "bla", {0}, {gd.ktot});
    std::rotate(bla.rbegin(), bla.rbegin()+gd.kstart, bla.rend());

    // Prepare statistics.
    const std::string group_name = "default";

    if (sw_timedep)
        stats.add_prof("bla", "Bla bla!", "kg K-1", "z" , group_name);
    else
        stats.add_fixed_prof("bla", "Bla bla!", "kg K-1", "z" , group_name, bla);
}

#ifndef USECUDA
template <typename TF>
void Aerosol<TF>::exec(Thermo<TF>& thermo)
{
    if (!sw_aerosol)
        return;

    auto& gd = grid.get_grid_data();

    std::cout << "exec() aerosols" << std::endl;
}

template <typename TF>
void Aerosol<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (!sw_aerosol)
        return;

    if (sw_timedep)
    {
        std::cout << "update_time_dependent() aerosols" << std::endl;

        throw std::runtime_error("Time dependent aerosols are not (yet) supported!\n");
    }
}
#endif

template<typename TF>
void Aerosol<TF>::exec_stats(Stats<TF>& stats)
{
    if (!sw_aerosol)
        return;

    auto& gd = grid.get_grid_data();

    if (sw_timedep)
    {
        std::cout << "exec_stats() aerosols" << std::endl;

        stats.set_prof("bla", bla);
    }
}

template class Aerosol<double>;
template class Aerosol<float>;
