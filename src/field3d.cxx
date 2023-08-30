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

#include <cstdio>
#include <iostream>
#include <cmath>
#include "master.h"
#include "grid.h"
#include "field3d.h"
#include "defines.h"

template<typename TF>
Field3d<TF>::Field3d(
        Master& masterin, Grid<TF>& gridin,
        const std::string& namein, const std::string& longnamein,
        const std::string& unitin, const std::string& groupin,
        const std::array<int,3>& locin) :
    master(masterin),
    grid(gridin)
{
    name     = namein;
    longname = longnamein;
    unit     = unitin;
    group    = groupin;
    loc      = locin;
}

template<typename TF>
Field3d<TF>::~Field3d()
{
}

template<typename TF>
int Field3d<TF>::init()
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    // Calculate the total field memory size
    const long long field_memory_size = (gd.ncells + 6*gd.ijcells + gd.kcells)*sizeof(TF);

    // Keep track of the total memory in fields
    static long long total_memory_size = 0;
    int nerror = 0;
    try
    {
        total_memory_size += field_memory_size;

        // Allocate all fields belonging to the 3d field
        fld     .resize(gd.ncells);
        fld_bot .resize(gd.ijcells);
        fld_top .resize(gd.ijcells);
        fld_mean.resize(gd.kcells);
        grad_bot.resize(gd.ijcells);
        grad_top.resize(gd.ijcells);
        flux_bot.resize(gd.ijcells);
        flux_top.resize(gd.ijcells);
    }
    catch (std::exception &e)
    {
        std::string msg = "Field " + name + " cannot be allocated, total fields memsize " + std::to_string(total_memory_size) + "is too large";
        master.print_message(msg);
        nerror++;
    }
    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("In Field3d::init");

    // set all values to zero
    for (int n=0; n<gd.ncells; ++n)
        fld[n] = 0.;

    for (int n=0; n<gd.kcells; ++n)
        fld_mean[n] = 0.;

    for (int n=0; n<gd.ijcells; ++n)
    {
        fld_bot [n] = 0.;
        fld_top [n] = 0.;
        grad_bot[n] = 0.;
        grad_top[n] = 0.;
        flux_bot[n] = 0.;
        flux_top[n] = 0.;
    }

    return 0;
}


#ifdef FLOAT_SINGLE
template class Field3d<float>;
#else
template class Field3d<double>;
#endif
