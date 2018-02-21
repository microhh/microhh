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

#include <cstdio>
#include <iostream>
#include <cmath>
#include "master.h"
#include "grid.h"
#include "field3d.h"
#include "defines.h"

template<typename TF>
Field3d<TF>::Field3d(Master& masterin, Grid<TF>& gridin, std::string namein, std::string longnamein, std::string unitin) :
    master(masterin),
    grid(gridin)
{
    name     = namein;
    longname = longnamein;
    unit     = unitin;
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
        master.print_error("Field %s cannot be allocated, total fields memsize %lu is too large\n", name.c_str(), total_memory_size);
        throw;
    }

    // set all values to zero
    for (int n=0; n<gd.ncells; ++n)
        fld[n] = 0.;

    for (int n=0; n<gd.kcells; ++n)
        fld_mean[n] = 0.;

    for (int n=0; n<gd.ijcells; ++n)
    {
        fld_bot    [n] = 0.;
        fld_top    [n] = 0.;
        grad_bot[n] = 0.;
        grad_top[n] = 0.;
        flux_bot[n] = 0.;
        flux_top[n] = 0.;
    }

    return 0;
}

#ifndef USECUDA
template<typename TF>
void Field3d<TF>::calc_mean_profile(TF * tmp)
{
    const auto& gd = grid.get_grid_data();

    for (int k=0; k<gd.kcells; ++k)
    {
        fld_mean[k] = 0.;
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk  = i + j*gd.icells + k*gd.ijcells;
                fld_mean[k] += fld[ijk];
            }
    }

    master.sum(fld_mean.data(), gd.kcells);

    const double n = gd.itot * gd.jtot;

    for (int k=0; k<gd.kcells; ++k)
        fld_mean[k] /= n;
}

// Calculate the volume weighted total mean
// BvS: for now only at full levels
template<typename TF>
TF Field3d<TF>::calc_mean(TF * tmp)
{
    const auto& gd = grid.get_grid_data();

    TF sum = 0;

    for (int k=gd.kstart; k<gd.kend; ++k)
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk  = i + j*gd.icells + k*gd.ijcells;
                sum += fld[ijk] * gd.dz[k];
            }

    master.sum(&sum, 1);
    const TF mean = sum / (gd.itot * gd.jtot * gd.zsize);

    return mean;
}
#endif

template class Field3d<double>;
template class Field3d<float>;
