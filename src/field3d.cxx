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

#ifndef USECUDA
template<typename TF>
Field3d<TF>::~Field3d()
{
}

template<typename TF>
int Field3d<TF>::init()
{
    // Calculate the total field memory size
    const long long field_memory_size = (grid->ncells + 6*grid->ijcells + grid->kcells)*sizeof(double);

    // Keep track of the total memory in fields
    static long long total_memory_size = 0;
    try
    {
        total_memory_size += field_memory_size;

        // Allocate all fields belonging to the 3d field
        data       .resize(grid->ncells);
        databot    .resize(grid->ijcells);
        datatop    .resize(grid->ijcells);
        datamean   .resize(grid->kcells);
        datagradbot.resize(grid->ijcells);
        datagradtop.resize(grid->ijcells);
        datafluxbot.resize(grid->ijcells);
        datafluxtop.resize(grid->ijcells);
    }
    catch (std::exception &e)
    {
        master.print_error("Field %s cannot be allocated, total fields memsize %lu is too large\n", name.c_str(), total_memory_size);
        throw;
    }

    // set all values to zero
    for (int n=0; n<grid->ncells; ++n)
        data[n] = 0.;

    for (int n=0; n<grid->kcells; ++n)
        datamean[n] = 0.;

    for (int n=0; n<grid->ijcells; ++n)
    {
        databot    [n] = 0.;
        datatop    [n] = 0.;
        datagradbot[n] = 0.;
        datagradtop[n] = 0.;
        datafluxbot[n] = 0.;
        datafluxtop[n] = 0.;
    }

    return 0;
}
#endif
