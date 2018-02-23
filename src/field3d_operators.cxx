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
#include "fields.h"
#include "defines.h"
#include "field3d_operators.h"


template<typename TF>
Field3d_operators<TF>::Field3d_operators(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin) :
    master(masterin),
    grid(gridin),
    fields(fieldsin)
{
}

template<typename TF>
Field3d_operators<TF>::~Field3d_operators()
{
}

#ifndef USECUDA
template<typename TF>
void Field3d_operators<TF>::calc_mean_profile(Field3d<TF>* fld)
{
    const auto& gd = grid.get_grid_data();

    for (int k=0; k<gd.kcells; ++k)
    {
        fld->fld_mean[k] = 0.;
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk  = i + j*gd.icells + k*gd.ijcells;
                fld->fld_mean[k] += fld->fld[ijk];
            }
    }

    master.sum(fld->fld_mean.data(), gd.kcells);

    const double n = gd.itot * gd.jtot;

    for (int k=0; k<gd.kcells; ++k)
        fld->fld_mean[k] /= n;
}

// Calculate the volume weighted total mean
// BvS: for now only at full levels
template<typename TF>
TF Field3d_operators<TF>::calc_mean(Field3d<TF>* fld)
{
    const auto& gd = grid.get_grid_data();

    TF sum = 0;

    for (int k=gd.kstart; k<gd.kend; ++k)
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk  = i + j*gd.icells + k*gd.ijcells;
                sum += fld->fld[ijk] * gd.dz[k];
            }

    master.sum(&sum, 1);
    const TF mean = sum / (gd.itot * gd.jtot * gd.zsize);

    return mean;
}
#endif

template class Field3d_operators<double>;
template class Field3d_operators<float>;
