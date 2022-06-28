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
#include "fields.h"
#include "defines.h"
#include "field3d_operators.h"


template<typename TF>
Field3d_operators<TF>::Field3d_operators(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
}

template<typename TF>
Field3d_operators<TF>::~Field3d_operators()
{
}

template<typename TF>
void Field3d_operators<TF>::calc_mean_profile(TF* const restrict prof, const TF* const restrict fld)
{
    const auto& gd = grid.get_grid_data();
    const double n = gd.itot * gd.jtot;

    #pragma omp parallel for
    for (int k=0; k<gd.kcells; ++k)
    {
        double tmp = 0.;
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk  = i + j*gd.icells + k*gd.ijcells;
                tmp += fld[ijk];
            }
        prof[k] = tmp / n;
    }

    master.sum(prof, gd.kcells);
}

template<typename TF>
TF Field3d_operators<TF>::calc_mean_2d(const TF* const restrict fld)
{
    const auto& gd = grid.get_grid_data();
    const double n = gd.itot * gd.jtot;

    double value = 0.;

    #pragma omp parallel for
    for (int j=gd.jstart; j<gd.jend; ++j)
        #pragma ivdep
        for (int i=gd.istart; i<gd.iend; ++i)
        {
            const int ij  = i + j*gd.icells;
            value += fld[ij];
        }

    value /= n;
    master.sum(&value, 1);

    return value;
}

template<typename TF>
void Field3d_operators<TF>::calc_mean_profile_nogc(
        TF* const restrict prof, const TF* const restrict fld, bool is_hlf)
{
    const auto& gd = grid.get_grid_data();
    const double n = gd.itot*gd.jtot;

    #pragma omp parallel for
    for (int k=0; k<gd.ktot+is_hlf; ++k)
    {
        double tmp = 0.;
        for (int j=0; j<gd.jmax; ++j)
            #pragma ivdep
            for (int i=0; i<gd.imax; ++i)
            {
                const int ijk  = i + j*gd.imax + k*gd.imax*gd.jmax;
                tmp += fld[ijk];
            }
        prof[k] = tmp / n;
    }
    master.sum(prof, gd.ktot+is_hlf);
}

template<typename TF>
void Field3d_operators<TF>::subtract_mean_profile(TF* const restrict fld, const TF* const restrict prof)
{
    const auto& gd = grid.get_grid_data();

    #pragma omp parallel for
    for (int k=0; k<gd.kcells; ++k)
    {
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk  = i + j*gd.icells + k*gd.ijcells;
                fld[ijk] -= prof[k];
            }
    }
}

// Calculate the volume weighted total mean
// BvS: for now only at full levels
template<typename TF>
TF Field3d_operators<TF>::calc_mean(const TF* const restrict fld)
{
    const auto& gd = grid.get_grid_data();

    double sum = 0;

    #pragma omp parallel for
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

template class Field3d_operators<double>;
template class Field3d_operators<float>;
