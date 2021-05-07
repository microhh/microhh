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
#include <cmath>
#include <vector>
#include <array>
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "source.h"
#include "defines.h"
#include "fast_math.h"

namespace
{
    template<typename TF>
    std::vector<int> calc_shape(
            const TF* restrict x, const TF x0, const TF sigma_x, const TF line_x, int istart, int iend)
    {
        std::vector<int> range(2);

        int i = istart;
        range[0] = iend;

        for (; i<iend; ++i)
        {
            if ( x[i]-x0 + TF(4)*sigma_x > TF(0) )
            {
                range[0] = i;
                break;
            }
        }

        i = istart;
        for (; i<iend; ++i)
        {
            range[1] = iend;

            if ( x[i]-x0-line_x - TF(4)*sigma_x > TF(0) )
            {
                range[1] = i;
                break;
            }
        }

        return range;
    }
}

// Constructor: read values from ini file that do not need info from other classes
template<typename TF>
Source<TF>::Source(Master& master, Grid<TF>& grid, Fields<TF>& fields, Input& input) :
    master(master), grid(grid), fields(fields)
{
    swsource = input.get_item<std::string>("source", "swsource", "", "0");

    if (swsource == "1")
    {
        sourcelist = input.get_list<std::string>("source", "sourcelist", "");

        source_x0 = input.get_list<TF>("source", "source_x0", "");
        source_y0 = input.get_list<TF>("source", "source_y0", "");
        source_z0 = input.get_list<TF>("source", "source_z0", "");
        sigma_x   = input.get_list<TF>("source", "sigma_x"  , "");
        sigma_y   = input.get_list<TF>("source", "sigma_y"  , "");
        sigma_z   = input.get_list<TF>("source", "sigma_z"  , "");
        strength  = input.get_list<TF>("source", "strength" , "");
        line_x    = input.get_list<TF>("source", "line_x"   , "");
        line_y    = input.get_list<TF>("source", "line_y"   , "");
        line_z    = input.get_list<TF>("source", "line_z"   , "");
    }
}

template<typename TF>
Source<TF>::~Source()
{
}

// Init function: allocate memory that you need
template<typename TF>
void Source<TF>::init()
{
    if (swsource == "1")
    {
        shape.resize(source_x0.size());
        norm.resize(source_x0.size());
    }
}

// Create function: read information from ini file that does need info from other class.
template<typename TF>
void Source<TF>::create(Input& input)
{
    auto& gd = grid.get_grid_data();

    for (int n=0; n<source_x0.size(); ++n)
    {
        // Shape of the source in each direction
        shape[n].range_x = calc_shape(gd.x.data(), source_x0[n], sigma_x[n], line_x[n], gd.istart, gd.iend);
        shape[n].range_y = calc_shape(gd.y.data(), source_y0[n], sigma_y[n], line_y[n], gd.jstart, gd.jend);
        shape[n].range_z = calc_shape(gd.z.data(), source_z0[n], sigma_z[n], line_z[n], gd.kstart, gd.kend);

        norm[n] = calc_norm(
                gd.x.data(), source_x0[n], sigma_x[n], line_x[n],
                gd.y.data(), source_y0[n], sigma_y[n], line_y[n],
                gd.z.data(), source_z0[n], sigma_z[n], line_z[n],
                shape[n].range_x, shape[n].range_y, shape[n].range_z);
    }
}

template<typename TF>
TF Source<TF>::calc_norm(
        const TF* const restrict x, const TF x0, const TF sigma_x, const TF line_x,
        const TF* const restrict y, const TF y0, const TF sigma_y, const TF line_y,
        const TF* const restrict z, const TF z0, const TF sigma_z, const TF line_z,
        std::vector<int> range_x, std::vector<int>range_y, std::vector<int> range_z)
{
    namespace fm = Fast_math;
    auto& gd = grid.get_grid_data();

    TF sum = 0.;
    TF blob_norm = 0.;
    const int jj = gd.icells;
    const int kk = gd.ijcells;

    for (int k=gd.kstart; k<gd.kend; ++k)
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
                        blob_norm = exp(
                                - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                }
                else
                    blob_norm = TF(0);

                sum += blob_norm*gd.dx*gd.dy*gd.dz[k];
            }

    master.sum(&sum, 1);

    return sum;
}

// Add the source to the fields. This function is called in the main time loop.
#ifndef USECUDA
template<typename TF>
void Source<TF>::exec()
{
    auto& gd = grid.get_grid_data();

    for (int n=0; n<sourcelist.size(); ++n)
        calc_source(
                fields.st[sourcelist[n]]->fld.data(),
                gd.x.data(), source_x0[n], sigma_x[n], line_x[n],
                gd.y.data(), source_y0[n], sigma_y[n], line_y[n],
                gd.z.data(), source_z0[n], sigma_z[n], line_z[n],
                shape[n].range_x, shape[n].range_y, shape[n].range_z,
                strength[n], norm[n]);
}

template<typename TF>
void Source<TF>::calc_source(
        TF* const restrict st,
        const TF* const restrict x, const TF x0, const TF sigma_x, const TF line_x,
        const TF* const restrict y, const TF y0, const TF sigma_y, const TF line_y,
        const TF* const restrict z, const TF z0, const TF sigma_z, const TF line_z,
        std::vector<int> range_x, std::vector<int> range_y, std::vector<int> range_z,
        const TF strength, TF norm)
{
    namespace fm = Fast_math;
    auto& gd = grid.get_grid_data();

    // Return if one of the blob dimensions is zero:
    if (range_x[0] == range_x[1] || range_y[0] == range_y[1] || range_z[0] == range_z[1])
        return;

    for(int k = range_z[0]; k<range_z[1]; ++k)
        for(int j = range_y[0]; j<range_y[1]; ++j)
            for(int i = range_x[0]; i<range_x[1]; ++i)
            {
                const int ijk = i + j*gd.icells + k*gd.ijcells;

                if (line_x != 0)
                {
                    if (x[i] >= x0+line_x)
                        st[ijk] += strength/norm*exp(
                                - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                    else if (x[i] <= x0)
                        st[ijk] += strength/norm*exp(
                                - fm::pow2(x[i]-x0)       /fm::pow2(sigma_x)
                                - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                    else
                        st[ijk] += strength/norm*exp(
                                - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                }
                else if (line_y != 0)
                {
                    if (y[j] >= y0+line_y)
                        st[ijk] += strength/norm*exp(
                                - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                    else if (y[j] <= y0)
                        st[ijk] += strength/norm*exp(
                                - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                - fm::pow2(y[j]-y0)       /fm::pow2(sigma_y)
                                - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                    else
                        st[ijk] += strength/norm*exp(
                                - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                }
                else if (line_z != 0)
                {
                    if (z[k] >= z0+line_z)
                        st[ijk] += strength/norm*exp(
                                - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
                    else if (z[k] <= z0)
                        st[ijk] += strength/norm*exp(
                                - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                                - fm::pow2(z[k]-z0)       /fm::pow2(sigma_z));
                    else
                        st[ijk] += strength/norm*exp(
                                - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                                - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y));
                }
                else
                    st[ijk] += strength/norm*exp(
                            - fm::pow2(x[i]-x0-line_x)/fm::pow2(sigma_x)
                            - fm::pow2(y[j]-y0-line_y)/fm::pow2(sigma_y)
                            - fm::pow2(z[k]-z0-line_z)/fm::pow2(sigma_z));
            }
}
#endif

template class Source<double>;
template class Source<float>;

