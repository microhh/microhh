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

#include "tools.h"
#include "data_block.h"
#include "timedep.h"


template<typename TF>
Timedep<TF>::Timedep(Master& masterin, Grid<TF>& gridin, std::string varname, std::bool is_timedep) :
    master(masterin), grid(gridin)
{
    sw = is_timedep;
    varname = varname;
}

template<typename TF>
Timedep<TF>::~Timedep()
{
    #ifdef USECUDA
        cuda_safe_call(cudaFree(data_g));
    #endif
}

template <typename TF>
void Timedep<TF>::create_timedep_prof()
{
    if (sw == false)
        return;

    auto& gd = grid.get_grid_data();
    Data_block data_block(master, varname+".time");

    std::vector<std::string> headers = data_block.get_headers();
    // Sort the times
    std::sort(headers.begin()+1, headers.end());

    for (auto& it : headers)
    {
        time.push_back(std::stod(it));
        data_block.get_vector(tmp, it, gd.kmax, 0, gd.kstart);
        data.insert(data[var].end(),tmp.begin(),tmp.end());
    }
}

template <typename TF>
void Timedep<TF>::create_timedep()
{
    if (sw == false)
        return;

    auto& gd = grid.get_grid_data();
    Data_block data_block(master, varname+".time");

    int length = data_block.get_vector_length("time");
    time.resize(length);
    data.resize(length);
    data_block.get_vector(time, "time", length, 0, 0);
    data_block.get_vector(time, varname, length, 0, 0);
}

template <typename TF>
void Timedep<TF>::update_time_dependent_prof(std::vector<TF> prof, Timeloop<TF>& timeloop)
{
    if (sw == false)
        return;

    auto& gd = grid.get_grid_data();
    const int kk = gd.kmax;
    const int kgc = gd.kgc;

    // Get/calculate the interpolation indexes/factors
    int index0 = 0, index1 = 0;
    TF fac0 = 0., fac1 = 0.;
    timeloop.get_interpolation_factors(index0, index1, fac0, fac1, time);

    // Calculate the new vertical profile
    for (int k=0; k<gd.kmax; ++k)
        prof[k+kgc] = fac0 * data[index0*kk+k] + fac1 * data[index1*kk+k];
}
template <typename TF>
TF Timedep<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (sw == false)
        return;

    // Get/calculate the interpolation indexes/factors
    int index0 = 0, index1 = 0;
    TF fac0 = 0., fac1 = 0.;
    timeloop.get_interpolation_factors(index0, index1, fac0, fac1, time);

    return fac0 * data + fac1 * data[name];
}

template class Timedep<double>;
template class Timedep<float>;
