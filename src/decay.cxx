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

//#include <cstdio>
//#include <algorithm>
#include <iostream>
//#include <math.h>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "decay.h"

namespace
{
    template<typename TF>
    void enforce_exponential_decay(TF* restrict tend, const TF* var, const TF decaytime, const TF dt, const int istart, const int iend, const int jstart, const int jend,
                            const int kstart, const int kend, const int jj, const int kk)
    {
        const TF rate = 1./(std::max(decaytime, dt));

        for (int k=kstart; k<kend; ++k)
           for (int j=jstart; j<jend; ++j)
               #pragma ivdep
               for (int i=istart; i<iend; ++i)
               {
                   const int ijk = i + j*jj + k*kk;
                   tend[ijk] -= rate * var[ijk];
               }
    }
}

template<typename TF>
Decay<TF>::Decay(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    // Read the switches from the input
    std::string swdecay_in = inputin.get_item<std::string>("decay", "swdecay", "", "0");

    // Set the internal switches and read other required input
    if (swdecay_in == "0")
        swdecay = Decay_type::disabled;
    else if (swdecay_in == "exponential")
        swdecay = Decay_type::enabled;
    else
        throw std::runtime_error("Invalid option for \"swdecay\"");
}

template <typename TF>
Decay<TF>::~Decay()
{
}

template <typename TF>
void Decay<TF>::init(Input& inputin)
{
    if (swdecay == Decay_type::enabled)
    {
        std::string type;
        for (auto& it : fields.st)
        {
            type = inputin.get_item<std::string>("decay", "type", it.first,"0");
            if (type == "0")
            {
            }
            else if (type == "exponential")
            {
                dmap[it.first].type = Decay_type::exponential;
                dmap[it.first].timescale = inputin.get_item<TF>("decay", "timescale", it.first);
            }
            else
                throw std::runtime_error("Invalid option for \"decay type\"");
        }

    }

}

template <typename TF>
void Decay<TF>::create(Input& inputin)
{
}

//#ifndef USECUDA
template <typename TF>
void Decay<TF>::exec(double dt)
{
    auto& gd = grid.get_grid_data();

    if (swdecay == Decay_type::enabled)
    {
        for (auto& it : dmap)
        {
            if (it.second.type == Decay_type::exponential)
            {
                enforce_exponential_decay<TF>(fields.st.at(it.first)->fld.data(), fields.sp.at(it.first)->fld.data(), it.second.timescale, dt, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
            }
        }

    }

}
//#endif


template class Decay<double>;
template class Decay<float>;
