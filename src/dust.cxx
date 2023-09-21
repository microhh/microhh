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

#include <stdexcept>
#include <iostream>

#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "dust.h"
#include "constants.h"
#include "timeloop.h"

namespace
{
    template<typename TF>
    void settle_dust(
            TF* const restrict st,
            const TF* const restrict s,
            const TF* const dzhi,
            const TF ws,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        // Simple upwind advection, like in subsidence.
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jstride + k*kstride;
                    st[ijk] -= ws * (s[ijk+kstride]-s[ijk])*dzhi[k+1];
                }
    }
}


template<typename TF>
Dust<TF>::Dust(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    sw_dust = inputin.get_item<bool>("dust", "swdust", "", false);

    if (sw_dust)
    {
        std::vector<std::string> scalars = inputin.get_list<std::string>("dust", "dustlist", "", std::vector<std::string>());

        // Read gravitational settling velocities.
        for (auto& scalar : scalars)
        {
            ws.emplace(scalar, inputin.get_item<TF>("dust", "ws", scalar));

            // Raise error if any of the velocities is positive.
            if (ws.at(scalar) > 0)
                throw std::runtime_error("Gravitational settling velocities need to be negative!");
        }

        // Constraint on time stepping.
        cfl_max = inputin.get_item<TF>("dust", "cfl_max", "", 1.2);
    }
}

template<typename TF>
Dust<TF>::~Dust()
{
}

template<typename TF>
void Dust<TF>::create(const double ifactor)
{
    auto& gd = grid.get_grid_data();

    // Find minimum vertical grid spacing.
    TF dz_min = TF(Constants::dbig);
    for (int k=gd.kstart; k<gd.kend; ++k)
       dz_min = std::min(dz_min, gd.dz[k]);

    // Find maximum gravitational settling velocity.
    TF w_max = -TF(Constants::dbig);
    for (auto& w : ws)
        w_max = std::max(w_max, std::abs(w.second));

    // Calculate maximum time step.
    const double dt_max = cfl_max / w_max * dz_min;
    idt_max = dt_max * ifactor;
}

template<typename TF>
unsigned long Dust<TF>::get_time_limit()
{
    return idt_max;
}

#ifndef USECUDA
template<typename TF>
void Dust<TF>::exec(Stats<TF>& stats)
{
    if (!sw_dust)
        return;

    auto& gd = grid.get_grid_data();

    for (auto& w : ws)
        settle_dust<TF>(
                fields.st.at(w.first)->fld.data(),
                fields.sp.at(w.first)->fld.data(),
                gd.dzhi.data(),
                w.second,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);
}
#endif

#ifdef FLOAT_SINGLE
template class Dust<float>;
#else
template class Dust<double>;
#endif
