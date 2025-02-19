/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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
#include "constants.h"
#include "timeloop.h"
#include "constants.h"

#include "particle_bin.h"

namespace
{
    template<typename TF>
    void settle_particles(
            TF* const restrict st,
            const TF* const restrict s,
            const TF* const dzhi,
            const TF w_particle,
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
                    st[ijk] -= w_particle * (s[ijk+kstride]-s[ijk])*dzhi[k+1];
                }
    }
}


template<typename TF>
Particle_bin<TF>::Particle_bin(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    sw_particle = inputin.get_item<bool>("particle_bin", "sw_particle", "", false);

    if (sw_particle)
    {
        std::vector<std::string> scalars = inputin.get_list<std::string>("particle_bin", "particle_list", "", std::vector<std::string>());

        // Read gravitational settling velocities.
        for (auto& scalar : scalars)
        {
            w_particle.emplace(scalar, inputin.get_item<TF>("particle_bin", "w_particle", scalar));

            // Raise error if any of the velocities is positive.
            if (w_particle.at(scalar) > 0)
                throw std::runtime_error("Gravitational settling velocities need to be negative!");
        }

        // Constraint on time stepping.
        cfl_max = inputin.get_item<TF>("particle_bin", "cfl_max", "", 1.2);
    }
}

template<typename TF>
Particle_bin<TF>::~Particle_bin()
{
}

template<typename TF>
void Particle_bin<TF>::create(Timeloop<TF>& timeloop)
{
    if (!sw_particle)
        return;

    auto& gd = grid.get_grid_data();

    // Find minimum vertical grid spacing.
    TF dz_min = TF(Constants::dbig);
    for (int k=gd.kstart; k<gd.kend; ++k)
       dz_min = std::min(dz_min, gd.dz[k]);

    // Find maximum gravitational settling velocity.
    TF w_max = -TF(Constants::dbig);
    for (auto& w : w_particle)
        w_max = std::max(w_max, std::abs(w.second));

    // Calculate maximum time step.
    const double dt_max = cfl_max / w_max * dz_min;

    idt_max = convert_to_itime(dt_max);
}

template<typename TF>
unsigned long Particle_bin<TF>::get_time_limit()
{
    if (!sw_particle)
        return Constants::ulhuge;

    return idt_max;
}

#ifndef USECUDA
template<typename TF>
void Particle_bin<TF>::exec(Stats<TF>& stats)
{
    if (!sw_particle)
        return;

    auto& gd = grid.get_grid_data();

    for (auto& w : w_particle)
        settle_particles<TF>(
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
template class Particle_bin<float>;
#else
template class Particle_bin<double>;
#endif
