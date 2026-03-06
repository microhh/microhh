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

#include <iostream>

#include "tools.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "timeloop.h"

#include "source_3d.h"
#include "source_3d_kernels.cuh"

namespace s3k = Source_3d_kernels_g;


#ifdef USECUDA
template<typename TF>
void Source_3d<TF>::exec(Thermo<TF>& thermo, Timeloop<TF>& timeloop)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, this->ktot);
    dim3 blockGPU(blocki, blockj, 1);

    for (auto& specie : sourcelist)
        s3k::add_source_tend_g<TF><<<gridGPU, blockGPU>>>(
            fields.st.at(specie)->fld_g,
            emission_g.at(specie),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kstart + this->ktot,
            gd.icells, gd.ijcells);
    cuda_check_error();

    if (sw_heat)
    {
        const TF subdti = TF(1) / timeloop.get_sub_time_step();

        auto tmp = fields.get_tmp_g();
        thermo.get_thermo_field_g(*tmp, "T", false);

        TF* exnref_g = thermo.get_basestate_fld_g("exner");

        // YIKES^3... Create a `thermo.get_temperature_var()` function?
        std::string th_var;
        if (thermo.get_switch() == Thermo_type::Dry)
            th_var = "th";
        else if (thermo.get_switch() == Thermo_type::Moist)
            th_var = "thl";
        else
            throw std::runtime_error("No temperature field found.");

        s3k::add_source_tend_heat_g<TF><<<gridGPU, blockGPU>>>(
            fields.st.at(th_var)->fld_g,
            emission_g.at("te"),
            emission_g.at("qe"),
            tmp->fld_g,
            gd.dz_g,
            exnref_g,
            gd.dx, gd.dy,
            subdti,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kstart + this->ktot,
            gd.icells, gd.ijcells);

        fields.release_tmp_g(tmp);
    }
}


template<typename TF>
void Source_3d<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (!sw_timedep)
        return;

    auto& gd = grid.get_grid_data();

    unsigned long itime = timeloop.get_itime();

    if (itime > iloadtime_next)
    {
        // Update two time states in memory.
        iloadtime_prev = iloadtime_next;
        iloadtime_next = iloadtime_prev + iloadfreq;

        unsigned long iiotimeprec = timeloop.get_iiotimeprec();
        const unsigned long iotime_next = int(iloadtime_next / iiotimeprec);

        auto swap_and_load = [&](const std::string& specie)
        {
            // Swap/update host fields.
            emission_prev.at(specie) = emission_next.at(specie);
            load_emission(emission_next.at(specie), specie, iotime_next);

            // Swap/update device fields.
            const int memsize = gd.ijcells * this->ktot * sizeof(TF);
            cuda_safe_call(cudaMemcpy(emission_prev_g.at(specie), emission_next_g.at(specie), memsize, cudaMemcpyDeviceToDevice));
            cuda_safe_call(cudaMemcpy(emission_next_g.at(specie), emission_next.at(specie).data(), memsize, cudaMemcpyHostToDevice));
        };

        // Swap next -> prev field, and read new time.
        for (auto& specie : sourcelist)
            swap_and_load(specie);

        if (sw_heat)
        {
            swap_and_load("qe");
            swap_and_load("te");
        }
    }

    // Interpolate emissions in time.
    const TF fac1 = TF(itime - iloadtime_prev) / TF(iloadtime_next - iloadtime_prev);
    const TF fac0 = TF(1) - fac1;

    // Emission fields don't have ghost cells.
    const int kstart = 0;
    const int kend = this->ktot;

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, this->ktot);
    dim3 blockGPU(blocki, blockj, 1);

    auto interpolate = [&](const std::string& specie)
    {
        s3k::interpolate_emission_g<TF><<<gridGPU, blockGPU>>>(
            emission_g.at(specie),
            emission_prev_g.at(specie),
            emission_next_g.at(specie),
            fac0,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            kstart, kend,
            gd.icells, gd.ijcells);
        cuda_check_error();
    };

    // Interpolate emissions linearly in time.
    for (auto& specie : sourcelist)
        interpolate(specie);

    if (sw_heat)
    {
        interpolate("qe");
        interpolate("te");
    }
}


template<typename TF>
void Source_3d<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();

    const int ncells = gd.ijcells * this->ktot;
    const int memsize = ncells * sizeof(TF);

    auto add_emission = [&](const std::string& specie)
    {
        emission_g.emplace(specie, cuda_vector<TF>(ncells));

        if (sw_timedep)
        {
            emission_prev_g.emplace(specie, cuda_vector<TF>(ncells));
            emission_next_g.emplace(specie, cuda_vector<TF>(ncells));

            cuda_safe_call(cudaMemcpy(emission_prev_g.at(specie), emission_prev.at(specie).data(), memsize, cudaMemcpyHostToDevice));
            cuda_safe_call(cudaMemcpy(emission_next_g.at(specie), emission_next.at(specie).data(), memsize, cudaMemcpyHostToDevice));
        }
        else
            cuda_safe_call(cudaMemcpy(emission_g.at(specie), emission.at(specie).data(), memsize, cudaMemcpyHostToDevice));
    };

    for (auto& specie : sourcelist)
        add_emission(specie);

    if (sw_heat)
    {
        add_emission("qe");
        add_emission("te");
    }
}
#endif


#ifdef FLOAT_SINGLE
template class Source_3d<float>;
#else
template class Source_3d<double>;
#endif
