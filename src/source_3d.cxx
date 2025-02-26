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

#include <cmath>
#include <algorithm>

#include "grid.h"
#include "fields.h"
#include "master.h"
#include "netcdf_interface.h"
#include "timeloop.h"

#include "source.h"
#include "source_3d.h"
#include "source_3d_kernels.h"

namespace s3k = Source_3d_kernels;

template<typename TF>
Source_3d<TF>::Source_3d(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Source<TF>(masterin, gridin, fieldsin, inputin)
{
    sourcelist = inputin.get_list<std::string>("source", "sourcelist", "");
    ktot = inputin.get_item<int>("source", "ktot", "");
    sw_timedep = inputin.get_item<bool>("source", "swtimedep", "", false);

    if (sw_timedep)
    {
        const int loadtime = inputin.get_item<int>("source", "loadtime", "");
        iloadfreq = convert_to_itime(loadtime);
    }
}


template<typename TF>
Source_3d<TF>::~Source_3d()
{
}


template<typename TF>
void Source_3d<TF>::init()
{
    auto& gd = grid.get_grid_data();

    // Size without vertical ghost cells.
    const int size = gd.ijcells * this->ktot;

    if (sw_timedep)
    {
        for (auto& specie : sourcelist)
        {
            emission_prev.emplace(specie, std::vector<TF>(size));
            emission_next.emplace(specie, std::vector<TF>(size));
        }
    }

    for (auto& specie : sourcelist)
        emission.emplace(specie, std::vector<TF>(size));
}


template<typename TF>
void Source_3d<TF>::create(Input& input, Timeloop<TF>& timeloop, Netcdf_handle& input_nc)
{
    const int itime = 0;

    if (sw_timedep)
    {
        // Read emissions around current time, for interpolation.
        const unsigned long itime = timeloop.get_itime();
        const unsigned long iiotimeprec = timeloop.get_iiotimeprec();

        // Previous and next load times.
        iloadtime_prev = itime/iloadfreq * iloadfreq;
        iloadtime_next = iloadtime_prev + iloadfreq;

        // IO time accounting for iotimeprec.
        const unsigned long iotime_prev = int(iloadtime_prev / iiotimeprec);
        const unsigned long iotime_next = int(iloadtime_next / iiotimeprec);

        // Read previous and next emissions.
        for (auto& specie : sourcelist)
        {
            load_emission(emission_prev.at(specie), specie, iotime_prev);
            load_emission(emission_next.at(specie), specie, iotime_next);
        }
    }
    else
    {
        // Read emissions which are constant in time.
        for (auto& specie : sourcelist)
            load_emission(emission.at(specie), specie, itime);
    }
}


#ifndef USECUDA
template<typename TF>
void Source_3d<TF>::exec()
{
    auto& gd = grid.get_grid_data();

    // Add emission to scalar tendencies.
    for (auto& specie : sourcelist)
        s3k::add_source_tend(
            fields.st.at(specie)->fld.data(),
            emission.at(specie).data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kstart + this->ktot,
            gd.icells, gd.ijcells);
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

        // Swap next -> prev field, and read new time.
        for (auto& specie : sourcelist)
        {
            emission_prev.at(specie) = emission_next.at(specie);
            load_emission(emission_next.at(specie), specie, iotime_next);
        }
    }

    // Interpolate emissions in time.
    const TF fac1 = TF(itime - iloadtime_prev) / TF(iloadtime_next - iloadtime_prev);
    const TF fac0 = TF(1) - fac1;

    // Emission fields don't have ghost cells.
    const int kstart = 0;
    const int kend = this->ktot;

    // Interpolate emissions linearly in time.
    for (auto& specie : sourcelist)
        s3k::interpolate_emission(
            emission.at(specie).data(),
            emission_prev.at(specie).data(),
            emission_next.at(specie).data(),
            fac0,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            kstart, kend,
            gd.icells, gd.ijcells);
}
#endif


template<typename TF>
void Source_3d<TF>::load_emission(std::vector<TF>& fld, const std::string& name, const int itime)
{
    auto& gd = grid.get_grid_data();

    // Create Field3_io instance to read the binary file.
    Field3d_io<TF> field3d_io(master, grid);

    auto tmp1 = fields.get_tmp();
    auto tmp2 = fields.get_tmp();

    char filename[256];
    std::snprintf(filename, 256, "%s_emission.%07d", name.c_str(), itime);
    master.print_message("Loading \"%s\" ... ", filename);

    // Storage arrays don't have vertical ghost cells.
    const int kstart = 0;
    const int kend = this->ktot;

    const TF no_offset = 0;

    if (field3d_io.load_field3d(
            fld.data(),
            tmp1->fld.data(),
            tmp2->fld.data(),
            filename,
            no_offset,
            kstart,
            kend))
    {
        master.print_message("FAILED\n");
        throw std::runtime_error("Reading binary failed.");
    }
    else
        master.print_message("OK\n");

    fields.release_tmp(tmp1);
    fields.release_tmp(tmp2);
}


#ifdef FLOAT_SINGLE
template class Source_3d<float>;
#else
template class Source_3d<double>;
#endif
