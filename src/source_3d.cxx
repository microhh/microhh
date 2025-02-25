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
    sw_timedep = inputin.get_item<bool>("source", "sw_timedep", "", false);
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

    for (auto& specie : sourcelist)
        emission[specie] = std::vector<TF>(size);
}


template<typename TF>
void Source_3d<TF>::create(Input& input, Netcdf_handle& input_nc)
{
    auto& gd = grid.get_grid_data();

    // Create Field3_io instance to read the binary files.
    Field3d_io<TF> field3d_io(master, grid);

    auto tmp1 = fields.get_tmp();
    auto tmp2 = fields.get_tmp();

    const int itime = 0;
    int nerror = 0;
    const TF no_offset = 0;

    auto load_3d_field = [&](
            TF* const restrict field, const std::string& name, const int itime)
    {
        char filename[256];
        std::snprintf(filename, 256, "%s.%07d", name.c_str(), itime);
        master.print_message("Loading \"%s\" ... ", filename);

        // Storage arrays don't have vertical ghost cells.
        const int kstart = 0;
        const int kend = this->ktot;

        if (field3d_io.load_field3d(
                field,
                tmp1->fld.data(),
                tmp2->fld.data(),
                filename,
                no_offset,
                kstart,
                kend))
        {
            master.print_message("FAILED\n");
            nerror += 1;
        }
        else
            master.print_message("OK\n");
    };

    for (auto& specie : sourcelist)
        load_3d_field(emission[specie].data(), specie, itime);
}


#ifndef USECUDA
template<typename TF>
void Source_3d<TF>::exec()
{
    auto& gd = grid.get_grid_data();

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
    throw std::runtime_error("Time dependent 3D emissions not (yet) implemented.");
}
#endif


#ifdef FLOAT_SINGLE
template class Source_3d<float>;
#else
template class Source_3d<double>;
#endif
