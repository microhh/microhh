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
#include "boundary_lateral.h"
#include "netcdf_interface.h"
#include "grid.h"
#include "input.h"
#include "master.h"

namespace
{
}

template<typename TF>
Boundary_lateral<TF>::Boundary_lateral(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
        master(masterin), grid(gridin), fields(fieldsin)
{
    sw_inoutflow = inputin.get_item<bool>("boundary", "swinoutflow", "", false);

    if (sw_inoutflow)
    {
        sw_inoutflow_u = inputin.get_item<bool>("boundary", "swinoutflow_u", "", true);
        sw_inoutflow_v = inputin.get_item<bool>("boundary", "swinoutflow_v", "", true);

        inoutflow_s = inputin.get_list<std::string>("boundary", "inoutflow_s", "", std::vector<std::string>());
    }
}

template <typename TF>
Boundary_lateral<TF>::~Boundary_lateral()
{
}

template <typename TF>
void Boundary_lateral<TF>::init()
{
    if (!sw_inoutflow)
        return;

    auto& gd = grid.get_grid_data();

    if (sw_inoutflow_u)
    {
        lbc_w.emplace("u", std::vector<TF>(gd.kcells*gd.jcells));
        lbc_e.emplace("u", std::vector<TF>(gd.kcells*gd.jcells));
    }

    if (sw_inoutflow_v)
    {
        lbc_s.emplace("v", std::vector<TF>(gd.kcells*gd.icells));
        lbc_n.emplace("v", std::vector<TF>(gd.kcells*gd.icells));
    }

    for (auto& fld : inoutflow_s)
    {
        lbc_w.emplace(fld, std::vector<TF>(gd.kcells*gd.jcells));
        lbc_e.emplace(fld, std::vector<TF>(gd.kcells*gd.jcells));
        lbc_s.emplace(fld, std::vector<TF>(gd.kcells*gd.icells));
        lbc_n.emplace(fld, std::vector<TF>(gd.kcells*gd.icells));
    }
}

template <typename TF>
void Boundary_lateral<TF>::create(Input& inputin, const std::string& sim_name)
{
    if (!sw_inoutflow)
        return;

    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    // Read NetCDF file with boundary data.
    Netcdf_file input_nc = Netcdf_file(master, sim_name + "_lbc_input.nc", Netcdf_mode::Read);

    const int ntime = input_nc.get_dimension_size("time");
    time_in = input_nc.get_variable<TF>("time", {ntime});

    auto copy_boundary = [&](
            std::map<std::string, std::vector<TF>>& map_out,
            std::vector<TF>& fld_in,
            const int istart, const int iend,
            const int igc, const int imax,
            const int istride_in, const int istride_out,
            const int mpicoord, const std::string name)
    {
        std::vector<TF> fld_out = std::vector<TF>(ntime * gd.kcells * istride_out);

        for (int t=0; t<ntime; t++)
            for (int k=gd.kstart; k<gd.kend; k++)
                for (int i=istart; i<iend; i++)
                {
                    const int kk = k-gd.kgc;
                    const int ii = i-igc;

                    const int ikt_out = i + k*istride_out + t*istride_out*gd.kcells;
                    const int ikt_in = (ii+mpicoord*imax) + kk*istride_in + t*istride_in*gd.ktot;

                    fld_out[ikt_out] = fld_in[ikt_in];
                }

        map_out.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(name),
                std::forward_as_tuple(std::move(fld_out)));
    };

    if (sw_inoutflow_u)
    {
        // TODO
    }

    if (sw_inoutflow_v)
    {
        // TODO
    }

    for (auto& fld : inoutflow_s)
    {
        // Boundaries for entire domain (exluding ghost cells).
        std::vector<TF> lbc_w_full = input_nc.get_variable<TF>(fld + "_west", {ntime, gd.ktot, gd.jtot});
        std::vector<TF> lbc_e_full = input_nc.get_variable<TF>(fld + "_east", {ntime, gd.ktot, gd.jtot});
        std::vector<TF> lbc_s_full = input_nc.get_variable<TF>(fld + "_south", {ntime, gd.ktot, gd.itot});
        std::vector<TF> lbc_n_full = input_nc.get_variable<TF>(fld + "_north", {ntime, gd.ktot, gd.itot});

        if (md.mpicoordx == 0)
            copy_boundary(
                lbc_w_in, lbc_w_full,
                gd.jstart, gd.jend, gd.jgc, gd.jmax,
                gd.jtot, gd.jcells, md.mpicoordy, fld);

        if (md.mpicoordx == md.npx-1)
            copy_boundary(
                    lbc_e_in, lbc_e_full,
                    gd.jstart, gd.jend, gd.jgc, gd.jmax,
                    gd.jtot, gd.jcells, md.mpicoordy, fld);

        if (md.mpicoordy == 0)
            copy_boundary(
                    lbc_s_in, lbc_s_full,
                    gd.istart, gd.iend, gd.igc, gd.imax,
                    gd.itot, gd.icells, md.mpicoordx, fld);

        if (md.mpicoordy == md.npy-1)
            copy_boundary(
                    lbc_n_in, lbc_n_full,
                    gd.istart, gd.iend, gd.igc, gd.imax,
                    gd.itot, gd.icells, md.mpicoordx, fld);
    }
}

template <typename TF>
void Boundary_lateral<TF>::set_ghost_cells()
{
    if (!sw_inoutflow)
        return;

    auto& gd = grid.get_grid_data();
}

template <typename TF>
void Boundary_lateral<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (!sw_inoutflow)
        return;
}

template class Boundary_lateral<double>;
template class Boundary_lateral<float>;
