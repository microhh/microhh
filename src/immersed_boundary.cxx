/*
 * MicroHH
 * Copyright (c) 2011-2019 Chiel van Heerwaarden
 * Copyright (c) 2011-2019 Thijs Heus
 * Copyright (c) 2014-2019 Bart van Stratum
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

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "input.h"

#include "immersed_boundary.h"

namespace
{
    /* Bi-linear interpolation of the 2D IB DEM
     * onto the requested (`x_req`, `y_req`) location */
    template<typename TF>
    TF interp2_dem(const TF x_req, const TF y_req,
            const TF* const restrict x, const TF* const restrict y,
            const TF* const restrict dem,
            const TF dx, const TF dy,
            const int igc, const int jgc, const int icells,
            const int imax, const int jmax,
            const int mpicoordx, const int mpicoordy)
    {
        const int ii = 1;
        const int jj = icells;

        const int i0 = (x_req - TF(0.5)*dx) / dx - mpicoordx*imax + igc;
        const int j0 = (y_req - TF(0.5)*dy) / dy - mpicoordy*jmax + jgc;
        const int ij = i0 + j0*jj;

        const TF f1x = (x_req - x[i0]) / dx;
        const TF f1y = (y_req - y[j0]) / dy;
        const TF f0x = TF(1) - f1x;
        const TF f0y = TF(1) - f1y;

        const TF z = f0y * (f0x * dem[ij   ] + f1x * dem[ij+ii   ]) +
                     f1y * (f0x * dem[ij+jj] + f1x * dem[ij+ii+jj]);

        return z;
    }


    template<typename TF>
    bool is_ghost_cell(std::vector<TF>& dem, std::vector<TF>& z,
            const int i, const int j, const int k, const int icells)
    {
        const int ij = i + j*icells;

        if (z[k] <= dem[ij])
        {
            for (int dk = -1; dk <= 1; ++dk)
                for (int dj = -1; dj <= 1; ++dj)
                    for (int di = -1; di <= 1; ++di)
                    {
                        const int ij2 = (i + di) + (j + dj)*icells;
                        if (z[k + dk] > dem[ij2])
                            return true;
                    }
        }

        return false;
    }


    template<typename TF>
    void calc_ghost_cells(std::vector<int>& ghost_i, std::vector<int>& ghost_j, std::vector<int>& ghost_k,
            std::vector<TF>& dem, std::vector<TF> x, std::vector<TF> y, std::vector<TF> z,
            const int istart, const int jstart, const int kstart,
            const int iend,   const int jend,   const int kend,
            const int icells)
    {
        // 1. Find the IB ghost cells
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                    for (int i=istart; i<iend; ++i)
                        if (is_ghost_cell(dem, z, i, j, k, icells))
                        {
                            ghost_i.push_back(i);
                            ghost_j.push_back(j);
                            ghost_k.push_back(k);
                        }

        // 2. Find nearest location on IB
        //for (int i=0; i<ghost_i.size(); ++i)
        //{
        //
        //
        //}
    }

    void print_statistics(std::vector<int>& ghost_i, std::string name, Master& master)
    {
        int nghost = ghost_i.size();
        master.sum(&nghost, 1);

        if (master.get_mpiid() == 0)
        {
            std::string message = "Found: " + std::to_string(nghost) + " IB ghost cells at the " + name + " location";
            master.print_message(message);
        }
    }

}

template<typename TF>
Immersed_boundary<TF>::Immersed_boundary(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin),
    field3d_io(masterin, gridin), boundary_cyclic(masterin, gridin)
{
    // Read IB switch from namelist, and set internal `sw_ib` switch
    std::string sw_ib_str = inputin.get_item<std::string>("IB", "sw_immersed_boundary", "", "0");

    if (sw_ib_str == "0")
        sw_ib = IB_type::Disabled;
    else if (sw_ib_str == "dem")
        sw_ib = IB_type::DEM;
    else
    {
        std::string error = "\"" + sw_ib_str + "\" is an illegal value for \"sw_ib\"";
        throw std::runtime_error(error);
    }

    if (sw_ib != IB_type::Disabled)
    {
        // Set a minimum of 2 ghost cells in the horizontal
        const int ijgc = 2;
        grid.set_minimum_ghost_cells(ijgc, ijgc, 0);

        // Read additional settings
        n_idw_points = inputin.get_item<int>("IB", "n_idw_points", "");
    }
}

template <typename TF>
Immersed_boundary<TF>::~Immersed_boundary()
{
}

template <typename TF>
void Immersed_boundary<TF>::init(Input& inputin)
{
    auto& gd = grid.get_grid_data();

    if (sw_ib == IB_type::DEM)
        dem.resize(gd.ijcells);
}

template <typename TF>
void Immersed_boundary<TF>::create()
{
    if (sw_ib == IB_type::Disabled)
        return;

    // Init the toolbox classes.
    boundary_cyclic.init();

    // Get grid information
    auto& gd = grid.get_grid_data();

    if (sw_ib == IB_type::DEM)
    {
        // Read the IB height (DEM) map
        char filename[256] = "dem.0000000";
        auto tmp = fields.get_tmp();
        master.print_message("Loading \"%s\" ... ", filename);

        if (field3d_io.load_xy_slice(dem.data(), tmp->fld.data(), filename))
        {
            master.print_message("FAILED\n");
            throw std::runtime_error("Reading input DEM field failed");
        }
        else
        {
            master.print_message("OK\n");
        }

        fields.release_tmp(tmp);
        boundary_cyclic.exec_2d(dem.data());

        // Find ghost cells (grid points inside IB, which have at least one
        // neighbouring grid point outside of IB). Different for each
        // location on staggered grid.
        calc_ghost_cells(ghost_u.i, ghost_u.j, ghost_u.k, dem, gd.xh, gd.y, gd.z,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend, gd.icells);

        calc_ghost_cells(ghost_v.i, ghost_v.j, ghost_v.k, dem, gd.x, gd.yh, gd.z,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend, gd.icells);

        calc_ghost_cells(ghost_w.i, ghost_w.j, ghost_w.k, dem, gd.x, gd.y, gd.zh,
                gd.istart, gd.jstart, gd.kstart+1,
                gd.iend,   gd.jend,   gd.kend, gd.icells);

        // Print some statistics (number of ghost cells)
        print_statistics(ghost_u.i, std::string("u"), master);
        print_statistics(ghost_v.i, std::string("v"), master);
        print_statistics(ghost_w.i, std::string("w"), master);

        if (fields.sp.size() > 0) {
            calc_ghost_cells(ghost_s.i, ghost_s.j, ghost_s.k, dem, gd.x, gd.y, gd.z,
                             gd.istart, gd.jstart, gd.kstart,
                             gd.iend, gd.jend, gd.kend, gd.icells);
            print_statistics(ghost_s.i, std::string("s"), master);
        }


    }
}

template class Immersed_boundary<double>;
template class Immersed_boundary<float>;
