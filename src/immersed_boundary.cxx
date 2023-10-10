/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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
#include <algorithm>

#include <constants.h>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "input.h"
#include "immersed_boundary.h"
#include "fast_math.h"
#include "stats.h"
#include "cross.h"

namespace
{
    template<typename TF>
    TF absolute_distance(
            const TF x1, const TF y1, const TF z1,
            const TF x2, const TF y2, const TF z2)
    {
        return std::pow(Fast_math::pow2(x2-x1) + Fast_math::pow2(y2-y1) + Fast_math::pow2(z2-z1), TF(0.5));
    }

    // Help function for sorting std::vector with Neighbour points
    template<typename TF>
    bool compare_value(const Neighbour<TF>& a, const Neighbour<TF>& b)
    {
        return a.distance < b.distance;
    }

    bool has_ending(const std::string& full_string, const std::string& ending)
    {
        if (full_string.length() >= ending.length())
            return (0 == full_string.compare(full_string.length() - ending.length(), ending.length(), ending));
        else
            return false;
    };

    /* Bi-linear interpolation of the 2D IB DEM
     * onto the requested (`x_goal`, `y_goal`) location */
    template<typename TF>
    TF interp2_dem(
            const TF x_goal, const TF y_goal,
            const std::vector<TF>& x, const std::vector<TF>& y, const std::vector<TF>& dem,
            const TF dx, const TF dy,
            const int icells, const int jcells,
            const int mpi_offset_x, const int mpi_offset_y)
    {
        const int ii = 1;
        const int jj = icells;

        // Indices west and south of `x_goal`, `y_goal`
        int i0 = (x_goal - TF(0.5)*dx) / dx + mpi_offset_x;
        int j0 = (y_goal - TF(0.5)*dy) / dy + mpi_offset_y;

        // Account for interpolation in last ghost cell (east and north)
        if (i0 == icells-1)
            i0 -= 1;
        if (j0 == jcells-1)
            j0 -= 1;

        const int ij = i0 + j0*jj;

        // Bounds check...
        if (i0 < 0 or i0 >= icells-1 or j0 < 0 or j0 >= jcells-1)
        {
            std::string error = "IB dem interpolation out of bounds!";
            throw std::runtime_error(error);
        }

        // Interpolation factors
        const TF f1x = (x_goal - x[i0]) / dx;
        const TF f1y = (y_goal - y[j0]) / dy;
        const TF f0x = TF(1) - f1x;
        const TF f0y = TF(1) - f1y;

        const TF z = f0y * (f0x * dem[ij   ] + f1x * dem[ij+ii   ]) +
                     f1y * (f0x * dem[ij+jj] + f1x * dem[ij+ii+jj]);
        return z;
    }


    template<typename TF>
    bool is_ghost_cell(
            const std::vector<TF>& dem,
            const std::vector<TF>& x, const std::vector<TF>& y, const std::vector<TF>& z,
            const TF dx, const TF dy,
            const int i, const int j, const int k,
            const int icells, const int jcells,
            const int mpi_offset_x, const int mpi_offset_y)
    {
        const TF zdem = interp2_dem(
                x[i], y[j], x, y, dem, dx, dy,
                icells, jcells, mpi_offset_x, mpi_offset_y);

        // Check if grid point is below IB. If so; check if
        // one of the neighbouring grid points is outside.
        if (z[k] <= zdem)
        {
            // OLD METHOD:
            //if (z[k+1] > zdem)
            //    return true;

            //for (int dj = -1; dj <= 1; ++dj)
            //{
            //    const TF zdem = interp2_dem(x[i], y[j+dj], x, y, dem, dx, dy,
            //                                icells, mpi_offset_x, mpi_offset_y);
            //    if (z[k] > zdem)
            //        return true;
            //}

            //for (int di = -1; di <= 1; ++di)
            //{
            //    // Interpolate DEM to account for half-level locations x,y
            //    const TF zdem = interp2_dem(x[i+di], y[j], x, y, dem, dx, dy,
            //                                icells, mpi_offset_x, mpi_offset_y);
            //    if (z[k] > zdem)
            //        return true;
            //}

            //// NEW METHOD
            //for (int dj = -1; dj <= 1; ++dj)
            //{
            //    // Interpolate DEM to account for half-level locations x,y
            //    const TF zdem = interp2_dem(
            //            x[i], y[j+dj], x, y, dem, dx, dy,
            //            icells, jcells, mpi_offset_x, mpi_offset_y);

            //    for (int dk = -1; dk <= 1; ++dk)
            //        if (z[k + dk] > zdem)
            //            return true;
            //}

            //for (int di = -1; di <= 1; ++di)
            //{
            //    // Interpolate DEM to account for half-level locations x,y
            //    const TF zdem = interp2_dem(
            //            x[i+di], y[j], x, y, dem, dx, dy,
            //            icells, jcells, mpi_offset_x, mpi_offset_y);

            //    for (int dk = -1; dk <= 1; ++dk)
            //        if (z[k + dk] > zdem)
            //            return true;
            //}

            // NEW METHOD
            for (int dj = -1; dj <= 1; ++dj)
                for (int di = -1; di <= 1; ++di)
                {
                    // Interpolate DEM to account for half-level locations x,y
                    const TF zdem = interp2_dem(
                            x[i+di], y[j+dj], x, y, dem, dx, dy,
                            icells, jcells, mpi_offset_x, mpi_offset_y);

                    for (int dk = -1; dk <= 1; ++dk)
                        if (z[k + dk] > zdem)
                            return true;
                }
        }

        return false;
    }

    template<typename TF>
    void find_nearest_location_wall(
            TF& xb, TF& yb, TF& zb,
            const std::vector<TF>& x, const std::vector<TF>& y, const std::vector<TF>& dem,
            const TF x0, const TF y0, const TF z0,
            const TF dx, const TF dy,
            const int icells, const int jcells,
            const int mpi_offset_x, const int mpi_offset_y)
    {
        TF d_min = 1e12;
        TF x_min, y_min, z_min;
        const int n = 40;

        for (int ii = -n/2; ii < n/2+1; ++ii)
            for (int jj = -n/2; jj < n/2+1; ++jj)
            {
                const TF xc = x0 + 2 * ii / (double) n * dx;
                const TF yc = y0 + 2 * jj / (double) n * dy;
                const TF zc = interp2_dem(xc, yc, x, y, dem, dx, dy, icells, jcells, mpi_offset_x, mpi_offset_y);
                const TF d  = absolute_distance(x0, y0, z0, xc, yc, zc);

                if (d < d_min)
                {
                    d_min = d;
                    x_min = xc;
                    y_min = yc;
                    z_min = zc;
                }
            }

        xb   = x_min;
        yb   = y_min;
        zb   = z_min;
    }

    template<typename TF>
    void find_interpolation_points(
            std::vector<int>& ip_i, std::vector<int>& ip_j, std::vector<int>& ip_k,
            std::vector<TF>& ip_d, std::vector<TF>& c_idw, 
            const int index, const int n_idw,
            const std::vector<TF>& x, const std::vector<TF>& y, const std::vector<TF>& z, const std::vector<TF>& dem,
            const TF d_lim, const TF dx, const TF dy,
            const int i, const int j, const int k,
            const int kstart, const int icells, const int jcells, const int ijcells,
            const int mpi_offset_x, const int mpi_offset_y)
    {
        // Vectors including all neighbours outside IB
        std::vector<Neighbour<TF>> neighbours;

        // Limit vertical stencil near surface
        const int dk0 = std::max(-2, kstart-k);

        // Find neighbouring grid points outside IB
        for (int dk=dk0; dk<6; ++dk)
            for (int dj=-1; dj<2; ++dj)
                for (int di=-1; di<2; ++di)
                {
                    const TF zd = interp2_dem(x[i+di], y[j+dj], x, y, dem, dx, dy, icells, jcells, mpi_offset_x, mpi_offset_y);

                    // Check if grid point is outside IB
                    if (z[k+dk] > zd)
                    {
                        // Calculate distance to IB
                        //TF xb, yb, zb;
                        //find_nearest_location_wall(
                        //        xb, yb, zb, x, y, dem, x[i+di], y[j+dj], z[k+dk],
                        //        dx, dy, icells, jcells, mpi_offset_x, mpi_offset_y);
                        //const TF dist = absolute_distance(xb, yb, zb, x[i+di], y[j+dj], z[k+dk]);

                        //// Exclude if grid point is too close to the IB
                        //if (dist > d_lim)
                        //{
                            const TF distance = absolute_distance(x[i], y[j], z[k], x[i+di], y[j+dj], z[k+dk]);
                            Neighbour<TF> tmp_neighbour = {i+di, j+dj, k+dk, distance};
                            neighbours.push_back(tmp_neighbour);
                        //}
                    }
                }

        // Sort them on distance
        std::sort(neighbours.begin(), neighbours.end(), compare_value<TF>);

        if (neighbours.size() < n_idw)
        {
           std::cout << "ERROR: only found " << neighbours.size() << " interpolation points @ ";
           std::cout << "i=" << i << ", j=" << j << ", k=" << k << std::endl; 
           throw 1;
        }

        // Save `n_idw` nearest neighbours
        for (int ii=0; ii<n_idw; ++ii)
        {
            const int in = ii + index*n_idw;
            ip_i[in] = neighbours[ii].i;
            ip_j[in] = neighbours[ii].j;
            ip_k[in] = neighbours[ii].k;
            ip_d[in] = neighbours[ii].distance;
        }
    }

    template<typename TF>
    void precalculate_idw(std::vector<TF>& c_idw, std::vector<TF>& c_idw_sum,
                          const std::vector<int>& ip_i, const std::vector<int>& ip_j, const std::vector<int>& ip_k,
                          const std::vector<TF>& x, const std::vector<TF>& y, const std::vector<TF>& z,
                          const TF xi, const TF yi, const TF zi,
                          const TF xb, const TF yb, const TF zb,
                          Boundary_type bc, const int index, const int n_idw)
    {
        // Dirichlet BCs use one interpolation point less, and include
        // the boundary value as an interpolation point
        const int n = (bc == Boundary_type::Dirichlet_type) ? n_idw-1 : n_idw;

        TF dist_max = TF(0);

        // Temp vector for calculations
        std::vector<TF> tmp(n_idw);

        // Calculate distances interpolation points -> image point
        for (int i=0; i<n; ++i)
        {
            const int ii  = i + index * n_idw;
            const int ipi = ip_i[ii];
            const int ipj = ip_j[ii];
            const int ipk = ip_k[ii];

            tmp[i] = absolute_distance(xi, yi, zi, x[ipi], y[ipj], z[ipk]);
            dist_max = std::max(dist_max, tmp[i]);
        }

        // For Dirichlet, add distance image point to IB
        if (bc == Boundary_type::Dirichlet_type)
        {
            tmp[n_idw-1] = std::max(absolute_distance(xi, yi, zi, xb, yb, zb), TF(1e-9));
            dist_max = std::max(dist_max, tmp[n_idw-1]);
        }

        // Calculate interpolation coefficients
        for (int i=0; i<n_idw; ++i)
        {
            const int ii = i + index * n_idw;

            c_idw[ii] = std::pow((dist_max - tmp[i]) / (dist_max * tmp[i]), 0.5) + TF(1e-9);
            c_idw_sum[index] += c_idw[ii];
        }
    }

    template<typename TF>
    void calc_ghost_cells(
            Ghost_cells<TF>& ghost, const std::vector<TF>& dem, 
            const std::vector<TF>& x, const std::vector<TF>& y, const std::vector<TF>& z,
            Boundary_type bc, const TF dx, const TF dy, const std::vector<TF>& dz, 
            const int n_idw,
            const int istart, const int jstart, const int kstart,
            const int iend,   const int jend,   const int kend,
            const int icells, const int jcells, const int ijcells,
            const int mpi_offset_x, const int mpi_offset_y)
    {
        // 1. Find the IB ghost cells
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                    if (is_ghost_cell(dem, x, y, z, dx, dy, i, j, k,
                                      icells, jcells, mpi_offset_x, mpi_offset_y))
                    {
                        ghost.i.push_back(i);
                        ghost.j.push_back(j);
                        ghost.k.push_back(k);
                    }

        const int nghost = ghost.i.size();
        ghost.nghost = nghost;

        // 2. For each ghost cell, find the nearest location on the wall,
        // the image point (ghost cell mirrored across the IB),
        // and distance between image point and ghost cell.
        ghost.xb.resize(nghost);
        ghost.yb.resize(nghost);
        ghost.zb.resize(nghost);

        ghost.xi.resize(nghost);
        ghost.yi.resize(nghost);
        ghost.zi.resize(nghost);

        ghost.di.resize(nghost);

        for (int n=0; n<nghost; ++n)
        {
            // Indices ghost cell in 3D field
            const int i = ghost.i[n];
            const int j = ghost.j[n];
            const int k = ghost.k[n];

            find_nearest_location_wall(
                    ghost.xb[n], ghost.yb[n], ghost.zb[n],
                    x, y, dem, x[i], y[j], z[k],
                    dx, dy, icells, jcells, mpi_offset_x, mpi_offset_y);

            // Image point
            ghost.xi[n] = 2*ghost.xb[n] - x[i];
            ghost.yi[n] = 2*ghost.yb[n] - y[j];
            ghost.zi[n] = 2*ghost.zb[n] - z[k];

            // Distance image point -> ghost cell
            ghost.di[n] = absolute_distance(ghost.xi[n], ghost.yi[n], ghost.zi[n], x[i], y[j], z[k]);
        }

        // 3. Find N interpolation points outside of IB
        ghost.ip_i .resize(nghost*n_idw);
        ghost.ip_j .resize(nghost*n_idw);
        ghost.ip_k .resize(nghost*n_idw);
        ghost.ip_d .resize(nghost*n_idw);
        ghost.c_idw.resize(nghost*n_idw);

        for (int n=0; n<nghost; ++n)
        {
            // Exclude interpolation points closer than `d_lim` to IB
            const TF dist_lim = 0.1 * std::min(std::min(dx, dy), dz[ghost.k[n]]);

            find_interpolation_points(
                    ghost.ip_i, ghost.ip_j, ghost.ip_k, ghost.ip_d, ghost.c_idw,
                    n, n_idw, x, y, z, dem, dist_lim, dx, dy,
                    ghost.i[n], ghost.j[n], ghost.k[n], kstart,
                    icells, jcells, ijcells,
                    mpi_offset_x, mpi_offset_y);
        }

        // 4. Calculate interpolation coefficients
        ghost.c_idw.resize(nghost*n_idw);
        ghost.c_idw_sum.resize(nghost);

        for (int n=0; n<nghost; ++n)
        {
            precalculate_idw(
                    ghost.c_idw, ghost.c_idw_sum,
                    ghost.ip_i, ghost.ip_j, ghost.ip_k, x, y, z,
                    ghost.xi[n], ghost.yi[n], ghost.zi[n],
                    ghost.xb[n], ghost.yb[n], ghost.zb[n],
                    bc, n, n_idw);
        }
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

    template<typename TF>
    void set_ghost_cells(
            TF* const restrict fld, const TF* const restrict boundary_value,
            const TF* const restrict c_idw, const TF* const restrict c_idw_sum,
            const TF* const restrict di,
            const int* const restrict gi, const int* const restrict gj, const int* const restrict gk,
            const int* const restrict ipi, const int* const restrict ipj, const int* const restrict ipk,
            Boundary_type bc, const TF visc, const int n_ghostcells, const int n_idw,
            const int icells, const int ijcells)
    {
        const int n_idw_loc = (bc == Boundary_type::Dirichlet_type) ? n_idw-1 : n_idw;

        for (int n=0; n<n_ghostcells; ++n)
        {
            const int ijkg = gi[n] + gj[n]*icells + gk[n]*ijcells;

            // Sum the IDW coefficient times the value at the neighbouring grid points
            TF vI = TF(0);
            for (int i=0; i<n_idw_loc; ++i)
            {
                const int ii = i + n*n_idw;
                const int ijki = ipi[ii] + ipj[ii]*icells + ipk[ii]*ijcells;
                vI += c_idw[ii] * fld[ijki];
            }

            // For Dirichlet BCs, add the boundary value
            if (bc == Boundary_type::Dirichlet_type)
            {
                const int ii = n_idw-1 + n*n_idw;
                vI += c_idw[ii] * boundary_value[n];
            }

            vI /= c_idw_sum[n];

            // Set the ghost cells, depending on the IB boundary conditions
            if (bc == Boundary_type::Dirichlet_type)
                fld[ijkg] = 2*boundary_value[n] - vI;       // Image value reflected across IB
            else if (bc == Boundary_type::Neumann_type)
                fld[ijkg] = vI - boundary_value[n] * di[n]; // Image value minus gradient times distance
            else if (bc == Boundary_type::Flux_type)
            {
                const TF grad = -boundary_value[n] / visc;
                fld[ijkg] = vI - grad * di[n];              // Image value minus gradient times distance
            }
        }
    }


    template<typename TF>
    void calc_mask(
            TF* const restrict mask, TF* const restrict maskh,
            const TF* const restrict z_dem,
            const TF* const restrict z, const TF* const restrict zh,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ij  = i + j*jj;

                    const int is_not_ib   = z [k] > z_dem[ij];
                    const int is_not_ib_h = zh[k] > z_dem[ij];

                    mask[ijk]  = static_cast<TF>(is_not_ib);
                    maskh[ijk] = static_cast<TF>(is_not_ib_h);
                }
    }


    template<typename TF>
    void find_k_dem(
            unsigned int* const restrict k_dem,
            const TF* const restrict dem,
            const TF* const restrict z,
            const int istart, const int iend, 
            const int jstart, const int jend, 
            const int kstart, const int kend,
            const int jj)
    {
        for (int j=jstart; j<jend; ++j)
            for (int i=istart; i<iend; ++i)
            {
                const int ij = i + j*jj;
                int k = -1;
                for (k=kstart; k<kend; ++k)
                {
                    if (z[k] > dem[ij])
                        break;
                }
                k_dem[ij] = k;
            }
    }


    template<typename TF>
    void calc_fluxes(
            TF* const restrict flux,
            const unsigned int* const restrict k_dem,
            const TF* restrict s,
            const TF dx, const TF dy, const TF* const restrict dz,
            const TF dxi, const TF dyi, const TF* const restrict dzhi,
            const TF svisc,
            const int istart, const int iend, 
            const int jstart, const int jend, 
            const int kstart, const int kend,
            const int jj, const int kk)
    {
        const int ii = 1;

        for (int j=jstart; j<jend; ++j)
            for (int i=istart; i<iend; ++i)
            {
                // Weight all fluxes by the respective area of the face through which they go.
                // Fluxes are only exchanged with neighbors that are a ghost cell. This requires
                // a drawing...

                // Add the vertical flux.
                const int ij  = i + j*jj;
                {
                    const int ijk = i + j*jj + k_dem[ij]*kk;
                    flux[ij] = -svisc*(s[ijk]-s[ijk-kk])*dzhi[k_dem[ij]] * dx*dy;
                }

                // West flux.
                for (int k=k_dem[ij]; k<k_dem[ij-ii]; ++k)
                {
                    const int ijk = i + j*jj + k*kk;
                    flux[ij] += -svisc*(s[ijk]-s[ijk-ii])*dxi * dy*dz[k];
                }
                // East flux.
                for (int k=k_dem[ij]; k<k_dem[ij+ii]; ++k)
                {
                    const int ijk = i + j*jj + k*kk;
                    flux[ij] += svisc*(s[ijk+ii]-s[ijk])*dxi * dy*dz[k];
                }
                // South flux.
                for (int k=k_dem[ij]; k<k_dem[ij-jj]; ++k)
                {
                    const int ijk = i + j*jj + k*kk;
                    flux[ij] += -svisc*(s[ijk]-s[ijk-jj])*dyi * dx*dz[k];
                }
                // North flux.
                for (int k=k_dem[ij]; k<k_dem[ij+jj]; ++k)
                {
                    const int ijk = i + j*jj + k*kk;
                    flux[ij] += svisc*(s[ijk+jj]-s[ijk])*dyi * dx*dz[k];
                }

                // Normalize the fluxes back to the correct units.
                flux[ij] /= dx*dy;
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

        // Set available masks
        available_masks.insert(available_masks.end(), {"ib"});
    }
}

template <typename TF>
Immersed_boundary<TF>::~Immersed_boundary()
{
}

#ifndef USECUDA
template <typename TF>
void Immersed_boundary<TF>::exec_momentum()
{
    if (sw_ib == IB_type::Disabled)
        return;

    auto& gd = grid.get_grid_data();

    set_ghost_cells(
            fields.mp.at("u")->fld.data(), ghost.at("u").mbot.data(),
            ghost.at("u").c_idw.data(), ghost.at("u").c_idw_sum.data(), ghost.at("u").di.data(),
            ghost.at("u").i.data(), ghost.at("u").j.data(), ghost.at("u").k.data(),
            ghost.at("u").ip_i.data(), ghost.at("u").ip_j.data(), ghost.at("u").ip_k.data(),
            Boundary_type::Dirichlet_type, fields.visc, ghost.at("u").i.size(), n_idw_points,
            gd.icells, gd.ijcells);

    set_ghost_cells(
            fields.mp.at("v")->fld.data(), ghost.at("v").mbot.data(),
            ghost.at("v").c_idw.data(), ghost.at("v").c_idw_sum.data(), ghost.at("v").di.data(),
            ghost.at("v").i.data(), ghost.at("v").j.data(), ghost.at("v").k.data(),
            ghost.at("v").ip_i.data(), ghost.at("v").ip_j.data(), ghost.at("v").ip_k.data(),
            Boundary_type::Dirichlet_type, fields.visc, ghost.at("v").i.size(), n_idw_points,
            gd.icells, gd.ijcells);

    set_ghost_cells(
            fields.mp.at("w")->fld.data(), ghost.at("w").mbot.data(),
            ghost.at("w").c_idw.data(), ghost.at("w").c_idw_sum.data(), ghost.at("w").di.data(),
            ghost.at("w").i.data(), ghost.at("w").j.data(), ghost.at("w").k.data(),
            ghost.at("w").ip_i.data(), ghost.at("w").ip_j.data(), ghost.at("w").ip_k.data(),
            Boundary_type::Dirichlet_type, fields.visc, ghost.at("w").i.size(), n_idw_points,
            gd.icells, gd.ijcells);

    boundary_cyclic.exec(fields.mp.at("u")->fld.data());
    boundary_cyclic.exec(fields.mp.at("v")->fld.data());
    boundary_cyclic.exec(fields.mp.at("w")->fld.data());
}

template <typename TF>
void Immersed_boundary<TF>::exec_scalars()
{
    if (sw_ib == IB_type::Disabled)
        return;

    auto& gd = grid.get_grid_data();

    for (auto& it : fields.sp)
    {
        set_ghost_cells(
                it.second->fld.data(), ghost.at("s").sbot.at(it.first).data(),
                ghost.at("s").c_idw.data(), ghost.at("s").c_idw_sum.data(), ghost.at("s").di.data(),
                ghost.at("s").i.data(), ghost.at("s").j.data(), ghost.at("s").k.data(),
                ghost.at("s").ip_i.data(), ghost.at("s").ip_j.data(), ghost.at("s").ip_k.data(),
                sbcbot, it.second->visc, ghost.at("s").i.size(), n_idw_points,
                gd.icells, gd.ijcells);

        boundary_cyclic.exec(it.second->fld.data());
    }
}
#endif

template <typename TF>
void Immersed_boundary<TF>::init(Input& inputin, Cross<TF>& cross)
{
    auto& gd = grid.get_grid_data();

    if (sw_ib == IB_type::Disabled)
        return;
    else if (sw_ib == IB_type::DEM)
    {
        dem  .resize(gd.ijcells);
        k_dem.resize(gd.ijcells);
    }

    // Process the boundary conditions for scalars
    if (fields.sp.size() > 0)
    {
        // All scalars have the same boundary type (for now)
        std::string swbot = inputin.get_item<std::string>("IB", "sbcbot", "");

        if (swbot == "flux")
            sbcbot = Boundary_type::Flux_type;
        else if (swbot == "dirichlet")
            sbcbot = Boundary_type::Dirichlet_type;
        else if (swbot == "neumann")
            sbcbot = Boundary_type::Neumann_type;
        else
        {
            std::string error = "IB sbcbot=" + swbot + " is not a valid choice (options: dirichlet, neumann, flux)";
            throw std::runtime_error(error);
        }

        // Process boundary values per scalar
        for (auto& it : fields.sp)
            sbc.emplace(it.first, inputin.get_item<TF>("IB", "sbot", it.first));

        // Read the scalars with spatial patterns
        sbot_spatial_list = inputin.get_list<std::string>("IB", "sbot_spatial", "", std::vector<std::string>());
    }

    // Check input list of cross variables (crosslist)
    std::vector<std::string>& crosslist_global = cross.get_crosslist();
    std::vector<std::string>::iterator it = crosslist_global.begin();
    while (it != crosslist_global.end())
    {
        const std::string fluxbot_ib_string = "fluxbot_ib";
        if (has_ending(*it, fluxbot_ib_string))
        {
            // Strip the ending.
            std::string scalar = *it;
            scalar.erase(it->length() - fluxbot_ib_string.length());

            // Check if array is exists, else cycle.
            if (fields.sp.find(scalar) != fields.sp.end())
            {
                // Remove variable from global list, put in local list
                crosslist.push_back(*it);
                crosslist_global.erase(it); // erase() returns iterator of next element..
            }
            else
                ++it;
        }
        else
            ++it;
    }
}

template <typename TF>
void Immersed_boundary<TF>::create()
{
    if (sw_ib == IB_type::Disabled)
        return;

    // Init the toolbox classes.
    boundary_cyclic.init();

    // Get grid and MPI information
    auto& gd  = grid.get_grid_data();
    auto& mpi = master.get_MPI_data();

    if (sw_ib == IB_type::DEM)
    {
        // Offsets used in the 2D DEM interpolation
        const int mpi_offset_x = -mpi.mpicoordx * gd.imax + gd.igc;
        const int mpi_offset_y = -mpi.mpicoordy * gd.jmax + gd.jgc;

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
        ghost.emplace("u", Ghost_cells<TF>());
        ghost.emplace("v", Ghost_cells<TF>());
        ghost.emplace("w", Ghost_cells<TF>());

        master.print_message("Calculating ghost cells u\n");
        calc_ghost_cells(
                ghost.at("u"), dem, gd.xh, gd.y, gd.z,
                Boundary_type::Dirichlet_type,
                gd.dx, gd.dy, gd.dz, n_idw_points,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.jcells, gd.ijcells,
                mpi_offset_x, mpi_offset_y);

        master.print_message("Calculating ghost cells v\n");
        calc_ghost_cells(
                ghost.at("v"), dem, gd.x, gd.yh, gd.z,
                Boundary_type::Dirichlet_type,
                gd.dx, gd.dy, gd.dz, n_idw_points,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.jcells, gd.ijcells,
                mpi_offset_x, mpi_offset_y);

        master.print_message("Calculating ghost cells w\n");
        calc_ghost_cells(
                ghost.at("w"), dem, gd.x, gd.y, gd.zh,
                Boundary_type::Dirichlet_type,
                gd.dx, gd.dy, gd.dzh, n_idw_points,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.jcells, gd.ijcells,
                mpi_offset_x, mpi_offset_y);

        // Print some statistics (number of ghost cells)
        print_statistics(ghost.at("u").i, std::string("u"), master);
        print_statistics(ghost.at("v").i, std::string("v"), master);
        print_statistics(ghost.at("w").i, std::string("w"), master);

        // Momentum boundary condition
        ghost.at("u").mbot.resize(ghost.at("u").nghost);
        ghost.at("v").mbot.resize(ghost.at("v").nghost);
        ghost.at("w").mbot.resize(ghost.at("w").nghost);

        std::fill(ghost.at("u").mbot.begin(), ghost.at("u").mbot.begin(), 0.);
        std::fill(ghost.at("v").mbot.begin(), ghost.at("v").mbot.begin(), 0.);
        std::fill(ghost.at("w").mbot.begin(), ghost.at("w").mbot.begin(), 0.);

        if (fields.sp.size() > 0)
        {
            ghost.emplace("s", Ghost_cells<TF>());

            master.print_message("Calculating ghost cells s\n");
            calc_ghost_cells(
                    ghost.at("s"), dem, gd.x, gd.y, gd.z, sbcbot,
                    gd.dx, gd.dy, gd.dz, n_idw_points,
                    gd.istart, gd.jstart, gd.kstart,
                    gd.iend,   gd.jend,   gd.kend,
                    gd.icells, gd.jcells, gd.ijcells,
                    mpi_offset_x, mpi_offset_y);

            print_statistics(ghost.at("s").i, std::string("s"), master);

            // Read spatially varying boundary conditions (if necessary)
            for (auto& scalar : fields.sp)
            {
                ghost.at("s").sbot.emplace(scalar.first, std::vector<TF>(ghost.at("s").nghost));

                if (std::find(sbot_spatial_list.begin(), sbot_spatial_list.end(), scalar.first) != sbot_spatial_list.end())
                {
                    // Read 2D sbot into tmp field
                    auto tmp = fields.get_tmp();

                    std::string sbot_file = scalar.first + "_sbot.0000000";
                    master.print_message("Loading \"%s\" ... ", sbot_file.c_str());

                    if (field3d_io.load_xy_slice(tmp->fld_bot.data(), tmp->fld.data(), sbot_file.c_str()))
                    {
                        master.print_message("FAILED\n");
                        throw std::runtime_error("Reading input sbot field failed");
                    }
                    else
                        master.print_message("OK\n");

                    // Interpolate 2D sbot onto the ghost cell boundary locations

                    for (int i=0; i<ghost.at("s").nghost; ++i)
                    {
                        ghost.at("s").sbot.at(scalar.first)[i] =
                            interp2_dem(ghost.at("s").xb[i], ghost.at("s").yb[i],
                                   gd.x, gd.y, tmp->fld_bot, gd.dx, gd.dy,
                                   gd.icells, gd.jcells, mpi_offset_x, mpi_offset_y);
                    }
                }
                else
                {
                    for (int i=0; i<ghost.at("s").nghost; ++i)
                        ghost.at("s").sbot.at(scalar.first)[i] = sbc.at(scalar.first);
                }
            }

        }

        // Create the array with vertical indices that give the first cell above the DEM. 
        find_k_dem(
                k_dem.data(), dem.data(), gd.z.data(),
                gd.istart, gd.iend, 
                gd.jstart, gd.jend, 
                gd.kstart, gd.kend,
                gd.icells);

        boundary_cyclic.exec_2d(k_dem.data());
    }
}

template<typename TF>
bool Immersed_boundary<TF>::has_mask(std::string mask_name)
{
    if (std::find(available_masks.begin(), available_masks.end(), mask_name) != available_masks.end())
        return true;
    else
        return false;
}

template<typename TF>
void Immersed_boundary<TF>::get_mask(Stats<TF>& stats, std::string mask_name)
{
    auto& gd = grid.get_grid_data();

    auto mask  = fields.get_tmp();
    auto maskh = fields.get_tmp();

    calc_mask(
            mask->fld.data(), maskh->fld.data(), dem.data(), gd.z.data(), gd.zh.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    stats.set_mask_thres("ib", *mask, *maskh, TF(0.5), Stats_mask_type::Plus);

    fields.release_tmp(mask );
    fields.release_tmp(maskh);
}


template<typename TF>
void Immersed_boundary<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    auto& gd = grid.get_grid_data();
    TF no_offset = 0.;

    if (cross.get_switch())
    {
        for (auto& s : crosslist)
        {
            const std::string fluxbot_ib_string = "fluxbot_ib";
            if (has_ending(s, fluxbot_ib_string))
            {
                // Strip the scalar from the fluxbot_ib
                std::string scalar = s;
                scalar.erase(s.length() - fluxbot_ib_string.length());

                auto tmp = fields.get_tmp();

                calc_fluxes(
                        tmp->flux_bot.data(), k_dem.data(),
                        fields.sp.at(scalar)->fld.data(),
                        gd.dx,  gd.dy,  gd.dz.data(),
                        gd.dxi, gd.dyi, gd.dzhi.data(),
                        fields.sp.at(scalar)->visc,
                        gd.istart, gd.iend, 
                        gd.jstart, gd.jend, 
                        gd.kstart, gd.kend,
                        gd.icells, gd.ijcells);

                cross.cross_plane(tmp->flux_bot.data(), no_offset, scalar+"_fluxbot_ib", iotime);

                fields.release_tmp(tmp);
            }
        }
    }
}


#ifdef FLOAT_SINGLE
template class Immersed_boundary<float>;
#else
template class Immersed_boundary<double>;
#endif
