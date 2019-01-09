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

#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>
#include "master.h"
#include "model.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "stats.h"
#include "cross.h"
#include "constants.h"
#include "finite_difference.h"
#include "immersed_boundary.h"
#include "input.h"
#include "diff.h"

namespace
{
    // CvH TEMPORARY UGLINESS NEEDS PERMANENT SOLUTION
    std::vector<double> sbot_ib;
    // CvH END

    // Overloaded method to calculate the absolute distance in 3D
    inline double abs_distance(const double x1, const double x2, const double y1, const double y2, const double z1, const double z2)
    {
        return std::pow(std::pow(x1-x2, 2) + std::pow(y1-y2, 2) + std::pow(z1-z2, 2), 0.5);
    }

    // Overloaded method to calculate the absolute distance in 2D
    inline double abs_distance(const double x1, const double x2, const double z1, const double z2)
    {
        return std::pow(std::pow(x1-x2, 2) + std::pow(z1-z2, 2), 0.5);
    }

    // Help function for sorting std::vector with Neighbour points
    bool compare_value(const Neighbour& a, const Neighbour& b)
    {
        return a.distance < b.distance;
    }

    //  Pre-calculate the coefficients used in the inverse distance weighted interpolation
    void precalculate_idw(Ghost_cell& ghost_cell, const double* const restrict x, const double* const restrict y, const double* const restrict z, Immersed_boundary::Boundary_type bc)
    {
        const int n = ghost_cell.neighbours.size();

        // Calculate distances between image and interpolation points
        for (std::vector<Neighbour>::const_iterator it=ghost_cell.neighbours.begin(); it<ghost_cell.neighbours.end(); ++it)
        {
            const double dist = std::max(abs_distance(ghost_cell.xI, x[it->i],
                                                      ghost_cell.yI, y[it->j],
                                                      ghost_cell.zI, z[it->k]), Constants::dsmall);
            ghost_cell.c_idw.push_back(dist);
        }

        // For dirichlet BCs, add the distance between boundary and image point
        if (bc == Immersed_boundary::Dirichlet_type)
        {
            const double dist = std::max(abs_distance(ghost_cell.xI, ghost_cell.xB,
                                                      ghost_cell.yI, ghost_cell.yB,
                                                      ghost_cell.zI, ghost_cell.zB), Constants::dsmall);

            ghost_cell.c_idw.push_back(dist);
        }

        // Save maximum distance
        const double max_distance = *std::max_element(ghost_cell.c_idw.begin(), ghost_cell.c_idw.end());

        // Calculate IDW coefficients
        ghost_cell.c_idw_sum = 0;

        for (int l=0; l<ghost_cell.c_idw.size(); ++l)
        {
            ghost_cell.c_idw[l] = std::pow((max_distance - ghost_cell.c_idw[l]) / (max_distance * ghost_cell.c_idw[l]), 0.5) + Constants::dsmall;
            ghost_cell.c_idw_sum += ghost_cell.c_idw[l];
        }
    }

    // Set the ghost cells according to the chosen boundary conditions
    void set_ghost_cells(std::vector<Ghost_cell> const &ghost_cells, double* const restrict field, const double boundary_value,
                         const double* const restrict x, const double* const restrict y, const double* const restrict z,
                         const int n_idw, const int ii, const int jj, const int kk,
                         Immersed_boundary::Boundary_type bc, const double visc, bool print=false)
    {
        const int n = (bc == Immersed_boundary::Boundary_type::Dirichlet_type) ? n_idw-1 : n_idw;

        for (std::vector<Ghost_cell>::const_iterator it=ghost_cells.begin(); it<ghost_cells.end(); ++it)
        {
            // Sum the IDW coefficient times the value at the neighbouring grid points
            double vI = 0;
            for (int i=0; i<n; ++i)
                vI += it->c_idw[i] * field[it->neighbours[i].ijk];

            // For Dirichlet BCs, add the boundary value
            if (bc == Immersed_boundary::Boundary_type::Dirichlet_type)
                vI += it->c_idw[n] * boundary_value;

            vI /= it->c_idw_sum;

            // Set the correct BC in the ghost cell
            if (bc == Immersed_boundary::Boundary_type::Dirichlet_type)
                field[it->ijk] = 2*boundary_value - vI;         // Image value reflected across IB
            else if (bc == Immersed_boundary::Boundary_type::Neumann_type)
                field[it->ijk] = vI - boundary_value * it->dI;  // Image value minus gradient times distance
            else if (bc == Immersed_boundary::Boundary_type::Flux_type)
            {
                const double grad = -boundary_value / visc;
                field[it->ijk] = vI - grad * it->dI;            // Image value minus gradient times distance
            }
        }
    }

    // Set the tendencies in the ghost cells to zero
    void zero_ghost_tendency(std::vector<Ghost_cell> const &ghost_cells, double* const restrict ft)
    {

        for (std::vector<Ghost_cell>::const_iterator it=ghost_cells.begin(); it<ghost_cells.end(); ++it)
            ft[it->ijk] = 0;
    }

    // Interpolate (bi-linear) the IB DEM to the requested location
    double interp2_dem(const double x_req, const double y_req,
                       const double* const restrict x, const double* const restrict y,
                       const double* const restrict dem,
                       const double dx, const double dy,
                       const int igc, const int jgc, const int icells,
                       const int imax, const int jmax,
                       const int mpicoordx, const int mpicoordy)
    {
        const int ii = 1;
        const int jj = icells;

        const int i0 = (x_req - 0.5*dx) / dx - mpicoordx*imax + igc;
        const int j0 = (y_req - 0.5*dy) / dy - mpicoordy*jmax + jgc;
        const int ij = i0 + j0*jj;

        const double f1x = (x_req - x[i0]) / dx;
        const double f1y = (y_req - y[j0]) / dy;
        const double f0x = 1 - f1x;
        const double f0y = 1 - f1y;

        const double z = f0y * (f0x * dem[ij   ] + f1x * dem[ij+ii   ]) +
                         f1y * (f0x * dem[ij+jj] + f1x * dem[ij+ii+jj]);

        return z;
    }

    void find_k_dem(
            int* const restrict k_dem,
            const double* const restrict dem,
            const double* const restrict z,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
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

    void calc_fluxes(
            double* const restrict flux,
            const int* const restrict k_dem,
            const double* restrict s,
            const double dx, const double dy, const double* const restrict dz,
            const double dxi, const double dyi, const double* const restrict dzhi,
            const double svisc,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
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

    bool has_ending(const std::string& full_string, const std::string& ending)
    {
        if (full_string.length() >= ending.length())
            return (0 == full_string.compare(full_string.length() - ending.length(), ending.length(), ending));
        else
            return false;
    };
}

Immersed_boundary::Immersed_boundary(Model* modelin, Input* inputin)
{
    model  = modelin;
    fields = model->fields;
    grid   = model->grid;

    int nerror = 0;

    // Get the IB switch, default to sw_ib = false
    inputin->get_item(&sw_ib, "IB", "sw_ib", "", "0");

    // Get the IB type, and required parameters for that type
    if (sw_ib == "0")
        ib_type = None_type;
    else if (sw_ib == "poly")
        ib_type = Poly_type;
    else if (sw_ib == "dem")
        ib_type = Dem_type;
    else if (sw_ib == "sine")
    {
        ib_type = Sine_type;
        nerror += inputin->get_item(&xy_dims,      "IB", "xy_dims",      "", 1 );
        nerror += inputin->get_item(&amplitude,    "IB", "amplitude",    ""    );
        nerror += inputin->get_item(&wavelength_x, "IB", "wavelength_x", ""    );
        nerror += inputin->get_item(&wavelength_y, "IB", "wavelength_y", "", -1);
        nerror += inputin->get_item(&z_offset,     "IB", "z_offset",     "", 0 );
    }
    else if (sw_ib == "gaussian" || sw_ib == "agnesi")
    {
        if (sw_ib == "gaussian")
            ib_type = Gaus_type;
        else if (sw_ib == "agnesi")
            ib_type = Agnesi_type;
        nerror += inputin->get_item(&xy_dims,      "IB", "xy_dims",      "", 1 );
        nerror += inputin->get_item(&amplitude,    "IB", "amplitude",    ""    );
        nerror += inputin->get_item(&x0_hill,      "IB", "x0_hill",      ""    );
        nerror += inputin->get_item(&sigma_x_hill, "IB", "sigma_x_hill", ""    );
        nerror += inputin->get_item(&y0_hill,      "IB", "y0_hill",      "", -1);
        nerror += inputin->get_item(&sigma_y_hill, "IB", "sigma_y_hill", "", -1);
        nerror += inputin->get_item(&z_offset,     "IB", "z_offset",     "", 0 );
    }
    else if (sw_ib == "flat")
    {
        ib_type = Flat_type;
        nerror += inputin->get_item(&z_offset,     "IB", "z_offset",     "", 0 );
    }
    else
    {
        model->master->print_error("sw_ib = \"%s\" not (yet) supported\n", sw_ib.c_str());
        ++nerror;
    }

    if (ib_type != None_type)
    {
        nerror += inputin->get_item(&n_idw,           "IB", "n_idw",           "");          // Number of grid points used in interpolation
        nerror += inputin->get_item(&zero_ghost_tend, "IB", "zero_ghost_tend", "", false);   // Zero tendencies in ghost cells

        // Fixed viscosity at the wall, in case of LES
        if (model->diff->get_switch() == "smag2")
            nerror += inputin->get_item(&visc_wall, "IB", "visc_wall", "");

        grid->set_minimum_ghost_cells(2, 2, 1);  // Set at leat two ghost cells in the horizontal

        // Get the scalar boundary conditions. For now, the momentum BC is hardcoded at no-slip.
        if (fields->sp.size() > 0)
        {
            std::string swbot;
            for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
            {
                sbc[it->first] = new Field3dBc;
                nerror += inputin->get_item(&swbot, "IB", "sbcbot", it->first);
                nerror += inputin->get_item(&sbc[it->first]->bot, "IB", "sbot", it->first);

                // Set the bottom BC
                if (swbot == "dirichlet")
                    sbc[it->first]->bcbot = Dirichlet_type;
                else if (swbot == "neumann")
                    sbc[it->first]->bcbot = Neumann_type;
                else if (swbot == "flux")
                    sbc[it->first]->bcbot = Flux_type;
                else
                {
                    model->master->print_error("%s is illegal value for sbcbot\n", swbot.c_str());
                    nerror++;
                }
            }

            // Check whether all scalar BCs are the same type
            Boundary_type bc = sbc[fields->sp.begin()->first]->bcbot;

            for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
            {
                if (sbc[it->first]->bcbot != bc)
                {
                    model->master->print_error("All IB \"sbcbot\" need to be of the same type!\n");
                    nerror++;
                }
            }
        }
    }

    if (nerror > 0)
        throw 1;
}

Immersed_boundary::~Immersed_boundary()
{
}

void Immersed_boundary::init()
{
    stats = model->stats;

    if (ib_type == None_type)
        return;

    if (ib_type == Dem_type)
    {
        // Resize the 2D DEM field
        dem.resize(grid->ijcells);
        k_dem.resize(grid->ijcells);

        // CvH TEMPORARY
        sbot_ib.resize(grid->ijcells);
        // CvH

        std::vector<std::string>& crosslist_global = model->cross->get_crosslist();

        // Check input list of cross variables (crosslist)
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
                if (fields->sp.find(scalar) != fields->sp.end())
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
}

void Immersed_boundary::create()
{
    if (ib_type == None_type)
        return;

    if (ib_type == Poly_type)
    {
        std::string file_name;

        file_name = "u.ib_input";
        read_ghost_cells(ghost_cells_u, file_name, grid->xh, grid->y, grid->z, Dirichlet_type);

        file_name = "v.ib_input";
        read_ghost_cells(ghost_cells_v, file_name, grid->x, grid->yh, grid->z, Dirichlet_type);

        file_name = "w.ib_input";
        read_ghost_cells(ghost_cells_w, file_name, grid->x, grid->y, grid->zh, Dirichlet_type);

        if (fields->sp.size() > 0)
        {
            Boundary_type bc = sbc[fields->sp.begin()->first]->bcbot;
            file_name = "s.ib_input";
            read_ghost_cells(ghost_cells_s, file_name, grid->x, grid->y, grid->z, bc);
        }
    }
    else
    {
        Boundary_type bc;
        if (fields->sp.size() > 0)
            bc = sbc[fields->sp.begin()->first]->bcbot;

        // Find the IB ghost cells
        if (ib_type == Dem_type)
        {
            // 1. Read input DEM
            char filename[256] = "dem.0000000";
            model->master->print_message("Loading \"%s\" ... ", filename);
            if(grid->load_xy_slice(dem.data(), fields->atmp["tmp1"]->data, filename))
            {
                model->master->print_message("FAILED\n");
                throw 1;
            }
            else
                model->master->print_message("OK\n");
            grid->boundary_cyclic_2d(dem.data());

            // Find ghost cells
            find_ghost_cells<Dem_type, 1>(ghost_cells_u, grid->xh, grid->y,  grid->z,  grid->kstart,   Dirichlet_type);
            find_ghost_cells<Dem_type, 1>(ghost_cells_v, grid->x,  grid->yh, grid->z,  grid->kstart,   Dirichlet_type);
            find_ghost_cells<Dem_type, 1>(ghost_cells_w, grid->x,  grid->y,  grid->zh, grid->kstart+1, Dirichlet_type);

            if (fields->sp.size() > 0)
                find_ghost_cells<Dem_type, 1>(ghost_cells_s, grid->x,  grid->y,  grid->z, grid->kstart,  bc);

            // Create the array with vertical indices that give the first cell above the DEM.
            // Calculate which neighbors are DEM.
            find_k_dem(
                    k_dem.data(), dem.data(), grid->z,
                    grid->istart, grid->iend, grid->jstart, grid->jend, grid->kstart, grid->kend,
                    grid->icells);

            char filename2[256] = "sbot_ib.0000000";
            model->master->print_message("Loading \"%s\" ... ", filename2);
            if (grid->load_xy_slice(sbot_ib.data(), fields->atmp["tmp1"]->data, filename2))
            {
                model->master->print_message("FAILED\n");
                throw 1;
            }
            else
                model->master->print_message("OK\n");
            grid->boundary_cyclic_2d(sbot_ib.data());
        }

        if (ib_type == Flat_type)
        {
            find_ghost_cells<Flat_type, 1>(ghost_cells_u, grid->xh, grid->y,  grid->z,  grid->kstart,   Dirichlet_type);
            find_ghost_cells<Flat_type, 1>(ghost_cells_v, grid->x,  grid->yh, grid->z,  grid->kstart,   Dirichlet_type);
            find_ghost_cells<Flat_type, 1>(ghost_cells_w, grid->x,  grid->y,  grid->zh, grid->kstart+1, Dirichlet_type);

            if (fields->sp.size() > 0)
                find_ghost_cells<Flat_type, 1>(ghost_cells_s, grid->x,  grid->y,  grid->z, grid->kstart,  bc);
        }
        else if (xy_dims == 1)
        {
            if (ib_type == Sine_type)
            {
                find_ghost_cells<Sine_type, 1>(ghost_cells_u, grid->xh, grid->y,  grid->z,  grid->kstart,   Dirichlet_type);
                find_ghost_cells<Sine_type, 1>(ghost_cells_v, grid->x,  grid->yh, grid->z,  grid->kstart,   Dirichlet_type);
                find_ghost_cells<Sine_type, 1>(ghost_cells_w, grid->x,  grid->y,  grid->zh, grid->kstart+1, Dirichlet_type);
                if (fields->sp.size() > 0)
                    find_ghost_cells<Sine_type, 1>(ghost_cells_s, grid->x,  grid->y,  grid->z, grid->kstart,  bc);
            }
            else if (ib_type == Gaus_type)
            {
                find_ghost_cells<Gaus_type, 1>(ghost_cells_u, grid->xh, grid->y,  grid->z,  grid->kstart,   Dirichlet_type);
                find_ghost_cells<Gaus_type, 1>(ghost_cells_v, grid->x,  grid->yh, grid->z,  grid->kstart,   Dirichlet_type);
                find_ghost_cells<Gaus_type, 1>(ghost_cells_w, grid->x,  grid->y,  grid->zh, grid->kstart+1, Dirichlet_type);
                if (fields->sp.size() > 0)
                    find_ghost_cells<Gaus_type, 1>(ghost_cells_s, grid->x,  grid->y,  grid->z, grid->kstart,  bc);
            }
            else if (ib_type == Agnesi_type)
            {
                find_ghost_cells<Agnesi_type, 1>(ghost_cells_u, grid->xh, grid->y,  grid->z,  grid->kstart,   Dirichlet_type);
                find_ghost_cells<Agnesi_type, 1>(ghost_cells_v, grid->x,  grid->yh, grid->z,  grid->kstart,   Dirichlet_type);
                find_ghost_cells<Agnesi_type, 1>(ghost_cells_w, grid->x,  grid->y,  grid->zh, grid->kstart+1, Dirichlet_type);
                if (fields->sp.size() > 0)
                    find_ghost_cells<Agnesi_type, 1>(ghost_cells_s, grid->x,  grid->y,  grid->z, grid->kstart,  bc);
            }
        }
        else if (xy_dims == 2)
        {
            if (ib_type == Sine_type)
            {
                find_ghost_cells<Sine_type, 2>(ghost_cells_u, grid->xh, grid->y,  grid->z,  grid->kstart,   Dirichlet_type);
                find_ghost_cells<Sine_type, 2>(ghost_cells_v, grid->x,  grid->yh, grid->z,  grid->kstart,   Dirichlet_type);
                find_ghost_cells<Sine_type, 2>(ghost_cells_w, grid->x,  grid->y,  grid->zh, grid->kstart+1, Dirichlet_type);
                if (fields->sp.size() > 0)
                    find_ghost_cells<Sine_type, 2>(ghost_cells_s, grid->x,  grid->y,  grid->z, grid->kstart,  bc);
            }
            else if (ib_type == Gaus_type)
            {
                find_ghost_cells<Gaus_type, 2>(ghost_cells_u, grid->xh, grid->y,  grid->z,  grid->kstart,   Dirichlet_type);
                find_ghost_cells<Gaus_type, 2>(ghost_cells_v, grid->x,  grid->yh, grid->z,  grid->kstart,   Dirichlet_type);
                find_ghost_cells<Gaus_type, 2>(ghost_cells_w, grid->x,  grid->y,  grid->zh, grid->kstart+1, Dirichlet_type);
                if (fields->sp.size() > 0)
                    find_ghost_cells<Gaus_type, 2>(ghost_cells_s, grid->x,  grid->y,  grid->z, grid->kstart,  bc);
            }
            else if (ib_type == Agnesi_type)
            {
                find_ghost_cells<Agnesi_type, 2>(ghost_cells_u, grid->xh, grid->y,  grid->z,  grid->kstart,   Dirichlet_type);
                find_ghost_cells<Agnesi_type, 2>(ghost_cells_v, grid->x,  grid->yh, grid->z,  grid->kstart,   Dirichlet_type);
                find_ghost_cells<Agnesi_type, 2>(ghost_cells_w, grid->x,  grid->y,  grid->zh, grid->kstart+1, Dirichlet_type);
                if (fields->sp.size() > 0)
                    find_ghost_cells<Agnesi_type, 2>(ghost_cells_s, grid->x,  grid->y,  grid->z, grid->kstart,  bc);
            }
        }
    }

    // Print some debugging output
    int n_ghost_u = ghost_cells_u.size();
    int n_ghost_v = ghost_cells_v.size();
    int n_ghost_w = ghost_cells_w.size();

    model->master->sum(&n_ghost_u, 1);
    model->master->sum(&n_ghost_v, 1);
    model->master->sum(&n_ghost_w, 1);

    model->master->print_message("Found %i IB[u] ghost cells \n", n_ghost_u);
    model->master->print_message("Found %i IB[v] ghost cells \n", n_ghost_v);
    model->master->print_message("Found %i IB[w] ghost cells \n", n_ghost_w);

    if (fields->sp.size() > 0)
    {
        int n_ghost_s = ghost_cells_s.size();
        model->master->sum(&n_ghost_s, 1);
        model->master->print_message("Found %i IB[s] ghost cells \n", n_ghost_s);
    }
}

// Return the height of the IB as a function of x,y position
template <Immersed_boundary::IB_type sw, int dims>
double Immersed_boundary::boundary_function(const double x, const double y)
{
    if (sw == Flat_type)
        return z_offset;
    else if (sw == Dem_type)
        return interp2_dem(x, y, grid->x, grid->y, dem.data(), grid->dx, grid->dy,
                           grid->igc, grid->jgc, grid->icells, grid->imax, grid->jmax,
                           model->master->mpicoordx, model->master->mpicoordy);
    else if (dims == 1)
    {
        if (sw == Sine_type)
            return z_offset + amplitude + amplitude * std::sin(2*Constants::pi*x/wavelength_x);
        else if (sw == Gaus_type)
            return z_offset + amplitude * std::exp(-pow((x-x0_hill)/(2*sigma_x_hill), 2));
        else if (sw == Agnesi_type)
            return z_offset + amplitude / (1. + pow((x-x0_hill)/sigma_x_hill, 2));
    }
    else if (dims == 2)
    {
        if (sw == Sine_type)
            return z_offset + amplitude + amplitude * std::sin(2*Constants::pi*x/wavelength_x)
                                                    * std::sin(2*Constants::pi*y/wavelength_y);
        else if (sw == Gaus_type)
            return z_offset + amplitude * std::exp(-pow((x-x0_hill)/(2*sigma_x_hill), 2))
                                        * std::exp(-pow((y-y0_hill)/(2*sigma_y_hill), 2));
        else if (sw == Agnesi_type)
            return z_offset + amplitude / (1. + pow((x-x0_hill)/sigma_x_hill, 2)
                                              + pow((y-y0_hill)/sigma_y_hill, 2));
    }

    return 0;   // Otherwise the compiler complains
}

template <Immersed_boundary::IB_type sw, int dims>
bool Immersed_boundary::is_ghost_cell(const double* const restrict x, const double* const restrict y, const double* const restrict z,
                                      const int i, const int j, const int k)
{
    if (z[k] <= boundary_function<sw, dims>(x[i], y[j]))  // Inside IB
    {
        // Check if one of the neighbouring grid cells is outside the IB
        for (int dk=-1; dk<2; ++dk)
            if (z[k+dk] > boundary_function<sw, dims>(x[i], y[j]))
                return true;
        for (int dj=-1; dj<2; ++dj)
            if (z[k] > boundary_function<sw, dims>(x[i], y[j+dj]))
                return true;
        for (int di=-1; di<2; ++di)
            if (z[k] > boundary_function<sw, dims>(x[i+di], y[j]))
                return true;
    }
    return false;
}

template <Immersed_boundary::IB_type sw, int dims>
void Immersed_boundary::find_nearest_location_wall(double& x_min, double& y_min, double& z_min, double& d_min,
                                                   const double x, const double y, const double z,
                                                   const int i, const int j, const int k)
{
    const double dx = grid->dx;
    const double dy = grid->dy;

    d_min = Constants::dbig;
    const int n  = 40;

    for (int ii = -n/2; ii < n/2+1; ++ii)
        for (int jj = -n/2; jj < n/2+1; ++jj)
        {
            const double xc = x + 2*ii/(double)n*dx;
            const double yc = y + 2*jj/(double)n*dy;
            const double zc = boundary_function<sw, dims>(xc, yc);
            const double d  = abs_distance(x, xc, y, yc, z, zc);

            if (d < d_min)
            {
                d_min = d;
                x_min = xc;
                y_min = yc;
                z_min = zc;
            }
        }
}

template <Immersed_boundary::IB_type sw, int dims>
void Immersed_boundary::find_interpolation_points(Ghost_cell& ghost_cell,
                                                  const double* const restrict x, const double* const restrict y, const double* const restrict z,
                                                  const int i, const int j, const int k, Boundary_type bc)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    double x_min, y_min, z_min, d_min;  // x, y, z location and distance to wall

    // Minimal distance that a grid point used in the interpolation has to be
    // away from the IB, to prevent getting extreme extrapolated values in the ghost cells
    const double d_lim = 0.1 * std::min(std::min(grid->dx, grid->dy), grid->dz[k]);

    // Limit vertical stencil near surface
    double dk0 = -2;
    if (k+dk0 < grid->kstart)
        dk0 = grid->kstart - k;

    // Find the neighbouring grid points outside the IB
    for (int dk=dk0; dk<3; ++dk)
        for (int dj=-2; dj<3; ++dj)
            for (int di=-2; di<3; ++di)
            {
                // Check if grid point is outside IB
                if (z[k+dk] > boundary_function<sw, dims>(x[i+di], y[j+dj]))
                {
                    // Calculate distance (d_min) of current grid point to the IB
                    find_nearest_location_wall<sw, dims>(x_min, y_min, z_min, d_min, x[i+di], y[j+dj], z[k+dk], i+di, j+dj, k+dk);

                    // As described above; exclude grid points which are close to the IB
                    if (d_min > d_lim)
                    {
                        const int ijk = (i+di) + (j+dj)*jj + (k+dk)*kk;
                        Neighbour tmp_neighbour = {i+di, j+dj, k+dk, ijk, abs_distance(x[i], x[i+di], y[j], y[j+dj], z[k], z[k+dk])};
                        ghost_cell.neighbours.push_back(tmp_neighbour);
                    }
                }
            }

    // Sort them on distance:
    std::sort(ghost_cell.neighbours.begin(), ghost_cell.neighbours.end(), compare_value);

    // Abort if there are insufficient neighbouring fluid points
    if (ghost_cell.neighbours.size() < n_idw)
    {
        model->master->print_error("Only found %i of n_idw=%i neighbour points\n", ghost_cell.neighbours.size(), n_idw);
        throw 1;
    }

    // For Dirichlet BCs, one interpolation point is the boundary value at the IB,
    // so we only require n_idw-1 interpolation points outside the IB
    const int n = (bc == Dirichlet_type) ? n_idw-1 : n_idw;

    // Only keep the N nearest neighbour fluid points
    ghost_cell.neighbours.erase(ghost_cell.neighbours.begin()+n, ghost_cell.neighbours.end());

    // Pre-calculate the inverse distance weighting coefficients
    precalculate_idw(ghost_cell, x, y, z, bc);
}


void Immersed_boundary::read_ghost_cells(std::vector<Ghost_cell> &ghost_cells, std::string file_name,
                                         const double* const restrict x, const double* const restrict y, const double* const restrict z,
                                         Boundary_type bc)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    bool is_optional = false;
    Data_map input;

    // Read the input data into a std::map< header_name, std::vector<data> >
    const int nerror = model->input->read_data_file(&input, file_name, is_optional);
    if (nerror > 0)
        throw 1;

    const int mpioffsx = model->master->mpicoordx * grid->imax;
    const int mpioffsy = model->master->mpicoordy * grid->jmax;

    const int i0 = mpioffsx + grid->igc;
    const int i1 = mpioffsx + grid->igc + grid->imax;
    const int j0 = mpioffsy + grid->jgc;
    const int j1 = mpioffsy + grid->jgc + grid->jmax;

    for (int n=0; n<input["i"].size(); ++n)
    {
        int i = input["i"][n] + grid->igc;
        int j = input["j"][n] + grid->jgc;
        int k = input["k"][n] + grid->kgc;

        // Check if grid point on this MPI task
        if (i >= i0 && i < i1 && j >= j0 && j < j1)
        {
            i -= mpioffsx;
            j -= mpioffsy;

            const int ijk = i + j*jj + k*kk;

            Ghost_cell tmp_ghost = {i, j, k, ijk};

            // Location on boundary
            tmp_ghost.xB = input["xb"][n];
            tmp_ghost.yB = input["yb"][n];
            tmp_ghost.zB = input["zb"][n];

            // Location image point
            tmp_ghost.xI = 2*tmp_ghost.xB - x[i];
            tmp_ghost.yI = 2*tmp_ghost.yB - y[j];
            tmp_ghost.zI = 2*tmp_ghost.zB - z[k];

            // Neighbours
            const int n_neighbours = (bc == Immersed_boundary::Boundary_type::Dirichlet_type) ? n_idw-1 : n_idw;

            for (int nn=0; nn<n_neighbours; ++nn)
            {
                // Note BvS: the cast to long long is only necessary as intel currently doesn't have the full c++11 set of to_string() implemented
                const int in = input["i"+std::to_string(static_cast<long long>(nn))][n] + grid->istart - mpioffsx;
                const int jn = input["j"+std::to_string(static_cast<long long>(nn))][n] + grid->jstart - mpioffsy;
                const int kn = input["k"+std::to_string(static_cast<long long>(nn))][n] + grid->kstart;

                const int ijk = in + jn*jj + kn*kk;

                Neighbour tmp_neighbour = {in, jn, kn, ijk, abs_distance(x[i], x[in], y[j], y[jn], z[k], z[kn])};
                tmp_ghost.neighbours.push_back(tmp_neighbour);
            }

            precalculate_idw(tmp_ghost, x, y, z, Dirichlet_type);

            ghost_cells.push_back(tmp_ghost);
        }
    }
}

template <Immersed_boundary::IB_type sw, int dims>
void Immersed_boundary::find_ghost_cells(std::vector<Ghost_cell> &ghost_cells,
                                         const double* const restrict x, const double* const restrict y, const double* const restrict z,
                                         const int kstart_search, Boundary_type bc)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    double d_wall;  // Distance to IB (m)

    for (int k=kstart_search; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                // 1. Check if this is a ghost cell, i.e. inside the IB, with a neighbour outside the IB
                if (is_ghost_cell<sw, dims>(x, y, z, i, j, k))
                {
                    const int ijk = i + j*jj + k*kk;

                    Ghost_cell tmp_ghost = {i, j, k, ijk};

                    // 2. Find the closest location on the IB
                    find_nearest_location_wall<sw, dims>(tmp_ghost.xB, tmp_ghost.yB, tmp_ghost.zB, d_wall, x[i], y[j], z[k], i, j, k);

                    // 2.1 Location image point
                    tmp_ghost.xI = 2*tmp_ghost.xB - x[i];
                    tmp_ghost.yI = 2*tmp_ghost.yB - y[j];
                    tmp_ghost.zI = 2*tmp_ghost.zB - z[k];

                    // 2.2 Save distance ghost cell to image point
                    tmp_ghost.dI = abs_distance(tmp_ghost.xI, x[i], tmp_ghost.yI, y[j], tmp_ghost.zI, z[k]);

                    // 3. Find the closest `n_idw` grid points outside the IB
                    find_interpolation_points<sw, dims>(tmp_ghost, x, y, z, i, j, k, bc);

                    // 4. Add to collection (list) of ghost cells
                    ghost_cells.push_back(tmp_ghost);
                }
            }
}

void Immersed_boundary::exec_stats(Mask *m)
{
    if (ib_type == None_type)
        return;
    else
        return;  // ...
}

void Immersed_boundary::exec_cross()
{
    if (ib_type == None_type)
        return;

    else if (ib_type == Dem_type)
    {
        for (const auto& s : crosslist)
        {
            const std::string fluxbot_ib_string = "fluxbot_ib";
            if (has_ending(s, fluxbot_ib_string))
            {
                // Strip the scalar from the fluxbot_ib
                std::string scalar = s;
                scalar.erase(s.length() - fluxbot_ib_string.length());

                calc_fluxes(
                        fields->atmp["tmp1"]->datafluxbot,
                        k_dem.data(),
                        fields->sp[scalar]->data,
                        grid->dx, grid->dy, grid->dz,
                        grid->dxi, grid->dyi, grid->dzhi,
                        fields->sp[scalar]->visc,
                        grid->istart, grid->iend, grid->jstart, grid->jend, grid->kstart, grid->kend,
                        grid->icells, grid->ijcells);

                // Calculate the scalar fluxes.
                int nerror = 0;
                nerror += model->cross->cross_plane(fields->atmp["tmp1"]->datafluxbot, fields->atmp["tmp1"]->data, scalar + "fluxbot_ib");

                if (nerror > 0)
                    throw 1;
            }
        }
    }
    else
        return;  // ...
}


template <Immersed_boundary::IB_type sw, int dims>
void Immersed_boundary::calc_mask(double* const restrict mask, double* const restrict maskh, double* const restrict maskbot,
                                  int* const restrict nmask, int* const restrict nmaskh, int* const restrict nmaskbot,
                                  const double* const restrict x, const double* const restrict y,
                                  const double* const restrict z, const double* const restrict zh)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    // Temporarily store IB heigh in maskbot
    double* const restrict zb = maskbot;

    for (int j=grid->jstart; j<grid->jend; j++)
        #pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ij  = i + j*jj;
            zb[ij] = boundary_function<sw,dims>(x[i], y[j]);
        }

    // Set the mask for outside (1) or inside (0) IB
    for (int k=grid->kstart; k<grid->kend; ++k)
    {
        nmask[k]  = 0;
        nmaskh[k] = 0;
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                const int ij  = i + j*jj;

                const int is_not_ib    = z [k] > zb[ij];
                const int is_not_ib_h  = zh[k] > zb[ij];

                mask[ijk]  = static_cast<double>(is_not_ib);
                maskh[ijk] = static_cast<double>(is_not_ib_h);

                nmask[k]  += is_not_ib;
                nmaskh[k] += is_not_ib_h;
            }
    }

    // Mask for surface projected quantities
    for (int j=grid->jstart; j<grid->jend; j++)
        #pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            maskbot[ij] = maskh[ijk];
        }

    grid->boundary_cyclic(mask);
    grid->boundary_cyclic(maskh);
    grid->boundary_cyclic_2d(maskbot);

    model->master->sum(nmask , grid->kcells);
    model->master->sum(nmaskh, grid->kcells);
    *nmaskbot = nmaskh[grid->kstart];
}


void Immersed_boundary::get_mask(Field3d *mfield, Field3d *mfieldh)
{
    // Mask is currently only implemented for boundaries set up from f(x,y)
    if ((ib_type != Sine_type) && (ib_type != Gaus_type) && (ib_type != Agnesi_type) && (ib_type != Flat_type) && (ib_type != Dem_type))
    {
        model->master->print_error("Get_mask() not yet implemented for chosen IB type\n");
        throw 1;
    }

    if (ib_type == Dem_type)
        calc_mask<Dem_type, 1>(mfield->data, mfieldh->data, mfieldh->databot, stats->nmask, stats->nmaskh, &stats->nmaskbot, grid->x, grid->y, grid->z, grid->zh);
    else if (ib_type == Flat_type)
        calc_mask<Flat_type, 1>(mfield->data, mfieldh->data, mfieldh->databot, stats->nmask, stats->nmaskh, &stats->nmaskbot, grid->x, grid->y, grid->z, grid->zh);
    else if (xy_dims == 1)
    {
        if (ib_type == Sine_type)
            calc_mask<Sine_type, 1>(mfield->data, mfieldh->data, mfieldh->databot, stats->nmask, stats->nmaskh, &stats->nmaskbot, grid->x, grid->y, grid->z, grid->zh);
        else if (ib_type == Gaus_type)
            calc_mask<Gaus_type, 1>(mfield->data, mfieldh->data, mfieldh->databot, stats->nmask, stats->nmaskh, &stats->nmaskbot, grid->x, grid->y, grid->z, grid->zh);
        else if (ib_type == Agnesi_type)
            calc_mask<Agnesi_type, 1>(mfield->data, mfieldh->data, mfieldh->databot, stats->nmask, stats->nmaskh, &stats->nmaskbot, grid->x, grid->y, grid->z, grid->zh);
    }
    else if (xy_dims == 2)
    {
        if (ib_type == Sine_type)
            calc_mask<Sine_type, 2>(mfield->data, mfieldh->data, mfieldh->databot, stats->nmask, stats->nmaskh, &stats->nmaskbot, grid->x, grid->y, grid->z, grid->zh);
        else if (ib_type == Gaus_type)
            calc_mask<Gaus_type, 2>(mfield->data, mfieldh->data, mfieldh->databot, stats->nmask, stats->nmaskh, &stats->nmaskbot, grid->x, grid->y, grid->z, grid->zh);
        else if (ib_type == Agnesi_type)
            calc_mask<Agnesi_type, 2>(mfield->data, mfieldh->data, mfieldh->databot, stats->nmask, stats->nmaskh, &stats->nmaskbot, grid->x, grid->y, grid->z, grid->zh);
    }
}

template <Immersed_boundary::IB_type sw, int dims>
void Immersed_boundary::zero_ib_tendency(double* const restrict tend, double* const restrict tmp2d,
                                         const double* const restrict x, const double* const restrict y, const double* const restrict z)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    // Calculate IB height in 2D tmp field
    for (int j=grid->jstart; j<grid->jend; j++)
        #pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ij  = i + j*jj;
            tmp2d[ij] = boundary_function<sw,dims>(x[i], y[j]);
        }


    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + k*kk;

                if (z[k] < tmp2d[ij])
                    tend[ijk] = 0.;
            }
}

void Immersed_boundary::exec_momentum()
{
    if (ib_type == None_type)
        return;

    const int ii = 1;
    const double no_slip = 0.;

    // Set the immersed boundary ghost cells to enforce a no-slip BC for the velocity components
    set_ghost_cells(ghost_cells_u, fields->u->data, no_slip, grid->xh, grid->y,  grid->z,  n_idw, ii, grid->icells, grid->ijcells, Dirichlet_type, fields->visc);
    set_ghost_cells(ghost_cells_v, fields->v->data, no_slip, grid->x,  grid->yh, grid->z,  n_idw, ii, grid->icells, grid->ijcells, Dirichlet_type, fields->visc);
    set_ghost_cells(ghost_cells_w, fields->w->data, no_slip, grid->x,  grid->y,  grid->zh, n_idw, ii, grid->icells, grid->ijcells, Dirichlet_type, fields->visc);

    grid->boundary_cyclic(fields->u->data);
    grid->boundary_cyclic(fields->v->data);
    grid->boundary_cyclic(fields->w->data);

    // Zero ghost cell tendencies
    if (zero_ghost_tend)
    {
        zero_ghost_tendency(ghost_cells_u, fields->at["u"]->data);
        zero_ghost_tendency(ghost_cells_v, fields->at["v"]->data);
        zero_ghost_tendency(ghost_cells_w, fields->at["w"]->data);

        grid->boundary_cyclic(fields->at["u"]->data);
        grid->boundary_cyclic(fields->at["v"]->data);
        grid->boundary_cyclic(fields->at["w"]->data);

        for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
        {
            zero_ghost_tendency(ghost_cells_s, fields->st[it->first]->data);
            grid->boundary_cyclic(fields->at[it->first]->data);
        }
    }
}

void Immersed_boundary::exec_scalars()
{
    if (ib_type == None_type)
        return;

    const int ii = 1;

    // For LES, prescribe a fixed viscosity at the IB wall
    if (model->diff->get_switch() == "smag2")
        set_ghost_cells(ghost_cells_s, fields->sd["evisc"]->data, visc_wall, grid->x, grid->y,  grid->z,  n_idw, ii, grid->icells, grid->ijcells, Dirichlet_type, fields->visc);

    // Set the ghost cells for scalars, depending on their BC
    for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
        if (model->diff->get_switch() == "smag2")
            set_ghost_cells(ghost_cells_s, fields->sp[it->first]->data, sbc[it->first]->bot,
                            grid->x, grid->y, grid->z, n_idw, ii, grid->icells, grid->ijcells, sbc[it->first]->bcbot, visc_wall, false);
        else
            set_ghost_cells(ghost_cells_s, fields->sp[it->first]->data, sbc[it->first]->bot,
                            grid->x, grid->y, grid->z, n_idw, ii, grid->icells, grid->ijcells, sbc[it->first]->bcbot, fields->sp[it->first]->visc, false);

        grid->boundary_cyclic(fields->ap[it->first]->data);
    }

    // Zero ghost cell tendencies
    if (zero_ghost_tend)
    {
        for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
        {
            zero_ghost_tendency(ghost_cells_s, fields->st[it->first]->data);
            grid->boundary_cyclic(fields->at[it->first]->data);
        }
    }
}
