/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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
#include "master.h"
#include "model.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "stats.h"
#include "constants.h"
#include "finite_difference.h"
#include "immersed_boundary.h"

#include <fstream>  // TMP BvS
#include <sstream>

namespace
{
    void write_debug(std::vector<Ghost_cell> ghost_cells, std::string name)
    {
        std::ofstream debugfile;
        debugfile.open(name);

        for (std::vector<Ghost_cell>::iterator it=ghost_cells.begin(); it<ghost_cells.end(); ++it)
        {
            debugfile << it->i                << ", " <<
                         it->j                << ", " <<
                         it->k                << ", " <<
                         it->xB               << ", " <<
                         it->yB               << ", " <<
                         it->zB               << ", " <<
                         it->neighbours[0].i  << ", " <<
                         it->neighbours[0].j  << ", " <<
                         it->neighbours[0].k  << ", " <<
                         it->neighbours[1].i  << ", " <<
                         it->neighbours[1].j  << ", " <<
                         it->neighbours[1].k  << ", " <<
                         it->neighbours[2].i  << ", " <<
                         it->neighbours[2].j  << ", " <<
                         it->neighbours[2].k  << std::endl;
        }
        debugfile.close();
    }

    inline double abs_distance(const double x1, const double x2, const double y1, const double y2, const double z1, const double z2)
    {
        return std::pow(std::pow(x1-x2, 2) + std::pow(y1-y2, 2) + std::pow(z1-z2, 2), 0.5);
    }

    inline double abs_distance(const double x1, const double x2, const double z1, const double z2)
    {
        return std::pow(std::pow(x1-x2, 2) + std::pow(z1-z2, 2), 0.5);
    }

    // Help function for sorting std::vector with Neighbour points
    bool compare_value(const Neighbour& a, const Neighbour& b)
    {
        return a.distance < b.distance;
    }

    void print_mat(std::vector<std::vector<double> > m)
    {
        for (int i=0; i<m.size(); ++i)
        {
            for (int j=0; j<m.size(); ++j)
            {
                std::cout << std::fixed << std::setw(6) << std::setprecision(3) << m[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }

    void inverse_mat_4x4(std::vector<std::vector<double> >& result, std::vector<std::vector<double> > input)
    {
        result[0][0] =  input[1][1] * input[2][2] * input[3][3] - 
                        input[1][1] * input[3][2] * input[2][3] - 
                        input[1][2] * input[2][1] * input[3][3] + 
                        input[1][2] * input[3][1] * input[2][3] +
                        input[1][3] * input[2][1] * input[3][2] - 
                        input[1][3] * input[3][1] * input[2][2];

        result[0][1] = -input[0][1] * input[2][2] * input[3][3] + 
                        input[0][1] * input[3][2] * input[2][3] + 
                        input[0][2] * input[2][1] * input[3][3] - 
                        input[0][2] * input[3][1] * input[2][3] - 
                        input[0][3] * input[2][1] * input[3][2] + 
                        input[0][3] * input[3][1] * input[2][2];

        result[0][2] =  input[0][1] * input[1][2] * input[3][3] - 
                        input[0][1] * input[3][2] * input[1][3] - 
                        input[0][2] * input[1][1] * input[3][3] + 
                        input[0][2] * input[3][1] * input[1][3] + 
                        input[0][3] * input[1][1] * input[3][2] - 
                        input[0][3] * input[3][1] * input[1][2];

        result[0][3] = -input[0][1] * input[1][2] * input[2][3] + 
                        input[0][1] * input[2][2] * input[1][3] +
                        input[0][2] * input[1][1] * input[2][3] - 
                        input[0][2] * input[2][1] * input[1][3] - 
                        input[0][3] * input[1][1] * input[2][2] + 
                        input[0][3] * input[2][1] * input[1][2];

        result[1][0] = -input[1][0] * input[2][2] * input[3][3] + 
                        input[1][0] * input[3][2] * input[2][3] + 
                        input[1][2] * input[2][0] * input[3][3] - 
                        input[1][2] * input[3][0] * input[2][3] - 
                        input[1][3] * input[2][0] * input[3][2] + 
                        input[1][3] * input[3][0] * input[2][2];

        result[1][1] =  input[0][0] * input[2][2] * input[3][3] - 
                        input[0][0] * input[3][2] * input[2][3] - 
                        input[0][2] * input[2][0] * input[3][3] + 
                        input[0][2] * input[3][0] * input[2][3] + 
                        input[0][3] * input[2][0] * input[3][2] - 
                        input[0][3] * input[3][0] * input[2][2];

        result[1][2] = -input[0][0] * input[1][2] * input[3][3] + 
                        input[0][0] * input[3][2] * input[1][3] + 
                        input[0][2] * input[1][0] * input[3][3] - 
                        input[0][2] * input[3][0] * input[1][3] - 
                        input[0][3] * input[1][0] * input[3][2] + 
                        input[0][3] * input[3][0] * input[1][2];

        result[1][3] =  input[0][0] * input[1][2] * input[2][3] - 
                        input[0][0] * input[2][2] * input[1][3] - 
                        input[0][2] * input[1][0] * input[2][3] + 
                        input[0][2] * input[2][0] * input[1][3] + 
                        input[0][3] * input[1][0] * input[2][2] - 
                        input[0][3] * input[2][0] * input[1][2];

        result[2][0] =  input[1][0] * input[2][1] * input[3][3] - 
                        input[1][0] * input[3][1] * input[2][3] - 
                        input[1][1] * input[2][0] * input[3][3] + 
                        input[1][1] * input[3][0] * input[2][3] + 
                        input[1][3] * input[2][0] * input[3][1] - 
                        input[1][3] * input[3][0] * input[2][1];

        result[2][1] = -input[0][0] * input[2][1] * input[3][3] + 
                        input[0][0] * input[3][1] * input[2][3] + 
                        input[0][1] * input[2][0] * input[3][3] - 
                        input[0][1] * input[3][0] * input[2][3] - 
                        input[0][3] * input[2][0] * input[3][1] + 
                        input[0][3] * input[3][0] * input[2][1];

        result[2][2] =  input[0][0] * input[1][1] * input[3][3] - 
                        input[0][0] * input[3][1] * input[1][3] - 
                        input[0][1] * input[1][0] * input[3][3] + 
                        input[0][1] * input[3][0] * input[1][3] + 
                        input[0][3] * input[1][0] * input[3][1] - 
                        input[0][3] * input[3][0] * input[1][1];

        result[2][3] = -input[0][0] * input[1][1] * input[2][3] + 
                        input[0][0] * input[2][1] * input[1][3] + 
                        input[0][1] * input[1][0] * input[2][3] - 
                        input[0][1] * input[2][0] * input[1][3] - 
                        input[0][3] * input[1][0] * input[2][1] + 
                        input[0][3] * input[2][0] * input[1][1];

        result[3][0] = -input[1][0] * input[2][1] * input[3][2] + 
                        input[1][0] * input[3][1] * input[2][2] + 
                        input[1][1] * input[2][0] * input[3][2] - 
                        input[1][1] * input[3][0] * input[2][2] - 
                        input[1][2] * input[2][0] * input[3][1] + 
                        input[1][2] * input[3][0] * input[2][1];

        result[3][1] =  input[0][0] * input[2][1] * input[3][2] - 
                        input[0][0] * input[3][1] * input[2][2] - 
                        input[0][1] * input[2][0] * input[3][2] + 
                        input[0][1] * input[3][0] * input[2][2] + 
                        input[0][2] * input[2][0] * input[3][1] - 
                        input[0][2] * input[3][0] * input[2][1];

        result[3][2] = -input[0][0] * input[1][1] * input[3][2] + 
                        input[0][0] * input[3][1] * input[1][2] + 
                        input[0][1] * input[1][0] * input[3][2] - 
                        input[0][1] * input[3][0] * input[1][2] - 
                        input[0][2] * input[1][0] * input[3][1] + 
                        input[0][2] * input[3][0] * input[1][1];

        result[3][3] =  input[0][0] * input[1][1] * input[2][2] - 
                        input[0][0] * input[2][1] * input[1][2] - 
                        input[0][1] * input[1][0] * input[2][2] + 
                        input[0][1] * input[2][0] * input[1][2] + 
                        input[0][2] * input[1][0] * input[2][1] - 
                        input[0][2] * input[2][0] * input[1][1];

        double det = input[0][0] * result[0][0] + input[1][0] * result[0][1] + input[2][0] * result[0][2] + input[3][0] * result[0][3];

        if (det == 0)
            std::cout << "DET == 0!" << std::endl;

        det = 1.0 / det;

        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                result[i][j] *= det;
    }

    void set_ghost_cells(std::vector<Ghost_cell> ghost_cells, double* const restrict field, const double boundary_value,
                         const double* const restrict x, const double* const restrict y, const double* const restrict z,
                         const int nn, const int ii, const int jj, const int kk)
    {
        for (std::vector<Ghost_cell>::iterator it=ghost_cells.begin(); it<ghost_cells.end(); ++it)
        {
            // Calculate interpolant at image point
            double vI = 0;
            for (int i=0; i<nn; ++i)
            {
                const int ijk = it->neighbours[i].i + it->neighbours[i].j*jj + it->neighbours[i].k*kk;
                vI += it->c_idw[i] * field[ijk];
            }
            vI += it->c_idw[nn] * boundary_value; 
            vI /= it->c_idw_sum; 

            // Reflect across boundary
            const int ijk = it->i + it->j*jj + it->k*kk;
            field[ijk] = 2*boundary_value - vI;
        }
    }
}

Immersed_boundary::Immersed_boundary(Model* modelin, Input* inputin)
{
    model  = modelin;
    fields = model->fields;
    grid   = model->grid;

    int nerror = 0;

    inputin->get_item(&sw_ib,   "immersed_boundary", "sw_ib",   "", "0");

    // Get the IB type, and required parameters for that type
    if (sw_ib == "0")
    {
        ib_type = None_type;
    }
    else if (sw_ib == "sine")
    {
        ib_type = Sine_type;
        nerror += inputin->get_item(&xy_dims,      "immersed_boundary", "xy_dims",      "", 1 );
        nerror += inputin->get_item(&amplitude,    "immersed_boundary", "amplitude",    ""    );
        nerror += inputin->get_item(&wavelength_x, "immersed_boundary", "wavelength_x", ""    );
        nerror += inputin->get_item(&wavelength_y, "immersed_boundary", "wavelength_y", "", -1);
        nerror += inputin->get_item(&z_offset,     "immersed_boundary", "z_offset",     "", 0 );
    }
    else if (sw_ib == "gaussian")
    {
        ib_type = Gaus_type; 
        nerror += inputin->get_item(&xy_dims,      "immersed_boundary", "xy_dims",      "", 1 );
        nerror += inputin->get_item(&amplitude,    "immersed_boundary", "amplitude",    ""    );
        nerror += inputin->get_item(&x0_hill,      "immersed_boundary", "x0_hill",      ""    );
        nerror += inputin->get_item(&sigma_x_hill, "immersed_boundary", "sigma_x_hill", ""    );
        nerror += inputin->get_item(&y0_hill,      "immersed_boundary", "y0_hill",      "", -1);
        nerror += inputin->get_item(&sigma_y_hill, "immersed_boundary", "sigma_y_hill", "", -1);
        nerror += inputin->get_item(&z_offset,     "immersed_boundary", "z_offset",     "", 0 );
    }
    else
    {
        model->master->print_error("sw_ib = \"%s\" not (yet) supported\n", sw_ib.c_str());
        throw 1;
    }

    // Interpolation method
    if (ib_type != None_type)
    {
        grid->set_minimum_ghost_cells(2, 2, 1);
        inputin->get_item(&sw_interpol, "immersed_boundary", "sw_interpolation", "");

        if (sw_interpol == "linear")
            interpol_type = Linear_type;
        else if (sw_interpol == "inverse_distance")
            interpol_type = Distance_type;
        else
        {
           model->master->print_error("sw_interpolate = \"%s\" not (yet) supported\n", sw_interpol.c_str());
           throw 1;
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
}

void Immersed_boundary::create()
{
    if (ib_type == None_type)
        return;

    if (xy_dims == 1)
    {
        if (ib_type == Sine_type)
        {
            find_ghost_cells<Sine_type, 1>(&ghost_cells_u, grid->xh, grid->y,  grid->z );
            find_ghost_cells<Sine_type, 1>(&ghost_cells_v, grid->x,  grid->yh, grid->z );
            find_ghost_cells<Sine_type, 1>(&ghost_cells_w, grid->x,  grid->y,  grid->zh);
            find_ghost_cells<Sine_type, 1>(&ghost_cells_s, grid->x,  grid->y,  grid->z );
        }
        else if (ib_type == Gaus_type)
        {
            find_ghost_cells<Gaus_type, 1>(&ghost_cells_u, grid->xh, grid->y,  grid->z );
            find_ghost_cells<Gaus_type, 1>(&ghost_cells_v, grid->x,  grid->yh, grid->z );
            find_ghost_cells<Gaus_type, 1>(&ghost_cells_w, grid->x,  grid->y,  grid->zh);
            find_ghost_cells<Gaus_type, 1>(&ghost_cells_s, grid->x,  grid->y,  grid->z );
        }
    }
    else if (xy_dims == 2)
    {
        if (ib_type == Sine_type)
        {
            find_ghost_cells<Sine_type, 2>(&ghost_cells_u, grid->xh, grid->y,  grid->z );
            find_ghost_cells<Sine_type, 2>(&ghost_cells_v, grid->x,  grid->yh, grid->z );
            find_ghost_cells<Sine_type, 2>(&ghost_cells_w, grid->x,  grid->y,  grid->zh);
            find_ghost_cells<Sine_type, 2>(&ghost_cells_s, grid->x,  grid->y,  grid->z );
        }
        else if (ib_type == Gaus_type)
        {
            find_ghost_cells<Gaus_type, 2>(&ghost_cells_u, grid->xh, grid->y,  grid->z );
            find_ghost_cells<Gaus_type, 2>(&ghost_cells_v, grid->x,  grid->yh, grid->z );
            find_ghost_cells<Gaus_type, 2>(&ghost_cells_w, grid->x,  grid->y,  grid->zh);
            find_ghost_cells<Gaus_type, 2>(&ghost_cells_s, grid->x,  grid->y,  grid->z );
        }
    }

    model->master->print_message("Found %i IB[u] ghost cells \n", ghost_cells_u.size());
    model->master->print_message("Found %i IB[v] ghost cells \n", ghost_cells_v.size());
    model->master->print_message("Found %i IB[w] ghost cells \n", ghost_cells_w.size());
    model->master->print_message("Found %i IB[s] ghost cells \n", ghost_cells_s.size());

    // Debug.....
    std::string name;
    name = "debug_u.txt"; write_debug(ghost_cells_u, name);
    name = "debug_v.txt"; write_debug(ghost_cells_v, name);
    name = "debug_w.txt"; write_debug(ghost_cells_w, name);
    name = "debug_s.txt"; write_debug(ghost_cells_s, name);
}

template<Immersed_boundary::IB_type sw, int dims>
double Immersed_boundary::boundary_function(const double x, const double y)
{
    if (dims == 1)
    {
        if (sw == Sine_type)
        {
            const double pi = std::acos((double)-1.); 
            return z_offset + amplitude + amplitude * std::sin(2*pi*x/wavelength_x);
        }
        else if (sw == Gaus_type)
        {
            return z_offset + amplitude * std::exp(-pow((x-x0_hill)/(2*sigma_x_hill), 2));
        }
    }
    else if (dims == 2)
    {
        if (sw == Sine_type)
        {
            const double pi = std::acos((double)-1.); 
            return z_offset + amplitude + amplitude * std::sin(2*pi*x/wavelength_x) * std::sin(2*pi*y/wavelength_y);
        }
        else if (sw == Gaus_type)
        {
            return z_offset + amplitude * std::exp(-pow((x-x0_hill)/(2*sigma_x_hill), 2))
                                        * std::exp(-pow((y-y0_hill)/(2*sigma_y_hill), 2));
        }
    }
}

template<Immersed_boundary::IB_type sw, int dims>
bool Immersed_boundary::is_ghost_cell(const double* const restrict x, const double* const restrict y, const double* const restrict z,
                                      const int i, const int j, const int k)
{
    if (z[k] < boundary_function<sw, dims>(x[i], y[j]))  // Inside IB
    {
        // Check if one of the neighbouring grid cells is outside the IB
        for (int dk=-1; dk<2; ++dk) // BvS Fix this
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

template<Immersed_boundary::IB_type sw, int dims>
void Immersed_boundary::find_nearest_location_wall(double& x_min, double& y_min, double& z_min, double& d_min,
                                                   const double x, const double y, const double z, 
                                                   const int i, const int j, const int k)
{
    const double dx = grid->dx;
    const double dy = grid->dy;

    d_min = Constants::dbig;
    const int n  = 20;

    for (int ii = -n/2; ii < n/2+1; ++ii)
        for (int jj = -n/2; jj < n/2+1; ++jj)
        {
            const double xc = x + 2*ii/(double)n*dx;
            const double yc = y + 2*jj/(double)n*dy;
            const double zc = boundary_function<sw,dims>(xc, yc);
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

template<Immersed_boundary::IB_type sw, int dims>
void Immersed_boundary::find_interpolation_points(Ghost_cell& ghost_cell, 
                                                  const double* const restrict x, const double* const restrict y, const double* const restrict z,
                                                  const int i, const int j, const int k, const int n)
{
    double x_min, y_min, z_min, d_min;  // x, y, z location and distance to wall 

    // Minimal distance that a grid point used in the interpolation has to be
    // away from the IB, to prevent getting extreme extrapolated values in the ghost cells
    const double d_lim = 0.25 * std::min(std::min(grid->dx, grid->dy), grid->dz[k]);

    // Find the neighbouring grid points outside the IB
    for (int dk=-1; dk<3; ++dk)
        for (int dj=-2; dj<3; ++dj)
            for (int di=-2; di<3; ++di)
            {
                // Check if grid point is outside IB
                if (z[k+dk] > boundary_function<sw,dims>(x[i+di], y[j+dj]))
                {
                    // Calculate distance (d_min) of current grid point to the IB
                    find_nearest_location_wall<sw,dims>(x_min, y_min, z_min, d_min, x[i+di], y[j+dj], z[k+dk], i+di, j+dj, k+dk); 

                    // As described above; exclude grid points which are close to the IB
                    if (d_min > d_lim)
                    {
                        Neighbour tmp_neighbour = {i+di, j+dj, k+dk, abs_distance(x[i], x[i+di], y[j], y[j+dj], z[k], z[k+dk])};
                        ghost_cell.neighbours.push_back(tmp_neighbour);
                    }
                }
            }

    // Sort them on distance:
    std::sort(ghost_cell.neighbours.begin(), ghost_cell.neighbours.end(), compare_value);

    // Abort if there are insufficient neighbouring fluid points
    if (ghost_cell.neighbours.size() < n)
    {
        model->master->print_error("Only found %i of %i neighbour points\n", ghost_cell.neighbours.size(), n);
        throw 1;
    }

    // Only keep the N nearest neighbour fluid points
    ghost_cell.neighbours.erase(ghost_cell.neighbours.begin()+n, ghost_cell.neighbours.end());

    // Pre-calculate the coefficient for the inverse distance weighted interpolation
    ghost_cell.c_idw.resize(n+1); // n neighbour points + location on wall

    // Calculate distance interpolation points to image point
    for (int l=0; l<n; ++l)
        ghost_cell.c_idw[l] = abs_distance(ghost_cell.xI, x[ghost_cell.neighbours[l].i], 
                                           ghost_cell.yI, y[ghost_cell.neighbours[l].j],
                                           ghost_cell.zI, z[ghost_cell.neighbours[l].k]);

    ghost_cell.c_idw[n] = abs_distance(ghost_cell.xI, ghost_cell.xB, 
                                       ghost_cell.yI, ghost_cell.yB,
                                       ghost_cell.zI, ghost_cell.zB);

    // Save maximum distance
    const double max_distance = *std::max_element(std::begin(ghost_cell.c_idw), std::end(ghost_cell.c_idw));

    // Calculate coefficients
    ghost_cell.c_idw_sum = 0;
    for (int l=0; l<n+1; ++l)
    {
        ghost_cell.c_idw[l] = std::pow((max_distance - ghost_cell.c_idw[l]) / (max_distance * ghost_cell.c_idw[l]), 0.5);
        ghost_cell.c_idw_sum += ghost_cell.c_idw[l];
    }
}

//void Immersed_boundary::define_distance_matrix(Ghost_cell& ghost_cell, 
//                                               const double* const restrict x, const double* const restrict y, const double* const restrict z)
//{
//    const int n = 4;
//
//    // Temporary matrix
//    std::vector< std::vector<double> > tmp(n, std::vector<double>(n,0));
//    for (int i=0; i<n; ++i)
//        tmp[i][0] = 1;
//
//    tmp[0][1] = ghost_cell.xb;                   // x-location of given value on boundary
//    tmp[1][1] = x[ghost_cell.neighbours[0].i];   // x-location of fluid point #1
//    tmp[2][1] = x[ghost_cell.neighbours[1].i];   // x-location of fluid point #2
//    tmp[3][1] = x[ghost_cell.neighbours[2].i];   // x-location of fluid point #3
//
//    tmp[0][2] = ghost_cell.yb;                   // y-location of given value on boundary
//    tmp[1][2] = y[ghost_cell.neighbours[0].j];   // x-location of fluid point #1
//    tmp[2][2] = y[ghost_cell.neighbours[1].j];   // x-location of fluid point #2
//    tmp[3][2] = y[ghost_cell.neighbours[2].j];   // x-location of fluid point #3
//    
//    tmp[0][3] = ghost_cell.zb;                   // z-location of given value on boundary
//    tmp[1][3] = z[ghost_cell.neighbours[0].k];   // z-location of fluid point #1
//    tmp[2][3] = z[ghost_cell.neighbours[1].k];   // z-location of fluid point #1
//    tmp[3][3] = z[ghost_cell.neighbours[2].k];   // z-location of fluid point #1
//
//    // Resize matrix in ghost_cell to 4x4 
//    ghost_cell.B.resize(n, std::vector<double>(n));
//    
//    // Save the inverse of the matrix
//    inverse_mat_4x4(ghost_cell.B, tmp);
//}
   
template<Immersed_boundary::IB_type sw, int dims>
void Immersed_boundary::find_ghost_cells(std::vector<Ghost_cell>* ghost_cells, 
                                         const double* const restrict x, const double* const restrict y, const double* const restrict z)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    double d_wall;  // Distance to IB (m)

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                // 1. Check if this is a ghost cell, i.e. inside the IB, with a neighbour outside the IB
                if (is_ghost_cell<sw, dims>(x, y, z, i, j, k))
                {
                    Ghost_cell tmp_ghost = {i, j, k};

                    // 2. Find the closest location on the IB
                    find_nearest_location_wall<sw,dims>(tmp_ghost.xB, tmp_ghost.yB, tmp_ghost.zB, d_wall, x[i], y[j], z[k], i, j, k);

                    // 2.1 Location image point
                    tmp_ghost.xI = 2*tmp_ghost.xB - x[i];
                    tmp_ghost.yI = 2*tmp_ghost.yB - y[j];
                    tmp_ghost.zI = 2*tmp_ghost.zB - z[k];

                    // 3. Find the closest 3 grid points outside the IB
                    find_interpolation_points<sw,dims>(tmp_ghost, x, y, z, i, j, k, 5);

                    // 4. Define the matrix with distances to the boundary and each grid point used in the interpolation
                    //define_distance_matrix(tmp_ghost, x, y, z);

                    ghost_cells->push_back(tmp_ghost);
                }
            }
}

void Immersed_boundary::exec_stats(Mask *m)
{
    if (ib_type == None_type)
        return;
}

void Immersed_boundary::exec()
{
    if (ib_type == None_type)
        return;

    const int ii = 1;
    const double boundary_value = 0.;

    set_ghost_cells(ghost_cells_u, fields->u->data, boundary_value, grid->xh, grid->y, grid->z, 5, ii, grid->icells, grid->ijcells);
    set_ghost_cells(ghost_cells_v, fields->v->data, boundary_value, grid->x, grid->yh, grid->z, 5, ii, grid->icells, grid->ijcells);
    set_ghost_cells(ghost_cells_w, fields->w->data, boundary_value, grid->x, grid->y, grid->zh, 5, ii, grid->icells, grid->ijcells);
}
