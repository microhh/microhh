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
#include <cmath>
#include <algorithm>
#include "master.h"
#include "model.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "stats.h"
#include "finite_difference.h"
#include "immersed_boundary.h"

#include <fstream>  // TMP BvS

Immersed_boundary::Immersed_boundary(Model* modelin, Input* inputin)
{
    model  = modelin;
    fields = model->fields;
    grid   = model->grid;

    inputin->get_item(&sw_ib,   "boundary", "sw_ib",   "", "0");
    inputin->get_item(&x0_hill, "boundary", "x0_hill", "",  0);
    inputin->get_item(&lz_hill, "boundary", "lz_hill", "",  0);
    inputin->get_item(&lx_hill, "boundary", "lx_hill", "",  0);
}

Immersed_boundary::~Immersed_boundary()
{
}

namespace
{
    double gaussian_bump(const double x, const double x0, const double height, const double width)
    {
        return height * std::exp(-pow((x-x0)/(2*width), 2));
    }

    double gaussian_bump_i(const double y, const double x0, const double height, const double width, const int mode)
    {
        const double diff = 2*width * pow(std::log(height/y), 2);
        if (mode == 0)
            return x0-diff;
        else
            return x0+diff;
    }

    double abs_distance(const double dx, const double dy, const double dz)
    {
        return std::pow(std::pow(dx, 2) + std::pow(dy, 2) + std::pow(dz, 2), 0.5);
    }

    double abs_distance(const double dx, const double dz)
    {
        return std::pow(std::pow(dx, 2) + std::pow(dz, 2), 0.5);
    }

    void inverse_mat_3x3(std::vector<std::vector<double> >& result, std::vector<std::vector<double> > input)
    {
        const double det = +input[0][0]*(input[1][1]*input[2][2]-input[2][1]*input[1][2])
                           -input[0][1]*(input[1][0]*input[2][2]-input[1][2]*input[2][0])
                           +input[0][2]*(input[1][0]*input[2][1]-input[1][1]*input[2][0]);
        const double invdet = 1./det;
        result[0][0]     =  (input[1][1]*input[2][2]-input[2][1]*input[1][2])*invdet;
        result[0][1]     = -(input[0][1]*input[2][2]-input[0][2]*input[2][1])*invdet;
        result[0][2]     =  (input[0][1]*input[1][2]-input[0][2]*input[1][1])*invdet;
        result[1][0]     = -(input[1][0]*input[2][2]-input[1][2]*input[2][0])*invdet;
        result[1][1]     =  (input[0][0]*input[2][2]-input[0][2]*input[2][0])*invdet;
        result[1][2]     = -(input[0][0]*input[1][2]-input[1][0]*input[0][2])*invdet;
        result[2][0]     =  (input[1][0]*input[2][1]-input[2][0]*input[1][1])*invdet;
        result[2][1]     = -(input[0][0]*input[2][1]-input[2][0]*input[0][1])*invdet;
        result[2][2]     =  (input[0][0]*input[1][1]-input[1][0]*input[0][1])*invdet;
    }

    // Help function for sorting std::vector
    bool compare_value(const Neighbour& a, const Neighbour& b)
    {
        return a.distance < b.distance;
    }

    void matvecmul(std::vector<double>& vec_out, std::vector<std::vector<double> > mat_in, std::vector<double> vec_in)
    {
        const int size = vec_in.size();
        for (int i=0; i<size; ++i)
        {
            vec_out[i] = 0;
            for (int j=0; j<size; ++j)
                vec_out[i] += mat_in[i][j] * vec_in[j];
        }
    }

    bool is_ghost_cell(const double* const x, const double* const z,
                       const int i, const int k, 
                       const double x0, const double z0, const double L)
    {
        // Gaussian hill; x0 = center hill, z0 = height hill, L = std.dev hill
        if (false)
        {
            if ( z[k] < gaussian_bump(x[i], x0, z0, L) && (
                        z[k+1] >= gaussian_bump(x[i  ], x0, z0, L) ||
                        z[k  ] >= gaussian_bump(x[i-1], x0, z0, L) ||
                        z[k  ] >= gaussian_bump(x[i+1], x0, z0, L) ) )
                return true;
            else
                return false;
        }   

        // Circle; {x0,z0} = center, L = radius
        if(true)
        {
            const double r2 = std::pow(L, 2);
            if( pow(x[i  ]-x0, 2) + pow(z[k  ]-z0, 2) < r2 && (
                pow(x[i+1]-x0, 2) + pow(z[k  ]-z0, 2) > r2 ||
                pow(x[i-1]-x0, 2) + pow(z[k  ]-z0, 2) > r2 ||
                pow(x[i  ]-x0, 2) + pow(z[k+1]-z0, 2) > r2 ||
                pow(x[i  ]-x0, 2) + pow(z[k-1]-z0, 2) > r2    ) )
                return true;
            else
                return false;  
        }
    }

    void find_min_distance_wall(double& x_min, double& z_min,
                                const double x, const double z, const double dx,
                                const double x0, const double z0, const double L)
    {
        double d_min = 1e12;
        const int n  = 100;

        // Gaussian hill; x0 = center hill, z0 = height hill, L = std.dev hill
        if (false)
        {
            for (int ii=-n/2; ii<n/2+1; ++ii)
            {
                const double xc = x + 2*ii/(double)n*dx;
                const double zc = gaussian_bump(xc, x0, z0, L);
                const double dc = abs_distance(x-xc, z-zc);

                if (dc < d_min)
                {
                    d_min = dc;
                    x_min = xc;
                    z_min = zc;
                }
            }
        }

        // Circle; {x0,z0} = center, L = radius
        if(true)
        {
            const double alpha = std::atan2(z-z0, x-x0);
            const double dz    = std::sin(alpha) * L;
            const double dx    = std::cos(alpha) * L;

            x_min = x0 + std::cos(alpha) * L;
            z_min = z0 + std::sin(alpha) * L;
        } 
    }

    void find_neighbour_points(Ghost_cell& tmp_ghost, const int n_neighbours,
                               const double* const x, const double* const z,
                               const double dx, const int i, const int j, const int k,
                               const double x0, const double z0, const double L)
    {
        // Gaussian hill; x0 = center hill, z0 = height hill, L = std.dev hill
        if (false)
        {
            for (int ii = -2; ii < 3; ++ii)
            {
                const double zb = gaussian_bump(x[i+ii], x0, z0, L);
                for (int kk = -2; kk < 3; ++kk)
                    if(z[k+kk] > zb)
                    {
                        // Check distance of current grid point to the wall
                        int mode;
                        x[i+ii] < x0 ? mode = 0 : mode = 1;

                        const double dz_wall = std::abs(z[k+kk]-gaussian_bump  (x[i+ii], x0, z0, L));
                        const double dx_wall = std::abs(x[i+ii]-gaussian_bump_i(z[k+kk], x0, z0, L, mode));
                        const double d_wall = 0.5 * abs_distance(dx_wall, dz_wall);

                        if(std::min(dz_wall, dx_wall) > 0.1*dx)
                        {
                            Neighbour nb = {i+ii, j, k+kk, abs_distance(tmp_ghost.xb-x[i+ii], tmp_ghost.zb-z[k+kk])};
                            tmp_ghost.fluid_neighbours.push_back(nb);
                        }
                    }
            }
        }

        // Circle; {x0,z0} = center, L = radius
        if (true)
        {
            for (int ii = -2; ii < 3; ++ii)
                for (int kk = -2; kk < 3; ++kk)
                {
                    const double r = std::pow(std::pow(x[i+ii]-x0, 2) + std::pow(z[k+kk]-z0, 2), 0.5);
                    if (r > L)
                    {
                        double x_min, z_min;
                        find_min_distance_wall(x_min, z_min, x[i+ii], z[k+kk], dx, x0, z0, L);
                        const double dist = abs_distance(x[i+ii]-x_min, z[k+kk]-z_min);

                        if (dist > 0.5 * dx)
                        {
                            Neighbour nb = {i+ii, j, k+kk, dist};
                            tmp_ghost.fluid_neighbours.push_back(nb);
                        }
                    }
                }
        }

        // Sort the neighbouring points on distance
        //std::sort(tmp_ghost.fluid_neighbours.begin(), tmp_ghost.fluid_neighbours.end(), [](Neighbour& a, Neighbour& b)
        //    { return a.distance < b.distance; }); // Doesn't work with CLANG?!
        std::sort(tmp_ghost.fluid_neighbours.begin(), tmp_ghost.fluid_neighbours.end(), compare_value);

        // Abort if there are insufficient neighbouring fluid points
        if (tmp_ghost.fluid_neighbours.size() < n_neighbours)
        {
            std::cout << "NEIN" << std::endl;
            throw 1;
        }

        // Only keep the N nearest neighbour fluid points
        tmp_ghost.fluid_neighbours.erase(tmp_ghost.fluid_neighbours.begin()+n_neighbours, tmp_ghost.fluid_neighbours.end());
    }

    void determine_ghost_cells(std::vector<Ghost_cell>* ghost_cells, std::vector<Cell>* boundary_cells,
                               const double* const x, const double* const y, const double* const z,
                               const int n_neighbours, const double x0_hill, const double lz_hill, const double lx_hill,
                               const double dx, const double dy,
                               const int istart, const int iend, const int jstart, const int jend,
                               const int kstart, const int kend, const int jj, const int kk)
    {
        const int ii = 1;

        std::vector< std::vector<double> > tmp_matrix(n_neighbours+1, std::vector<double>(n_neighbours+1,0));
        tmp_matrix[0][0] = 1;
        tmp_matrix[1][0] = 1;
        tmp_matrix[2][0] = 1;

        // Determine the ghost cells
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    if (is_ghost_cell(x, z, i, k, x0_hill, lz_hill, lx_hill))
                    {
                        Ghost_cell tmp_ghost;
                        tmp_ghost.i = i;
                        tmp_ghost.j = j;
                        tmp_ghost.k = k;

                        // Find the minimum distance from ghost cell to the wall
                        find_min_distance_wall(tmp_ghost.xb, tmp_ghost.zb, x[i], z[k], dx, x0_hill, lz_hill, lx_hill);

                        // Find the nearest neighbours outside of the hill
                        find_neighbour_points(tmp_ghost, n_neighbours, x, z, dx, i, j, k, x0_hill, lz_hill, lx_hill);

                        // Define the distance matrix for the interpolations
                        tmp_ghost.B.resize(n_neighbours+1, std::vector<double>(n_neighbours+1));

                        tmp_matrix[0][1] = tmp_ghost.xb;                         // x-location of given value on boundary
                        tmp_matrix[1][1] = x[tmp_ghost.fluid_neighbours[0].i];   // x-location of fluid point #1
                        tmp_matrix[2][1] = x[tmp_ghost.fluid_neighbours[1].i];   // x-location of fluid point #2

                        tmp_matrix[0][2] = tmp_ghost.zb;                         // z-location of given value on boundary
                        tmp_matrix[1][2] = z[tmp_ghost.fluid_neighbours[0].k];   // z-location of fluid point #1
                        tmp_matrix[2][2] = z[tmp_ghost.fluid_neighbours[1].k];   // z-location of fluid point #1

                        // Save the inverse of the matrix
                        inverse_mat_3x3(tmp_ghost.B, tmp_matrix);

                        // Add to vector with ghost_cells
                        ghost_cells->push_back(tmp_ghost);
                    }
                    //else if(z[k] < z_hill)
                    //{
                    //    Cell tmp = {i, j, k};
                    //    boundary_cells->push_back(tmp);
                    //}
                }
    }

    void set_ghost_cells(std::vector<Ghost_cell> ghost_cells, double* const restrict field, const double boundary_value,
                         const double* const restrict x, const double* const restrict z, const int ii, const int jj, const int kk)
    {
        std::vector<double> coefficients(3);
        std::vector<double> known_values(3);

        for (std::vector<Ghost_cell>::iterator it=ghost_cells.begin(); it<ghost_cells.end(); ++it)
        {
            // Values at the nearest fluid points
            const double v1 = field[it->fluid_neighbours[0].i + it->fluid_neighbours[0].j*jj + it->fluid_neighbours[0].k*kk];
            const double v2 = field[it->fluid_neighbours[1].i + it->fluid_neighbours[1].j*jj + it->fluid_neighbours[1].k*kk];

            // Populate vector with values
            known_values[0] = boundary_value;
            known_values[1] = v1;
            known_values[2] = v2;

            // Solve coefficients
            matvecmul(coefficients, it->B, known_values);

            // Location of image point
            const double xi = 2*it->xb - x[it->i];
            const double zi = 2*it->zb - z[it->k];

            // Value at image point
            double image_value = coefficients[0] + coefficients[1]*xi + coefficients[2]*zi;

            // Set ghost cell value
            const int ijk = it->i + it->j*jj + it->k*kk;
            field[ijk] = 2*boundary_value - image_value;
        }
    }

    void set_boundary_cells(std::vector<Cell> boundary_cells, double* const restrict tendency,
                            const int ii, const int jj, const int kk)
    {

        for (std::vector<Cell>::iterator it=boundary_cells.begin(); it<boundary_cells.end(); ++it)
        {
            const int ijk = it->i + it->j*jj + it->k*kk;
            tendency[ijk] = 0;
        }
    }

    void print_debug(std::vector<Ghost_cell> ghost_cells)
    {
        std::ofstream debugfile;
        debugfile.open("debug.txt");

        for (std::vector<Ghost_cell>::iterator it=ghost_cells.begin(); it<ghost_cells.end(); ++it)
        {
            debugfile << it->i                          << ", " <<
                         it->k                          << ", " <<
                         it->xb                         << ", " <<
                         it->zb                         << ", " <<
                         it->fluid_neighbours[0].i      << ", " <<
                         it->fluid_neighbours[0].k      << ", " <<
                         it->fluid_neighbours[1].i      << ", " <<
                         it->fluid_neighbours[1].k      << std::endl;
        }

        debugfile.close();
    }
}

void Immersed_boundary::init()
{
    stats = model->stats;

    if (sw_ib != "1")
        return;
}

void Immersed_boundary::create()
{
    if (sw_ib != "1")
        return;

    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const int n_neighbours = 2;     // 2 for 2D linear

    // Determine ghost cells for u-component
    determine_ghost_cells(&ghost_cells_u, &boundary_cells_u, grid->xh, grid->y, grid->z, n_neighbours, x0_hill, lz_hill, lx_hill,
                          grid->dx, grid->dy, grid->istart, grid->iend, grid->jstart, grid->jend, grid->kstart, grid->kend, grid->icells, grid->ijcells);

    // Determine ghost cells for w-component
    determine_ghost_cells(&ghost_cells_w, &boundary_cells_w, grid->x, grid->y, grid->zh, n_neighbours, x0_hill, lz_hill, lx_hill,
                          grid->dx, grid->dy, grid->istart, grid->iend, grid->jstart, grid->jend, grid->kstart+1, grid->kend, grid->icells, grid->ijcells);

    // Determine ghost cells for scalar location
    determine_ghost_cells(&ghost_cells_s, &boundary_cells_s, grid->x, grid->y, grid->z, n_neighbours, x0_hill, lz_hill, lx_hill,
                          grid->dx, grid->dy, grid->istart, grid->iend, grid->jstart, grid->jend, grid->kstart, grid->kend, grid->icells, grid->ijcells);

    model->master->print_message("Found %i IB[u] ghost cells \n",ghost_cells_u.size());
    model->master->print_message("Found %i IB[w] ghost cells \n",ghost_cells_w.size());
    model->master->print_message("Found %i IB[s] ghost cells \n",ghost_cells_s.size());

    print_debug(ghost_cells_s);
}

void Immersed_boundary::exec_stats(Mask *m)
{
//    if (sw_ib != "1")
//        return;
//
//    // Temporary hack (BvS) to get wall velocities in statistics, doesn't account for mask!!
//    m->tseries["IB_u_max"].data = boundary_u_max;
//    m->tseries["IB_w_max"].data = boundary_w_max;
}

void Immersed_boundary::exec()
{
    if (sw_ib != "1")
        return;

    const double boundary_value = 0;
    const int ii = 1;

    set_ghost_cells(ghost_cells_u, fields->u->data, boundary_value, grid->xh, grid->z, ii, grid->icells, grid->ijcells);
    set_ghost_cells(ghost_cells_w, fields->w->data, boundary_value, grid->x, grid->zh, ii, grid->icells, grid->ijcells);

    //const double sboundary_value = 10;
    //for (FieldMap::const_iterator it = fields->sp.begin(); it!=fields->sp.end(); it++)
    //    set_ghost_cells(ghost_cells_s, it->second->data, sboundary_value, grid->x, grid->z, ii, grid->icells, grid->ijcells);
}

void Immersed_boundary::exec_tend()
{
    if (sw_ib != "1")
        return;

    const int ii = 1;

    set_boundary_cells(boundary_cells_u, fields->ut->data, ii, grid->icells, grid->ijcells);
    set_boundary_cells(boundary_cells_w, fields->wt->data, ii, grid->icells, grid->ijcells);

    for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
        set_boundary_cells(boundary_cells_s, it->second->data, ii, grid->icells, grid->ijcells);
}

