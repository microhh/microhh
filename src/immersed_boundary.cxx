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

Immersed_boundary::Immersed_boundary(Model* modelin, Input* inputin)
{
    model  = modelin;
    fields = model->fields;
    grid   = model->grid;

    inputin->get_item(&sw_ib, "boundary", "sw_ib", "", "0");
}

Immersed_boundary::~Immersed_boundary()
{
}

namespace
{
    double gaussian_bump(const double x, const double x0, const double height, const double width)
    {
        return height * std::exp(-pow((x-x0)/width, 2));    
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

    void determine_ghost_cells(std::vector<Ghost_cell>* ghost_cells, const double* const x, const double* const y, const double* const z,
                               const int n_neighbours, const double x0_hill, const double lz_hill, const double lx_hill, const double dx, const double dy,
                               const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int jj, const int kk)
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
                    
                    const double z_hill    = gaussian_bump(x[i],   x0_hill, lz_hill, lx_hill);
                    const double z_hill_im = gaussian_bump(x[i-1], x0_hill, lz_hill, lx_hill);
                    const double z_hill_ip = gaussian_bump(x[i+1], x0_hill, lz_hill, lx_hill);

                    // If this grid point is inside the hill, and one of its neighbours is not, this is a ghost cell
                    if (z[k] < z_hill && (z[k+1] >= z_hill || z[k] >= z_hill_im || z[k] >= z_hill_ip))
                    {
                        Ghost_cell tmp_ghost;
                        tmp_ghost.i = i;
                        tmp_ghost.j = j;
                        tmp_ghost.k = k;

                        // Find the two nearest neighbours outside of the hill
                        for (int ii = -1; ii < 2; ++ii)
                            for (int kk = -1; kk < 2; ++kk)
                                if(z[k+kk] > gaussian_bump(x[i+ii], x0_hill, lz_hill, lx_hill))
                                {
                                    Neighbour nb = {i+ii, j, k+kk, abs_distance(x[i]-x[i+ii], z[k]-z[k+kk])};
                                    tmp_ghost.fluid_neighbours.push_back(nb);
                                }

                        // Sort the neighbouring points on distance
                        //std::sort(tmp_ghost.fluid_neighbours.begin(), tmp_ghost.fluid_neighbours.end(), [](Neighbour& a, Neighbour& b)
                        //    { return a.distance < b.distance; });
                        std::sort(tmp_ghost.fluid_neighbours.begin(), tmp_ghost.fluid_neighbours.end(), compare_value);
                   
                        // Abort if there are insufficient neighbouring fluid points 
                        if (tmp_ghost.fluid_neighbours.size() < n_neighbours)
                        {
                            std::cout << "NEIN" << std::endl;
                            throw 1;
                        }

                        // Only keep the N nearest neighbour fluid points
                        tmp_ghost.fluid_neighbours.erase(tmp_ghost.fluid_neighbours.begin()+n_neighbours, tmp_ghost.fluid_neighbours.end());
                    
                        // Find closest point on boundary
                        double d_min = 1e12;
                        double x_min = 0;
                        double z_min = 0;
                        const int n  = 100;
                        for (int ii=-n/2; ii<n/2+1; ++ii)
                        {
                            const double xc = x[i] + 2*ii/(double)n*dx; 
                            const double zc = gaussian_bump(xc, x0_hill, lz_hill, lx_hill);
                            const double dc = abs_distance(x[i]-xc, z[k]-zc);

                            if (dc < d_min)
                            {
                                d_min = dc;
                                x_min = xc;
                                z_min = zc;
                            }
                        }

                        // Define the distance matrix for the interpolations
                        //std::vector< std::vector<double> > tmp_ghost.B_u(n+1, std::vector<double>(n+1,0));
                        tmp_ghost.B.resize(n+1, std::vector<double>(n+1));

                        tmp_matrix[0][1] = x_min;                                // x-location of given value on boundary 
                        tmp_matrix[1][1] = x[tmp_ghost.fluid_neighbours[0].i];   // x-location of fluid point #1
                        tmp_matrix[2][1] = x[tmp_ghost.fluid_neighbours[1].i];   // x-location of fluid point #2

                        tmp_matrix[0][2] = z_min;                                // z-location of given value on boundary
                        tmp_matrix[1][2] = z[tmp_ghost.fluid_neighbours[0].k];   // z-location of fluid point #1
                        tmp_matrix[2][2] = z[tmp_ghost.fluid_neighbours[1].k];   // z-location of fluid point #1
                  
                        // Save the inverse of the matrix
                        inverse_mat_3x3(tmp_ghost.B, tmp_matrix);

                        // Add to vector with ghost_cells 
                        ghost_cells->push_back(tmp_ghost);
                    }
                }
    }

    void set_ghost_cells(std::vector<Ghost_cell> ghost_cells, double* const restrict field, const double boundary_value,
                         const double* const restrict x, const double* const restrict z, const int ii, const int jj, const int kk)
    {
        std::vector<double> a(3);
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

            // Solve for ghost cell value
            matvecmul(a, it->B, known_values);

            double ghost_value = a[0] + a[1]*x[it->i] + a[2]*z[it->k];

            if(std::isnan(ghost_value))
            {
                std::cout << "shit -> fan" << std::endl;
                ghost_value = 0;        
            }

            const int ijk = it->i + it->j*jj + it->k*kk;
            field[ijk] = ghost_value;
        }
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

    // Tmp BvS: setting Gaussian hill
    const double x0_hill = grid->xsize / 2.;
    const double lz_hill = grid->zsize / 4.;
    const double lx_hill = grid->xsize / 16.;

    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    
    const int n_neighbours = 2;     // 2 for 2D linear

    // Determine ghost cells for u-component
    determine_ghost_cells(&ghost_cells_u, grid->xh, grid->y, grid->z, n_neighbours, x0_hill, lz_hill, lx_hill,
                          grid->dx, grid->dy, grid->istart, grid->iend, grid->jstart, grid->jend, grid->kstart, grid->kend, grid->icells, grid->ijcells);

    // Determine ghost cells for w-component
    determine_ghost_cells(&ghost_cells_w, grid->x, grid->y, grid->zh, n_neighbours, x0_hill, lz_hill, lx_hill,
                          grid->dx, grid->dy, grid->istart, grid->iend, grid->jstart, grid->jend, grid->kstart+1, grid->kend, grid->icells, grid->ijcells);

    // Determine ghost cells for scalar location
    determine_ghost_cells(&ghost_cells_s, grid->x, grid->y, grid->z, n_neighbours, x0_hill, lz_hill, lx_hill,
                          grid->dx, grid->dy, grid->istart, grid->iend, grid->jstart, grid->jend, grid->kstart, grid->kend, grid->icells, grid->ijcells);

    //for (std::vector<Ghost_cell>::iterator it=ghost_cells_u.begin(); it<ghost_cells_u.end(); ++it)
    //{
    //    std::cout << "------------------------" << std::endl;
    //    std::cout << it->i << " - " << it->j << " - " << it->k << std::endl;
    //    std::cout << "------------" << std::endl;
    //    for (std::vector<Neighbour>::iterator it2=it->fluid_neighbours.begin(); it2 < it->fluid_neighbours.end(); ++it2)
    //        std::cout << "     " << it2->i << " - " << it2->j << " - " << it2->k << " -- " << it2->distance << std::endl;
    //}

    //for (std::vector<Ghost_cell>::iterator it=ghost_cells_w.begin(); it<ghost_cells_w.end(); ++it)
    //{
    //    std::cout << "------------------------" << std::endl;
    //    std::cout << it->i << " - " << it->j << " - " << it->k << std::endl;
    //    std::cout << "------------" << std::endl;
    //    for (std::vector<Neighbour>::iterator it2=it->fluid_neighbours.begin(); it2 < it->fluid_neighbours.end(); ++it2)
    //        std::cout << "     " << it2->i << " - " << it2->j << " - " << it2->k << " -- " << it2->distance << std::endl;
    //}
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

    for (FieldMap::const_iterator it = fields->sp.begin(); it!=fields->sp.end(); it++)
        set_ghost_cells(ghost_cells_s, it->second->data, boundary_value, grid->x, grid->z, ii, grid->icells, grid->ijcells);
}


