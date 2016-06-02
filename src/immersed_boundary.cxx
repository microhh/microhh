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
        return height * std::exp(pow(-(x-x0)/width, 2));    
    }

    double abs_distance(const double dx, const double dy, const double dz)
    {
        return std::pow(std::pow(dx, 2) + std::pow(dy, 2) + std::pow(dz, 2), 0.5);
    }

    double abs_distance(const double dx, const double dz)
    {
        return std::pow(std::pow(dx, 2) + std::pow(dz, 2), 0.5);
    }

    void determine_ghost_cells(std::vector<Ghost_cell>* ghost_cells, const double* const x, const double* const y, const double* const z,
                               const int n_neighbours, const double x0_hill, const double lz_hill, const double lx_hill, const double dx, const double dy,
                               const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int jj, const int kk)
    {
        const int ii = 1;

        // Determine the ghost cells 
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    
                    // u-location
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
                        std::sort(tmp_ghost.fluid_neighbours.begin(), tmp_ghost.fluid_neighbours.end(), [](Neighbour& a, Neighbour& b)
                            { return a.distance < b.distance; });
                   
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
                        std::vector< std::vector<double> > B_u(n+1, std::vector<double>(n+2,0));
                       
                        B_u[0][0] = 1;
                        B_u[1][0] = 1;
                        B_u[2][0] = 1;

                        B_u[0][1] = x_min;                                // x-location of given value on boundary 
                        B_u[1][1] = x[tmp_ghost.fluid_neighbours[0].i];   // x-location of fluid point #1
                        B_u[2][1] = x[tmp_ghost.fluid_neighbours[1].i];   // x-location of fluid point #2

                        B_u[0][2] = z_min;                                // z-location of given value on boundary
                        B_u[1][2] = z[tmp_ghost.fluid_neighbours[0].k];   // z-location of fluid point #1
                        B_u[2][2] = z[tmp_ghost.fluid_neighbours[1].k];   // z-location of fluid point #1
                    
                        // Add to vector with ghost_cells 
                        ghost_cells->push_back(tmp_ghost);
                    }
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
    const double lx_hill = grid->xsize / 8.;

    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    
    const int n_neighbours = 2;     // 2 for 2D linear

    // Determine ghost cells for u-component
    determine_ghost_cells(&ghost_cells_u, grid->xh, grid->y, grid->z, n_neighbours, x0_hill, lz_hill, lx_hill,
                          grid->dx, grid->dy, grid->istart, grid->iend, grid->jstart, grid->jend, grid->kstart, grid->kend, grid->icells, grid->ijcells);

    // Determine ghost cells for w-component
    determine_ghost_cells(&ghost_cells_w, grid->x, grid->y, grid->zh, n_neighbours, x0_hill, lz_hill, lx_hill,
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
//    if (sw_ib != "1")
//        return;
//
//    if (stats->get_switch() == "1")
//        check_no_penetration(&boundary_u_max, &boundary_w_max, fields->u->data, fields->w->data, IB_cells, grid->icells, grid->ijcells);
//
//    set_no_penetration_new(fields->ut->data, fields->vt->data, fields->wt->data,
//                           fields->u->data, fields->v->data, fields->w->data, 
//                           IB_cells, grid->icells, grid->ijcells);
//
//    set_no_slip_new(fields->ut->data, fields->wt->data,
//                    fields->u->data, fields->w->data,
//                    fields->rhoref, fields->rhorefh,
//                    grid->dzi, grid->dzhi,
//                    grid->dxi,
//                    fields->visc,
//                    IB_cells,
//                    grid->icells, grid->ijcells);
//
//    for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
//        set_scalar_new(it->second->data, fields->sp[it->first]->data,
//                       fields->u->data, fields->w->data,
//                       fields->rhoref, fields->rhorefh,
//                       grid->dzi, grid->dzhi,
//                       grid->dxi,
//                       fields->sp[it->first]->visc,
//                       IB_cells,
//                       grid->icells, grid->ijcells);

    //set_no_penetration(fields->ut->data, fields->wt->data,
    //                   fields->u->data, fields->w->data,
    //                   grid->istart, grid->iend,
    //                   grid->jstart, grid->jend,
    //                   grid->kstart, grid->kend,
    //                   grid->icells, grid->ijcells);

    //set_no_slip(fields->ut->data, fields->wt->data,
    //            fields->u->data, fields->w->data,
    //            fields->rhoref, fields->rhorefh,
    //            grid->dzi, grid->dzhi,
    //            grid->dxi,
    //            fields->visc,
    //            grid->istart, grid->iend,
    //            grid->jstart, grid->jend,
    //            grid->kstart, grid->kend,
    //            grid->icells, grid->ijcells);

    //for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    //    set_scalar(it->second->data, fields->sp[it->first]->data,
    //               fields->u->data, fields->w->data,
    //               fields->rhoref, fields->rhorefh,
    //               grid->dzi, grid->dzhi,
    //               grid->dxi,
    //               fields->sp[it->first]->visc,
    //               grid->istart, grid->iend,
    //               grid->jstart, grid->jend,
    //               grid->kstart, grid->kend,
    //               grid->icells, grid->ijcells);
}


