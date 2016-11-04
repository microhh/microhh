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
#include <vector>
#include "master.h"
#include "model.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "stats.h"
#include "constants.h"
#include "finite_difference.h"
#include "immersed_boundary.h"
#include "input.h"

//#include <fstream>  // TMP BvS
//#include <sstream>

namespace
{
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
    void precalculate_idw(Ghost_cell& ghost_cell, const double* const restrict x, const double* const restrict y, const double* const restrict z)
    {
        const int n = ghost_cell.neighbours.size();

        // Pre-calculate the coefficient for the inverse distance weighted interpolation
        ghost_cell.c_idw.resize(n+1); // n neighbour points + location on wall

        // Calculate distance interpolation points to image point
        for (int l=0; l<n; ++l)
            ghost_cell.c_idw[l] = std::max(abs_distance(ghost_cell.xI, x[ghost_cell.neighbours[l].i],
                                                        ghost_cell.yI, y[ghost_cell.neighbours[l].j],
                                                        ghost_cell.zI, z[ghost_cell.neighbours[l].k]), Constants::dsmall);

        // Add the distance between boundary and image point 
        ghost_cell.c_idw[n] = std::max(abs_distance(ghost_cell.xI, ghost_cell.xB,
                                                    ghost_cell.yI, ghost_cell.yB,
                                                    ghost_cell.zI, ghost_cell.zB), Constants::dsmall);

        // Save maximum distance
        const double max_distance = *std::max_element(ghost_cell.c_idw.begin(), ghost_cell.c_idw.end());

        // Calculate coefficients
        ghost_cell.c_idw_sum = 0;
        for (int l=0; l<n+1; ++l)
        {
            ghost_cell.c_idw[l] = std::pow((max_distance - ghost_cell.c_idw[l]) / (max_distance * ghost_cell.c_idw[l]), 0.5);
            ghost_cell.c_idw_sum += ghost_cell.c_idw[l];
        }
    }

    // Set the ghost cells according to the chosen boundary conditions
    void set_ghost_cells(std::vector<Ghost_cell> const &ghost_cells, double* const restrict field, const double boundary_value,
                         const double* const restrict x, const double* const restrict y, const double* const restrict z,
                         const int nn, const int ii, const int jj, const int kk)
    {
        for (std::vector<Ghost_cell>::const_iterator it=ghost_cells.begin(); it<ghost_cells.end(); ++it)
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

    inputin->get_item(&sw_ib, "IB", "sw_ib", "", "0");

    // Get the IB type, and required parameters for that type
    if (sw_ib == "0")
    {
        ib_type = None_type;
    }
    else if (sw_ib == "user")
    {
        ib_type = User_type;
    }
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
        throw 1;
    }

    if (ib_type != None_type)
    {
        nerror += inputin->get_item(&n_idw, "IB", "n_idw", "");  // Number of grid points used in interpolation

        grid->set_minimum_ghost_cells(2, 2, 1);  // Set at leat two ghost cells in the horizontal

        // Get the scalar boundary conditions. For now, the momentum BC is hardcoded at a no-slip condition.
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

    if (ib_type == User_type)
    {
        std::string file_name;

        file_name = "u.ib_input";
        read_ghost_cells(ghost_cells_u, file_name, grid->xh, grid->y, grid->z);

        file_name = "v.ib_input";
        read_ghost_cells(ghost_cells_v, file_name, grid->x, grid->yh, grid->z);

        file_name = "w.ib_input";
        read_ghost_cells(ghost_cells_w, file_name, grid->x, grid->y, grid->zh);

        file_name = "s.ib_input";
        read_ghost_cells(ghost_cells_s, file_name, grid->x, grid->y, grid->z);
    }
    else 
    {
        find_ghost_cells(ghost_cells_u, grid->xh, grid->y,  grid->z );
        find_ghost_cells(ghost_cells_v, grid->x,  grid->yh, grid->z );
        find_ghost_cells(ghost_cells_w, grid->x,  grid->y,  grid->zh);
        find_ghost_cells(ghost_cells_s, grid->x,  grid->y,  grid->z );
    }

    model->master->print_message("Found %i IB[u] ghost cells \n", ghost_cells_u.size());
    model->master->print_message("Found %i IB[v] ghost cells \n", ghost_cells_v.size());
    model->master->print_message("Found %i IB[w] ghost cells \n", ghost_cells_w.size());
    model->master->print_message("Found %i IB[s] ghost cells \n", ghost_cells_s.size());
}

// Return the height of the IB as a function of x,y position
double Immersed_boundary::boundary_function(const double x, const double y)
{
    if (ib_type == Flat_type)
        return z_offset;
    else if (xy_dims == 1)
    {
        if (ib_type == Sine_type)
        {
            const double pi = std::acos((double)-1.);
            return z_offset + amplitude + amplitude * std::sin(2*pi*x/wavelength_x);
        }
        else if (ib_type == Gaus_type)
        {
            return z_offset + amplitude * std::exp(-pow((x-x0_hill)/(2*sigma_x_hill), 2));
        }
        else if (ib_type == Agnesi_type)
        {
            return z_offset + amplitude / (1. + pow((x-x0_hill)/sigma_x_hill, 2));
        }
    }
    else if (xy_dims == 2)
    {
        if (ib_type == Sine_type)
        {
            const double pi = std::acos((double)-1.);
            return z_offset + amplitude + amplitude * std::sin(2*pi*x/wavelength_x) * std::sin(2*pi*y/wavelength_y);
        }
        else if (ib_type == Gaus_type)
        {
            return z_offset + amplitude * std::exp(-pow((x-x0_hill)/(2*sigma_x_hill), 2))
                                        * std::exp(-pow((y-y0_hill)/(2*sigma_y_hill), 2));
        }
        else if (ib_type == Agnesi_type)
        {
            return z_offset + amplitude / (1. + pow((x-x0_hill)/sigma_x_hill, 2)
                                              + pow((y-y0_hill)/sigma_y_hill, 2));
        }
    }

    return 0;   // Otherwise the compiler complains
}

bool Immersed_boundary::is_ghost_cell(const double* const restrict x, const double* const restrict y, const double* const restrict z,
                                      const int i, const int j, const int k)
{
    if (z[k] <= boundary_function(x[i], y[j]))  // Inside IB
    {
        // Check if one of the neighbouring grid cells is outside the IB
        for (int dk=-1; dk<2; ++dk) // BvS Fix this
            if (z[k+dk] > boundary_function(x[i], y[j]))
                return true;
        for (int dj=-1; dj<2; ++dj)
            if (z[k] > boundary_function(x[i], y[j+dj]))
                return true;
        for (int di=-1; di<2; ++di)
            if (z[k] > boundary_function(x[i+di], y[j]))
                return true;
    }
    return false;
}

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
            const double zc = boundary_function(xc, yc);
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

void Immersed_boundary::find_interpolation_points(Ghost_cell& ghost_cell,
                                                  const double* const restrict x, const double* const restrict y, const double* const restrict z,
                                                  const int i, const int j, const int k, const int n_idw)
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
                if (z[k+dk] > boundary_function(x[i+di], y[j+dj]))
                {
                    // Calculate distance (d_min) of current grid point to the IB
                    find_nearest_location_wall(x_min, y_min, z_min, d_min, x[i+di], y[j+dj], z[k+dk], i+di, j+dj, k+dk);

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
    if (ghost_cell.neighbours.size() < n_idw)
    {
        model->master->print_error("Only found %i of n_idw=%i neighbour points\n", ghost_cell.neighbours.size(), n_idw);
        throw 1;
    }

    // Only keep the N nearest neighbour fluid points
    ghost_cell.neighbours.erase(ghost_cell.neighbours.begin()+n_idw, ghost_cell.neighbours.end());

    // Pre-calculate the inverse distance weighting coefficients
    precalculate_idw(ghost_cell, x, y, z);
}


void Immersed_boundary::read_ghost_cells(std::vector<Ghost_cell> &ghost_cells, std::string file_name,
                                         const double* const restrict x, const double* const restrict y, const double* const restrict z)
{
    bool is_optional = false;
    Data_map input;

    // Read the input data into a std::map< header_name, std::vecor<data> >
    model->input->read_data_file(&input, file_name, is_optional);

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

            Ghost_cell tmp_ghost = {i, j, k};

            // Location on boundary
            tmp_ghost.xB = input["xb"][n];
            tmp_ghost.yB = input["yb"][n];
            tmp_ghost.zB = input["zb"][n];

            // Location image point
            tmp_ghost.xI = 2*tmp_ghost.xB - x[i];
            tmp_ghost.yI = 2*tmp_ghost.yB - y[j];
            tmp_ghost.zI = 2*tmp_ghost.zB - z[k];

            // Neighbours
            const int n_neighbours = input["nn"][n];

            for (int nn=0; nn<n_neighbours; ++nn)
            {
                // Note BvS: the cast to long long is only necessary as intel currently doesn't have the full c++11 set of to_string() implemented
                const int in = input["i"+std::to_string(static_cast<long long>(nn))][n] + grid->istart - mpioffsx;
                const int jn = input["j"+std::to_string(static_cast<long long>(nn))][n] + grid->jstart - mpioffsy;
                const int kn = input["k"+std::to_string(static_cast<long long>(nn))][n] + grid->kstart;

                Neighbour tmp_neighbour = {in, jn, kn, abs_distance(x[i], x[in], y[j], y[jn], z[k], z[kn])};
                tmp_ghost.neighbours.push_back(tmp_neighbour);
            }

            precalculate_idw(tmp_ghost, x, y, z);

            ghost_cells.push_back(tmp_ghost);
        }
    }
}

void Immersed_boundary::find_ghost_cells(std::vector<Ghost_cell> &ghost_cells,
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
                if (is_ghost_cell(x, y, z, i, j, k))
                {
                    Ghost_cell tmp_ghost = {i, j, k};

                    // 2. Find the closest location on the IB
                    find_nearest_location_wall(tmp_ghost.xB, tmp_ghost.yB, tmp_ghost.zB, d_wall, x[i], y[j], z[k], i, j, k);

                    // 2.1 Location image point
                    tmp_ghost.xI = 2*tmp_ghost.xB - x[i];
                    tmp_ghost.yI = 2*tmp_ghost.yB - y[j];
                    tmp_ghost.zI = 2*tmp_ghost.zB - z[k];

                    // 3. Find the closest `n_idw` grid points outside the IB
                    find_interpolation_points(tmp_ghost, x, y, z, i, j, k, n_idw);

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

void Immersed_boundary::calc_mask(double* const restrict mask, double* const restrict maskh, double* const restrict maskbot,
                                  int* const restrict nmask, int* const restrict nmaskh, int* const restrict nmaskbot,
                                  const double* const restrict x, const double* const restrict y,
                                  const double* const restrict z, const double* const restrict zh)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

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

                const double zb = boundary_function(x[i], y[j]);
                const int is_not_ib    = z [k] > zb;
                const int is_not_ib_h  = zh[k] > zb;

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
    if ((ib_type != Sine_type) && (ib_type != Gaus_type) && (ib_type != Agnesi_type) && (ib_type != Flat_type))
    {
        model->master->print_error("Get_mask() not yet implemented for chosen IB type\n");
        throw 1;
    }

    calc_mask(mfield->data, mfieldh->data, mfieldh->databot, stats->nmask, stats->nmaskh, &stats->nmaskbot,
              grid->x, grid->y, grid->z, grid->zh);
}

void Immersed_boundary::exec()
{
    if (ib_type == None_type)
        return;

    const int ii = 1;
    const double boundary_value = 0.;

    // Set the immersed boundary ghost cells to enforce a no-slip BC for the velocity components
    set_ghost_cells(ghost_cells_u, fields->u->data, boundary_value, grid->xh, grid->y, grid->z, 4, ii, grid->icells, grid->ijcells);
    set_ghost_cells(ghost_cells_v, fields->v->data, boundary_value, grid->x, grid->yh, grid->z, 4, ii, grid->icells, grid->ijcells);
    set_ghost_cells(ghost_cells_w, fields->w->data, boundary_value, grid->x, grid->y, grid->zh, 4, ii, grid->icells, grid->ijcells);

    grid->boundary_cyclic(fields->u->data);
    grid->boundary_cyclic(fields->v->data);
    grid->boundary_cyclic(fields->w->data);

    // Set the ghost cells for scalars, depending on their BC
    for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
        set_ghost_cells(ghost_cells_s, fields->ap[it->first]->data, sbc[it->first]->bot, grid->x, grid->y, grid->z, 4, ii, grid->icells, grid->ijcells);
        grid->boundary_cyclic(fields->ap[it->first]->data);
    }
}





// "Archive"....
/*
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

    void inverse_mat_4x4(std::vector<std::vector<double> > &result, std::vector<std::vector<double> > const &input)
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

void Immersed_boundary::define_distance_matrix(Ghost_cell& ghost_cell,
                                               const double* const restrict x, const double* const restrict y, const double* const restrict z)
{
    const int n = 4;

    // Temporary matrix
    std::vector< std::vector<double> > tmp(n, std::vector<double>(n,0));
    for (int i=0; i<n; ++i)
        tmp[i][0] = 1;

    tmp[0][1] = ghost_cell.xb;                   // x-location of given value on boundary
    tmp[1][1] = x[ghost_cell.neighbours[0].i];   // x-location of fluid point #1
    tmp[2][1] = x[ghost_cell.neighbours[1].i];   // x-location of fluid point #2
    tmp[3][1] = x[ghost_cell.neighbours[2].i];   // x-location of fluid point #3

    tmp[0][2] = ghost_cell.yb;                   // y-location of given value on boundary
    tmp[1][2] = y[ghost_cell.neighbours[0].j];   // x-location of fluid point #1
    tmp[2][2] = y[ghost_cell.neighbours[1].j];   // x-location of fluid point #2
    tmp[3][2] = y[ghost_cell.neighbours[2].j];   // x-location of fluid point #3

    tmp[0][3] = ghost_cell.zb;                   // z-location of given value on boundary
    tmp[1][3] = z[ghost_cell.neighbours[0].k];   // z-location of fluid point #1
    tmp[2][3] = z[ghost_cell.neighbours[1].k];   // z-location of fluid point #1
    tmp[3][3] = z[ghost_cell.neighbours[2].k];   // z-location of fluid point #1

    // Resize matrix in ghost_cell to 4x4
    ghost_cell.B.resize(n, std::vector<double>(n));

    // Save the inverse of the matrix
    inverse_mat_4x4(ghost_cell.B, tmp);
}


 */

