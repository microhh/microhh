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

#ifndef IMMERSED_BOUNDARY_H
#define IMMERSED_BOUNDARY_H

#include "boundary.h"
#include "boundary_cyclic.h"

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;
template<typename> class Cross;

enum class IB_type {Disabled, DEM, User};

// Ghost cell info on staggered grid
template<typename TF>
struct Ghost_cells
{
    int nghost;

    //
    // CPU
    //

    // Indices of IB ghost cells:
    std::vector<int> i;     // size = number of ghost cells
    std::vector<int> j;
    std::vector<int> k;

    // Nearest location of IB to ghost cell:
    std::vector<TF> xb;     // size = number of ghost cells
    std::vector<TF> yb;
    std::vector<TF> zb;

    // Location of interpolation point outside IB:
    std::vector<TF> xi;     // size = number of ghost cells
    std::vector<TF> yi;
    std::vector<TF> zi;

    std::vector<TF> di; // Distance ghost cell to interpolation point

    // Points outside IB used for IDW interpolation:
    std::vector<int> ip_i;     // size = number of ghost cells x n_idw_points
    std::vector<int> ip_j;
    std::vector<int> ip_k;
    std::vector<TF>  ip_d;  // Distance to interpolation point

    // Interpolation coefficients
    std::vector<TF> c_idw;     // size = number of ghost cells x n_idw_points
    std::vector<TF> c_idw_sum; // size = number of ghost cells

    // Spatially varying scalar (and momentum..) boundary conditions
    std::map<std::string, std::vector<TF>> sbot;
    std::vector<TF> mbot;
    
    //
    // GPU 
    //

    // Indices of IB ghost cells:
    int* i_g;
    int* j_g;
    int* k_g;

    // Nearest location of IB to ghost cell:
    TF* xb_g;
    TF* yb_g;
    TF* zb_g;

    // Location of interpolation point outside IB:
    TF* xi_g;
    TF* yi_g;
    TF* zi_g;

    TF* di_g;  // Distance ghost cell to interpolation point

    // Points outside IB used for IDW interpolation:
    int* ip_i_g;
    int* ip_j_g;
    int* ip_k_g;
    TF* ip_d_g;

    // Interpolation coefficients
    TF* c_idw_g;
    TF* c_idw_sum_g;

    // Spatially varying scalar (and momentum..) boundary conditions
    std::map<std::string, TF*> sbot_g;
    TF* mbot_g;
};

// Convenience struct to simplify sorting
template<typename TF>
struct Neighbour
{
    int i;
    int j;
    int k;
    TF distance;
};

template<typename TF>
class Immersed_boundary
{
    public:
        Immersed_boundary(Master&, Grid<TF>&, Fields<TF>&, Input&);
        ~Immersed_boundary();

        void init(Input&, Cross<TF>&);
        void create();

        void exec_momentum();
        void exec_scalars();

        void exec_cross(Cross<TF>&, unsigned long);

        bool has_mask(std::string);
        void get_mask(Stats<TF>&, std::string);

        void prepare_device();
        void clear_device();

        IB_type get_switch() const { return sw_ib; }


    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Field3d_io<TF> field3d_io;
        Boundary_cyclic<TF> boundary_cyclic;

        IB_type sw_ib;

        int n_idw_points;       // Number of interpolation points in IDW interpolation

        // Boundary conditions for scalars
        Boundary_type sbcbot;
        std::map<std::string, TF> sbc;
        std::vector<std::string> sbot_spatial_list;

        // IB input from DEM
        std::vector<TF> dem;
        std::vector<unsigned int> k_dem;

        // All ghost cell properties
        std::map<std::string, Ghost_cells<TF>> ghost;

        // Statistics
        std::vector<std::string> available_masks;

        // Cross-sections
        std::vector<std::string> crosslist;
};

#endif
