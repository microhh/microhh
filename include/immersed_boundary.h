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

#ifndef IMMERSED_BOUNDARY
#define IMMERSED_BOUNDARY

#include "boundary.h"
#include "boundary_cyclic.h"

class Master;
class Input;
template<typename> class Grid;
template<typename> class Fields;

enum class IB_type {Disabled, DEM, User};

// Ghost cell info on staggered grid
template<typename TF>
struct Ghost_cells
{
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

        void init(Input&);
        void create();
        void exec_momentum();

    private:
        Master& master;
        Grid<TF>& grid;
        Fields<TF>& fields;
        Field3d_io<TF> field3d_io;
        Boundary_cyclic<TF> boundary_cyclic;

        IB_type sw_ib;

        int n_idw_points;       // Number of interpolation points in IDW interpolation
        Boundary_type sbcbot;   // Boundary type for scalars

        // IB input from DEM
        std::vector<TF> dem;

        // Structs with the ghost cell properties
        Ghost_cells<TF> ghost_u;
        Ghost_cells<TF> ghost_v;
        Ghost_cells<TF> ghost_w;
        Ghost_cells<TF> ghost_s;
};

#endif
