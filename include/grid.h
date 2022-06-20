/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
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

#ifndef GRID_H
#define GRID_H

#ifdef USEMPI
#include <mpi.h>
#endif

#include <string>
#include <vector>
#include <array>
#include "defines.h"
#include "transpose.h"
#include "cuda_buffer.h"

class Master;
class Input;
class Netcdf_handle;

enum class Grid_order { Second, Fourth };

template<typename TF>
struct Grid_data
{
    int itot; // Total number of grid cells in the x-direction.
    int jtot; // Total number of grid cells in the y-direction.
    int ktot; // Total number of grid cells in the z-direction.
    int ntot; // Total number of grid cells.

    int imax; // Number of grid cells in the x-direction for one process.
    int jmax; // Number of grid cells in the y-direction for one process.
    int kmax; // Number of grid cells in the z-direction for one process.
    int nmax; // Number of grid cells for one process.

    int iblock; // Number of grid cells in the x-direction for calculation blocks in transposes.
    int jblock; // Number of grid cells in the y-direction for calculation blocks in transposes.
    int kblock; // Number of grid cells in the z-direction for calculation blocks in transposes.

    int igc; // Number of ghost cells in the x-direction.
    int jgc; // Number of ghost cells in the y-direction.
    int kgc; // Number of ghost cells in the z-direction.

    int icells;  // Number of grid cells in the x-direction including ghost cells for one process.
    int jcells;  // Number of grid cells in the y-direction including ghost cells for one process.
    int ijcells; // Number of grid cells in the xy-plane including ghost cells for one process.
    int kcells;  // Number of grid cells in the z-direction including ghost cells for one process.
    int ncells;  // Total number of grid cells for one process including ghost cells.
    int istart;  // Index of the first grid point in the x-direction.
    int jstart;  // Index of the first grid point in the y-direction.
    int kstart;  // Index of the first grid point in the z-direction.
    int iend;    // Index of the last gridpoint+1 in the x-direction.
    int jend;    // Index of the last gridpoint+1 in the y-direction.
    int kend;    // Index of the last gridpoint+1 in the z-direction.

    TF xsize; // Size of the domain in the x-direction.
    TF ysize; // Size of the domain in the y-direction.
    TF zsize; // Size of the domain in the z-direction.

    TF dx;    // Distance between the center of two grid cell in the x-direction.
    TF dy;    // Distance between the center of two grid cell in the y-direction.
    TF dxi;   // Reciprocal of dx.
    TF dyi;   // Reciprocal of dy.

    std::vector<TF> dz;    // Distance between the faces of two grid cells in the z-direction.
    std::vector<TF> dzh;   // Distance between the centers of two grid cells in the z-direction.
    std::vector<TF> dzi;   // Reciprocal of dz.
    std::vector<TF> dzhi;  // Reciprocal of dzh.
    std::vector<TF> dzi4;  // Fourth order gradient of the distance between cell faces to be used in 4th-order schemes.
    std::vector<TF> dzhi4; // Fourth order gradient of the distance between cell centers to be used in 4th-order schemes.

    TF dzhi4bot;
    TF dzhi4top;

    std::vector<TF> x;  // Grid coordinate of cell center in x-direction.
    std::vector<TF> y;  // Grid coordinate of cell center in y-direction.
    std::vector<TF> z;  // Grid coordinate of cell center in z-direction.
    std::vector<TF> xh; // Grid coordinate of cell faces in x-direction.
    std::vector<TF> yh; // Grid coordinate of cell faces in x-direction.
    std::vector<TF> zh; // Grid coordinate of cell faces in x-direction.

    // GPU fields and settings
    int ithread_block; // Number of grid cells in the x-direction for GPU thread block.
    int jthread_block; // Number of grid cells in the y-direction for GPU thread block.

    cuda_buffer<TF> x_g;
    cuda_buffer<TF> y_g;
    cuda_buffer<TF> z_g;
    cuda_buffer<TF> zh_g;
    cuda_buffer<TF> dz_g;
    cuda_buffer<TF> dzh_g;
    cuda_buffer<TF> dzi_g;
    cuda_buffer<TF> dzhi_g;
    cuda_buffer<TF> dzi4_g;
    cuda_buffer<TF> dzhi4_g;

    const std::array<int,3> uloc  = {{1,0,0}}; // Location of the u-velocity on the staggered grid
    const std::array<int,3> vloc  = {{0,1,0}}; // Location of the v-velocity on the staggered grid
    const std::array<int,3> wloc  = {{0,0,1}}; // Location of the w-velocity on the staggered grid
    const std::array<int,3> sloc  = {{0,0,0}}; // Location of the cell center on the staggered grid
    const std::array<int,3> uwloc = {{1,0,1}}; // Location of the u-flux on the staggered grid
    const std::array<int,3> vwloc = {{0,1,1}}; // Location of the v-flux on the staggered grid
};

/**
 * Class for the grid settings and operators.
 * This class contains the grid properties, such as dimensions and resolution.
 * The public funcions in this class contain operations that are called from many routines
 * in order to interpolate, transpose and save data. The MPI operations that work over multiple
 * processes on the entire grid are contained in this class.
 */
template<typename TF>
class Grid
{
    public:
        Grid(Master&, Input&); // Constructor of the grid class.
        ~Grid();               // Destructor of the grid class.

        void init();              // Initialization of the grid arrays.
        void create(Netcdf_handle&); // Creation of the grid data.
        void save();              // Saves grid data to file.
        void load();              // Loads grid data to file.

        TF utrans; // Galilean transformation velocity in x-direction.
        TF vtrans; // Galilean transformation velocity in y-direction.

        const Grid_data<TF>& get_grid_data();
        Grid_order get_spatial_order() const { return spatial_order; }

        void set_minimum_ghost_cells(int, int, int);

        // MPI functions
        void init_mpi(); // Creates the MPI data types used in grid operations.
        void exit_mpi(); // Destructs the MPI data types used in grid operations.

        // interpolation functions
        void interpolate_2nd(TF*, const TF*, const int[3], const int[3]); // Second order interpolation
        void interpolate_4th(TF*, const TF*, const int[3], const int[3]); // Fourth order interpolation

        // GPU functions
        void prepare_device(); // Load the arrays onto the GPU
        void clear_device();   // Deallocate the arrays onto the GPU

    private:
        Master& master; // Reference to master class.
        Transpose<TF> transpose;

        Grid_order spatial_order; // Default spatial order of the operators to be used on this grid.

        bool mpitypes;  // Boolean to check whether MPI datatypes are created.

        void calculate(); // Computation of dimensions, faces and ghost cells.
        void check_ghost_cells(); // Check whether slice thickness is at least equal to number of ghost cells.

        void save_grid(); // Save the grid properties.
        void load_grid(); // Load the grid properties.

        Grid_data<TF> gd;

        #ifdef USEMPI
        MPI_Datatype subi; // MPI datatype containing a subset of the entire x-axis.
        MPI_Datatype subj; // MPI datatype containing a subset of the entire y-axis.
        #endif
};
#endif
