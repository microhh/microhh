/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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

#ifndef SOIL_GRID
#define SOIL_GRID

#include <vector>
#include "defines.h"

class Master;
class Input;
class Netcdf_handle;
template<typename> class Grid;

template<typename TF>
struct Soil_grid_data
{
    bool is_enabled;  // Is the soil grid active/initialised?

    int ktot;   // Total number of full level grid cells in the z-direction.
    int kmax;   // Number of full level grid cells in the z-direction for one process.
    int kmaxh;  // Number of half level grid point in the z-direction for one process.
    int kgc;    // Number of ghost cells in the z-direction.

    int kcells;   // Number of grid cells in the z-direction including ghost cells for one process.
    int kcellsh;  // Number of grid cells in the z-direction including ghost cells for one process.

    int ncells;   // Total number of grid cells for one process including ghost cells.
    int ncellsh;  // Total number of grid cells for one process including ghost cells.

    int kstart;  // Index of the first grid point in the z-direction.
    int kend;    // Index of the last gridpoint+1 in the z-direction.

    TF zsize; // Size of the domain in the z-direction.

    std::vector<TF> dz;    // Distance between the faces of two grid cells in the z-direction.
    std::vector<TF> dzh;   // Distance between the centers of two grid cells in the z-direction.
    std::vector<TF> dzi;   // Reciprocal of dz.
    std::vector<TF> dzhi;  // Reciprocal of dzh.

    std::vector<TF> z;  // Grid coordinate of cell center in z-direction.
    std::vector<TF> zh; // Grid coordinate of cell faces in x-direction.
};

/**
 * Class for the soil grid settings and operators.
 * This class contains the grid properties, such as dimensions and resolution.
 */
template<typename TF>
class Soil_grid
{
    public:
        Soil_grid(Master&, Grid<TF>&, Input&);  // Constructor of the grid class.
        ~Soil_grid();  // Destructor of the grid class.

        void init();                  // Initialization of the grid arrays.
        void create(Netcdf_handle&);  // Creation of the grid data.

        const Soil_grid_data<TF>& get_grid_data(); // Function to get grid data struct

    private:
        Master& master;    // Reference to master class.
        Grid<TF>& grid;    // Reference to atmospheric grid

        bool sw_land_surface;
        Soil_grid_data<TF> gd;  // Struct holding the grid data
};
#endif
