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

#ifndef IMMERSED_BOUNDARY
#define IMMERSED_BOUNDARY

#include <vector>
#include <bitset>

class Model;
class Grid;
class Fields;
class Input;
class Stats;
struct Mask;

struct Neighbour
{
    int i;
    int j;
    int k;
    double distance;
};

struct Ghost_cell
{
    int i;  ///< i-location of ghost cell
    int j;  ///< j-location of ghost cell
    int k;  ///< k-location of ghost cell

    double xb;  ///< Nearest location at boundary
    double zb;  ///< Nearest location at boundary

    std::vector< std::vector<double> > B;
    std::vector<Neighbour> fluid_neighbours; ///< Neighbouring fluid points used in interpolation
};

struct Cell
{
    int i;  ///< i-location of cell
    int j;  ///< j-location of cell
    int k;  ///< k-location of cell
};

class Immersed_boundary
{
    public:
        Immersed_boundary(Model*, Input*); ///< Constructor of the class.
        ~Immersed_boundary();              ///< Destructor of the class.

        void init();
        void create();
        void exec();
        void exec_tend();

        void exec_stats(Mask*); ///< Execute statistics of immersed boundaries

    private:
        Model*  model;  ///< Pointer to model class.
        Fields* fields; ///< Pointer to fields class.
        Grid*   grid;   ///< Pointer to grid class.
        Stats*  stats;  ///< Pointer to grid class.

        std::vector<Cell> boundary_cells_u;  ///< Vector holding info on all the ghost cells within the boundary
        std::vector<Cell> boundary_cells_w;  ///< Vector holding info on all the ghost cells within the boundary
        std::vector<Cell> boundary_cells_s;  ///< Vector holding info on all the ghost cells within the boundary

        std::vector<Ghost_cell> ghost_cells_u;  ///< Vector holding info on all the ghost cells within the boundary
        std::vector<Ghost_cell> ghost_cells_w;  ///< Vector holding info on all the ghost cells within the boundary
        std::vector<Ghost_cell> ghost_cells_s;  ///< Vector holding info on all the ghost cells within the boundary

        std::string sw_ib; ///< Immersed boundary switch

        double x0_hill;     ///< Middle of Gaussian hill
        double lz_hill;     ///< Height of  "   "    "
        double lx_hill;     ///< Length of  "    "   "
};
#endif
