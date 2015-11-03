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

struct IB_cell
{
    int i;
    int j;
    int k;
    std::bitset<6> sw_outer_wall;
};

enum IB_edge {Bottom_edge, Top_edge, West_edge, East_edge, South_edge, North_edge};

class Immersed_boundary
{
    public:
        Immersed_boundary(Model*, Input*); ///< Constructor of the class.
        ~Immersed_boundary();              ///< Destructor of the class.

        void init();
        void create();
        void exec();

        void exec_stats(Mask*); ///< Execute statistics of immersed boundaries

    private:
        Model*  model;  ///< Pointer to model class.
        Fields* fields; ///< Pointer to fields class.
        Grid*   grid;   ///< Pointer to grid class.
        Stats*  stats;  ///< Pointer to grid class.

        std::vector<IB_cell> IB_cells;

        // Some (temporary?) statistics walls
        double boundary_u_max;
        double boundary_w_max;
};
#endif
