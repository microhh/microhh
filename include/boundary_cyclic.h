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

#ifndef BOUNDARY_CYCLIC
#define BOUNDARY_CYCLIC

class Master;
template<typename> class Grid;

template<typename TF>
class Boundary_cyclic
{
    public:
        Boundary_cyclic(Master&, Grid<TF>&); // Constuctor of the boundary class.
        ~Boundary_cyclic();                  // Destructor of the boundary class.

        void init();   // Initialize the fields.
        void create(); // Create the fields.

    private:
        Master& master; // Reference to master class.
        Grid<TF>& grid; // Reference to grid class.

        #ifdef USEMPI
        MPI_Datatype eastwestedge;     ///< MPI datatype containing the ghostcells at the east-west sides.
        MPI_Datatype northsouthedge;   ///< MPI datatype containing the ghostcells at the north-south sides.
        #endif
};
#endif
