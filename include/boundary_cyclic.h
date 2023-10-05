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

#ifndef BOUNDARY_CYCLIC_H
#define BOUNDARY_CYCLIC_H

#ifdef USEMPI
#include <mpi.h>
#endif

class Master;
template<typename> class Grid;

enum class Edge {East_west_edge, North_south_edge, Both_edges};

template<typename TF>
class Boundary_cyclic
{
    public:
        Boundary_cyclic(Master&, Grid<TF>&); // Constuctor of the boundary class.
        ~Boundary_cyclic();                  // Destructor of the boundary class.

        void init();   // Initialize the fields.
        void exec(TF* const restrict, Edge=Edge::Both_edges); // Fills the ghost cells in the periodic directions.
        void exec_2d(TF* const restrict); // Fills the ghost cells of one slice in the periodic direction.

        void exec(unsigned int* const restrict, Edge=Edge::Both_edges); // Fills the ghost cells in the periodic directions.
        void exec_2d(unsigned int* const restrict); // Fills the ghost cells of one slice in the periodic direction.

        void exec_g(TF*);   // Fills the ghost cells in the periodic directions.
        void exec_2d_g(TF*); // Fills the ghost cells of one slice in the periodic directions.

    private:
        Master& master; // Reference to master class.
        Grid<TF>& grid; // Reference to grid class.

        void init_mpi();
        void exit_mpi();
        bool mpi_types_allocated;

        #ifdef USEMPI
        MPI_Datatype eastwestedge;     ///< MPI datatype containing the ghostcells at the east-west sides.
        MPI_Datatype northsouthedge;   ///< MPI datatype containing the ghostcells at the north-south sides.
        MPI_Datatype eastwestedge2d;   ///< MPI datatype containing the ghostcells for one slice at the east-west sides.
        MPI_Datatype northsouthedge2d; ///< MPI datatype containing the ghostcells for one slice at the north-south sides.
        MPI_Datatype eastwestedge_uint;     ///< MPI datatype containing the ghostcells at the east-west sides.
        MPI_Datatype northsouthedge_uint;   ///< MPI datatype containing the ghostcells at the north-south sides.
        MPI_Datatype eastwestedge2d_uint;   ///< MPI datatype containing the ghostcells for one slice at the east-west sides.
        MPI_Datatype northsouthedge2d_uint; ///< MPI datatype containing the ghostcells for one slice at the north-south sides.
        #endif
};
#endif
