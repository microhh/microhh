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

#ifndef TRANSPOSE
#define TRANSPOSE

#ifdef USEMPI
#include <mpi.h>
#endif

class Master;
template<typename> class Grid;

template<typename TF>
class Transpose
{
    public:
        Transpose(Master&, Grid<TF>&);
        ~Transpose();

        void init();

    private:
        Master& master;
        Grid<TF>& grid;

        void init_mpi();
        void exit_mpi();
        bool mpi_types_allocated;

        #ifdef USEMPI
        MPI_Datatype transposez;  ///< MPI datatype containing base blocks for z-orientation in zx-transpose.
        MPI_Datatype transposez2; ///< MPI datatype containing base blocks for z-orientation in zy-transpose.
        MPI_Datatype transposex;  ///< MPI datatype containing base blocks for x-orientation in zx-transpose.
        MPI_Datatype transposex2; ///< MPI datatype containing base blocks for x-orientation in xy-transpose.
        MPI_Datatype transposey;  ///< MPI datatype containing base blocks for y-orientation in xy-transpose.
        MPI_Datatype transposey2; ///< MPI datatype containing base blocks for y-orientation in zy-transpose.
        #endif
};
#endif
