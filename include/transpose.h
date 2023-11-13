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

#ifndef TRANSPOSE_H
#define TRANSPOSE_H

#ifdef USEMPI
#include <mpi.h>
#endif

#include "defines.h"

class Master;
template<typename> class Grid;

template<typename TF>
class Transpose
{
    public:
        Transpose(Master&, Grid<TF>&);
        ~Transpose();

        void init();

        void exec_zx(TF* const restrict, TF* const restrict); ///< Changes the transpose orientation from z to x.
        void exec_xz(TF* const restrict, TF* const restrict); ///< Changes the transpose orientation from x to z.
        void exec_xy(TF* const restrict, TF* const restrict); ///< changes the transpose orientation from x to y.
        void exec_yx(TF* const restrict, TF* const restrict); ///< Changes the transpose orientation from y to x.
        void exec_yz(TF* const restrict, TF* const restrict); ///< Changes the transpose orientation from y to z.
        void exec_zy(TF* const restrict, TF* const restrict); ///< Changes the transpose orientation from z to y.

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
