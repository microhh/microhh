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

#include "master.h"
#include "grid.h"
#include "transpose.h"

template<typename TF>
Transpose<TF>::Transpose(Master& masterin, Grid<TF>& gridin) :
    master(masterin),
    grid(gridin),
    mpi_types_allocated(false)
{
}

template<typename TF>
Transpose<TF>::~Transpose()
{
    exit_mpi();
}

template<typename TF>
void Transpose<TF>::init()
{
    init_mpi();
}

#ifdef USEMPI
namespace
{
    template<typename TF> MPI_Datatype mpi_fp_type();
    template<> MPI_Datatype mpi_fp_type<double>() { return MPI_DOUBLE; }
    template<> MPI_Datatype mpi_fp_type<float>() { return MPI_FLOAT; }
}

template<typename TF>
void Transpose<TF>::init_mpi()
{
    auto& gd = grid.get_grid_data();

    int datacount, datablock, datastride;

    mpi_types_allocated = true;
}

template<typename TF>
void Transpose<TF>::exit_mpi()
{
    if (mpi_types_allocated)
    {
    }
}


#else

template<typename TF>
void Transpose<TF>::init_mpi()
{
}

template<typename TF>
void Transpose<TF>::exit_mpi()
{
}
#endif

template class Transpose<double>;
template class Transpose<float>;
