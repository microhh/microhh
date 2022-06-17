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

#include <iostream>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "column.h"
#include "tools.h"
#include "netcdf_interface.h"

#ifdef USECUDA
template<typename TF>
void Column<TF>::calc_column(
        std::string profname, const TF* const restrict data, const TF offset, const bool copy_from_gpu)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    for (auto& col : columns)
    {
        if ( (col.coord[0] / gd.imax == md.mpicoordx ) && (col.coord[1] / gd.jmax == md.mpicoordy ) )
        {
            const int i_col = col.coord[0] % gd.imax + gd.istart;
            const int j_col = col.coord[1] % gd.jmax + gd.jstart;

            if (copy_from_gpu)
            {
                const int kbeg  = i_col + j_col * gd.icells;

                cuda_safe_call(cudaMemcpy2D(
                            col.profs.at(profname).data.data(),
                            sizeof(TF), &data[kbeg], gd.ijcells*sizeof(TF), sizeof(TF), gd.kcells, cudaMemcpyDeviceToHost));

                for (int k=0; k<gd.kcells; ++k)
                    col.profs.at(profname).data[k] += offset;
            }
            else
            {
                for (int k=0; k<gd.kcells; k++)
                {
                    const int ijk = i_col + j_col*gd.icells + k*gd.ijcells;
                    col.profs.at(profname).data[k] = (data[ijk] + offset);
                }
            }
        }
    }
}

template<typename TF>
void Column<TF>::calc_time_series(
        std::string name, const TF* const restrict data, const TF offset)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    for (auto& col : columns)
    {
        // Check if coordinate is in range.
        if ( (col.coord[0] / gd.imax == md.mpicoordx ) && (col.coord[1] / gd.jmax == md.mpicoordy ) )
        {
            const int i_col = col.coord[0] % gd.imax + gd.istart;
            const int j_col = col.coord[1] % gd.jmax + gd.jstart;
            const int ij = i_col + j_col*gd.icells;

            cuda_safe_call(cudaMemcpy(
                    &col.time_series.at(name).data, &data[ij], sizeof(TF), cudaMemcpyDeviceToHost));

            col.time_series.at(name).data += offset;
        }
    }
}

template<typename TF>
int* Column<TF>::get_column_location_g(const std::string& dim)
{
    if (dim == "i")
        return col_i_g;
    else if (dim == "j")
        return col_j_g;
    else
        throw std::runtime_error("Can not get column locations for dimension " + dim);
}

template<typename TF>
void Column<TF>::prepare_device()
{
    // Init individual column statistics.
    std::vector<int> col_i;
    std::vector<int> col_j;
    get_column_locations(col_i, col_j);

    if (col_i.size() > 0)
    {
        const int colsize = col_i.size() * sizeof(int);

        cuda_safe_call(cudaMalloc(&col_i_g, colsize));
        cuda_safe_call(cudaMalloc(&col_j_g, colsize));

        cuda_safe_call(cudaMemcpy(col_i_g, col_i.data(), colsize, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(col_j_g, col_j.data(), colsize, cudaMemcpyHostToDevice));
    }
}

template<typename TF>
void Column<TF>::clear_device()
{
    int n = get_n_columns();

    if (n > 0)
    {
        cuda_safe_call(cudaFree(col_i_g));
        cuda_safe_call(cudaFree(col_j_g));
    }
}
#endif

template class Column<double>;
template class Column<float>;
