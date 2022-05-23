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

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "column.h"
#include "tools.h"
#include "netcdf_interface.h"

#ifdef USECUDA
template<typename TF>
void Column<TF>::calc_column(
        std::string profname, const TF* const restrict data, const TF offset)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    for (auto& col : columns)
    {
        if ( (col.coord[0] / gd.imax == md.mpicoordx ) && (col.coord[1] / gd.jmax == md.mpicoordy ) )
        {
            const int i_col = col.coord[0] % gd.imax + gd.istart;
            const int j_col = col.coord[1] % gd.jmax + gd.jstart;
            const int kbeg  = i_col + j_col * gd.icells;

            cuda_safe_call(cudaMemcpy2D(
                        col.profs.at(profname).data.data(),
                        sizeof(TF), &data[kbeg], gd.ijcells*sizeof(TF), sizeof(TF), gd.kcells, cudaMemcpyDeviceToHost));

            for (int k=0; k<gd.kcells; ++k)
                col.profs.at(profname).data[k] += offset;
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
void Column<TF>::set_individual_column(
        std::string profname, const TF* const restrict prof,
        const TF offset, const int i_col, const int j_col)
{
}
#endif

template class Column<double>;
template class Column<float>;
