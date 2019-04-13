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

    // CvH : This does not work in case CUDA is combined with MPI.
    for (auto& col : columns)
    {
        const int kbeg = col.coord[0] + col.coord[1]*gd.icells;

        cuda_safe_call(cudaMemcpy2D(
                    col.profs.at(profname).data.data(),
                    sizeof(TF), &data[kbeg], gd.ijcells*sizeof(TF), sizeof(TF), gd.kcells, cudaMemcpyDeviceToHost));

        for (int k=0; k<gd.kcells; ++k)
            col.profs.at(profname).data[k] += offset;
    }
}
#endif

template class Column<double>;
template class Column<float>;
