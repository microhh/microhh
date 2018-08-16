/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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


#ifdef USECUDA
template<typename TF>
void Column<TF>::calc_column(std::string profname, const TF* const restrict data,
                      const TF offset)
{
    auto& gd = grid.get_grid_data();

    for(auto& it: columns)
    {
        const int kbeg = it.coord[0] + it.coord[1] * gd.icells + gd.ijcells * gd.kstart;

        cuda_safe_call(cudaMemcpy2D(&it.profs.at(profname).data.data()[gd.kstart], sizeof(TF),&data[kbeg], gd.ijcells*sizeof(TF), sizeof(TF), gd.kmax, cudaMemcpyDeviceToHost));

        for (int k=gd.kstart; k<gd.kend; k++)
            it.profs.at(profname).data.data()[k] += offset;
    }
}
#endif

template class Column<double>;
template class Column<float>;
