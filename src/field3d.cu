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

#include <map>
#include <vector>
#include "field3d.h"
#include "grid.h"
#include "master.h"
#include "tools.h"
#include <iostream>

template<typename TF>
void Field3d<TF>::init_device()
{
    const Grid_data<TF>& gd = grid.get_grid_data();
    fld_g.resize(gd.ncells);
    fld_bot_g.resize(gd.ijcells);
    fld_top_g.resize(gd.ijcells);
    grad_bot_g.resize(gd.ijcells);
    grad_top_g.resize(gd.ijcells);
    flux_bot_g.resize(gd.ijcells);
    flux_top_g.resize(gd.ijcells);
    fld_mean_g.resize(gd.kcells);
}

template<typename TF>
void Field3d<TF>::clear_device()
{
    fld_g.free();
    fld_bot_g.free();
    fld_top_g.free();
    grad_bot_g.free();
    grad_top_g.free();
    flux_bot_g.free();
    flux_top_g.free();
    fld_mean_g.free();
}


#ifdef FLOAT_SINGLE
template class Field3d<float>;
#else
template class Field3d<double>;
#endif
