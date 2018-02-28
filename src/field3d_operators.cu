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

#include <cstdio>
#include <iostream>
#include <cmath>
#include "master.h"
#include "grid.h"
#include "field3d.h"
#include "fields.h"
#include "defines.h"
#include "tools.h"
#include "field3d_operators.h"

#ifdef USECUDA
template<typename TF>
void Field3d_operators<TF>::calc_mean_profile(TF* const restrict prof, const TF* const restrict fld)
{
    using namespace Tools_g;

    const Grid_data<TF>& gd = grid.get_grid_data();
    const TF scalefac = 1./(gd.itot*gd.jtot);
    auto tmp = fields.get_tmp_g();
    // Reduce 3D field excluding ghost cells and padding to jtot*kcells values
    reduce_interior<TF>(&fld[gd.memoffset], &tmp->fld_g[gd.memoffset], gd.itot, gd.istart, gd.iend, gd.jtot, gd.jstart, gd.jend, gd.kcells, 0, gd.icellsp, gd.ijcellsp, Sum_type);
    // Reduce jtot*kcells to kcells values
    reduce_all<TF>     (&tmp->fld_g[gd.memoffset], prof, gd.jtot*gd.kcells, gd.kcells, gd.jtot, Sum_type, scalefac);
    fields.release_tmp_g(tmp);
}

template<typename TF>
TF Field3d_operators<TF>::calc_mean(const TF* const restrict fld)
{
    using namespace Tools_g;

    const Grid_data<TF>& gd = grid.get_grid_data();
    const TF scalefac = 1./(gd.itot*gd.jtot*gd.ktot);
    TF mean_value;

    auto tmp = fields.get_tmp_g();
    // Reduce 3D field excluding ghost cells and padding to jtot*ktot values
    reduce_interior<TF>(&fld[gd.memoffset], &tmp->fld_g[gd.memoffset], gd.itot, gd.istart, gd.iend, gd.jtot, gd.jstart, gd.jend, gd.ktot, gd.kstart, gd.icellsp, gd.ijcellsp, Sum_type);
    // Reduce jtot*ktot to ktot values
    reduce_all<TF>     (&tmp->fld_g[gd.memoffset], &tmp->fld_g[gd.memoffset+gd.jtot*gd.ktot], gd.jtot*gd.ktot, gd.ktot, gd.jtot, Sum_type, 1.);
    // Reduce ktot values to a single value
    reduce_all<TF>     (&tmp->fld_g[gd.memoffset+gd.jtot*gd.ktot], &tmp->fld_g[gd.memoffset], gd.ktot, 1, gd.ktot, Sum_type, scalefac);
    // Copy back result from GPU
    cuda_safe_call(cudaMemcpy(&mean_value, &tmp->fld_g[gd.memoffset], sizeof(TF), cudaMemcpyDeviceToHost));
    fields.release_tmp_g(tmp);
    return mean_value;
}

template<typename TF>
TF Field3d_operators<TF>::calc_max(const TF* const restrict fld)
{
    using namespace Tools_g;

    const Grid_data<TF>& gd = grid.get_grid_data();
    const TF scalefac = 1.;
    TF max_value;

    auto tmp = fields.get_tmp_g();
    // Reduce 3D field excluding ghost cells and padding to jtot*ktot values
    reduce_interior<TF>(&fld[gd.memoffset], &tmp->fld_g[gd.memoffset], gd.itot, gd.istart, gd.iend, gd.jtot, gd.jstart, gd.jend, gd.ktot, gd.kstart, gd.icellsp, gd.ijcellsp, Max_type);
    // Reduce jtot*ktot to ktot values
    reduce_all<TF>     (&tmp->fld_g[gd.memoffset], &tmp->fld_g[gd.memoffset+gd.jtot*gd.ktot], gd.jtot*gd.ktot, gd.ktot, gd.jtot, Max_type, scalefac);
    // Reduce ktot values to a single value
    reduce_all<TF>     (&tmp->fld_g[gd.memoffset+gd.jtot*gd.ktot], &tmp->fld_g[gd.memoffset], gd.ktot, 1, gd.ktot, Max_type, scalefac);
    // Copy back result from GPU
    cuda_safe_call(cudaMemcpy(&max_value, &tmp->fld_g[gd.memoffset], sizeof(TF), cudaMemcpyDeviceToHost));
    fields.release_tmp_g(tmp);

    return max_value;
}
#endif

template class Field3d_operators<double>;
template class Field3d_operators<float>;
