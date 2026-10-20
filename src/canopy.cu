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

#include <algorithm>

#include "canopy.h"
#include "grid.h"
#include "fields.h"
#include "field3d_operators.h"
#include "field3d_io.h"
#include "tools.h"

#include "canopy_kernels.cuh"
#include "cuda_launcher.h"
#include "cuda_tiling.h"

#ifdef USECUDA
namespace
{
    template<typename TF, bool sw_pad_3d>
    void canopy_drag(
            Fields<TF>& fields,
            const cuda_vector<TF>& pad_g,
            const Grid_layout& grid_layout_uv,
            const Grid_layout& grid_layout_w,
            const TF utrans,
            const TF vtrans,
            const TF cd)
    {
        launch_grid_kernel<Canopy_kernels::canopy_drag_u_g<TF, sw_pad_3d>>(
                grid_layout_uv,
                fields.mt.at("u")->fld_g.view(),
                fields.mp.at("u")->fld_g,
                fields.mp.at("v")->fld_g,
                fields.mp.at("w")->fld_g,
                pad_g,
                utrans,
                vtrans,
                cd);

        launch_grid_kernel<Canopy_kernels::canopy_drag_v_g<TF, sw_pad_3d>>(
                grid_layout_uv,
                fields.mt.at("v")->fld_g.view(),
                fields.mp.at("u")->fld_g,
                fields.mp.at("v")->fld_g,
                fields.mp.at("w")->fld_g,
                pad_g,
                utrans,
                vtrans,
                cd);

        launch_grid_kernel<Canopy_kernels::canopy_drag_w_g<TF, sw_pad_3d>>(
                grid_layout_w,
                fields.mt.at("w")->fld_g.view(),
                fields.mp.at("u")->fld_g,
                fields.mp.at("v")->fld_g,
                fields.mp.at("w")->fld_g,
                pad_g,
                utrans,
                vtrans,
                cd);
    }
}

template <typename TF>
void Canopy<TF>::exec()
{
    if (!sw_canopy)
        return;

    auto& gd = grid.get_grid_data();

    // No canopy present, prevent launching kernels with an empty grid.
    if (kend_canopy == gd.kstart)
        return;

    Grid_layout grid_layout_uv = {
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, kend_canopy,
            gd.istride,
            gd.jstride,
            gd.kstride};

    // Include the half level at the top of the canopy, where `pad` interpolated to `w` is non-zero.
    const int kend_canopy_w = std::min(kend_canopy+1, gd.kend);

    Grid_layout grid_layout_w = {
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart+1, kend_canopy_w,
            gd.istride,
            gd.jstride,
            gd.kstride};

    if (sw_3d_pad)
        canopy_drag<TF, true>(fields, pad_g, grid_layout_uv, grid_layout_w, gd.utrans, gd.vtrans, cd);
    else
        canopy_drag<TF, false>(fields, pad_g, grid_layout_uv, grid_layout_w, gd.utrans, gd.vtrans, cd);
}

template <typename TF>
void Canopy<TF>::prepare_device()
{
    if (!sw_canopy)
        return;

    // `pad` is either the 1D profile or the 3D field, including the ghost cells filled in `create()`.
    pad_g.allocate(pad.size());
    cuda_safe_call(cudaMemcpy(pad_g, pad.data(), pad.size()*sizeof(TF), cudaMemcpyHostToDevice));
}

template <typename TF>
void Canopy<TF>::clear_device()
{
    // No need for cleanup with `cuda_vector`'s.
}
#endif

#ifdef FLOAT_SINGLE
template class Canopy<float>;
#else
template class Canopy<double>;
#endif
