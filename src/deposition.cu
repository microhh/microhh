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

#include "deposition.h"
#include "tools.h"

namespace
{
}

#ifdef USECUDA
template <typename TF>
void Deposition<TF>::update_time_dependent(
        Timeloop<TF>& timeloop,
        Boundary<TF>& boundary,
        TF* restrict vdo3,
        TF* restrict vdno,
        TF* restrict vdno2,
        TF* restrict vdhno3,
        TF* restrict vdh2o2,
        TF* restrict vdrooh,
        TF* restrict vdhcho
        )
{
    if (!sw_deposition)
        return;

    auto& gd = grid.get_grid_data();

    // Get information from the land-surface model:
    auto& tiles = boundary.get_tiles();
    int* water_mask_g = boundary.get_water_mask_g();
    TF* lai_g = boundary.get_lai_g();
    TF* c_veg_g = boundary.get_c_veg_g();

    // TODO
}

template <typename TF>
void Deposition<TF>::prepare_device()
{
    if (!sw_deposition)
        return;

    auto& gd = grid.get_grid_data();

    // Allocate GPU arrays.
    for (auto& tile : deposition_tiles)
    {
        tile.second.vdo3_g.allocate(gd.ijcells);
        tile.second.vdno_g.allocate(gd.ijcells);
        tile.second.vdno2_g.allocate(gd.ijcells);
        tile.second.vdhno3_g.allocate(gd.ijcells);
        tile.second.vdh2o2_g.allocate(gd.ijcells);
        tile.second.vdrooh_g.allocate(gd.ijcells);
        tile.second.vdhcho_g.allocate(gd.ijcells);
    }

    const int size = rmes.size();
    rmes_g.allocate(size);
    rsoil_g.allocate(size);
    rcut_g.allocate(size);
    rws_g.allocate(size);
    rwat_g.allocate(size);
    diff_g.allocate(size);
    diff_scl_g.allocate(size);
    henry_g.allocate(size);
    f0_g.allocate(size);

    // Copy data to device.
    const int memsize_ij  = gd.ijcells * sizeof(TF);
    const int memsize_res = size * sizeof(TF);

    for (auto& tile : deposition_tiles)
    {
        cuda_safe_call(cudaMemcpy(tile.second.vdo3_g,  tile.second.vdo3.data(),  memsize_ij, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.second.vdno_g,  tile.second.vdno.data(),  memsize_ij, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.second.vdno2_g, tile.second.vdno2.data(), memsize_ij, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.second.vdhno3_g,tile.second.vdhno3.data(),memsize_ij, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.second.vdh2o2_g,tile.second.vdh2o2.data(),memsize_ij, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.second.vdrooh_g,tile.second.vdrooh.data(),memsize_ij, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.second.vdhcho_g,tile.second.vdhcho.data(),memsize_ij, cudaMemcpyHostToDevice));
    }

    cuda_safe_call(cudaMemcpy(rmes_g,    rmes.data(),    memsize_res, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(rsoil_g,   rsoil.data(),   memsize_res, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(rcut_g,    rcut.data(),    memsize_res, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(rws_g,     rws.data(),     memsize_res, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(rwat_g,    rwat.data(),    memsize_res, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(diff_g,    diff.data(),    memsize_res, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(diff_scl_g,diff_scl.data(),memsize_res, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(henry_g,   henry.data(),   memsize_res, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(f0_g,      f0.data(),      memsize_res, cudaMemcpyHostToDevice));
}
#endif

#ifdef FLOAT_SINGLE
template class Deposition<float>;
#else
template class Deposition<double>;
#endif