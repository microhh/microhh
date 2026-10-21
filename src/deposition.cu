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

#include "tools.h"
#include "constants.h"
#include "land_surface_kernels_gpu.h"

#include "deposition.h"
#include "deposition_kernels.cuh"

namespace dkg  = Deposition_kernels_g;
namespace lskg = Land_surface_kernels_g;


#ifdef USECUDA
template <typename TF>
void Deposition<TF>::update_time_dependent(
        Timeloop<TF>& timeloop,
        Boundary<TF>& boundary,
        TF* restrict vdo3_g,
        TF* restrict vdno_g,
        TF* restrict vdno2_g,
        TF* restrict vdhno3_g,
        TF* restrict vdh2o2_g,
        TF* restrict vdrooh_g,
        TF* restrict vdhcho_g
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

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);
    dim3 grid_gpu (gridi,  gridj,  1);
    dim3 block_gpu(blocki, blockj, 1);

    auto& dep_veg  = deposition_tiles.at("veg");
    auto& dep_soil = deposition_tiles.at("soil");
    auto& dep_wet  = deposition_tiles.at("wet");

    dkg::calc_deposition_veg<TF><<<grid_gpu, block_gpu>>>(
            dep_veg.vd_g.at("o3"),
            dep_veg.vd_g.at("no"),
            dep_veg.vd_g.at("no2"),
            dep_veg.vd_g.at("hno3"),
            dep_veg.vd_g.at("h2o2"),
            dep_veg.vd_g.at("rooh"),
            dep_veg.vd_g.at("hcho"),
            lai_g,
            tiles.at("veg").rs_g,
            tiles.at("veg").ra_g,
            tiles.at("veg").ustar_g,
            tiles.at("veg").fraction_g,
            rmes_g,
            rsoil_g,
            rcut_g,
            diff_scl_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
    cuda_check_error();

    dkg::calc_deposition_soil<TF><<<grid_gpu, block_gpu>>>(
            dep_soil.vd_g.at("o3"),
            dep_soil.vd_g.at("no"),
            dep_soil.vd_g.at("no2"),
            dep_soil.vd_g.at("hno3"),
            dep_soil.vd_g.at("h2o2"),
            dep_soil.vd_g.at("rooh"),
            dep_soil.vd_g.at("hcho"),
            tiles.at("soil").ra_g,
            tiles.at("soil").ustar_g,
            tiles.at("soil").fraction_g,
            rsoil_g,
            diff_scl_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
    cuda_check_error();

    dkg::calc_deposition_wet<TF><<<grid_gpu, block_gpu>>>(
            dep_wet.vd_g.at("o3"),
            dep_wet.vd_g.at("no"),
            dep_wet.vd_g.at("no2"),
            dep_wet.vd_g.at("hno3"),
            dep_wet.vd_g.at("h2o2"),
            dep_wet.vd_g.at("rooh"),
            dep_wet.vd_g.at("hcho"),
            lai_g,
            c_veg_g,
            tiles.at("wet").rs_g,
            tiles.at("veg").rs_g,
            tiles.at("wet").ra_g,
            tiles.at("wet").ustar_g,
            tiles.at("wet").fraction_g,
            rmes_g,
            rsoil_g,
            rws_g,
            diff_scl_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
    cuda_check_error();

    auto calc_vd_g = [&](TF* vd, const std::string& sp)
    {
        lskg::calc_tiled_mean_g<TF><<<grid_gpu, block_gpu>>>(
                vd,
                dep_veg.vd_g.at(sp),
                dep_soil.vd_g.at(sp),
                dep_wet.vd_g.at(sp),
                tiles.at("veg").fraction_g,
                tiles.at("soil").fraction_g,
                tiles.at("wet").fraction_g,
                TF(1),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells);
        cuda_check_error();

        // Use wet-tile u* and ra: calculated in lsm with f_wet = 100%.
        const int s = species_idx.at(sp);
        dkg::calc_vd_water<TF><<<grid_gpu, block_gpu>>>(
                vd,
                tiles.at("wet").ra_g,
                tiles.at("wet").ustar_g,
                water_mask_g,
                diff_scl[s], rwat[s],
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells);
        cuda_check_error();

        // TODO: spatial_avg_vd.
    };

    calc_vd_g(vdo3_g,  "o3");
    calc_vd_g(vdno_g,  "no");
    calc_vd_g(vdno2_g, "no2");
    calc_vd_g(vdhno3_g,"hno3");
    calc_vd_g(vdh2o2_g,"h2o2");
    calc_vd_g(vdrooh_g,"rooh");
    calc_vd_g(vdhcho_g,"hcho");
}

template <typename TF>
void Deposition<TF>::prepare_device()
{
    if (!sw_deposition)
        return;

    auto& gd = grid.get_grid_data();

    // Allocate GPU arrays.
    for (auto& tile : deposition_tiles)
        for (auto& [sp, idx] : species_idx)
            tile.second.vd_g[sp].allocate(gd.ijcells);

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
        for (auto& [sp, idx] : species_idx)
            cuda_safe_call(cudaMemcpy(tile.second.vd_g.at(sp), tile.second.vd.at(sp).data(), memsize_ij, cudaMemcpyHostToDevice));

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

template <typename TF>
void Deposition<TF>::backward_device()
{
    if (!sw_deposition)
        return;

    auto& gd = grid.get_grid_data();

    for (auto& tile : deposition_tiles)
        for (auto& [sp, idx] : species_idx)
            cuda_safe_call(cudaMemcpy(tile.second.vd.at(sp).data(), tile.second.vd_g.at(sp), gd.ijcells * sizeof(TF), cudaMemcpyDeviceToHost));
}
#endif

#ifdef FLOAT_SINGLE
template class Deposition<float>;
#else
template class Deposition<double>;
#endif