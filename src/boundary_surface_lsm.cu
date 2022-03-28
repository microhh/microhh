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

#include "boundary_surface_lsm.h"
#include "boundary.h"
#include "land_surface_kernels_gpu.h"
#include "tools.h"
#include "grid.h"
#include "soil_grid.h"

namespace
{
    namespace lsmk = Land_surface_kernels_g;
}


#ifdef USECUDA
template<typename TF>
void Boundary_surface_lsm<TF>::exec(
        Thermo<TF>& thermo, Radiation<TF>& radiation,
        Microphys<TF>& microphys, Timeloop<TF>& timeloop)
{
}

template<typename TF>
void Boundary_surface_lsm<TF>::exec_column(Column<TF>& column)
{
}

template<typename TF>
void Boundary_surface_lsm<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    const int tf_memsize_ij  = gd.ijcells*sizeof(TF);
    const int int_memsize_ij = gd.ijcells*sizeof(int);
    const int float_memsize_mo_lut = nzL_lut*sizeof(float);

    // Surface layer / Monin-Obukhov:
    cuda_safe_call(cudaMalloc(&obuk_g,  tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&ustar_g, tf_memsize_ij));

    cuda_safe_call(cudaMalloc(&z0m_g,   tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&z0h_g,   tf_memsize_ij));

    cuda_safe_call(cudaMalloc(&dudz_mo_g, tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&dvdz_mo_g, tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&dbdz_mo_g, tf_memsize_ij));

    cuda_safe_call(cudaMalloc(&nobuk_g, int_memsize_ij));
    cuda_safe_call(cudaMalloc(&zL_sl_g, float_memsize_mo_lut));
    cuda_safe_call(cudaMalloc(&f_sl_g,  float_memsize_mo_lut));

    // Land-surface:
    // 1. Init tiles:
    for (auto& tile : tiles)
        lsmk::init_tile(tile.second, gd.ijcells);

    // 2. Init 2D surface properties:
    cuda_safe_call(cudaMalloc(&gD_coeff_g, tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&c_veg_g, tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&lai_g, tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&rs_veg_min_g, tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&rs_soil_min_g, tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&lambda_stable_g, tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&lambda_unstable_g, tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&cs_veg_g, tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&water_mask_g, int_memsize_ij));
    cuda_safe_call(cudaMalloc(&t_bot_water_g, tf_memsize_ij));

    cuda_safe_call(cudaMalloc(&interception_g, tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&throughfall_g, tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&infiltration_g, tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&runoff_g, tf_memsize_ij));

    // 3. Init 3D soil properties:
    const int tf_memsize_ijk  = sgd.ncells*sizeof(TF);
    const int int_memsize_ijk = sgd.ncells*sizeof(int);

    cuda_safe_call(cudaMalloc(&soil_index_g, int_memsize_ijk));
    cuda_safe_call(cudaMalloc(&diffusivity_g, tf_memsize_ijk));
    cuda_safe_call(cudaMalloc(&diffusivity_h_g, tf_memsize_ijk));
    cuda_safe_call(cudaMalloc(&conductivity_g, tf_memsize_ijk));
    cuda_safe_call(cudaMalloc(&conductivity_h_g, tf_memsize_ijk));
    cuda_safe_call(cudaMalloc(&source_g, tf_memsize_ijk));
    cuda_safe_call(cudaMalloc(&root_fraction_g, tf_memsize_ijk));

    // 4. Init lookup table with van Genuchten parameters:
    const int memsize_vg_lut = theta_res.size() * sizeof(TF);

    cuda_safe_call(cudaMalloc(&theta_res_g, memsize_vg_lut));
    cuda_safe_call(cudaMalloc(&theta_wp_g, memsize_vg_lut));
    cuda_safe_call(cudaMalloc(&theta_fc_g, memsize_vg_lut));
    cuda_safe_call(cudaMalloc(&theta_sat_g, memsize_vg_lut));

    cuda_safe_call(cudaMalloc(&gamma_theta_sat_g, memsize_vg_lut));

    cuda_safe_call(cudaMalloc(&vg_a_g, memsize_vg_lut));
    cuda_safe_call(cudaMalloc(&vg_l_g, memsize_vg_lut));
    cuda_safe_call(cudaMalloc(&vg_n_g, memsize_vg_lut));
    cuda_safe_call(cudaMalloc(&vg_m_g, memsize_vg_lut));

    cuda_safe_call(cudaMalloc(&kappa_theta_max_g, memsize_vg_lut));
    cuda_safe_call(cudaMalloc(&kappa_theta_min_g, memsize_vg_lut));
    cuda_safe_call(cudaMalloc(&gamma_theta_max_g, memsize_vg_lut));
    cuda_safe_call(cudaMalloc(&gamma_theta_min_g, memsize_vg_lut));
    cuda_safe_call(cudaMalloc(&gamma_T_dry_g, memsize_vg_lut));
    cuda_safe_call(cudaMalloc(&rho_C_g, memsize_vg_lut));

    // Copy data from host to device
    forward_device();
}

template<typename TF>
void Boundary_surface_lsm<TF>::forward_device()
{
    auto& gd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    const int tf_memsize_ij  = gd.ijcells*sizeof(TF);
    const int int_memsize_ij = gd.ijcells*sizeof(int);
    const int float_memsize_lut = nzL_lut*sizeof(float);

    // Surface layer / Monin-Obukhov:
    cuda_safe_call(cudaMemcpy(obuk_g,  obuk.data(),  tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(ustar_g, ustar.data(), tf_memsize_ij, cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy(z0m_g, z0m.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(z0h_g, z0h.data(), tf_memsize_ij, cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy(dudz_mo_g, dudz_mo.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(dvdz_mo_g, dvdz_mo.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(dbdz_mo_g, dbdz_mo.data(), tf_memsize_ij, cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy(nobuk_g, nobuk.data(), int_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(zL_sl_g, zL_sl.data(), float_memsize_lut, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(f_sl_g,  f_sl.data(),  float_memsize_lut, cudaMemcpyHostToDevice));

    // Land-surface:
    // 1. Copy tiles:
    for (auto& tile : tiles)
        lsmk::forward_device_tile(tile.second, gd.ijcells);

    // 2. Copy 2D surface properties:
    cuda_safe_call(cudaMemcpy(gD_coeff_g, gD_coeff.data(), int_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(c_veg_g, c_veg.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(lai_g, lai.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(rs_veg_min_g, rs_veg_min.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(rs_soil_min_g, rs_soil_min.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(lambda_stable_g, lambda_stable.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(lambda_unstable_g, lambda_unstable.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(cs_veg_g, cs_veg.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(water_mask_g, water_mask.data(), int_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(t_bot_water_g, t_bot_water.data(), tf_memsize_ij, cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy(interception_g, interception.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(throughfall_g, throughfall.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(infiltration_g, infiltration.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(runoff_g, runoff.data(), tf_memsize_ij, cudaMemcpyHostToDevice));

    // 3. Copy 3D soil properties:
    const int tf_memsize_ijk  = sgd.ncells*sizeof(TF);
    const int int_memsize_ijk = sgd.ncells*sizeof(int);

    cuda_safe_call(cudaMemcpy(soil_index_g, soil_index.data(), int_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(diffusivity_g, diffusivity.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(diffusivity_h_g, diffusivity_h.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(conductivity_g, conductivity.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(conductivity_h_g, conductivity_h.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(source_g, source.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(root_fraction_g, root_fraction.data(), tf_memsize_ij, cudaMemcpyHostToDevice));

    // 4. Copy lookup table with van Genuchten parameters:
    const int memsize_vg_lut = theta_res.size() * sizeof(TF);

    cuda_safe_call(cudaMemcpy(theta_res_g, theta_res.data(), memsize_vg_lut, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(theta_wp_g, theta_wp.data(), memsize_vg_lut, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(theta_fc_g, theta_fc.data(), memsize_vg_lut, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(theta_sat_g, theta_sat.data(), memsize_vg_lut, cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy(gamma_theta_sat_g, gamma_theta_sat.data(), memsize_vg_lut, cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy(vg_a_g, vg_a.data(), memsize_vg_lut, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vg_l_g, vg_l.data(), memsize_vg_lut, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vg_n_g, vg_n.data(), memsize_vg_lut, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vg_m_g, vg_m.data(), memsize_vg_lut, cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy(kappa_theta_max_g, kappa_theta_max.data(), memsize_vg_lut, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(kappa_theta_min_g, kappa_theta_min.data(), memsize_vg_lut, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gamma_theta_max_g, gamma_theta_max.data(), memsize_vg_lut, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gamma_theta_min_g, gamma_theta_min.data(), memsize_vg_lut, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(gamma_T_dry_g, gamma_T_dry.data(), memsize_vg_lut, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(rho_C_g, rho_C.data(), memsize_vg_lut, cudaMemcpyHostToDevice));
}

template<typename TF>
void Boundary_surface_lsm<TF>::backward_device()
{
    auto& gd = grid.get_grid_data();

    const int tf_memsize_ij  = gd.ijcells*sizeof(TF);
    const int int_memsize_ij = gd.ijcells*sizeof(int);

    // NOTE: only copy back the required/useful data...

    cuda_safe_call(cudaMemcpy(obuk.data(),  obuk_g,  tf_memsize_ij, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(ustar.data(), ustar_g, tf_memsize_ij, cudaMemcpyDeviceToHost));

    cuda_safe_call(cudaMemcpy(dudz_mo.data(), dudz_mo_g, tf_memsize_ij, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(dvdz_mo.data(), dvdz_mo_g, tf_memsize_ij, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(dbdz_mo.data(), dbdz_mo_g, tf_memsize_ij, cudaMemcpyDeviceToHost));

    cuda_safe_call(cudaMemcpy(nobuk.data(), nobuk_g, int_memsize_ij, cudaMemcpyDeviceToHost));

    // TODO: which fields are needed from the land-surface?
}

template<typename TF>
void Boundary_surface_lsm<TF>::clear_device()
{
    //
    // De-llocate fields on GPU
    //
    // Monin-Obukhov stuff:
    cuda_safe_call(cudaFree(obuk_g));
    cuda_safe_call(cudaFree(ustar_g));

    cuda_safe_call(cudaFree(z0m_g));
    cuda_safe_call(cudaFree(z0h_g));

    cuda_safe_call(cudaFree(dudz_mo_g));
    cuda_safe_call(cudaFree(dvdz_mo_g));
    cuda_safe_call(cudaFree(dbdz_mo_g));

    cuda_safe_call(cudaFree(nobuk_g));
    cuda_safe_call(cudaFree(zL_sl_g));
    cuda_safe_call(cudaFree(f_sl_g));

    // Land-surface stuff:
}
#endif

template class Boundary_surface_lsm<double>;
template class Boundary_surface_lsm<float>;
