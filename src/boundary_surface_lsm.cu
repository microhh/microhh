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

#include "boundary.h"
#include "boundary_surface_lsm.h"
#include "land_surface_kernels_gpu.h"
#include "boundary_surface_kernels_gpu.h"
#include "soil_kernels_gpu.h"
#include "tools.h"
#include "grid.h"
#include "soil_grid.h"
#include "fields.h"
#include "soil_field3d.h"
#include "radiation.h"
#include "thermo.h"
#include "microphys.h"

namespace
{
    namespace lsmk = Land_surface_kernels_g;
    namespace bsk = Boundary_surface_kernels_g;
    namespace sk = Soil_kernels_g;
}


#ifdef USECUDA
template<typename TF>
void Boundary_surface_lsm<TF>::exec(
        Thermo<TF>& thermo, Radiation<TF>& radiation,
        Microphys<TF>& microphys, Timeloop<TF>& timeloop)
{
    auto& gd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;

    // For 2D field excluding ghost cells
    int gridi = gd.imax/blocki + (gd.imax%blocki > 0);
    int gridj = gd.jmax/blockj + (gd.jmax%blockj > 0);
    dim3 gridGPU (gridi,  gridj,  1);
    dim3 blockGPU(blocki, blockj, 1);

    // For 2D field including ghost cells
    gridi = gd.icells/blocki + (gd.icells%blocki > 0);
    gridj = gd.jcells/blockj + (gd.jcells%blockj > 0);
    dim3 gridGPU2 (gridi,  gridj,  1);
    dim3 blockGPU2(blocki, blockj, 1);

    // Calculate filtered wind speed difference surface-atmosphere.
    auto tmp1 = fields.get_tmp_g();
    // Aarrghh, TODO: replace with `get_tmp_xy_g()......`.
    TF* du_tot = tmp1->fld_bot_g;

    bsk::calc_dutot_g<<<gridGPU, blockGPU>>>(
        du_tot,
        fields.mp.at("u")->fld_g,
        fields.mp.at("v")->fld_g,
        fields.mp.at("u")->fld_bot_g,
        fields.mp.at("v")->fld_bot_g,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart,
        gd.icells, gd.ijcells);
    cuda_check_error();

    // 2D cyclic boundaries on dutot
    boundary_cyclic.exec_2d_g(du_tot);

    //
    // Retrieve necessary data from other classes.
    //
    // Get references to surface radiation fluxes
    TF* sw_dn = radiation.get_surface_radiation_g("sw_down");
    TF* sw_up = radiation.get_surface_radiation_g("sw_up");
    TF* lw_dn = radiation.get_surface_radiation_g("lw_down");
    TF* lw_up = radiation.get_surface_radiation_g("lw_up");

    // Get (near-) surface thermo.
    // Aarrghh, TODO: replace with `get_tmp_xy_g()......`.
    TF* T_bot = tmp1->flux_bot_g;
    TF* T_a = tmp1->grad_bot_g;
    TF* vpd = tmp1->fld_top_g;
    TF* qsat_bot = tmp1->flux_top_g;
    TF* dqsatdT_bot = tmp1->grad_top_g;

    thermo.get_land_surface_fields_g(
        T_bot, T_a, vpd, qsat_bot, dqsatdT_bot);

    // Get (near-) surface buoyancy.
    auto buoy = fields.get_tmp_g();
    thermo.get_buoyancy_surf_g(*buoy);
    const TF db_ref = thermo.get_db_ref();

    // Get basestate vectors.
    TF* rhorefh = thermo.get_basestate_fld_g("rhoh");
    TF* thvrefh = thermo.get_basestate_fld_g("thvh");
    TF* exnrefh = thermo.get_basestate_fld_g("exnerh");
    TF* prefh   = thermo.get_basestate_fld_g("prefh");

    // Get surface precipitation (positive downwards, kg m-2 s-1 = mm s-1)
    auto tmp2 = fields.get_tmp_g();
    TF* rain_rate = tmp2->fld_bot_g;
    microphys.get_surface_rain_rate_g(rain_rate);

    // XY tmp fields for intermediate calculations
    // Aarrghh, TODO: replace with `get_tmp_xy_g()......`.
    TF* f1  = tmp2->flux_bot_g;
    TF* f2  = tmp2->grad_bot_g;
    TF* f2b = tmp2->fld_top_g;
    TF* f3  = tmp2->flux_top_g;
    TF* theta_mean_n = tmp2->grad_top_g;

    const double subdt = timeloop.get_sub_time_step();

    const int iter = timeloop.get_iteration();
    const int subs = timeloop.get_substep();
    const int mpiid = master.get_mpiid();

    //
    // LSM calculations
    //
    lsmk::calc_tile_fractions_g<<<gridGPU, blockGPU>>>(
            tiles.at("veg").fraction_g,
            tiles.at("soil").fraction_g,
            tiles.at("wet").fraction_g,
            fields.ap2d.at("wl")->fld_g,
            c_veg_g, lai_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
    cuda_check_error();

    // Calculate root fraction weighted mean soil water content
    sk::calc_root_weighted_mean_theta_g<<<gridGPU, blockGPU>>>(
            theta_mean_n,
            fields.sps.at("theta")->fld_g,
            soil_index_g,
            root_fraction_g,
            theta_wp_g,
            theta_fc_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            sgd.kstart, sgd.kend,
            gd.icells, gd.ijcells);
    cuda_check_error();

    // Calculate vegetation/soil resistance functions `f`.
    lsmk::calc_resistance_functions_g<<<gridGPU, blockGPU>>>(
            f1, f2, f2b, f3,
            sw_dn,
            fields.sps.at("theta")->fld_g,
            theta_mean_n, vpd,
            gD_coeff_g,
            c_veg_g,
            theta_wp_g,
            theta_fc_g,
            theta_res_g,
            soil_index_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            sgd.kend,
            gd.icells, gd.ijcells);
    cuda_check_error();

    // Calculate canopy resistance for veg and soil tiles.
    lsmk::calc_canopy_resistance_g<<<gridGPU, blockGPU>>>(
            tiles.at("veg").rs_g,
            rs_veg_min_g, lai_g,
            f1, f2, f3,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
    cuda_check_error();

    lsmk::calc_soil_resistance_g<<<gridGPU, blockGPU>>>(
            tiles.at("soil").rs_g,
            rs_soil_min_g, f2b,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
    cuda_check_error();

    // Loop over tiles, and calculate tile properties and fluxes
    for (auto& tile : tiles)
    {
        bool use_cs_veg = (tile.first == "veg");

        //
        // 1) Calculate obuk/ustar/ra using thl_bot and qt_bot
        // from previous time step (= old method, similar to DALES).
        // 2) Calculate new thl_bot such that SEB closes.
        //
        thermo.get_buoyancy_surf_g(
                buoy->fld_bot_g,
                tile.second.thl_bot_g,
                tile.second.qt_bot_g);

        // Calculate Obuk, ustar, and ra.
        if (sw_constant_z0)
            lsmk::calc_stability_g<TF, true><<<gridGPU2, blockGPU2>>>(
                    tile.second.ustar_g,
                    tile.second.obuk_g,
                    tile.second.bfluxbot_g,
                    tile.second.ra_g,
                    tile.second.nobuk_g,
                    du_tot,
                    buoy->fld_g,
                    buoy->fld_bot_g,
                    z0m_g, z0h_g,
                    zL_sl_g,
                    f_sl_g,
                    db_ref,
                    gd.z[gd.kstart],
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart,
                    gd.icells, gd.jcells,
                    gd.ijcells);
        //else
        //    lsmk::calc_stability<TF, false>(
        //            tile.second.ustar.data(),
        //            tile.second.obuk.data(),
        //            tile.second.bfluxbot.data(),
        //            tile.second.ra.data(),
        //            tile.second.nobuk.data(),
        //            (*dutot).data(),
        //            buoy->fld.data(),
        //            buoy->fld_bot.data(),
        //            z0m.data(), z0h.data(),
        //            zL_sl.data(),
        //            f_sl.data(),
        //            db_ref,
        //            gd.z[gd.kstart],
        //            gd.istart, gd.iend,
        //            gd.jstart, gd.jend,
        //            gd.kstart,
        //            gd.icells, gd.jcells,
        //            gd.ijcells);

        //std::cout << "from GPU:" << std::endl;
        //print_ij(tile.second.ustar_g);
        //throw 1;
    }




    fields.release_tmp_g(buoy);
    fields.release_tmp_g(tmp1);
    fields.release_tmp_g(tmp2);
}

template<typename TF>
void Boundary_surface_lsm<TF>::exec_column(Column<TF>& column)
{
}

template<typename TF>
void Boundary_surface_lsm<TF>::print_ij(
    const TF* const __restrict__ fld_g)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;

    int gridi = gd.imax/blocki + (gd.imax%blocki > 0);
    int gridj = gd.jmax/blockj + (gd.jmax%blockj > 0);
    dim3 gridGPU (gridi,  gridj,  1);
    dim3 blockGPU(blocki, blockj, 1);

    lsmk::print_ij<<<gridGPU, blockGPU>>>(
        fld_g,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.icells);
    cuda_check_error();
    cudaDeviceSynchronize();
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

    if (sw_water)
    {
        cuda_safe_call(cudaMalloc(&water_mask_g, int_memsize_ij));
        cuda_safe_call(cudaMalloc(&t_bot_water_g, tf_memsize_ij));
    }

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

    if (sw_water)
    {
        cuda_safe_call(cudaMemcpy(water_mask_g, water_mask.data(), int_memsize_ij, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(t_bot_water_g, t_bot_water.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    }

    cuda_safe_call(cudaMemcpy(interception_g, interception.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(throughfall_g, throughfall.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(infiltration_g, infiltration.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(runoff_g, runoff.data(), tf_memsize_ij, cudaMemcpyHostToDevice));

    // 3. Copy 3D soil properties:
    const int tf_memsize_ijk  = sgd.ncells*sizeof(TF);
    const int int_memsize_ijk = sgd.ncells*sizeof(int);

    cuda_safe_call(cudaMemcpy(soil_index_g, soil_index.data(), int_memsize_ijk, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(diffusivity_g, diffusivity.data(), tf_memsize_ijk, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(diffusivity_h_g, diffusivity_h.data(), tf_memsize_ijk, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(conductivity_g, conductivity.data(), tf_memsize_ijk, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(conductivity_h_g, conductivity_h.data(), tf_memsize_ijk, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(source_g, source.data(), tf_memsize_ijk, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(root_fraction_g, root_fraction.data(), tf_memsize_ijk, cudaMemcpyHostToDevice));

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
