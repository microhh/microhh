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
#include "column.h"

namespace
{
    namespace lsmk = Land_surface_kernels_g;
    namespace bsk = Boundary_surface_kernels_g;
    namespace sk = Soil_kernels_g;
}

namespace
{
    template<typename TF>
    void dump_field(
        TF* const restrict fld,
        TF* const restrict tmp,
        std::string name,
        const int size)
    {
        std::cout << "Saving: " << name << std::endl;
        cuda_safe_call(cudaMemcpy(tmp, fld, size*sizeof(TF), cudaMemcpyDeviceToHost));
        FILE *pFile;
        pFile = fopen(name.c_str(), "wb");
        if (pFile == NULL)
            std::cout << "Error opening file" << std::endl;
        fwrite(tmp, sizeof(TF), size, pFile);
        fclose(pFile);
    }
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
    dim3 grid_gpu_2d (gridi,  gridj,  1);
    dim3 block_gpu_2d(blocki, blockj, 1);

    // For 2D field including ghost cells
    gridi = gd.icells/blocki + (gd.icells%blocki > 0);
    gridj = gd.jcells/blockj + (gd.jcells%blockj > 0);
    dim3 grid_gpu_2d_gc (gridi,  gridj,  1);
    dim3 block_gpu_2d_gc(blocki, blockj, 1);

    // Calculate filtered wind speed difference surface-atmosphere.
    auto tmp1 = fields.get_tmp_g();
    // Aarrghh, TODO: replace with `get_tmp_xy_g()......`.
    TF* du_tot = tmp1->fld_bot_g;

    bsk::calc_dutot_g<<<grid_gpu_2d, block_gpu_2d>>>(
        du_tot,
        fields.mp.at("u")->fld_g,
        fields.mp.at("v")->fld_g,
        fields.mp.at("u")->fld_bot_g,
        fields.mp.at("v")->fld_bot_g,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.kstart,
        gd.icells, gd.ijcells);
    boundary_cyclic.exec_2d_g(du_tot);
    cuda_check_error();

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
    lsmk::calc_tile_fractions_g<<<grid_gpu_2d, block_gpu_2d>>>(
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
    sk::calc_root_weighted_mean_theta_g<<<grid_gpu_2d, block_gpu_2d>>>(
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
    lsmk::calc_resistance_functions_g<<<grid_gpu_2d, block_gpu_2d>>>(
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
    lsmk::calc_canopy_resistance_g<<<grid_gpu_2d, block_gpu_2d>>>(
            tiles.at("veg").rs_g,
            rs_veg_min_g, lai_g,
            f1, f2, f3,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
    cuda_check_error();

    lsmk::calc_soil_resistance_g<<<grid_gpu_2d, block_gpu_2d>>>(
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
            lsmk::calc_stability_g<TF, true><<<grid_gpu_2d_gc, block_gpu_2d_gc>>>(
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
        else
            lsmk::calc_stability_g<TF, false><<<grid_gpu_2d_gc, block_gpu_2d_gc>>>(
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
        cuda_check_error();

        //auto tmp_cpu = fields.get_tmp();
        //dump_field(tile.second.ustar_g, tmp_cpu->fld_bot.data(), "dump_gpu", gd.ijcells);
        //fields.release_tmp(tmp_cpu);
        //cudaDeviceSynchronize();
        //throw 1;

        // Calculate surface fluxes
        lsmk::calc_fluxes_g<<<grid_gpu_2d, block_gpu_2d>>>(
                tile.second.H_g,
                tile.second.LE_g,
                tile.second.G_g,
                tile.second.S_g,
                tile.second.thl_bot_g,
                tile.second.qt_bot_g,
                T_a,
                fields.sp.at("qt")->fld_g,
                fields.sps.at("t")->fld_g,
                qsat_bot, dqsatdT_bot,
                tile.second.ra_g,
                tile.second.rs_g,
                lambda_stable_g,
                lambda_unstable_g,
                cs_veg_g,
                sw_dn,
                sw_up,
                lw_dn,
                lw_up,
                buoy->fld_g,
                buoy->fld_bot_g,
                rhorefh,
                exnrefh,
                db_ref, emis_sfc,
                TF(subdt),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, sgd.kend,
                gd.icells, gd.ijcells,
                use_cs_veg);
        cuda_check_error();
    }

    // Override grid point with water
    if (sw_water)
    {
        // Set BCs for water grid points
        lsmk::set_water_tiles<<<grid_gpu_2d, block_gpu_2d>>>(
                tiles.at("veg").fraction_g,
                tiles.at("soil").fraction_g,
                tiles.at("wet").fraction_g,
                tiles.at("veg").H_g,
                tiles.at("soil").H_g,
                tiles.at("wet").H_g,
                tiles.at("veg").LE_g,
                tiles.at("soil").LE_g,
                tiles.at("wet").LE_g,
                tiles.at("veg").G_g,
                tiles.at("soil").G_g,
                tiles.at("wet").G_g,
                tiles.at("veg").rs_g,
                tiles.at("soil").rs_g,
                tiles.at("wet").rs_g,
                tiles.at("wet").thl_bot_g,
                tiles.at("wet").qt_bot_g,
                water_mask_g,
                t_bot_water_g,
                fields.sp.at("thl")->fld_g,
                fields.sp.at("qt")->fld_g,
                fields.sp.at("thl")->fld_bot_g,
                fields.sp.at("qt")->fld_bot_g,
                tiles.at("wet").ra_g,
                rhorefh,
                prefh,
                exnrefh,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart,
                gd.icells, gd.ijcells);
        cuda_check_error();
    }

    // Calculate tile averaged surface fluxes and values.
    const TF rhoref_bot = thermo.get_basestate_vector("rhoh")[gd.kstart];
    const TF rhocpi = TF(1) / (rhoref_bot * Constants::cp<TF>);
    const TF rholvi = TF(1) / (rhoref_bot * Constants::Lv<TF>);
    const TF no_scaling = TF(1);

    // Surface fluxes.
    get_tiled_mean_g(fields.sp.at("thl")->flux_bot_g, "H", rhocpi);
    get_tiled_mean_g(fields.sp.at("qt")->flux_bot_g, "LE", rholvi);
    get_tiled_mean_g(ustar_g, "ustar", no_scaling);
    get_tiled_mean_g(buoy->flux_bot_g, "bfluxbot", no_scaling);

    // Surface values.
    get_tiled_mean_g(fields.sp.at("thl")->fld_bot_g, "thl_bot", TF(1));
    get_tiled_mean_g(fields.sp.at("qt")->fld_bot_g, "qt_bot", TF(1));

    // Set ghost cells `thl_bot`, `qt_bot`, needed for surface scheme
    boundary_cyclic.exec_2d_g(fields.sp.at("thl")->fld_bot_g);
    boundary_cyclic.exec_2d_g(fields.sp.at("qt")->fld_bot_g);

    // Calculate bulk Obukhov length.
    lsmk::calc_bulk_obuk_g<<<grid_gpu_2d, block_gpu_2d>>>(
            obuk_g,
            buoy->flux_bot_g,
            ustar_g,
            gd.z[gd.kstart],
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
    cuda_check_error();

    boundary_cyclic.exec_2d_g(ustar_g);
    boundary_cyclic.exec_2d_g(obuk_g);

    // Redistribute ustar over `uw` and `vw`.
    lsmk::set_bcs_momentum_g<<<grid_gpu_2d, block_gpu_2d>>>(
            fields.mp.at("u")->flux_bot_g,
            fields.mp.at("v")->flux_bot_g,
            fields.mp.at("u")->grad_bot_g,
            fields.mp.at("v")->grad_bot_g,
            ustar_g,
            fields.mp.at("u")->fld_g,
            fields.mp.at("u")->fld_bot_g,
            fields.mp.at("v")->fld_g,
            fields.mp.at("v")->fld_bot_g,
            z0m_g, gd.z[gd.kstart],
            gd.istart, gd.iend,
            gd.jstart, gd.jend, gd.kstart,
            gd.icells, gd.jcells, gd.ijcells);
    cuda_check_error();

    boundary_cyclic.exec_2d_g(fields.mp.at("u")->flux_bot_g);
    boundary_cyclic.exec_2d_g(fields.mp.at("v")->flux_bot_g);
    boundary_cyclic.exec_2d_g(fields.mp.at("u")->grad_bot_g);
    boundary_cyclic.exec_2d_g(fields.mp.at("v")->grad_bot_g);

    // Set BCs (gradients) thl + qt
    lsmk::set_bcs_thl_qt_g<<<grid_gpu_2d_gc, block_gpu_2d_gc>>>(
            fields.sp.at("thl")->grad_bot_g,
            fields.sp.at("qt")->grad_bot_g,
            fields.sp.at("thl")->fld_g,
            fields.sp.at("qt")->fld_g,
            fields.sp.at("thl")->fld_bot_g,
            fields.sp.at("qt")->fld_bot_g,
            gd.z[gd.kstart], gd.kstart,
            gd.icells, gd.jcells, gd.ijcells);
    cuda_check_error();

    // Set BCs other scalars
    for (auto& it : fields.sp)
        if (it.first != "thl" and it.first != "qt")
        {
            if (sbc.at(it.first).bcbot == Boundary_type::Dirichlet_type)
                lsmk::set_bcs_scalars_dirichlet_g<<<grid_gpu_2d_gc, block_gpu_2d_gc>>>(
                    it.second->fld_bot_g,
                    it.second->grad_bot_g,
                    it.second->flux_bot_g,
                    ustar_g, obuk_g,
                    it.second->fld_g, z0h_g,
                    gd.z[gd.kstart],
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend, gd.kstart,
                    gd.icells, gd.jcells, gd.ijcells);

            else if (sbc.at(it.first).bcbot == Boundary_type::Flux_type)
                lsmk::set_bcs_scalars_flux_g<<<grid_gpu_2d_gc, block_gpu_2d_gc>>>(
                    it.second->fld_bot_g,
                    it.second->grad_bot_g,
                    it.second->flux_bot_g,
                    ustar_g, obuk_g,
                    it.second->fld_g, z0h_g,
                    gd.z[gd.kstart],
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend, gd.kstart,
                    gd.icells, gd.jcells, gd.ijcells);
            cuda_check_error();
        }

    // Calc MO gradients, for subgrid scheme
    bsk::calc_duvdz_mo_g<<<grid_gpu_2d_gc, block_gpu_2d_gc>>>(
            dudz_mo_g, dvdz_mo_g,
            fields.mp.at("u")->fld_g,
            fields.mp.at("v")->fld_g,
            fields.mp.at("u")->fld_bot_g,
            fields.mp.at("v")->fld_bot_g,
            fields.mp.at("u")->flux_bot_g,
            fields.mp.at("v")->flux_bot_g,
            ustar_g, obuk_g, z0m_g,
            gd.z[gd.kstart],
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart,
            gd.icells, gd.ijcells);
    cuda_check_error();

    bsk::calc_dbdz_mo_g<<<grid_gpu_2d_gc, block_gpu_2d_gc>>>(
            dbdz_mo_g, buoy->flux_bot_g,
            ustar_g, obuk_g,
            gd.z[gd.kstart],
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
    cuda_check_error();

    // Calculate changes in the liquid water reservoir
    lsmk::calc_liquid_water_reservoir_g<<<grid_gpu_2d, block_gpu_2d>>>(
            fields.at2d.at("wl")->fld_g,
            interception_g,
            throughfall_g,
            fields.ap2d.at("wl")->fld_g,
            tiles.at("veg").LE_g,
            tiles.at("soil").LE_g,
            tiles.at("wet").LE_g,
            tiles.at("veg").fraction_g,
            tiles.at("soil").fraction_g,
            tiles.at("wet").fraction_g,
            rain_rate,
            c_veg_g,
            lai_g, subdt,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
    cuda_check_error();

    fields.release_tmp_g(buoy);
    fields.release_tmp_g(tmp2);

    if (sw_homogenize_sfc)
    {
        const int blockGPU = 256;
        const int gridGPU = gd.ijcells/blockGPU + (gd.ijcells%blockGPU > 0);

        auto homogenize = [&](TF* const __restrict__ field)
        {
            const TF mean_value = field3d_operators.calc_mean_2d_g(field);
            Tools_g::set_to_val<<<gridGPU, blockGPU>>>(field, gd.ijcells, mean_value);
        };

        // Homogenize the surface fields which interact with the atmosphere.
        homogenize(fields.sp.at("thl")->flux_bot_g);
        homogenize(fields.sp.at("qt")->flux_bot_g);

        homogenize(fields.mp.at("u")->flux_bot_g),
        homogenize(fields.mp.at("v")->flux_bot_g),

        homogenize(dudz_mo_g);
        homogenize(dvdz_mo_g);
        homogenize(dbdz_mo_g);
    }

    //
    // Calculate soil tendencies
    //
    // Only soil moisture has a source and conductivity term
    const bool sw_source_term_t = false;
    const bool sw_conductivity_term_t = false;
    const bool sw_source_term_theta = true;
    const bool sw_conductivity_term_theta = true;

    // Soil GPU grid without ghost cells.
    gridi = gd.imax/blocki + (gd.imax%blocki > 0);
    gridj = gd.jmax/blockj + (gd.jmax%blockj > 0);
    dim3 grid_gpu_3d (gridi,  gridj,  sgd.kmax);
    dim3 block_gpu_3d(blocki, blockj, 1);

    //
    // Soil temperature
    //
    // Calculate the thermal diffusivity at full levels
    sk::calc_thermal_properties_g<<<grid_gpu_3d, block_gpu_3d>>>(
            diffusivity_g,
            conductivity_g,
            soil_index_g,
            fields.sps.at("theta")->fld_g,
            theta_sat_g,
            gamma_T_dry_g,
            rho_C_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            sgd.kstart, sgd.kend,
            gd.icells, gd.ijcells);
    cuda_check_error();

    // Linear interpolation diffusivity to half levels
    sk::interp_2_vertical_g<TF, Soil_interpolation_type::Harmonic_mean><<<grid_gpu_3d, block_gpu_3d>>>(
            diffusivity_h_g,
            diffusivity_g,
            sgd.dz_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            sgd.kstart, sgd.kend,
            gd.icells, gd.ijcells);
    cuda_check_error();

    // Set flux boundary conditions at top and bottom of soil column
    // Top = soil heat flux (G) averaged over all tiles, bottom = zero flux.
    get_tiled_mean_g(tmp1->fld_bot_g, "G", TF(1));

    sk::set_bcs_temperature_g<<<grid_gpu_2d, block_gpu_2d>>>(
            fields.sps.at("t")->flux_top_g,
            fields.sps.at("t")->flux_bot_g,
            tmp1->fld_bot_g,
            rho_C_g,
            soil_index_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            sgd.kend,
            gd.icells, gd.ijcells);
    cuda_check_error();

    // Calculate diffusive tendency
    sk::diff_explicit_g<TF, sw_source_term_t, sw_conductivity_term_t><<<grid_gpu_3d, block_gpu_3d>>>(
            fields.sts.at("t")->fld_g,
            fields.sps.at("t")->fld_g,
            diffusivity_h_g,
            conductivity_h_g,
            source_g,
            fields.sps.at("t")->flux_top_g,
            fields.sps.at("t")->flux_bot_g,
            sgd.dzi_g, sgd.dzhi_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            sgd.kstart, sgd.kend,
            gd.icells, gd.ijcells);
    cuda_check_error();

    //
    // Soil moisture
    //
    // Calculate the hydraulic diffusivity and conductivity at full levels
    sk::calc_hydraulic_properties_g<<<grid_gpu_3d, block_gpu_3d>>>(
            diffusivity_g,
            conductivity_g,
            soil_index_g,
            fields.sps.at("theta")->fld_g,
            theta_sat_g,
            theta_res_g,
            vg_a_g,
            vg_l_g,
            vg_m_g,
            gamma_theta_sat_g,
            gamma_theta_min_g,
            gamma_theta_max_g,
            kappa_theta_min_g,
            kappa_theta_max_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            sgd.kstart, sgd.kend,
            gd.icells, gd.ijcells);
    cuda_check_error();

    // Interpolation diffusivity and conductivity to half levels,
    // using the IFS method, which uses the max value from the
    // two surrounding grid points.
    sk::interp_2_vertical_g<TF, Soil_interpolation_type::Max><<<grid_gpu_3d, block_gpu_3d>>>(
            diffusivity_h_g,
            diffusivity_g,
            sgd.dz_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            sgd.kstart, sgd.kend,
            gd.icells, gd.ijcells);
    cuda_check_error();

    sk::interp_2_vertical_g<TF, Soil_interpolation_type::Max><<<grid_gpu_3d, block_gpu_3d>>>(
            conductivity_h_g,
            conductivity_g,
            sgd.dz_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            sgd.kstart, sgd.kend,
            gd.icells, gd.ijcells);
    cuda_check_error();

    // Calculate infiltration/runoff
    sk::calc_infiltration_g<<<grid_gpu_2d, block_gpu_2d>>>(
            infiltration_g,
            runoff_g,
            throughfall_g,
            fields.sps.at("theta")->fld_g,
            theta_sat_g,
            kappa_theta_max_g,
            gamma_theta_max_g,
            sgd.dz_g,
            soil_index_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            sgd.kend,
            gd.icells, gd.ijcells);
    cuda_check_error();

    // Set the boundary conditions.
    // Top = evaporation from bare soil tile.
    // Bottom = optionally free drainage (or else closed)
    sk::set_bcs_moisture_g<<<grid_gpu_2d, block_gpu_2d>>>(
            fields.sps.at("theta")->flux_top_g,
            fields.sps.at("theta")->flux_bot_g,
            conductivity_h_g,
            tiles.at("soil").LE_g,
            tiles.at("soil").fraction_g,
            infiltration_g,
            sw_free_drainage,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            sgd.kstart, sgd.kend,
            gd.icells, gd.ijcells);
    cuda_check_error();

    // Calculate root water extraction
    lsmk::scale_tile_with_fraction_g<<<grid_gpu_2d, block_gpu_2d>>>(
            tmp1->fld_bot_g,
            tiles.at("veg").LE_g,
            tiles.at("veg").fraction_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
    cuda_check_error();

    sk::calc_root_water_extraction_g<<<grid_gpu_2d, block_gpu_2d>>>(
            source_g,
            tmp1->fld_top_g,
            fields.sps.at("theta")->fld_g,
            root_fraction_g,
            tmp1->fld_bot_g,
            sgd.dzi_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            sgd.kstart, sgd.kend,
            gd.icells, gd.ijcells);
    cuda_check_error();

    // Calculate diffusive tendency
    sk::diff_explicit_g<TF, sw_source_term_theta, sw_conductivity_term_theta><<<grid_gpu_3d, block_gpu_3d>>>(
            fields.sts.at("theta")->fld_g,
            fields.sps.at("theta")->fld_g,
            diffusivity_h_g,
            conductivity_h_g,
            source_g,
            fields.sps.at("theta")->flux_top_g,
            fields.sps.at("theta")->flux_bot_g,
            sgd.dzi_g, sgd.dzhi_g,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            sgd.kstart, sgd.kend,
            gd.icells, gd.ijcells);

    fields.release_tmp_g(tmp1);
}

template<typename TF>
void Boundary_surface_lsm<TF>::exec_column(Column<TF>& column)
{
    const TF no_offset = 0.;

    auto tmp = fields.get_tmp_g();

    column.calc_time_series("obuk", obuk_g, no_offset);
    column.calc_time_series("ustar", ustar_g, no_offset);
    column.calc_time_series("wl", fields.ap2d.at("wl")->fld_g, no_offset);

    get_tiled_mean_g(tmp->fld_bot_g, "H", TF(1));
    column.calc_time_series("H", tmp->fld_bot_g, no_offset);

    get_tiled_mean_g(tmp->fld_bot_g, "LE", TF(1));
    column.calc_time_series("LE", tmp->fld_bot_g, no_offset);

    get_tiled_mean_g(tmp->fld_bot_g, "G", TF(1));
    column.calc_time_series("G", tmp->fld_bot_g, no_offset);

    get_tiled_mean_g(tmp->fld_bot_g, "S", TF(1));
    column.calc_time_series("S", tmp->fld_bot_g, no_offset);

    if (sw_tile_stats_col)
        for (auto& tile : tiles)
        {
            column.calc_time_series("c_"+tile.first, tile.second.fraction_g, no_offset);

            column.calc_time_series("ustar_"+tile.first, tile.second.ustar_g, no_offset);
            column.calc_time_series("obuk_"+tile.first, tile.second.obuk_g, no_offset);

            column.calc_time_series("rs_"+tile.first, tile.second.rs_g, no_offset);
            column.calc_time_series("ra_"+tile.first, tile.second.ra_g, no_offset);

            column.calc_time_series("thl_bot_"+tile.first, tile.second.thl_bot_g, no_offset);
            column.calc_time_series("qt_bot_"+tile.first, tile.second.qt_bot_g, no_offset);

            column.calc_time_series("H_"+tile.first, tile.second.H_g, no_offset);
            column.calc_time_series("LE_"+tile.first, tile.second.LE_g, no_offset);
            column.calc_time_series("G_"+tile.first, tile.second.G_g, no_offset);
            column.calc_time_series("S_"+tile.first, tile.second.S_g, no_offset);
        }

    fields.release_tmp_g(tmp);
}

template<typename TF>
void Boundary_surface_lsm<TF>::get_tiled_mean_g(
    TF* const __restrict__ fld_out, std::string name, const TF scale_factor)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;

    // For 2D field excluding ghost cells
    int gridi = gd.imax/blocki + (gd.imax%blocki > 0);
    int gridj = gd.jmax/blockj + (gd.jmax%blockj > 0);
    dim3 gridGPU (gridi,  gridj,  1);
    dim3 blockGPU(blocki, blockj, 1);

    TF* fld_veg;
    TF* fld_soil;
    TF* fld_wet;

    // Yikes..
    if (name == "H")
    {
        fld_veg  = tiles.at("veg").H_g;
        fld_soil = tiles.at("soil").H_g;
        fld_wet  = tiles.at("wet").H_g;
    }
    else if (name == "LE")
    {
        fld_veg  = tiles.at("veg").LE_g;
        fld_soil = tiles.at("soil").LE_g;
        fld_wet  = tiles.at("wet").LE_g;
    }
    else if (name == "G")
    {
        fld_veg  = tiles.at("veg").G_g;
        fld_soil = tiles.at("soil").G_g;
        fld_wet  = tiles.at("wet").G_g;
    }
    else if (name == "S")
    {
        fld_veg  = tiles.at("veg").S_g;
        fld_soil = tiles.at("soil").S_g;
        fld_wet  = tiles.at("wet").S_g;
    }
    else if (name == "bfluxbot")
    {
        fld_veg  = tiles.at("veg").bfluxbot_g;
        fld_soil = tiles.at("soil").bfluxbot_g;
        fld_wet  = tiles.at("wet").bfluxbot_g;
    }
    else if (name == "ustar")
    {
        fld_veg  = tiles.at("veg").ustar_g;
        fld_soil = tiles.at("soil").ustar_g;
        fld_wet  = tiles.at("wet").ustar_g;
    }
    else if (name == "thl_bot")
    {
        fld_veg  = tiles.at("veg").thl_bot_g;
        fld_soil = tiles.at("soil").thl_bot_g;
        fld_wet  = tiles.at("wet").thl_bot_g;
    }
    else if (name == "qt_bot")
    {
        fld_veg  = tiles.at("veg").qt_bot_g;
        fld_soil = tiles.at("soil").qt_bot_g;
        fld_wet  = tiles.at("wet").qt_bot_g;
    }
    else
        throw std::runtime_error("Cannot calculate tiled mean for variable \"" + name + "\"\\n");

    lsmk::calc_tiled_mean_g<<<gridGPU, blockGPU>>>(
            fld_out,
            tiles.at("veg").fraction_g,
            tiles.at("soil").fraction_g,
            tiles.at("wet").fraction_g,
            fld_veg,
            fld_soil,
            fld_wet,
            scale_factor,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
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

    // Prepare base boundary, for inflow profiles.
    Boundary<TF>::prepare_device();

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

    if (sw_constant_z0)
    {
        cuda_safe_call(cudaMalloc(&nobuk_g, int_memsize_ij));
        cuda_safe_call(cudaMalloc(&zL_sl_g, float_memsize_mo_lut));
        cuda_safe_call(cudaMalloc(&f_sl_g,  float_memsize_mo_lut));
    }

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
    const int tf_memsizeh_ijk = sgd.ncellsh*sizeof(TF);
    const int int_memsize_ijk = sgd.ncells*sizeof(int);

    cuda_safe_call(cudaMalloc(&soil_index_g, int_memsize_ijk));
    cuda_safe_call(cudaMalloc(&diffusivity_g, tf_memsize_ijk));
    cuda_safe_call(cudaMalloc(&diffusivity_h_g, tf_memsizeh_ijk));
    cuda_safe_call(cudaMalloc(&conductivity_g, tf_memsize_ijk));
    cuda_safe_call(cudaMalloc(&conductivity_h_g, tf_memsizeh_ijk));
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

    if (sw_constant_z0)
    {
        cuda_safe_call(cudaMemcpy(nobuk_g, nobuk.data(), int_memsize_ij, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(zL_sl_g, zL_sl.data(), float_memsize_lut, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(f_sl_g,  f_sl.data(),  float_memsize_lut, cudaMemcpyHostToDevice));
    }

    // Land-surface:
    // 1. Copy tiles:
    for (auto& tile : tiles)
        lsmk::forward_device_tile(tile.second, gd.ijcells);

    // 2. Copy 2D surface properties:
    cuda_safe_call(cudaMemcpy(gD_coeff_g, gD_coeff.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
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

    // TODO: which fields are needed from the land-surface?
    // Nearly all tile fields are used in the statistics:
    for (auto& tile : tiles)
        lsmk::backward_device_tile(tile.second, gd.ijcells);
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

    if (sw_constant_z0)
    {
        cuda_safe_call(cudaFree(nobuk_g));
        cuda_safe_call(cudaFree(zL_sl_g));
        cuda_safe_call(cudaFree(f_sl_g));
    }
    // Land-surface stuff:
}
#endif

template class Boundary_surface_lsm<double>;
template class Boundary_surface_lsm<float>;
