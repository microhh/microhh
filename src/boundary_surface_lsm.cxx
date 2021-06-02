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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "boundary_surface_lsm.h"
#include "boundary.h"

#include "master.h"
#include "input.h"
#include "grid.h"
#include "soil_grid.h"
#include "fields.h"
#include "soil_field3d.h"
#include "diff.h"
#include "defines.h"
#include "constants.h"
#include "thermo.h"
#include "master.h"
#include "cross.h"
#include "column.h"
#include "monin_obukhov.h"
#include "fast_math.h"
#include "netcdf_interface.h"
#include "radiation.h"
#include "microphys.h"

#include "boundary_surface_kernels.h"
#include "land_surface_kernels.h"
#include "soil_kernels.h"

namespace
{
    // Make a shortcut in the file scope.
    namespace most = Monin_obukhov;
    namespace fm = Fast_math;
    namespace bsk = Boundary_surface_kernels;
    namespace lsmk = Land_surface_kernels;
    namespace sk = Soil_kernels;
}

namespace
{
    template<typename TF>
    void calc_dutot(
            TF* const restrict dutot,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict ubot,
            const TF* const restrict vbot,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int ijcells)
    {
        const int ii = 1;
        const int ii2 = 2;
        const int jj = icells;
        const int jj2 = 2*icells;
        const int kk = ijcells;

        // Calculate total wind.
        const TF minval = 1.e-1;

        // Interpolate the wind (3x3 mean) to the scalar location.
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;

                const TF u_filtered = TF(1./9) *
                    ( TF(0.5)*u[ijk-ii-jj] + u[ijk-jj] + u[ijk+ii-jj] + TF(0.5)*u[ijk+ii2-jj]
                    + TF(0.5)*u[ijk-ii   ] + u[ijk   ] + u[ijk+ii   ] + TF(0.5)*u[ijk+ii2   ]
                    + TF(0.5)*u[ijk-ii+jj] + u[ijk+jj] + u[ijk+ii+jj] + TF(0.5)*u[ijk+ii2+jj] );

                const TF v_filtered = TF(1./9) *
                    ( TF(0.5)*v[ijk-ii-jj] + v[ijk-ii] + v[ijk-ii+jj] + TF(0.5)*v[ijk-ii+jj2]
                    + TF(0.5)*v[ijk   -jj] + v[ijk   ] + v[ijk   +jj] + TF(0.5)*v[ijk   +jj2]
                    + TF(0.5)*v[ijk+ii-jj] + v[ijk+ii] + v[ijk+ii+jj] + TF(0.5)*v[ijk+ii+jj2] );

                const TF du2 = fm::pow2(u_filtered - TF(0.5)*(ubot[ij] + ubot[ij+ii]))
                             + fm::pow2(v_filtered - TF(0.5)*(vbot[ij] + vbot[ij+jj]));

                // Prevent the absolute wind gradient from reaching values less than 0.01 m/s,
                // otherwise evisc at k = kstart blows up
                dutot[ij] = std::max(std::pow(du2, TF(0.5)), minval);
            }
    }

    template<typename TF>
    void calc_tiled_mean(
            TF* const restrict fld,
            const TF* const restrict c_veg,
            const TF* const restrict c_soil,
            const TF* const restrict c_wet,
            const TF* const restrict fld_veg,
            const TF* const restrict fld_soil,
            const TF* const restrict fld_wet,
            const TF fac,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;

                fld[ij] = (
                    c_veg [ij] * fld_veg [ij] +
                    c_soil[ij] * fld_soil[ij] +
                    c_wet [ij] * fld_wet [ij] ) * fac;
            }
    }

    template<typename TF>
    void calc_bulk_obuk(
            TF* const restrict obuk,
            const TF* const restrict bfluxbot,
            const TF* const restrict ustar,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                obuk[ij] = -fm::pow3(ustar[ij]) / (Constants::kappa<TF> * bfluxbot[ij]);
                obuk[ij] = zsl/std::min(std::max(zsl/obuk[ij], bsk::zL_min<TF>), bsk::zL_max<TF>);
            }
    }

    template<typename TF>
    void set_bcs_momentum(
            TF* const restrict ufluxbot,
            TF* const restrict vfluxbot,
            TF* const restrict ugradbot,
            TF* const restrict vgradbot,
            const TF* const restrict ustar,
            const TF* const restrict u,
            const TF* const restrict ubot,
            const TF* const restrict v,
            const TF* const restrict vbot,
            const TF* const restrict z0m,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int jcells,
            const int ijcells,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int ii = 1;
        const int jj = icells;

        const TF minval = 1.e-2;

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*ijcells;

                const TF vonu2 = std::max(minval, TF(0.25)*(
                            fm::pow2(v[ijk-ii]-vbot[ij-ii]) + fm::pow2(v[ijk-ii+jj]-vbot[ij-ii+jj])
                          + fm::pow2(v[ijk   ]-vbot[ij   ]) + fm::pow2(v[ijk   +jj]-vbot[ij   +jj])) );
                const TF uonv2 = std::max(minval, TF(0.25)*(
                            fm::pow2(u[ijk-jj]-ubot[ij-jj]) + fm::pow2(u[ijk+ii-jj]-ubot[ij+ii-jj])
                          + fm::pow2(u[ijk   ]-ubot[ij   ]) + fm::pow2(u[ijk+ii   ]-ubot[ij+ii   ])) );

                const TF u2 = std::max(minval, fm::pow2(u[ijk]-ubot[ij]) );
                const TF v2 = std::max(minval, fm::pow2(v[ijk]-vbot[ij]) );

                const TF ustaronu4 = TF(0.5)*(fm::pow4(ustar[ij-ii]) + fm::pow4(ustar[ij]));
                const TF ustaronv4 = TF(0.5)*(fm::pow4(ustar[ij-jj]) + fm::pow4(ustar[ij]));

                ufluxbot[ij] = -copysign(TF(1), u[ijk]-ubot[ij]) * std::pow(ustaronu4 / (TF(1) + vonu2 / u2), TF(0.5));
                vfluxbot[ij] = -copysign(TF(1), v[ijk]-vbot[ij]) * std::pow(ustaronv4 / (TF(1) + uonv2 / v2), TF(0.5));
            }

        boundary_cyclic.exec_2d(ufluxbot);
        boundary_cyclic.exec_2d(vfluxbot);

        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*ijcells;

                // Use the linearly interpolated grad, rather than the MO grad,
                // to prevent giving unresolvable gradients to advection schemes
                ugradbot[ij] = (u[ijk]-ubot[ij])/zsl;
                vgradbot[ij] = (v[ijk]-vbot[ij])/zsl;
            }
    }

    template<typename TF>
    void set_bcs_thl_qt(
            TF* const restrict thl_gradbot,
            TF* const restrict qt_gradbot,
            const TF* const restrict thl,
            const TF* const restrict qt,
            const TF* const restrict thl_bot,
            const TF* const restrict qt_bot,
            const TF zsl, const int kstart,
            const int icells, const int jcells,
            const int ijcells)
    {
        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ij  = i + j*icells;
                const int ijk = i + j*icells + kstart*ijcells;

                // Use the linearly interpolated grad, rather than the MO grad,
                // to prevent giving unresolvable gradients to advection schemes
                thl_gradbot[ij] = (thl[ijk]-thl_bot[ij])/zsl;
                qt_gradbot[ij]  = (qt [ijk]-qt_bot [ij])/zsl;
            }
    }
}

template<typename TF>
Boundary_surface_lsm<TF>::Boundary_surface_lsm(
        Master& masterin, Grid<TF>& gridin, Soil_grid<TF>& soilgridin,
        Fields<TF>& fieldsin, Input& inputin) :
        Boundary<TF>(masterin, gridin, soilgridin, fieldsin, inputin)
{
    swboundary = "surface_lsm";

    // Read .ini settings:
    sw_constant_z0   = inputin.get_item<bool>("boundary", "swconstantz0", "", true);
    sw_homogeneous   = inputin.get_item<bool>("land_surface", "swhomogeneous", "", true);
    sw_free_drainage = inputin.get_item<bool>("land_surface", "swfreedrainage", "", true);
    sw_water         = inputin.get_item<bool>("land_surface", "swwater", "", false);
    sw_tile_stats    = inputin.get_item<bool>("land_surface", "swtilestats", "", false);

    if (sw_water)
        tskin_water = inputin.get_item<TF>("land_surface", "tskin_water", "");

    // Create prognostic 2D and 3D fields;
    fields.init_prognostic_soil_field("t", "Soil temperature", "K");
    fields.init_prognostic_soil_field("theta", "Soil volumetric water content", "m3 m-3");
    fields.init_prognostic_2d_field("wl");

    // Create surface tiles:
    for (auto& name : tile_names)
        tiles.emplace(name, Surface_tile<TF>{});

    // Open NetCDF file with soil lookup table:
    nc_lookup_table =
        std::make_shared<Netcdf_file>(master, "van_genuchten_parameters.nc", Netcdf_mode::Read);

    // Checks:
    if (sw_homogeneous && sw_water)
        throw std::runtime_error("Homogeneous land-surface with water is not supported!\n");

    //#ifdef USECUDA
    //ustar_g = 0;
    //obuk_g  = 0;
    //nobuk_g = 0;
    //zL_sl_g = 0;
    //f_sl_g  = 0;
    //#endif
}

template<typename TF>
Boundary_surface_lsm<TF>::~Boundary_surface_lsm()
{
    //#ifdef USECUDA
    //clear_device();
    //#endif
}

#ifndef USECUDA
template<typename TF>
void Boundary_surface_lsm<TF>::exec(
        Thermo<TF>& thermo, Radiation<TF>& radiation,
        Microphys<TF>& microphys, Timeloop<TF>& timeloop)
{
    auto& gd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    //
    // Calculate tile independant properties
    //
    auto dutot = fields.get_tmp_xy();

    calc_dutot(
            (*dutot).data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("u")->fld_bot.data(),
            fields.mp.at("v")->fld_bot.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart,
            gd.icells, gd.ijcells);

    boundary_cyclic.exec_2d((*dutot).data());

    //
    // Retrieve necessary data from other classes.
    //
    // Get references to surface radiation fluxes
    std::vector<TF>& sw_dn = radiation.get_surface_radiation("sw_down");
    std::vector<TF>& sw_up = radiation.get_surface_radiation("sw_up");
    std::vector<TF>& lw_dn = radiation.get_surface_radiation("lw_down");
    std::vector<TF>& lw_up = radiation.get_surface_radiation("lw_up");

    // Get (near-) surface thermo
    auto T_bot = fields.get_tmp_xy();
    auto T_a = fields.get_tmp_xy();
    auto vpd = fields.get_tmp_xy();
    auto qsat_bot = fields.get_tmp_xy();
    auto dqsatdT_bot = fields.get_tmp_xy();

    thermo.get_land_surface_fields(
        *T_bot, *T_a, *vpd, *qsat_bot, *dqsatdT_bot);

    // NOTE: `get_buoyancy_surf` calculates the first model level buoyancy only,
    //       but since this is written at `kstart`, we can't use a 2D slice...
    //auto b_a = fields.get_tmp_xy();
    auto buoy = fields.get_tmp();
    auto b_bot = fields.get_tmp_xy();

    thermo.get_buoyancy_surf(buoy->fld, *b_bot, false);
    const TF db_ref = thermo.get_db_ref();

    const std::vector<TF>& rhorefh = thermo.get_basestate_vector("rhoh");
    const std::vector<TF>& thvrefh = thermo.get_basestate_vector("thvh");
    const std::vector<TF>& exnrefh = thermo.get_basestate_vector("exnerh");
    const std::vector<TF>& prefh = thermo.get_basestate_vector("ph");

    // Get surface precipitation (positive downwards, kg m-2 s-1 = mm s-1)
    auto rain_rate = fields.get_tmp_xy();
    microphys.get_surface_rain_rate(*rain_rate);

    // XY tmp fields for intermediate calculations
    auto f1  = fields.get_tmp_xy();
    auto f2  = fields.get_tmp_xy();
    auto f2b = fields.get_tmp_xy();
    auto f3  = fields.get_tmp_xy();
    auto theta_mean_n = fields.get_tmp_xy();
    auto ra = fields.get_tmp_xy();

    const double subdt = timeloop.get_sub_time_step();

    //
    // LSM calculations
    // Calculate dynamic tile fractions
    lsmk::calc_tile_fractions(
            tiles.at("veg").fraction.data(),
            tiles.at("soil").fraction.data(),
            tiles.at("wet").fraction.data(),
            fields.ap2d.at("wl").data(),
            c_veg.data(),
            lai.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

    // Calculate root fraction weighted mean soil water content
    sk::calc_root_weighted_mean_theta(
            (*theta_mean_n).data(),
            fields.sps.at("theta")->fld.data(),
            soil_index.data(),
            root_fraction.data(),
            theta_wp.data(),
            theta_fc.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            sgd.kstart, sgd.kend,
            gd.icells, gd.ijcells);

    // Calculate vegetation/soil resistance functions `f`
    lsmk::calc_resistance_functions(
            (*f1).data(), (*f2).data(),
            (*f2b).data(), (*f3).data(),
            sw_dn.data(),
            fields.sps.at("theta")->fld.data(),
            (*theta_mean_n).data(),
            (*vpd).data(),
            gD_coeff.data(),
            c_veg.data(),
            theta_wp.data(),
            theta_fc.data(),
            theta_res.data(),
            soil_index.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            sgd.kend,
            gd.icells, gd.ijcells);

    // Calculate canopy resistance per tile
    lsmk::calc_canopy_resistance(
            tiles.at("veg").rs.data(),
            rs_veg_min.data(), lai.data(),
            (*f1).data(), (*f2).data(), (*f3).data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

    lsmk::calc_soil_resistance(
            tiles.at("soil").rs.data(),
            rs_soil_min.data(), (*f2b).data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

    // Loop over tiles
    for (auto& tile : tiles)
    {
        bool use_cs_veg = (tile.first == "veg");

        // Solve SEB and SL combined
        lsmk::calc_stability_and_fluxes(
                tile.second.H.data(),
                tile.second.LE.data(),
                tile.second.G.data(),
                tile.second.S.data(),
                tile.second.thl_bot.data(),
                tile.second.qt_bot.data(),
                tile.second.ustar.data(),
                tile.second.obuk.data(),
                tile.second.bfluxbot.data(),
                sw_dn.data(), sw_up.data(),
                lw_dn.data(), lw_up.data(),
                (*dutot).data(),
                (*T_a).data(),
                fields.sp.at("qt")->fld.data(),
                buoy->fld.data(),
                fields.sps.at("t")->fld.data(),
                tile.second.rs.data(),
                z0m.data(), z0h.data(),
                lambda_stable.data(),
                lambda_unstable.data(),
                cs_veg.data(),
                rhorefh.data(),
                exnrefh.data(),
                thvrefh.data(),
                prefh.data(),
                db_ref, gd.z[gd.kstart],
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, sgd.kend,
                gd.icells, gd.ijcells,
                use_cs_veg,
                tile.first);
    }

    // Calculate tile averaged surface fluxes and values.
    const TF rhocpi = TF(1) / (rhorefh[gd.kstart] * Constants::cp<TF>);
    const TF rholvi = TF(1) / (rhorefh[gd.kstart] * Constants::Lv<TF>);

    // Surface fluxes.
    calc_tiled_mean(
            fields.sp.at("thl")->flux_bot.data(),
            tiles.at("veg").fraction.data(),
            tiles.at("soil").fraction.data(),
            tiles.at("wet").fraction.data(),
            tiles.at("veg").H.data(),
            tiles.at("soil").H.data(),
            tiles.at("wet").H.data(),
            rhocpi,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

    calc_tiled_mean(
            fields.sp.at("qt")->flux_bot.data(),
            tiles.at("veg").fraction.data(),
            tiles.at("soil").fraction.data(),
            tiles.at("wet").fraction.data(),
            tiles.at("veg").LE.data(),
            tiles.at("soil").LE.data(),
            tiles.at("wet").LE.data(),
            rholvi,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

    calc_tiled_mean(
            ustar.data(),
            tiles.at("veg").fraction.data(),
            tiles.at("soil").fraction.data(),
            tiles.at("wet").fraction.data(),
            tiles.at("veg").ustar.data(),
            tiles.at("soil").ustar.data(),
            tiles.at("wet").ustar.data(),
            TF(1),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

    calc_tiled_mean(
            buoy->flux_bot.data(),
            tiles.at("veg").fraction.data(),
            tiles.at("soil").fraction.data(),
            tiles.at("wet").fraction.data(),
            tiles.at("veg").bfluxbot.data(),
            tiles.at("soil").bfluxbot.data(),
            tiles.at("wet").bfluxbot.data(),
            TF(1),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

    // Surface values.
    calc_tiled_mean(
            fields.sp.at("thl")->fld_bot.data(),
            tiles.at("veg").fraction.data(),
            tiles.at("soil").fraction.data(),
            tiles.at("wet").fraction.data(),
            tiles.at("veg").thl_bot.data(),
            tiles.at("soil").thl_bot.data(),
            tiles.at("wet").thl_bot.data(),
            TF(1),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

    calc_tiled_mean(
            fields.sp.at("qt")->fld_bot.data(),
            tiles.at("veg").fraction.data(),
            tiles.at("soil").fraction.data(),
            tiles.at("wet").fraction.data(),
            tiles.at("veg").qt_bot.data(),
            tiles.at("soil").qt_bot.data(),
            tiles.at("wet").qt_bot.data(),
            TF(1),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

    // Calculate bulk Obukhov length.
    calc_bulk_obuk(
            obuk.data(),
            buoy->flux_bot.data(),
            ustar.data(),
            gd.z[gd.kstart],
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

    // Redistribute ustar over `uw` and `vw`.
    set_bcs_momentum(
            fields.mp.at("u")->flux_bot.data(),
            fields.mp.at("v")->flux_bot.data(),
            fields.mp.at("u")->grad_bot.data(),
            fields.mp.at("v")->grad_bot.data(),
            ustar.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("u")->fld_bot.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("v")->fld_bot.data(),
            z0m.data(), gd.z[gd.kstart],
            gd.istart, gd.iend,
            gd.jstart, gd.jend, gd.kstart,
            gd.icells, gd.jcells, gd.ijcells,
            boundary_cyclic);

    // Set BCs (gradients) thl + qt
    set_bcs_thl_qt(
            fields.sp.at("thl")->grad_bot.data(),
            fields.sp.at("qt")->grad_bot.data(),
            fields.sp.at("thl")->fld.data(),
            fields.sp.at("qt")->fld.data(),
            fields.sp.at("thl")->fld_bot.data(),
            fields.sp.at("qt")->fld_bot.data(),
            gd.z[gd.kstart], gd.kstart,
            gd.icells, gd.jcells, gd.ijcells);


    fields.release_tmp_xy(dutot);

    fields.release_tmp_xy(T_bot);
    fields.release_tmp_xy(T_a);
    fields.release_tmp_xy(vpd);
    fields.release_tmp_xy(qsat_bot);
    fields.release_tmp_xy(dqsatdT_bot);

    fields.release_tmp(buoy);
    fields.release_tmp_xy(b_bot);

    fields.release_tmp_xy(rain_rate);

    fields.release_tmp_xy(f1);
    fields.release_tmp_xy(f2);
    fields.release_tmp_xy(f2b);
    fields.release_tmp_xy(f3);
    fields.release_tmp_xy(theta_mean_n);
    fields.release_tmp_xy(ra);
}
#endif

template<typename TF>
void Boundary_surface_lsm<TF>::init(Input& inputin, Thermo<TF>& thermo)
{
    // Process the boundary conditions now all fields are registered.
    process_bcs(inputin);

    // Read and check the boundary_surface specific settings.
    process_input(inputin, thermo);

    // Allocate and initialize the 2D surface fields.
    init_surface_layer(inputin);
    init_land_surface();

    // Initialize the boundary cyclic.
    boundary_cyclic.init();
}

template<typename TF>
void Boundary_surface_lsm<TF>::process_input(Input& inputin, Thermo<TF>& thermo)
{
    // Check whether the prognostic thermo vars are of the same type
    std::vector<std::string> thermolist;
    thermo.get_prog_vars(thermolist);

    // Boundary_surface_lsm only supports Dirichlet BCs
    if (mbcbot != Boundary_type::Dirichlet_type)
        throw std::runtime_error("swboundary=surface_lsm requires mbcbot=noslip");

    for (auto& it : sbc)
        if (it.second.bcbot != Boundary_type::Dirichlet_type)
            throw std::runtime_error("swboundary_surface_lsm requires sbcbot=dirichlet");

    thermobc = Boundary_type::Dirichlet_type;
}

template<typename TF>
void Boundary_surface_lsm<TF>::init_surface_layer(Input& input)
{
    auto& gd = grid.get_grid_data();

    obuk.resize(gd.ijcells);
    ustar.resize(gd.ijcells);

    dudz_mo.resize(gd.ijcells);
    dvdz_mo.resize(gd.ijcells);
    dbdz_mo.resize(gd.ijcells);

    z0m.resize(gd.ijcells);
    z0h.resize(gd.ijcells);

    if (sw_constant_z0)
    {
        nobuk.resize(gd.ijcells);

        const TF z0m_hom = input.get_item<TF>("boundary", "z0m", "");
        const TF z0h_hom = input.get_item<TF>("boundary", "z0h", "");

        std::fill(z0m.begin(), z0m.end(), z0m_hom);
        std::fill(z0h.begin(), z0h.end(), z0h_hom);
        std::fill(nobuk.begin(), nobuk.end(), 0);
    }
    // else: z0m and z0h are read from 2D input files in `load()`.

    // Initialize the obukhov length on a small number.
    std::fill(obuk.begin(), obuk.end(), Constants::dsmall);

    // Also initialise ustar at small number, to prevent div/0
    // in calculation surface gradients during cold start.
    std::fill(ustar.begin(), ustar.end(), Constants::dsmall);
}

template<typename TF>
void Boundary_surface_lsm<TF>::init_land_surface()
{
    auto& gd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    // Allocate the surface tiles
    for (auto& tile : tiles)
        lsmk::init_tile(tile.second, gd.ijcells);
    tiles.at("veg" ).long_name = "vegetation";
    tiles.at("soil").long_name = "bare soil";
    tiles.at("wet" ).long_name = "wet skin";

    gD_coeff.resize(gd.ijcells);
    c_veg.resize(gd.ijcells);
    lai.resize(gd.ijcells);
    rs_veg_min.resize(gd.ijcells);
    rs_soil_min.resize(gd.ijcells);
    lambda_stable.resize(gd.ijcells);
    lambda_unstable.resize(gd.ijcells);
    cs_veg.resize(gd.ijcells);

    if (sw_water)
        water_mask.resize(gd.ijcells);

    interception.resize(gd.ijcells);
    throughfall.resize(gd.ijcells);
    infiltration.resize(gd.ijcells);
    runoff.resize(gd.ijcells);

    // Resize the vectors which contain the soil properties
    soil_index.resize(sgd.ncells);
    diffusivity.resize(sgd.ncells);
    diffusivity_h.resize(sgd.ncellsh);
    conductivity.resize(sgd.ncells);
    conductivity_h.resize(sgd.ncellsh);
    source.resize(sgd.ncells);
    root_fraction.resize(sgd.ncells);

    // Resize the lookup table with van Genuchten parameters
    const int size = nc_lookup_table->get_dimension_size("index");

    theta_res.resize(size);
    theta_wp.resize(size);
    theta_fc.resize(size);
    theta_sat.resize(size);

    gamma_theta_sat.resize(size);
    vg_a.resize(size);
    vg_l.resize(size);
    vg_n.resize(size);
    vg_m.resize(size);

    kappa_theta_max.resize(size);
    kappa_theta_min.resize(size);
    gamma_theta_max.resize(size);
    gamma_theta_min.resize(size);

    gamma_T_dry.resize(size);
    rho_C.resize(size);
}

template<typename TF>
void Boundary_surface_lsm<TF>::create_cold_start(Netcdf_handle& input_nc)
{
    auto& agd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    Netcdf_group& soil_group = input_nc.get_group("soil");
    Netcdf_group& init_group = input_nc.get_group("init");

    // Init the soil variables
    if (sw_homogeneous)
    {
        // Read initial profiles from input NetCDF file
        std::vector<TF> t_prof(sgd.ktot);
        std::vector<TF> theta_prof(sgd.ktot);

        soil_group.get_variable(t_prof, "t_soil", {0}, {sgd.ktot});
        soil_group.get_variable(theta_prof, "theta_soil", {0}, {sgd.ktot});

        // Initialise soil as spatially homogeneous
        sk::init_soil_homogeneous(
                fields.sps.at("t")->fld.data(), t_prof.data(),
                agd.istart, agd.iend,
                agd.jstart, agd.jend,
                sgd.kstart, sgd.kend,
                agd.icells, agd.ijcells);

        sk::init_soil_homogeneous(
                fields.sps.at("theta")->fld.data(), theta_prof.data(),
                agd.istart, agd.iend,
                agd.jstart, agd.jend,
                sgd.kstart, sgd.kend,
                agd.icells, agd.ijcells);
    }
    // else: these fields will be provided by the user as binary input files.

    // Initialise the prognostic surface variables, and/or
    // variables which are needed for consistent restarts.
    std::fill(fields.ap2d.at("wl").begin(), fields.ap2d.at("wl").end(), TF(0));

    // Set initial surface potential temperature and humidity to the atmospheric values (...)
    std::vector<TF> thl_1(1);
    std::vector<TF> qt_1(1);

    init_group.get_variable(thl_1, "thl", {0}, {1});
    init_group.get_variable(qt_1,  "qt",  {0}, {1});

    std::fill(
            fields.sp.at("thl")->fld_bot.begin(),
            fields.sp.at("thl")->fld_bot.end(), thl_1[0]);
    std::fill(
            fields.sp.at("qt")->fld_bot.begin(),
            fields.sp.at("qt")->fld_bot.end(), qt_1[0]);

    // Init surface temperature tiles
    for (auto& tile : tiles)
    {
        std::fill(
                tile.second.thl_bot.begin(),
                tile.second.thl_bot.end(), thl_1[0]);

        std::fill(
                tile.second.qt_bot.begin(),
                tile.second.qt_bot.end(), qt_1[0]);
    }

    // Init surface fluxes to some small non-zero value
    std::fill(
            fields.sp.at("thl")->flux_bot.begin(),
            fields.sp.at("thl")->flux_bot.end(), Constants::dsmall);
    std::fill(
            fields.sp.at("qt")->flux_bot.begin(),
            fields.sp.at("qt")->flux_bot.end(), Constants::dsmall);
}

template<typename TF>
void Boundary_surface_lsm<TF>::create(
        Input& input, Netcdf_handle& input_nc,
        Stats<TF>& stats, Column<TF>& column,
        Cross<TF>& cross, Timeloop<TF>& timeloop)
{
    auto& agd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    // Setup statiscics, cross-sections and column statistics
    create_stats(stats, column, cross);

    // Init soil properties
    if (sw_homogeneous)
    {
        Netcdf_group& soil_group = input_nc.get_group("soil");

        // Soil index
        std::vector<int> soil_index_prof(sgd.ktot);
        soil_group.get_variable<int>(soil_index_prof, "index_soil", {0}, {sgd.ktot});

        sk::init_soil_homogeneous<int>(
                soil_index.data(), soil_index_prof.data(),
                agd.istart, agd.iend,
                agd.jstart, agd.jend,
                sgd.kstart, sgd.kend,
                agd.icells, agd.ijcells);

        // Root fraction
        std::vector<TF> root_frac_prof(sgd.ktot);
        soil_group.get_variable<TF>(root_frac_prof, "root_frac", {0}, {sgd.ktot});

        sk::init_soil_homogeneous<TF>(
                root_fraction.data(), root_frac_prof.data(),
                agd.istart, agd.iend,
                agd.jstart, agd.jend,
                sgd.kstart, sgd.kend,
                agd.icells, agd.ijcells);

        // Lambda function to read namelist value, and init 2D field homogeneous.
        auto init_homogeneous = [&](std::vector<TF>& field, std::string name)
        {
            const TF value = input.get_item<TF>("land_surface", name.c_str(), "");
            std::fill(field.begin(), field.begin()+agd.ijcells, value);
        };

        // Land-surface properties
        init_homogeneous(gD_coeff, "gD");
        init_homogeneous(c_veg, "c_veg");
        init_homogeneous(lai, "lai");
        init_homogeneous(rs_veg_min, "rs_veg_min");
        init_homogeneous(rs_soil_min, "rs_soil_min");
        init_homogeneous(lambda_stable, "lambda_stable");
        init_homogeneous(lambda_unstable, "lambda_unstable");
        init_homogeneous(cs_veg, "cs_veg");
    }
    // else: these fields are read from 2D input files in `boundary->load()`.

    // Set the canopy resistance of the liquid water tile at zero
    std::fill(tiles.at("wet").rs.begin(), tiles.at("wet").rs.begin()+agd.ijcells, 0.);

    // Read the lookup table with soil properties
    const int size = nc_lookup_table->get_dimension_size("index");
    nc_lookup_table->get_variable<TF>(theta_res, "theta_res", {0}, {size});
    nc_lookup_table->get_variable<TF>(theta_wp,  "theta_wp",  {0}, {size});
    nc_lookup_table->get_variable<TF>(theta_fc,  "theta_fc",  {0}, {size});
    nc_lookup_table->get_variable<TF>(theta_sat, "theta_sat", {0}, {size});

    nc_lookup_table->get_variable<TF>(gamma_theta_sat, "gamma_sat", {0}, {size});

    nc_lookup_table->get_variable<TF>(vg_a, "alpha", {0}, {size});
    nc_lookup_table->get_variable<TF>(vg_l, "l",     {0}, {size});
    nc_lookup_table->get_variable<TF>(vg_n, "n",     {0}, {size});

    // Calculate derived properties of the lookup table
    sk::calc_soil_properties(
            kappa_theta_min.data(), kappa_theta_max.data(),
            gamma_theta_min.data(), gamma_theta_max.data(), vg_m.data(),
            gamma_T_dry.data(), rho_C.data(),
            vg_a.data(), vg_l.data(), vg_n.data(), gamma_theta_sat.data(),
            theta_res.data(), theta_sat.data(), theta_fc.data(), size);
}

template<typename TF>
void Boundary_surface_lsm<TF>::create_stats(
        Stats<TF>& stats, Column<TF>& column, Cross<TF>& cross)
{
    auto& agd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    const std::string group_name = "land_surface";
    const std::string group_name_tiles = "land_surface_tiles";

    // add variables to the statistics
    if (stats.get_switch())
    {
        // Surface layer
        stats.add_time_series("ustar", "Surface friction velocity", "m s-1", group_name);
        stats.add_time_series("obuk", "Obukhov length", "m", group_name);

        // LSM
        stats.add_time_series("H", "Surface sensible heat flux", "W m-2", group_name);
        stats.add_time_series("LE", "Surface latent heat flux", "W m-2", group_name);
        stats.add_time_series("G", "Surface soil heat flux", "W m-2", group_name);

        stats.add_time_series("rs_veg", "Canopy resistance", "s m-1", group_name);
        stats.add_time_series("rs_soil", "Soil resistance", "s m-1", group_name);

        if (sw_tile_stats)
            for (auto& tile : tiles)
            {
                stats.add_time_series("c_"+tile.first, "Subgrid fraction "+tile.second.long_name, "-", group_name_tiles);

                stats.add_time_series("ustar_"+tile.first, "Surface friction velocity "+tile.second.long_name, "m s-1", group_name_tiles);
                stats.add_time_series("obuk_"+tile.first, "Obukhov length "+tile.second.long_name, "m", group_name_tiles);

                stats.add_time_series("thl_bot_"+tile.first, "Surface potential temperature "+tile.second.long_name, "K", group_name_tiles);
                stats.add_time_series("qt_bot_"+tile.first, "Obukhov specific humidity "+tile.second.long_name, "kg kg-1", group_name_tiles);

                stats.add_time_series("H_"+tile.first, "Surface sensible heat flux "+tile.second.long_name, "W m-2", group_name_tiles);
                stats.add_time_series("LE_"+tile.first, "Surface latent heat flux "+tile.second.long_name, "W m-2", group_name_tiles);
                stats.add_time_series("G_"+tile.first, "Surface soil heat flux "+tile.second.long_name, "W m-2", group_name_tiles);
            }
    }

    if (column.get_switch())
    {
        column.add_time_series("ustar", "Surface friction velocity", "m s-1");
        column.add_time_series("obuk", "Obukhov length", "m");
    }

    if (cross.get_switch())
    {
        const std::vector<std::string> allowed_crossvars = {"ustar", "obuk", "ra"};
        cross_list = cross.get_enabled_variables(allowed_crossvars);
    }
}

template<typename TF>
void Boundary_surface_lsm<TF>::load(const int iotime)
{
    auto& sgd = soil_grid.get_grid_data();

    auto tmp1 = fields.get_tmp();
    auto tmp2 = fields.get_tmp();

    int nerror = 0;
    const TF no_offset = TF(0);

    // Lambda function to load 2D fields.
    auto load_2d_field = [&](
            TF* const restrict field, const std::string& name, const int itime)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", name.c_str(), itime);
        master.print_message("Loading \"%s\" ... ", filename);

        if (field3d_io.load_xy_slice(
                field, tmp1->fld.data(),
                filename))
        {
            master.print_message("FAILED\n");
            nerror += 1;
        }
        else
            master.print_message("OK\n");

        boundary_cyclic.exec_2d(field);
    };

    // Lambda function to load 3D soil fields.
    auto load_3d_field = [&](
            TF* const restrict field, const std::string& name, const int itime)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", name.c_str(), itime);
        master.print_message("Loading \"%s\" ... ", filename);

        if (field3d_io.load_field3d(
                field,
                tmp1->fld.data(), tmp2->fld.data(),
                filename, no_offset,
                sgd.kstart, sgd.kend))
        {
            master.print_message("FAILED\n");
            nerror += 1;
        }
        else
            master.print_message("OK\n");
    };

    // MO gradients are always needed, as the calculation of the
    // eddy viscosity uses the gradients from the previous time step.
    load_2d_field(dudz_mo.data(), "dudz_mo", iotime);
    load_2d_field(dvdz_mo.data(), "dvdz_mo", iotime);
    load_2d_field(dbdz_mo.data(), "dbdz_mo", iotime);

    // Obukhov length restart files are only needed for the iterative solver
    if (!sw_constant_z0)
    {
        // Read Obukhov length
        load_2d_field(obuk.data(), "obuk", iotime);

        // Read spatial z0 fields
        load_2d_field(z0m.data(), "z0m", 0);
        load_2d_field(z0h.data(), "z0h", 0);
    }

    // Load the 3D soil temperature and moisture fields.
    load_3d_field(fields.sps.at("t")->fld.data(), "t_soil", iotime);
    load_3d_field(fields.sps.at("theta")->fld.data(), "theta_soil", iotime);

    // Load the surface temperature, humidity and liquid water content.
    load_2d_field(fields.ap2d.at("wl").data(), "wl_skin", iotime);

    //load_2d_field(fields.sp.at("thl")->fld_bot.data(), "thl_bot", iotime);
    //load_2d_field(fields.sp.at("qt")->fld_bot.data(), "qt_bot", iotime);

    //load_2d_field(fields.sp.at("thl")->flux_bot.data(), "thl_fluxbot", iotime);
    //load_2d_field(fields.sp.at("qt")->flux_bot.data(), "qt_fluxbot", iotime);

    for (auto& tile : tiles)
    {
        load_2d_field(tile.second.thl_bot.data(), "thl_bot_" + tile.first, iotime);
        load_2d_field(tile.second.qt_bot.data(),  "qt_bot_"  + tile.first, iotime);
    }

    // In case of heterogeneous land-surface, read spatial properties.
    //if (!sw_homogeneous)
    //{
    //    // 3D (soil) fields
    //    // BvS: yikes.. Read soil index as float, round/cast to int... TO-DO: fix this!
    //    load_3d_field(tmp3->fld.data(), "index_soil", 0);
    //    for (int i=0; i<sgd.ncells; ++i)
    //        soil_index[i] = std::round(tmp3->fld[i]);

    //    if (sw_water)
    //    {
    //        // More yikes.. Read water mask as float, cast to bool
    //        load_2d_field(tmp3->fld.data(), "water_mask", 0);
    //        for (int i=0; i<agd.ijcells; ++i)
    //            water_mask[i] = tmp3->fld[i] > TF(0.5) ? 1 : 0;
    //    }

    //    load_3d_field(root_fraction.data(), "root_frac", 0);

    //    // 2D (surface) fields
    //    load_2d_field(gD_coeff.data(), "gD", 0);
    //    load_2d_field(c_veg.data(), "c_veg", 0);
    //    load_2d_field(lai.data(), "lai", 0);
    //    load_2d_field(rs_veg_min.data(), "rs_veg_min", 0);
    //    load_2d_field(rs_soil_min.data(), "rs_soil_min", 0);
    //    load_2d_field(lambda_stable.data(), "lambda_stable", 0);
    //    load_2d_field(lambda_unstable.data(), "lambda_unstable", 0);
    //    load_2d_field(cs_veg.data(), "cs_veg", 0);
    //}

    // Check for any failures.
    master.sum(&nerror, 1);
    if (nerror)
        throw std::runtime_error("Error loading field(s)");

    fields.release_tmp(tmp1);
    fields.release_tmp(tmp2);
}

template<typename TF>
void Boundary_surface_lsm<TF>::save(const int iotime)
{
    auto& sgd = soil_grid.get_grid_data();

    auto tmp1 = fields.get_tmp();
    auto tmp2 = fields.get_tmp();

    int nerror = 0;
    const TF no_offset = TF(0);

    // Lambda function to save 2D fields.
    auto save_2d_field = [&](
            TF* const restrict field, const std::string& name)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", name.c_str(), iotime);
        master.print_message("Saving \"%s\" ... ", filename);

        const int kslice = 0;
        if (field3d_io.save_xy_slice(
                field, tmp1->fld.data(), filename, kslice))
        {
            master.print_message("FAILED\n");
            nerror += 1;
        }
        else
            master.print_message("OK\n");
    };

    // Lambda function to save the 3D soil fields.
    auto save_3d_field = [&](TF* const restrict field, const std::string& name)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", name.c_str(), iotime);
        master.print_message("Saving \"%s\" ... ", filename);

        if (field3d_io.save_field3d(
                field, tmp1->fld.data(), tmp2->fld.data(),
                filename, no_offset,
                sgd.kstart, sgd.kend))
        {
            master.print_message("FAILED\n");
            nerror += 1;
        }
        else
            master.print_message("OK\n");
    };

    // MO gradients are always needed, as the calculation of the
    // eddy viscosity use the gradients from the previous time step.
    save_2d_field(dudz_mo.data(), "dudz_mo");
    save_2d_field(dvdz_mo.data(), "dvdz_mo");
    save_2d_field(dbdz_mo.data(), "dbdz_mo");

    // Obukhov length restart files are only needed for the iterative solver.
    if (!sw_constant_z0)
        save_2d_field(obuk.data(), "obuk");

    // Don't save the initial soil temperature/moisture for heterogeneous runs.
    if (sw_homogeneous || iotime > 0)
    {
        save_3d_field(fields.sps.at("t")->fld.data(), "t_soil");
        save_3d_field(fields.sps.at("theta")->fld.data(), "theta_soil");
    }

    // Surface fields.
    save_2d_field(fields.ap2d.at("wl").data(), "wl_skin");

    //save_2d_field(fields.sp.at("thl")->fld_bot.data(), "thl_bot");
    //save_2d_field(fields.sp.at("qt")->fld_bot.data(), "qt_bot");

    //save_2d_field(fields.sp.at("thl")->flux_bot.data(), "thl_fluxbot");
    //save_2d_field(fields.sp.at("qt")->flux_bot.data(), "qt_fluxbot");

    for (auto& tile : tiles)
    {
        save_2d_field(tile.second.thl_bot.data(), "thl_bot_" + tile.first);
        save_2d_field(tile.second.qt_bot.data(),  "qt_bot_"  + tile.first);
    }

    // Check for any failures.
    master.sum(&nerror, 1);
    if (nerror)
        throw std::runtime_error("Error saving field(s)");

    fields.release_tmp(tmp1);
    fields.release_tmp(tmp2);
}

template<typename TF>
void Boundary_surface_lsm<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    auto& gd = grid.get_grid_data();
    //auto tmp1 = fields.get_tmp();
    //
    //for (auto& it : cross_list)
    //{
    //    if (it == "ustar")
    //        cross.cross_plane(ustar.data(), "ustar", iotime);
    //    else if (it == "obuk")
    //        cross.cross_plane(obuk.data(), "obuk", iotime);
    //    else if (it == "ra")
    //    {
    //        bsk::calc_ra(
    //                tmp1->flux_bot.data(), ustar.data(), obuk.data(),
    //                z0h.data(), gd.z[gd.kstart], gd.istart,
    //                gd.iend, gd.jstart, gd.jend, gd.icells);
    //        cross.cross_plane(tmp1->flux_bot.data(), "ra", iotime);
    //    }
    //}

    //fields.release_tmp(tmp1);
}

template<typename TF>
void Boundary_surface_lsm<TF>::exec_stats(Stats<TF>& stats)
{
    const TF no_offset = 0.;

    auto fld_mean = fields.get_tmp_xy();

    stats.calc_stats_2d("obuk", obuk, no_offset);
    stats.calc_stats_2d("ustar", ustar, no_offset);

    get_tiled_mean(*fld_mean, "H");
    stats.calc_stats_2d("H", *fld_mean, no_offset);

    get_tiled_mean(*fld_mean, "LE");
    stats.calc_stats_2d("LE", *fld_mean, no_offset);

    get_tiled_mean(*fld_mean, "G");
    stats.calc_stats_2d("G", *fld_mean, no_offset);

    stats.calc_stats_2d("rs_veg", tiles.at("veg").rs, no_offset);
    stats.calc_stats_2d("rs_soil", tiles.at("soil").rs, no_offset);

    if (sw_tile_stats)
        for (auto& tile : tiles)
        {
            stats.calc_stats_2d("c_"+tile.first, tile.second.fraction, no_offset);

            stats.calc_stats_2d("ustar_"+tile.first, tile.second.ustar, no_offset);
            stats.calc_stats_2d("obuk_"+tile.first, tile.second.obuk, no_offset);

            stats.calc_stats_2d("thl_bot_"+tile.first, tile.second.thl_bot, no_offset);
            stats.calc_stats_2d("qt_bot_"+tile.first, tile.second.qt_bot, no_offset);

            stats.calc_stats_2d("H_"+tile.first, tile.second.H, no_offset);
            stats.calc_stats_2d("LE_"+tile.first, tile.second.LE, no_offset);
            stats.calc_stats_2d("G_"+tile.first, tile.second.G, no_offset);
        }

    fields.release_tmp_xy(fld_mean);
}

#ifndef USECUDA
template<typename TF>
void Boundary_surface_lsm<TF>::exec_column(Column<TF>& column)
{
    const TF no_offset = 0.;
    //column.calc_time_series("obuk", obuk.data(), no_offset);
    //column.calc_time_series("ustar", ustar.data(), no_offset);
}
#endif

template<typename TF>
void Boundary_surface_lsm<TF>::set_values()
{
    auto& gd = grid.get_grid_data();

    // Call the base class function.
    Boundary<TF>::set_values();

    // Override the boundary settings in order to enforce dirichlet BC for surface model.
    bsk::set_bc<TF>(
            fields.mp.at("u")->fld_bot.data(),
            fields.mp.at("u")->grad_bot.data(),
            fields.mp.at("u")->flux_bot.data(),
            Boundary_type::Dirichlet_type, ubot,
            fields.visc, grid.utrans,
            gd.icells, gd.jcells);

    bsk::set_bc<TF>(
            fields.mp.at("v")->fld_bot.data(),
            fields.mp.at("v")->grad_bot.data(),
            fields.mp.at("v")->flux_bot.data(),
            Boundary_type::Dirichlet_type, vbot,
            fields.visc, grid.vtrans,
            gd.icells, gd.jcells);

    // Prepare the lookup table for the surface solver
    if (sw_constant_z0)
        init_solver();
}

// Prepare the surface layer solver.
template<typename TF>
void Boundary_surface_lsm<TF>::init_solver()
{
    auto& gd = grid.get_grid_data();

    zL_sl.resize(nzL_lut);
    f_sl.resize(nzL_lut);

    bsk::prepare_lut(
        zL_sl.data(),
        f_sl.data(),
        z0m[0], z0h[0],
        gd.z[gd.kstart], nzL_lut,
        mbcbot, thermobc);
}


template<typename TF>
void Boundary_surface_lsm<TF>::update_slave_bcs()
{
    // This function does nothing when the surface model is enabled, because
    // the fields are computed by the surface model in update_bcs.
}


template<typename TF>
void Boundary_surface_lsm<TF>::get_tiled_mean(std::vector<TF>& fld_out, std::string name)
{
    auto& gd = grid.get_grid_data();

    TF* fld_veg;
    TF* fld_soil;
    TF* fld_wet;

    if (name == "H")
    {
        fld_veg  = tiles.at("veg").H.data();
        fld_soil = tiles.at("soil").H.data();
        fld_wet  = tiles.at("wet").H.data();
    }
    else if (name == "LE")
    {
        fld_veg  = tiles.at("veg").LE.data();
        fld_soil = tiles.at("soil").LE.data();
        fld_wet  = tiles.at("wet").LE.data();
    }
    else if (name == "G")
    {
        fld_veg  = tiles.at("veg").G.data();
        fld_soil = tiles.at("soil").G.data();
        fld_wet  = tiles.at("wet").G.data();
    }
    else
        throw std::runtime_error("Cannot calculate tiled mean for variable \"" + name + "\"\\n");

    calc_tiled_mean(
            fld_out.data(),
            tiles.at("veg").fraction.data(),
            tiles.at("soil").fraction.data(),
            tiles.at("wet").fraction.data(),
            fld_veg,
            fld_soil,
            fld_wet,
            TF(1),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
}


template class Boundary_surface_lsm<double>;
template class Boundary_surface_lsm<float>;
