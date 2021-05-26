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

namespace bs
{
    template<typename TF, bool sw_constant_z0>
    void stability(
            TF* const restrict ustar,
            TF* const restrict obuk,
            const TF* const restrict bfluxbot,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict b,
            const TF* const restrict ubot,
            const TF* const restrict vbot,
            const TF* const restrict bbot,
            TF* const restrict dutot,
            const TF* const restrict z,
            const TF* const restrict z0m,
            const TF* const restrict z0h,
            const float* const zL_sl,
            const float* const f_sl,
            int* const nobuk,
            const TF db_ref,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int jcells, const int kk,
            Boundary_type mbcbot, Boundary_type thermobc,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int ii = 1;
        const int jj = icells;

        // Calculate total wind.
        const TF minval = 1.e-1;

        // First, interpolate the wind to the scalar location.
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;
                const TF du2 = fm::pow2(TF(0.5)*(u[ijk] + u[ijk+ii]) - TF(0.5)*(ubot[ij] + ubot[ij+ii]))
                             + fm::pow2(TF(0.5)*(v[ijk] + v[ijk+jj]) - TF(0.5)*(vbot[ij] + vbot[ij+jj]));
                // Prevent the absolute wind gradient from reaching values less than 0.01 m/s,
                // otherwise evisc at k = kstart blows up
                dutot[ij] = std::max(std::pow(du2, TF(0.5)), minval);
            }

        boundary_cyclic.exec_2d(dutot);

        // Calculate Obukhov length
        // Case 1: fixed buoyancy flux and fixed ustar
        if (mbcbot == Boundary_type::Ustar_type && thermobc == Boundary_type::Flux_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    obuk[ij] = -fm::pow3(ustar[ij]) / (Constants::kappa<TF>*bfluxbot[ij]);
                }
        }
        // Case 2: fixed buoyancy surface value and free ustar
        else if (mbcbot == Boundary_type::Dirichlet_type && thermobc == Boundary_type::Flux_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;

                    // Switch between the iterative and lookup solver
                    if (sw_constant_z0)
                        obuk[ij] = bsk::calc_obuk_noslip_flux_lookup(
                                zL_sl, f_sl, nobuk[ij], dutot[ij], bfluxbot[ij], z[kstart]);
                    else
                        obuk[ij] = bsk::calc_obuk_noslip_flux_iterative(
                                obuk[ij], dutot[ij], bfluxbot[ij], z[kstart], z0m[ij]);

                    ustar[ij] = dutot[ij] * most::fm(z[kstart], z0m[ij], obuk[ij]);
                }
        }
        else if (mbcbot == Boundary_type::Dirichlet_type && thermobc == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    const TF db = b[ijk] - bbot[ij] + db_ref;

                    // Switch between the iterative and lookup solver
                    if (sw_constant_z0)
                        obuk[ij] = bsk::calc_obuk_noslip_dirichlet_lookup(
                                zL_sl, f_sl, nobuk[ij], dutot[ij], db, z[kstart]);
                    else
                        obuk[ij] = bsk::calc_obuk_noslip_dirichlet_iterative(
                                obuk[ij], dutot[ij], db, z[kstart], z0m[ij], z0h[ij]);

                    ustar[ij] = dutot[ij] * most::fm(z[kstart], z0m[ij], obuk[ij]);
                }
        }
    }

    template<typename TF>
    void stability_neutral(
            TF* const restrict ustar,
            TF* const restrict obuk,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict ubot,
            const TF* const restrict vbot,
            TF* const restrict dutot,
            const TF* const restrict z,
            const TF* const restrict z0m,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int jcells, const int kk,
            Boundary_type mbcbot,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int ii = 1;
        const int jj = icells;

        // calculate total wind
        TF du2;
        const TF minval = 1.e-1;

        // first, interpolate the wind to the scalar location
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;
                du2 = fm::pow2(TF(0.5)*(u[ijk] + u[ijk+ii]) - TF(0.5)*(ubot[ij] + ubot[ij+ii]))
                    + fm::pow2(TF(0.5)*(v[ijk] + v[ijk+jj]) - TF(0.5)*(vbot[ij] + vbot[ij+jj]));
                // prevent the absolute wind gradient from reaching values less than 0.01 m/s,
                // otherwise evisc at k = kstart blows up
                dutot[ij] = std::max(std::pow(du2, TF(0.5)), minval);
            }

        boundary_cyclic.exec_2d(dutot);

        // set the Obukhov length to a very large negative number
        // case 1: fixed buoyancy flux and fixed ustar
        if (mbcbot == Boundary_type::Ustar_type)
        {
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;
                    obuk[ij] = -Constants::dbig;
                }
        }
        // case 2: free ustar
        else if (mbcbot == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    obuk [ij] = -Constants::dbig;
                    ustar[ij] = dutot[ij] * most::fm(z[kstart], z0m[ij], obuk[ij]);
                }
        }
    }

    template<typename TF>
    void surfm(
            TF* const restrict ufluxbot,
            TF* const restrict vfluxbot,
            TF* const restrict ugradbot,
            TF* const restrict vgradbot,
            const TF* const restrict ustar,
            const TF* const restrict obuk,
            const TF* const restrict u, const TF* const restrict ubot,
            const TF* const restrict v, const TF* const restrict vbot,
            const TF* const restrict z0m,
            const TF zsl, const Boundary_type bcbot,
            const int istart, const int iend,
            const int jstart, const int jend, const int kstart,
            const int icells, const int jcells, const int kk,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int ii = 1;
        const int jj = icells;

        // the surface value is known, calculate the flux and gradient
        if (bcbot == Boundary_type::Dirichlet_type)
        {
            // first calculate the surface value
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;

                    // interpolate the whole stability function rather than ustar or obuk
                    ufluxbot[ij] = -(u[ijk]-ubot[ij])*TF(0.5)*
                        (ustar[ij-ii]*most::fm(zsl, z0m[ij-ii], obuk[ij-ii]) + ustar[ij]*most::fm(zsl, z0m[ij], obuk[ij]));
                    vfluxbot[ij] = -(v[ijk]-vbot[ij])*TF(0.5)*
                        (ustar[ij-jj]*most::fm(zsl, z0m[ij-jj], obuk[ij-jj]) + ustar[ij]*most::fm(zsl, z0m[ij], obuk[ij]));
                }

            boundary_cyclic.exec_2d(ufluxbot);
            boundary_cyclic.exec_2d(vfluxbot);
        }

        // the flux is known, calculate the surface value and gradient
        else if (bcbot == Boundary_type::Ustar_type)
        {
            // first redistribute ustar over the two flux components
            const TF minval = 1.e-2;

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;

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
        }

        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;
                // use the linearly interpolated grad, rather than the MO grad,
                // to prevent giving unresolvable gradients to advection schemes
                // vargradbot[ij] = -varfluxbot[ij] / (kappa*z0m*ustar[ij]) * phih(zsl/obuk[ij]);
                ugradbot[ij] = (u[ijk]-ubot[ij])/zsl;
                vgradbot[ij] = (v[ijk]-vbot[ij])/zsl;
            }
    }

    template<typename TF>
    void surfs(
            TF* const restrict varbot,
            TF* const restrict vargradbot,
            TF* const restrict varfluxbot,
            const TF* const restrict ustar,
            const TF* const restrict obuk,
            const TF* const restrict var,
            const TF* const restrict z0h,
            const TF zsl, const Boundary_type bcbot,
            const int istart, const int iend,
            const int jstart, const int jend, const int kstart,
            const int icells, const int jcells, const int kk,
            Boundary_cyclic<TF>& boundary_cyclic)
    {
        const int jj = icells;

        // the surface value is known, calculate the flux and gradient
        if (bcbot == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    varfluxbot[ij] = -(var[ijk]-varbot[ij])*ustar[ij]*most::fh(zsl, z0h[ij], obuk[ij]);
                    // vargradbot[ij] = -varfluxbot[ij] / (kappa*z0h*ustar[ij]) * phih(zsl/obuk[ij]);
                    // use the linearly interpolated grad, rather than the MO grad,
                    // to prevent giving unresolvable gradients to advection schemes
                    vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
                }
        }
        else if (bcbot == Boundary_type::Flux_type)
        {
            // the flux is known, calculate the surface value and gradient
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    varbot[ij] = varfluxbot[ij] / (ustar[ij]*most::fh(zsl, z0h[ij], obuk[ij])) + var[ijk];
                    // vargradbot[ij] = -varfluxbot[ij] / (kappa*z0h*ustar[ij]) * phih(zsl/obuk[ij]);
                    // use the linearly interpolated grad, rather than the MO grad,
                    // to prevent giving unresolvable gradients to advection schemes
                    vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
                }
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
    #ifdef USECUDA
    clear_device();
    #endif
}

template<typename TF>
void Boundary_surface_lsm<TF>::create(
        Input& input, Netcdf_handle& input_nc,
        Stats<TF>& stats, Column<TF>& column, Cross<TF>& cross)
{
    const std::string group_name = "default";

    //// add variables to the statistics
    //if (stats.get_switch())
    //{
    //    stats.add_time_series("ustar", "Surface friction velocity", "m s-1", group_name);
    //    stats.add_time_series("obuk", "Obukhov length", "m", group_name);
    //}

    //if (column.get_switch())
    //{
    //    column.add_time_series("ustar", "Surface friction velocity", "m s-1");
    //    column.add_time_series("obuk", "Obukhov length", "m");
    //}

    //if (cross.get_switch())
    //{
    //    const std::vector<std::string> allowed_crossvars = {"ustar", "obuk", "ra"};
    //    cross_list = cross.get_enabled_variables(allowed_crossvars);
    //}
}

template<typename TF>
void Boundary_surface_lsm<TF>::init(Input& inputin, Thermo<TF>& thermo)
{
    // Surface-layer
    // Process the boundary conditions now all fields are registered.
    process_bcs(inputin);

    // Read and check the boundary_surface specific settings.
    process_input(inputin, thermo);

    // Allocate and initialize the 2D surface fields.
    init_surface(inputin);
    init_lsm();

    // Initialize the boundary cyclic.
    boundary_cyclic.init();
}

template<typename TF>
void Boundary_surface_lsm<TF>::process_input(Input& inputin, Thermo<TF>& thermo)
{
    // check whether the prognostic thermo vars are of the same type
    std::vector<std::string> thermolist;
    thermo.get_prog_vars(thermolist);

    // save the bc of the first thermo field in case thermo is enabled
    auto it = thermolist.begin();
    if (it != thermolist.end())
        thermobc = sbc[*it].bcbot;
    else
        // Set the thermobc to Flux_type to avoid uninitialized errors.
        thermobc = Boundary_type::Flux_type;
}

template<typename TF>
void Boundary_surface_lsm<TF>::init_surface(Input& input)
{
    auto& gd = grid.get_grid_data();

    obuk.resize(gd.ijcells);
    ustar.resize(gd.ijcells);

    dudz_mo.resize(gd.ijcells);
    dvdz_mo.resize(gd.ijcells);
    dbdz_mo.resize(gd.ijcells);

    if (sw_constant_z0)
        nobuk.resize(gd.ijcells);

    z0m.resize(gd.ijcells);
    z0h.resize(gd.ijcells);

    if (sw_constant_z0)
    {
        const TF z0m_hom = input.get_item<TF>("boundary", "z0m", "");
        const TF z0h_hom = input.get_item<TF>("boundary", "z0h", "");

        std::fill(z0m.begin(), z0m.end(), z0m_hom);
        std::fill(z0h.begin(), z0h.end(), z0h_hom);
    }

    // Initialize the obukhov length on a small number.
    std::fill(obuk.begin(), obuk.end(), Constants::dsmall);

    // Also initialise ustar at small number, to prevent div/0
    // in calculation surface gradients during cold start.
    std::fill(ustar.begin(), ustar.end(), Constants::dsmall);

    if (sw_constant_z0)
        std::fill(nobuk.begin(), nobuk.end(), 0);
}

template<typename TF>
void Boundary_surface_lsm<TF>::init_lsm()
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

    diffusivity.resize   (sgd.ncells);
    diffusivity_h.resize (sgd.ncellsh);
    conductivity.resize  (sgd.ncells);
    conductivity_h.resize(sgd.ncellsh);
    source.resize        (sgd.ncells);

    root_fraction.resize(sgd.ncells);

    // Resize the lookup table
    lookup_table_size = nc_lookup_table->get_dimension_size("index");

    theta_res.resize(lookup_table_size);
    theta_wp.resize(lookup_table_size);
    theta_fc.resize(lookup_table_size);
    theta_sat.resize(lookup_table_size);

    gamma_theta_sat.resize(lookup_table_size);
    vg_a.resize(lookup_table_size);
    vg_l.resize(lookup_table_size);
    vg_n.resize(lookup_table_size);

    vg_m.resize(lookup_table_size);
    kappa_theta_max.resize(lookup_table_size);
    kappa_theta_min.resize(lookup_table_size);
    gamma_theta_max.resize(lookup_table_size);
    gamma_theta_min.resize(lookup_table_size);

    gamma_T_dry.resize(lookup_table_size);
    rho_C.resize(lookup_table_size);
}

template<typename TF>
void Boundary_surface_lsm<TF>::load(const int iotime)
{
    auto tmp1 = fields.get_tmp();
    int nerror = 0;

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

    // MO gradients are always needed, as the calculation of the
    // eddy viscosity use the gradients from the previous time step.
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

    // Check for any failures.
    master.sum(&nerror, 1);
    if (nerror)
        throw std::runtime_error("Error loading field(s)");

    fields.release_tmp(tmp1);
}

template<typename TF>
void Boundary_surface_lsm<TF>::save(const int iotime)
{
    auto tmp1 = fields.get_tmp();
    int nerror = 0;

    auto save_2d_field = [&](
            TF* const restrict field, const std::string& name, const int itime)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", name.c_str(), itime);
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

    // MO gradients are always needed, as the calculation of the
    // eddy viscosity use the gradients from the previous time step.
    save_2d_field(dudz_mo.data(), "dudz_mo", iotime);
    save_2d_field(dvdz_mo.data(), "dvdz_mo", iotime);
    save_2d_field(dbdz_mo.data(), "dbdz_mo", iotime);

    // Obukhov length restart files are only needed for the iterative solver
    if (!sw_constant_z0)
        save_2d_field(obuk.data(), "obuk", iotime);

    // Check for any failures.
    master.sum(&nerror, 1);
    if (nerror)
        throw std::runtime_error("Error saving field(s)");

    fields.release_tmp(tmp1);
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
    //stats.calc_stats_2d("obuk", obuk, no_offset);
    //stats.calc_stats_2d("ustar", ustar, no_offset);
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

#ifndef USECUDA
template<typename TF>
void Boundary_surface_lsm<TF>::exec(Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();

    // Start with retrieving the stability information.
    if (thermo.get_switch() == "0")
    {
        auto dutot = fields.get_tmp();
        bs::stability_neutral(
                ustar.data(), obuk.data(),
                fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(),
                fields.mp.at("u")->fld_bot.data(), fields.mp.at("v")->fld_bot.data(),
                dutot->fld.data(), gd.z.data(), z0m.data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend, gd.kstart,
                gd.icells, gd.jcells, gd.ijcells,
                mbcbot, boundary_cyclic);
        fields.release_tmp(dutot);
    }
    else
    {
        auto buoy = fields.get_tmp();
        auto tmp = fields.get_tmp();

        thermo.get_buoyancy_surf(*buoy, false);
        const TF db_ref = thermo.get_db_ref();

        if (sw_constant_z0)
            bs::stability<TF, true>(
                    ustar.data(), obuk.data(),
                    buoy->flux_bot.data(),
                    fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(),
                    buoy->fld.data(), fields.mp.at("u")->fld_bot.data(),
                    fields.mp.at("v")->fld_bot.data(), buoy->fld_bot.data(),
                    tmp->fld.data(), gd.z.data(),
                    z0m.data(), z0h.data(),
                    zL_sl.data(), f_sl.data(), nobuk.data(), db_ref,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend, gd.kstart,
                    gd.icells, gd.jcells, gd.ijcells,
                    mbcbot, thermobc, boundary_cyclic);
        else
            bs::stability<TF, false>(
                    ustar.data(), obuk.data(),
                    buoy->flux_bot.data(),
                    fields.mp.at("u")->fld.data(), fields.mp.at("v")->fld.data(), buoy->fld.data(),
                    fields.mp.at("u")->fld_bot.data(), fields.mp.at("v")->fld_bot.data(), buoy->fld_bot.data(),
                    tmp->fld.data(), gd.z.data(),
                    z0m.data(), z0h.data(),
                    zL_sl.data(), f_sl.data(), nobuk.data(), db_ref,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart,
                    gd.icells, gd.jcells, gd.ijcells,
                    mbcbot, thermobc, boundary_cyclic);

        fields.release_tmp(buoy);
        fields.release_tmp(tmp);
    }

    // Calculate the surface value, gradient and flux depending on the chosen boundary condition.
    // Momentum:
    bs::surfm(fields.mp.at("u")->flux_bot.data(),
          fields.mp.at("v")->flux_bot.data(),
          fields.mp.at("u")->grad_bot.data(),
          fields.mp.at("v")->grad_bot.data(),
          ustar.data(), obuk.data(),
          fields.mp.at("u")->fld.data(), fields.mp.at("u")->fld_bot.data(),
          fields.mp.at("v")->fld.data(), fields.mp.at("v")->fld_bot.data(),
          z0m.data(), gd.z[gd.kstart], mbcbot,
          gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart,
          gd.icells, gd.jcells, gd.ijcells,
          boundary_cyclic);

    // Calculate MO gradients
    bsk::calc_duvdz(
            dudz_mo.data(), dvdz_mo.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("u")->fld_bot.data(),
            fields.mp.at("v")->fld_bot.data(),
            fields.mp.at("u")->flux_bot.data(),
            fields.mp.at("v")->flux_bot.data(),
            ustar.data(), obuk.data(), z0m.data(),
            gd.z[gd.kstart],
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart,
            gd.icells, gd.ijcells);

    // Scalars:
    for (auto& it : fields.sp)
        bs::surfs(it.second->fld_bot.data(),
              it.second->grad_bot.data(),
              it.second->flux_bot.data(),
              ustar.data(), obuk.data(),
              it.second->fld.data(), z0h.data(),
              gd.z[gd.kstart], sbc.at(it.first).bcbot,
              gd.istart, gd.iend,
              gd.jstart, gd.jend, gd.kstart,
              gd.icells, gd.jcells, gd.ijcells,
              boundary_cyclic);

    // Calculate MO gradients
    auto buoy = fields.get_tmp();
    thermo.get_buoyancy_fluxbot(*buoy, false);

    bsk::calc_dbdz(
            dbdz_mo.data(), buoy->flux_bot.data(),
            ustar.data(), obuk.data(),
            gd.z[gd.kstart],
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

    fields.release_tmp(buoy);
}
#endif

template<typename TF>
void Boundary_surface_lsm<TF>::update_slave_bcs()
{
    // This function does nothing when the surface model is enabled, because
    // the fields are computed by the surface model in update_bcs.
}

template class Boundary_surface_lsm<double>;
template class Boundary_surface_lsm<float>;
