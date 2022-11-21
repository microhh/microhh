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

#include "master.h"
#include "input.h"
#include "grid.h"
#include "soil_grid.h"
#include "fields.h"
#include "diff.h"
#include "boundary.h"
#include "boundary_surface.h"
#include "boundary_surface_kernels.h"
#include "defines.h"
#include "constants.h"
#include "thermo.h"
#include "master.h"
#include "cross.h"
#include "column.h"
#include "monin_obukhov.h"
#include "fast_math.h"

namespace
{
    // Make a shortcut in the file scope.
    namespace most = Monin_obukhov;
    namespace fm = Fast_math;
    namespace bsk = Boundary_surface_kernels;

    template<typename TF, bool sw_constant_z0>
    void stability(
            TF* const restrict ustar,
            TF* const restrict obuk,
            const TF* const restrict bfluxbot,
            const TF* const restrict b,
            const TF* const restrict bbot,
            const TF* const restrict dutot,
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
        const int ii2 = 2;
        const int jj = icells;
        const int jj2 = 2*icells;

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
            const TF* const restrict dutot,
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

            // CvH: I think that the problem is not closed, since both the fluxes and the surface values
            // of u and v are unknown. You have to assume a no slip in order to get the fluxes and therefore
            // should not update the surface values with those that belong to the flux. This procedure needs
            // to be checked more carefully.
            /*
            // calculate the surface values
            for (int j=grid->jstart; j<grid->jend; ++j)
                #pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                ij  = i + j*jj;
                ijk = i + j*jj + kstart*kk;
                // interpolate the whole stability function rather than ustar or obuk
                ubot[ij] = 0.;// ufluxbot[ij] / (0.5*(ustar[ij-ii]*fm(zsl, z0m, obuk[ij-ii]) + ustar[ij]*fm(zsl, z0m, obuk[ij]))) + u[ijk];
                vbot[ij] = 0.;// vfluxbot[ij] / (0.5*(ustar[ij-jj]*fm(zsl, z0m, obuk[ij-jj]) + ustar[ij]*fm(zsl, z0m, obuk[ij]))) + v[ijk];
            }

            grid->boundary_cyclic_2d(ubot);
            grid->boundary_cyclic_2d(vbot);
            */
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

    template<typename TF>
    void calc_z0_charnock(
            TF* const restrict z0m,
            TF* const restrict z0h,
            const TF* const restrict ustar,
            const TF alpha_m, const TF alpha_ch, const TF alpha_h,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        const TF visc = TF(1.5e-5);
        const TF gi = TF(1)/Constants::grav<TF>;
        const TF min_ustar = TF(1e-8);

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;

                // Limit u* to prevent div/0:
                const TF ustar_lim = std::max(ustar[ij], min_ustar);

                // Roughness lengths, like IFS:
                z0m[ij] = alpha_m * visc/ustar_lim + alpha_ch * fm::pow2(ustar_lim) * gi;
                z0h[ij] = alpha_h * visc/ustar_lim;
            }
    }
}

template<typename TF>
Boundary_surface<TF>::Boundary_surface(
        Master& masterin, Grid<TF>& gridin, Soil_grid<TF>& soilgridin,
        Fields<TF>& fieldsin, Input& inputin) :
        Boundary<TF>(masterin, gridin, soilgridin, fieldsin, inputin)
{
    swboundary = "surface";

    #ifdef USECUDA
    ustar_g = 0;
    obuk_g  = 0;
    nobuk_g = 0;
    zL_sl_g = 0;
    f_sl_g  = 0;
    #endif
}

template<typename TF>
Boundary_surface<TF>::~Boundary_surface()
{
}

template<typename TF>
void Boundary_surface<TF>::create(
        Input& input, Netcdf_handle& input_nc,
        Stats<TF>& stats, Column<TF>& column,
        Cross<TF>& cross, Timeloop<TF>& timeloop)
{
    const std::string group_name = "default";
    Boundary<TF>::process_time_dependent(input, input_nc, timeloop);
    Boundary<TF>::process_inflow(input, input_nc);

    // add variables to the statistics
    if (stats.get_switch())
    {
        stats.add_time_series("ustar", "Surface friction velocity", "m s-1", group_name);
        stats.add_time_series("obuk", "Obukhov length", "m", group_name);

        if (sw_charnock)
        {
            stats.add_time_series("z0m", "Roughness length momentum", "m", group_name);
            stats.add_time_series("z0h", "Roughness length heat", "m", group_name);
        }
    }

    if (column.get_switch())
    {
        column.add_time_series("ustar", "Surface friction velocity", "m s-1");
        column.add_time_series("obuk", "Obukhov length", "m");

        if (sw_charnock)
        {
            column.add_time_series("z0m", "Roughness length momentum", "m");
            column.add_time_series("z0h", "Roughness length heat", "m");
        }
    }

    if (cross.get_switch())
    {
        std::vector<std::string> allowed_crossvars = {"ustar", "obuk", "ra"};

        if (sw_charnock)
        {
            allowed_crossvars.push_back("z0m");
            allowed_crossvars.push_back("z0h");
        }

        cross_list = cross.get_enabled_variables(allowed_crossvars);
    }
}

template<typename TF>
void Boundary_surface<TF>::create_cold_start(Netcdf_handle& input_nc)
{
}

template<typename TF>
void Boundary_surface<TF>::init(Input& inputin, Thermo<TF>& thermo)
{
    // 1. Process the boundary conditions now all fields are registered.
    process_bcs(inputin);

    // 2. Read and check the boundary_surface specific settings.
    process_input(inputin, thermo);

    // 3. Allocate and initialize the 2D surface fields.
    init_surface(inputin, thermo);

    // 4. Initialize the boundary cyclic.
    boundary_cyclic.init();
}

template<typename TF>
void Boundary_surface<TF>::process_input(Input& inputin, Thermo<TF>& thermo)
{
    // Switch between heterogeneous and homogeneous z0's
    sw_constant_z0 = inputin.get_item<bool>("boundary", "swconstantz0", "", true);

    #ifdef USECUDA
    if (!sw_constant_z0)
        throw std::runtime_error("\"boundary_surface\" with heterogeneous z0s is not (yet) supported");
    #endif

    // Switch for z0 as function of u* (Charnock relation)
    sw_charnock = inputin.get_item<bool>("boundary", "swcharnock", "", false);

    if (sw_charnock && sw_constant_z0)
        throw std::runtime_error("\"swcharnock=true\" requires \"swconstantz0=false\"");

    if (sw_charnock)
    {
        alpha_m  = inputin.get_item<TF>("boundary", "alpha_m", "");
        alpha_ch = inputin.get_item<TF>("boundary", "alpha_ch", "");
        alpha_h  = inputin.get_item<TF>("boundary", "alpha_h", "");
    }

    // crash in case fixed gradient is prescribed
    if (mbcbot == Boundary_type::Neumann_type)
    {
        std::string msg = "Neumann bc is not supported in surface model";
        throw std::runtime_error(msg);
    }

    // read the ustar value only if fixed fluxes are prescribed
    else if (mbcbot == Boundary_type::Ustar_type)
        ustarin = inputin.get_item<TF>("boundary", "ustar", "");

    // process the scalars
    for (auto& it : sbc)
    {
        // crash in case fixed gradient is prescribed
        if (it.second.bcbot == Boundary_type::Neumann_type)
        {
            std::string msg = "Fixed Gradient bc is not supported in surface model";
            throw std::runtime_error(msg);
        }

        // crash in case of fixed momentum flux and dirichlet bc for scalar
        if (it.second.bcbot == Boundary_type::Dirichlet_type && mbcbot == Boundary_type::Ustar_type)
        {
            std::string msg = "Fixed Ustar bc in combination with Dirichlet bc for scalars is not supported";
            throw std::runtime_error(msg);
        }
    }

    // check whether the prognostic thermo vars are of the same type
    std::vector<std::string> thermolist;
    thermo.get_prog_vars(thermolist);

    auto it = thermolist.begin();

    // save the bc of the first thermo field in case thermo is enabled
    if (it != thermolist.end())
        thermobc = sbc[*it].bcbot;
    else
        // Set the thermobc to Flux_type to avoid ininitialized errors.
        thermobc = Boundary_type::Flux_type;

    while (it != thermolist.end())
    {
        if (sbc[*it].bcbot != thermobc)
        {

            std::string msg = "All thermo variables need to have the same bc type";
            throw std::runtime_error(msg);
        }
        ++it;
    }
}

template<typename TF>
void Boundary_surface<TF>::init_surface(Input& input, Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();

    obuk.resize(gd.ijcells);
    ustar.resize(gd.ijcells);

    dudz_mo.resize(gd.ijcells);
    dvdz_mo.resize(gd.ijcells);
    if (thermo.get_switch() != "0")
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
    //std::fill(ustar.begin(), ustar.end(), Constants::dsmall);
    std::fill(ustar.begin(), ustar.end(), 1e-2);

    if (sw_constant_z0)
        std::fill(nobuk.begin(), nobuk.end(), 0);
}

template<typename TF>
void Boundary_surface<TF>::load(const int iotime, Thermo<TF>& thermo)
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
    if (thermo.get_switch() != "0")
        load_2d_field(dbdz_mo.data(), "dbdz_mo", iotime);

    // The `fld->gradbot`'s are only needed for flux BCs, required by
    // `set_ghost_cells()` at the start of the time loop.
    for (auto& it : fields.sp)
        if (sbc.at(it.first).bcbot == Boundary_type::Flux_type)
            load_2d_field(it.second->grad_bot.data(), it.first + "_gradbot", iotime);

    // Obukhov length restart files are only needed for the iterative solver
    if (!sw_constant_z0)
    {
        // Read Obukhov length
        load_2d_field(obuk.data(), "obuk", iotime);

        // Read spatial z0 fields
        if (!sw_charnock)
        {
            load_2d_field(z0m.data(), "z0m", 0);
            load_2d_field(z0h.data(), "z0h", 0);
        }
    }

    // Check for any failures.
    master.sum(&nerror, 1);
    if (nerror)
        throw std::runtime_error("Error loading field(s)");

    fields.release_tmp(tmp1);
}

template<typename TF>
void Boundary_surface<TF>::save(const int iotime, Thermo<TF>& thermo)
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
    if (thermo.get_switch() != "0")
        save_2d_field(dbdz_mo.data(), "dbdz_mo", iotime);

    // The `fld->gradbot`'s are only needed for flux BCs, required by
    // `set_ghost_cells()` at the start of the time loop.
    for (auto& it : fields.sp)
        if (sbc.at(it.first).bcbot == Boundary_type::Flux_type)
            save_2d_field(it.second->grad_bot.data(), it.first + "_gradbot", iotime);

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
void Boundary_surface<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    auto& gd = grid.get_grid_data();
    auto tmp1 = fields.get_tmp();

    for (auto& it : cross_list)
    {
        if (it == "ustar")
            cross.cross_plane(ustar.data(), "ustar", iotime);
        else if (it == "obuk")
            cross.cross_plane(obuk.data(), "obuk", iotime);
        else if (it == "z0m")
            cross.cross_plane(z0m.data(), "z0m", iotime);
        else if (it == "z0h")
            cross.cross_plane(z0h.data(), "z0h", iotime);
        else if (it == "ra")
        {
            bsk::calc_ra(
                    tmp1->flux_bot.data(), ustar.data(), obuk.data(),
                    z0h.data(), gd.z[gd.kstart], gd.istart,
                    gd.iend, gd.jstart, gd.jend, gd.icells);
            cross.cross_plane(tmp1->flux_bot.data(), "ra", iotime);
        }
    }

    fields.release_tmp(tmp1);
}

template<typename TF>
void Boundary_surface<TF>::exec_stats(Stats<TF>& stats)
{
    const TF no_offset = 0.;
    stats.calc_stats_2d("obuk", obuk, no_offset);
    stats.calc_stats_2d("ustar", ustar, no_offset);

    if (sw_charnock)
    {
        stats.calc_stats_2d("z0m", z0m, no_offset);
        stats.calc_stats_2d("z0h", z0h, no_offset);
    }
}

#ifndef USECUDA
template<typename TF>
void Boundary_surface<TF>::exec_column(Column<TF>& column)
{
    const TF no_offset = 0.;
    column.calc_time_series("obuk", obuk.data(), no_offset);
    column.calc_time_series("ustar", ustar.data(), no_offset);

    if (sw_charnock)
    {
        column.calc_time_series("z0m", z0m.data(), no_offset);
        column.calc_time_series("z0h", z0h.data(), no_offset);
    }
}
#endif

template<typename TF>
void Boundary_surface<TF>::set_values()
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

    // in case the momentum has a fixed ustar, set the value to that of the input
    if (mbcbot == Boundary_type::Ustar_type)
        set_ustar();

    // Prepare the lookup table for the surface solver
    if (sw_constant_z0)
        init_solver();
}

template<typename TF>
void Boundary_surface<TF>::set_ustar()
{
    auto& gd = grid.get_grid_data();
    const int jj = gd.icells;

    bsk::set_bc<TF>(
            fields.mp.at("u")->fld_bot.data(),
            fields.mp.at("u")->grad_bot.data(),
            fields.mp.at("u")->flux_bot.data(),
            mbcbot, ubot, fields.visc, grid.utrans,
            gd.icells, gd.jcells);

    bsk::set_bc<TF>(
            fields.mp.at("v")->fld_bot.data(),
            fields.mp.at("v")->grad_bot.data(),
            fields.mp.at("v")->flux_bot.data(),
            mbcbot, vbot, fields.visc, grid.vtrans,
            gd.icells, gd.jcells);

    for (int j=0; j<gd.jcells; ++j)
        #pragma ivdep
        for (int i=0; i<gd.icells; ++i)
        {
            const int ij = i + j*jj;
            // Limit ustar at 1e-4 to avoid zero divisions.
            ustar[ij] = std::max(static_cast<TF>(0.0001), ustarin);
        }
}

// Prepare the surface layer solver.
template<typename TF>
void Boundary_surface<TF>::init_solver()
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
void Boundary_surface<TF>::exec(
        Thermo<TF>& thermo, Radiation<TF>& radiation,
        Microphys<TF>& microphys, Timeloop<TF>& timeloop)
{
    auto& gd = grid.get_grid_data();

    // Update roughness lengths when Charnock relation is used,
    // using friction velocity from previous time step (?).
    if (sw_charnock)
    {
        calc_z0_charnock(
                z0m.data(), z0h.data(),
                ustar.data(),
                alpha_m, alpha_ch, alpha_h,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells);

        boundary_cyclic.exec_2d(z0m.data());
        boundary_cyclic.exec_2d(z0h.data());
    }

    // Calculate (limited and filtered) total wind speed difference surface-atmosphere:
    auto dutot = fields.get_tmp();

    bsk::calc_dutot(
            dutot->fld.data(),
            fields.mp.at("u")->fld.data(),
            fields.mp.at("v")->fld.data(),
            fields.mp.at("u")->fld_bot.data(),
            fields.mp.at("v")->fld_bot.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart,
            gd.icells, gd.jcells, gd.ijcells,
            boundary_cyclic);

    // Start with retrieving the stability information.
    if (thermo.get_switch() == "0")
    {
        stability_neutral(
                ustar.data(),
                obuk.data(),
                dutot->fld.data(),
                gd.z.data(), z0m.data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart,
                gd.icells, gd.jcells,
                gd.ijcells, mbcbot,
                boundary_cyclic);
    }
    else
    {
        auto buoy = fields.get_tmp();

        thermo.get_buoyancy_surf(buoy->fld, buoy->fld_bot, false);
        thermo.get_buoyancy_fluxbot(buoy->flux_bot, false);

        const TF db_ref = thermo.get_db_ref();

        if (sw_constant_z0)
            stability<TF, true>(
                    ustar.data(), obuk.data(),
                    buoy->flux_bot.data(),
                    buoy->fld.data(),
                    buoy->fld_bot.data(),
                    dutot->fld.data(),
                    gd.z.data(),
                    z0m.data(), z0h.data(),
                    zL_sl.data(), f_sl.data(), nobuk.data(), db_ref,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend, gd.kstart,
                    gd.icells, gd.jcells, gd.ijcells,
                    mbcbot, thermobc, boundary_cyclic);
        else
            stability<TF, false>(
                    ustar.data(), obuk.data(),
                    buoy->flux_bot.data(),
                    buoy->fld.data(),
                    buoy->fld_bot.data(),
                    dutot->fld.data(), gd.z.data(),
                    z0m.data(), z0h.data(),
                    zL_sl.data(), f_sl.data(), nobuk.data(), db_ref,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart,
                    gd.icells, gd.jcells, gd.ijcells,
                    mbcbot, thermobc, boundary_cyclic);

        fields.release_tmp(buoy);
    }

    fields.release_tmp(dutot);

    // Calculate the surface value, gradient and flux depending on the chosen boundary condition.
    // Momentum:
    surfm(fields.mp.at("u")->flux_bot.data(),
          fields.mp.at("v")->flux_bot.data(),
          fields.mp.at("u")->grad_bot.data(),
          fields.mp.at("v")->grad_bot.data(),
          ustar.data(), obuk.data(),
          fields.mp.at("u")->fld.data(),
          fields.mp.at("u")->fld_bot.data(),
          fields.mp.at("v")->fld.data(),
          fields.mp.at("v")->fld_bot.data(),
          z0m.data(), gd.z[gd.kstart], mbcbot,
          gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart,
          gd.icells, gd.jcells, gd.ijcells,
          boundary_cyclic);

    // Scalars:
    for (auto& it : fields.sp)
        surfs(it.second->fld_bot.data(),
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
    bsk::calc_duvdz_mo(
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

    if (thermo.get_switch() != "0")
    {
        auto buoy = fields.get_tmp();
        thermo.get_buoyancy_fluxbot(buoy->flux_bot, false);

        bsk::calc_dbdz_mo(
                dbdz_mo.data(), buoy->flux_bot.data(),
                ustar.data(), obuk.data(),
                gd.z[gd.kstart],
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells);

        fields.release_tmp(buoy);
    }
}
#endif

template<typename TF>
void Boundary_surface<TF>::update_slave_bcs()
{
    // This function does nothing when the surface model is enabled, because
    // the fields are computed by the surface model in update_bcs.
}

template class Boundary_surface<double>;
template class Boundary_surface<float>;
