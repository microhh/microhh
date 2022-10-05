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

#include <cmath>
#include "fast_math.h"
#include "master.h"
#include "input.h"
#include "grid.h"
#include "soil_grid.h"
#include "fields.h"
#include "thermo.h"
#include "boundary_surface_bulk.h"
#include "boundary_surface_kernels.h"
#include "monin_obukhov.h"
#include "constants.h"

namespace
{
    namespace fm = Fast_math;
    namespace most = Monin_obukhov;
    namespace bsk = Boundary_surface_kernels;

    template<typename TF>
    void momentum_fluxgrad(
            TF* restrict ufluxbot, TF* restrict vfluxbot,
            TF* restrict ugradbot, TF* restrict vgradbot,
            TF* restrict u, TF* restrict v,
            TF* restrict ubot, TF* restrict vbot,
            TF* restrict dutot, const TF Cm, const TF zsl,
            const int istart, const int iend, const int jstart, const int jend, const int kstart,
            const int jj, const int kk)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;

                ufluxbot[ij] = -Cm * dutot[ij] * (u[ijk]-ubot[ij]);
                vfluxbot[ij] = -Cm * dutot[ij] * (v[ijk]-vbot[ij]);
                ugradbot[ij] = (u[ijk]-ubot[ij])/zsl;
                vgradbot[ij] = (v[ijk]-vbot[ij])/zsl;
            }
    }

    template<typename TF>
    void scalar_fluxgrad(
            TF* restrict sfluxbot, TF* restrict sgradbot, TF* restrict s, TF* restrict sbot,
            TF* restrict dutot, const TF Cs, const TF zsl,
            const int istart, const int iend, const int jstart, const int jend, const int kstart,
            const int jj, const int kk)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + kstart*kk;

                sfluxbot[ij] = -Cs * dutot[ij] * (s[ijk]-sbot[ij]);
                sgradbot[ij] = (s[ijk]-sbot[ij])/zsl;
            }
    }

    template<typename TF>
    void surface_scaling(
            TF* restrict ustar, TF* restrict obuk, TF* restrict dutot, TF* restrict bfluxbot, const TF Cm,
            const int istart, const int iend, const int jstart, const int jend, const int jj)
    {
        const double sqrt_Cm = sqrt(Cm);

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*jj;

                ustar[ij] = sqrt_Cm * dutot[ij];
                obuk[ij] = -fm::pow3(ustar[ij]) / (Constants::kappa<TF> * bfluxbot[ij]);
            }
    }
}

template<typename TF>
Boundary_surface_bulk<TF>::Boundary_surface_bulk(
        Master& masterin, Grid<TF>& gridin, Soil_grid<TF>& soilgridin,
        Fields<TF>& fieldsin, Input& inputin) :
    Boundary<TF>(masterin, gridin, soilgridin, fieldsin, inputin)
{
    swboundary = "surface_bulk";

    #ifdef USECUDA
    ustar_g = 0;
    obuk_g  = 0;
    #endif
}

template<typename TF>
Boundary_surface_bulk<TF>::~Boundary_surface_bulk()
{
}

template<typename TF>
void Boundary_surface_bulk<TF>::create(
        Input& input, Netcdf_handle& input_nc,
        Stats<TF>& stats, Column<TF>& column,
        Cross<TF>& cross, Timeloop<TF>& timeloop)
{
    const std::string group_name = "default";

    Boundary<TF>::process_time_dependent(input, input_nc, timeloop);

    // add variables to the statistics
    if (stats.get_switch())
    {
        stats.add_time_series("ustar", "Surface friction velocity", "m s-1", group_name);
        stats.add_time_series("obuk", "Obukhov length", "m", group_name);
    }
}

template<typename TF>
void Boundary_surface_bulk<TF>::create_cold_start(Netcdf_handle& input_nc)
{
}

template<typename TF>
void Boundary_surface_bulk<TF>::init(Input& inputin, Thermo<TF>& thermo)
{
    // 1. Process the boundary conditions now all fields are registered.
    process_bcs(inputin);

    // 2. Read and check the boundary_surface specific settings.
    process_input(inputin, thermo);

    // 3. Allocate and initialize the 2D surface fields.
    init_surface(inputin);

    // 4. Initialize the boundary cyclic.
    boundary_cyclic.init();
}

template<typename TF>
void Boundary_surface_bulk<TF>::process_input(Input& inputin, Thermo<TF>& thermo)
{

    // crash in case fixed gradient is prescribed
    if (mbcbot != Boundary_type::Dirichlet_type)
    {
        std::string msg = "Only \"noslip\" is allowed as mbcbot with swboundary=\"bulk\"";
        throw std::runtime_error(msg);
    }

    bulk_cm = inputin.get_item<TF>("boundary", "bulk_cm", "");

    // process the scalars
    for (auto& it : sbc)
    {
        if (it.second.bcbot != Boundary_type::Dirichlet_type)
        {
            std::string msg = "Only \"noslip\" is allowed as mbcbot with swboundary=\"bulk\"";
            throw std::runtime_error(msg);
        }
        bulk_cs[it.first] = inputin.get_item<TF>("boundary", "bulk_cs", it.first);
    }
}

template<typename TF>
void Boundary_surface_bulk<TF>::init_surface(Input& input)
{
    auto& gd = grid.get_grid_data();

    obuk.resize(gd.ijcells);
    ustar.resize(gd.ijcells);

    dudz_mo.resize(gd.ijcells);
    dvdz_mo.resize(gd.ijcells);
    dbdz_mo.resize(gd.ijcells);

    z0m.resize(gd.ijcells);
    z0h.resize(gd.ijcells);

    // BvS TMP
    const TF z0m_hom = input.get_item<TF>("boundary", "z0m", "");
    const TF z0h_hom = input.get_item<TF>("boundary", "z0h", "");

    std::fill(z0m.begin(), z0m.end(), z0m_hom);
    std::fill(z0h.begin(), z0h.end(), z0h_hom);

    // Initialize the obukhov length on a small number.
    std::fill(obuk.begin(), obuk.end(), Constants::dsmall);
}

template<typename TF>
void Boundary_surface_bulk<TF>::load(const int iotime, Thermo<TF>& thermo)
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

    // Check for any failures.
    master.sum(&nerror, 1);
    if (nerror)
        throw std::runtime_error("Error loading field(s)");

    fields.release_tmp(tmp1);
}

template<typename TF>
void Boundary_surface_bulk<TF>::save(const int iotime, Thermo<TF>& thermo)
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

    // Check for any failures.
    master.sum(&nerror, 1);
    if (nerror)
        throw std::runtime_error("Error saving field(s)");

    fields.release_tmp(tmp1);
}

template<typename TF>
void Boundary_surface_bulk<TF>::exec_stats(Stats<TF>& stats)
{
    const TF no_offset = 0.;
    stats.calc_stats_2d("obuk", obuk, no_offset);
    stats.calc_stats_2d("ustar", ustar, no_offset);
}

template<typename TF>
void Boundary_surface_bulk<TF>::exec_column(Column<TF>& column)
{
}

template<typename TF>
void Boundary_surface_bulk<TF>::set_values()
{
    // Call the base class function.
    Boundary<TF>::set_values();
}

#ifndef USECUDA
template<typename TF>
void Boundary_surface_bulk<TF>::exec(
        Thermo<TF>& thermo, Radiation<TF>& radiation,
        Microphys<TF>& microphys, Timeloop<TF>& timeloop)
{
    auto& gd = grid.get_grid_data();
    const TF zsl = gd.z[gd.kstart];

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

    // Calculate surface momentum fluxes and gradients
    momentum_fluxgrad(
            fields.mp.at("u")->flux_bot.data(),fields.mp.at("v")->flux_bot.data(),
            fields.mp.at("u")->grad_bot.data(),fields.mp.at("v")->grad_bot.data(),
            fields.mp.at("u")->fld.data(),fields.mp.at("v")->fld.data(),
            fields.mp.at("u")->fld_bot.data(),fields.mp.at("v")->fld_bot.data(),
            dutot->fld.data(), bulk_cm, zsl,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.icells, gd.ijcells);

    boundary_cyclic.exec_2d(fields.mp.at("u")->flux_bot.data());
    boundary_cyclic.exec_2d(fields.mp.at("v")->flux_bot.data());
    boundary_cyclic.exec_2d(fields.mp.at("u")->grad_bot.data());
    boundary_cyclic.exec_2d(fields.mp.at("v")->grad_bot.data());

    // Calculate surface scalar fluxes and gradients
    for (auto& it : fields.sp)
    {
        scalar_fluxgrad(it.second->flux_bot.data(), it.second->grad_bot.data(),
                        it.second->fld.data(), it.second->fld_bot.data(),
                        dutot->fld.data(), bulk_cs.at(it.first), zsl,
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.icells, gd.ijcells);
        boundary_cyclic.exec_2d(it.second->flux_bot.data());
        boundary_cyclic.exec_2d(it.second->grad_bot.data());
    }

    auto b= fields.get_tmp();
    thermo.get_buoyancy_fluxbot(b->flux_bot, false);
    surface_scaling(
            ustar.data(), obuk.data(),
            dutot->fld.data(), b->flux_bot.data(),
            bulk_cm,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

    fields.release_tmp(b);
    fields.release_tmp(dutot);

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
#endif

template<typename TF>
void Boundary_surface_bulk<TF>::update_slave_bcs()
{
    // This function does nothing when the surface model is enabled, because
    // the fields are computed by the surface model in update_bcs.
}

template class Boundary_surface_bulk<double>;
template class Boundary_surface_bulk<float>;
