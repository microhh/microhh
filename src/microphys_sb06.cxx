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
#include <cmath>
#include <vector>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "diff.h"
#include "stats.h"
#include "cross.h"
#include "column.h"
#include "thermo.h"
#include "thermo_moist_functions.h"
#include "constants.h"

#include "microphys.h"
#include "microphys_sb06.h"

namespace
{
    using namespace Constants;
    namespace fm = Fast_math;

    // For simplicity, define constants here for now. These should probably move to the header.
    template<typename TF> constexpr TF ql_min   = 1.e-6;               // Min cloud liquid water for which calculations are performed
    template<typename TF> constexpr TF qr_min   = 1.e-15;              // Min rain liquid water for which calculations are performed

    template<typename TF>
    void autoconversion(
            TF* const restrict qrt,
            TF* const restrict nrt,
            TF* const restrict qtt,
            TF* const restrict thlt,
            const TF* const restrict qr,
            const TF* const restrict ql,
            const TF* const restrict rho,
            const TF* const restrict exner,
            const TF nc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const TF x_star = 2.6e-10;       // SB06, list of symbols, same as UCLA-LES

        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;

                    if (ql[ijk] > ql_min<TF>)
                    {
                        const TF au_tend = TF(1e-12);

                        qrt[ijk]  += au_tend;
                        nrt[ijk]  += au_tend * rho[k] / x_star;
                        qtt[ijk]  -= au_tend;
                        thlt[ijk] += Lv<TF> / (cp<TF> * exner[k]) * au_tend;
                    }
                }
    }
}

template<typename TF>
Microphys_sb06<TF>::Microphys_sb06(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
        Microphys<TF>(masterin, gridin, fieldsin, inputin)
{
    auto& gd = grid.get_grid_data();
    swmicrophys = Microphys_type::SB06;

    // Read microphysics switches and settings
    cflmax = inputin.get_item<TF>("micro", "cflmax", "", 1.2);
    Nc0 = inputin.get_item<TF>("micro", "Nc0", "");

    // Initialize the qr (rain water specific humidity) and nr (droplot number concentration) fields
    const std::string group_name = "thermo";
    fields.init_prognostic_field("qi", "Ice specific humidity", "kg kg-1", group_name, gd.sloc);
    fields.init_prognostic_field("qr", "Rain water specific humidity", "kg kg-1", group_name, gd.sloc);
    fields.init_prognostic_field("qs", "Snow specific humidity", "kg kg-1", group_name, gd.sloc);
    fields.init_prognostic_field("qg", "Graupel specific humidity", "kg kg-1", group_name, gd.sloc);
    fields.init_prognostic_field("qh", "Hail specific humidity", "kg kg-1", group_name, gd.sloc);

    fields.init_prognostic_field("ni", "Number density ice", "m-3", group_name, gd.sloc);
    fields.init_prognostic_field("nr", "Number density rain", "m-3", group_name, gd.sloc);
    fields.init_prognostic_field("ns", "Number density snow", "m-3", group_name, gd.sloc);
    fields.init_prognostic_field("ng", "Number density graupel", "m-3", group_name, gd.sloc);
    fields.init_prognostic_field("nh", "Number density hail", "m-3", group_name, gd.sloc);

    // Load the viscosity for both fields.
    fields.sp.at("qr")->visc = inputin.get_item<TF>("fields", "svisc", "qr");
    fields.sp.at("qg")->visc = inputin.get_item<TF>("fields", "svisc", "qg");
    fields.sp.at("qs")->visc = inputin.get_item<TF>("fields", "svisc", "qs");
    fields.sp.at("qh")->visc = inputin.get_item<TF>("fields", "svisc", "qh");

    fields.sp.at("nr")->visc = inputin.get_item<TF>("fields", "svisc", "nr");
    fields.sp.at("ng")->visc = inputin.get_item<TF>("fields", "svisc", "ng");
    fields.sp.at("ns")->visc = inputin.get_item<TF>("fields", "svisc", "ns");
    fields.sp.at("nh")->visc = inputin.get_item<TF>("fields", "svisc", "nh");
}

template<typename TF>
Microphys_sb06<TF>::~Microphys_sb06()
{
}

template<typename TF>
void Microphys_sb06<TF>::init()
{
    auto& gd = grid.get_grid_data();

    rr_bot.resize(gd.ijcells);
    rs_bot.resize(gd.ijcells);
    rg_bot.resize(gd.ijcells);
}

template<typename TF>
void Microphys_sb06<TF>::create(
        Input& inputin, Netcdf_handle& input_nc,
        Stats<TF>& stats, Cross<TF>& cross, Dump<TF>& dump, Column<TF>& column)
{
    const std::string group_name = "thermo";

    // Add variables to the statistics
    if (stats.get_switch())
    {
        // Time series
        stats.add_time_series("rr", "Mean surface rain rate", "kg m-2 s-1", group_name);
        stats.add_time_series("rs", "Mean surface snow rate", "kg m-2 s-1", group_name);
        stats.add_time_series("rg", "Mean surface graupel rate", "kg m-2 s-1", group_name);

        stats.add_tendency(*fields.st.at("thl"), "z", tend_name, tend_longname);
        stats.add_tendency(*fields.st.at("qt") , "z", tend_name, tend_longname);
        stats.add_tendency(*fields.st.at("qr") , "z", tend_name, tend_longname);
        stats.add_tendency(*fields.st.at("qs") , "z", tend_name, tend_longname);
        stats.add_tendency(*fields.st.at("qg") , "z", tend_name, tend_longname);
    }

    if (column.get_switch())
    {
        column.add_time_series("rr", "Surface rain rate", "kg m-2 s-1");
        column.add_time_series("rs", "Surface snow rate", "kg m-2 s-1");
        column.add_time_series("rg", "Surface graupel rate", "kg m-2 s-1");
    }

    // Create cross sections
    // 1. Variables that this class can calculate/provide:
    const std::vector<std::string> allowed_crossvars = {"rr_bot", "rs_bot", "rg_bot"};

    // 2. Cross-reference with the variables requested in the .ini file:
    crosslist = cross.get_enabled_variables(allowed_crossvars);
}

#ifndef USECUDA
template<typename TF>
void Microphys_sb06<TF>::exec(Thermo<TF>& thermo, const double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    // Get liquid water, ice and pressure variables before starting.
    auto ql = fields.get_tmp();

    thermo.get_thermo_field(*ql, "ql", false, false);

    const std::vector<TF>& p = thermo.get_basestate_vector("p");
    const std::vector<TF>& exner = thermo.get_basestate_vector("exner");

    for (int k=gd.kend-1; k<=gd.kstart; --k)
    {
        // Sedimentation
        // TODO

        // Autoconversion; formation of rain drop by coagulating cloud droplets.
        autoconversion(
            fields.st.at("qr")->fld.data(),
            fields.st.at("nr")->fld.data(),
            fields.st.at("qt")->fld.data(),
            fields.st.at("thl")->fld.data(),
            fields.sp.at("qr")->fld.data(),
            ql->fld.data(),
            fields.rhoref.data(),
            exner.data(), Nc0,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            k, k+1,
            gd.icells, gd.ijcells);

        // Sedimentation
        // TODO
    }

    // Release temporary fields.
    fields.release_tmp(ql);

    // Calculate tendencies.
    stats.calc_tend(*fields.st.at("thl"), tend_name);
    stats.calc_tend(*fields.st.at("qt" ), tend_name);
    stats.calc_tend(*fields.st.at("qi" ), tend_name);
    stats.calc_tend(*fields.st.at("qr" ), tend_name);
    stats.calc_tend(*fields.st.at("qs" ), tend_name);
    stats.calc_tend(*fields.st.at("qg" ), tend_name);
    stats.calc_tend(*fields.st.at("qh" ), tend_name);
}
#endif

template<typename TF>
void Microphys_sb06<TF>::exec_stats(Stats<TF>& stats, Thermo<TF>& thermo, const double dt)
{
    // Time series
    const TF no_offset = 0.;
    stats.calc_stats_2d("rr", rr_bot, no_offset);
    stats.calc_stats_2d("rs", rs_bot, no_offset);
    stats.calc_stats_2d("rg", rg_bot, no_offset);
}

#ifndef USECUDA
template<typename TF>
void Microphys_sb06<TF>::exec_column(Column<TF>& column)
{
    const TF no_offset = 0.;
    column.calc_time_series("rr", rr_bot.data(), no_offset);
    column.calc_time_series("rs", rs_bot.data(), no_offset);
    column.calc_time_series("rg", rg_bot.data(), no_offset);
}
#endif

template<typename TF>
void Microphys_sb06<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    if (cross.get_switch())
    {
        for (auto& it : crosslist)
        {
            if (it == "rr_bot")
                cross.cross_plane(rr_bot.data(), "rr_bot", iotime);

            if (it == "rs_bot")
                cross.cross_plane(rs_bot.data(), "rs_bot", iotime);

            if (it == "rg_bot")
                cross.cross_plane(rg_bot.data(), "rg_bot", iotime);
        }
    }
}

#ifndef USECUDA
template<typename TF>
unsigned long Microphys_sb06<TF>::get_time_limit(unsigned long idt, const double dt)
{
    auto& gd = grid.get_grid_data();
    auto tmp = fields.get_tmp();

    double cfl = 0.;

    // TO-DO

    // Get maximum CFL across all MPI tasks
    master.max(&cfl, 1);
    fields.release_tmp(tmp);

    // Prevent zero division.
    cfl = std::max(cfl, 1.e-5);

    return idt * this->cflmax / cfl;
}
#endif

template<typename TF>
bool Microphys_sb06<TF>::has_mask(std::string name)
{
    return false;
}

template<typename TF>
void Microphys_sb06<TF>::get_mask(Stats<TF>& stats, std::string mask_name)
{
    std::string message = "SB06 microphysics scheme can not provide mask: \"" + mask_name +"\"";
    throw std::runtime_error(message);
}

template<typename TF>
void Microphys_sb06<TF>::get_surface_rain_rate(std::vector<TF>& field)
{
    // Make a hard copy of the surface rain precipitation field
    field = rr_bot;

    // Add snow and graupel surface precipitation
    std::transform(field.begin(), field.end(), rs_bot.begin(), field.begin(), std::plus<TF>());
    std::transform(field.begin(), field.end(), rg_bot.begin(), field.begin(), std::plus<TF>());
}

template class Microphys_sb06<double>;
template class Microphys_sb06<float>;
