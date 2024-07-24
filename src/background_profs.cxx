/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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
#include <algorithm>

#include "background_profs.h"

#include "timeloop.h"
#include "input.h"
#include "grid.h"
#include "netcdf_interface.h"
#include "stats.h"
#include "thermo.h"
#include "fields.h"
#include "timedep.h"
#include "Array.h"
#include "Gas_concs.h"

using Aerosol_concs = Gas_concs;

template<typename TF>
Background<TF>::Background(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
        master(masterin), grid(gridin), fields(fieldsin)
{
    // Read `.ini` settings.
    swtimedep_background = inputin.get_item<bool>("radiation", "swtimedep_background", "", false);

    if (!swtimedep_background)
        return;

    sw_aerosol = inputin.get_item<bool>("aerosol", "swaerosol", "", false);
    swtimedep_aerosol = inputin.get_item<bool>("aerosol", "swtimedep", "", false);
    dt_rad = inputin.get_item<double>("radiation", "dt_rad", "");
    gaslist = inputin.get_list<std::string>("radiation", "timedeplist_gas", "", std::vector<std::string>());

    const std::vector<std::string> possible_gases = {
            "h2o", "co2" ,"o3", "n2o", "co", "ch4", "o2", "n2",
            "ccl4", "cfc11", "cfc12", "cfc22",
            "hfc143a", "hfc125", "hfc23", "hfc32", "hfc134a",
            "cf4", "no2" };

    for (auto& it : gaslist)
    {
        if (std::find(possible_gases.begin(), possible_gases.end(), it) != possible_gases.end())
        {
            tdep_gases.emplace(it, new Timedep<TF>(master, grid, it+"_bg", swtimedep_background));
        }
        else
        {
            std::cout << "Unsupported gas \"" + it+"_bg" + "\" in timedeplist_gas" << std::endl;
        }
    }
}

template <typename TF>
Background<TF>::~Background()
{
}

template <typename TF>
void Background<TF>::init(Netcdf_handle& input_nc)
{
    // Always get dimensions background levels, if radiation group is present.
    if (input_nc.group_exists("radiation"))
    {
        n_lay = -1;
        n_lev = -1;

        Netcdf_handle& rad_nc = input_nc.get_group("radiation");

        if (rad_nc.dimension_exists("lay"))
            n_lay = rad_nc.get_dimension_size("lay");
        if (rad_nc.dimension_exists("lev"))
            n_lev = rad_nc.get_dimension_size("lev");
    }

    if (!swtimedep_background)
        return;

    idt_rad = convert_to_itime(dt_rad);

    // temperature, pressure and moisture
    t_lay.resize(n_lay);
    t_lev.resize(n_lev);
    p_lay.resize(n_lay);
    p_lev.resize(n_lev);
    h2o.resize(n_lay);

    for (auto& it : gaslist)
        gasprofs[it] = std::vector<TF>(n_lay);

    // aerosols
    aermr01.resize(n_lay);
    aermr02.resize(n_lay);
    aermr03.resize(n_lay);
    aermr04.resize(n_lay);
    aermr05.resize(n_lay);
    aermr06.resize(n_lay);
    aermr07.resize(n_lay);
    aermr08.resize(n_lay);
    aermr09.resize(n_lay);
    aermr10.resize(n_lay);
    aermr11.resize(n_lay);
}

template <typename TF>
void Background<TF>::create(Input& inputin, Netcdf_handle& input_nc, Stats<TF>& stats)
{
    // Read input from NetCDF and prepare statistics output.
    if (!swtimedep_background)
        return;

    // create time dependent profiles
    const TF offset = 0;
    std::string timedep_dim = "time_rad";

    // temperature, pressure and moisture
    tdep_t_lay = std::make_unique<Timedep<TF>>(master, grid, "t_lay", swtimedep_background);
    tdep_t_lay->create_timedep_prof(input_nc, offset, timedep_dim, n_lay);
    tdep_t_lev = std::make_unique<Timedep<TF>>(master, grid, "t_lev", swtimedep_background);
    tdep_t_lev->create_timedep_prof(input_nc, offset, timedep_dim, n_lev);
    tdep_p_lay = std::make_unique<Timedep<TF>>(master, grid, "p_lay", swtimedep_background);
    tdep_p_lay->create_timedep_prof(input_nc, offset, timedep_dim, n_lay);
    tdep_p_lev = std::make_unique<Timedep<TF>>(master, grid, "p_lev", swtimedep_background);
    tdep_p_lev->create_timedep_prof(input_nc, offset, timedep_dim, n_lev);
    tdep_h2o = std::make_unique<Timedep<TF>>(master, grid, "h2o_bg", swtimedep_background);
    tdep_h2o->create_timedep_prof(input_nc, offset, timedep_dim, n_lay);

    // gasses
    for (auto& it : tdep_gases)
        it.second->create_timedep_prof(input_nc, offset, timedep_dim, n_lay);

    //aerosols
    if (sw_aerosol && swtimedep_aerosol)
    {
        tdep_aermr01 = std::make_unique<Timedep<TF>>(master, grid, "aermr01_bg", swtimedep_background);
        tdep_aermr01->create_timedep_prof(input_nc, offset, timedep_dim, n_lay);
        tdep_aermr02 = std::make_unique<Timedep<TF>>(master, grid, "aermr02_bg", swtimedep_background);
        tdep_aermr02->create_timedep_prof(input_nc, offset, timedep_dim, n_lay);
        tdep_aermr03 = std::make_unique<Timedep<TF>>(master, grid, "aermr03_bg", swtimedep_background);
        tdep_aermr03->create_timedep_prof(input_nc, offset, timedep_dim, n_lay);
        tdep_aermr04 = std::make_unique<Timedep<TF>>(master, grid, "aermr04_bg", swtimedep_background);
        tdep_aermr04->create_timedep_prof(input_nc, offset, timedep_dim, n_lay);
        tdep_aermr05 = std::make_unique<Timedep<TF>>(master, grid, "aermr05_bg", swtimedep_background);
        tdep_aermr05->create_timedep_prof(input_nc, offset, timedep_dim, n_lay);
        tdep_aermr06 = std::make_unique<Timedep<TF>>(master, grid, "aermr06_bg", swtimedep_background);
        tdep_aermr06->create_timedep_prof(input_nc, offset, timedep_dim, n_lay);
        tdep_aermr07 = std::make_unique<Timedep<TF>>(master, grid, "aermr07_bg", swtimedep_background);
        tdep_aermr07->create_timedep_prof(input_nc, offset, timedep_dim, n_lay);
        tdep_aermr08 = std::make_unique<Timedep<TF>>(master, grid, "aermr08_bg", swtimedep_background);
        tdep_aermr08->create_timedep_prof(input_nc, offset, timedep_dim, n_lay);
        tdep_aermr09 = std::make_unique<Timedep<TF>>(master, grid, "aermr09_bg", swtimedep_background);
        tdep_aermr09->create_timedep_prof(input_nc, offset, timedep_dim, n_lay);
        tdep_aermr10 = std::make_unique<Timedep<TF>>(master, grid, "aermr10_bg", swtimedep_background);
        tdep_aermr10->create_timedep_prof(input_nc, offset, timedep_dim, n_lay);
        tdep_aermr11 = std::make_unique<Timedep<TF>>(master, grid, "aermr11_bg", swtimedep_background);
        tdep_aermr11->create_timedep_prof(input_nc, offset, timedep_dim, n_lay);
    }

    // Prepare statistics.
    const std::string group_name = "radiation";
    stats.add_dimension("lay", n_lay);
    stats.add_dimension("lev", n_lev);

    // temperature, pressure and moisture
    stats.add_prof("t_lay_bg", "temperature at model layers of background profile", "K", "lay", group_name);
    stats.add_prof("t_lev_bg", "temperature at model levels of background profile", "K", "lev", group_name);
    stats.add_prof("p_lay_bg", "pressure at model layers of background profile", "Pa", "lay", group_name);
    stats.add_prof("p_lev_bg", "pressure at model levels of background profile", "Pa", "lev", group_name);
    stats.add_prof("h2o_bg", "h2o", "kg Kg-1", "lay", group_name);

    // gasses
    if (std::find(gaslist.begin(), gaslist.end(), "o3") != gaslist.end())
        stats.add_prof("o3_bg", "o3", "kg Kg-1", "lay", group_name);

    // aerosols
    if (sw_aerosol && swtimedep_aerosol)
    {
        stats.add_prof("aermr01_bg", "Sea salt (0.03 - 0.5 um) mixing ratio", "kg Kg-1", "lay", group_name);
        stats.add_prof("aermr02_bg", "Sea salt (0.5 - 5 um) mixing ratio", "kg Kg-1", "lay", group_name);
        stats.add_prof("aermr03_bg", "Sea salt (5 - 20 um) mixing ratio", "kg Kg-1", "lay", group_name);
        stats.add_prof("aermr04_bg", "Dust (0.03 - 0.55 um) mixing ratio", "kg Kg-1", "lay", group_name);
        stats.add_prof("aermr05_bg", "Dust (0.55 - 0.9 um) mixing ratio", "kg Kg-1", "lay", group_name);
        stats.add_prof("aermr06_bg", "Dust (0.9 - 20 um) mixing ratio", "kg Kg-1", "lay", group_name);
        stats.add_prof("aermr07_bg", "Organic matter (hydrophilic) mixing ratio", "kg Kg-1", "lay", group_name);
        stats.add_prof("aermr08_bg", "Organic matter (hydrophobic) mixing ratio", "kg Kg-1", "lay", group_name);
        stats.add_prof("aermr09_bg", "Black carbon (hydrophilic) mixing ratio", "kg Kg-1", "lay", group_name);
        stats.add_prof("aermr10_bg", "Black carbon (hydrophobic) mixing ratio", "kg Kg-1", "lay", group_name);
        stats.add_prof("aermr11_bg", "Sulfates mixing ratio", "kg Kg-1", "lay", group_name);
    }
}

template <typename TF>
void Background<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (!swtimedep_background)
        return;

    const bool do_radiation = ((timeloop.get_itime() % idt_rad == 0) && !timeloop.in_substep()) ;

    if (do_radiation)
    {
        // temperature, pressure and moisture
        tdep_t_lay   ->update_time_dependent_prof(t_lay, timeloop, n_lay);
        tdep_t_lev   ->update_time_dependent_prof(t_lev, timeloop, n_lev);
        tdep_p_lay   ->update_time_dependent_prof(p_lay, timeloop, n_lay);
        tdep_p_lev   ->update_time_dependent_prof(p_lev, timeloop, n_lev);
        tdep_h2o     ->update_time_dependent_prof(h2o, timeloop, n_lay);

        // gasses
        for (auto& it : tdep_gases)
            it.second->update_time_dependent_prof(gasprofs.at(it.first), timeloop, n_lay);

        // aerosols
        if (sw_aerosol && swtimedep_aerosol)
        {
            tdep_aermr01 ->update_time_dependent_prof(aermr01, timeloop, n_lay);
            tdep_aermr02 ->update_time_dependent_prof(aermr02, timeloop, n_lay);
            tdep_aermr03 ->update_time_dependent_prof(aermr03, timeloop, n_lay);
            tdep_aermr04 ->update_time_dependent_prof(aermr04, timeloop, n_lay);
            tdep_aermr05 ->update_time_dependent_prof(aermr05, timeloop, n_lay);
            tdep_aermr06 ->update_time_dependent_prof(aermr06, timeloop, n_lay);
            tdep_aermr07 ->update_time_dependent_prof(aermr07, timeloop, n_lay);
            tdep_aermr08 ->update_time_dependent_prof(aermr08, timeloop, n_lay);
            tdep_aermr09 ->update_time_dependent_prof(aermr09, timeloop, n_lay);
            tdep_aermr10 ->update_time_dependent_prof(aermr10, timeloop, n_lay);
            tdep_aermr11 ->update_time_dependent_prof(aermr11, timeloop, n_lay);
        }
    }
}


template<typename TF>
void Background<TF>::exec_stats(Stats<TF>& stats)
{
    if (!swtimedep_background)
        return;

    auto& gd = grid.get_grid_data();

    // temperature, pressure and moisture
    stats.set_prof_background("t_lay_bg", t_lay);
    stats.set_prof_background("t_lev_bg", t_lev);
    stats.set_prof_background("p_lay_bg", p_lay);
    stats.set_prof_background("p_lev_bg", p_lev);
    stats.set_prof_background("h2o_bg", h2o);

    // gasses
    if (std::find(gaslist.begin(), gaslist.end(), "o3") != gaslist.end())
        stats.set_prof_background("o3_bg", gasprofs.at("o3"));

    //aerosols
    if (sw_aerosol && swtimedep_aerosol)
    {
        stats.set_prof_background("aermr01_bg", aermr01);
        stats.set_prof_background("aermr02_bg", aermr02);
        stats.set_prof_background("aermr03_bg", aermr03);
        stats.set_prof_background("aermr04_bg", aermr04);
        stats.set_prof_background("aermr05_bg", aermr05);
        stats.set_prof_background("aermr06_bg", aermr06);
        stats.set_prof_background("aermr07_bg", aermr07);
        stats.set_prof_background("aermr08_bg", aermr08);
        stats.set_prof_background("aermr09_bg", aermr09);
        stats.set_prof_background("aermr10_bg", aermr10);
        stats.set_prof_background("aermr11_bg", aermr11);
    }
}

template<typename TF>
void Background<TF>::get_tpm(Array<Float,2>& t_lay_col, Array<Float,2>& t_lev_col,
                            Array<Float,2>& p_lay_col, Array<Float,2>& p_lev_col,
                             Gas_concs& gas_concs_col)
{
    for (int k=0; k<n_lay; ++k)
    {
        t_lay_col({1, k+1}) = t_lay[k];
        p_lay_col({1, k+1}) = p_lay[k];
    }
    for (int k=0; k<n_lev; ++k)
    {
        t_lev_col({1, k+1}) = t_lev[k];
        p_lev_col({1, k+1}) = p_lev[k];
    }

    Array<Float,2> h2o_bg_a({1, n_lay});
    for (int k=0; k<n_lay; ++k)
    {
        h2o_bg_a({1, k+1}) = h2o[k];
    }
    gas_concs_col.set_vmr("h2o", h2o_bg_a);
}

template<typename TF>
void Background<TF>::get_gasses(Gas_concs& gas_concs_col)
{
    for (auto& it : tdep_gases)
    {
        Array<Float,2> tmp_array({1, n_lay});
        for (int k=0; k<n_lay; ++k)
        {
            tmp_array({1, k+1}) = gasprofs.at(it.first)[k];
        }
        gas_concs_col.set_vmr(it.first, tmp_array);
    }

}

template<typename TF>
void Background<TF>::get_aerosols(Aerosol_concs& aerosol_concs_col)
{
    Array<Float,2> aermr01_bg_a({1, n_lay});
    Array<Float,2> aermr02_bg_a({1, n_lay});
    Array<Float,2> aermr03_bg_a({1, n_lay});
    Array<Float,2> aermr04_bg_a({1, n_lay});
    Array<Float,2> aermr05_bg_a({1, n_lay});
    Array<Float,2> aermr06_bg_a({1, n_lay});
    Array<Float,2> aermr07_bg_a({1, n_lay});
    Array<Float,2> aermr08_bg_a({1, n_lay});
    Array<Float,2> aermr09_bg_a({1, n_lay});
    Array<Float,2> aermr10_bg_a({1, n_lay});
    Array<Float,2> aermr11_bg_a({1, n_lay});

    for (int k=0; k<n_lay; ++k)
    {
        aermr01_bg_a({1, k+1}) = aermr01[k];
        aermr02_bg_a({1, k+1}) = aermr02[k];
        aermr03_bg_a({1, k+1}) = aermr03[k];
        aermr04_bg_a({1, k+1}) = aermr04[k];
        aermr05_bg_a({1, k+1}) = aermr05[k];
        aermr06_bg_a({1, k+1}) = aermr06[k];
        aermr07_bg_a({1, k+1}) = aermr07[k];
        aermr08_bg_a({1, k+1}) = aermr08[k];
        aermr09_bg_a({1, k+1}) = aermr09[k];
        aermr10_bg_a({1, k+1}) = aermr10[k];
        aermr11_bg_a({1, k+1}) = aermr11[k];
    }

    aerosol_concs_col.set_vmr("aermr01", aermr01_bg_a);
    aerosol_concs_col.set_vmr("aermr02", aermr02_bg_a);
    aerosol_concs_col.set_vmr("aermr03", aermr03_bg_a);
    aerosol_concs_col.set_vmr("aermr04", aermr04_bg_a);
    aerosol_concs_col.set_vmr("aermr05", aermr05_bg_a);
    aerosol_concs_col.set_vmr("aermr06", aermr06_bg_a);
    aerosol_concs_col.set_vmr("aermr07", aermr07_bg_a);
    aerosol_concs_col.set_vmr("aermr08", aermr08_bg_a);
    aerosol_concs_col.set_vmr("aermr09", aermr09_bg_a);
    aerosol_concs_col.set_vmr("aermr10", aermr10_bg_a);
    aerosol_concs_col.set_vmr("aermr11", aermr11_bg_a);
}


#ifdef FLOAT_SINGLE
template class Background<float>;
#else
template class Background<double>;
#endif
