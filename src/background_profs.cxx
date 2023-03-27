//
// Created by Mirjam Tijhuis on 27/03/2023.
//

#include <iostream>
#include <algorithm>

#include "background_profs.h"

#include "timeloop.h"
#include "input.h"
#include "grid.h"
#include "netcdf_interface.h"
#include "stats.h"
//#include "constants.h"
#include "thermo.h"
#include "fields.h"
#include "timedep.h"

namespace
{
    // Kernels...
}

template<typename TF>
Background<TF>::Background(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
        master(masterin), grid(gridin), fields(fieldsin)
{
    // Read `.ini` settings.
    sw_update_background = inputin.get_item<bool>("radiation", "swupdatecolumn", "", false);
    sw_aerosol = inputin.get_item<bool>("aerosol", "swaerosol", "", false);
}

template <typename TF>
Background<TF>::~Background()
{
}

template <typename TF>
void Background<TF>::init()
{
    // Allocate (`.resize`) arrays.
    if (!sw_update_background)
        return;

    const TF n_era_layers = 136;

    std::cout << "init() timedepent background profiles" << std::endl;

    h2o.resize(n_era_layers);
    o3.resize(n_era_layers);
    aermr01.resize(n_era_layers);
    aermr02.resize(n_era_layers);
    aermr03.resize(n_era_layers);
    aermr04.resize(n_era_layers);
    aermr05.resize(n_era_layers);
    aermr06.resize(n_era_layers);
    aermr07.resize(n_era_layers);
    aermr08.resize(n_era_layers);
    aermr09.resize(n_era_layers);
    aermr10.resize(n_era_layers);
    aermr11.resize(n_era_layers);

}

template <typename TF>
void Background<TF>::create(Input& inputin, Netcdf_handle& input_nc, Stats<TF>& stats)
{
    // Read input from NetCDF and prepare statistics output.
    if (!sw_update_background)
        return;

    std::cout << "create() timedependent background profiles" << std::endl;

    // create time dependent profiles
    const TF offset = 0;
    std::string timedep_dim_ls = "time_ls";
    tdep_h2o = std::make_unique<Timedep<TF>>(master, grid, "h2o", inputin.get_item<bool>("radiation", "swupdatecolumn", "", false));
    tdep_h2o->create_timedep_background_prof(input_nc, offset, timedep_dim_ls, n_era_layers);
    tdep_o3 = std::make_unique<Timedep<TF>>(master, grid, "o3", inputin.get_item<bool>("radiation", "swupdatecolumn", "", false));
    tdep_o3->create_timedep_background_prof(input_nc, offset, timedep_dim_ls, n_era_layers);

    if (sw_aerosol)
    {
        std::string timedep_dim = "time_aerosols";
        tdep_aermr01 = std::make_unique<Timedep<TF>>(master, grid, "aermr01_bg", inputin.get_item<bool>("radiation", "swupdatecolumn", "", false));
        tdep_aermr01->create_timedep_background_prof(input_nc, offset, timedep_dim, n_era_layers);
        tdep_aermr02 = std::make_unique<Timedep<TF>>(master, grid, "aermr02_bg", inputin.get_item<bool>("radiation", "swupdatecolumn", "", false));
        tdep_aermr02->create_timedep_background_prof(input_nc, offset, timedep_dim, n_era_layers);
        tdep_aermr03 = std::make_unique<Timedep<TF>>(master, grid, "aermr03_bg", inputin.get_item<bool>("radiation", "swupdatecolumn", "", false));
        tdep_aermr03->create_timedep_background_prof(input_nc, offset, timedep_dim, n_era_layers);
        tdep_aermr04 = std::make_unique<Timedep<TF>>(master, grid, "aermr04_bg", inputin.get_item<bool>("radiation", "swupdatecolumn", "", false));
        tdep_aermr04->create_timedep_background_prof(input_nc, offset, timedep_dim, n_era_layers);
        tdep_aermr05 = std::make_unique<Timedep<TF>>(master, grid, "aermr05_bg", inputin.get_item<bool>("radiation", "swupdatecolumn", "", false));
        tdep_aermr05->create_timedep_background_prof(input_nc, offset, timedep_dim, n_era_layers);
        tdep_aermr06 = std::make_unique<Timedep<TF>>(master, grid, "aermr06_bg", inputin.get_item<bool>("radiation", "swupdatecolumn", "", false));
        tdep_aermr06->create_timedep_background_prof(input_nc, offset, timedep_dim, n_era_layers);
        tdep_aermr07 = std::make_unique<Timedep<TF>>(master, grid, "aermr07_bg", inputin.get_item<bool>("radiation", "swupdatecolumn", "", false));
        tdep_aermr07->create_timedep_background_prof(input_nc, offset, timedep_dim, n_era_layers);
        tdep_aermr08 = std::make_unique<Timedep<TF>>(master, grid, "aermr08_bg", inputin.get_item<bool>("radiation", "swupdatecolumn", "", false));
        tdep_aermr08->create_timedep_background_prof(input_nc, offset, timedep_dim, n_era_layers);
        tdep_aermr09 = std::make_unique<Timedep<TF>>(master, grid, "aermr09_bg", inputin.get_item<bool>("radiation", "swupdatecolumn", "", false));
        tdep_aermr09->create_timedep_background_prof(input_nc, offset, timedep_dim, n_era_layers);
        tdep_aermr10 = std::make_unique<Timedep<TF>>(master, grid, "aermr10_bg", inputin.get_item<bool>("radiation", "swupdatecolumn", "", false));
        tdep_aermr10->create_timedep_background_prof(input_nc, offset, timedep_dim, n_era_layers);
        tdep_aermr11 = std::make_unique<Timedep<TF>>(master, grid, "aermr11_bg", inputin.get_item<bool>("radiation", "swupdatecolumn", "", false));
        tdep_aermr11->create_timedep_background_prof(input_nc, offset, timedep_dim, n_era_layers);
    }

    // Prepare statistics.
    const std::string group_name = "default";
    stats.add_dimension("era_layers", n_era_layers);

    stats.add_prof("h2o_bg", "h2o", "kg Kg-1", "era_layers", group_name);
    stats.add_prof("o3_bg", "o3", "kg Kg-1", "era_layers", group_name);

    if (sw_aerosol)
    {
        stats.add_prof("aermr01_bg", "Sea salt (0.03 - 0.5 um) mixing ratio", "kg Kg-1", "era_layers", group_name);
        stats.add_prof("aermr02_bg", "Sea salt (0.5 - 5 um) mixing ratio", "kg Kg-1", "era_layers", group_name);
        stats.add_prof("aermr03_bg", "Sea salt (5 - 20 um) mixing ratio", "kg Kg-1", "era_layers", group_name);
        stats.add_prof("aermr04_bg", "Dust (0.03 - 0.55 um) mixing ratio", "kg Kg-1", "era_layers", group_name);
        stats.add_prof("aermr05_bg", "Dust (0.55 - 0.9 um) mixing ratio", "kg Kg-1", "era_layers", group_name);
        stats.add_prof("aermr06_bg", "Dust (0.9 - 20 um) mixing ratio", "kg Kg-1", "era_layers", group_name);
        stats.add_prof("aermr07_bg", "Organic matter (hydrophilic) mixing ratio", "kg Kg-1", "era_layers", group_name);
        stats.add_prof("aermr08_bg", "Organic matter (hydrophobic) mixing ratio", "kg Kg-1", "era_layers", group_name);
        stats.add_prof("aermr09_bg", "Black carbon (hydrophilic) mixing ratio", "kg Kg-1", "era_layers", group_name);
        stats.add_prof("aermr10_bg", "Black carbon (hydrophobic) mixing ratio", "kg Kg-1", "era_layers", group_name);
        stats.add_prof("aermr11_bg", "Sulfates mixing ratio", "kg Kg-1", "era_layers", group_name);
    }

}

#ifndef USECUDA
template <typename TF>
void Background<TF>::exec(Thermo<TF>& thermo)
{
    if (!sw_update_background)
        return;

    auto& gd = grid.get_grid_data();
}

template <typename TF>
void Background<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (!sw_update_background)
        return;

    tdep_h2o     ->update_time_dependent_background_prof(h2o, timeloop, n_era_layers);
    tdep_o3      ->update_time_dependent_background_prof(o3, timeloop, n_era_layers);

    if (sw_aerosol)
    {
        tdep_aermr01 ->update_time_dependent_background_prof(aermr01, timeloop, n_era_layers);
        tdep_aermr02 ->update_time_dependent_background_prof(aermr02, timeloop, n_era_layers);
        tdep_aermr03 ->update_time_dependent_background_prof(aermr03, timeloop, n_era_layers);
        tdep_aermr04 ->update_time_dependent_background_prof(aermr04, timeloop, n_era_layers);
        tdep_aermr05 ->update_time_dependent_background_prof(aermr05, timeloop, n_era_layers);
        tdep_aermr06 ->update_time_dependent_background_prof(aermr06, timeloop, n_era_layers);
        tdep_aermr07 ->update_time_dependent_background_prof(aermr07, timeloop, n_era_layers);
        tdep_aermr08 ->update_time_dependent_background_prof(aermr08, timeloop, n_era_layers);
        tdep_aermr09 ->update_time_dependent_background_prof(aermr09, timeloop, n_era_layers);
        tdep_aermr10 ->update_time_dependent_background_prof(aermr10, timeloop, n_era_layers);
        tdep_aermr11 ->update_time_dependent_background_prof(aermr11, timeloop, n_era_layers);
    }

}
#endif

template<typename TF>
void Background<TF>::exec_stats(Stats<TF>& stats)
{
    if (!sw_update_background)
        return;

    auto& gd = grid.get_grid_data();

    stats.set_prof_background("h2o_bg", h2o);
    stats.set_prof_background("o3_bg", o3);

    if (sw_aerosol)
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
void Background<TF>::get_radiation_fields(std::vector<TF>& h2o_bg, std::vector<TF>& o3_bg,
                                          std::vector<TF>& mr01, std::vector<TF>& mr02, std::vector<TF>& mr03,
                                          std::vector<TF>& mr04, std::vector<TF>& mr05, std::vector<TF>& mr06,
                                          std::vector<TF>& mr07, std::vector<TF>& mr08, std::vector<TF>& mr09,
                                          std::vector<TF>& mr10, std::vector<TF>& mr11)
{
    h2o_bg = h2o;
    o3_bg = o3;
    mr01 = aermr01;
    mr02 = aermr02;
    mr03 = aermr03;
    mr04 = aermr04;
    mr05 = aermr05;
    mr06 = aermr06;
    mr07 = aermr07;
    mr08 = aermr08;
    mr09 = aermr09;
    mr10 = aermr10;
    mr11 = aermr11;
}

template class Background<double>;
template class Background<float>;
