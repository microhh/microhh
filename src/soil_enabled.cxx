/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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
#include <iostream>
#include <cmath>

#include "master.h"
#include "grid.h"
#include "soil_grid.h"
#include "fields.h"
#include "soil_field3d.h"
#include "stats.h"
#include "cross.h"
#include "constants.h"
#include "netcdf_interface.h"
#include "constants.h"

#include "soil.h"
#include "soil_enabled.h"

using namespace Constants;

namespace
{
    template<typename TF>
    void init_soil_homogeneous(
            TF* const restrict soil_fld, const TF* const restrict soil_prof,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int isize, const int ijsize)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*isize + k*ijsize;
                    soil_fld[ijk] = soil_prof[k-kstart];
                }
    }

    template<typename TF>
    inline TF calc_diffusivity_vg(
            const TF vg_a, const TF vg_l, const TF vg_m, const TF gamma_sat,
            const TF theta_res, const TF theta_sat, const TF theta_norm)
    {
        const TF vg_mi = TF(1) / vg_m;

        return (TF(1)-vg_m)*gamma_sat / (vg_a * vg_m * (theta_sat-theta_res)) * pow(theta_norm, (vg_l-vg_mi)) *
             (  pow((TF(1)-pow(theta_norm, vg_mi)), -vg_m) + pow((TF(1)-pow(theta_norm, vg_mi)), vg_m) - TF(2));
    }

    template<typename TF>
    void calc_soil_properties(
            TF* const restrict kappa_theta_min, TF* const restrict kappa_theta_max,
            TF* const restrict gamma_theta_min, TF* const restrict gamma_theta_max,
            TF* const restrict vg_m,
            TF* const restrict gamma_T_dry, TF* const restrict rho_C,
            const TF* const restrict vg_a, const TF* const restrict vg_l, const TF* const restrict vg_n,
            const TF* const restrict gamma_theta_sat,
            const TF* const restrict theta_res, const TF* const restrict theta_sat,
            const TF* const restrict theta_fc,
            const int table_size)
    {
        for (int i = 0; i<table_size; ++i)
        {
            // van Genuchten parameter `m`
            vg_m[i] = (TF(1) - (TF(1) / vg_n[i]));

            // Min/max values diffusivity soil moisture
            const TF theta_norm_min = (TF(1.001) * theta_res[i] - theta_res[i]) / (theta_sat[i] - theta_res[i]);
            const TF theta_norm_max = (TF(0.999) * theta_sat[i] - theta_res[i]) / (theta_sat[i] - theta_res[i]);

            kappa_theta_min[i] = calc_diffusivity_vg(
                    vg_a[i], vg_l[i], vg_m[i], gamma_theta_sat[i], theta_res[i], theta_sat[i], theta_norm_min);
            kappa_theta_max[i] = calc_diffusivity_vg(
                    vg_a[i], vg_l[i], vg_m[i], gamma_theta_sat[i], theta_res[i], theta_sat[i], theta_norm_max);

            // Min/max values conductivity soil moisture
            gamma_theta_min[i] = TF(0);
            gamma_theta_max[i] = gamma_theta_sat[i];

            // Conductivity temperature
            const TF rho_solid = TF(2700);                          // Density of dry solid soil (kg m-3); PL98, eq. 6
            const TF rho_dry = (TF(1) - theta_sat[i]) * rho_solid;      // Density of soil (kg m-3)

            gamma_T_dry[i] = (TF(0.135) * rho_dry + TF(64.7)) / (rho_solid - TF(0.947) * rho_dry);
            rho_C[i] = (TF(1) - theta_sat[i]) * Constants::rho_C_matrix<TF> + theta_fc[i] * Constants::rho_C_water<TF>;
        }
    }
}


template<typename TF>
Soil_enabled<TF>::Soil_enabled(
        Master& masterin, Grid<TF>& gridin, Soil_grid<TF>& soilgridin,
        Fields<TF>& fieldsin, Input& inputin) :
    Soil<TF>(masterin, gridin, soilgridin, fieldsin, inputin)
{
    sw_soil = Soil_type::Enabled;

    sw_interactive = inputin.get_item<bool>("soil", "sw_interactive", "", false);
    sw_homogeneous = inputin.get_item<bool>("soil", "sw_homogeneous", "", true);

    // Checks on input & limitations
    //if (sw_interactive)
    //    throw std::runtime_error("Interactive soil not (yet) implemented");
    if (!sw_homogeneous)
        throw std::runtime_error("Heterogeneous soil input not (yet) implemented");

    // Create soil fields (temperature and volumetric water content)
    fields.init_prognostic_soil_field("t_soil",     "Soil temperature", "K");
    fields.init_prognostic_soil_field("theta_soil", "Soil volumetric water content", "m3 m-3");

    // Open NetCDF file with soil lookup table
    nc_lookup_table = std::make_shared<Netcdf_file>(master, "van_genuchten_parameters.nc", Netcdf_mode::Read);
}

template<typename TF>
Soil_enabled<TF>::~Soil_enabled()
{
}

template<typename TF>
void Soil_enabled<TF>::init()
{
    /*
       Allocate/resize the soil fields, properties, and grid definition.
    */
    auto& sgd = soil_grid.get_grid_data();

    // Resize the vectors which contain the soil properties
    soil_index.resize(sgd.ncells);

    if (sw_interactive)
    {
        diffusivity.resize   (sgd.ncells);
        diffusivity_h.resize (sgd.ncellsh);
        conductivity.resize  (sgd.ncells);
        conductivity_h.resize(sgd.ncellsh);
        source.resize        (sgd.ncells);

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
}

template<typename TF>
void Soil_enabled<TF>::create_cold_start(Input& input, Netcdf_handle& input_nc)
{
    /*
       Create the prognostic soil fields, initialised either
       homogeneous from the input NetCDF file, or heterogeneous
       from "other" (yet to be defined..) sources.
       This routine is only called in the `init` phase of the model (from model.cxx),
       in the `run` phase these fields are read from the restart files.
     */
    auto& agd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    // Init the soil variables
    if (sw_homogeneous)
    {
        // Read initial profiles from input NetCDF file
        Netcdf_group& soil_group = input_nc.get_group("soil");

        std::vector<TF> t_prof(sgd.ktot);
        std::vector<TF> theta_prof(sgd.ktot);

        soil_group.get_variable(t_prof, "t_soil", {0}, {sgd.ktot});
        soil_group.get_variable(theta_prof, "theta_soil", {0}, {sgd.ktot});

        // Initialise soil as spatially homogeneous
        init_soil_homogeneous(
                fields.sps.at("t_soil")->fld.data(), t_prof.data(),
                agd.istart, agd.iend,
                agd.jstart, agd.jend,
                sgd.kstart, sgd.kend,
                agd.icells, agd.ijcells);

        init_soil_homogeneous(
                fields.sps.at("theta_soil")->fld.data(), theta_prof.data(),
                agd.istart, agd.iend,
                agd.jstart, agd.jend,
                sgd.kstart, sgd.kend,
                agd.icells, agd.ijcells);
    }
}

template<typename TF>
void Soil_enabled<TF>::create_fields_grid_stats(
        Input& input, Netcdf_handle& input_nc, Stats<TF>& stats, Cross<TF>& cross)
{
    /*
       Create/set the non-prognostic fields (soil type, ...) from the input files,
       calculate/define the soil grid, and init the soil statistics and cross-sections.
    */
    auto& agd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    // Init soil properties
    if (sw_homogeneous)
    {
        Netcdf_group& soil_group = input_nc.get_group("soil");
        std::vector<int> soil_index_prof(sgd.ktot);

        soil_group.get_variable<int>(soil_index_prof, "soil_index", {0}, {sgd.ktot});

        init_soil_homogeneous<int>(
                soil_index.data(), soil_index_prof.data(),
                agd.istart, agd.iend,
                agd.jstart, agd.jend,
                sgd.kstart, sgd.kend,
                agd.icells, agd.ijcells);
    }

    // Read lookup table soil
    nc_lookup_table->get_variable<TF>(theta_res, "theta_res", {0}, {lookup_table_size});
    nc_lookup_table->get_variable<TF>(theta_wp,  "theta_wp",  {0}, {lookup_table_size});
    nc_lookup_table->get_variable<TF>(theta_fc,  "theta_fc",  {0}, {lookup_table_size});
    nc_lookup_table->get_variable<TF>(theta_sat, "theta_sat", {0}, {lookup_table_size});

    nc_lookup_table->get_variable<TF>(gamma_theta_sat, "gamma_sat", {0}, {lookup_table_size});

    nc_lookup_table->get_variable<TF>(vg_a, "alpha", {0}, {lookup_table_size});
    nc_lookup_table->get_variable<TF>(vg_l, "l",     {0}, {lookup_table_size});
    nc_lookup_table->get_variable<TF>(vg_n, "n",     {0}, {lookup_table_size});

    // Calculate derived properties
    calc_soil_properties(
            kappa_theta_min.data(), kappa_theta_max.data(),
            gamma_theta_min.data(), gamma_theta_max.data(), vg_m.data(),
            gamma_T_dry.data(), rho_C.data(),
            vg_a.data(), vg_l.data(), vg_n.data(), gamma_theta_sat.data(),
            theta_res.data(), theta_sat.data(), theta_fc.data(), lookup_table_size);


    // Init the soil statistics
    if (stats.get_switch())
    {
        std::string group_name = "soil";

        // Add soil dimensions to each of the statistics masks
        auto& masks = stats.get_masks();
        for (auto& mask : masks)
        {
            auto& m = mask.second;

            // Add dimensions to NetCDF file
            m.data_file->add_dimension("zs",  sgd.ktot);
            m.data_file->add_dimension("zsh", sgd.ktot+1);

            // Write the attributes
            Netcdf_variable<TF> zs_var = m.data_file->template add_variable<TF>("zs", {"zs"});
            zs_var.add_attribute("units", "m");
            zs_var.add_attribute("long_name", "Full level soil height");

            Netcdf_variable<TF> zsh_var = m.data_file->template add_variable<TF>("zsh", {"zsh"});
            zsh_var.add_attribute("units", "m");
            zsh_var.add_attribute("long_name", "Half level soil height");

            // Write the grid levels
            zs_var .insert(sgd.z,  {0});
            zsh_var.insert(sgd.zh, {0});

            m.data_file->sync();
        }

        // Add the statistics variables
        stats.add_prof("t_soil", "Soil temperature", "K", "zs", group_name);
        stats.add_prof("theta_soil", "Soil volumetric water content", "-", "zs", group_name);
    }

    // Init the soil cross-sections
    if (cross.get_switch())
    {
        std::vector<std::string> allowed_crossvars = {"t_soil", "theta_soil"};
        crosslist = cross.get_enabled_variables(allowed_crossvars);
    }
}

template<typename TF>
void Soil_enabled<TF>::exec_stats(Stats<TF>& stats)
{
    const TF offset = 0;

    stats.calc_stats_soil("t_soil",     fields.sps.at("t_soil")->fld,     offset);
    stats.calc_stats_soil("theta_soil", fields.sps.at("theta_soil")->fld, offset);
}

template<typename TF>
void Soil_enabled<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    for (auto& it : crosslist)
    {
        if (it == "t_soil")
            cross.cross_soil(fields.sps.at("t_soil")->fld.data(), it, iotime);
        else if (it == "theta_soil")
            cross.cross_soil(fields.sps.at("t_soil")->fld.data(), it, iotime);
    }
}

template<typename TF>
void Soil_enabled<TF>::save_prognostic_fields(const int itime)
{
    auto field3d_io = Field3d_io<TF>(master, grid);
    auto& sgd = soil_grid.get_grid_data();

    const TF no_offset = 0.;
    int nerror = 0;

    auto tmp1 = fields.get_tmp();
    auto tmp2 = fields.get_tmp();

    // Soil temperature
    char filename[256];
    std::sprintf(filename, "%s.%07d", "t_soil", itime);
    master.print_message("Saving \"%s\" ... ", filename);

    if (field3d_io.save_field3d(
                fields.sps.at("t_soil")->fld.data(),
                tmp1->fld.data(), tmp2->fld.data(),
                filename, no_offset,
                sgd.kstart, sgd.kend))
    {
        master.print_message("FAILED\n");
        ++nerror;
    }
    else
        master.print_message("OK\n");

    // Soil moisture
    std::sprintf(filename, "%s.%07d", "theta_soil", itime);
    master.print_message("Saving \"%s\" ... ", filename);

    if (field3d_io.save_field3d(
                fields.sps.at("theta_soil")->fld.data(),
                tmp1->fld.data(), tmp2->fld.data(),
                filename, no_offset,
                sgd.kstart, sgd.kend))
    {
        master.print_message("FAILED\n");
        ++nerror;
    }
    else
        master.print_message("OK\n");

    fields.release_tmp(tmp1);
    fields.release_tmp(tmp2);

    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("Error saving soil fields");
}

template<typename TF>
void Soil_enabled<TF>::load_prognostic_fields(const int itime)
{
    auto field3d_io = Field3d_io<TF>(master, grid);
    auto& sgd = soil_grid.get_grid_data();

    const TF no_offset = 0.;
    int nerror = 0;

    auto tmp1 = fields.get_tmp();
    auto tmp2 = fields.get_tmp();

    // Soil temperature
    char filename[256];
    std::sprintf(filename, "%s.%07d", "t_soil", itime);
    master.print_message("Loading \"%s\" ... ", filename);

    if (field3d_io.load_field3d(
                fields.sps.at("t_soil")->fld.data(),
                tmp1->fld.data(), tmp2->fld.data(),
                filename, no_offset,
                sgd.kstart, sgd.kend))
    {
        master.print_message("FAILED\n");
        ++nerror;
    }
    else
        master.print_message("OK\n");

    // Soil moisture
    std::sprintf(filename, "%s.%07d", "theta_soil", itime);
    master.print_message("Loading \"%s\" ... ", filename);

    if (field3d_io.load_field3d(
                fields.sps.at("theta_soil")->fld.data(),
                tmp1->fld.data(), tmp2->fld.data(),
                filename, no_offset,
                sgd.kstart, sgd.kend))
    {
        master.print_message("FAILED\n");
        ++nerror;
    }
    else
        master.print_message("OK\n");

    fields.release_tmp(tmp1);
    fields.release_tmp(tmp2);

    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("Error loading soil fields");
}

template class Soil_enabled<double>;
template class Soil_enabled<float>;
