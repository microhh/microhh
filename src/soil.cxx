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
#include <algorithm>
#include <cmath>
#include <math.h>

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

using namespace Constants;

namespace
{
    template<typename TF>
    inline TF calc_diffusivity_vg(
            const TF vg_a, const TF vg_l, const TF vg_m, const TF gamma_sat,
            const TF theta_res, const TF theta_sat, const TF theta_norm)
    {
        const TF vg_mi = TF(1) / vg_m;

        return (TF(1) - vg_m) * gamma_sat / (vg_a * vg_m * (theta_sat - theta_res)) * pow(theta_norm, (vg_l - vg_mi)) *
               (pow((TF(1) - pow(theta_norm, vg_mi)), -vg_m) + pow((TF(1) - pow(theta_norm, vg_mi)), vg_m) - TF(2));
    }

    template<typename TF>
    inline TF calc_conductivity_vg(
            const TF theta_norm, const TF vg_l, const TF vg_m, const TF gamma_sat)
    {
        return gamma_sat * pow(theta_norm, vg_l) * pow((TF(1) - pow((TF(1) - pow(theta_norm, (1. / vg_m))), vg_m)), 2);
    }

    template<typename TF>
    void init_soil_homogeneous(
            TF* const restrict soil_fld,
            const TF* const restrict soil_prof,
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
                    const int ijk = i+j * isize + k*ijsize;
                    soil_fld[ijk] = soil_prof[k-kstart];
                }
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
        for (int i=0; i<table_size; ++i)
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
            const TF rho_solid = TF(2700);  // Density of dry solid soil (kg m-3); PL98, eq. 6
            const TF rho_dry = (TF(1) - theta_sat[i]) * rho_solid;  // Density of soil (kg m-3)

            gamma_T_dry[i] = (TF(0.135) * rho_dry + TF(64.7)) / (rho_solid - TF(0.947) * rho_dry);
            rho_C[i] = (TF(1) - theta_sat[i]) * Constants::rho_C_matrix<TF> + theta_fc[i] * Constants::rho_C_water<TF>;
        }
    }

    template<typename TF>
    void calc_thermal_properties(
            TF* const restrict kappa,
            TF* const restrict gamma,
            const int* const restrict soil_index,
            const TF* const restrict theta,
            const TF* const restrict theta_sat,
            const TF* const restrict gamma_dry,
            const TF* const restrict rho_C,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    const int si = soil_index[ijk];

                    // Heat conductivity at saturation (from IFS code..)
                    const TF lambda_T_sat = pow(Constants::gamma_T_matrix<TF>, (TF(1) - theta_sat[si]))
                                            * pow(Constants::gamma_T_water<TF>, theta[ijk])
                                            * pow(TF(2.2), (theta_sat[si] - theta[ijk]));

                    // Kersten number for fine soils [IFS eq 8.64] (-)
                    const TF kersten = log10(std::max(TF(0.1), theta[ijk] / theta_sat[si])) + TF(1);

                    // Heat conductivity soil [IFS eq 8.62] (W m-1 K-1)
                    gamma[ijk] = kersten * (lambda_T_sat - gamma_dry[si]) + gamma_dry[si];

                    // Heat diffusivity (m2 s-1)
                    kappa[ijk] = gamma[ijk] / rho_C[si];
                }
    }

    template<typename TF>
    void calc_hydraulic_properties(
            TF* const restrict kappa,
            TF* const restrict gamma,
            const int* const restrict soil_index,
            const TF* const restrict theta,
            const TF* const restrict theta_sat,
            const TF* const restrict theta_res,
            const TF* const restrict vg_a,
            const TF* const restrict vg_l,
            const TF* const restrict vg_m,
            const TF* const restrict gamma_sat,
            const TF* const restrict gamma_min,
            const TF* const restrict gamma_max,
            const TF* const restrict kappa_min,
            const TF* const restrict kappa_max,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    const int si = soil_index[ijk];

                    // Limit soil moisture just above the residual soil moisture content
                    const TF theta_lim = std::max(theta[ijk], TF(1.001) * theta_res[si]);

                    // Dimensionless soil water content
                    const TF theta_norm = (theta_lim - theta_res[si]) / (theta_sat[si] - theta_res[si]);

                    // Calculate & limit the diffusivity
                    kappa[ijk] = calc_diffusivity_vg(
                            vg_a[si], vg_l[si], vg_m[si], gamma_sat[si],
                            theta_res[si], theta_sat[si], theta_norm);
                    kappa[ijk] = std::max(std::min(kappa_max[si], kappa[ijk]), kappa_min[si]);

                    // Calculate & limit the conductivity
                    gamma[ijk] = calc_conductivity_vg(
                            theta_norm, vg_l[si], vg_m[si], gamma_sat[si]);
                    gamma[ijk] = std::max(std::min(gamma_max[si], gamma[ijk]), gamma_min[si]);
                }
    }

    template<typename TF, Soil_interpolation_type interpolation_type>
    void interp_2_vertical(
            TF* const restrict fldh,
            const TF* const restrict fld,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int kk = ijcells;

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;

                    if (interpolation_type == Soil_interpolation_type::Mean)
                        fldh[ijk] = TF(0.5) * (fld[ijk] + fld[ijk-kk]);
                    else if(interpolation_type == Soil_interpolation_type::Max)
                        fldh[ijk] = std::max(fld[ijk], fld[ijk-kk]);
                }
    }

    template<typename TF>
    void set_bcs_temperature(
            TF* const restrict flux_top,
            TF* const restrict flux_bot,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij = i + j*icells;
                flux_top[ij] = TF(0);    // Eventually: G/rho
                flux_bot[ij] = TF(0);
            }
    }

    template<typename TF, bool sw_free_drainage>
    void set_bcs_moisture(
            TF* const restrict flux_top,
            TF* const restrict flux_bot,
            TF* const restrict conductivity_h,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int kk = ijcells;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij = i + j*icells;
                flux_top[ij] = TF(0);    // Eventually: LE & infiltration
                flux_bot[ij] = TF(0);

                // Set free drainage bottom BC:
                const int ijk = ij + kstart*ijcells;
                if (sw_free_drainage)
                    conductivity_h[ijk] = conductivity_h[ijk+kk];
                else
                    conductivity_h[ijk] = TF(0);
            }
    }

    template<typename TF, bool sw_source_term, bool sw_conductivity_term>
    void diff_explicit(
            TF* const restrict tend,
            const TF* const restrict fld,
            const TF* const restrict kappa_h,
            const TF* const restrict gamma_h,
            const TF* const restrict source,
            const TF* const restrict flux_top,
            const TF* const restrict flux_bot,
            const TF* const restrict dzi,
            const TF* const restrict dzhi,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int kk = ijcells;
        int k;

        // Bottom soil level
        k = kstart;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                const int ijk = ij + k*ijcells;

                tend[ijk] += ((kappa_h[ijk+kk] * (fld[ijk+kk] - fld[ijk]) * dzhi[k+1]) + flux_bot[ij])*dzi[k];

                if (sw_conductivity_term)
                    tend[ijk] += (gamma_h[ijk+kk] - gamma_h[ijk]) * dzi[k];
                if (sw_source_term)
                    tend[ijk] += source[ijk];
            }

        // Top soil level
        k = kend-1;
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                const int ijk = ij + k*ijcells;

                tend[ijk] += (-flux_top[ij] - (kappa_h[ijk] * (fld[ijk] - fld[ijk-kk]) * dzhi[k]))*dzi[k];

                if (sw_conductivity_term)
                    tend[ijk] -= gamma_h[ijk] * dzi[k];
                if (sw_source_term)
                    tend[ijk] += source[ijk];
            }

        // Interior
        for (int k=kstart+1; k<kend-1; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;

                    tend[ijk] += ((kappa_h[ijk+kk] * (fld[ijk+kk] - fld[ijk   ]) * dzhi[k+1])
                               - (kappa_h[ijk   ] * (fld[ijk   ] - fld[ijk-kk]) * dzhi[k  ])) * dzi[k];

                    if (sw_conductivity_term)
                        tend[ijk] += (gamma_h[ijk+kk] - gamma_h[ijk]) * dzi[k];
                    if (sw_source_term)
                        tend[ijk] += source[ijk];
                }
    }
}

template<typename TF>
Soil<TF>::Soil(
        Master& masterin, Grid<TF>& gridin, Soil_grid<TF>& soilgridin,
        Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), soil_grid(soilgridin), fields(fieldsin) 
{
    sw_soil = inputin.get_item<bool>("soil", "sw_soil", "", false);
    
    if (sw_soil)
    {
        sw_homogeneous   = inputin.get_item<bool>("soil", "sw_homogeneous", "", true);
        sw_free_drainage = inputin.get_item<bool>("soil", "sw_free_drainage", "", true);

        // Checks on input & limitations
        if (!sw_homogeneous)
            throw std::runtime_error("Heterogeneous soil input not (yet) implemented");

        // Create soil fields (temperature and volumetric water content)
        fields.init_prognostic_soil_field("t",     "Soil temperature", "K");
        fields.init_prognostic_soil_field("theta", "Soil volumetric water content", "m3 m-3");

        // Open NetCDF file with soil lookup table
        nc_lookup_table = std::make_shared<Netcdf_file>(master, "van_genuchten_parameters.nc", Netcdf_mode::Read);
    }
}

template<typename TF>
Soil<TF>::~Soil()
{
}

template<typename TF>
void Soil<TF>::init()
{
    /*
       Allocate/resize the soil fields, properties, and grid definition.
    */
    if (!sw_soil)
        return;

    auto& sgd = soil_grid.get_grid_data();

    // Resize the vectors which contain the soil properties
    soil_index.resize(sgd.ncells);

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

template<typename TF>
void Soil<TF>::create_cold_start(Input& input, Netcdf_handle& input_nc)
{
    /*
       Create the prognostic soil fields, initialised either
       homogeneous from the input NetCDF file, or heterogeneous
       from "other" (yet to be defined..) sources.
       This routine is only called in the `init` phase of the model (from model.cxx),
       in the `run` phase these fields are read from the restart files.
     */
    if (!sw_soil)
        return;

    auto& agd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    // Init the soil variables
    if (sw_homogeneous)
    {
        // Read initial profiles from input NetCDF file
        Netcdf_group& soil_group = input_nc.get_group("soil");

        std::vector<TF> t_prof(sgd.ktot);
        std::vector<TF> theta_prof(sgd.ktot);

        soil_group.get_variable(t_prof, "t", {0}, {sgd.ktot});
        soil_group.get_variable(theta_prof, "theta", {0}, {sgd.ktot});

        // Initialise soil as spatially homogeneous
        init_soil_homogeneous(
                fields.sps.at("t")->fld.data(), t_prof.data(),
                agd.istart, agd.iend,
                agd.jstart, agd.jend,
                sgd.kstart, sgd.kend,
                agd.icells, agd.ijcells);

        init_soil_homogeneous(
                fields.sps.at("theta")->fld.data(), theta_prof.data(),
                agd.istart, agd.iend,
                agd.jstart, agd.jend,
                sgd.kstart, sgd.kend,
                agd.icells, agd.ijcells);
    }
}

template<typename TF>
void Soil<TF>::create_fields_grid_stats(
        Input& input, Netcdf_handle& input_nc, Stats<TF>& stats, Cross<TF>& cross)
{
    /*
       Create/set the non-prognostic fields (soil type, ...) from the input files,
       calculate/define the soil grid, and init the soil statistics and cross-sections.
    */
    if (!sw_soil)
        return;

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
        stats.add_prof("t", "Soil temperature", "K", "zs", group_name);
        stats.add_prof("theta", "Soil volumetric water content", "-", "zs", group_name);
    }

    // Init the soil cross-sections
    if (cross.get_switch())
    {
        std::vector<std::string> allowed_crossvars = {"t_soil", "theta_soil"};
        crosslist = cross.get_enabled_variables(allowed_crossvars);
    }
}

template<typename TF>
void Soil<TF>::calc_tendencies()
{
    if (!sw_soil)
        return;

    auto& agd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    // Only soil moisture has a source and conductivity term
    const bool sw_source_term_t = false;
    const bool sw_conductivity_term_t = false;
    const bool sw_source_term_theta = true;
    const bool sw_conductivity_term_theta = true;

    //
    // Soil temperature
    //
    // Calculate the thermal diffusivity at full levels
    calc_thermal_properties(
            diffusivity.data(),
            conductivity.data(),
            soil_index.data(),
            fields.sps.at("theta")->fld.data(),
            theta_sat.data(),
            gamma_T_dry.data(),
            rho_C.data(),
            agd.istart, agd.iend,
            agd.jstart, agd.jend,
            sgd.kstart, sgd.kend,
            agd.icells, agd.ijcells);

    // Linear interpolation diffusivity to half levels
    interp_2_vertical<TF, Soil_interpolation_type::Mean>(
            diffusivity_h.data(),
            diffusivity.data(),
            agd.istart, agd.iend,
            agd.jstart, agd.jend,
            sgd.kstart, sgd.kend,
            agd.icells, agd.ijcells);

    // Set flux boundary conditions at top and bottom of soil column
    set_bcs_temperature(
            fields.sps.at("t")->flux_top.data(),
            fields.sps.at("t")->flux_bot.data(),
            agd.istart, agd.iend,
            agd.jstart, agd.jend,
            sgd.kstart, sgd.kend,
            agd.icells, agd.ijcells);

    // Calculate diffusive tendency
    diff_explicit<TF, sw_source_term_t, sw_conductivity_term_t>(
            fields.sts.at("t")->fld.data(),
            fields.sps.at("t")->fld.data(),
            diffusivity_h.data(),
            conductivity_h.data(),
            source.data(),
            fields.sps.at("t")->flux_top.data(),
            fields.sps.at("t")->flux_bot.data(),
            sgd.dzi.data(), sgd.dzhi.data(),
            agd.istart, agd.iend,
            agd.jstart, agd.jend,
            sgd.kstart, sgd.kend,
            agd.icells, agd.ijcells);

    //
    // Soil moisture
    //
    // Calculate the hydraulic diffusivity and conductivity at full levels
    calc_hydraulic_properties(
            diffusivity.data(),
            conductivity.data(),
            soil_index.data(),
            fields.sps.at("theta")->fld.data(),
            theta_sat.data(),
            theta_res.data(),
            vg_a.data(),
            vg_l.data(),
            vg_m.data(),
            gamma_theta_sat.data(),
            gamma_theta_min.data(),
            gamma_theta_max.data(),
            kappa_theta_min.data(),
            kappa_theta_max.data(),
            agd.istart, agd.iend,
            agd.jstart, agd.jend,
            sgd.kstart, sgd.kend,
            agd.icells, agd.ijcells);

    // Interpolation diffusivity and conductivity to half levels,
    // using the IFS method, which uses the max value from the
    // two surrounding grid points.
    interp_2_vertical<TF, Soil_interpolation_type::Max>(
            diffusivity_h.data(),
            diffusivity.data(),
            agd.istart, agd.iend,
            agd.jstart, agd.jend,
            sgd.kstart, sgd.kend,
            agd.icells, agd.ijcells);

    interp_2_vertical<TF, Soil_interpolation_type::Max>(
            conductivity_h.data(),
            conductivity.data(),
            agd.istart, agd.iend,
            agd.jstart, agd.jend,
            sgd.kstart, sgd.kend,
            agd.icells, agd.ijcells);

    // Set the flux boundary conditions at the top and bottom
    // of the soil layer, and a free drainage conditions at the bottom.
    if (sw_free_drainage)
        set_bcs_moisture<TF, true>(
                fields.sps.at("theta")->flux_top.data(),
                fields.sps.at("theta")->flux_bot.data(),
                conductivity_h.data(),
                agd.istart, agd.iend,
                agd.jstart, agd.jend,
                sgd.kstart, sgd.kend,
                agd.icells, agd.ijcells);
    else
        set_bcs_moisture<TF, false>(
                fields.sps.at("theta")->flux_top.data(),
                fields.sps.at("theta")->flux_bot.data(),
                conductivity_h.data(),
                agd.istart, agd.iend,
                agd.jstart, agd.jend,
                sgd.kstart, sgd.kend,
                agd.icells, agd.ijcells);

    // Calculate diffusive tendency
    diff_explicit<TF, sw_source_term_theta, sw_conductivity_term_theta>(
            fields.sts.at("theta")->fld.data(),
            fields.sps.at("theta")->fld.data(),
            diffusivity_h.data(),
            conductivity_h.data(),
            source.data(),
            fields.sps.at("theta")->flux_top.data(),
            fields.sps.at("theta")->flux_bot.data(),
            sgd.dzi.data(), sgd.dzhi.data(),
            agd.istart, agd.iend,
            agd.jstart, agd.jend,
            sgd.kstart, sgd.kend,
            agd.icells, agd.ijcells);

}

template<typename TF>
void Soil<TF>::exec_stats(Stats<TF>& stats)
{
    if (!sw_soil)
        return;

    const TF offset = 0;

    // Soil prognostic fields
    stats.calc_stats_soil("t",     fields.sps.at("t")->fld,     offset);
    stats.calc_stats_soil("theta", fields.sps.at("theta")->fld, offset);
}

template<typename TF>
void Soil<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    if (!sw_soil)
        return;

    for (auto& it : crosslist)
    {
        if (it == "t_soil")
            cross.cross_soil(fields.sps.at("t")->fld.data(), it, iotime);
        else if (it == "theta_soil")
            cross.cross_soil(fields.sps.at("theta")->fld.data(), it, iotime);
    }
}

template<typename TF>
void Soil<TF>::save_prognostic_fields(const int itime)
{
    if (!sw_soil)
        return;

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
                fields.sps.at("t")->fld.data(),
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
                fields.sps.at("theta")->fld.data(),
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
void Soil<TF>::load_prognostic_fields(const int itime)
{
    if (!sw_soil)
        return;

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
                fields.sps.at("t")->fld.data(),
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
                fields.sps.at("theta")->fld.data(),
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

template class Soil<double>;
template class Soil<float>;
