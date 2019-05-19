/*
 * MicroHH
 * Copyright (c) 2011-2019 Chiel van Heerwaarden
 * Copyright (c) 2011-2019 Thijs Heus
 * Copyright (c) 2014-2019 Bart van Stratum
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

#include <boost/algorithm/string.hpp>
#include <numeric>

#include "radiation_rrtmgp.h"
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "input.h"
#include "netcdf_interface.h"

#include "Array.h"
#include "Optical_props.h"
#include "Gas_optics.h"
#include "Gas_concs.h"
#include "Fluxes.h"
#include "Rte_lw.h"
#include "Rte_sw.h"
#include "Source_functions.h"

namespace
{
    std::vector<std::string> get_variable_string(
            const std::string& var_name,
            std::vector<int> i_count,
            Netcdf_handle& input_nc,
            const int string_len,
            bool trim=true)
    {
        // Multiply all elements in i_count.
        int total_count = std::accumulate(i_count.begin(), i_count.end(), 1, std::multiplies<>());

        // Add the string length as the rightmost dimension.
        i_count.push_back(string_len);

        // Multiply all elements in i_count.
        // int total_count_char = std::accumulate(i_count.begin(), i_count.end(), 1, std::multiplies<>());

        // Read the entire char array;
        std::vector<char> var_char;
        var_char = input_nc.get_variable<char>(var_name, i_count);

        std::vector<std::string> var;

        for (int n=0; n<total_count; ++n)
        {
            std::string s(var_char.begin()+n*string_len, var_char.begin()+(n+1)*string_len);
            if (trim)
                boost::trim(s);
            var.push_back(s);
        }

        return var;
    }

    template<typename TF>
    void load_gas_concs(Gas_concs<TF>& gas_concs, Netcdf_handle& input_nc)
    {
        // This part is contained in the create
        Netcdf_group rad_nc = input_nc.get_group("radiation");

        const int n_lay = rad_nc.get_dimension_size("lay");

        gas_concs.set_vmr("h2o",
                Array<TF,1>(rad_nc.get_variable<TF>("h2o", {n_lay}), {n_lay}));
        gas_concs.set_vmr("co2",
                rad_nc.get_variable<TF>("co2"));
        gas_concs.set_vmr("o3",
                Array<TF,1>(rad_nc.get_variable<TF>("o3", {n_lay}), {n_lay}));
        gas_concs.set_vmr("n2o",
                rad_nc.get_variable<TF>("n2o"));
        gas_concs.set_vmr("ch4",
                rad_nc.get_variable<TF>("ch4"));
        gas_concs.set_vmr("o2",
                rad_nc.get_variable<TF>("o2"));
        gas_concs.set_vmr("n2",
                rad_nc.get_variable<TF>("n2"));
    }

    Gas_optics<double> load_and_init_gas_optics(
            Master& master,
            const Gas_concs<double>& gas_concs,
            const std::string& coef_file)
    {
        // READ THE COEFFICIENTS FOR THE OPTICAL SOLVER.
        Netcdf_file coef_nc(master, coef_file, Netcdf_mode::Read);

        // Read k-distribution information.
        int n_temps = coef_nc.get_dimension_size("temperature");
        int n_press = coef_nc.get_dimension_size("pressure");
        int n_absorbers = coef_nc.get_dimension_size("absorber");
        int n_char = coef_nc.get_dimension_size("string_len");
        int n_minorabsorbers = coef_nc.get_dimension_size("minor_absorber");
        int n_extabsorbers = coef_nc.get_dimension_size("absorber_ext");
        int n_mixingfracs = coef_nc.get_dimension_size("mixing_fraction");
        int n_layers = coef_nc.get_dimension_size("atmos_layer");
        int n_bnds = coef_nc.get_dimension_size("bnd");
        int n_gpts = coef_nc.get_dimension_size("gpt");
        int n_pairs = coef_nc.get_dimension_size("pair");
        int n_minor_absorber_intervals_lower = coef_nc.get_dimension_size("minor_absorber_intervals_lower");
        int n_minor_absorber_intervals_upper = coef_nc.get_dimension_size("minor_absorber_intervals_upper");
        int n_contributors_lower = coef_nc.get_dimension_size("contributors_lower");
        int n_contributors_upper = coef_nc.get_dimension_size("contributors_upper");

        // Read gas names.
        Array<std::string,1> gas_names(
                get_variable_string("gas_names", {n_absorbers}, coef_nc, n_char, true), {n_absorbers});

        Array<int,3> key_species(
                coef_nc.get_variable<int>("key_species", {n_bnds, n_layers, 2}),
                {2, n_layers, n_bnds});
        Array<double,2> band_lims(coef_nc.get_variable<double>("bnd_limits_wavenumber", {n_bnds, 2}), {2, n_bnds});
        Array<int,2> band2gpt(coef_nc.get_variable<int>("bnd_limits_gpt", {n_bnds, 2}), {2, n_bnds});
        Array<double,1> press_ref(coef_nc.get_variable<double>("press_ref", {n_press}), {n_press});
        Array<double,1> temp_ref(coef_nc.get_variable<double>("temp_ref", {n_temps}), {n_temps});

        double temp_ref_p = coef_nc.get_variable<double>("absorption_coefficient_ref_P");
        double temp_ref_t = coef_nc.get_variable<double>("absorption_coefficient_ref_T");
        double press_ref_trop = coef_nc.get_variable<double>("press_ref_trop");

        Array<double,3> kminor_lower(
                coef_nc.get_variable<double>("kminor_lower", {n_temps, n_mixingfracs, n_contributors_lower}),
                {n_contributors_lower, n_mixingfracs, n_temps});
        Array<double,3> kminor_upper(
                coef_nc.get_variable<double>("kminor_upper", {n_temps, n_mixingfracs, n_contributors_upper}),
                {n_contributors_upper, n_mixingfracs, n_temps});

        Array<std::string,1> gas_minor(get_variable_string("gas_minor", {n_minorabsorbers}, coef_nc, n_char),
                {n_minorabsorbers});

        Array<std::string,1> identifier_minor(
                get_variable_string("identifier_minor", {n_minorabsorbers}, coef_nc, n_char), {n_minorabsorbers});

        Array<std::string,1> minor_gases_lower(
                get_variable_string("minor_gases_lower", {n_minor_absorber_intervals_lower}, coef_nc, n_char),
                {n_minor_absorber_intervals_lower});
        Array<std::string,1> minor_gases_upper(
                get_variable_string("minor_gases_upper", {n_minor_absorber_intervals_upper}, coef_nc, n_char),
                {n_minor_absorber_intervals_upper});

        Array<int,2> minor_limits_gpt_lower(
                coef_nc.get_variable<int>("minor_limits_gpt_lower", {n_minor_absorber_intervals_lower, n_pairs}),
                {n_pairs, n_minor_absorber_intervals_lower});
        Array<int,2> minor_limits_gpt_upper(
                coef_nc.get_variable<int>("minor_limits_gpt_upper", {n_minor_absorber_intervals_upper, n_pairs}),
                {n_pairs, n_minor_absorber_intervals_upper});

        Array<int,1> minor_scales_with_density_lower(
                coef_nc.get_variable<int>("minor_scales_with_density_lower", {n_minor_absorber_intervals_lower}),
                {n_minor_absorber_intervals_lower});
        Array<int,1> minor_scales_with_density_upper(
                coef_nc.get_variable<int>("minor_scales_with_density_upper", {n_minor_absorber_intervals_upper}),
                {n_minor_absorber_intervals_upper});

        Array<int,1> scale_by_complement_lower(
                coef_nc.get_variable<int>("scale_by_complement_lower", {n_minor_absorber_intervals_lower}),
                {n_minor_absorber_intervals_lower});
        Array<int,1> scale_by_complement_upper(
                coef_nc.get_variable<int>("scale_by_complement_upper", {n_minor_absorber_intervals_upper}),
                {n_minor_absorber_intervals_upper});

        Array<std::string,1> scaling_gas_lower(
                get_variable_string("scaling_gas_lower", {n_minor_absorber_intervals_lower}, coef_nc, n_char),
                {n_minor_absorber_intervals_lower});
        Array<std::string,1> scaling_gas_upper(
                get_variable_string("scaling_gas_upper", {n_minor_absorber_intervals_upper}, coef_nc, n_char),
                {n_minor_absorber_intervals_upper});

        Array<int,1> kminor_start_lower(
                coef_nc.get_variable<int>("kminor_start_lower", {n_minor_absorber_intervals_lower}),
                {n_minor_absorber_intervals_lower});
        Array<int,1> kminor_start_upper(
                coef_nc.get_variable<int>("kminor_start_upper", {n_minor_absorber_intervals_upper}),
                {n_minor_absorber_intervals_upper});

        Array<double,3> vmr_ref(
                coef_nc.get_variable<double>("vmr_ref", {n_temps, n_extabsorbers, n_layers}),
                {n_layers, n_extabsorbers, n_temps});

        Array<double,4> kmajor(
                coef_nc.get_variable<double>("kmajor", {n_temps, n_press+1, n_mixingfracs, n_gpts}),
                {n_gpts, n_mixingfracs, n_press+1, n_temps});

        // Keep the size at zero, if it does not exist.
        Array<double,3> rayl_lower;
        Array<double,3> rayl_upper;

        if (coef_nc.variable_exists("rayl_lower"))
        {
            rayl_lower.set_dims({n_gpts, n_mixingfracs, n_temps});
            rayl_upper.set_dims({n_gpts, n_mixingfracs, n_temps});
            rayl_lower = coef_nc.get_variable<double>("rayl_lower", {n_temps, n_mixingfracs, n_gpts});
            rayl_upper = coef_nc.get_variable<double>("rayl_upper", {n_temps, n_mixingfracs, n_gpts});
        }

        // Is it really LW if so read these variables as well.
        if (coef_nc.variable_exists("totplnk"))
        {
            int n_internal_sourcetemps = coef_nc.get_dimension_size("temperature_Planck");

            Array<double,2> totplnk(
                    coef_nc.get_variable<double>( "totplnk", {n_bnds, n_internal_sourcetemps}),
                    {n_internal_sourcetemps, n_bnds});
            Array<double,4> planck_frac(
                    coef_nc.get_variable<double>("plank_fraction", {n_temps, n_press+1, n_mixingfracs, n_gpts}),
                    {n_gpts, n_mixingfracs, n_press+1, n_temps});

            // Construct the k-distribution.
            return Gas_optics<double>(
                    gas_concs,
                    gas_names,
                    key_species,
                    band2gpt,
                    band_lims,
                    press_ref,
                    press_ref_trop,
                    temp_ref,
                    temp_ref_p,
                    temp_ref_t,
                    vmr_ref,
                    kmajor,
                    kminor_lower,
                    kminor_upper,
                    gas_minor,
                    identifier_minor,
                    minor_gases_lower,
                    minor_gases_upper,
                    minor_limits_gpt_lower,
                    minor_limits_gpt_upper,
                    minor_scales_with_density_lower,
                    minor_scales_with_density_upper,
                    scaling_gas_lower,
                    scaling_gas_upper,
                    scale_by_complement_lower,
                    scale_by_complement_upper,
                    kminor_start_lower,
                    kminor_start_upper,
                    totplnk,
                    planck_frac,
                    rayl_lower,
                    rayl_upper);
        }
        else
        {
            Array<double,1> solar_src(
                    coef_nc.get_variable<double>("solar_source", {n_gpts}), {n_gpts});

            return Gas_optics<double>(
                    gas_concs,
                    gas_names,
                    key_species,
                    band2gpt,
                    band_lims,
                    press_ref,
                    press_ref_trop,
                    temp_ref,
                    temp_ref_p,
                    temp_ref_t,
                    vmr_ref,
                    kmajor,
                    kminor_lower,
                    kminor_upper,
                    gas_minor,
                    identifier_minor,
                    minor_gases_lower,
                    minor_gases_upper,
                    minor_limits_gpt_lower,
                    minor_limits_gpt_upper,
                    minor_scales_with_density_lower,
                    minor_scales_with_density_upper,
                    scaling_gas_lower,
                    scaling_gas_upper,
                    scale_by_complement_lower,
                    scale_by_complement_upper,
                    kminor_start_lower,
                    kminor_start_upper,
                    solar_src,
                    rayl_lower,
                    rayl_upper);
        }
        // End reading of k-distribution.
    }
    template<typename TF>
    void solve_longwave(
            std::unique_ptr<Optical_props_arry<TF>>& optical_props,
            Array<TF,2>& flux_up,
            Array<TF,2>& flux_dn,
            Array<TF,2>& flux_net,
            const Gas_concs<TF>& gas_concs,
            const std::unique_ptr<Gas_optics<TF>>& kdist_lw,
            const std::unique_ptr<Source_func_lw<TF>>& sources,
            const Array<TF,2>& col_dry,
            const Array<TF,2>& p_lay, const Array<TF,2>& p_lev,
            const Array<TF,2>& t_lay, const Array<TF,2>& t_lev,
            const Array<TF,1>& t_sfc, const Array<TF,2>& emis_sfc,
            const int n_col, const int n_col_block, const int n_lay)
    {
        const int n_lev = n_lay + 1;
        const int n_blocks = n_col / n_col_block;
        const int n_col_block_left = n_col % n_col_block;

        // Store the number of bands in a variable.
        const int n_bnd = kdist_lw->get_nband();

        // Set the number of angles to 1.
        const int n_ang = 1;

        // Check the dimension ordering.
        const int top_at_1 = p_lay({1, 1}) < p_lay({1, n_lay});

        // Only subset if the amount of columns is larger than the actual block;
        if (n_col > n_col_block)
        {
            std::unique_ptr<Optical_props_arry<TF>> optical_props_subset =
                    std::make_unique<Optical_props_1scl<TF>>(n_col_block, n_lay, *kdist_lw);
            std::unique_ptr<Source_func_lw<TF>> sources_subset =
                    std::make_unique<Source_func_lw<TF>>(n_col_block, n_lay, *kdist_lw);

            std::unique_ptr<Optical_props_arry<TF>> optical_props_left =
                    std::make_unique<Optical_props_1scl<TF>>(n_col_block_left, n_lay, *kdist_lw);
            std::unique_ptr<Source_func_lw<TF>> sources_left =
                    std::make_unique<Source_func_lw<TF>>(n_col_block_left, n_lay, *kdist_lw);

            // Lambda function for solving optical properties subset.
            auto calc_optical_props_subset = [&](
                    const int col_s_in, const int col_e_in,
                    std::unique_ptr<Optical_props_arry<TF>>& optical_props_subset_in,
                    Source_func_lw<TF>& sources_subset_in)
            {
                const int n_col_in = col_e_in - col_s_in + 1;
                Gas_concs<TF> gas_concs_subset(gas_concs, col_s_in, n_col_in);

                kdist_lw->gas_optics(
                        p_lay.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}),
                        p_lev.subset({{ {col_s_in, col_e_in}, {1, n_lev} }}),
                        t_lay.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}),
                        t_sfc.subset({{ {col_s_in, col_e_in} }}),
                        gas_concs_subset,
                        optical_props_subset_in,
                        sources_subset_in,
                        col_dry.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}),
                        t_lev.subset  ({{ {col_s_in, col_e_in}, {1, n_lev} }}) );

                optical_props->set_subset(optical_props_subset_in, col_s_in, col_e_in);
                sources->set_subset(sources_subset_in, col_s_in, col_e_in);
            };

            // Lambda function for solving fluxes of a subset.
            auto calc_fluxes_subset = [&](
                    const int col_s_in, const int col_e_in,
                    const std::unique_ptr<Optical_props_arry<TF>>& optical_props_subset_in,
                    const Source_func_lw<TF>& sources_subset_in,
                    const Array<TF,2> emis_sfc_subset_in,
                    std::unique_ptr<Fluxes_broadband<TF>>& fluxes)
            {
                const int n_col_block_subset = col_e_in - col_s_in + 1;

                Rte_lw<TF>::rte_lw(
                        optical_props_subset_in,
                        top_at_1,
                        sources_subset_in,
                        emis_sfc_subset_in,
                        fluxes,
                        n_ang);

                // Copy the data to the output.
                for (int ilev=1; ilev<=n_lev; ++ilev)
                    for (int icol=1; icol<=n_col_block_subset; ++icol)
                    {
                        flux_up ({icol+col_s_in-1, ilev}) = fluxes->get_flux_up ()({icol, ilev});
                        flux_dn ({icol+col_s_in-1, ilev}) = fluxes->get_flux_dn ()({icol, ilev});
                        flux_net({icol+col_s_in-1, ilev}) = fluxes->get_flux_net()({icol, ilev});
                    }
            };

            for (int b=1; b<=n_blocks; ++b)
            {
                const int col_s = (b-1) * n_col_block + 1;
                const int col_e =  b    * n_col_block;

                calc_optical_props_subset(
                        col_s, col_e,
                        optical_props_subset,
                        *sources_subset);
            }

            if (n_col_block_left > 0)
            {
                const int col_s = n_col - n_col_block_left + 1;
                const int col_e = n_col;

                calc_optical_props_subset(
                        col_s, col_e,
                        optical_props_left,
                        *sources_left);
            }

            for (int b=1; b<=n_blocks; ++b)
            {
                const int col_s = (b-1) * n_col_block + 1;
                const int col_e =  b    * n_col_block;

                optical_props_subset->get_subset(optical_props, col_s, col_e);
                sources_subset->get_subset(*sources, col_s, col_e);

                Array<TF,2> emis_sfc_subset = emis_sfc.subset({{ {1, n_bnd}, {col_s, col_e} }});

                std::unique_ptr<Fluxes_broadband<TF>> fluxes_subset =
                        std::make_unique<Fluxes_byband<TF>>(n_col_block, n_lev, n_bnd);

                calc_fluxes_subset(
                        col_s, col_e,
                        optical_props_subset,
                        *sources_subset,
                        emis_sfc_subset,
                        fluxes_subset);
            }

            if (n_col_block_left > 0)
            {
                const int col_s = n_col - n_col_block_left + 1;
                const int col_e = n_col;

                optical_props_left->get_subset(optical_props, col_s, col_e);
                sources_left->get_subset(*sources, col_s, col_e);

                Array<TF,2> emis_sfc_left = emis_sfc.subset({{ {1, n_bnd}, {col_s, col_e} }});

                std::unique_ptr<Fluxes_broadband<TF>> fluxes_left =
                        std::make_unique<Fluxes_byband<TF>>(n_col_block_left, n_lev, n_bnd);

                calc_fluxes_subset(
                        col_s, col_e,
                        optical_props_left,
                        *sources_left,
                        emis_sfc_left,
                        fluxes_left);
            }
        }
        else
        {
            // Solve a single block, this does not require subsetting.
            kdist_lw->gas_optics(
                    p_lay,
                    p_lev,
                    t_lay,
                    t_sfc,
                    gas_concs,
                    optical_props,
                    *sources,
                    col_dry,
                    t_lev);

            std::unique_ptr<Fluxes_broadband<double>> fluxes =
                    std::make_unique<Fluxes_broadband<double>>(n_col, n_lev);

            const int top_at_1 = p_lay({1, 1}) < p_lay({1, n_lay});

            Rte_lw<double>::rte_lw(
                    optical_props,
                    top_at_1,
                    *sources,
                    emis_sfc,
                    fluxes,
                    n_ang);

            // Copy the data to the output.
            for (int ilev=1; ilev<=n_lev; ++ilev)
                for (int icol=1; icol<=n_col; ++icol)
                {
                    flux_up ({icol, ilev}) = fluxes->get_flux_up ()({icol, ilev});
                    flux_dn ({icol, ilev}) = fluxes->get_flux_dn ()({icol, ilev});
                    flux_net({icol, ilev}) = fluxes->get_flux_net()({icol, ilev});
                }
        }
    }

    template<typename TF>
    void solve_longwave_column(
            std::unique_ptr<Optical_props_arry<TF>>& optical_props,
            Array<TF,2>& flux_up, Array<TF,2>& flux_dn, Array<TF,2>& flux_net,
            Array<TF,2>& flux_dn_inc, const TF p_top,
            const Gas_concs<TF>& gas_concs,
            const std::unique_ptr<Gas_optics<TF>>& kdist_lw,
            const std::unique_ptr<Source_func_lw<TF>>& sources,
            const Array<TF,2>& col_dry,
            const Array<TF,2>& p_lay, const Array<TF,2>& p_lev,
            const Array<TF,2>& t_lay, const Array<TF,2>& t_lev,
            const Array<TF,1>& t_sfc, const Array<TF,2>& emis_sfc,
            const int n_lay)
    {
        const int n_col = 1;
        const int n_lev = n_lay + 1;

        // Set the number of angles to 1.
        const int n_ang = 1;

        // Check the dimension ordering.
        const int top_at_1 = p_lay({1, 1}) < p_lay({1, n_lay});

        // Solve a single block, this does not require subsetting.
        kdist_lw->gas_optics(
                p_lay,
                p_lev,
                t_lay,
                t_sfc,
                gas_concs,
                optical_props,
                *sources,
                col_dry,
                t_lev);

        std::unique_ptr<Fluxes_broadband<double>> fluxes =
                std::make_unique<Fluxes_broadband<double>>(n_col, n_lev);

        const int n_gpt = kdist_lw->get_ngpt();
        Array<double,3> gpt_flux_up({n_col, n_lev, n_gpt});
        Array<double,3> gpt_flux_dn({n_col, n_lev, n_gpt});

        Rte_lw<double>::rte_lw(
                optical_props,
                top_at_1,
                *sources,
                emis_sfc,
                gpt_flux_up,
                gpt_flux_dn,
                n_ang);

        fluxes->reduce(gpt_flux_up, gpt_flux_dn, optical_props, top_at_1);

        // Copy the data to the output.
        for (int ilev=1; ilev<=n_lev; ++ilev)
            for (int icol=1; icol<=n_col; ++icol)
            {
                flux_up ({icol, ilev}) = fluxes->get_flux_up ()({icol, ilev});
                flux_dn ({icol, ilev}) = fluxes->get_flux_dn ()({icol, ilev});
                flux_net({icol, ilev}) = fluxes->get_flux_net()({icol, ilev});
            }
    }

    template<typename TF>
    void solve_shortwave_column(
            std::unique_ptr<Optical_props_arry<TF>>& optical_props,
            Array<TF,2>& flux_up, Array<TF,2>& flux_dn,
            Array<TF,2>& flux_dn_dir, Array<TF,2>& flux_net,
            Array<TF,2>& flux_dn_inc, Array<TF,2>& flux_dn_dir_inc, const TF p_top,
            const Gas_concs<TF>& gas_concs,
            const Gas_optics<TF>& kdist_sw,
            const Array<TF,2>& col_dry,
            const Array<TF,2>& p_lay, const Array<TF,2>& p_lev,
            const Array<TF,2>& t_lay, const Array<TF,2>& t_lev,
            const Array<TF,1>& mu0,
            const Array<TF,2>& sfc_alb_dir, const Array<TF,2>& sfc_alb_dif,
            const TF tsi_scaling,
            const int n_lay)
    {
        const int n_col = 1;
        const int n_lev = n_lay + 1;

        // Check the dimension ordering.
        const int top_at_1 = p_lay({1, 1}) < p_lay({1, n_lay});

        // Create the field for the top of atmosphere source.
        const int n_gpt = kdist_sw.get_ngpt();
        Array<TF,2> toa_src({n_col, n_gpt});

        kdist_sw.gas_optics(
                p_lay,
                p_lev,
                t_lay,
                gas_concs,
                optical_props,
                toa_src,
                col_dry);

        if (tsi_scaling >= 0)
            for (int igpt=1; igpt<=n_gpt; ++igpt)
                toa_src({1, igpt}) *= tsi_scaling;

        std::unique_ptr<Fluxes_broadband<TF>> fluxes =
                std::make_unique<Fluxes_broadband<TF>>(n_col, n_lev);

        Array<double,3> gpt_flux_up    ({n_col, n_lev, n_gpt});
        Array<double,3> gpt_flux_dn    ({n_col, n_lev, n_gpt});
        Array<double,3> gpt_flux_dn_dir({n_col, n_lev, n_gpt});

        Rte_sw<TF>::rte_sw(
                optical_props,
                top_at_1,
                mu0,
                toa_src,
                sfc_alb_dir,
                sfc_alb_dif,
                gpt_flux_up,
                gpt_flux_dn,
                gpt_flux_dn_dir);

        fluxes->reduce(
                gpt_flux_up, gpt_flux_dn, gpt_flux_dn_dir,
                optical_props, top_at_1);

        // Copy the data to the output.
        for (int ilev=1; ilev<=n_lev; ++ilev)
            for (int icol=1; icol<=n_col; ++icol)
            {
                flux_up    ({icol, ilev}) = fluxes->get_flux_up    ()({icol, ilev});
                flux_dn    ({icol, ilev}) = fluxes->get_flux_dn    ()({icol, ilev});
                flux_dn_dir({icol, ilev}) = fluxes->get_flux_dn_dir()({icol, ilev});
                flux_net   ({icol, ilev}) = fluxes->get_flux_net   ()({icol, ilev});
            }
    }

    template<typename TF>
    void solve_shortwave(
            std::unique_ptr<Optical_props_arry<TF>>& optical_props,
            Array<TF,2>& flux_up, Array<TF,2>& flux_dn,
            Array<TF,2>& flux_dn_dir, Array<TF,2>& flux_net,
            const Gas_concs<TF>& gas_concs,
            const std::unique_ptr<Gas_optics<TF>>& kdist_sw,
            const Array<TF,2>& col_dry,
            const Array<TF,2>& p_lay, const Array<TF,2>& p_lev,
            const Array<TF,2>& t_lay, const Array<TF,2>& t_lev,
            const Array<TF,1>& mu0,
            const Array<TF,2>& sfc_alb_dir, const Array<TF,2>& sfc_alb_dif,
            const TF tsi_scaling,
            const int n_col, const int n_col_block, const int n_lay)
    {
        const int n_lev = n_lay + 1;
        const int n_blocks = n_col / n_col_block;
        const int n_col_block_left = n_col % n_col_block;

        // Check the dimension ordering.
        const int top_at_1 = p_lay({1, 1}) < p_lay({1, n_lay});

        // Create the field for the top of atmosphere source.
        const int n_gpt = kdist_sw->get_ngpt();
        Array<TF,2> toa_src({n_col, n_gpt});

        // Store the number of bands in a variable.
        const int n_bnd = kdist_sw->get_nband();

        std::unique_ptr<Optical_props_arry<TF>> optical_props_subset =
                std::make_unique<Optical_props_2str<TF>>(n_col_block, n_lay, *kdist_sw);

        std::unique_ptr<Optical_props_arry<TF>> optical_props_left =
                std::make_unique<Optical_props_2str<TF>>(n_col_block_left, n_lay, *kdist_sw);

        // Lambda function for solving optical properties subset.
        auto calc_optical_props_subset = [&](
                const int col_s_in, const int col_e_in,
                std::unique_ptr<Optical_props_arry<TF>>& optical_props_subset_in)
        {
            const int n_col_in = col_e_in - col_s_in + 1;

            Gas_concs<TF> gas_concs_subset(gas_concs, col_s_in, n_col_in);
            Array<TF,2> toa_src_subset({n_col_in, n_gpt});

            kdist_sw->gas_optics(
                    p_lay.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}),
                    p_lev.subset({{ {col_s_in, col_e_in}, {1, n_lev} }}),
                    t_lay.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}),
                    gas_concs_subset,
                    optical_props_subset_in,
                    toa_src_subset,
                    col_dry.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}) );

            optical_props->set_subset(optical_props_subset_in, col_s_in, col_e_in);

            // Copy the data to the output.
            for (int igpt=1; igpt<=n_gpt; ++igpt)
                for (int icol=1; icol<=n_col_in; ++icol)
                    toa_src({icol+col_s_in-1, igpt}) = toa_src_subset({icol, igpt});
        };

        // Lambda function for solving fluxes of a subset.
        auto calc_fluxes_subset = [&](
                const int col_s_in, const int col_e_in,
                const std::unique_ptr<Optical_props_arry<TF>>& optical_props_subset_in,
                const Array<TF,1>& mu0_subset_in,
                const Array<TF,2>& toa_src_subset_in,
                const Array<TF,2>& sfc_alb_dir_subset_in,
                const Array<TF,2>& sfc_alb_dif_subset_in,
                std::unique_ptr<Fluxes_broadband<TF>>& fluxes)
        {
            const int n_col_block_subset = col_e_in - col_s_in + 1;

            Rte_sw<TF>::rte_sw(
                    optical_props_subset_in,
                    top_at_1,
                    mu0_subset_in,
                    toa_src_subset_in,
                    sfc_alb_dir_subset_in,
                    sfc_alb_dif_subset_in,
                    fluxes);

            // Copy the data to the output.
            for (int ilev=1; ilev<=n_lev; ++ilev)
                for (int icol=1; icol<=n_col_block_subset; ++icol)
                {
                    flux_up    ({icol+col_s_in-1, ilev}) = fluxes->get_flux_up    ()({icol, ilev});
                    flux_dn    ({icol+col_s_in-1, ilev}) = fluxes->get_flux_dn    ()({icol, ilev});
                    flux_dn_dir({icol+col_s_in-1, ilev}) = fluxes->get_flux_dn_dir()({icol, ilev});
                    flux_net   ({icol+col_s_in-1, ilev}) = fluxes->get_flux_net   ()({icol, ilev});
                }
        };

        for (int b=1; b<=n_blocks; ++b)
        {
            const int col_s = (b-1) * n_col_block + 1;
            const int col_e =  b    * n_col_block;

            calc_optical_props_subset(
                    col_s, col_e,
                    optical_props_subset);
        }

        if (n_col_block_left > 0)
        {
            const int col_s = n_col - n_col_block_left + 1;
            const int col_e = n_col;

            calc_optical_props_subset(
                    col_s, col_e,
                    optical_props_left);
        }

        for (int b=1; b<=n_blocks; ++b)
        {
            const int col_s = (b-1) * n_col_block + 1;
            const int col_e =  b    * n_col_block;

            optical_props_subset->get_subset(optical_props, col_s, col_e);

            Array<TF,1> mu0_subset = mu0.subset({{ {col_s, col_e} }});
            Array<TF,2> toa_src_subset = toa_src.subset({{ {col_s, col_e}, {1, n_gpt} }});
            Array<TF,2> sfc_alb_dir_subset = sfc_alb_dir.subset({{ {1, n_bnd}, {col_s, col_e} }});
            Array<TF,2> sfc_alb_dif_subset = sfc_alb_dif.subset({{ {1, n_bnd}, {col_s, col_e} }});

            std::unique_ptr<Fluxes_broadband<TF>> fluxes_subset =
                    std::make_unique<Fluxes_byband<TF>>(n_col_block, n_lev, n_bnd);

            calc_fluxes_subset(
                    col_s, col_e,
                    optical_props_subset,
                    mu0_subset,
                    toa_src_subset,
                    sfc_alb_dir_subset,
                    sfc_alb_dif_subset,
                    fluxes_subset);
        }

        if (n_col_block_left > 0)
        {
            const int col_s = n_col - n_col_block_left + 1;
            const int col_e = n_col;

            optical_props_left->get_subset(optical_props, col_s, col_e);

            Array<TF,1> mu0_left = mu0.subset({{ {col_s, col_e} }});
            Array<TF,2> toa_src_left = toa_src.subset({{ {col_s, col_e}, {1, n_gpt} }});
            Array<TF,2> sfc_alb_dir_left = sfc_alb_dir.subset({{ {1, n_bnd}, {col_s, col_e} }});
            Array<TF,2> sfc_alb_dif_left = sfc_alb_dif.subset({{ {1, n_bnd}, {col_s, col_e} }});

            std::unique_ptr<Fluxes_broadband<TF>> fluxes_left =
                    std::make_unique<Fluxes_broadband<TF>>(n_col_block_left, n_lev);

            calc_fluxes_subset(
                    col_s, col_e,
                    optical_props_left,
                    mu0_left,
                    toa_src_left,
                    sfc_alb_dir_left,
                    sfc_alb_dif_left,
                    fluxes_left);
        }
    }
}

template<typename TF>
Radiation_rrtmgp<TF>::Radiation_rrtmgp(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
	Radiation<TF>(masterin, gridin, fieldsin, inputin)
{
	swradiation = "rrtmgp";
}

template<typename TF>
void Radiation_rrtmgp<TF>::init()
{
}

template<typename TF>
void Radiation_rrtmgp<TF>::create(
        Input& input, Netcdf_handle& input_nc, Thermo<TF>& thermo,
        Stats<TF>& stats, Column<TF>& column, Cross<TF>& cross, Dump<TF>& dump)
{
    auto& gd = grid.get_grid_data();

    // 1. Load the available gas concentrations from the group of the netcdf file.
    load_gas_concs<double>(gas_concs, input_nc);

    // 2. Set up the gas optics classes for long and shortwave.
    kdist_lw = std::make_unique<Gas_optics<double>>(
            load_and_init_gas_optics(master, gas_concs, "coefficients_lw.nc"));
    kdist_sw = std::make_unique<Gas_optics<double>>(
            load_and_init_gas_optics(master, gas_concs, "coefficients_sw.nc"));

    // 3. Read the atmospheric pressure and temperature.
    const int n_col = 1;

    Netcdf_group rad_nc = input_nc.get_group("radiation");

    const int n_lay = rad_nc.get_dimension_size("lay");
    const int n_lev = rad_nc.get_dimension_size("lev");

    Array<double,2> p_lay(rad_nc.get_variable<double>("p_lay", {n_lay, n_col}), {n_col, n_lay});
    Array<double,2> t_lay(rad_nc.get_variable<double>("t_lay", {n_lay, n_col}), {n_col, n_lay});
    Array<double,2> p_lev(rad_nc.get_variable<double>("p_lev", {n_lev, n_col}), {n_col, n_lev});
    Array<double,2> t_lev(rad_nc.get_variable<double>("t_lev", {n_lev, n_col}), {n_col, n_lev});

    Array<double,2> col_dry({n_col, n_lay});
    if (input_nc.variable_exists("col_dry"))
        col_dry = rad_nc.get_variable<double>("col_dry", {n_lay, n_col});
    else
    {
        kdist_lw->get_col_dry(col_dry, gas_concs.get_vmr("h2o"), p_lev);
        kdist_sw->get_col_dry(col_dry, gas_concs.get_vmr("h2o"), p_lev);
    }

    const double tsi_scaling = input.get_item<double>("radiation", "tsi_scaling", "", -999.);

    // 4. Read the boundary conditions for the longwave and shortwave solver.
    // Set the surface temperature and emissivity.
    // CvH: communicate with surface scheme.
    const double t_sfc_in = input.get_item<double>("radiation", "t_sfc", "");
    Array<double,1> t_sfc({1});
    t_sfc({1}) = t_sfc_in;

    const double emis_sfc_in = input.get_item<double>("radiation", "emis_sfc", "");
    const int n_bnd = kdist_lw->get_nband();
    Array<double,2> emis_sfc({n_bnd, 1});
    for (int ibnd=1; ibnd<=n_bnd; ++ibnd)
        emis_sfc({ibnd, 1}) = emis_sfc_in;

    // Set the solar zenith angle and albedo.
    Array<double,1> sza({n_col});
    Array<double,2> sfc_alb_dir({n_bnd, n_col});
    Array<double,2> sfc_alb_dif({n_bnd, n_col});

    const double sza_in = input.get_item<double>("radiation", "sza", "");
    sza({1}) = sza_in;

    const double sfc_alb_dir_in = input.get_item<double>("radiation", "sfc_alb_dir", "");
    const double sfc_alb_dif_in = input.get_item<double>("radiation", "sfc_alb_dif", "");
    for (int ibnd=1; ibnd<=n_bnd; ++ibnd)
    {
        sfc_alb_dir({ibnd, 1}) = sfc_alb_dir_in;
        sfc_alb_dif({ibnd, 1}) = sfc_alb_dif_in;
    }

    Array<double,1> mu0({n_col});
    mu0({1}) = std::cos(sza({1}));

    // Compute the longwave for the reference profile.
    std::unique_ptr<Source_func_lw<double>> sources_lw =
            std::make_unique<Source_func_lw<double>>(n_col, n_lay, *kdist_lw);

    std::unique_ptr<Optical_props_arry<double>> optical_props_lw =
            std::make_unique<Optical_props_1scl<double>>(n_col, n_lay, *kdist_lw);

    Array<double,2> lw_flux_up ({n_col, n_lev});
    Array<double,2> lw_flux_dn ({n_col, n_lev});
    Array<double,2> lw_flux_net({n_col, n_lev});

    const int n_gpt = kdist_sw->get_ngpt();
    Array<double,2> lw_flux_dn_inc({n_col, n_gpt});

    solve_longwave_column<double>(
            optical_props_lw,
            lw_flux_up, lw_flux_dn, lw_flux_net,
            lw_flux_dn_inc, thermo.get_ph_vector()[gd.kstart],
            gas_concs,
            kdist_lw,
            sources_lw,
            col_dry,
            p_lay, p_lev,
            t_lay, t_lev,
            t_sfc, emis_sfc,
            n_lay);

    std::unique_ptr<Optical_props_arry<double>> optical_props_sw =
            std::make_unique<Optical_props_2str<double>>(n_col, n_lay, *kdist_sw);

    Array<double,2> sw_flux_up    ({n_col, n_lev});
    Array<double,2> sw_flux_dn    ({n_col, n_lev});
    Array<double,2> sw_flux_dn_dir({n_col, n_lev});
    Array<double,2> sw_flux_net   ({n_col, n_lev});

    Array<double,2> sw_flux_dn_inc    ({n_col, n_gpt});
    Array<double,2> sw_flux_dn_dir_inc({n_col, n_gpt});

    solve_shortwave_column<double>(
            optical_props_sw,
            sw_flux_up, sw_flux_dn, sw_flux_dn_dir, sw_flux_net,
            sw_flux_dn_inc, sw_flux_dn_dir_inc, thermo.get_ph_vector()[gd.kstart],
            gas_concs,
            *kdist_sw,
            col_dry,
            p_lay, p_lev,
            t_lay, t_lev,
            mu0,
            sfc_alb_dir, sfc_alb_dif,
            tsi_scaling,
            n_lay);

    // Save the reference radiation fluxes in the stats.
}

template<typename TF>
void Radiation_rrtmgp<TF>::exec(
        Thermo<TF>& thermo, const double time, Timeloop<TF>& timeloop, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();
}

template class Radiation_rrtmgp<double>;
template class Radiation_rrtmgp<float>;
