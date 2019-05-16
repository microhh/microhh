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

    /*
    // READ THE ATMOSPHERIC DATA.
    int n_lay = input_nc.get_dimension_size("lay");
    int n_lev = input_nc.get_dimension_size("lev");
    int n_col = input_nc.get_dimension_size("col");

    Array<TF,2> p_lay(input_nc.get_variable<TF>("p_lay", {n_lay, n_col}), {n_col, n_lay});
    Array<TF,2> t_lay(input_nc.get_variable<TF>("t_lay", {n_lay, n_col}), {n_col, n_lay});
    Array<TF,2> p_lev(input_nc.get_variable<TF>("p_lev", {n_lev, n_col}), {n_col, n_lev});
    Array<TF,2> t_lev(input_nc.get_variable<TF>("t_lev", {n_lev, n_col}), {n_col, n_lev});

    const int top_at_1 = p_lay({1, 1}) < p_lay({1, n_lay});

    gas_concs.set_vmr("h2o",
            Array<TF,2>(input_nc.get_variable<TF>("vmr_h2o", {n_lay, n_col}), {n_col, n_lay}));
    gas_concs.set_vmr("co2",
            Array<TF,2>(input_nc.get_variable<TF>("vmr_co2", {n_lay, n_col}), {n_col, n_lay}));
    gas_concs.set_vmr("o3",
            Array<TF,2>(input_nc.get_variable<TF>("vmr_o3", {n_lay, n_col}), {n_col, n_lay}));
    gas_concs.set_vmr("n2o",
            Array<TF,2>(input_nc.get_variable<TF>("vmr_n2o", {n_lay, n_col}), {n_col, n_lay}));
    gas_concs.set_vmr("co",
            Array<TF,2>(input_nc.get_variable<TF>("vmr_co", {n_lay, n_col}), {n_col, n_lay}));
    gas_concs.set_vmr("ch4",
            Array<TF,2>(input_nc.get_variable<TF>("vmr_ch4", {n_lay, n_col}), {n_col, n_lay}));
    gas_concs.set_vmr("o2",
            Array<TF,2>(input_nc.get_variable<TF>("vmr_o2", {n_lay, n_col}), {n_col, n_lay}));
    gas_concs.set_vmr("n2",
            Array<TF,2>(input_nc.get_variable<TF>("vmr_n2", {n_lay, n_col}), {n_col, n_lay}));

    // CvH: does this one need to be present?
    Array<TF,2> col_dry(input_nc.get_variable<TF>("col_dry", {n_lay, n_col}), {n_col, n_lay});
    */

    // Construct the gas optics class, always in double precision.
    kdist_lw = std::make_unique<Gas_optics<double>>(
            load_and_init_gas_optics(master, gas_concs, "coefficients_lw.nc"));
    kdist_sw = std::make_unique<Gas_optics<double>>(
            load_and_init_gas_optics(master, gas_concs, "coefficients_sw.nc"));

    const int n_bnd = kdist_lw->get_nband();
}

template<typename TF>
void Radiation_rrtmgp<TF>::exec(
        Thermo<TF>& thermo, const double time, Timeloop<TF>& timeloop, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    const int n_col = gd.imax;

    ////// SOLVING THE OPTICAL PROPERTIES FOR LONGWAVE RADIATION //////
    const int n_col_block = 4;

    // Download surface boundary conditions for long wave.
    Array<double,2> emis_sfc_tmp(
            input_nc.get_variable<double>(
                    "emis_sfc", {n_col, n_bnd}), {n_bnd, n_col});
    Array<double,1> t_sfc_tmp(
            input_nc.get_variable<double>(
                    "t_sfc", {n_col}), {n_col});

    emis_sfc = emis_sfc_tmp;
    t_sfc = t_sfc_tmp;

    // Read the sources and create containers for the substeps.
    int n_blocks = n_col / n_col_block;
    int n_col_block_left = n_col % n_col_block;

    optical_props        = std::make_unique<Optical_props_1scl<double>>(n_col      , n_lay, kdist);
    optical_props_subset = std::make_unique<Optical_props_1scl<double>>(n_col_block, n_lay, kdist);

    Source_func_lw<double> sources       (n_col      , n_lay, kdist);
    Source_func_lw<double> sources_subset(n_col_block, n_lay, kdist);

    auto calc_optical_props_subset = [&](
            const int col_s_in, const int col_e_in,
            std::unique_ptr<Optical_props_arry<double>>& optical_props_subset_in,
            Source_func_lw<double>& sources_subset_in)
    {
        const int n_col_in = col_e_in - col_s_in + 1;
        Gas_concs<double> gas_concs_subset(gas_concs, col_s_in, n_col_in);

        kdist.gas_optics(
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
        sources.set_subset(sources_subset_in, col_s_in, col_e_in);
    };

    for (int b=1; b<=n_blocks; ++b)
    {
        const int col_s = (b-1) * n_col_block + 1;
        const int col_e =  b    * n_col_block;

        calc_optical_props_subset(
                col_s, col_e,
                optical_props_subset,
                sources_subset);
    }

    if (n_col_block_left > 0)
    {
        optical_props_left = std::make_unique<Optical_props_1scl<double>>(n_col_block_left, n_lay, kdist);
        Source_func_lw<double> sources_left(n_col_block_left, n_lay, kdist);

        const int col_s = n_col - n_col_block_left + 1;
        const int col_e = n_col;

        calc_optical_props_subset(
                col_s, col_e,
                optical_props_left,
                sources_left);
    }
    
    ////// SOLVING THE FLUXES FOR LONGWAVE RADIATION //////
    master.print_message("STEP 2: Computing the longwave radiation fluxes.\n");

    const int n_ang = input_nc.get_variable<double>("angle");

    Array<double,2> flux_up ({n_col, n_lev});
    Array<double,2> flux_dn ({n_col, n_lev});
    Array<double,2> flux_net({n_col, n_lev});

    auto calc_fluxes_subset = [&](
            const int col_s_in, const int col_e_in,
            const std::unique_ptr<Optical_props_arry<double>>& optical_props_subset_in,
            const Source_func_lw<double>& sources_subset_in,
            const Array<double,2> emis_sfc_subset_in,
            std::unique_ptr<Fluxes_broadband<double>>& fluxes)
    {
        const int n_col_block_subset = col_e_in - col_s_in + 1;

        // CvH: I removed the pointer assignments of the fluxes, as this is unportable Fortran code.
        Rte_lw<double>::rte_lw(
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

        optical_props_subset->get_subset(optical_props, col_s, col_e);
        sources_subset.get_subset(sources, col_s, col_e);

        Array<double,2> emis_sfc_subset = emis_sfc.subset({{ {1, n_bnd}, {col_s, col_e} }});

        std::unique_ptr<Fluxes_broadband<double>> fluxes_subset =
                std::make_unique<Fluxes_broadband<double>>(n_col_block, n_lev, n_bnd);

        calc_fluxes_subset(
                col_s, col_e,
                optical_props_subset,
                sources_subset,
                emis_sfc_subset,
                fluxes_subset);
    }

    if (n_col_block_left > 0)
    {
        const int col_s = n_col - n_col_block_left + 1;
        const int col_e = n_col;

        // CvH, check for reuse of this field.
        Source_func_lw<double> sources_left(n_col_block_left, n_lay, kdist);

        optical_props_left->get_subset(optical_props, col_s, col_e);
        sources_left.get_subset(sources, col_s, col_e);

        Array<double,2> emis_sfc_left = emis_sfc.subset({{ {1, n_bnd}, {col_s, col_e} }});

        std::unique_ptr<Fluxes_broadband<double>> fluxes_left =
                std::make_unique<Fluxes_broadband<double>>(n_col_block_left, n_lev, n_bnd);

        calc_fluxes_subset(
                col_s, col_e,
                optical_props_left,
                sources_left,
                emis_sfc_left,
                fluxes_left);
    }
}

template class Radiation_rrtmgp<double>;
template class Radiation_rrtmgp<float>;
