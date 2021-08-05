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

#include <boost/algorithm/string.hpp>
#include <numeric>
#include <string>
#include <cmath>

#include "radiation_rrtmgp.h"
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "input.h"
#include "netcdf_interface.h"
#include "stats.h"
#include "cross.h"
#include "column.h"
#include "constants.h"
#include "timeloop.h"

// RRTMGP headers.
#include "Array.h"
#include "Optical_props.h"
#include "Gas_optics_rrtmgp.h"
#include "Gas_concs.h"
#include "Fluxes.h"
#include "Rte_lw.h"
#include "Rte_sw.h"
#include "Source_functions.h"
#include "Cloud_optics.h"

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
    void load_gas_concs(
            Gas_concs<TF>& gas_concs, Netcdf_handle& input_nc, const std::string& dim_name)
    {
        const int n_lay = input_nc.get_dimension_size(dim_name);

        const std::vector<std::string> possible_gases = {
                "h2o", "co2" ,"o3", "n2o", "co", "ch4", "o2", "n2",
                "ccl4", "cfc11", "cfc12", "cfc22",
                "hfc143a", "hfc125", "hfc23", "hfc32", "hfc134a",
                "cf4", "no2" };

        for (const std::string& gas_name : possible_gases)
        {
            if (input_nc.variable_exists(gas_name))
            {
                std::map<std::string, int> dims = input_nc.get_variable_dimensions(gas_name);
                const int n_dims = dims.size();

                if (n_dims == 0)
                {
                    gas_concs.set_vmr(gas_name, input_nc.get_variable<TF>(gas_name));
                }
                else if (n_dims == 1)
                {
                    if (dims.at(dim_name) == n_lay)
                        gas_concs.set_vmr(gas_name,
                                Array<TF,1>(input_nc.get_variable<TF>(gas_name, {n_lay}), {n_lay}));
                    else
                        throw std::runtime_error("Illegal dimensions of gas \"" + gas_name + "\" in input");
                }
                else
                {
                    throw std::runtime_error("Illegal dimensions of gas \"" + gas_name + "\" in input");
                }
            }
        };
    }

    template<typename TF>
    Gas_optics_rrtmgp<TF> load_and_init_gas_optics(
            Master& master,
            const Gas_concs<TF>& gas_concs,
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
        Array<TF,2> band_lims(coef_nc.get_variable<TF>("bnd_limits_wavenumber", {n_bnds, 2}), {2, n_bnds});
        Array<int,2> band2gpt(coef_nc.get_variable<int>("bnd_limits_gpt", {n_bnds, 2}), {2, n_bnds});
        Array<TF,1> press_ref(coef_nc.get_variable<TF>("press_ref", {n_press}), {n_press});
        Array<TF,1> temp_ref(coef_nc.get_variable<TF>("temp_ref", {n_temps}), {n_temps});

        TF temp_ref_p = coef_nc.get_variable<TF>("absorption_coefficient_ref_P");
        TF temp_ref_t = coef_nc.get_variable<TF>("absorption_coefficient_ref_T");
        TF press_ref_trop = coef_nc.get_variable<TF>("press_ref_trop");

        Array<TF,3> kminor_lower(
                coef_nc.get_variable<TF>("kminor_lower", {n_temps, n_mixingfracs, n_contributors_lower}),
                {n_contributors_lower, n_mixingfracs, n_temps});
        Array<TF,3> kminor_upper(
                coef_nc.get_variable<TF>("kminor_upper", {n_temps, n_mixingfracs, n_contributors_upper}),
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

        Array<BOOL_TYPE,1> minor_scales_with_density_lower(
                coef_nc.get_variable<BOOL_TYPE>("minor_scales_with_density_lower", {n_minor_absorber_intervals_lower}),
                {n_minor_absorber_intervals_lower});
        Array<BOOL_TYPE,1> minor_scales_with_density_upper(
                coef_nc.get_variable<BOOL_TYPE>("minor_scales_with_density_upper", {n_minor_absorber_intervals_upper}),
                {n_minor_absorber_intervals_upper});

        Array<BOOL_TYPE,1> scale_by_complement_lower(
                coef_nc.get_variable<BOOL_TYPE>("scale_by_complement_lower", {n_minor_absorber_intervals_lower}),
                {n_minor_absorber_intervals_lower});
        Array<BOOL_TYPE,1> scale_by_complement_upper(
                coef_nc.get_variable<BOOL_TYPE>("scale_by_complement_upper", {n_minor_absorber_intervals_upper}),
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

        Array<TF,3> vmr_ref(
                coef_nc.get_variable<TF>("vmr_ref", {n_temps, n_extabsorbers, n_layers}),
                {n_layers, n_extabsorbers, n_temps});

        Array<TF,4> kmajor(
                coef_nc.get_variable<TF>("kmajor", {n_temps, n_press+1, n_mixingfracs, n_gpts}),
                {n_gpts, n_mixingfracs, n_press+1, n_temps});

        // Keep the size at zero, if it does not exist.
        Array<TF,3> rayl_lower;
        Array<TF,3> rayl_upper;

        if (coef_nc.variable_exists("rayl_lower"))
        {
            rayl_lower.set_dims({n_gpts, n_mixingfracs, n_temps});
            rayl_upper.set_dims({n_gpts, n_mixingfracs, n_temps});
            rayl_lower = coef_nc.get_variable<TF>("rayl_lower", {n_temps, n_mixingfracs, n_gpts});
            rayl_upper = coef_nc.get_variable<TF>("rayl_upper", {n_temps, n_mixingfracs, n_gpts});
        }

        // Is it really LW if so read these variables as well.
        if (coef_nc.variable_exists("totplnk"))
        {
            int n_internal_sourcetemps = coef_nc.get_dimension_size("temperature_Planck");

            Array<TF,2> totplnk(
                    coef_nc.get_variable<TF>( "totplnk", {n_bnds, n_internal_sourcetemps}),
                    {n_internal_sourcetemps, n_bnds});
            Array<TF,4> planck_frac(
                    coef_nc.get_variable<TF>("plank_fraction", {n_temps, n_press+1, n_mixingfracs, n_gpts}),
                    {n_gpts, n_mixingfracs, n_press+1, n_temps});

            // Construct the k-distribution.
            return Gas_optics_rrtmgp<TF>(
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
            Array<TF,1> solar_src_quiet(
                    coef_nc.get_variable<TF>("solar_source_quiet", {n_gpts}), {n_gpts});
            Array<TF,1> solar_src_facular(
                    coef_nc.get_variable<TF>("solar_source_facular", {n_gpts}), {n_gpts});
            Array<TF,1> solar_src_sunspot(
                    coef_nc.get_variable<TF>("solar_source_sunspot", {n_gpts}), {n_gpts});

            TF tsi = coef_nc.get_variable<TF>("tsi_default");
            TF mg_index = coef_nc.get_variable<TF>("mg_default");
            TF sb_index = coef_nc.get_variable<TF>("sb_default");

            return Gas_optics_rrtmgp<TF>(
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
                    solar_src_quiet,
                    solar_src_facular,
                    solar_src_sunspot,
                    tsi,
                    mg_index,
                    sb_index,
                    rayl_lower,
                    rayl_upper);
        }
        // End reading of k-distribution.
    }

    template<typename TF>
    Cloud_optics<TF> load_and_init_cloud_optics(
            Master& master,
            const std::string& coef_file)
    {
        // READ THE COEFFICIENTS FOR THE OPTICAL SOLVER.
        Netcdf_file coef_nc(master, coef_file, Netcdf_mode::Read);

        // Read look-up table coefficient dimensions
        int n_band     = coef_nc.get_dimension_size("nband");
        int n_rghice   = coef_nc.get_dimension_size("nrghice");
        int n_size_liq = coef_nc.get_dimension_size("nsize_liq");
        int n_size_ice = coef_nc.get_dimension_size("nsize_ice");

        Array<TF,2> band_lims_wvn(coef_nc.get_variable<TF>("bnd_limits_wavenumber", {n_band, 2}), {2, n_band});

        // Read look-up table constants.
        TF radliq_lwr = coef_nc.get_variable<TF>("radliq_lwr");
        TF radliq_upr = coef_nc.get_variable<TF>("radliq_upr");
        TF radliq_fac = coef_nc.get_variable<TF>("radliq_fac");

        TF radice_lwr = coef_nc.get_variable<TF>("radice_lwr");
        TF radice_upr = coef_nc.get_variable<TF>("radice_upr");
        TF radice_fac = coef_nc.get_variable<TF>("radice_fac");

        Array<TF,2> lut_extliq(
                coef_nc.get_variable<TF>("lut_extliq", {n_band, n_size_liq}), {n_size_liq, n_band});
        Array<TF,2> lut_ssaliq(
                coef_nc.get_variable<TF>("lut_ssaliq", {n_band, n_size_liq}), {n_size_liq, n_band});
        Array<TF,2> lut_asyliq(
                coef_nc.get_variable<TF>("lut_asyliq", {n_band, n_size_liq}), {n_size_liq, n_band});

        Array<TF,3> lut_extice(
                coef_nc.get_variable<TF>("lut_extice", {n_rghice, n_band, n_size_ice}), {n_size_ice, n_band, n_rghice});
        Array<TF,3> lut_ssaice(
                coef_nc.get_variable<TF>("lut_ssaice", {n_rghice, n_band, n_size_ice}), {n_size_ice, n_band, n_rghice});
        Array<TF,3> lut_asyice(
                coef_nc.get_variable<TF>("lut_asyice", {n_rghice, n_band, n_size_ice}), {n_size_ice, n_band, n_rghice});

        return Cloud_optics<TF>(
                band_lims_wvn,
                radliq_lwr, radliq_upr, radliq_fac,
                radice_lwr, radice_upr, radice_fac,
                lut_extliq, lut_ssaliq, lut_asyliq,
                lut_extice, lut_ssaice, lut_asyice);
    }

    template<typename TF>
    void calc_tendency(
            TF* restrict thlt_rad,
            const TF* restrict flux_up, const TF* restrict flux_dn,
            const TF* restrict rho, const TF* exner, const TF* dz,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int igc, const int jgc, const int kgc,
            const int jj, const int kk,
            const int jj_nogc, const int kk_nogc)
    {
        for (int k=kstart; k<kend; ++k)
        {
            // Conversion from energy to temperature.
            const TF fac = TF(1.) / (rho[k]*Constants::cp<TF>*exner[k]*dz[k]);

            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ijk_nogc = (i-igc) + (j-jgc)*jj_nogc + (k-kgc)*kk_nogc;

                    thlt_rad[ijk] -= fac *
                        ( flux_up[ijk_nogc+kk_nogc] - flux_up[ijk_nogc]
                        - flux_dn[ijk_nogc+kk_nogc] + flux_dn[ijk_nogc] );
                }
        }
    }

    template<typename TF>
    void add_tendency(
            TF* restrict thlt, const TF* restrict thlt_rad,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    thlt[ijk] += thlt_rad[ijk];
                }
    }

    template<typename TF>
    void store_surface_fluxes(
            TF* restrict flux_up_sfc, TF* restrict flux_dn_sfc,
            const TF* restrict flux_up, const TF* restrict flux_dn,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int igc, const int jgc,
            const int jj, const int kk,
            const int jj_nogc)
    {
        for (int j=jstart; j<jend; ++j)
            for (int i=istart; i<iend; ++i)
            {
                const int ij = i + j*jj;
                const int ijk_nogc = (i-igc) + (j-jgc)*jj_nogc;

                flux_up_sfc[ij] = flux_up[ijk_nogc];
                flux_dn_sfc[ij] = flux_dn[ijk_nogc];
            }
    }

    template<typename TF>
    void solve_longwave_column(
            std::unique_ptr<Optical_props_arry<TF>>& optical_props,
            Array<TF,2>& flux_up, Array<TF,2>& flux_dn, Array<TF,2>& flux_net,
            Array<TF,2>& flux_dn_inc, const TF p_top,
            const Gas_concs<TF>& gas_concs,
            const std::unique_ptr<Gas_optics_rrtmgp<TF>>& kdist_lw,
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

        std::unique_ptr<Fluxes_broadband<TF>> fluxes =
                std::make_unique<Fluxes_broadband<TF>>(n_col, n_lev);

        const int n_gpt = kdist_lw->get_ngpt();
        Array<TF,3> gpt_flux_up({n_col, n_lev, n_gpt});
        Array<TF,3> gpt_flux_dn({n_col, n_lev, n_gpt});

        Rte_lw<TF>::rte_lw(
                optical_props,
                top_at_1,
                *sources,
                emis_sfc,
                Array<TF,2>(),
                gpt_flux_up,
                gpt_flux_dn,
                n_ang);

        fluxes->reduce(gpt_flux_up, gpt_flux_dn, optical_props, top_at_1);

        // Find the index where p_lev exceeds p_top.
        int idx_top=1;
        for (; idx_top<=n_lev; ++idx_top)
        {
            if (p_lev({1, idx_top}) < p_top)
                break;
        }

        // Calculate the interpolation factors.
        const int idx_bot = idx_top - 1;
        const TF fac_bot = (p_top - p_lev({1, idx_top})) / (p_lev({1, idx_bot}) - p_lev({1, idx_top}));
        const TF fac_top = 1. - fac_bot;

        // Interpolate the top boundary conditions.
        for (int igpt=1; igpt<=n_gpt; ++igpt)
            flux_dn_inc({1, igpt}) = fac_bot * gpt_flux_dn({1, idx_bot, igpt}) + fac_top * gpt_flux_dn({1, idx_top, igpt});

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
            Array<TF,2>& flux_dn_dir_inc, Array<TF,2>& flux_dn_dif_inc, const TF p_top,
            const Gas_concs<TF>& gas_concs,
            const Gas_optics_rrtmgp<TF>& kdist_sw,
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

        Array<TF,3> gpt_flux_up    ({n_col, n_lev, n_gpt});
        Array<TF,3> gpt_flux_dn    ({n_col, n_lev, n_gpt});
        Array<TF,3> gpt_flux_dn_dir({n_col, n_lev, n_gpt});

        Rte_sw<TF>::rte_sw(
                optical_props,
                top_at_1,
                mu0,
                toa_src,
                sfc_alb_dir,
                sfc_alb_dif,
                Array<TF,2>(),
                gpt_flux_up,
                gpt_flux_dn,
                gpt_flux_dn_dir);

        fluxes->reduce(
                gpt_flux_up, gpt_flux_dn, gpt_flux_dn_dir,
                optical_props, top_at_1);

        // Find the index where p_lev exceeds p_top.
        int idx_top=1;
        for (; idx_top<=n_lev; ++idx_top)
        {
            if (p_lev({1, idx_top}) < p_top)
                break;
        }

        // Calculate the interpolation factors.
        const int idx_bot = idx_top - 1;
        const TF fac_bot = (p_top - p_lev({1, idx_top})) / (p_lev({1, idx_bot}) - p_lev({1, idx_top}));
        const TF fac_top = 1. - fac_bot;

        // Interpolate the top boundary conditions.
        for (int igpt=1; igpt<=n_gpt; ++igpt)
        {
            const TF flux_dn_tot = fac_bot * gpt_flux_dn    ({1, idx_bot, igpt}) + fac_top * gpt_flux_dn    ({1, idx_top, igpt});
            const TF flux_dn_dir = fac_bot * gpt_flux_dn_dir({1, idx_bot, igpt}) + fac_top * gpt_flux_dn_dir({1, idx_top, igpt});
            // Divide out the cosine of the solar zenith angle.
            flux_dn_dir_inc({1, igpt}) = flux_dn_dir / mu0({1});
            flux_dn_dif_inc({1, igpt}) = flux_dn_tot - flux_dn_dir;
        }

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
    TF deg_to_rad(const TF deg)
    {
        return TF(2.*M_PI/360. * deg);
    }

    template<typename TF>
    TF calc_cos_zenith_angle(
            const TF lat, const TF lon, const int day_of_year,
            const TF seconds_since_midnight, const int year)
    {
        /* Based on: Paltridge, G. W. and Platt, C. M. R. (1976).
                     Radiative Processes in Meteorology and Climatology.
                     Elsevier, New York, 318 pp. */

        const TF pi = TF(M_PI);

        // Account for leap year
        int days_per_year;
        if ((year%4 == 0) && ((year%100 != 0) || (year%400 == 0)))
            days_per_year = 366;
        else
            days_per_year = 365;

        // DOY in time calculations are zero based:
        const int doy = day_of_year-1;

        // Lat/lon in radians
        const TF radlat = lat * pi/TF(180);
        const TF radlon = lon * pi/TF(180);

        // DOY in range (0,2*pi)
        const TF doy_pi = TF(2)*pi*doy/days_per_year;

        // Solar declination angle
        const TF declination_angle = \
                TF(0.006918) - TF(0.399912) * std::cos(doy_pi) + TF(0.070257) * std::sin(doy_pi)
              - TF(0.006758) * std::cos(TF(2)*doy_pi) + TF(0.000907) * std::sin(TF(2)*doy_pi)
              - TF(0.002697) * std::cos(TF(3)*doy_pi) + TF(0.00148)  * std::sin(TF(3)*doy_pi);

        // Hour angle in radians, using true solar time
        const TF a1 = (TF(1.00554) * doy - TF( 6.28306)) * pi/TF(180);
        const TF a2 = (TF(1.93946) * doy + TF(23.35089)) * pi/TF(180);
        const TF a3 = (TF(7.67825) * std::sin(a1) + TF(10.09176) * std::sin(a2)) / TF(60);

        const TF hour_solar_time = (seconds_since_midnight/TF(3600)) - a3 + radlon * (TF(180.)/pi/TF(15.0));
        const TF hour_angle = (hour_solar_time-TF(12))*TF(15.0)*(pi/TF(180));

        // Cosine of solar zenith angle
        const TF cos_zenith = std::sin(radlat) * std::sin(declination_angle) +
                              std::cos(radlat) * std::cos(declination_angle) * std::cos(hour_angle);

        return cos_zenith;
    }

    template<typename TF>
    TF calc_sun_distance_factor(const TF frac_doy)
    {
        // Based on: An Introduction to Atmospheric Radiation, Liou, Eq. 2.2.9.
        constexpr TF an [] = {1.000110, 0.034221, 0.000719};
        constexpr TF bn [] = {0,        0.001280, 0.000077};

        const TF pi = TF(M_PI);

        // NOTE: DOY in time calculations are zero based
        const TF t = TF(2)*pi*(frac_doy-TF(1))/TF(365);

        TF factor = TF(0);
        for (int n=0; n<3; ++n)
            factor += an[n]*std::cos(TF(n)*t) + bn[n]*std::sin(TF(n)*t);

        return factor;
    }
}

template<typename TF>
Radiation_rrtmgp<TF>::Radiation_rrtmgp(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
        Radiation<TF>(masterin, gridin, fieldsin, inputin)
{
    swradiation = "rrtmgp";

    sw_longwave  = inputin.get_item<bool>("radiation", "swlongwave" , "", true);
    sw_shortwave = inputin.get_item<bool>("radiation", "swshortwave", "", true);
    sw_fixed_sza = inputin.get_item<bool>("radiation", "swfixedsza", "", true);

    sw_clear_sky_stats = inputin.get_item<bool>("radiation", "swclearskystats", "", false);

    dt_rad = inputin.get_item<double>("radiation", "dt_rad", "");

    t_sfc       = inputin.get_item<TF>("radiation", "t_sfc"      , "");
    emis_sfc    = inputin.get_item<TF>("radiation", "emis_sfc"   , "");
    sfc_alb_dir = inputin.get_item<TF>("radiation", "sfc_alb_dir", "");
    sfc_alb_dif = inputin.get_item<TF>("radiation", "sfc_alb_dif", "");
    tsi_scaling = inputin.get_item<TF>("radiation", "tsi_scaling", "", -999.);

    if (sw_fixed_sza)
    {
        const TF sza = inputin.get_item<TF>("radiation", "sza", "");
        mu0 = std::cos(sza);
    }
    else
    {
        lat = inputin.get_item<TF>("radiation", "lat", "");
        lon = inputin.get_item<TF>("radiation", "lon", "");
    }

    auto& gd = grid.get_grid_data();
    fields.init_diagnostic_field("thlt_rad", "Tendency by radiation", "K s-1", "radiation", gd.sloc);
}

template<typename TF>
void Radiation_rrtmgp<TF>::init(Timeloop<TF>& timeloop)
{
    auto& gd = grid.get_grid_data();

    idt_rad = static_cast<unsigned long>(timeloop.get_ifactor() * dt_rad + 0.5);

    // Check if restarttime is dividable by dt_rad
    if (timeloop.get_isavetime() % idt_rad != 0)
        throw std::runtime_error("Restart \"savetime\" is not an (integer) multiple of \"dt_rad\"");

    // Resize surface radiation fields
    lw_flux_dn_sfc.resize(gd.ijcells);
    lw_flux_up_sfc.resize(gd.ijcells);

    sw_flux_dn_sfc.resize(gd.ijcells);
    sw_flux_up_sfc.resize(gd.ijcells);
}

template<typename TF>
unsigned long Radiation_rrtmgp<TF>::get_time_limit(unsigned long itime)
{
    unsigned long idtlim = idt_rad - itime % idt_rad;
    return idtlim;
}

template<typename TF>
void Radiation_rrtmgp<TF>::create(
        Input& input, Netcdf_handle& input_nc, Thermo<TF>& thermo,
        Stats<TF>& stats, Column<TF>& column, Cross<TF>& cross, Dump<TF>& dump)
{
    // Check if the thermo supports the radiation.
    if (thermo.get_switch() != "moist")
    {
        const std::string error = "Radiation does not support thermo mode " + thermo.get_switch();
        throw std::runtime_error(error);
    }

    // Initialize the tendency if the radiation is used.
    if (stats.get_switch() && (sw_longwave || sw_shortwave))
        stats.add_tendency(*fields.st.at("thl"), "z", tend_name, tend_longname);

    // Create the gas optics solver that is needed for the column and model solver.
    // The solver is in a try catch, because in theory it could crash on a single core,
    // which requires an MPI_Abort in order to prevent a deadlock.
    try
    {
        create_solver(input, input_nc, thermo, stats, column);

        // Solve the reference column to compute upper boundary conditions.
        create_column(input, input_nc, thermo, stats);
    }
    catch (std::exception& e)
    {
        #ifdef USEMPI
        std::cout << "SINGLE PROCESS EXCEPTION: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        #else
        throw;
        #endif
    }

    if (stats.get_switch() && sw_shortwave)
    {
        const std::string group_name = "radiation";
        stats.add_time_series("sza", "solar zenith angle", "rad", group_name);
        stats.add_time_series("sw_flux_dn_toa", "shortwave downwelling flux at toa", "W m-2", group_name);
    }

    // Get the allowed cross sections from the cross list
    std::vector<std::string> allowed_crossvars_radiation;

    if (sw_shortwave)
    {
        allowed_crossvars_radiation.push_back("sw_flux_up");
        allowed_crossvars_radiation.push_back("sw_flux_dn");
        allowed_crossvars_radiation.push_back("sw_flux_dn_dir");

        if (sw_clear_sky_stats)
        {
            allowed_crossvars_radiation.push_back("sw_flux_up_clear");
            allowed_crossvars_radiation.push_back("sw_flux_dn_clear");
            allowed_crossvars_radiation.push_back("sw_flux_dn_dir_clear");
        }
    }

    if (sw_longwave)
    {
        allowed_crossvars_radiation.push_back("lw_flux_up");
        allowed_crossvars_radiation.push_back("lw_flux_dn");

        if (sw_clear_sky_stats)
        {
            allowed_crossvars_radiation.push_back("lw_flux_up_clear");
            allowed_crossvars_radiation.push_back("lw_flux_dn_clear");
        }
    }

    crosslist = cross.get_enabled_variables(allowed_crossvars_radiation);
}

template<typename TF>
void Radiation_rrtmgp<TF>::create_column(
        Input& input, Netcdf_handle& input_nc, Thermo<TF>& thermo, Stats<TF>& stats)
{
    // 1. Load the available gas concentrations from the group of the netcdf file.
    Netcdf_handle& rad_nc = input_nc.get_group("radiation");

    load_gas_concs<TF>(gas_concs_col, rad_nc, "lay");

    // 2. Set the coordinate for the reference profiles in the stats, before calling the other creates.
    if (stats.get_switch() && (sw_longwave || sw_shortwave))
    {
        const int n_col = 1;
        const int n_lev = rad_nc.get_dimension_size("lev");
        Array<TF,2> p_lev(rad_nc.get_variable<TF>("p_lev", {n_lev, n_col}), {n_col, n_lev});

        stats.add_dimension("p_rad", n_lev);

        const std::string group_name = "radiation";
        const std::string root_group= "";

        stats.add_fixed_prof_raw(
                "p_rad",
                "Pressure of radiation reference column",
                "Pa", "p_rad", root_group,
                p_lev.v());
    }

    // 3. Read background profiles on pressure levels
    read_background_profiles(rad_nc, gas_concs_col);

    // 4. Call the column solvers for longwave and shortwave.
    if (sw_longwave)
        create_column_longwave (input, rad_nc, thermo, stats, gas_concs_col);
    if (sw_shortwave)
        create_column_shortwave(input, rad_nc, thermo, stats, gas_concs_col);
}

template<typename TF>
void Radiation_rrtmgp<TF>::read_background_profiles(
        Netcdf_handle& rad_nc,
        const Gas_concs<TF>& gas_concs_col)
{
    // Read the atmospheric pressure and temperature.
    n_col = 1;
    n_lay_col = rad_nc.get_dimension_size("lay");
    n_lev_col = rad_nc.get_dimension_size("lev");

    p_lay_col.set_dims({n_col, n_lay_col});
    t_lay_col.set_dims({n_col, n_lay_col});
    p_lev_col.set_dims({n_col, n_lev_col});
    t_lev_col.set_dims({n_col, n_lev_col});
    col_dry.  set_dims({n_col, n_lay_col});

    p_lay_col = rad_nc.get_variable<TF>("p_lay", {n_lay_col, n_col});
    t_lay_col = rad_nc.get_variable<TF>("t_lay", {n_lay_col, n_col});
    p_lev_col = rad_nc.get_variable<TF>("p_lev", {n_lev_col, n_col});
    t_lev_col = rad_nc.get_variable<TF>("t_lev", {n_lev_col, n_col});

    if (rad_nc.variable_exists("col_dry"))
        col_dry = rad_nc.get_variable<TF>("col_dry", {n_lay_col, n_col});
    else
        Gas_optics_rrtmgp<TF>::get_col_dry(col_dry, gas_concs_col.get_vmr("h2o"), p_lev_col);
}

template<typename TF>
void Radiation_rrtmgp<TF>::create_column_longwave(
        Input& input, Netcdf_handle& rad_nc, Thermo<TF>& thermo, Stats<TF>& stats,
        const Gas_concs<TF>& gas_concs_col)
{
    auto& gd = grid.get_grid_data();

    // Read the boundary conditions.
    // Set the surface temperature and emissivity.
    Array<TF,1> t_sfc({1});
    t_sfc({1}) = t_lev_col({1,1});

    const int n_bnd = kdist_lw->get_nband();
    Array<TF,2> emis_sfc({n_bnd, 1});
    for (int ibnd=1; ibnd<=n_bnd; ++ibnd)
        emis_sfc({ibnd, 1}) = this->emis_sfc;

    // Compute the longwave for the reference profile.
    sources_lw = std::make_unique<Source_func_lw<TF>>(n_col, n_lay_col, *kdist_lw);
    optical_props_lw = std::make_unique<Optical_props_1scl<TF>>(n_col, n_lay_col, *kdist_lw);

    lw_flux_up_col .set_dims({n_col, n_lev_col});
    lw_flux_dn_col .set_dims({n_col, n_lev_col});
    lw_flux_net_col.set_dims({n_col, n_lev_col});

    const int n_gpt = kdist_lw->get_ngpt();
    lw_flux_dn_inc.set_dims({n_col, n_gpt});

    solve_longwave_column<TF>(
            optical_props_lw,
            lw_flux_up_col, lw_flux_dn_col, lw_flux_net_col,
            lw_flux_dn_inc, thermo.get_basestate_vector("ph")[gd.kend],
            gas_concs_col,
            kdist_lw,
            sources_lw,
            col_dry,
            p_lay_col, p_lev_col,
            t_lay_col, t_lev_col,
            t_sfc, emis_sfc,
            n_lay_col);

    // Save the reference profile fluxes in the stats.
    if (stats.get_switch())
    {
        const std::string group_name = "radiation";

        stats.add_fixed_prof_raw(
                "lw_flux_up_ref",
                "Longwave upwelling flux of reference column",
                "W m-2", "p_rad", group_name,
                lw_flux_up_col.v());

        stats.add_fixed_prof_raw(
                "lw_flux_dn_ref",
                "Longwave downwelling flux of reference column",
                "W m-2", "p_rad", group_name,
                lw_flux_dn_col.v());
    }
}

template<typename TF>
void Radiation_rrtmgp<TF>::create_column_shortwave(
        Input& input, Netcdf_handle& rad_nc, Thermo<TF>& thermo, Stats<TF>& stats,
        const Gas_concs<TF>& gas_concs_col)
{
    auto& gd = grid.get_grid_data();

    // Read the boundary conditions.
    const int n_bnd = kdist_sw->get_nband();

    optical_props_sw = std::make_unique<Optical_props_2str<TF>>(n_col, n_lay_col, *kdist_sw);

    sw_flux_up_col    .set_dims({n_col, n_lev_col});
    sw_flux_dn_col    .set_dims({n_col, n_lev_col});
    sw_flux_dn_dir_col.set_dims({n_col, n_lev_col});
    sw_flux_net_col   .set_dims({n_col, n_lev_col});

    const int n_gpt = kdist_sw->get_ngpt();
    sw_flux_dn_dir_inc.set_dims({n_col, n_gpt});
    sw_flux_dn_dif_inc.set_dims({n_col, n_gpt});

    if (sw_fixed_sza)
    {
        // Set the solar zenith angle and albedo.
        Array<TF,2> sfc_alb_dir({n_bnd, n_col});
        Array<TF,2> sfc_alb_dif({n_bnd, n_col});

        for (int ibnd=1; ibnd<=n_bnd; ++ibnd)
        {
            sfc_alb_dir({ibnd, 1}) = this->sfc_alb_dir;
            sfc_alb_dif({ibnd, 1}) = this->sfc_alb_dif;
        }

        Array<TF,1> mu0({n_col});
        mu0({1}) = std::max(this->mu0, Constants::mu0_min<TF>);

        solve_shortwave_column<TF>(
                optical_props_sw,
                sw_flux_up_col, sw_flux_dn_col, sw_flux_dn_dir_col, sw_flux_net_col,
                sw_flux_dn_dir_inc, sw_flux_dn_dif_inc, thermo.get_basestate_vector("ph")[gd.kend],
                gas_concs_col,
                *kdist_sw,
                col_dry,
                p_lay_col, p_lev_col,
                t_lay_col, t_lev_col,
                mu0,
                sfc_alb_dir, sfc_alb_dif,
                tsi_scaling,
                n_lay_col);

        // Save the reference profile fluxes in the stats.
        if (stats.get_switch())
        {
            const std::string group_name = "radiation";

            stats.add_fixed_prof_raw(
                    "sw_flux_up_ref",
                    "Shortwave upwelling flux of reference column",
                    "W m-2", "p_rad", group_name,
                    sw_flux_up_col.v());
            stats.add_fixed_prof_raw(
                    "sw_flux_dn_ref",
                    "Shortwave downwelling flux of reference column",
                    "W m-2", "p_rad", group_name,
                    sw_flux_dn_col.v());
            stats.add_fixed_prof_raw(
                    "sw_flux_dn_dir_ref",
                    "Shortwave direct downwelling flux of reference column",
                    "W m-2", "p_rad", group_name,
                    sw_flux_dn_dir_col.v());
        }
    }
}

template<typename TF>
void Radiation_rrtmgp<TF>::create_solver(
        Input& input, Netcdf_handle& input_nc, Thermo<TF>& thermo,
        Stats<TF>& stats, Column<TF>& column)
{
    // 1. Load the available gas concentrations from the group of the netcdf file.
    Netcdf_handle& rad_input_nc = input_nc.get_group("init");
    load_gas_concs<TF>(gas_concs, rad_input_nc, "z");

    // 2. Pass the gas concentrations to the solver initializers.
    if (sw_longwave)
        create_solver_longwave(input, input_nc, thermo, stats, column, gas_concs);
    if (sw_shortwave)
        create_solver_shortwave(input, input_nc, thermo, stats, column, gas_concs);
}

template<typename TF>
void Radiation_rrtmgp<TF>::create_solver_longwave(
        Input& input, Netcdf_handle& input_nc, Thermo<TF>& thermo,
        Stats<TF>& stats, Column<TF>& column,
        const Gas_concs<TF>& gas_concs)
{
    const std::string group_name = "radiation";

    // Set up the gas optics classes for long and shortwave.
    kdist_lw = std::make_unique<Gas_optics_rrtmgp<TF>>(
            load_and_init_gas_optics<TF>(master, gas_concs, "coefficients_lw.nc"));

    cloud_lw = std::make_unique<Cloud_optics<TF>>(
            load_and_init_cloud_optics<TF>(master, "cloud_coefficients_lw.nc"));

    // Set up the statistics.
    if (stats.get_switch())
    {
        stats.add_prof("lw_flux_up", "Longwave upwelling flux"  , "W m-2", "zh", group_name);
        stats.add_prof("lw_flux_dn", "Longwave downwelling flux", "W m-2", "zh", group_name);

        if (sw_clear_sky_stats)
        {
            stats.add_prof("lw_flux_up_clear", "Clear-sky longwave upwelling flux"  , "W m-2", "zh", group_name);
            stats.add_prof("lw_flux_dn_clear", "Clear-sky longwave downwelling flux", "W m-2", "zh", group_name);
        }
    }

    // Set up the column statistics
    if (column.get_switch())
    {
        column.add_prof("lw_flux_up", "Longwave upwelling flux"  , "W m-2", "zh");
        column.add_prof("lw_flux_dn", "Longwave downwelling flux", "W m-2", "zh");

        if (sw_clear_sky_stats)
        {
            column.add_prof("lw_flux_up_clear", "Clear-sky longwave upwelling flux"  , "W m-2", "zh");
            column.add_prof("lw_flux_dn_clear", "Clear-sky longwave downwelling flux", "W m-2", "zh");
        }
    }
}

template<typename TF>
void Radiation_rrtmgp<TF>::create_solver_shortwave(
        Input& input, Netcdf_handle& input_nc, Thermo<TF>& thermo,
        Stats<TF>& stats, Column<TF>& column,
        const Gas_concs<TF>& gas_concs)
{
    const std::string group_name = "radiation";

    // Set up the gas optics classes for long and shortwave.
    kdist_sw = std::make_unique<Gas_optics_rrtmgp<TF>>(
            load_and_init_gas_optics<TF>(master, gas_concs, "coefficients_sw.nc"));

    cloud_sw = std::make_unique<Cloud_optics<TF>>(
            load_and_init_cloud_optics<TF>(master, "cloud_coefficients_sw.nc"));

    // Set up the statistics.
    if (stats.get_switch())
    {
        stats.add_prof("sw_flux_up"    , "Shortwave upwelling flux"         , "W m-2", "zh", group_name);
        stats.add_prof("sw_flux_dn"    , "Shortwave downwelling flux"       , "W m-2", "zh", group_name);
        stats.add_prof("sw_flux_dn_dir", "Shortwave direct downwelling flux", "W m-2", "zh", group_name);

        if (sw_clear_sky_stats)
        {
            stats.add_prof("sw_flux_up_clear"    , "Clear-sky shortwave upwelling flux"         , "W m-2", "zh", group_name);
            stats.add_prof("sw_flux_dn_clear"    , "Clear-sky shortwave downwelling flux"       , "W m-2", "zh", group_name);
            stats.add_prof("sw_flux_dn_dir_clear", "Clear-sky shortwave direct downwelling flux", "W m-2", "zh", group_name);
        }
    }

    // Set up the column statistics
    if (column.get_switch())
    {
        column.add_prof("sw_flux_up"    , "Shortwave upwelling flux"         , "W m-2", "zh");
        column.add_prof("sw_flux_dn"    , "Shortwave downwelling flux"       , "W m-2", "zh");
        column.add_prof("sw_flux_dn_dir", "Shortwave direct downwelling flux", "W m-2", "zh");

        if (sw_clear_sky_stats)
        {
            column.add_prof("sw_flux_up_clear"    , "Clear-sky shortwave upwelling flux"         , "W m-2", "zh");
            column.add_prof("sw_flux_dn_clear"    , "Clear-sky shortwave downwelling flux"       , "W m-2", "zh");
            column.add_prof("sw_flux_dn_dir_clear", "Clear-sky shortwave direct downwelling flux", "W m-2", "zh");
        }
    }
}

#ifndef USECUDA
template<typename TF>
void Radiation_rrtmgp<TF>::exec(
        Thermo<TF>& thermo, const double time, Timeloop<TF>& timeloop, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    const bool do_radiation = ((timeloop.get_itime() % idt_rad == 0) && !timeloop.in_substep()) ;

    if (do_radiation)
    {
        // Set the tendency to zero.
        std::fill(fields.sd.at("thlt_rad")->fld.begin(), fields.sd.at("thlt_rad")->fld.end(), TF(0.));

        auto t_lay = fields.get_tmp();
        auto t_lev = fields.get_tmp();
        auto h2o   = fields.get_tmp(); // This is the volume mixing ratio, not the specific humidity of vapor.
        auto clwp  = fields.get_tmp();
        auto ciwp  = fields.get_tmp();

        // Set the input to the radiation on a 3D grid without ghost cells.
        thermo.get_radiation_fields(*t_lay, *t_lev, *h2o, *clwp, *ciwp);

        const int nmaxh = gd.imax*gd.jmax*(gd.ktot+1);
        const int ijmax = gd.imax*gd.jmax;

        Array<TF,2> t_lay_a(t_lay->fld, {gd.imax*gd.jmax, gd.ktot});
        Array<TF,2> t_lev_a(t_lev->fld, {gd.imax*gd.jmax, gd.ktot+1});
        Array<TF,1> t_sfc_a(t_lev->fld_bot, {gd.imax*gd.jmax});
        Array<TF,2> h2o_a(h2o->fld, {gd.imax*gd.jmax, gd.ktot});
        Array<TF,2> clwp_a(clwp->fld, {gd.imax*gd.jmax, gd.ktot});
        Array<TF,2> ciwp_a(ciwp->fld, {gd.imax*gd.jmax, gd.ktot});

        Array<TF,2> flux_up ({gd.imax*gd.jmax, gd.ktot+1});
        Array<TF,2> flux_dn ({gd.imax*gd.jmax, gd.ktot+1});
        Array<TF,2> flux_net({gd.imax*gd.jmax, gd.ktot+1});

        const bool compute_clouds = true;

        try
        {
            if (sw_longwave)
            {
                exec_longwave(
                        thermo, timeloop, stats,
                        flux_up, flux_dn, flux_net,
                        t_lay_a, t_lev_a, t_sfc_a, h2o_a, clwp_a, ciwp_a,
                        compute_clouds);

                calc_tendency(
                        fields.sd.at("thlt_rad")->fld.data(),
                        flux_up.ptr(), flux_dn.ptr(),
                        fields.rhoref.data(), thermo.get_basestate_vector("exner").data(),
                        gd.dz.data(),
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                        gd.igc, gd.jgc, gd.kgc,
                        gd.icells, gd.ijcells,
                        gd.imax, gd.imax*gd.jmax);

                store_surface_fluxes(
                        lw_flux_up_sfc.data(), lw_flux_dn_sfc.data(),
                        flux_up.ptr(), flux_dn.ptr(),
                        gd.istart, gd.iend,
                        gd.jstart, gd.jend,
                        gd.igc, gd.jgc,
                        gd.icells, gd.ijcells,
                        gd.imax);
            }

            if (sw_shortwave)
            {
                if (!sw_fixed_sza)
                {
                    // Update the solar zenith angle, and calculate new shortwave reference column
                    const int day_of_year = int(timeloop.calc_day_of_year());
                    const int year = timeloop.get_year();
                    const TF seconds_after_midnight = TF(timeloop.calc_hour_of_day()*3600);
                    this->mu0 = calc_cos_zenith_angle(lat, lon, day_of_year, seconds_after_midnight, year);

                    // Calculate correction factor for impact Sun's distance on the solar "constant"
                    const TF frac_day_of_year = TF(day_of_year) + seconds_after_midnight / TF(86400);
                    this->tsi_scaling = calc_sun_distance_factor(frac_day_of_year);

                    if (is_day(this->mu0))
                    {
                        const int n_bnd = kdist_sw->get_nband();

                        // Set the solar zenith angle and albedo.
                        Array<TF,2> sfc_alb_dir({n_bnd, n_col});
                        Array<TF,2> sfc_alb_dif({n_bnd, n_col});

                        for (int ibnd=1; ibnd<=n_bnd; ++ibnd)
                        {
                            sfc_alb_dir({ibnd, 1}) = this->sfc_alb_dir;
                            sfc_alb_dif({ibnd, 1}) = this->sfc_alb_dif;
                        }

                        Array<TF,1> mu0({n_col});
                        mu0({1}) = this->mu0;

                        solve_shortwave_column<TF>(
                                optical_props_sw,
                                sw_flux_up_col, sw_flux_dn_col, sw_flux_dn_dir_col, sw_flux_net_col,
                                sw_flux_dn_dir_inc, sw_flux_dn_dif_inc, thermo.get_basestate_vector("ph")[gd.kend],
                                gas_concs_col,
                                *kdist_sw,
                                col_dry,
                                p_lay_col, p_lev_col,
                                t_lay_col, t_lev_col,
                                mu0,
                                sfc_alb_dir, sfc_alb_dif,
                                tsi_scaling,
                                n_lay_col);
                    }
                }

                if (is_day(this->mu0))
                {
                    Array<TF,2> flux_dn_dir({gd.imax*gd.jmax, gd.ktot+1});

                    exec_shortwave(
                            thermo, timeloop, stats,
                            flux_up, flux_dn, flux_dn_dir, flux_net,
                            t_lay_a, t_lev_a, h2o_a, clwp_a, ciwp_a,
                            compute_clouds);

                    calc_tendency(
                            fields.sd.at("thlt_rad")->fld.data(),
                            flux_up.ptr(), flux_dn.ptr(),
                            fields.rhoref.data(), thermo.get_basestate_vector("exner").data(),
                            gd.dz.data(),
                            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                            gd.igc, gd.jgc, gd.kgc,
                            gd.icells, gd.ijcells,
                            gd.imax, gd.imax*gd.jmax);

                    store_surface_fluxes(
                            sw_flux_up_sfc.data(), sw_flux_dn_sfc.data(),
                            flux_up.ptr(), flux_dn.ptr(),
                            gd.istart, gd.iend,
                            gd.jstart, gd.jend,
                            gd.igc, gd.jgc,
                            gd.icells, gd.ijcells,
                            gd.imax);
                }
                else
                {
                    // Set the surface fluxes to zero, for (e.g.) the land-surface model.
                    std::fill(sw_flux_up_sfc.begin(), sw_flux_up_sfc.end(), TF(0));
                    std::fill(sw_flux_dn_sfc.begin(), sw_flux_dn_sfc.end(), TF(0));
                }
            }
        } // End try block.
        catch (std::exception& e)
        {
            #ifdef USEMPI
            std::cout << "SINGLE PROCESS EXCEPTION: " << e.what() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            #else
            throw;
            #endif
        }

        fields.release_tmp(t_lay);
        fields.release_tmp(t_lev);
        fields.release_tmp(h2o);
        fields.release_tmp(clwp);
        fields.release_tmp(ciwp);
    }

    // Always add the tendency.
    add_tendency(
            fields.st.at("thl")->fld.data(),
            fields.sd.at("thlt_rad")->fld.data(),
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    stats.calc_tend(*fields.st.at("thl"), tend_name);
}


template<typename TF>
std::vector<TF>& Radiation_rrtmgp<TF>::get_surface_radiation(std::string name)
{
    if (name == "sw_down")
        return sw_flux_dn_sfc;
    else if (name == "sw_up")
        return sw_flux_up_sfc;
    else if (name == "lw_down")
        return lw_flux_dn_sfc;
    else if (name == "lw_up")
        return lw_flux_up_sfc;
    else
    {
        std::string error = "Variable \"" + name + "\" is not a valid surface radiation field";
        throw std::runtime_error(error);
    }
}
#endif

namespace
{
    template<typename TF>
    void add_ghost_cells(
            TF* restrict out, const TF* restrict in,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend_field,
            const int igc, const int jgc, const int kgc,
            const int jj, const int kk,
            const int jj_nogc, const int kk_nogc)
    {
        // Value of kend_field is either kend or kend+1.
        #pragma omp parallel for
        for (int k=kstart; k<kend_field; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int ijk_nogc = (i-igc) + (j-jgc)*jj_nogc + (k-kgc)*kk_nogc;
                    out[ijk] = in[ijk_nogc];
                }
    }
}

template<typename TF>
void Radiation_rrtmgp<TF>::exec_all_stats(
        Stats<TF>& stats, Cross<TF>& cross,
        Dump<TF>& dump, Column<TF>& column,
        Thermo<TF>& thermo, Timeloop<TF>& timeloop,
        const unsigned long itime, const int iotime)
{
    const bool do_stats  = stats.do_statistics(itime);
    const bool do_cross  = cross.do_cross(itime) && crosslist.size() > 0;
    const bool do_column = column.do_column(itime);

    // Return in case of no stats or cross section.
    if ( !(do_stats || do_cross || do_column) )
        return;

    const TF no_offset = 0.;
    const TF no_threshold = 0.;

    // CvH: lots of code repetition with exec()
    auto& gd = grid.get_grid_data();

    auto t_lay = fields.get_tmp();
    auto t_lev = fields.get_tmp();
    auto h2o   = fields.get_tmp(); // This is the volume mixing ratio, not the specific humidity of vapor.
    auto clwp  = fields.get_tmp();
    auto ciwp  = fields.get_tmp();

    // Set the input to the radiation on a 3D grid without ghost cells.
    thermo.get_radiation_fields(*t_lay, *t_lev, *h2o, *clwp, *ciwp);

    const int nmaxh = gd.imax*gd.jmax*(gd.ktot+1);
    const int ijmax = gd.imax*gd.jmax;

    Array<TF,2> t_lay_a(t_lay->fld, {gd.imax*gd.jmax, gd.ktot});
    Array<TF,2> t_lev_a(t_lev->fld, {gd.imax*gd.jmax, gd.ktot+1});
    Array<TF,1> t_sfc_a(t_lev->fld_bot, {gd.imax*gd.jmax});
    Array<TF,2> h2o_a(h2o->fld, {gd.imax*gd.jmax, gd.ktot});
    Array<TF,2> clwp_a(clwp->fld, {gd.imax*gd.jmax, gd.ktot});
    Array<TF,2> ciwp_a(ciwp->fld, {gd.imax*gd.jmax, gd.ktot});

    fields.release_tmp(t_lay);
    fields.release_tmp(t_lev);
    fields.release_tmp(h2o);
    fields.release_tmp(clwp);
    fields.release_tmp(ciwp);

    Array<TF,2> flux_up ({gd.imax*gd.jmax, gd.ktot+1});
    Array<TF,2> flux_dn ({gd.imax*gd.jmax, gd.ktot+1});
    Array<TF,2> flux_net({gd.imax*gd.jmax, gd.ktot+1});

    auto tmp = fields.get_tmp();
    tmp->loc = gd.wloc;

    const bool compute_clouds = true;

    // Use a lambda function to avoid code repetition.
    auto save_stats_and_cross = [&](
            const Array<TF,2>& array, const std::string& name, const std::array<int,3>& loc)
    {
        if (do_stats || do_cross || do_column)
        {
            // Make sure that the top boundary is taken into account in case of fluxes.
            const int kend = gd.kstart + array.dim(2);
            add_ghost_cells(
                    tmp->fld.data(), array.ptr(),
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, kend,
                    gd.igc, gd.jgc, gd.kgc,
                    gd.icells, gd.ijcells,
                    gd.imax, gd.imax*gd.jmax);
        }

        if (do_stats)
            stats.calc_stats(name, *tmp, no_offset, no_threshold);

        if (do_cross)
        {
            if (std::find(crosslist.begin(), crosslist.end(), name) != crosslist.end())
                cross.cross_simple(tmp->fld.data(), name, iotime, loc);
        }

        if (do_column)
            column.calc_column(name, tmp->fld.data(), no_offset);
    };

    try
    {
        if (sw_longwave)
        {
            exec_longwave(
                    thermo, timeloop, stats,
                    flux_up, flux_dn, flux_net,
                    t_lay_a, t_lev_a, t_sfc_a, h2o_a, clwp_a, ciwp_a,
                    compute_clouds);

            save_stats_and_cross(flux_up, "lw_flux_up", gd.wloc);
            save_stats_and_cross(flux_dn, "lw_flux_dn", gd.wloc);

            if (sw_clear_sky_stats)
            {
                exec_longwave(
                        thermo, timeloop, stats,
                        flux_up, flux_dn, flux_net,
                        t_lay_a, t_lev_a, t_sfc_a, h2o_a, clwp_a, ciwp_a,
                        !compute_clouds);

                save_stats_and_cross(flux_up, "lw_flux_up_clear", gd.wloc);
                save_stats_and_cross(flux_dn, "lw_flux_dn_clear", gd.wloc);
            }
        }

        if (sw_shortwave)
        {
            Array<TF,2> flux_dn_dir({gd.imax*gd.jmax, gd.ktot+1});

            if (!is_day(mu0))
            {
                flux_up.fill(0.);
                flux_dn.fill(0.);
                flux_dn_dir.fill(0.);
                flux_net.fill(0.);
            }

            if (is_day(mu0))
            {
                exec_shortwave(
                        thermo, timeloop, stats,
                        flux_up, flux_dn, flux_dn_dir, flux_net,
                        t_lay_a, t_lev_a, h2o_a, clwp_a, ciwp_a,
                        compute_clouds);
            }

            save_stats_and_cross(flux_up,     "sw_flux_up"    , gd.wloc);
            save_stats_and_cross(flux_dn,     "sw_flux_dn"    , gd.wloc);
            save_stats_and_cross(flux_dn_dir, "sw_flux_dn_dir", gd.wloc);

            if (sw_clear_sky_stats)
            {
                if (is_day(mu0))
                {
                    exec_shortwave(
                            thermo, timeloop, stats,
                            flux_up, flux_dn, flux_dn_dir, flux_net,
                            t_lay_a, t_lev_a, h2o_a, clwp_a, ciwp_a,
                            !compute_clouds);
                }

                save_stats_and_cross(flux_up,     "sw_flux_up_clear"    , gd.wloc);
                save_stats_and_cross(flux_dn,     "sw_flux_dn_clear"    , gd.wloc);
                save_stats_and_cross(flux_dn_dir, "sw_flux_dn_dir_clear", gd.wloc);
            }

            stats.set_time_series("sza", std::acos(mu0));
            stats.set_time_series("sw_flux_dn_toa", sw_flux_dn_col({1,n_lev_col}));
        }
    }
    catch (std::exception& e)
    {
        #ifdef USEMPI
        std::cout << "SINGLE PROCESS EXCEPTION: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        #else
        throw;
        #endif
    }

    fields.release_tmp(tmp);
}

template<typename TF>
void Radiation_rrtmgp<TF>::exec_longwave(
        Thermo<TF>& thermo, Timeloop<TF>& timeloop, Stats<TF>& stats,
        Array<TF,2>& flux_up, Array<TF,2>& flux_dn, Array<TF,2>& flux_net,
        const Array<TF,2>& t_lay, const Array<TF,2>& t_lev, const Array<TF,1>& t_sfc,
        const Array<TF,2>& h2o, const Array<TF,2>& clwp, const Array<TF,2>& ciwp,
        const bool compute_clouds)
{
    // How many profiles are solved simultaneously?
    constexpr int n_col_block = 4;

    auto& gd = grid.get_grid_data();

    const int n_lay = gd.ktot;
    const int n_lev = gd.ktot+1;
    const int n_col = gd.imax*gd.jmax;

    const int n_blocks = n_col / n_col_block;
    const int n_col_block_left = n_col % n_col_block;

    // Store the number of bands and gpt in a variable.
    const int n_bnd = kdist_lw->get_nband();
    const int n_gpt = kdist_lw->get_ngpt();

    // Set the number of angles to 1.
    const int n_ang = 1;

    // Check the dimension ordering. The top is not at 1 in MicroHH, but the surface is.
    const int top_at_1 = 0;

    // Define the pointers for the subsetting.
    std::unique_ptr<Optical_props_arry<TF>> optical_props_subset =
            std::make_unique<Optical_props_1scl<TF>>(n_col_block, n_lay, *kdist_lw);
    std::unique_ptr<Source_func_lw<TF>> sources_subset =
            std::make_unique<Source_func_lw<TF>>(n_col_block, n_lay, *kdist_lw);
    std::unique_ptr<Optical_props_1scl<TF>> cloud_optical_props_subset =
            std::make_unique<Optical_props_1scl<TF>>(n_col_block, n_lay, *cloud_lw);

    std::unique_ptr<Optical_props_arry<TF>> optical_props_left =
            std::make_unique<Optical_props_1scl<TF>>(n_col_block_left, n_lay, *kdist_lw);
    std::unique_ptr<Source_func_lw<TF>> sources_left =
            std::make_unique<Source_func_lw<TF>>(n_col_block_left, n_lay, *kdist_lw);
    std::unique_ptr<Optical_props_1scl<TF>> cloud_optical_props_left =
            std::make_unique<Optical_props_1scl<TF>>(n_col_block_left, n_lay, *cloud_lw);

    // Define the arrays that contain the subsets.
    const std::vector<TF>& p  = thermo.get_basestate_vector("p");
    const std::vector<TF>& ph = thermo.get_basestate_vector("ph");
    Array<TF,2> p_lay(std::vector<TF>(p.begin()  + gd.kstart, p.begin()  + gd.kend  ), {1, n_lay});
    Array<TF,2> p_lev(std::vector<TF>(ph.begin() + gd.kstart, ph.begin() + gd.kend+1), {1, n_lev});

    Array<TF,2> emis_sfc(std::vector<TF>(n_bnd, this->emis_sfc), {n_bnd, 1});

    gas_concs.set_vmr("h2o", h2o);
    Array<TF,2> col_dry({n_col, n_lay});
    Gas_optics_rrtmgp<TF>::get_col_dry(col_dry, gas_concs.get_vmr("h2o"), p_lev.subset({{ {1, n_col}, {1, n_lev} }}));

    // Lambda function for solving optical properties subset.
    auto call_kernels = [&](
            const int col_s_in, const int col_e_in,
            std::unique_ptr<Optical_props_arry<TF>>& optical_props_subset_in,
            std::unique_ptr<Optical_props_1scl<TF>>& cloud_optical_props_in,
            Source_func_lw<TF>& sources_subset_in,
            const Array<TF,2>& emis_sfc_subset_in,
            const Array<TF,2>& lw_flux_dn_inc_subset_in,
            std::unique_ptr<Fluxes_broadband<TF>>& fluxes)
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
                t_lev  .subset({{ {col_s_in, col_e_in}, {1, n_lev} }}) );

        // 2. Solve the cloud optical properties.
        if (compute_clouds)
        {
            Array<TF,2> clwp_subset(clwp.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}));
            Array<TF,2> ciwp_subset(ciwp.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}));

            // Compute the effective droplet radius.
            Array<TF,2> rel({n_col_in, n_lay});
            Array<TF,2> rei({n_col_in, n_lay});

            const TF sig_g = 1.34;
            const TF fac = std::exp(std::log(sig_g)*std::log(sig_g)); // no conversion to micron yet.

            // CvH: Numbers according to RCEMIP.
            const TF Nc0 = 100.e6;
            const TF Ni0 = 1.e5;

            const TF four_third_pi_Nc0_rho_w = (4./3.)*M_PI*Nc0*Constants::rho_w<TF>;
            const TF four_third_pi_Ni0_rho_i = (4./3.)*M_PI*Ni0*Constants::rho_i<TF>;

            for (int ilay=1; ilay<=n_lay; ++ilay)
            {
                // const TF layer_mass = (p_lev({1, ilay}) - p_lev({1, ilay+1})) / Constants::grav<TF>;
                const TF layer_thickness = gd.dz[ilay + gd.kstart - 1];

                for (int icol=1; icol<=n_col_in; ++icol)
                {
                    // Parametrization according to Martin et al., 1994 JAS. Fac multiplication taken from DALES.
                    // CvH: Potentially better using moments from microphysics.
                    TF rel_value = clwp_subset({icol, ilay}) > TF(0.) ? 
                        1.e6 * fac * std::pow((clwp_subset({icol, ilay})/layer_thickness) / four_third_pi_Nc0_rho_w, (1./3.)) : TF(0.);

                    // Limit the values between 2.5 and 21.5 (limits of cloud optics lookup table).
                    rel({icol, ilay}) = std::max(TF(2.5), std::min(rel_value, TF(21.5)));

                    // Calculate the effective radius of ice from the mass and the number concentration.
                    TF rei_value = ciwp_subset({icol, ilay}) > TF(0.) ? 
                        1.e6 * std::pow((ciwp_subset({icol, ilay})/layer_thickness) / four_third_pi_Ni0_rho_i, (1./3.)) : TF(0.);

                    // Limit the values between 10. and 180 (limits of cloud optics lookup table).
                    rei({icol, ilay}) = std::max(TF(10.), std::min(rei_value, TF(180.)));
                }
            }

            // Convert to g/m2.
            for (int i=0; i<clwp_subset.size(); ++i)
                clwp_subset.v()[i] *= 1e3;

            for (int i=0; i<ciwp_subset.size(); ++i)
                ciwp_subset.v()[i] *= 1e3;

            cloud_lw->cloud_optics(
                    clwp_subset, ciwp_subset,
                    rel, rei,
                    *cloud_optical_props_in);

            // Add the cloud optical props to the gas optical properties.
            add_to(
                    dynamic_cast<Optical_props_1scl<TF>&>(*optical_props_subset_in),
                    dynamic_cast<Optical_props_1scl<TF>&>(*cloud_optical_props_in));
        }

        Array<TF,3> gpt_flux_up({n_col_in, n_lev, n_gpt});
        Array<TF,3> gpt_flux_dn({n_col_in, n_lev, n_gpt});

        Rte_lw<TF>::rte_lw(
                optical_props_subset_in,
                top_at_1,
                sources_subset_in,
                emis_sfc_subset_in,
                lw_flux_dn_inc_subset_in,
                gpt_flux_up, gpt_flux_dn,
                n_ang);

        fluxes->reduce(gpt_flux_up, gpt_flux_dn, optical_props_subset_in, top_at_1);

        // Copy the data to the output.
        for (int ilev=1; ilev<=n_lev; ++ilev)
            for (int icol=1; icol<=n_col_in; ++icol)
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

        Array<TF,2> emis_sfc_subset = emis_sfc.subset({{ {1, n_bnd}, {col_s, col_e} }});
        Array<TF,2> lw_flux_dn_inc_subset = lw_flux_dn_inc.subset({{ {col_s, col_e}, {1, n_gpt} }});
        std::unique_ptr<Fluxes_broadband<TF>> fluxes_subset =
                std::make_unique<Fluxes_broadband<TF>>(n_col_block, n_lev);

        call_kernels(
                col_s, col_e,
                optical_props_subset,
                cloud_optical_props_subset,
                *sources_subset,
                emis_sfc_subset,
                lw_flux_dn_inc_subset,
                fluxes_subset);
    }

    if (n_col_block_left > 0)
    {
        const int col_s = n_col - n_col_block_left + 1;
        const int col_e = n_col;

        Array<TF,2> emis_sfc_left = emis_sfc.subset({{ {1, n_bnd}, {col_s, col_e} }});
        Array<TF,2> lw_flux_dn_inc_left = lw_flux_dn_inc.subset({{ {col_s, col_e}, {1, n_gpt} }});
        std::unique_ptr<Fluxes_broadband<TF>> fluxes_left =
                std::make_unique<Fluxes_broadband<TF>>(n_col_block_left, n_lev);

        call_kernels(
                col_s, col_e,
                optical_props_left,
                cloud_optical_props_left,
                *sources_left,
                emis_sfc_left,
                lw_flux_dn_inc_left,
                fluxes_left);
    }
}

template<typename TF>
void Radiation_rrtmgp<TF>::exec_shortwave(
        Thermo<TF>& thermo, Timeloop<TF>& timeloop, Stats<TF>& stats,
        Array<TF,2>& flux_up, Array<TF,2>& flux_dn, Array<TF,2>& flux_dn_dir, Array<TF,2>& flux_net,
        const Array<TF,2>& t_lay, const Array<TF,2>& t_lev,
        const Array<TF,2>& h2o, const Array<TF,2>& clwp, const Array<TF,2>& ciwp,
        const bool compute_clouds)
{
    // How many profiles are solved simultaneously?
    constexpr int n_col_block = 4;

    auto& gd = grid.get_grid_data();

    const int n_lay = gd.ktot;
    const int n_lev = gd.ktot+1;
    const int n_col = gd.imax*gd.jmax;

    const int n_blocks = n_col / n_col_block;
    const int n_col_block_left = n_col % n_col_block;

    // Store the number of bands and gpt in a variable.
    const int n_bnd = kdist_sw->get_nband();
    const int n_gpt = kdist_sw->get_ngpt();

    // Check the dimension ordering. The top is not at 1 in MicroHH, but the surface is.
    const int top_at_1 = 0;

    // Define the pointers for the subsetting.
    std::unique_ptr<Optical_props_arry<TF>> optical_props_subset =
            std::make_unique<Optical_props_2str<TF>>(n_col_block, n_lay, *kdist_sw);
    std::unique_ptr<Optical_props_2str<TF>> cloud_optical_props_subset =
            std::make_unique<Optical_props_2str<TF>>(n_col_block, n_lay, *cloud_sw);

    std::unique_ptr<Optical_props_arry<TF>> optical_props_left =
            std::make_unique<Optical_props_2str<TF>>(n_col_block_left, n_lay, *kdist_sw);
    std::unique_ptr<Optical_props_2str<TF>> cloud_optical_props_left =
            std::make_unique<Optical_props_2str<TF>>(n_col_block_left, n_lay, *cloud_sw);

    // Define the arrays that contain the subsets.
    std::vector<TF> p  = thermo.get_basestate_vector("p");
    std::vector<TF> ph = thermo.get_basestate_vector("ph");
    Array<TF,2> p_lay(std::vector<TF>(p.begin()  + gd.kstart, p.begin()  + gd.kend  ), {1, n_lay});
    Array<TF,2> p_lev(std::vector<TF>(ph.begin() + gd.kstart, ph.begin() + gd.kend+1), {1, n_lev});

    // Create the boundary conditions
    Array<TF,1> mu0(std::vector<TF>(1, this->mu0), {1});
    Array<TF,2> sfc_alb_dir(std::vector<TF>(n_bnd, this->sfc_alb_dir), {n_bnd, 1});
    Array<TF,2> sfc_alb_dif(std::vector<TF>(n_bnd, this->sfc_alb_dif), {n_bnd, 1});

    // Create the field for the top of atmosphere source.
    Array<TF,2> toa_src({n_col, n_gpt});

    gas_concs.set_vmr("h2o", h2o);
    Array<TF,2> col_dry({n_col, n_lay});
    Gas_optics_rrtmgp<TF>::get_col_dry(col_dry, gas_concs.get_vmr("h2o"), p_lev.subset({{ {1, n_col}, {1, n_lev} }}));

    // Lambda function for solving optical properties subset.
    auto call_kernels = [&](
            const int col_s_in, const int col_e_in,
            std::unique_ptr<Optical_props_arry<TF>>& optical_props_subset_in,
            std::unique_ptr<Optical_props_2str<TF>>& cloud_optical_props_in,
            const Array<TF,1>& mu0_subset_in,
            const Array<TF,2>& toa_src_subset_in,
            const Array<TF,2>& sfc_alb_dir_subset_in,
            const Array<TF,2>& sfc_alb_dif_subset_in,
            const Array<TF,2>& sw_flux_dn_dif_inc_subset_in,
            std::unique_ptr<Fluxes_broadband<TF>>& fluxes)
    {
        const int n_col_in = col_e_in - col_s_in + 1;

        Gas_concs<TF> gas_concs_subset(gas_concs, col_s_in, n_col_in);
        Array<TF,2> toa_src_dummy({n_col_in, n_gpt});

        // 1. Solve the gas optical properties.
        kdist_sw->gas_optics(
                p_lay.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}),
                p_lev.subset({{ {col_s_in, col_e_in}, {1, n_lev} }}),
                t_lay.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}),
                gas_concs_subset,
                optical_props_subset_in,
                toa_src_dummy,
                col_dry.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}) );

        // 2. Solve the cloud optical properties.
        if (compute_clouds)
        {
            Array<TF,2> clwp_subset(clwp.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}));
            Array<TF,2> ciwp_subset(ciwp.subset({{ {col_s_in, col_e_in}, {1, n_lay} }}));

            // Compute the effective droplet radius.
            Array<TF,2> rel({n_col_in, n_lay});
            Array<TF,2> rei({n_col_in, n_lay});

            const TF sig_g = 1.34;
            const TF fac = std::exp(std::log(sig_g)*std::log(sig_g)); // no conversion to micron yet.

            // CvH: Numbers according to RCEMIP.
            const TF Nc0 = 100.e6;
            const TF Ni0 = 1.e5;

            const TF four_third_pi_Nc0_rho_w = (4./3.)*M_PI*Nc0*Constants::rho_w<TF>;
            const TF four_third_pi_Ni0_rho_i = (4./3.)*M_PI*Ni0*Constants::rho_i<TF>;

            for (int ilay=1; ilay<=n_lay; ++ilay)
            {
                // const TF layer_mass = (p_lev({1, ilay}) - p_lev({1, ilay+1})) / Constants::grav<TF>;
                const TF layer_thickness = gd.dz[ilay + gd.kstart - 1];

                for (int icol=1; icol<=n_col_in; ++icol)
                {
                    // Parametrization according to Martin et al., 1994 JAS. Fac multiplication taken from DALES.
                    // CvH: Potentially better using moments from microphysics.
                    TF rel_value = clwp_subset({icol, ilay}) > TF(0.) ? 
                        1.e6 * fac * std::pow((clwp_subset({icol, ilay})/layer_thickness) / four_third_pi_Nc0_rho_w, (1./3.)) : TF(0.);

                    // Limit the values between 2.5 and 21.5 (limits of cloud optics lookup table).
                    rel({icol, ilay}) = std::max(TF(2.5), std::min(rel_value, TF(21.5)));

                    // Calculate the effective radius of ice from the mass and the number concentration.
                    TF rei_value = ciwp_subset({icol, ilay}) > TF(0.) ? 
                        1.e6 * std::pow((ciwp_subset({icol, ilay})/layer_thickness) / four_third_pi_Ni0_rho_i, (1./3.)) : TF(0.);

                    // Limit the values between 10. and 180 (limits of cloud optics lookup table).
                    rei({icol, ilay}) = std::max(TF(10.), std::min(rei_value, TF(180.)));
                }
            }

            // Convert to g/m2.
            for (int i=0; i<clwp_subset.size(); ++i)
                clwp_subset.v()[i] *= 1e3;

            for (int i=0; i<ciwp_subset.size(); ++i)
                ciwp_subset.v()[i] *= 1e3;

            cloud_sw->cloud_optics(
                    clwp_subset, ciwp_subset,
                    rel, rei,
                    *cloud_optical_props_in);

            cloud_optical_props_in->delta_scale();

            // Add the cloud optical props to the gas optical properties.
            add_to(
                    dynamic_cast<Optical_props_2str<TF>&>(*optical_props_subset_in),
                    dynamic_cast<Optical_props_2str<TF>&>(*cloud_optical_props_in));
        }

        // 3. Solve the fluxes.
        Array<TF,3> gpt_flux_up    ({n_col_in, n_lev, n_gpt});
        Array<TF,3> gpt_flux_dn    ({n_col_in, n_lev, n_gpt});
        Array<TF,3> gpt_flux_dn_dir({n_col_in, n_lev, n_gpt});

        Rte_sw<TF>::rte_sw(
                optical_props_subset_in,
                top_at_1,
                mu0_subset_in,
                toa_src_subset_in,
                sfc_alb_dir_subset_in,
                sfc_alb_dif_subset_in,
                sw_flux_dn_dif_inc_subset_in,
                gpt_flux_up,
                gpt_flux_dn,
                gpt_flux_dn_dir);

        // 4. Reduce the fluxes to the needed information.
        fluxes->reduce(
                gpt_flux_up, gpt_flux_dn, gpt_flux_dn_dir,
                optical_props_subset_in, top_at_1);

        // Copy the data to the output.
        for (int ilev=1; ilev<=n_lev; ++ilev)
            for (int icol=1; icol<=n_col_in; ++icol)
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

        Array<TF,1> mu0_subset = mu0.subset({{ {col_s, col_e} }});
        Array<TF,2> toa_src_subset = sw_flux_dn_dir_inc.subset({{ {col_s, col_e}, {1, n_gpt} }});
        Array<TF,2> sfc_alb_dir_subset = sfc_alb_dir.subset({{ {1, n_bnd}, {col_s, col_e} }});
        Array<TF,2> sfc_alb_dif_subset = sfc_alb_dif.subset({{ {1, n_bnd}, {col_s, col_e} }});
        Array<TF,2> sw_flux_dn_dif_inc_subset = sw_flux_dn_dif_inc.subset({{ {col_s, col_e}, {1, n_gpt} }});

        std::unique_ptr<Fluxes_broadband<TF>> fluxes_subset =
                std::make_unique<Fluxes_broadband<TF>>(n_col_block, n_lev);

        call_kernels(
                col_s, col_e,
                optical_props_subset,
                cloud_optical_props_subset,
                mu0_subset,
                toa_src_subset,
                sfc_alb_dir_subset,
                sfc_alb_dif_subset,
                sw_flux_dn_dif_inc_subset,
                fluxes_subset);
    }

    if (n_col_block_left > 0)
    {
        const int col_s = n_col - n_col_block_left + 1;
        const int col_e = n_col;

        Array<TF,1> mu0_left = mu0.subset({{ {col_s, col_e} }});
        Array<TF,2> toa_src_left = sw_flux_dn_dir_inc.subset({{ {col_s, col_e}, {1, n_gpt} }});
        Array<TF,2> sfc_alb_dir_left = sfc_alb_dir.subset({{ {1, n_bnd}, {col_s, col_e} }});
        Array<TF,2> sfc_alb_dif_left = sfc_alb_dif.subset({{ {1, n_bnd}, {col_s, col_e} }});
        Array<TF,2> sw_flux_dn_dif_inc_left = sw_flux_dn_dif_inc.subset({{ {col_s, col_e}, {1, n_gpt} }});

        std::unique_ptr<Fluxes_broadband<TF>> fluxes_left =
                std::make_unique<Fluxes_broadband<TF>>(n_col_block_left, n_lev);

        call_kernels(
                col_s, col_e,
                optical_props_left,
                cloud_optical_props_left,
                mu0_left,
                toa_src_left,
                sfc_alb_dir_left,
                sfc_alb_dif_left,
                sw_flux_dn_dif_inc_left,
                fluxes_left);
    }
}


template<typename TF>
bool Radiation_rrtmgp<TF>::is_day(const TF mu0)
{
    if (mu0 > Constants::mu0_min<TF>)
        return true;

    return false;
}

#ifdef FLOAT_SINGLE
template class Radiation_rrtmgp<float>;
#else
template class Radiation_rrtmgp<double>;
#endif
