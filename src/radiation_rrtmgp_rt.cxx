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

#include "radiation_rrtmgp_rt.h"
#include "radiation_rrtmgp_functions.h"
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

// RRTMGP RT headers.
#include "Optical_props_rt.h"
#include "Gas_optics_rrtmgp_rt.h"
#include "Fluxes_rt.h"
#include "Rte_lw_rt.h"
#include "Rte_sw_rt.h"
#include "Source_functions_rt.h"
#include "Cloud_optics_rt.h"


// IMPORTANT: The RTE+RRTMGP code sets the precision using a compiler flag RTE_RRTMGP_SINGLE_PRECISION, which defines
// a type Float that is float or double depending on the flag. The type of Float is coupled to the TF switch in MicroHH.
// To avoid confusion, we limit the use of TF in the code to the class headers and use Float.
using namespace Radiation_rrtmgp_functions;
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


    void load_gas_concs(
            Gas_concs& gas_concs, Netcdf_handle& input_nc, const std::string& dim_name)
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
                    gas_concs.set_vmr(gas_name, input_nc.get_variable<Float>(gas_name));
                }
                else if (n_dims == 1)
                {
                    if (dims.at(dim_name) == n_lay)
                        gas_concs.set_vmr(gas_name,
                                Array<Float,1>(input_nc.get_variable<Float>(gas_name, {n_lay}), {n_lay}));
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


    Gas_optics_rrtmgp load_and_init_gas_optics(
            Master& master,
            const Gas_concs& gas_concs,
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
        Array<Float,2> band_lims(coef_nc.get_variable<Float>("bnd_limits_wavenumber", {n_bnds, 2}), {2, n_bnds});
        Array<int,2> band2gpt(coef_nc.get_variable<int>("bnd_limits_gpt", {n_bnds, 2}), {2, n_bnds});
        Array<Float,1> press_ref(coef_nc.get_variable<Float>("press_ref", {n_press}), {n_press});
        Array<Float,1> temp_ref(coef_nc.get_variable<Float>("temp_ref", {n_temps}), {n_temps});

        Float temp_ref_p = coef_nc.get_variable<Float>("absorption_coefficient_ref_P");
        Float temp_ref_t = coef_nc.get_variable<Float>("absorption_coefficient_ref_T");
        Float press_ref_trop = coef_nc.get_variable<Float>("press_ref_trop");

        Array<Float,3> kminor_lower(
                coef_nc.get_variable<Float>("kminor_lower", {n_temps, n_mixingfracs, n_contributors_lower}),
                {n_contributors_lower, n_mixingfracs, n_temps});
        Array<Float,3> kminor_upper(
                coef_nc.get_variable<Float>("kminor_upper", {n_temps, n_mixingfracs, n_contributors_upper}),
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

        Array<Bool,1> minor_scales_with_density_lower(
                coef_nc.get_variable<Bool>("minor_scales_with_density_lower", {n_minor_absorber_intervals_lower}),
                {n_minor_absorber_intervals_lower});
        Array<Bool,1> minor_scales_with_density_upper(
                coef_nc.get_variable<Bool>("minor_scales_with_density_upper", {n_minor_absorber_intervals_upper}),
                {n_minor_absorber_intervals_upper});

        Array<Bool,1> scale_by_complement_lower(
                coef_nc.get_variable<Bool>("scale_by_complement_lower", {n_minor_absorber_intervals_lower}),
                {n_minor_absorber_intervals_lower});
        Array<Bool,1> scale_by_complement_upper(
                coef_nc.get_variable<Bool>("scale_by_complement_upper", {n_minor_absorber_intervals_upper}),
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

        Array<Float,3> vmr_ref(
                coef_nc.get_variable<Float>("vmr_ref", {n_temps, n_extabsorbers, n_layers}),
                {n_layers, n_extabsorbers, n_temps});

        Array<Float,4> kmajor(
                coef_nc.get_variable<Float>("kmajor", {n_temps, n_press+1, n_mixingfracs, n_gpts}),
                {n_gpts, n_mixingfracs, n_press+1, n_temps});

        // Keep the size at zero, if it does not exist.
        Array<Float,3> rayl_lower;
        Array<Float,3> rayl_upper;

        if (coef_nc.variable_exists("rayl_lower"))
        {
            rayl_lower.set_dims({n_gpts, n_mixingfracs, n_temps});
            rayl_upper.set_dims({n_gpts, n_mixingfracs, n_temps});
            rayl_lower = coef_nc.get_variable<Float>("rayl_lower", {n_temps, n_mixingfracs, n_gpts});
            rayl_upper = coef_nc.get_variable<Float>("rayl_upper", {n_temps, n_mixingfracs, n_gpts});
        }

        // Is it really LW if so read these variables as well.
        if (coef_nc.variable_exists("totplnk"))
        {
            int n_internal_sourcetemps = coef_nc.get_dimension_size("temperature_Planck");

            Array<Float,2> totplnk(
                    coef_nc.get_variable<Float>( "totplnk", {n_bnds, n_internal_sourcetemps}),
                    {n_internal_sourcetemps, n_bnds});
            Array<Float,4> planck_frac(
                    coef_nc.get_variable<Float>("plank_fraction", {n_temps, n_press+1, n_mixingfracs, n_gpts}),
                    {n_gpts, n_mixingfracs, n_press+1, n_temps});

            // Construct the k-distribution.
            return Gas_optics_rrtmgp(
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
            Array<Float,1> solar_src_quiet(
                    coef_nc.get_variable<Float>("solar_source_quiet", {n_gpts}), {n_gpts});
            Array<Float,1> solar_src_facular(
                    coef_nc.get_variable<Float>("solar_source_facular", {n_gpts}), {n_gpts});
            Array<Float,1> solar_src_sunspot(
                    coef_nc.get_variable<Float>("solar_source_sunspot", {n_gpts}), {n_gpts});

            Float tsi = coef_nc.get_variable<Float>("tsi_default");
            Float mg_index = coef_nc.get_variable<Float>("mg_default");
            Float sb_index = coef_nc.get_variable<Float>("sb_default");

            return Gas_optics_rrtmgp(
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


    Cloud_optics load_and_init_cloud_optics(
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

        Array<Float,2> band_lims_wvn(coef_nc.get_variable<Float>("bnd_limits_wavenumber", {n_band, 2}), {2, n_band});

        // Read look-up table constants.
        Float radliq_lwr = coef_nc.get_variable<Float>("radliq_lwr");
        Float radliq_upr = coef_nc.get_variable<Float>("radliq_upr");
        Float radliq_fac = coef_nc.get_variable<Float>("radliq_fac");

        Float radice_lwr = coef_nc.get_variable<Float>("radice_lwr");
        Float radice_upr = coef_nc.get_variable<Float>("radice_upr");
        Float radice_fac = coef_nc.get_variable<Float>("radice_fac");

        Array<Float,2> lut_extliq(
                coef_nc.get_variable<Float>("lut_extliq", {n_band, n_size_liq}), {n_size_liq, n_band});
        Array<Float,2> lut_ssaliq(
                coef_nc.get_variable<Float>("lut_ssaliq", {n_band, n_size_liq}), {n_size_liq, n_band});
        Array<Float,2> lut_asyliq(
                coef_nc.get_variable<Float>("lut_asyliq", {n_band, n_size_liq}), {n_size_liq, n_band});

        Array<Float,3> lut_extice(
                coef_nc.get_variable<Float>("lut_extice", {n_rghice, n_band, n_size_ice}), {n_size_ice, n_band, n_rghice});
        Array<Float,3> lut_ssaice(
                coef_nc.get_variable<Float>("lut_ssaice", {n_rghice, n_band, n_size_ice}), {n_size_ice, n_band, n_rghice});
        Array<Float,3> lut_asyice(
                coef_nc.get_variable<Float>("lut_asyice", {n_rghice, n_band, n_size_ice}), {n_size_ice, n_band, n_rghice});

        return Cloud_optics(
                band_lims_wvn,
                radliq_lwr, radliq_upr, radliq_fac,
                radice_lwr, radice_upr, radice_fac,
                lut_extliq, lut_ssaliq, lut_asyliq,
                lut_extice, lut_ssaice, lut_asyice);
    }

}

template<typename TF>
Radiation_rrtmgp_rt<TF>::Radiation_rrtmgp_rt(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
        Radiation<TF>(masterin, gridin, fieldsin, inputin)
{
    swradiation = "rrtmgp_rt";

    sw_longwave  = inputin.get_item<bool>("radiation", "swlongwave" , "", true);
    sw_shortwave = inputin.get_item<bool>("radiation", "swshortwave", "", true);
    sw_fixed_sza = inputin.get_item<bool>("radiation", "swfixedsza", "", true);

    sw_clear_sky_stats = inputin.get_item<bool>("radiation", "swclearskystats", "", false);

    dt_rad = inputin.get_item<double>("radiation", "dt_rad", "");

    t_sfc       = inputin.get_item<Float>("radiation", "t_sfc"      , "");
    emis_sfc    = inputin.get_item<Float>("radiation", "emis_sfc"   , "");
    sfc_alb_dir = inputin.get_item<Float>("radiation", "sfc_alb_dir", "");
    sfc_alb_dif = inputin.get_item<Float>("radiation", "sfc_alb_dif", "");
    tsi_scaling = inputin.get_item<Float>("radiation", "tsi_scaling", "", -999.);

    if (sw_fixed_sza)
    {
        const Float sza = inputin.get_item<Float>("radiation", "sza", "");
        mu0 = std::cos(sza);
    }
    else
    {
        lat = inputin.get_item<Float>("radiation", "lat", "");
        lon = inputin.get_item<Float>("radiation", "lon", "");
    }

    auto& gd = grid.get_grid_data();
    fields.init_diagnostic_field("thlt_rad", "Tendency by radiation", "K s-1", "radiation", gd.sloc);

    fields.init_diagnostic_field("lw_flux_up", "Longwave upwelling flux", "W m-2", "radiation", gd.wloc);
    fields.init_diagnostic_field("lw_flux_dn", "Longwave downwelling flux", "W m-2", "radiation", gd.wloc);
    
    if (sw_clear_sky_stats)
    {
        fields.init_diagnostic_field("lw_flux_up_clear", "Clear-sky longwave upwelling flux", "W m-2", "radiation", gd.wloc);
        fields.init_diagnostic_field("lw_flux_dn_clear", "Clear-sky longwave downwelling flux", "W m-2", "radiation", gd.wloc);
    }
    
    fields.init_diagnostic_field("sw_flux_up", "Shortwave upwelling flux", "W m-2", "radiation", gd.wloc);
    fields.init_diagnostic_field("sw_flux_dn", "Shortwave downwelling flux", "W m-2", "radiation", gd.wloc);
    fields.init_diagnostic_field("sw_flux_dn_dir", "Shortwave direct downwelling flux", "W m-2", "radiation", gd.wloc);
    
    fields.init_diagnostic_field("sw_heat_dir_rt", "Heating rates from direct raytraced radiation", "K s-1", "radiation", gd.sloc);
    fields.init_diagnostic_field("sw_heat_dif_rt", "Heating rates from diffuse raytraced radiation", "k s-1", "radiation", gd.sloc);

    if (sw_clear_sky_stats)
    {
        fields.init_diagnostic_field("sw_flux_up_clear", "Clear-sky shortwave upwelling flux", "W m-2", "radiation", gd.wloc);
        fields.init_diagnostic_field("sw_flux_dn_clear", "Clear-sky shortwave downwelling flux", "W m-2", "radiation", gd.wloc);
        fields.init_diagnostic_field("sw_flux_dn_dir_clear", "Clear-sky shortwave direct downwelling flux", "W m-2", "radiation", gd.wloc);
    }
}


template<typename TF>
void Radiation_rrtmgp_rt<TF>::init(Timeloop<TF>& timeloop)
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

    sw_flux_sfc_dir_rt.resize(gd.ijcells);
    sw_flux_sfc_dif_rt.resize(gd.ijcells);
    sw_flux_sfc_up_rt.resize(gd.ijcells);
    sw_flux_tod_up_rt.resize(gd.ijcells);

}


template<typename TF>
unsigned long Radiation_rrtmgp_rt<TF>::get_time_limit(unsigned long itime)
{
    unsigned long idtlim = idt_rad - itime % idt_rad;
    return idtlim;
}


template<typename TF>
void Radiation_rrtmgp_rt<TF>::create(
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
        
        stats.add_time_series("sw_flux_sfc_dir_rt", "raytraced shortwave downwelling direct flux at the surface", "W m-2", group_name);
        stats.add_time_series("sw_flux_sfc_dif_rt", "raytraced shortwave downwelling diffuse flux at the surface", "W m-2", group_name);
        stats.add_time_series("sw_flux_sfc_up_rt",  "raytraced shortwave downwelling upwelling flux at the surface", "W m-2", group_name);
        stats.add_time_series("sw_flux_tod_up_rt",  "raytraced shortwave downwelling upwelling flux at toa", "W m-2", group_name);
    }

    // Get the allowed cross sections from the cross list
    std::vector<std::string> allowed_crossvars_radiation;

    if (sw_shortwave)
    {
        allowed_crossvars_radiation.push_back("sw_flux_up");
        allowed_crossvars_radiation.push_back("sw_flux_dn");
        allowed_crossvars_radiation.push_back("sw_flux_dn_dir");

        allowed_crossvars_radiation.push_back("sw_heat_dir_rt");
        allowed_crossvars_radiation.push_back("sw_heat_dif_rt");
        
        allowed_crossvars_radiation.push_back("sw_flux_sfc_dir_rt");
        allowed_crossvars_radiation.push_back("sw_flux_sfc_dif_rt");
        allowed_crossvars_radiation.push_back("sw_flux_sfc_up_rt");
        allowed_crossvars_radiation.push_back("sw_flux_tod_up_rt");
        
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
void Radiation_rrtmgp_rt<TF>::solve_longwave_column(
        std::unique_ptr<Optical_props_arry>& optical_props,
        Array<Float,2>& flux_up, Array<Float,2>& flux_dn, Array<Float,2>& flux_net,
        Array<Float,2>& flux_dn_inc, const Float p_top,
        const Gas_concs& gas_concs,
        const std::unique_ptr<Gas_optics_rrtmgp>& kdist_lw,
        const std::unique_ptr<Source_func_lw>& sources,
        const Array<Float,2>& col_dry,
        const Array<Float,2>& p_lay, const Array<Float,2>& p_lev,
        const Array<Float,2>& t_lay, const Array<Float,2>& t_lev,
        const Array<Float,1>& t_sfc, const Array<Float,2>& emis_sfc,
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

    std::unique_ptr<Fluxes_broadband> fluxes =
            std::make_unique<Fluxes_broadband>(n_col, n_lev);

    const int n_gpt = kdist_lw->get_ngpt();
    Array<Float,3> gpt_flux_up({n_col, n_lev, n_gpt});
    Array<Float,3> gpt_flux_dn({n_col, n_lev, n_gpt});

    Rte_lw::rte_lw(
            optical_props,
            top_at_1,
            *sources,
            emis_sfc,
            Array<Float,2>({n_col, n_gpt}),
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
    const Float fac_bot = (p_top - p_lev({1, idx_top})) / (p_lev({1, idx_bot}) - p_lev({1, idx_top}));
    const Float fac_top = 1. - fac_bot;

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
void Radiation_rrtmgp_rt<TF>::solve_shortwave_column(
        std::unique_ptr<Optical_props_arry>& optical_props,
        Array<Float,2>& flux_up, Array<Float,2>& flux_dn,
        Array<Float,2>& flux_dn_dir, Array<Float,2>& flux_net,
        Array<Float,2>& flux_dn_dir_inc, Array<Float,2>& flux_dn_dif_inc, const Float p_top,
        const Gas_concs& gas_concs,
        const Gas_optics_rrtmgp& kdist_sw,
        const Array<Float,2>& col_dry,
        const Array<Float,2>& p_lay, const Array<Float,2>& p_lev,
        const Array<Float,2>& t_lay, const Array<Float,2>& t_lev,
        const Array<Float,1>& mu0,
        const Array<Float,2>& sfc_alb_dir, const Array<Float,2>& sfc_alb_dif,
        const Float tsi_scaling,
        const int n_lay)
{
    const int n_col = 1;
    const int n_lev = n_lay + 1;

    // Check the dimension ordering.
    const int top_at_1 = p_lay({1, 1}) < p_lay({1, n_lay});

    // Create the field for the top of atmosphere source.
    const int n_gpt = kdist_sw.get_ngpt();
    Array<Float,2> toa_src({n_col, n_gpt});

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

    std::unique_ptr<Fluxes_broadband> fluxes =
            std::make_unique<Fluxes_broadband>(n_col, n_lev);

    Array<Float,3> gpt_flux_up    ({n_col, n_lev, n_gpt});
    Array<Float,3> gpt_flux_dn    ({n_col, n_lev, n_gpt});
    Array<Float,3> gpt_flux_dn_dir({n_col, n_lev, n_gpt});

    Rte_sw::rte_sw(
            optical_props,
            top_at_1,
            mu0,
            toa_src,
            sfc_alb_dir,
            sfc_alb_dif,
            Array<Float,2>(),
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
    const Float fac_bot = (p_top - p_lev({1, idx_top})) / (p_lev({1, idx_bot}) - p_lev({1, idx_top}));
    const Float fac_top = 1. - fac_bot;

    // Interpolate the top boundary conditions.
    for (int igpt=1; igpt<=n_gpt; ++igpt)
    {
        const Float flux_dn_tot = fac_bot * gpt_flux_dn    ({1, idx_bot, igpt}) + fac_top * gpt_flux_dn    ({1, idx_top, igpt});
        const Float flux_dn_dir = fac_bot * gpt_flux_dn_dir({1, idx_bot, igpt}) + fac_top * gpt_flux_dn_dir({1, idx_top, igpt});
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
void Radiation_rrtmgp_rt<TF>::create_column(
        Input& input, Netcdf_handle& input_nc, Thermo<TF>& thermo, Stats<TF>& stats)
{
    // 1. Load the available gas concentrations from the group of the netcdf file.
    Netcdf_handle& rad_nc = input_nc.get_group("radiation");

    load_gas_concs(gas_concs_col, rad_nc, "lay");

    // 2. Set the coordinate for the reference profiles in the stats, before calling the other creates.
    if (stats.get_switch() && (sw_longwave || sw_shortwave))
    {
        const int n_col = 1;
        const int n_lev = rad_nc.get_dimension_size("lev");
        Array<Float,2> p_lev(rad_nc.get_variable<Float>("p_lev", {n_lev, n_col}), {n_col, n_lev});

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
void Radiation_rrtmgp_rt<TF>::read_background_profiles(
        Netcdf_handle& rad_nc,
        const Gas_concs& gas_concs_col)
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

    p_lay_col = rad_nc.get_variable<Float>("p_lay", {n_lay_col, n_col});
    t_lay_col = rad_nc.get_variable<Float>("t_lay", {n_lay_col, n_col});
    p_lev_col = rad_nc.get_variable<Float>("p_lev", {n_lev_col, n_col});
    t_lev_col = rad_nc.get_variable<Float>("t_lev", {n_lev_col, n_col});

    if (rad_nc.variable_exists("col_dry"))
        col_dry = rad_nc.get_variable<Float>("col_dry", {n_lay_col, n_col});
    else
        Gas_optics_rrtmgp::get_col_dry(col_dry, gas_concs_col.get_vmr("h2o"), p_lev_col);
}


template<typename TF>
void Radiation_rrtmgp_rt<TF>::create_column_longwave(
        Input& input, Netcdf_handle& rad_nc, Thermo<TF>& thermo, Stats<TF>& stats,
        const Gas_concs& gas_concs_col)
{
    auto& gd = grid.get_grid_data();

    // Read the boundary conditions.
    // Set the surface temperature and emissivity.
    Array<Float,1> t_sfc({1});
    t_sfc({1}) = t_lev_col({1,1});

    const int n_bnd = kdist_lw->get_nband();
    Array<Float,2> emis_sfc({n_bnd, 1});
    for (int ibnd=1; ibnd<=n_bnd; ++ibnd)
        emis_sfc({ibnd, 1}) = this->emis_sfc;

    // Compute the longwave for the reference profile.
    sources_lw = std::make_unique<Source_func_lw>(n_col, n_lay_col, *kdist_lw);
    optical_props_lw = std::make_unique<Optical_props_1scl>(n_col, n_lay_col, *kdist_lw);

    lw_flux_up_col .set_dims({n_col, n_lev_col});
    lw_flux_dn_col .set_dims({n_col, n_lev_col});
    lw_flux_net_col.set_dims({n_col, n_lev_col});

    const int n_gpt = kdist_lw->get_ngpt();
    lw_flux_dn_inc.set_dims({n_col, n_gpt});

    solve_longwave_column(
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
void Radiation_rrtmgp_rt<TF>::create_column_shortwave(
        Input& input, Netcdf_handle& rad_nc, Thermo<TF>& thermo, Stats<TF>& stats,
        const Gas_concs& gas_concs_col)
{
    auto& gd = grid.get_grid_data();

    // Read the boundary conditions.
    const int n_bnd = kdist_sw->get_nband();

    optical_props_sw = std::make_unique<Optical_props_2str>(n_col, n_lay_col, *kdist_sw);

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
        Array<Float,2> sfc_alb_dir({n_bnd, n_col});
        Array<Float,2> sfc_alb_dif({n_bnd, n_col});

        for (int ibnd=1; ibnd<=n_bnd; ++ibnd)
        {
            sfc_alb_dir({ibnd, 1}) = this->sfc_alb_dir;
            sfc_alb_dif({ibnd, 1}) = this->sfc_alb_dif;
        }

        Array<Float,1> mu0({n_col});
        mu0({1}) = std::max(this->mu0, Constants::mu0_min<Float>);

        solve_shortwave_column(
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
void Radiation_rrtmgp_rt<TF>::create_solver(
        Input& input, Netcdf_handle& input_nc, Thermo<TF>& thermo,
        Stats<TF>& stats, Column<TF>& column)
{
    // 1. Load the available gas concentrations from the group of the netcdf file.
    Netcdf_handle& rad_input_nc = input_nc.get_group("init");
    load_gas_concs(gas_concs, rad_input_nc, "z");

    // 2. Pass the gas concentrations to the solver initializers.
    if (sw_longwave)
        create_solver_longwave(input, input_nc, thermo, stats, column, gas_concs);
    if (sw_shortwave)
        create_solver_shortwave(input, input_nc, thermo, stats, column, gas_concs);
}


template<typename TF>
void Radiation_rrtmgp_rt<TF>::create_solver_longwave(
        Input& input, Netcdf_handle& input_nc, Thermo<TF>& thermo,
        Stats<TF>& stats, Column<TF>& column,
        const Gas_concs& gas_concs)
{
    const std::string group_name = "radiation";

    // Set up the gas optics classes for long and shortwave.
    kdist_lw = std::make_unique<Gas_optics_rrtmgp>(
            load_and_init_gas_optics(master, gas_concs, "coefficients_lw.nc"));

    cloud_lw = std::make_unique<Cloud_optics>(
            load_and_init_cloud_optics(master, "cloud_coefficients_lw.nc"));

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
void Radiation_rrtmgp_rt<TF>::create_solver_shortwave(
        Input& input, Netcdf_handle& input_nc, Thermo<TF>& thermo,
        Stats<TF>& stats, Column<TF>& column,
        const Gas_concs& gas_concs)
{
    const std::string group_name = "radiation";

    // Set up the gas optics classes for long and shortwave.
    kdist_sw = std::make_unique<Gas_optics_rrtmgp>(
            load_and_init_gas_optics(master, gas_concs, "coefficients_sw.nc"));

    cloud_sw = std::make_unique<Cloud_optics>(
            load_and_init_cloud_optics(master, "cloud_coefficients_sw.nc"));

    // Set up the statistics.
    if (stats.get_switch())
    {
        stats.add_prof("sw_flux_up"    , "Shortwave upwelling flux"         , "W m-2", "zh", group_name);
        stats.add_prof("sw_flux_dn"    , "Shortwave downwelling flux"       , "W m-2", "zh", group_name);
        stats.add_prof("sw_flux_dn_dir", "Shortwave direct downwelling flux", "W m-2", "zh", group_name);

        stats.add_prof("sw_heat_dir_rt"    , "Raytraced heating rates from direct radiation"   , "K s-2", "z", group_name);
        stats.add_prof("sw_heat_dif_rt"    , "Raytraced heating rates from diffuse radiation"  , "K s-2", "z", group_name);
        
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

template<typename TF>
void Radiation_rrtmgp_rt<TF>::set_sun_location(Timeloop<TF>& timeloop)
{
    // Update the solar zenith angle.
    const int day_of_year = int(timeloop.calc_day_of_year());
    const int year = timeloop.get_year();
    const TF seconds_after_midnight = TF(timeloop.calc_hour_of_day()*3600);
    this->mu0 = calc_cos_zenith_angle(lat, lon, day_of_year, seconds_after_midnight, year);

    // Calculate correction factor for impact Sun's distance on the solar "constant"
    const TF frac_day_of_year = TF(day_of_year) + seconds_after_midnight / TF(86400);
    this->tsi_scaling = calc_sun_distance_factor(frac_day_of_year);
}

template<typename TF>
void Radiation_rrtmgp_rt<TF>::set_background_column_shortwave(Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();

    const int n_bnd = kdist_sw->get_nband();

    // Set the solar zenith angle and albedo.
    Array<Float,2> sfc_alb_dir({n_bnd, n_col});
    Array<Float,2> sfc_alb_dif({n_bnd, n_col});

    for (int ibnd=1; ibnd<=n_bnd; ++ibnd)
    {
        sfc_alb_dir({ibnd, 1}) = this->sfc_alb_dir;
        sfc_alb_dif({ibnd, 1}) = this->sfc_alb_dif;
    }

    Array<Float,1> mu0({n_col});
    mu0({1}) = this->mu0;

    solve_shortwave_column(
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

#ifndef USECUDA
template<typename TF>
void Radiation_rrtmgp_rt<TF>::exec(
        Thermo<TF>& thermo, const double time, Timeloop<TF>& timeloop, Stats<TF>& stats)
{
    throw std::runtime_error("no raytracing in CPU mode, sorry!");
}


template<typename TF>
std::vector<TF>& Radiation_rrtmgp_rt<TF>::get_surface_radiation(const std::string& name)
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


#ifndef USECUDA
template<typename TF>
void Radiation_rrtmgp_rt<TF>::exec_all_stats(
        Stats<TF>& stats, Cross<TF>& cross,
        Dump<TF>& dump, Column<TF>& column,
        Thermo<TF>& thermo, Timeloop<TF>& timeloop,
        const unsigned long itime, const int iotime)
{
    throw std::runtime_error("no no no, no raytracing in cpu!");
}
#endif


template<typename TF>
void Radiation_rrtmgp_rt<TF>::exec_individual_column_stats(
        Column<TF>& column, Thermo<TF>& thermo, Timeloop<TF>& timeloop, Stats<TF>& stats)
{
    throw std::runtime_error("We are not running column output in raytracing mode!");
}

template<typename TF>
bool Radiation_rrtmgp_rt<TF>::is_day(const Float mu0)
{
    if (mu0 > Constants::mu0_min<Float>)
        return true;

    return false;
}

#ifdef RTE_RRTMGP_SINGLE_PRECISION
template class Radiation_rrtmgp_rt<float>;
#else
template class Radiation_rrtmgp_rt<double>;
#endif
