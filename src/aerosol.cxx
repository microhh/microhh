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
#include <algorithm>

#include "aerosol.h"
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
Aerosol<TF>::Aerosol(
    Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    // Read `.ini` settings.
    sw_aerosol = inputin.get_item<bool>("aerosol", "swaerosol", "", 0);
    sw_timedep = inputin.get_item<bool>("aerosol", "swtimedep", "", 0);
}

template <typename TF>
Aerosol<TF>::~Aerosol()
{
}

template <typename TF>
void Aerosol<TF>::init()
{
    // Allocate (`.resize`) arrays.
    if (!sw_aerosol)
        return;

    auto& gd = grid.get_grid_data();

    std::cout << "init() aerosols" << std::endl;

    aermr01.resize(gd.kcells);
    aermr02.resize(gd.kcells);
    aermr03.resize(gd.kcells);
    aermr04.resize(gd.kcells);
    aermr05.resize(gd.kcells);
    aermr06.resize(gd.kcells);
    aermr07.resize(gd.kcells);
    aermr08.resize(gd.kcells);
    aermr09.resize(gd.kcells);
    aermr10.resize(gd.kcells);
    aermr11.resize(gd.kcells);

//    od_wl10.resize(gd.kcells);
//    ssa_wl10.resize(gd.kcells);
//    g_wl10.resize(gd.kcells);
//
//    const std::string& coef_file{"aerosol_optics.nc"};
//    Netcdf_file opt_prop_nc(master, coef_file, Netcdf_mode::Read);
//
//    nwavelengths = opt_prop_nc.get_dimension_size("band_sw");
//    nspecies_phobic = opt_prop_nc.get_dimension_size("hydrophobic");
//    nspecies_philic = opt_prop_nc.get_dimension_size("hydrophilic");
//    nrh = opt_prop_nc.get_dimension_size("relative_humidity");
//
//    mext_phobic.resize(nspecies_phobic*nwavelengths);
//    ssa_phobic.resize(nspecies_phobic*nwavelengths);
//    g_phobic.resize(nspecies_phobic*nwavelengths);
//    mext_philic.resize(nspecies_philic*nrh*nwavelengths);
//    ssa_philic.resize(nspecies_philic*nrh*nwavelengths);
//    g_philic.resize(nspecies_philic*nrh*nwavelengths);
//    rh_classes.resize(nrh);
//
//    mmr.resize(gd.kcells);
//
//    aod_ml.resize(gd.ncells*nwavelengths);
//    ssa_ml.resize(gd.ncells*nwavelengths);
//    g_ml.resize(gd.ncells*nwavelengths);

}

template <typename TF>
void Aerosol<TF>::create(Input& inputin, Netcdf_handle& input_nc, Stats<TF>& stats)
{
    // Read input from NetCDF and prepare statistics output.
    if (!sw_aerosol)
        return;

    auto& gd = grid.get_grid_data();

    std::cout << "create() aerosols" << std::endl;

//    // Read aerosol optical properties from lookup table
//    // Open file
//    const std::string& coef_file{"aerosol_optics.nc"};
//    Netcdf_file opt_prop_nc(master, coef_file, Netcdf_mode::Read);
//    // Load properties
//    rh_classes = opt_prop_nc.get_variable<TF>("relative_humidity2", {nrh});
//    mext_phobic = opt_prop_nc.get_variable<TF>("mass_ext_sw_hydrophobic", {nspecies_phobic, nwavelengths});
//    ssa_phobic = opt_prop_nc.get_variable<TF>("ssa_sw_hydrophobic", {nspecies_phobic, nwavelengths});
//    g_phobic = opt_prop_nc.get_variable<TF>("asymmetry_sw_hydrophobic", {nspecies_phobic, nwavelengths});
//    mext_philic = opt_prop_nc.get_variable<TF>("mass_ext_sw_hydrophilic", {nspecies_philic, nrh, nwavelengths});
//    ssa_philic = opt_prop_nc.get_variable<TF>("ssa_sw_hydrophilic", {nspecies_philic, nrh, nwavelengths});
//    g_philic = opt_prop_nc.get_variable<TF>("asymmetry_sw_hydrophilic", {nspecies_philic, nrh, nwavelengths});

    if (sw_timedep)
    {
        // create time dependent profiles
        const TF offset = 0;
        std::string timedep_dim = "time_aerosols";
        tdep_aermr01 = std::make_unique<Timedep<TF>>(master, grid, "aermr01", inputin.get_item<bool>("aerosol", "swtimedep", "", false));
        tdep_aermr01->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr02 = std::make_unique<Timedep<TF>>(master, grid, "aermr02", inputin.get_item<bool>("aerosol", "swtimedep", "", false));
        tdep_aermr02->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr03 = std::make_unique<Timedep<TF>>(master, grid, "aermr03", inputin.get_item<bool>("aerosol", "swtimedep", "", false));
        tdep_aermr03->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr04 = std::make_unique<Timedep<TF>>(master, grid, "aermr04", inputin.get_item<bool>("aerosol", "swtimedep", "", false));
        tdep_aermr04->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr05 = std::make_unique<Timedep<TF>>(master, grid, "aermr05", inputin.get_item<bool>("aerosol", "swtimedep", "", false));
        tdep_aermr05->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr06 = std::make_unique<Timedep<TF>>(master, grid, "aermr06", inputin.get_item<bool>("aerosol", "swtimedep", "", false));
        tdep_aermr06->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr07 = std::make_unique<Timedep<TF>>(master, grid, "aermr07", inputin.get_item<bool>("aerosol", "swtimedep", "", false));
        tdep_aermr07->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr08 = std::make_unique<Timedep<TF>>(master, grid, "aermr08", inputin.get_item<bool>("aerosol", "swtimedep", "", false));
        tdep_aermr08->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr09 = std::make_unique<Timedep<TF>>(master, grid, "aermr09", inputin.get_item<bool>("aerosol", "swtimedep", "", false));
        tdep_aermr09->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr10 = std::make_unique<Timedep<TF>>(master, grid, "aermr10", inputin.get_item<bool>("aerosol", "swtimedep", "", false));
        tdep_aermr10->create_timedep_prof(input_nc, offset, timedep_dim);
        tdep_aermr11 = std::make_unique<Timedep<TF>>(master, grid, "aermr11", inputin.get_item<bool>("aerosol", "swtimedep", "", false));
        tdep_aermr11->create_timedep_prof(input_nc, offset, timedep_dim);
    }
    else
    {
        // Read NetCDF input.
        Netcdf_group& group_nc = input_nc.get_group("init");
        group_nc.get_variable(aermr01, "aermr01", {0}, {gd.ktot});
        std::rotate(aermr01.rbegin(), aermr01.rbegin()+gd.kstart, aermr01.rend());
        group_nc.get_variable(aermr02, "aermr02", {0}, {gd.ktot});
        std::rotate(aermr02.rbegin(), aermr02.rbegin()+gd.kstart, aermr02.rend());
        group_nc.get_variable(aermr03, "aermr03", {0}, {gd.ktot});
        std::rotate(aermr03.rbegin(), aermr03.rbegin()+gd.kstart, aermr03.rend());
        group_nc.get_variable(aermr04, "aermr04", {0}, {gd.ktot});
        std::rotate(aermr04.rbegin(), aermr04.rbegin()+gd.kstart, aermr04.rend());
        group_nc.get_variable(aermr05, "aermr05", {0}, {gd.ktot});
        std::rotate(aermr05.rbegin(), aermr05.rbegin()+gd.kstart, aermr05.rend());
        group_nc.get_variable(aermr06, "aermr06", {0}, {gd.ktot});
        std::rotate(aermr06.rbegin(), aermr06.rbegin()+gd.kstart, aermr06.rend());
        group_nc.get_variable(aermr07, "aermr07", {0}, {gd.ktot});
        std::rotate(aermr07.rbegin(), aermr07.rbegin()+gd.kstart, aermr07.rend());
        group_nc.get_variable(aermr08, "aermr08", {0}, {gd.ktot});
        std::rotate(aermr08.rbegin(), aermr08.rbegin()+gd.kstart, aermr08.rend());
        group_nc.get_variable(aermr09, "aermr09", {0}, {gd.ktot});
        std::rotate(aermr09.rbegin(), aermr09.rbegin()+gd.kstart, aermr09.rend());
        group_nc.get_variable(aermr10, "aermr10", {0}, {gd.ktot});
        std::rotate(aermr10.rbegin(), aermr10.rbegin()+gd.kstart, aermr10.rend());
        group_nc.get_variable(aermr11, "aermr11", {0}, {gd.ktot});
        std::rotate(aermr11.rbegin(), aermr11.rbegin()+gd.kstart, aermr11.rend());
    }

    // Prepare statistics.
    const std::string group_name = "default";

    if (sw_timedep)
    {
        stats.add_prof("aermr01", "Sea salt (0.03 - 0.5 um) mixing ratio", "kg Kg-1", "z", group_name);
        stats.add_prof("aermr02", "Sea salt (0.5 - 5 um) mixing ratio", "kg Kg-1", "z", group_name);
        stats.add_prof("aermr03", "Sea salt (5 - 20 um) mixing ratio", "kg Kg-1", "z", group_name);
        stats.add_prof("aermr04", "Dust (0.03 - 0.55 um) mixing ratio", "kg Kg-1", "z", group_name);
        stats.add_prof("aermr05", "Dust (0.55 - 0.9 um) mixing ratio", "kg Kg-1", "z", group_name);
        stats.add_prof("aermr06", "Dust (0.9 - 20 um) mixing ratio", "kg Kg-1", "z", group_name);
        stats.add_prof("aermr07", "Organic matter (hydrophilic) mixing ratio", "kg Kg-1", "z", group_name);
        stats.add_prof("aermr08", "Organic matter (hydrophobic) mixing ratio", "kg Kg-1", "z", group_name);
        stats.add_prof("aermr09", "Black carbon (hydrophilic) mixing ratio", "kg Kg-1", "z", group_name);
        stats.add_prof("aermr10", "Black carbon (hydrophobic) mixing ratio", "kg Kg-1", "z", group_name);
        stats.add_prof("aermr11", "Sulfates mixing ratio", "kg Kg-1", "z", group_name);
    }
    else
    {
        stats.add_fixed_prof("aermr01", "Sea salt (0.03 - 0.5 um) mixing ratio", "kg Kg-1", "z", group_name, aermr01);
        stats.add_fixed_prof("aermr02", "Sea salt (0.5 - 5 um) mixing ratio", "kg Kg-1", "z", group_name, aermr02);
        stats.add_fixed_prof("aermr03", "Sea salt (5 - 20 um) mixing ratio", "kg Kg-1", "z", group_name, aermr03);
        stats.add_fixed_prof("aermr04", "Dust (0.03 - 0.55 um) mixing ratio", "kg Kg-1", "z", group_name, aermr04);
        stats.add_fixed_prof("aermr05", "Dust (0.55 - 0.9 um) mixing ratio", "kg Kg-1", "z", group_name, aermr05);
        stats.add_fixed_prof("aermr06", "Dust (0.9 - 20 um) mixing ratio", "kg Kg-1", "z", group_name, aermr06);
        stats.add_fixed_prof("aermr07", "Organic matter (hydrophilic) mixing ratio", "kg Kg-1", "z", group_name, aermr07);
        stats.add_fixed_prof("aermr08", "Organic matter (hydrophobic) mixing ratio", "kg Kg-1", "z", group_name, aermr08);
        stats.add_fixed_prof("aermr09", "Black carbon (hydrophilic) mixing ratio", "kg Kg-1", "z", group_name, aermr09);
        stats.add_fixed_prof("aermr10", "Black carbon (hydrophobic) mixing ratio", "kg Kg-1", "z", group_name, aermr10);
        stats.add_fixed_prof("aermr11", "Sulfates mixing ratio", "kg Kg-1", "z", group_name, aermr11);
    }

//    stats.add_prof("od_wl10", "Optical depth for wavelength band 10", "-", "z", group_name);
//    stats.add_prof("ssa_wl10", "single scattering albedo for wavelength band 10", "-", "z", group_name);
//    stats.add_prof("g_wl10", "asymmetry for wavelength band 10", "-", "z", group_name);
}

#ifndef USECUDA
template <typename TF>
void Aerosol<TF>::exec(Thermo<TF>& thermo)
{
    if (!sw_aerosol)
        return;

    auto& gd = grid.get_grid_data();

//    // get pressure and relative humidity from thermo
//    std::vector<TF> ph_thermo = thermo.get_basestate_vector("ph");
//    auto rh = fields.get_tmp();
//    thermo.get_thermo_field(*rh, "rh", false, false);
//
//    std::vector<TF> local_od;
//    std::vector<TF> od_ext_ml;
//    std::vector<TF> od_scat_ml;
//    std::vector<TF> scatg_ml;
//
//    local_od.resize(gd.kcells*nwavelengths*gd.ijcells);
//    od_ext_ml.resize(gd.kcells*nwavelengths*gd.ijcells);
//    od_scat_ml.resize(gd.kcells*nwavelengths*gd.ijcells);
//    scatg_ml.resize(gd.kcells*nwavelengths*gd.ijcells);
//
//    // Loop over all aerosol types
//    std::string list_aerosol_types[]{"SS1", "SS2", "SS3", "DU1", "DU2", "DU3", "OM1", "OM2", "BC1", "BC2", "SU"};
//    for (auto &&aerosol_class : list_aerosol_types)
//        // loop over all wavelengths
//        for (int n=0; n<nwavelengths; ++n)
//            // loop over all heights
//            for (int k=gd.kstart; k<gd.kend; ++k)
//                // loop over all y
//                for (int j=gd.jstart; j<gd.jend; ++j)
//                    // loop over all x
//                    for (int i=gd.istart; i<gd.iend; ++i)
//                    {
//                        int ijk = i + j*gd.icells + k*gd.ijcells;
//                        int ijkn = i + j*gd.icells + k*gd.ijcells + n*gd.ncells;
//
//                        // determine the relative humidity category
//                        rh_class(rh->fld.data()[ijk], rh_classes);
//
//                        // get the aerosol properties for the combination of aerosol type, wavelength and rh.
//                        // This is done in a separate function, because
//                        // 1. the indexing depends on whether the aerosol type is hydrophilic or hydrophobic
//                        // 2. each type of aerosols has its own set of properties in the lookup tables
//                        set_aerosol_properties(aerosol_class, n, ihum);
//
//                        // calculate the pressure difference divided by g
//                        float dp_g = (ph_thermo[k]-ph_thermo[k+1])/Constants::grav<TF>;
//
//                        // calculate the optical thicknesses and add to the total
//                        local_od[ijkn] = mmr[k] * dp_g * mext;
//                        od_ext_ml[ijkn] += local_od[ijkn];
//                        od_scat_ml[ijkn] += local_od[ijkn] * ssa;
//                        scatg_ml[ijkn] += local_od[ijkn] * ssa * g;
//                    }
//
//    fields.release_tmp(rh);
//
//    // loop over the od values to compute aod, g, and ssa
//    for (int idx = 0; idx<nwavelengths*gd.ncells; ++idx)
//    {
//        aod_ml[idx] = od_ext_ml[idx];
//        ssa_ml[idx] = od_scat_ml[idx] / od_ext_ml[idx];
//        g_ml[idx] = scatg_ml[idx] / od_scat_ml[idx];
//    }
//
//    std::cout << "exec() aerosols" << std::endl;
}

template <typename TF>
void Aerosol<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (!sw_aerosol)
        return;

    if (sw_timedep)
    {
        tdep_aermr01 ->update_time_dependent_prof(aermr01, timeloop);
        tdep_aermr02 ->update_time_dependent_prof(aermr02, timeloop);
        tdep_aermr03 ->update_time_dependent_prof(aermr03, timeloop);
        tdep_aermr04 ->update_time_dependent_prof(aermr04, timeloop);
        tdep_aermr05 ->update_time_dependent_prof(aermr05, timeloop);
        tdep_aermr06 ->update_time_dependent_prof(aermr06, timeloop);
        tdep_aermr07 ->update_time_dependent_prof(aermr07, timeloop);
        tdep_aermr08 ->update_time_dependent_prof(aermr08, timeloop);
        tdep_aermr09 ->update_time_dependent_prof(aermr09, timeloop);
        tdep_aermr10 ->update_time_dependent_prof(aermr10, timeloop);
        tdep_aermr11 ->update_time_dependent_prof(aermr11, timeloop);

//        std::cout << "update_time_dependent() aerosols" << std::endl;
//        throw std::runtime_error("Time dependent aerosols are not (yet) supported!\n");
    }
}
#endif

template<typename TF>
void Aerosol<TF>::exec_stats(Stats<TF>& stats)
{
    if (!sw_aerosol)
        return;

    auto& gd = grid.get_grid_data();

//    // determine the optical depth at one wavelength, one xy, all heights to write to output
//    int n = 9;
//    int j = gd.jstart + gd.jtot/2;
//    int i = gd.istart + gd.itot/2;
//    for (int k=gd.kstart; k<gd.kend; ++k)
//    {
//        int ijkn = i + j * gd.icells + k * gd.ijcells + n * gd.ncells;
//        od_wl10[k] = aod_ml[ijkn];
//        ssa_wl10[k] = ssa_ml[ijkn];
//        g_wl10[k] = g_ml[ijkn];
//    }
//
//    stats.set_prof("od_wl10", od_wl10);
//    stats.set_prof("ssa_wl10", ssa_wl10);
//    stats.set_prof("g_wl10", g_wl10);

    if (sw_timedep)
    {
        stats.set_prof("aermr01", aermr01);
        stats.set_prof("aermr02", aermr02);
        stats.set_prof("aermr03", aermr03);
        stats.set_prof("aermr04", aermr04);
        stats.set_prof("aermr05", aermr05);
        stats.set_prof("aermr06", aermr06);
        stats.set_prof("aermr07", aermr07);
        stats.set_prof("aermr08", aermr08);
        stats.set_prof("aermr09", aermr09);
        stats.set_prof("aermr10", aermr10);
        stats.set_prof("aermr11", aermr11);

        // std::cout << "exec_stats() aerosols" << std::endl;
    }
}

//template<typename TF>
//void Aerosol<TF>::set_aerosol_properties(const std::string& aerosol_class, int wavelength, int humidity)
//{
//    if (aerosol_class == "DU1")
//    {
//        set_hydrophobic(0, wavelength);
//        mmr = aermr04;
//    }
//    else if (aerosol_class == "DU2")
//    {
//        set_hydrophobic(7, wavelength);
//        mmr = aermr05;
//    }
//    else if (aerosol_class == "DU3")
//    {
//        set_hydrophobic(5, wavelength);
//        mmr = aermr06;
//    }
//    else if (aerosol_class == "BC1")
//    {
//        set_hydrophobic(10, wavelength);
//        mmr = aermr09;
//    }
//    else if (aerosol_class == "BC2")
//    {
//        set_hydrophobic(10, wavelength);
//        mmr = aermr10;
//    }
//    else if (aerosol_class == "SS1")
//    {
//        set_hydrophilic(0, humidity, wavelength);
//        mmr = aermr01;
//    }
//    else if (aerosol_class == "SS2")
//    {
//        set_hydrophilic(1, humidity, wavelength);
//        mmr = aermr02;
//    }
//    else if (aerosol_class == "SS3")
//    {
//        set_hydrophilic(2, humidity, wavelength);
//        mmr = aermr03;
//    }
//    else if (aerosol_class == "SU")
//    {
//        set_hydrophilic(4, humidity, wavelength);
//        mmr = aermr11;
//    }
//    else if (aerosol_class == "OM1")
//    {
//        set_hydrophobic(9, wavelength);
//        mmr = aermr08;
//    }
//    else if (aerosol_class == "OM2")
//    {
//        set_hydrophilic(3, humidity, wavelength);
//        mmr = aermr07;
//    }
//}
//
//template<typename TF>
//void Aerosol<TF>::set_hydrophobic(int i, int wavelength)
//{
//    mext = mext_phobic[i*nwavelengths+wavelength];
//    ssa = ssa_phobic[i*nwavelengths+wavelength];
//    g = g_phobic[i*nwavelengths+wavelength];
//}
//
//template<typename TF>
//void Aerosol<TF>::set_hydrophilic(int i, int humidity, int wavelength)
//{
//    mext = mext_philic[i*nwavelengths*nrh+humidity*nwavelengths+wavelength];
//    ssa = ssa_philic[i*nwavelengths*nrh+humidity*nwavelengths+wavelength];
//    g = g_philic[i*nwavelengths*nrh+humidity*nwavelengths+wavelength];
//}
//
//template<typename TF>
//void Aerosol<TF>::rh_class(float rel_hum, std::vector<TF> rh2)
//{
//    {
//        ihum = 0;
//        float rh_class = rh2[ihum];
//        while (rh_class < rel_hum)
//        {
//            ihum += 1;
//            rh_class = rh2[ihum];
//        }
//    }
//}

template<typename TF>
void Aerosol<TF>::get_radiation_fields(std::vector<TF>& mr01, std::vector<TF>& mr02, std::vector<TF>& mr03, std::vector<TF>& mr04,
                                       std::vector<TF>& mr05, std::vector<TF>& mr06, std::vector<TF>& mr07, std::vector<TF>& mr08,
                                       std::vector<TF>& mr09, std::vector<TF>& mr10, std::vector<TF>& mr11)
{
    auto& gd = grid.get_grid_data();

    for (int k=gd.kstart; k<gd.kend; ++k)
            {
                const int k_nogc = k-gd.kgc;

                mr01[k_nogc] = aermr01[k];
                mr02[k_nogc] = aermr02[k];
                mr03[k_nogc] = aermr03[k];
                mr04[k_nogc] = aermr04[k];
                mr05[k_nogc] = aermr05[k];
                mr06[k_nogc] = aermr06[k];
                mr07[k_nogc] = aermr07[k];
                mr08[k_nogc] = aermr08[k];
                mr09[k_nogc] = aermr09[k];
                mr10[k_nogc] = aermr10[k];
                mr11[k_nogc] = aermr11[k];
            }
}

template class Aerosol<double>;
template class Aerosol<float>;
