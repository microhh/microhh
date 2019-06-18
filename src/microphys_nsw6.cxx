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

#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "diff.h"
#include "stats.h"
#include "cross.h"
#include "thermo.h"
#include "thermo_moist_functions.h"

#include "constants.h"
#include "microphys.h"
#include "microphys_nsw6.h"


namespace
{
    using namespace Constants;
    using namespace Thermo_moist_functions;

    template<typename TF>
    void remove_negative_values(TF* const restrict field,
                                const int istart, const int jstart, const int kstart,
                                const int iend,   const int jend,   const int kend,
                                const int jj,     const int kk)
    {
        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;
                    field[ijk] = std::max(TF(0.), field[ijk]);
                }
    }

    template<typename TF>
    void zero_field(TF* const restrict field, const int ncells)
    {
        for (int n=0; n<ncells; n++)
            field[n] = TF(0.);
    }

}

template<typename TF>
Microphys_nsw6<TF>::Microphys_nsw6(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Microphys<TF>(masterin, gridin, fieldsin, inputin)
{
    auto& gd = grid.get_grid_data();
    swmicrophys = Microphys_type::Nsw6;

    // Read microphysics switches and settings
    // swmicrobudget = inputin.get_item<bool>("micro", "swmicrobudget", "", false);
    // cflmax        = inputin.get_item<TF>("micro", "cflmax", "", 2.);
    // Nc0<TF>       = inputin.get_item<TF>("micro", "Nc0", "", 70e6);

    // Initialize the qr (rain water specific humidity) and nr (droplot number concentration) fields
    fields.init_prognostic_field("qr", "Rain water specific humidity", "kg kg-1", gd.sloc);
    fields.init_prognostic_field("qg", "Graupel specific humidity", "kg kg-1", gd.sloc);
    fields.init_prognostic_field("qs", "Snow specific humidity", "kg kg-1", gd.sloc);

    // Load the viscosity for both fields.
    fields.sp.at("qr")->visc = inputin.get_item<TF>("fields", "svisc", "qr");
    fields.sp.at("qg")->visc = inputin.get_item<TF>("fields", "svisc", "qg");
    fields.sp.at("qs")->visc = inputin.get_item<TF>("fields", "svisc", "qs");
}

template<typename TF>
Microphys_nsw6<TF>::~Microphys_nsw6()
{
}

template<typename TF>
void Microphys_nsw6<TF>::init()
{
    auto& gd = grid.get_grid_data();
}

template<typename TF>
void Microphys_nsw6<TF>::create(Input& inputin, Netcdf_handle& input_nc, Stats<TF>& stats, Cross<TF>& cross, Dump<TF>& dump)
{
    const std::string group_name = "thermo";

    /*
    // Add variables to the statistics
    if (stats.get_switch())
    {
        // Time series
        stats.add_time_series("rr", "Mean surface rain rate", "kg m-2 s-1", group_name);
        stats.add_profs(*fields.sp.at("qr"), "z", {"frac", "path", "cover"}, group_name);

        if (swmicrobudget)
        {
            // Microphysics tendencies for qr, nr, thl and qt
            stats.add_prof("sed_qrt", "Sedimentation tendency of qr", "kg kg-1 s-1", "z", group_name);
            stats.add_prof("sed_nrt", "Sedimentation tendency of nr", "m3 s-1", "z", group_name);

            stats.add_prof("auto_qrt" , "Autoconversion tendency qr",  "kg kg-1 s-1", "z", group_name);
            stats.add_prof("auto_nrt" , "Autoconversion tendency nr",  "m-3 s-1", "z", group_name);
            stats.add_prof("auto_thlt", "Autoconversion tendency thl", "K s-1", "z", group_name);
            stats.add_prof("auto_qtt" , "Autoconversion tendency qt",  "kg kg-1 s-1", "z", group_name);

            stats.add_prof("evap_qrt" , "Evaporation tendency qr",  "kg kg-1 s-1", "z", group_name);
            stats.add_prof("evap_nrt" , "Evaporation tendency nr",  "m-3 s-1", "z", group_name);
            stats.add_prof("evap_thlt", "Evaporation tendency thl", "K s-1", "z", group_name);
            stats.add_prof("evap_qtt" , "Evaporation tendency qt",  "kg kg-1 s-1", "z", group_name);

            stats.add_prof("scbr_nrt" , "Selfcollection and breakup tendency nr", "m-3 s-1", "z", group_name);

            stats.add_prof("accr_qrt" , "Accretion tendency qr",  "kg kg-1 s-1", "z", group_name);
            stats.add_prof("accr_thlt", "Accretion tendency thl", "K s-1", "z", group_name);
            stats.add_prof("accr_qtt" , "Accretion tendency qt",  "kg kg-1 s-1", "z", group_name);
        }

        stats.add_tendency(*fields.st.at("thl"), "z", tend_name, tend_longname);
        stats.add_tendency(*fields.st.at("qt") , "z", tend_name, tend_longname);
        stats.add_tendency(*fields.st.at("qr") , "z", tend_name, tend_longname);
        stats.add_tendency(*fields.st.at("nr") , "z", tend_name, tend_longname);
    }

    // Create cross sections
    // 1. Variables that this class can calculate/provide:
    std::vector<std::string> allowed_crossvars = {"rr_bot"};
    // 2. Cross-reference with the variables requested in the .ini file:
    crosslist = cross.get_enabled_variables(allowed_crossvars);
    */
}

#ifndef USECUDA
template<typename TF>
void Microphys_nsw6<TF>::exec(Thermo<TF>& thermo, const double dt, Stats<TF>& stats)
{
    auto& gd = grid.get_grid_data();

    // Remove spurious negative values from qr and nr fields
    remove_negative_values(
            fields.sp.at("qr")->fld.data(), gd.istart, gd.jstart, gd.kstart,
            gd.iend, gd.jend, gd.kend, gd.icells, gd.ijcells);

    // Get cloud liquid water specific humidity from thermodynamics
    auto ql = fields.get_tmp();
    auto qi = fields.get_tmp();
    thermo.get_thermo_field(*ql, "ql", false, false);
    thermo.get_thermo_field(*qi, "qi", false, false);

    // Get pressure and exner function from thermodynamics
    const std::vector<TF>& p     = thermo.get_p_vector();
    const std::vector<TF>& exner = thermo.get_exner_vector();

    fields.release_tmp(ql);
    fields.release_tmp(qi);

    stats.calc_tend(*fields.st.at("thl"), tend_name);
    stats.calc_tend(*fields.st.at("qt"),  tend_name);
    stats.calc_tend(*fields.st.at("qr"),  tend_name);
    stats.calc_tend(*fields.st.at("qg"),  tend_name);
    stats.calc_tend(*fields.st.at("qs"),  tend_name);
}
#endif

template<typename TF>
void Microphys_nsw6<TF>::exec_stats(Stats<TF>& stats, Thermo<TF>& thermo, const double dt)
{
}

template<typename TF>
void Microphys_nsw6<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
}

#ifndef USECUDA
template<typename TF>
unsigned long Microphys_nsw6<TF>::get_time_limit(unsigned long idt, const double dt)
{
    // return idt * cflmax / cfl;
    return Constants::ulhuge;
}
#endif

template<typename TF>
bool Microphys_nsw6<TF>::has_mask(std::string name)
{
    return false;
}

template<typename TF>
void Microphys_nsw6<TF>::get_mask(Stats<TF>& stats, std::string mask_name)
{
    auto& gd = grid.get_grid_data();

    std::string message = "Double moment warm microphysics can not provide mask: \"" + mask_name +"\"";
    throw std::runtime_error(message);
}

template class Microphys_nsw6<double>;
template class Microphys_nsw6<float>;
