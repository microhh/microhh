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
#include <cmath>

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

// Constants, move out later.
namespace
{
    template<typename TF> constexpr TF N_0r = 8.e6; // Intercept parameter rain (m-4).
    template<typename TF> constexpr TF N_0s = 3.e6; // Intercept parameter snow (m-4).
    template<typename TF> constexpr TF N_0g = 4.e6; // Intercept parameter graupel (m-4).
}

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

    // Autoconversion.
    template<typename TF>
    void autoconversion(
            TF* const restrict qrt, TF* const restrict qst, TF* const restrict qgt,
            TF* const restrict qtt, TF* const restrict thlt,
            const TF* const restrict qr, const TF* const restrict qs, const TF* const restrict qg,
            const TF* const restrict qt, const TF* const restrict thl,
            const TF* const restrict ql, const TF* const restrict qi,
            const TF* const restrict rho, const TF* const restrict exner,
            const TF N_d,
            const int istart, const int jstart, const int kstart,
            const int iend, const int jend, const int kend,
            const int jj, const int kk)
    {
        using Fast_math::pow2;

        const TF D_d = TF(0.146) - TF(5.964e-2)*std::log(N_d / TF(2.e3));

        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*jj + k*kk;

                    // Compute the T out of the known values of ql and qi, this saves memory and sat_adjust.
                    const TF T = exner[k]*thl[ijk] + Lv<TF>/cp<TF>*ql[ijk] + Ls<TF>/cp<TF>*qi[ijk];

                    constexpr TF gamma_saut = TF(0.025);
                    constexpr TF gamma_gaut = TF(0.09);

                    constexpr TF q_icrt = TF(0.);
                    constexpr TF q_scrt = TF(6.e-4);

                    const TF beta_1 = TF(1.e-3)*std::exp(gamma_saut * (T - T0<TF>));
                    const TF beta_2 = TF(1.e-3)*std::exp(gamma_gaut * (T - T0<TF>));

                    // Calculate the three autoconversion rates.
                    const TF P_raut = TF(16.7)/rho[k] * pow2(rho[k]*ql[ijk]) / (TF(5.) + TF(3.6e-5)*N_d/(D_d*rho[k]*ql[ijk]));
                    const TF P_saut = beta_1*(qi[ijk] - q_icrt);
                    const TF P_gaut = beta_2*(qs[ijk] - q_scrt);

                    // Cloud to rain.
                    qtt[ijk] -= P_raut;
                    qrt[ijk] += P_raut;

                    // Ice to snow.
                    qtt[ijk] -= P_saut;
                    qst[ijk] += P_saut;

                    // Snow to graupel.
                    qst[ijk] -= P_gaut;
                    qgt[ijk] += P_gaut;
                }
    }

    // Accretion.
    /*
    template<typename TF>
    void accretion(
            TF* const restrict qrt, TF* const restrict qst, TF* const restrict qgt,
            TF* const restrict qtt, TF* const restrict thlt,
            const TF* const restrict qr, const TF* const restrict qs, const TF* const restrict qg,
            const TF* const restrict ql, const TF* const restrict qi,
            const TF* const restrict rho, const TF* const restrict exner, const TF nc,
            const int istart, const int jstart, const int kstart,
            const int iend, const int jend, const int kend,
            const int jj, const int kk)
    {
    }
    */
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
    N_d = inputin.get_item<TF>("micro", "Nd", "", 50); // CvH: cm-3 do we need conversion, or do we stick with Tomita?

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

    // Get liquid water, ice and pressure variables before starting.
    auto ql = fields.get_tmp();
    auto qi = fields.get_tmp();

    thermo.get_thermo_field(*ql, "ql", false, false);
    thermo.get_thermo_field(*qi, "qi", false, false);

    const std::vector<TF>& p = thermo.get_p_vector();
    const std::vector<TF>& exner = thermo.get_exner_vector();

    // CLOUD WATER -> RAIN
    autoconversion(
            fields.st.at("qr")->fld.data(), fields.st.at("qs")->fld.data(), fields.st.at("qg")->fld.data(),
            fields.st.at("qt")->fld.data(), fields.st.at("thl")->fld.data(),
            fields.sp.at("qr")->fld.data(), fields.sp.at("qs")->fld.data(), fields.sp.at("qg")->fld.data(),
            fields.sp.at("qt")->fld.data(), fields.sp.at("thl")->fld.data(),
            ql->fld.data(), qi->fld.data(),
            fields.rhoref.data(), exner.data(),
            this->N_d,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.ijcells);

    // CLOUD WATER -> SNOW
    // CLOUD WATER -> GRAUPEL
    // ICE -> SNOW
    // ICE -> GRAUPEL
    // RAIN -> VAPOR
    // RAIN -> SNOW
    // RAIN -> GRAUPEL
    // SNOW <-> VAPOR
    // SNOW -> RAIN
    // SNOW -> GRAUPEL
    // GRAUPEL -> RAIN
    // GRAUPEL <-> VAPOR

    // Release the temp fields and save statistics.
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
