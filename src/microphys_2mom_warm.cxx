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

#include <iostream>
#include <cstdio>
#include <cmath>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"

#include "constants.h"
#include "microphys.h"
#include "microphys_2mom_warm.h"

namespace
{
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
    TF* get_tmp_slice(std::vector<std::shared_ptr<Field3d<TF>>> &tmp_fields, int &slice_counter,
                      const int jcells, const int ikcells)
    {
        const int tmp_index   = slice_counter / jcells;     // Which tmp field in tmp_fields vector?
        const int fld_index   = slice_counter % jcells;     // Which slice in tmp field?
        const int slice_start = fld_index * ikcells;        // Start index of slice

        slice_counter++;

        return &(tmp_fields[tmp_index]->fld[slice_start]);
    }
}

// Microphysics calculated over entire 3D field
namespace mp3d
{
    // Autoconversion: formation of rain drop by coagulating cloud droplets
    template<typename TF>
    void autoconversion(TF* const restrict qrt, TF* const restrict nrt,
                        TF* const restrict qtt, TF* const restrict thlt,
                        const TF* const restrict qr,  const TF* const restrict ql,
                        const TF* const restrict rho, const TF* const restrict exner,
                        const int istart, const int jstart, const int kstart,
                        const int iend,   const int jend,   const int kend,
                        const int jj, const int kk,
                        Micro_2mom_warm_constants<TF> micro_constans)
    {
        const TF x_star = 2.6e-10;       // SB06, list of symbols, same as UCLA-LES
        const TF k_cc   = 9.44e9;        // UCLA-LES (Long, 1974), 4.44e9 in SB06, p48
        const TF nu_c   = 1;             // SB06, Table 1., same as UCLA-LES
        //const TF kccxs  = k_cc / (TF(20.) * x_star) * (nu_c+2)*(nu_c+4) / pow(nu_c+1, 2);

        //for (int k=kstart; k<kend; k++)
        //    for (int j=jstart; j<jend; j++)
        //        #pragma ivdep
        //        for (int i=istart; i<iend; i++)
        //        {
        //            const int ijk = i + j*jj + k*kk;
        //            if(ql[ijk] > ql_min<TF>)
        //            {
        //                const TF xc      = rho[k] * ql[ijk] / Nc0<TF>; // Mean mass of cloud drops [kg]
        //                const TF tau     = TF(1.) - ql[ijk] / (ql[ijk] + qr[ijk] + dsmall); // SB06, Eq 5
        //                const TF phi_au  = TF(600.) * pow(tau, TF(0.68)) * pow(TF(1.) - pow(tau, TF(0.68)), 3); // UCLA-LES
        //                //const TF phi_au  = 400. * pow(tau, 0.7) * pow(1. - pow(tau, 0.7), 3); // SB06, Eq 6
        //                const TF au_tend = rho[k] * kccxs * pow(ql[ijk], 2) * pow(xc, 2) *
        //                                       (TF(1.) + phi_au / pow(TF(1.)-tau, 2)); // SB06, eq 4

        //                qrt[ijk]  += au_tend;
        //                nrt[ijk]  += au_tend * rho[k] / x_star;
        //                qtt[ijk]  -= au_tend;
        //                thlt[ijk] += Lv / (cp * exner[k]) * au_tend;
        //            }
        //        }
    }
}

// Microphysics calculated over 2D (xz) slices
namespace mp2d
{

}


template<typename TF>
Microphys_2mom_warm<TF>::Microphys_2mom_warm(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    Microphys<TF>(masterin, gridin, fieldsin, inputin)
{
    swmicrophys = Microphys_type::Warm_2mom;

    #ifdef USECUDA
    throw std::runtime_error("swmicro = \"2mom_warm\" not (yet) implemented in CUDA\n");
    #endif

    // Read microphysics switches and settings
    swmicrobudget = inputin.get_item<bool>("micro", "swmicrobudget", "", false);
    cflmax        = inputin.get_item<TF>("micro", "cflmax", "", 2.);

    // Initialize the qr (rain water specific humidity) and nr (droplot number concentration) fields
    fields.init_prognostic_field("qr", "Rain water mixing ratio", "kg kg-1");
    fields.init_prognostic_field("nr", "Number density rain", "m-3");
}

template<typename TF>
Microphys_2mom_warm<TF>::~Microphys_2mom_warm()
{
}

template<typename TF>
void Microphys_2mom_warm<TF>::exec(Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();

    // Remove spurious negative values from qr and nr fields
    remove_negative_values(fields.sp.at("qr")->fld.data(), gd.istart, gd.jstart, gd.kstart,
                           gd.iend, gd.jend, gd.kend, gd.icells, gd.ijcells);
    remove_negative_values(fields.sp.at("nr")->fld.data(), gd.istart, gd.jstart, gd.kstart,
                           gd.iend, gd.jend, gd.kend, gd.icells, gd.ijcells);

    // Get cloud liquid water specific humidity from thermodynamics
    auto ql = fields.get_tmp();
    thermo.get_thermo_field(*ql, "ql", false, false);

    // Microphysics is handled in XZ slices, to
    // (1) limit the required number of tmp fields
    // (2) re-use some expensive calculations used in multiple microphysics routines.
    const int ikcells    = gd.icells * gd.kcells;                           // Size of XZ slice
    const int n_slices   = 12;                                              // Number of XZ slices required
    const int n_tmp_flds = std::ceil(static_cast<TF>(n_slices)/gd.jcells);  // Number of required tmp fields

    // Load the required number of tmp fields:
    std::vector<std::shared_ptr<Field3d<TF>>> tmp_fields;
    for (int n=0; n<n_tmp_flds; ++n)
        tmp_fields.push_back(fields.get_tmp());

    // Get pointers to slices in tmp fields:
    int slice_counter = 0;

    TF* w_qr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
    TF* w_nr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

    TF* c_qr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
    TF* c_nr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

    TF* slope_qr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
    TF* slope_nr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

    TF* flux_qr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
    TF* flux_nr = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

    TF* rain_mass = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
    TF* rain_diam = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

    TF* lambda_r = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);
    TF* mu_r     = get_tmp_slice<TF>(tmp_fields, slice_counter, gd.jcells, ikcells);

    // ---------------------------------
    // Calculate microphysics tendencies
    // ---------------------------------

    // Autoconversion; formation of rain drop by coagulating cloud droplets
    //autoconversion(fields.st.at("qr")->fld.data(), fields.st.at("nr")->fld.data(), fields.st.at("qt")->fld.data(), fields.st.at("thl")->fld.data(),
    //               fields.sp.at("qr")->fld.data(), ql->fld.data(), fields.rhoref.data(), bs.exnref.data(),
    //               gd.istart, gd.jstart, gd.kstart,
    //               gd.iend,   gd.jend,   gd.kend,
    //               gd.icells, gd.ijcells,
    //               micro_constants);

    // Accretion; growth of raindrops collecting cloud droplets
    //accretion(fields.st.at("qr")->fld.data(), fields.st.at("qt")->fld.data(), fields.st.at("thl")->fld.data(),
    //          fields.sp.at("qr")->fld.data(), ql->fld.data(), fields.rhoref.data(), bs.exnref.data(),
    //          gd.istart, gd.jstart, gd.kstart,
    //          gd.iend,   gd.jend,   gd.kend,
    //          gd.icells, gd.ijcells);


    // Release all local tmp fields in use
    for (auto& it: tmp_fields)
        fields.release_tmp(it);

    fields.release_tmp(ql);

    throw 1;
}

template<typename TF>
void Microphys_2mom_warm<TF>::exec_stats(Stats<TF>& stats, std::string mask_name,
                                         Field3d<TF>& mask_field, Field3d<TF>& mask_fieldh, const double dt)
{
}

template<typename TF>
void Microphys_2mom_warm<TF>::exec_dump(Dump<TF>& dump, unsigned long iotime)
{
}

template<typename TF>
void Microphys_2mom_warm<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
}

template<typename TF>
unsigned long Microphys_2mom_warm<TF>::get_time_limit(unsigned long idt, const double dt)
{
    return Constants::ulhuge;
}

template<typename TF>
bool Microphys_2mom_warm<TF>::has_mask(std::string name)
{
    return false;
}

template<typename TF>
void Microphys_2mom_warm<TF>::get_mask(Field3d<TF>& mfield, Field3d<TF>& mfieldh, Stats<TF>& stats, std::string mask_name)
{
}

template class Microphys_2mom_warm<double>;
template class Microphys_2mom_warm<float>;
