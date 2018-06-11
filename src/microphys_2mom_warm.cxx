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


    // =============
    // DO MICRO!
    // =============


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
