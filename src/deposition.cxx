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

//#include <cstdio>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <math.h>
#include <iomanip>
#include <utility>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "stats.h"
#include "netcdf_interface.h"
#include "constants.h"
#include "timeloop.h"
#include "deposition.h"
#include "boundary_surface_lsm.h"
#include "boundary.h"
#include "cross.h"


namespace
{
    template<typename TF>
    void calc_tiled_mean(
            TF* const restrict fld,
            const TF* const restrict f_veg,
            const TF* const restrict f_soil,
            const TF* const restrict f_wet,
            const TF* const restrict fld_veg,
            const TF* const restrict fld_soil,
            const TF* const restrict fld_wet,
            const TF fac,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;

                fld[ij] = (
                    f_veg [ij] * fld_veg [ij] +
                    f_soil[ij] * fld_soil[ij] +
                    f_wet [ij] * fld_wet [ij] ) * fac;
            }
    }


    template<typename TF>
    void calc_spatial_avg_deposition(
            TF* const restrict fld,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        //Calculate sum and count
        TF n_dep = (TF)0.0;
        TF sum_dep = (TF)0.0;

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij = i + j*icells;

                sum_dep += fld[ij];
                n_dep += 1.0;
            }

        // Calculate average
        TF avg_dep = sum_dep / n_dep;

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij = i + j*icells;

                fld[ij] = avg_dep;
            }
    }


    template<typename TF>
    void calc_vd_water(
            TF* const restrict fld,
            const TF* const restrict ra,
            const TF* const restrict ustar,
            const int* const restrict water_mask,
            const TF diff_scl,
            const TF rwat,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
    for (int j=jstart; j<jend; ++j)
        #pragma ivdep
        for (int i=istart; i<iend; ++i)
        {

            const int ij = i + j*icells;

            if (water_mask[ij] == 1)
            {
                const TF ckarman = 0.4;
                const TF rb = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl;

                fld[ij] = (TF)1.0 / (ra[ij] + rb + rwat);
            }

        }

    }

    template<typename TF>
    void calc_deposition_per_tile(
        const std::basic_string<char> lu_type,
        TF* restrict vdo3,
        TF* restrict vdno,
        TF* restrict vdno2,
        TF* restrict vdhno3,
        TF* restrict vdh2o2,
        TF* restrict vdrooh,
        TF* restrict vdhcho,
        const TF* const restrict lai,
        const TF* const restrict c_veg,
        const TF* const restrict rs,
        const TF* const restrict rs_veg,
        const TF* const restrict ra,
        const TF* const restrict ustar,
        const TF* const restrict fraction,
        const TF* const restrict rmes,
        const TF* const restrict rsoil,
        const TF* const restrict rcut,
        const TF* const restrict rws,
        const TF* const restrict rwat,
        const TF* const restrict diff_scl,
        const int istart, const int iend,
        const int jstart, const int jend,
        const int jj)
    {

        const int ntrac_vd = 7;
        const TF ckarman = (TF)0.4;
        const TF hc = (TF)10.0; // constant for now...

        if (lu_type == "veg")
        {

            // Note: I think memory-wise it's more efficient to first loop over ij and then over species,
            // because otherwise rb and rc vectors must be allocated for the entire grid instead of for
            // the number of tracers. Also, it avoids the use of if statements (e.g. "if (t==0) vdo3[ij] = ...")
            std::vector<TF> rmes_local = {rmes[0], rmes[1], rmes[2], rmes[3], rmes[4], rmes[5], rmes[6]};
            std::vector<TF> rb(ntrac_vd, (TF)0.0);
            std::vector<TF> rc(ntrac_vd, (TF)0.0);

            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i) {

                    const int ij = i + j*jj;

                    //Do not proceed in loop if tile fraction is small
                    if (fraction[ij] < (TF)1e-12)
                        continue;

                    //rmes for NO and NO2 requires multiplication with rs, according to Ganzeveld et al. (1995)
                    rmes_local[1] = rmes[1] * rs[ij];
                    rmes_local[2] = rmes[2] * rs[ij];
                    const TF ra_inc = (TF)14. * hc * lai[ij] / ustar[ij];

                    for (int t=0; t<ntrac_vd; ++t)
                    {
                        rb[t] = TF(2.0) / (ckarman * ustar[ij]) * diff_scl[t];
                        rc[t] = TF(1.0) / ((TF)1.0 / (diff_scl[t] + rs[ij] + rmes_local[t]) + (TF)1.0 / rcut[t] + (TF)1.0 / (ra_inc + rsoil[t]));
                    }

                    vdo3[ij]   = (TF)1.0 / (ra[ij] + rb[0] + rc[0]);
                    vdno[ij]   = (TF)1.0 / (ra[ij] + rb[1] + rc[1]);
                    vdno2[ij]  = (TF)1.0 / (ra[ij] + rb[2] + rc[2]);
                    vdhno3[ij] = (TF)1.0 / (ra[ij] + rb[3] + rc[3]);
                    vdh2o2[ij] = (TF)1.0 / (ra[ij] + rb[4] + rc[4]);
                    vdrooh[ij] = (TF)1.0 / (ra[ij] + rb[5] + rc[5]);
                    vdhcho[ij] = (TF)1.0 / (ra[ij] + rb[6] + rc[6]);
                }

        }
        else if (lu_type == "soil")
        {
            std::vector<TF> rb(ntrac_vd, (TF)0.0);

            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i) {

                    const int ij = i + j*jj;

                    //Do not proceed in loop if tile fraction is small
                    if (fraction[ij] < (TF)1e-12) continue;

                    for (int t=0; t<ntrac_vd; ++t)
                    {
                        rb[t] = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl[t];
                    }

                    vdo3[ij]   = (TF)1.0 / (ra[ij] + rb[0] + rsoil[0]);
                    vdno[ij]   = (TF)1.0 / (ra[ij] + rb[1] + rsoil[1]);
                    vdno2[ij]  = (TF)1.0 / (ra[ij] + rb[2] + rsoil[2]);
                    vdhno3[ij] = (TF)1.0 / (ra[ij] + rb[3] + rsoil[3]);
                    vdh2o2[ij] = (TF)1.0 / (ra[ij] + rb[4] + rsoil[4]);
                    vdrooh[ij] = (TF)1.0 / (ra[ij] + rb[5] + rsoil[5]);
                    vdhcho[ij] = (TF)1.0 / (ra[ij] + rb[6] + rsoil[6]);
                }
        }
        else if (lu_type == "wet")
        {
            std::vector<TF> rb_veg(ntrac_vd, (TF)0.0);
            std::vector<TF> rb_soil(ntrac_vd, (TF)0.0);
            std::vector<TF> rc(ntrac_vd, (TF)0.0);
            std::vector<TF> rmes_local = {rmes[0], rmes[1], rmes[2], rmes[3], rmes[4], rmes[5], rmes[6]};

            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jj;

                    // Do not proceed in loop if tile fraction is small
                    if (fraction[ij] < (TF)1e-12) continue;

                    // rmes for NO and NO2 requires multiplication with rs, according to Ganzeveld et al. (1995)
                    // rmes_local[0] =
                    rmes_local[1] = rmes[1] * rs[ij];
                    rmes_local[2] = rmes[2] * rs[ij];
                    const TF ra_inc = (TF)14. * hc * lai[ij] / ustar[ij];

                    //Note that in rc calculation, rcut is replaced by rws for calculating wet skin uptake
                    for (int t=0; t<ntrac_vd; ++t)
                    {
                        rb_veg[t] = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl[t];
                        rb_soil[t] = (TF)1.0 / (ckarman * ustar[ij]) * diff_scl[t];
                        rc[t] = TF(1.0) / ((TF)1.0 / (diff_scl[t] + rs_veg[ij] + rmes_local[t]) + (TF)1.0 / rws[t] + (TF)1.0 / (ra_inc + rsoil[t]));
                    }

                    // Calculate vd for wet skin tile as the weighted average of vd to wet soil and to wet vegetation
                    vdo3[ij]   = c_veg[ij] / (ra[ij] + rb_veg[0] + rc[0]) + ((TF)1.0 - c_veg[ij]) / (ra[ij] + rb_soil[0] + rsoil[0]);
                    vdno[ij]   = c_veg[ij] / (ra[ij] + rb_veg[1] + rc[1]) + ((TF)1.0 - c_veg[ij]) / (ra[ij] + rb_soil[1] + rsoil[1]);
                    vdno2[ij]  = c_veg[ij] / (ra[ij] + rb_veg[2] + rc[2]) + ((TF)1.0 - c_veg[ij]) / (ra[ij] + rb_soil[2] + rsoil[2]);
                    vdhno3[ij] = c_veg[ij] / (ra[ij] + rb_veg[3] + rc[3]) + ((TF)1.0 - c_veg[ij]) / (ra[ij] + rb_soil[3] + rsoil[3]);
                    vdh2o2[ij] = c_veg[ij] / (ra[ij] + rb_veg[4] + rc[4]) + ((TF)1.0 - c_veg[ij]) / (ra[ij] + rb_soil[4] + rsoil[4]);
                    vdrooh[ij] = c_veg[ij] / (ra[ij] + rb_veg[5] + rc[5]) + ((TF)1.0 - c_veg[ij]) / (ra[ij] + rb_soil[5] + rsoil[5]);
                    vdhcho[ij] = c_veg[ij] / (ra[ij] + rb_veg[6] + rc[6]) + ((TF)1.0 - c_veg[ij]) / (ra[ij] + rb_soil[6] + rsoil[6]);
                }
        }
    }
}

template<typename TF>
Deposition<TF>::Deposition(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    sw_deposition = inputin.get_item<bool>("deposition", "swdeposition", "", false);
}


template <typename TF>
Deposition<TF>::~Deposition()
{
}


template <typename TF>
void Deposition<TF>::init(Input& inputin)
{
    // Always read the default deposition velocities. They are needed by 
    // chemistry, even if deposition is disabled.
    vd_o3   = inputin.get_item<TF>("deposition", "vdo3", "", (TF)0.005);
    vd_no   = inputin.get_item<TF>("deposition", "vdno", "", (TF)0.002);
    vd_no2  = inputin.get_item<TF>("deposition", "vdno2", "", (TF)0.005);
    vd_hno3 = inputin.get_item<TF>("deposition", "vdhno3", "", (TF)0.040);
    vd_h2o2 = inputin.get_item<TF>("deposition", "vdh2o2", "", (TF)0.018);
    vd_rooh = inputin.get_item<TF>("deposition", "vdrooh", "", (TF)0.008);
    vd_hcho = inputin.get_item<TF>("deposition", "vdhcho", "", (TF)0.0033);

    if (!sw_deposition)
        return;

    auto& gd = grid.get_grid_data();

    // Create surface tiles for deposition:
    for (auto& name : deposition_tile_names)
        deposition_tiles.emplace(name, Deposition_tile<TF>{});

    for (auto& tile : deposition_tiles)
    {
        tile.second.vdo3.resize(gd.ijcells);
        tile.second.vdno.resize(gd.ijcells);
        tile.second.vdno2.resize(gd.ijcells);
        tile.second.vdno2.resize(gd.ijcells);
        tile.second.vdhno3.resize(gd.ijcells);
        tile.second.vdh2o2.resize(gd.ijcells);
        tile.second.vdrooh.resize(gd.ijcells);
        tile.second.vdhcho.resize(gd.ijcells);
    }

    deposition_tiles.at("veg" ).long_name = "vegetation";
    deposition_tiles.at("soil").long_name = "bare soil";
    deposition_tiles.at("wet" ).long_name = "wet skin";
    deposition_var = inputin.get_item<TF>("deposition", "deposition_var","", (TF)1e5);

    henry_so2 = inputin.get_item<TF>("deposition", "henry_so2", "", (TF)1e5);
    rsoil_so2 = inputin.get_item<TF>("deposition", "rsoil_so2", "", (TF)250.0);
    rwat_so2 = inputin.get_item<TF>("deposition", "rwat_so2", "", (TF)1.0);
    rws_so2 = inputin.get_item<TF>("deposition", "rws_so2", "", (TF)100.0);

    // Note: rmes for NO and NO2 (indices 1 and 2) will still be scaled with rs
    rmes     = {(TF)1.0, (TF)5.0, (TF)0.5, (TF)0.0, (TF)0.0, (TF)0.0, (TF)0.0};
    rsoil    = {(TF)400.0, (TF)1e5, (TF)600.0, (TF)0.0, (TF)0.0, (TF)0.0, (TF)0.0};
    rcut     = {(TF)1e5, (TF)1e5, (TF)1e5, (TF)0.0, (TF)0.0, (TF)0.0, (TF)0.0};
    rws      = {(TF)2000.0, (TF)1e5, (TF)1e5, (TF)0.0, (TF)0.0, (TF)0.0, (TF)0.0};
    rwat     = {(TF)2000.0, (TF)1e5, (TF)1e5, (TF)0.0, (TF)0.0, (TF)0.0, (TF)0.0};
    diff     = {(TF)0.13, (TF)0.16, (TF)0.13, (TF)0.11, (TF)0.15, (TF)0.13, (TF)0.16};
    diff_scl = {(TF)1.6, (TF)1.3, (TF)1.6, (TF)1.9, (TF)1.4, (TF)1.6, (TF)1.3};
    henry    = {(TF)0.01, (TF)2e-3, (TF)0.01, (TF)1e14, (TF)1e5, (TF)240., (TF)6e3};
    f0       = {(TF)1.0, (TF)0.0, (TF)0.1, (TF)0.0, (TF)1.0, (TF)0.1, (TF)0.0};

    // Define uninitialized resistance values by scaling with O3 and SO2 resistances (Wesely 1989)
    for (int i=3; i<7; i++)
    {
        rmes[i]  = (TF)1.0 / (henry[i] / (TF)3000.0 + (TF)100.0 * f0[i]);
        rsoil[i] = (TF)1.0 / (henry[i] / (henry_so2 + rsoil_so2) + f0[i] / rsoil[0]);
        rcut[i]  = (TF)1.0 / (henry[i] / henry_so2  + f0[i]) * rcut[0];
        rws[i]   = (TF)1.0 / (TF(1.0) / ((TF)3.0*rws_so2) + (TF)1e-7 * henry[i] + f0[i] / rws[0]);
        rwat[i]  = (TF)1.0 / (henry[i] / (henry_so2 + rwat_so2)  + f0[i] / rwat[0]);
    }

    // Change diff_scl to diff_scl^(2/3) for use in rb calculation
    for (int i=0; i<7; i++) diff_scl[i] = pow(diff_scl[i], (TF)2.0/(TF)3.0);

    for (auto& tile : deposition_tiles)
    {
        std::fill(tile.second.vdo3.begin(),tile.second.vdo3.end(), vd_o3);
        std::fill(tile.second.vdno.begin(),tile.second.vdno.end(), vd_no);
        std::fill(tile.second.vdno2.begin(),tile.second.vdno2.end(), vd_no2);
        std::fill(tile.second.vdhno3.begin(),tile.second.vdhno3.end(), vd_hno3);
        std::fill(tile.second.vdh2o2.begin(),tile.second.vdh2o2.end(), vd_h2o2);
        std::fill(tile.second.vdrooh.begin(),tile.second.vdrooh.end(), vd_rooh);
        std::fill(tile.second.vdhcho.begin(),tile.second.vdhcho.end(), vd_hcho);
    }
}


template <typename TF>
void Deposition<TF>::create(Stats<TF>& stats, Cross<TF>& cross)
{
    if (!sw_deposition)
        return;

    // add cross-sections
    if (cross.get_switch())
    {
        std::vector<std::string> allowed_crossvars = {
                "vdo3_soil", "vdno_soil", "vdno2_soil", "vdhno3_soil",
                "vdh2o2_soil", "vdrooh_soil", "vdhcho_soil",
                "vdo3_wet", "vdno_wet", "vdno2_wet", "vdhno3_wet",
                "vdh2o2_wet", "vdrooh_wet", "vdhcho_wet",
                "vdo3_veg", "vdno_veg", "vdno2_veg", "vdhno3_veg",
                "vdh2o2_veg", "vdrooh_veg", "vdhcho_veg"};
        cross_list = cross.get_enabled_variables(allowed_crossvars);
    }
}


template <typename TF>
void Deposition<TF>::update_time_dependent(
        Timeloop<TF>& timeloop,
        Boundary<TF>& boundary,
        TF* restrict vdo3,
        TF* restrict vdno,
        TF* restrict vdno2,
        TF* restrict vdhno3,
        TF* restrict vdh2o2,
        TF* restrict vdrooh,
        TF* restrict vdhcho
        )
{
    if (!sw_deposition)
        return;

    auto& gd = grid.get_grid_data();

    // get information from lsm:
    auto& tiles = boundary.get_tiles();
    auto& lai = boundary.get_lai();
    auto& water_mask = boundary.get_water_mask();
    auto& c_veg = boundary.get_c_veg();

    // calculate deposition per tile:
    for (auto& tile : tiles)
    {
        calc_deposition_per_tile(
                tile.first,
                deposition_tiles.at(tile.first).vdo3.data(),
                deposition_tiles.at(tile.first).vdno.data(),
                deposition_tiles.at(tile.first).vdno2.data(),
                deposition_tiles.at(tile.first).vdhno3.data(),
                deposition_tiles.at(tile.first).vdh2o2.data(),
                deposition_tiles.at(tile.first).vdrooh.data(),
                deposition_tiles.at(tile.first).vdhcho.data(),
                lai.data(),
                c_veg.data(),
                tile.second.rs.data(),
                tiles.at("veg").rs.data(),
                tile.second.ra.data(),
                tile.second.ustar.data(),
                tile.second.fraction.data(),
                rmes.data(), rsoil.data(), rcut.data(),
                rws.data(), rwat.data(), diff_scl.data(),   // pass as constant....
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells);
    }

    // Calculate tile-mean deposition for chemistry
    get_tiled_mean(vdo3,"o3",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    get_tiled_mean(vdno,"no",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    get_tiled_mean(vdno2,"no2",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    get_tiled_mean(vdhno3,"hno3",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    get_tiled_mean(vdh2o2,"h2o2",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    get_tiled_mean(vdrooh,"rooh",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());
    get_tiled_mean(vdhcho,"hcho",(TF) 1.0,tiles.at("veg").fraction.data(), tiles.at("soil").fraction.data(), tiles.at("wet").fraction.data());

    // cmk: we use the wet-tile info for u* and ra, since these are calculated in lsm with f_wet = 100%
    update_vd_water(vdo3,"o3",tiles.at("wet").ra.data(),tiles.at("wet").ustar.data(),water_mask.data(),diff_scl.data(),rwat.data());
    update_vd_water(vdno,"no",tiles.at("wet").ra.data(),tiles.at("wet").ustar.data(),water_mask.data(),diff_scl.data(),rwat.data());
    update_vd_water(vdno2,"no2",tiles.at("wet").ra.data(),tiles.at("wet").ustar.data(),water_mask.data(),diff_scl.data(),rwat.data());
    update_vd_water(vdhno3,"hno3",tiles.at("wet").ra.data(),tiles.at("wet").ustar.data(),water_mask.data(),diff_scl.data(),rwat.data());
    update_vd_water(vdh2o2,"h2o2",tiles.at("wet").ra.data(),tiles.at("wet").ustar.data(),water_mask.data(),diff_scl.data(),rwat.data());
    update_vd_water(vdrooh,"rooh",tiles.at("wet").ra.data(),tiles.at("wet").ustar.data(),water_mask.data(),diff_scl.data(),rwat.data());
    update_vd_water(vdhcho,"hcho",tiles.at("wet").ra.data(),tiles.at("wet").ustar.data(),water_mask.data(),diff_scl.data(),rwat.data());

    spatial_avg_vd(vdo3);
    spatial_avg_vd(vdno);
    spatial_avg_vd(vdno2);
    spatial_avg_vd(vdhno3);
    spatial_avg_vd(vdh2o2);
    spatial_avg_vd(vdrooh);
    spatial_avg_vd(vdhcho);

}


template<typename TF>
void Deposition<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    if (!sw_deposition)
        return;

    auto& gd = grid.get_grid_data();

    const TF no_offset = TF(0);

    for (auto& name : cross_list)
    {
        if (name == "vdo3_veg")
            cross.cross_plane(deposition_tiles.at("veg").vdo3.data(), no_offset, name, iotime);
        else if (name == "vdno_veg")
            cross.cross_plane(deposition_tiles.at("veg").vdno.data(), no_offset, name, iotime);
        else if (name == "vdno2_veg")
            cross.cross_plane(deposition_tiles.at("veg").vdno2.data(), no_offset, name, iotime);
        else if (name == "vdhno3_veg")
            cross.cross_plane(deposition_tiles.at("veg").vdhno3.data(), no_offset, name, iotime);
        else if (name == "vdh2o2_veg")
            cross.cross_plane(deposition_tiles.at("veg").vdh2o2.data(), no_offset, name, iotime);
        else if (name == "vdrooh_veg")
            cross.cross_plane(deposition_tiles.at("veg").vdrooh.data(), no_offset, name, iotime);
        else if (name == "vdhcho_veg")
            cross.cross_plane(deposition_tiles.at("veg").vdhcho.data(), no_offset, name, iotime);
        else if (name == "vdo3_soil")
            cross.cross_plane(deposition_tiles.at("soil").vdo3.data(), no_offset, name, iotime);
        else if (name == "vdno_soil")
            cross.cross_plane(deposition_tiles.at("soil").vdno.data(), no_offset, name, iotime);
        else if (name == "vdno2_soil")
            cross.cross_plane(deposition_tiles.at("soil").vdno2.data(), no_offset, name, iotime);
        else if (name == "vdhno3_soil")
            cross.cross_plane(deposition_tiles.at("soil").vdhno3.data(), no_offset, name, iotime);
        else if (name == "vdh2o2_soil")
            cross.cross_plane(deposition_tiles.at("soil").vdh2o2.data(), no_offset, name, iotime);
        else if (name == "vdrooh_soil")
            cross.cross_plane(deposition_tiles.at("soil").vdrooh.data(), no_offset, name, iotime);
        else if (name == "vdhcho_soil")
            cross.cross_plane(deposition_tiles.at("soil").vdhcho.data(), no_offset, name, iotime);
        else if (name == "vdo3_wet")
            cross.cross_plane(deposition_tiles.at("wet").vdo3.data(), no_offset, name, iotime);
        else if (name == "vdno_wet")
            cross.cross_plane(deposition_tiles.at("wet").vdno.data(), no_offset, name, iotime);
        else if (name == "vdno2_wet")
            cross.cross_plane(deposition_tiles.at("wet").vdno2.data(), no_offset, name, iotime);
        else if (name == "vdhno3_wet")
            cross.cross_plane(deposition_tiles.at("wet").vdhno3.data(), no_offset, name, iotime);
        else if (name == "vdh2o2_wet")
            cross.cross_plane(deposition_tiles.at("wet").vdh2o2.data(), no_offset, name, iotime);
        else if (name == "vdrooh_wet")
            cross.cross_plane(deposition_tiles.at("wet").vdrooh.data(), no_offset, name, iotime);
        else if (name == "vdhcho_wet")
            cross.cross_plane(deposition_tiles.at("wet").vdhcho.data(), no_offset, name, iotime);
    }
}


template<typename TF>
const TF Deposition<TF>::get_vd(const std::string& name) const
{
    if (name == "o3")
        return vd_o3;
    else if (name == "no")
        return vd_no;
    else if (name == "no2")
        return vd_no2;
    else if (name == "hno3")
        return vd_hno3;
    else if (name == "h2o2")
        return vd_h2o2;
    else if (name == "rooh")
        return vd_rooh;
    else if (name == "hcho")
        return vd_hcho;
    else
    {
        std::string error = "Deposition::get_vd() can't return \"" + name + "\"";
        throw std::runtime_error(error);
    }
}


template<typename TF>
void Deposition<TF>::get_tiled_mean(
    TF* restrict fld_out, std::string name, const TF fac,
    const TF* const restrict fveg,
    const TF* const restrict fsoil,
    const TF* const restrict fwet)
{
    auto& gd = grid.get_grid_data();

    TF* fld_veg;
    TF* fld_soil;
    TF* fld_wet;

    // Yikes..
    if (name == "o3")
    {
        fld_veg  = deposition_tiles.at("veg").vdo3.data();
        fld_soil = deposition_tiles.at("soil").vdo3.data();
        fld_wet  = deposition_tiles.at("wet").vdo3.data();
    }
    else if (name == "no")
    {
        fld_veg  = deposition_tiles.at("veg").vdno.data();
        fld_soil = deposition_tiles.at("soil").vdno.data();
        fld_wet  = deposition_tiles.at("wet").vdno.data();
    }
    else if (name == "no2")
    {
        fld_veg  = deposition_tiles.at("veg").vdno2.data();
        fld_soil = deposition_tiles.at("soil").vdno2.data();
        fld_wet  = deposition_tiles.at("wet").vdno2.data();
    }
    else if (name == "hno3")
    {
        fld_veg  = deposition_tiles.at("veg").vdhno3.data();
        fld_soil = deposition_tiles.at("soil").vdhno3.data();
        fld_wet  = deposition_tiles.at("wet").vdhno3.data();
    }
    else if (name == "h2o2")
    {
        fld_veg  = deposition_tiles.at("veg").vdh2o2.data();
        fld_soil = deposition_tiles.at("soil").vdh2o2.data();
        fld_wet  = deposition_tiles.at("wet").vdh2o2.data();
    }
    else if (name == "rooh")
    {
        fld_veg  = deposition_tiles.at("veg").vdrooh.data();
        fld_soil = deposition_tiles.at("soil").vdrooh.data();
        fld_wet  = deposition_tiles.at("wet").vdrooh.data();
    }
    else if (name == "hcho")
    {
        fld_veg  = deposition_tiles.at("veg").vdhcho.data();
        fld_soil = deposition_tiles.at("soil").vdhcho.data();
        fld_wet  = deposition_tiles.at("wet").vdhcho.data();
    }
    else
        throw std::runtime_error("Cannot calculate tiled mean for variable \"" + name + "\"\\n");

    calc_tiled_mean(
            fld_out,
            fveg,
            fsoil,
            fwet,
            fld_veg,
            fld_soil,
            fld_wet,
            fac,
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);
}


template<typename TF>
void Deposition<TF>::update_vd_water(
        TF* restrict fld_out, std::string name,
        const TF* const restrict ra,
        const TF* const restrict ustar,
        const int* const restrict water_mask,
        const TF* const restrict diff_scl,
        const TF* const restrict rwat)
{
    auto& gd = grid.get_grid_data();

    // TF* fld;
    TF diff_scl_val;
    TF rwat_val;

    // Yikes...
    if (name == "o3")
    {
        // fld = vd_o3.data();
        diff_scl_val = diff_scl[0];
        rwat_val = rwat[0];
    }
    else if (name == "no")
    {
        // fld = vd_no.data();
        diff_scl_val = diff_scl[1];
        rwat_val = rwat[1];
    }
    else if (name == "no2")
    {
        // fld = vd_no2.data();
        diff_scl_val = diff_scl[2];
        rwat_val = rwat[2];
    }
    else if (name == "hno3")
    {
        // fld = vd_hno3.data();
        diff_scl_val = diff_scl[3];
        rwat_val = rwat[3];
    }
    else if (name == "h2o2")
    {
        // fld = vd_h2o2.data();
        diff_scl_val = diff_scl[4];
        rwat_val = rwat[4];
    }
    else  if (name == "rooh")
    {
        // fld = vd_rooh.data();
        diff_scl_val = diff_scl[5];
        rwat_val = rwat[5];
    }
    else if (name == "hcho")
    {
        // fld = vd_hcho.data();
        diff_scl_val = diff_scl[6];
        rwat_val = rwat[6];
    }
    else
        throw std::runtime_error("Cannot update vd to water for variable \"" + name + "\"\\n");

    calc_vd_water(
        fld_out,
        ra,
        ustar,
        water_mask,
        diff_scl_val,
        rwat_val,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.icells);
}


template<typename TF>
void Deposition<TF>::spatial_avg_vd(
    TF* restrict fld_out)
{
    auto& gd = grid.get_grid_data();

    calc_spatial_avg_deposition(
        fld_out,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.icells);

}


template class Deposition<double>;
//:template class Chemistry<float>;
