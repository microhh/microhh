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

#include <cstdio>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
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
#include "boundary_surface_lsm.h"
#include "boundary.h"
#include "cross.h"

#include "deposition.h"
#include "deposition_kernels.h"

namespace dk = Deposition_kernels;


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
    vd_o3   = inputin.get_item<TF>("deposition", "vdo3", "", TF(0.005));
    vd_no   = inputin.get_item<TF>("deposition", "vdno", "", TF(0.002));
    vd_no2  = inputin.get_item<TF>("deposition", "vdno2", "", TF(0.005));
    vd_hno3 = inputin.get_item<TF>("deposition", "vdhno3", "", TF(0.040));
    vd_h2o2 = inputin.get_item<TF>("deposition", "vdh2o2", "", TF(0.018));
    vd_rooh = inputin.get_item<TF>("deposition", "vdrooh", "", TF(0.008));
    vd_hcho = inputin.get_item<TF>("deposition", "vdhcho", "", TF(0.0033));

    if (!sw_deposition)
        return;

    auto& gd = grid.get_grid_data();

    // Create surface tiles for deposition:
    for (auto& name : deposition_tile_names)
        deposition_tiles.emplace(name, Deposition_tile<TF>{});

    for (auto& tile : deposition_tiles)
        for (auto& [sp, idx] : species_idx)
            tile.second.vd[sp].resize(gd.ijcells);

    deposition_tiles.at("veg" ).long_name = "vegetation";
    deposition_tiles.at("soil").long_name = "bare soil";
    deposition_tiles.at("wet" ).long_name = "wet skin";
    deposition_var = inputin.get_item<TF>("deposition", "deposition_var","", TF(1e5));

    henry_so2 = inputin.get_item<TF>("deposition", "henry_so2", "", TF(1e5));
    rsoil_so2 = inputin.get_item<TF>("deposition", "rsoil_so2", "", TF(250));
    rwat_so2  = inputin.get_item<TF>("deposition", "rwat_so2", "", TF(1));
    rws_so2   = inputin.get_item<TF>("deposition", "rws_so2", "", TF(100));

    // Note: rmes for NO and NO2 (indices 1 and 2) will still be scaled with rs
    rmes     = {TF(1),    TF(5),   TF(0.5), TF(0),   TF(0),   TF(0),   TF(0)};
    rsoil    = {TF(400),  TF(1e5), TF(600), TF(0),   TF(0),   TF(0),   TF(0)};
    rcut     = {TF(1e5),  TF(1e5), TF(1e5), TF(0),   TF(0),   TF(0),   TF(0)};
    rws      = {TF(2000), TF(1e5), TF(1e5), TF(0),   TF(0),   TF(0),   TF(0)};
    rwat     = {TF(2000), TF(1e5), TF(1e5), TF(0),   TF(0),   TF(0),   TF(0)};
    diff     = {TF(0.13), TF(0.16),TF(0.13),TF(0.11),TF(0.15),TF(0.13),TF(0.16)};
    diff_scl = {TF(1.6),  TF(1.3), TF(1.6), TF(1.9), TF(1.4), TF(1.6), TF(1.3)};
    henry    = {TF(0.01), TF(2e-3),TF(0.01),TF(1e14),TF(1e5), TF(240), TF(6e3)};
    f0       = {TF(1),    TF(0),   TF(0.1), TF(0),   TF(1),   TF(0.1), TF(0)};

    // Define uninitialized resistance values by scaling with O3 and SO2 resistances (Wesely 1989)
    for (int i=3; i<7; i++)
    {
        rmes[i]  = TF(1) / (henry[i] / TF(3000) + TF(100) * f0[i]);
        rsoil[i] = TF(1) / (henry[i] / (henry_so2 + rsoil_so2) + f0[i] / rsoil[0]);
        rcut[i]  = TF(1) / (henry[i] / henry_so2 + f0[i]) * rcut[0];
        rws[i]   = TF(1) / (TF(1) / (TF(3) * rws_so2) + TF(1e-7) * henry[i] + f0[i] / rws[0]);
        rwat[i]  = TF(1) / (henry[i] / (henry_so2 + rwat_so2) + f0[i] / rwat[0]);
    }

    // Change diff_scl to diff_scl^(2/3) for use in rb calculation
    for (int i=0; i<7; i++)
        diff_scl[i] = pow(diff_scl[i], TF(2)/TF(3));

    const std::map<std::string, TF> vd_defaults = {
        {"o3", vd_o3}, {"no", vd_no}, {"no2", vd_no2}, {"hno3", vd_hno3},
        {"h2o2", vd_h2o2}, {"rooh", vd_rooh}, {"hcho", vd_hcho}};

    for (auto& tile : deposition_tiles)
        for (auto& [sp, vd_default] : vd_defaults)
            std::fill(tile.second.vd[sp].begin(), tile.second.vd[sp].end(), vd_default);
}


template <typename TF>
void Deposition<TF>::create(Stats<TF>& stats, Cross<TF>& cross)
{
    if (!sw_deposition)
        return;

    // Add cross-sections.
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


#ifndef USECUDA
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

    // Get information from the land-surface model:
    auto& tiles = boundary.get_tiles();
    auto& water_mask = boundary.get_water_mask();
    auto& lai = boundary.get_lai();
    auto& c_veg = boundary.get_c_veg();

    // Calculate deposition per tile.
    auto& dep_veg  = deposition_tiles.at("veg");
    auto& dep_soil = deposition_tiles.at("soil");
    auto& dep_wet  = deposition_tiles.at("wet");

    dk::calc_deposition_veg(
            dep_veg.vd.at("o3").data(),
            dep_veg.vd.at("no").data(),
            dep_veg.vd.at("no2").data(),
            dep_veg.vd.at("hno3").data(),
            dep_veg.vd.at("h2o2").data(),
            dep_veg.vd.at("rooh").data(),
            dep_veg.vd.at("hcho").data(),
            lai.data(),
            tiles.at("veg").rs.data(),
            tiles.at("veg").ra.data(),
            tiles.at("veg").ustar.data(),
            tiles.at("veg").fraction.data(),
            rmes.data(),
            rsoil.data(),
            rcut.data(),
            diff_scl.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

    dk::calc_deposition_soil(
            dep_soil.vd.at("o3").data(),
            dep_soil.vd.at("no").data(),
            dep_soil.vd.at("no2").data(),
            dep_soil.vd.at("hno3").data(),
            dep_soil.vd.at("h2o2").data(),
            dep_soil.vd.at("rooh").data(),
            dep_soil.vd.at("hcho").data(),
            tiles.at("soil").ra.data(),
            tiles.at("soil").ustar.data(),
            tiles.at("soil").fraction.data(),
            rsoil.data(),
            diff_scl.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

    dk::calc_deposition_wet(
            dep_wet.vd.at("o3").data(),
            dep_wet.vd.at("no").data(),
            dep_wet.vd.at("no2").data(),
            dep_wet.vd.at("hno3").data(),
            dep_wet.vd.at("h2o2").data(),
            dep_wet.vd.at("rooh").data(),
            dep_wet.vd.at("hcho").data(),
            lai.data(),
            c_veg.data(),
            tiles.at("wet").rs.data(),
            tiles.at("veg").rs.data(),
            tiles.at("wet").ra.data(),
            tiles.at("wet").ustar.data(),
            tiles.at("wet").fraction.data(),
            rmes.data(),
            rsoil.data(),
            rws.data(),
            diff_scl.data(),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.icells);

    // Calculate tile-mean deposition and apply water correction per species.
    auto calc_vd = [&](TF* vd, const std::string& name)
    {
        dk::calc_tiled_mean(
                vd,
                tiles.at("veg").fraction.data(),
                tiles.at("soil").fraction.data(),
                tiles.at("wet").fraction.data(),
                deposition_tiles.at("veg").vd.at(name).data(),
                deposition_tiles.at("soil").vd.at(name).data(),
                deposition_tiles.at("wet").vd.at(name).data(),
                TF(1),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells);

        // Use wet-tile u* and ra: calculated in lsm with f_wet = 100%.
        const int s = species_idx.at(name);
        dk::calc_vd_water(
                vd,
                tiles.at("wet").ra.data(),
                tiles.at("wet").ustar.data(),
                water_mask.data(),
                diff_scl[s], rwat[s],
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.icells);
    };

    calc_vd(vdo3,  "o3");
    calc_vd(vdno,  "no");
    calc_vd(vdno2, "no2");
    calc_vd(vdhno3,"hno3");
    calc_vd(vdh2o2,"h2o2");
    calc_vd(vdrooh,"rooh");
    calc_vd(vdhcho,"hcho");

    spatial_avg_vd(vdo3);
    spatial_avg_vd(vdno);
    spatial_avg_vd(vdno2);
    spatial_avg_vd(vdhno3);
    spatial_avg_vd(vdh2o2);
    spatial_avg_vd(vdrooh);
    spatial_avg_vd(vdhcho);
}
#endif


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
            cross.cross_plane(deposition_tiles.at("veg").vd.at("o3").data(), no_offset, name, iotime);
        else if (name == "vdno_veg")
            cross.cross_plane(deposition_tiles.at("veg").vd.at("no").data(), no_offset, name, iotime);
        else if (name == "vdno2_veg")
            cross.cross_plane(deposition_tiles.at("veg").vd.at("no2").data(), no_offset, name, iotime);
        else if (name == "vdhno3_veg")
            cross.cross_plane(deposition_tiles.at("veg").vd.at("hno3").data(), no_offset, name, iotime);
        else if (name == "vdh2o2_veg")
            cross.cross_plane(deposition_tiles.at("veg").vd.at("h2o2").data(), no_offset, name, iotime);
        else if (name == "vdrooh_veg")
            cross.cross_plane(deposition_tiles.at("veg").vd.at("rooh").data(), no_offset, name, iotime);
        else if (name == "vdhcho_veg")
            cross.cross_plane(deposition_tiles.at("veg").vd.at("hcho").data(), no_offset, name, iotime);
        else if (name == "vdo3_soil")
            cross.cross_plane(deposition_tiles.at("soil").vd.at("o3").data(), no_offset, name, iotime);
        else if (name == "vdno_soil")
            cross.cross_plane(deposition_tiles.at("soil").vd.at("no").data(), no_offset, name, iotime);
        else if (name == "vdno2_soil")
            cross.cross_plane(deposition_tiles.at("soil").vd.at("no2").data(), no_offset, name, iotime);
        else if (name == "vdhno3_soil")
            cross.cross_plane(deposition_tiles.at("soil").vd.at("hno3").data(), no_offset, name, iotime);
        else if (name == "vdh2o2_soil")
            cross.cross_plane(deposition_tiles.at("soil").vd.at("h2o2").data(), no_offset, name, iotime);
        else if (name == "vdrooh_soil")
            cross.cross_plane(deposition_tiles.at("soil").vd.at("rooh").data(), no_offset, name, iotime);
        else if (name == "vdhcho_soil")
            cross.cross_plane(deposition_tiles.at("soil").vd.at("hcho").data(), no_offset, name, iotime);
        else if (name == "vdo3_wet")
            cross.cross_plane(deposition_tiles.at("wet").vd.at("o3").data(), no_offset, name, iotime);
        else if (name == "vdno_wet")
            cross.cross_plane(deposition_tiles.at("wet").vd.at("no").data(), no_offset, name, iotime);
        else if (name == "vdno2_wet")
            cross.cross_plane(deposition_tiles.at("wet").vd.at("no2").data(), no_offset, name, iotime);
        else if (name == "vdhno3_wet")
            cross.cross_plane(deposition_tiles.at("wet").vd.at("hno3").data(), no_offset, name, iotime);
        else if (name == "vdh2o2_wet")
            cross.cross_plane(deposition_tiles.at("wet").vd.at("h2o2").data(), no_offset, name, iotime);
        else if (name == "vdrooh_wet")
            cross.cross_plane(deposition_tiles.at("wet").vd.at("rooh").data(), no_offset, name, iotime);
        else if (name == "vdhcho_wet")
            cross.cross_plane(deposition_tiles.at("wet").vd.at("hcho").data(), no_offset, name, iotime);
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
void Deposition<TF>::spatial_avg_vd(
    TF* restrict fld_out)
{
    auto& gd = grid.get_grid_data();

    dk::calc_spatial_avg_deposition(
        fld_out,
        gd.istart, gd.iend,
        gd.jstart, gd.jend,
        gd.icells);

}


#ifdef FLOAT_SINGLE
template class Deposition<float>;
#else
template class Deposition<double>;
#endif