
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
#include "boundary.h"
#include "cross.h"

#include "chemistry.h"
#include "chemistry_plume_kernels.h"

namespace cpk = Chemistry_plume_kernels;


template<typename TF>
Chemistry<TF>::Chemistry(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), field3d_operators(master, grid, fields)
{
    const std::string group_name = "default";
    auto& gd = grid.get_grid_data();

    sw_chemistry = inputin.get_item<bool>("chemistry", "swchemistry", "", false);
    //sw_scheme = inputin.get_item<bool>("chemistry", "swscheme", "", false);

    if (!sw_chemistry)
        return;

    deposition = std::make_shared<Deposition <TF>>(masterin, gridin, fieldsin, inputin);
}

template <typename TF>
Chemistry<TF>::~Chemistry()
{
}

#ifndef USECUDA
template <typename TF>
void Chemistry<TF>::update_time_dependent(Timeloop<TF>& timeloop, Boundary<TF>& boundary)
{
    if (!sw_chemistry)
        return;

    Interpolation_factors<TF> ifac = timeloop.get_interpolation_factors(time);

    jval[Jval::o31d]   = ifac.fac0 * jo31d[ifac.index0]   + ifac.fac1 * jo31d[ifac.index1];
    jval[Jval::h2o2]   = ifac.fac0 * jh2o2[ifac.index0]   + ifac.fac1 * jh2o2[ifac.index1];
    jval[Jval::no2]    = ifac.fac0 * jno2[ifac.index0]    + ifac.fac1 * jno2[ifac.index1];
    jval[Jval::no3]    = ifac.fac0 * jno3[ifac.index0]    + ifac.fac1 * jno3[ifac.index1];
    jval[Jval::n2o5]   = ifac.fac0 * jn2o5[ifac.index0]   + ifac.fac1 * jn2o5[ifac.index1];
    jval[Jval::ch2or]  = ifac.fac0 * jch2or[ifac.index0]  + ifac.fac1 * jch2or[ifac.index1];
    jval[Jval::ch2om]  = ifac.fac0 * jch2om[ifac.index0]  + ifac.fac1 * jch2om[ifac.index1];
    jval[Jval::ch3o2h] = ifac.fac0 * jch3o2h[ifac.index0] + ifac.fac1 * jch3o2h[ifac.index1];

    deposition->update_time_dependent(
            timeloop,
            boundary,
            vdo3.data(),
            vdno.data(),
            vdno2.data(),
            vdhno3.data(),
            vdh2o2.data(),
            vdrooh.data(),
            vdhcho.data());
}

template <typename TF>
void Chemistry<TF>::exec(Thermo<TF>& thermo, const double sdt, const double dt)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();

    // Calculate the mean temperature profile.
    auto temperature = fields.get_tmp();
    thermo.get_thermo_field(*temperature, "T", true, false);
    field3d_operators.calc_mean_profile(temperature->fld_mean.data(), temperature->fld.data());

    cpk::pss<TF>(
        fields.st.at("hno3")->fld.data(),
        fields.st.at("h2o2")->fld.data(),
        fields.st.at("co")  ->fld.data(),
        fields.st.at("hcho")->fld.data(),
        fields.st.at("rooh")->fld.data(),
        fields.st.at("c3h6")->fld.data(),
        fields.st.at("o3")  ->fld.data(),
        fields.st.at("no")  ->fld.data(),
        fields.st.at("no2") ->fld.data(),
        fields.sp.at("hno3")->fld.data(),
        fields.sp.at("h2o2")->fld.data(),
        fields.sp.at("co")  ->fld.data(),
        fields.sp.at("hcho")->fld.data(),
        fields.sp.at("rooh")->fld.data(),
        fields.sp.at("c3h6")->fld.data(),
        fields.sp.at("o3")  ->fld.data(),
        fields.sp.at("no")  ->fld.data(),
        fields.sp.at("no2") ->fld.data(),
        jval.data(),
        vdo3.data(),
        vdno.data(),
        vdno2.data(),
        vdhno3.data(),
        vdh2o2.data(),
        vdrooh.data(),
        vdhcho.data(),
        temperature->fld_mean.data(),
        fields.sp.at("qt")->fld_mean.data(),
        gd.dzi.data(),
        fields.rhoref.data(),
        sdt,
        gd.istart,
        gd.iend,
        gd.jstart,
        gd.jend,
        gd.kstart,
        gd.kend,
        gd.icells,
        gd.ijcells);

    fields.release_tmp(temperature);
}
#endif

template <typename TF>
void Chemistry<TF>::init(Input& inputin)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();

    // initialize 2D deposition arrays:
    vdo3.resize(gd.ijcells);
    vdno.resize(gd.ijcells);
    vdno2.resize(gd.ijcells);
    vdhno3.resize(gd.ijcells);
    vdh2o2.resize(gd.ijcells);
    vdrooh.resize(gd.ijcells);
    vdhcho.resize(gd.ijcells);

    // initialize deposition routine:
    deposition-> init(inputin);

    // fill deposition with standard values:
    std::fill(vdo3.begin(), vdo3.end(), deposition-> get_vd("o3"));
    std::fill(vdno.begin(), vdno.end(), deposition-> get_vd("no"));
    std::fill(vdno2.begin(), vdno2.end(), deposition-> get_vd("no2"));
    std::fill(vdhno3.begin(), vdhno3.end(), deposition-> get_vd("hno3"));
    std::fill(vdh2o2.begin(), vdh2o2.end(), deposition-> get_vd("h2o2"));
    std::fill(vdrooh.begin(), vdrooh.end(), deposition-> get_vd("rooh"));
    std::fill(vdhcho.begin(), vdhcho.end(), deposition-> get_vd("hcho"));

    master.print_message("Deposition arrays initialized, e.g. with vdo3 = %13.5e m/s \n", deposition-> get_vd("o3"));
}

template <typename TF>
void Chemistry<TF>::create(
        const Timeloop<TF>& timeloop, std::string sim_name, Netcdf_handle& input_nc,
        Stats<TF>& stats, Cross<TF>& cross)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();
    int iotime = timeloop.get_iotime();

    Netcdf_group& group_nc = input_nc.get_group("timedep_chem");
    int time_dim_length;
    std::string time_dim;

    for (std::string varname : jname)    // check dimensions:
    {
        std::map<std::string, int> dims = group_nc.get_variable_dimensions(varname);
        std::pair<std::string, int> unique_time = cpk::check_for_unique_time_dim(dims);
        time_dim = unique_time.first;
        time_dim_length = unique_time.second;
        time.resize(time_dim_length);
    }

    jval.resize(n_jval);
    jo31d.resize(time_dim_length);
    jh2o2.resize(time_dim_length);
    jno2.resize(time_dim_length);
    jno3.resize(time_dim_length);
    jn2o5.resize(time_dim_length);
    jch2or.resize(time_dim_length);
    jch2om.resize(time_dim_length);
    jch3o2h.resize(time_dim_length);

    group_nc.get_variable(time, time_dim, {0}, {time_dim_length});
    group_nc.get_variable(jo31d, jname[0],  {0}, {time_dim_length});
    group_nc.get_variable(jh2o2, jname[1],  {0}, {time_dim_length});
    group_nc.get_variable(jno2, jname[2],  {0}, {time_dim_length});
    group_nc.get_variable(jno3, jname[3],  {0}, {time_dim_length});
    group_nc.get_variable(jn2o5, jname[4],  {0}, {time_dim_length});
    group_nc.get_variable(jch2or, jname[5],  {0}, {time_dim_length});
    group_nc.get_variable(jch2om, jname[6],  {0}, {time_dim_length});
    group_nc.get_variable(jch3o2h, jname[7],  {0}, {time_dim_length});

    if (stats.get_switch())
    {
        // Stats:
        const std::string group_name = "default";
        const std::vector<std::string> stat_op_def = {"mean", "2", "3", "4", "w", "grad", "diff", "flux", "path"};
        const std::vector<std::string> stat_op_w = {"mean", "2", "3", "4"};
        const std::vector<std::string> stat_op_p = {"mean", "2", "w", "grad"};

        // add the deposition-velocity timeseries in deposition group statistics
        const std::string group_named = "deposition";

        // used in chemistry:
        stats.add_time_series("vdo3", "O3 deposition velocity", "m s-1", group_named);
        stats.add_time_series("vdno", "NO deposition velocity", "m s-1", group_named);
        stats.add_time_series("vdno2", "NO2 deposition velocity", "m s-1", group_named);
        stats.add_time_series("vdhno3", "HNO3 deposition velocity", "m s-1", group_named);
        stats.add_time_series("vdh2o2", "H2O2 deposition velocity", "m s-1", group_named);
        stats.add_time_series("vdrooh", "ROOH deposition velocity", "m s-1", group_named);
        stats.add_time_series("vdhcho", "HCHO deposition velocity", "m s-1", group_named);
    }

    // add cross-sections
    if (cross.get_switch())
    {
        std::vector<std::string> allowed_crossvars = {"vdo3", "vdno", "vdno2", "vdhno3", "vdh2o2", "vdrooh", "vdhcho"};
        cross_list = cross.get_enabled_variables(allowed_crossvars);

        // `deposition->create()` only creates cross-sections.
        deposition->create(stats, cross);
    }
}

template<typename TF>
void Chemistry<TF>::exec_stats(const int iteration, const double time, Stats<TF>& stats)
{
    if (!sw_chemistry or stats.get_switch())
        return;

    auto& gd = grid.get_grid_data();

    const TF no_offset = 0.;
    const TF no_threshold = 0.;

    if (iteration != 0)   // this does not make sense for first step = t=0.
    {
        // add deposition velocities to statistics:
        stats.calc_stats_2d("vdo3"   , vdo3,   no_offset);
        stats.calc_stats_2d("vdno"   , vdno,   no_offset);
        stats.calc_stats_2d("vdno2"  , vdno2,  no_offset);
        stats.calc_stats_2d("vdhno3" , vdhno3, no_offset);
        stats.calc_stats_2d("vdh2o2" , vdh2o2, no_offset);
        stats.calc_stats_2d("vdrooh" , vdrooh, no_offset);
        stats.calc_stats_2d("vdhcho" , vdhcho, no_offset);
    }
}


template<typename TF>
void Chemistry<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    if (!sw_chemistry)
        return;

    auto& gd = grid.get_grid_data();

    const TF no_offset = TF(0);

    for (auto& name : cross_list)
    {
        if (name == "vdo3")
            cross.cross_plane(vdo3.data(), no_offset, name, iotime);
        else if (name == "vdno")
            cross.cross_plane(vdno.data(), no_offset, name, iotime);
        else if (name == "vdno2")
            cross.cross_plane(vdno2.data(), no_offset, name, iotime);
        else if (name == "vdhno3")
            cross.cross_plane(vdhno3.data(), no_offset, name, iotime);
        else if (name == "vdh2o2")
            cross.cross_plane(vdh2o2.data(), no_offset, name, iotime);
        else if (name == "vdrooh")
            cross.cross_plane(vdrooh.data(), no_offset, name, iotime);
        else if (name == "vdhcho")
            cross.cross_plane(vdhcho.data(), no_offset, name, iotime);
    }

    // See if to write per tile:
    deposition->exec_cross(cross, iotime);
}



#ifdef FLOAT_SINGLE
template class Chemistry<float>;
#else
template class Chemistry<double>;
#endif