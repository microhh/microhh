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

#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "input.h"
#include "netcdf_interface.h"
#include "constants.h"
#include "timeloop.h"
#include "trajectory.h"

namespace
{
    template<typename TF> TF netcdf_fp_fillvalue();
    template<> double netcdf_fp_fillvalue<double>() { return NC_FILL_DOUBLE; }
    template<> float  netcdf_fp_fillvalue<float>()  { return NC_FILL_FLOAT; }
}

template<typename TF>
Trajectory<TF>::Trajectory(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    names = inputin.get_list<std::string>("trajectory", "trajectories", "", std::vector<std::string>());
    sw_trajectory = names.size() > 0;

    if (sw_trajectory)
    {
        variables = inputin.get_list<std::string>("trajectory", "variables", "");
        sampletime = inputin.get_item<double>("trajectory", "sampletime", "");
    }
}

template<typename TF>
Trajectory<TF>::~Trajectory()
{
}

template<typename TF>
void Trajectory<TF>::init()
{
    if (!sw_trajectory)
        return;

    isampletime = convert_to_itime(sampletime);
    statistics_counter = 0;
}

template<typename TF>
void Trajectory<TF>::create(
        Input& inputin, Netcdf_file& input_nc,
        Timeloop<TF>& timeloop, std::string sim_name)
{
    // Do not create file if column is disabled.
    if (!sw_trajectory)
        return;

    for (auto& name : names)
    {
        // Get time and location from input NetCDF file.
        Netcdf_group& group_nc = input_nc.get_group("trajectory_" + name);
        const int n_traj = group_nc.get_dimension_size("itraj");

        // Create new trajectory instance, and add to map.
        trajectories.emplace(name, Single_trajectory<TF>{});

        // Short-cut to trajectory.
        Single_trajectory<TF>& traj = trajectories.at(name);

        // Retrieve input coordinates (time/space) from `case_input.nc`.
        traj.time_in.resize(n_traj);
        traj.x_in.resize(n_traj);
        traj.y_in.resize(n_traj);
        traj.z_in.resize(n_traj);

        group_nc.get_variable(traj.time_in, "time", {0}, {n_traj});
        group_nc.get_variable(traj.x_in, "x", {0}, {n_traj});
        group_nc.get_variable(traj.y_in, "y", {0}, {n_traj});
        group_nc.get_variable(traj.z_in, "z", {0}, {n_traj});

        // Create a NetCDF file for the statistics.
        std::stringstream filename;
        filename << sim_name << "." << "trajectory" << "." << name << "."
                 << std::setfill('0') << std::setw(7) << timeloop.get_iotime() << ".nc";

        // Create new NetCDF file.
        traj.data_file = std::make_unique<Netcdf_file>(master, filename.str(), Netcdf_mode::Create);
        traj.data_file->add_dimension("time");

        // Create dimensions and variables.
        traj.time_var = std::make_unique<Netcdf_variable<TF>>(
                    traj.data_file->template add_variable<TF>("time", {"time"}));
        if (timeloop.has_utc_time())
            traj.time_var->add_attribute("units", "seconds since " + timeloop.get_datetime_utc_start_string());
        else
            traj.time_var->add_attribute("units", "seconds since start");
        traj.time_var->add_attribute("long_name", "Time");

        traj.x_var = std::make_unique<Netcdf_variable<TF>>(
                    traj.data_file->template add_variable<TF>("x", {"time"}));
        traj.x_var->add_attribute("units", "m");
        traj.x_var->add_attribute("long_name", "X location trajectory");

        traj.y_var = std::make_unique<Netcdf_variable<TF>>(
                    traj.data_file->template add_variable<TF>("y", {"time"}));
        traj.y_var->add_attribute("units", "m");
        traj.y_var->add_attribute("long_name", "Y location trajectory");

        traj.z_var = std::make_unique<Netcdf_variable<TF>>(
                    traj.data_file->template add_variable<TF>("z", {"time"}));
        traj.z_var->add_attribute("units", "m");
        traj.z_var->add_attribute("long_name", "Z location trajectory");

        for (auto& name : variables)
        {
            // Check if requested field is prognostic or diagnostic. If not, raise exception.
            if (!fields.a.count(name))
            {
                std::string error = "Trajectory variable \"" + name + "\" is not a prognostic/diagnostic field.";
                throw std::runtime_error(error);
            }

            Time_var<TF> var{traj.data_file->template add_variable<TF>(name, {"time"}), TF(0)};
            var.ncvar.add_attribute("units", fields.a.at(name)->unit);
            var.ncvar.add_attribute("long_name", fields.a.at(name)->longname);
            //var.ncvar.add_attribute("_FillValue", netcdf_fp_fillvalue<TF>());  // ?!??!!?

            traj.time_series.emplace(
                    std::piecewise_construct, std::forward_as_tuple(name), std::forward_as_tuple(std::move(var)));
        }

        traj.data_file->sync();
    }
}

template<typename TF>
unsigned long Trajectory<TF>::get_time_limit(unsigned long itime)
{
    if (!sw_trajectory)
        return Constants::ulhuge;

    return isampletime - itime % isampletime;
}

template<typename TF>
bool Trajectory<TF>::do_trajectory(unsigned long itime)
{
    // Check if column are enabled
    if (!sw_trajectory)
        return false;

    // Check if time for execution
    if (itime % isampletime != 0)
        return false;

    // Return true such that column are computed
    return true;
}

template<typename TF>
void Trajectory<TF>::exec(Timeloop<TF>& timeloop, double time, unsigned long itime)
{
    // Write message in trajectory stats is triggered.
    master.print_message("Saving trajectory for time %f\n", time);

    auto& gd=grid.get_grid_data();
    auto& md=master.get_MPI_data();

    for (auto& name : names)
    {
        // Short-cut to trajectory.
        Single_trajectory<TF>& traj = trajectories.at(name);

        bool do_trajectory = true;

        // 1. Check on time bounds. Trajectory is allowed to start/end during the simulation.
        if (time < traj.time_in[0] || time > traj.time_in[traj.time_in.size()-1])
            do_trajectory = false;

        if (do_trajectory)
        {
            // 2. Get (x,y,z) location by interpolation of the input coordinates in time.
            Interpolation_factors<TF> ifac = timeloop.get_interpolation_factors(traj.time_in);

            traj.x_loc = ifac.fac0 * traj.x_in[ifac.index0] + ifac.fac1 * traj.x_in[ifac.index1];
            traj.y_loc = ifac.fac0 * traj.y_in[ifac.index0] + ifac.fac1 * traj.y_in[ifac.index1];
            traj.z_loc = ifac.fac0 * traj.z_in[ifac.index0] + ifac.fac1 * traj.z_in[ifac.index1];

            // Bounds check on domain. NOTE: keep the `>=xsize` instead of `>xsize`...
            if (traj.x_loc < 0 || traj.x_loc >= gd.xsize ||
                traj.y_loc < 0 || traj.y_loc >= gd.ysize ||
                traj.z_loc < 0 || traj.z_loc >= gd.zsize)
                do_trajectory = false;
        }

        if (do_trajectory)
        {
            // 2. Get indexes and interpolation factors for full and half levels.
            auto get_index_fac=[&](const std::vector<TF>& x_vec, const TF x_val, const int size)
            {
                // 1. Get index in `x_vec` left of `x_val`.
                int n0;
                for (int n=0; n<size; ++n)
                {
                    if (x_vec[n]<=x_val && x_vec[n + 1]>x_val)
                    {
                        n0=n;
                        break;
                    }
                }

                // 2. Calculate interpolation factor.
                const TF fac=TF(1) - (x_val - x_vec[n0]) / (x_vec[n0 + 1] - x_vec[n0]);

                return std::pair(n0, fac);
            };

            // Calculate `mpiid` which contains the current trajectory location.
            // All MPI tasks need to know this ID for the `broadcast` below.
            const TF sub_xsize=gd.imax * gd.dx;
            const TF sub_ysize=gd.jmax * gd.dy;

            const int mpicoordx=traj.x_loc / sub_xsize;
            const int mpicoordy=traj.y_loc / sub_ysize;

            const int mpiid_of_traj=master.calc_mpiid(mpicoordx, mpicoordy);

            if (md.mpiid == mpiid_of_traj)
            {
                std::pair<int, TF> ipx=get_index_fac(gd.x, traj.x_loc, gd.icells);
                std::pair<int, TF> ipxh=get_index_fac(gd.xh, traj.x_loc, gd.icells);

                std::pair<int, TF> jpx=get_index_fac(gd.y, traj.y_loc, gd.jcells);
                std::pair<int, TF> jpxh=get_index_fac(gd.yh, traj.y_loc, gd.jcells);

                std::pair<int, TF> kpx=get_index_fac(gd.z, traj.z_loc, gd.kcells);
                std::pair<int, TF> kpxh=get_index_fac(gd.zh, traj.z_loc, gd.kcells);

                // Put in vector for easier access with `field3d->loc`.
                std::vector<std::pair<int, TF>> ix={ipx, ipxh};
                std::vector<std::pair<int, TF>> iy={jpx, jpxh};
                std::vector<std::pair<int, TF>> iz={kpx, kpxh};

                // Interpolate variables.
                auto ijk=[&](const int i, const int j, const int k)
                {
                    return i + j * gd.icells + k * gd.ijcells;
                };

                auto interpolate=[&](
                        const std::vector<TF>& fld,
                        std::pair<int, TF> fx,
                        std::pair<int, TF> fy,
                        std::pair<int, TF> fz)
                {
                    // Short-cuts.
                    const int i0=fx.first;
                    const int j0=fy.first;
                    const int k0=fz.first;

                    const TF fx0=fx.second;
                    const TF fx1=TF(1) - fx.second;
                    const TF fy0=fy.second;
                    const TF fy1=TF(1) - fy.second;
                    const TF fz0=fz.second;
                    const TF fz1=TF(1) - fz.second;

                    // Tri-linear interpolation onto requested location.
                    const TF value=
                            fx0 * fy0 * fz0 * fld[ijk(i0, j0, k0)] +
                            fx1 * fy0 * fz0 * fld[ijk(i0 + 1, j0, k0)] +
                            fx0 * fy1 * fz0 * fld[ijk(i0, j0 + 1, k0)] +
                            fx1 * fy1 * fz0 * fld[ijk(i0 + 1, j0 + 1, k0)] +
                            fx0 * fy0 * fz1 * fld[ijk(i0, j0, k0 + 1)] +
                            fx1 * fy0 * fz1 * fld[ijk(i0 + 1, j0, k0 + 1)] +
                            fx0 * fy1 * fz1 * fld[ijk(i0, j0 + 1, k0 + 1)] +
                            fx1 * fy1 * fz1 * fld[ijk(i0 + 1, j0 + 1, k0 + 1)];

                    return value;
                };

                for (auto& name: variables)
                {
                    const int loc_i=fields.a.at(name)->loc[0];
                    const int loc_j=fields.a.at(name)->loc[1];
                    const int loc_k=fields.a.at(name)->loc[2];

                    traj.time_series.at(name).data=interpolate(
                            fields.a.at(name)->fld, ix[loc_i], iy[loc_j], iz[loc_k]);
                }
            }

            // This is a bit wasteful (?), a specific send/recv should be enough,
            // but we are only sending one float/double per variable....
            for (auto& name : variables)
                master.broadcast(&traj.time_series.at(name).data, 1, mpiid_of_traj);
        }
        else
        {
            // Trajectory is out of bounds in time or space (domain).
            for (auto& name: variables)
                traj.time_series.at(name).data = netcdf_fp_fillvalue<TF>();

            traj.x_loc = netcdf_fp_fillvalue<TF>();
            traj.y_loc = netcdf_fp_fillvalue<TF>();
            traj.z_loc = netcdf_fp_fillvalue<TF>();
        }

        // Store value in NetCDF file.
        const std::vector<int> time_index{statistics_counter};

        traj.time_var->insert(time, time_index);
        traj.x_var->insert(traj.x_loc, time_index);
        traj.y_var->insert(traj.y_loc, time_index);
        traj.z_var->insert(traj.z_loc, time_index);

        for (auto& name : variables)
            traj.time_series.at(name).ncvar.insert(traj.time_series.at(name).data, time_index);

        // Synchronize the NetCDF file
        traj.data_file->sync();
    }

    ++statistics_counter;
}

#ifdef FLOAT_SINGLE
template class Trajectory<float>;
#else
template class Trajectory<double>;
#endif
