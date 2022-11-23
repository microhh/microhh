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
#include <cmath>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <tuple>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "column.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "timeloop.h"
#include "netcdf_interface.h"

template<typename TF>
Column<TF>::Column(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    swcolumn = inputin.get_item<bool>("column", "swcolumn", "", false);

    if (swcolumn)
        sampletime = inputin.get_item<double>("column", "sampletime", "");
}

template<typename TF>
Column<TF>::~Column()
{
}

template<typename TF>
void Column<TF>::init(double ifactor)
{
    if (!swcolumn)
        return;

    isampletime = static_cast<unsigned long>(ifactor * sampletime);
    statistics_counter = 0;
}

template<typename TF>
void Column<TF>::create(Input& inputin, Timeloop<TF>& timeloop, std::string sim_name)
{
    // Do not create file if column is disabled.
    if (!swcolumn)
        return;

    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    std::vector<TF> coordx = inputin.get_list<TF>("column", "coordinates", "x", std::vector<TF>());
    std::vector<TF> coordy = inputin.get_list<TF>("column", "coordinates", "y", std::vector<TF>());

    if (coordx.size() != coordy.size())
    {
        std::string msg = "Column error: X-coord array and Y-coord array do not match in size";
        throw std::runtime_error(msg);
    }

    for (size_t n=0; n<coordx.size(); ++n)
    {
        if ( (coordx[n] < 0) || (coordx[n] > gd.xsize) || (coordy[n] < 0) || (coordy[n] > gd.ysize) )
        {
            std::string error = "Column #" + std::to_string(n) + " is outside the domain!";
            throw std::runtime_error(error);
        }

        int i = static_cast<int>(std::floor(coordx[n]/gd.dx));
        int j = static_cast<int>(std::floor(coordy[n]/gd.dy));

        // If column is at `xsize` or `ysize`, put it at itot/jtot-1:
        i = std::min(i, gd.itot-1);
        j = std::min(j, gd.jtot-1);

        columns.emplace_back(Column_struct{});
        columns.back().coord = {i, j};
    }

    // Create a NetCDF file for the statistics.
    for (auto& col : columns)
    {
        std::stringstream filename;
        filename << sim_name << "." << "column" << "."
                 << std::setfill('0') << std::setw(5) << col.coord[0] << "."
                 << std::setfill('0') << std::setw(5) << col.coord[1] << "."
                 << std::setfill('0') << std::setw(7) << timeloop.get_iotime() << ".nc";

        // Create new NetCDF file.
        // 1. Find the mpiid of the column.
        int mpiid_column = 0;
        if ( (col.coord[0] / gd.imax == md.mpicoordx ) && (col.coord[1] / gd.jmax == md.mpicoordy ) )
            mpiid_column = master.get_mpiid();
        master.sum(&mpiid_column, 1);

        // 2. Make the NetCDF file.
        col.data_file = std::make_unique<Netcdf_file>(
                master, filename.str(), Netcdf_mode::Create, mpiid_column);
    }

    // Create dimensions.
    for (auto& col : columns)
    {
        col.data_file->add_dimension("z",  gd.kmax);
        col.data_file->add_dimension("zh", gd.kmax+1);
        col.data_file->add_dimension("time");

        Netcdf_variable<TF> z_var = col.data_file->template add_variable<TF>("z", {"z"});
        z_var.add_attribute("units", "m");
        z_var.add_attribute("long_name", "Full level height");

        Netcdf_variable<TF> zh_var = col.data_file->template add_variable<TF>("zh", {"zh"});
        zh_var.add_attribute("units", "m");
        zh_var.add_attribute("long_name", "Half level height");

        // Save the grid variables.
        std::vector<TF> z_nogc (gd.z. begin() + gd.kstart, gd.z. begin() + gd.kend  );
        std::vector<TF> zh_nogc(gd.zh.begin() + gd.kstart, gd.zh.begin() + gd.kend+1);
        z_var .insert( z_nogc, {0});
        zh_var.insert(zh_nogc, {0});

        // create variables belonging to dimensions
        col.iter_var = std::make_unique<Netcdf_variable<int>>(
                col.data_file->template add_variable<int>("iter", {"time"}));
        col.iter_var->add_attribute("units", "-");
        col.iter_var->add_attribute("long_name", "Iteration number");

        col.time_var = std::make_unique<Netcdf_variable<TF>>(
                col.data_file->template add_variable<TF>("time", {"time"}));
        if (timeloop.has_utc_time())
            col.time_var->add_attribute("units", "seconds since " + timeloop.get_datetime_utc_start_string());
        else
            col.time_var->add_attribute("units", "seconds since start");
        col.time_var->add_attribute("long_name", "Time");

        col.data_file->sync();
    }
}

template<typename TF>
unsigned long Column<TF>::get_time_limit(unsigned long itime)
{
    if (!swcolumn)
        return Constants::ulhuge;

    // If sampletime is negative, output column every timestep
    if (isampletime == 0)
        return Constants::ulhuge;

    unsigned long idtlim = isampletime - itime % isampletime;
    return idtlim;
}

template<typename TF>
bool Column<TF>::do_column(unsigned long itime)
{
    // check if column are enabled
    if (!swcolumn)
        return false;

    // If sampletime is negative, output column every timestep
    if (isampletime == 0)
        return true;

    // check if time for execution
    if (itime % isampletime != 0)
        return false;

    // return true such that column are computed
    return true;
}

template<typename TF>
void Column<TF>::exec(int iteration, double time, unsigned long itime)
{
    // Write message in case stats is triggered.
    if (isampletime > 0)
        master.print_message("Saving columns for time %f\n", time);

    auto& gd = grid.get_grid_data();

    // Put the data into the NetCDF file.
    for (auto& col : columns)
    {
        const std::vector<int> time_index{statistics_counter};

        // Write the time and iteration number.
        col.time_var->insert(time     , time_index);
        col.iter_var->insert(iteration, time_index);

        const std::vector<int> time_height_index = {statistics_counter, 0};

        for (auto& p : col.profs)
        {
            const int ksize = p.second.ncvar.get_dim_sizes()[1];
            std::vector<int> time_height_size  = {1, ksize};

            std::vector<TF> prof_nogc(
                    p.second.data.begin() + gd.kstart,
                    p.second.data.begin() + gd.kstart + ksize);

            p.second.ncvar.insert(prof_nogc, time_height_index, time_height_size);
        }

        for (auto& t : col.time_series)
        {
            t.second.ncvar.insert(t.second.data, time_index);
        }

        // Synchronize the NetCDF file
        col.data_file->sync();
    }

    ++statistics_counter;
}

template<typename TF>
void Column<TF>::add_prof(std::string name, std::string longname, std::string unit, std::string zloc)
{
    auto& gd = grid.get_grid_data();

    // Create the NetCDF variable.
    for (auto& col : columns)
    {
        Prof_var var{col.data_file->template add_variable<TF>(name, {"time", zloc}),
            std::vector<TF>(gd.kcells)};

        var.ncvar.add_attribute("units", unit);
        var.ncvar.add_attribute("long_name", longname);

        // Insert the variable into the container.
        col.profs.emplace(
                std::piecewise_construct, std::forward_as_tuple(name), std::forward_as_tuple(std::move(var)));

        col.data_file->sync();
    }
}

template<typename TF>
void Column<TF>::add_time_series(std::string name, std::string longname, std::string unit)
{
    auto& gd = grid.get_grid_data();

    // Create the NetCDF variable.
    for (auto& col : columns)
    {
        Time_var var{col.data_file->template add_variable<TF>(name, {"time"}), TF(0)};

        var.ncvar.add_attribute("units", unit);
        var.ncvar.add_attribute("long_name", longname);

        // Insert the variable into the container.
        col.time_series.emplace(
                std::piecewise_construct, std::forward_as_tuple(name), std::forward_as_tuple(std::move(var)));

        col.data_file->sync();
    }
}

template<typename TF>
void Column<TF>::get_column_locations(std::vector<int>& i, std::vector<int>& j)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    for (auto& col : columns)
        if ( (col.coord[0] / gd.imax == md.mpicoordx ) && (col.coord[1] / gd.jmax == md.mpicoordy ) )
        {
            i.push_back(col.coord[0] % gd.imax + gd.istart);
            j.push_back(col.coord[1] % gd.jmax + gd.jstart);
        }
}

#ifndef USECUDA
template<typename TF>
void Column<TF>::calc_column(
        std::string profname, const TF* const restrict data, const TF offset, const bool copy_from_gpu)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    const int jj = gd.icells;
    const int kk = gd.ijcells;

    for (auto& col : columns)
    {
        // Check if coordinate is in range.
        if ( (col.coord[0] / gd.imax == md.mpicoordx ) && (col.coord[1] / gd.jmax == md.mpicoordy) )
        {
            const int i_col = col.coord[0] % gd.imax + gd.istart;
            const int j_col = col.coord[1] % gd.jmax + gd.jstart;

            for (int k=0; k<gd.kcells; k++)
            {
                const int ijk = i_col + j_col*jj + k*kk;
                col.profs.at(profname).data[k] = (data[ijk] + offset);
            }
        }
    }
}
#endif

template<typename TF>
void Column<TF>::set_individual_column(
        std::string profname, const TF* const restrict prof,
        const TF offset, const int i_col, const int j_col)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    for (auto& col : columns)
    {
        // Check if coordinate is in range.
        if ( (col.coord[0] / gd.imax == md.mpicoordx ) && (col.coord[1] / gd.jmax == md.mpicoordy) )
        {
            if (i_col == col.coord[0] % gd.imax + gd.istart &&
                j_col == col.coord[1] % gd.jmax + gd.jstart)
            {
                for (int k=0; k<gd.kcells; k++)
                    col.profs.at(profname).data[k] = (prof[k] + offset);
                return;
            }
        }
    }
    throw std::runtime_error("Cant set individual column for i,j=" + std::to_string(i_col) + "," + std::to_string(j_col));
}

#ifndef USECUDA
template<typename TF>
void Column<TF>::calc_time_series(
        std::string name, const TF* const restrict data, const TF offset)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    const int jj = gd.icells;
    const int kk = gd.ijcells;

    for (auto& col : columns)
    {
        // Check if coordinate is in range.
        if ( (col.coord[0] / gd.imax == md.mpicoordx ) && (col.coord[1] / gd.jmax == md.mpicoordy ) )
        {
            const int i_col = col.coord[0] % gd.imax + gd.istart;
            const int j_col = col.coord[1] % gd.jmax + gd.jstart;

            const int ij = i_col + j_col*jj;
            col.time_series.at(name).data = data[ij] + offset;
        }
    }
}
#endif

template class Column<double>;
template class Column<float>;
