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

#include <cstdio>
#include <cmath>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "column.h"
//#include "thermo_moist.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
//#include "diff_smag2.h"
#include "timeloop.h"

#include <netcdf>       // C++
#include <netcdf.h>     // C, for sync() using older netCDF-C++ versions
using namespace netCDF;
using namespace netCDF::exceptions;


namespace
{
    // Help functions to switch between the different NetCDF data types
    template<typename TF> NcType netcdf_fp_type();
    template<> NcType netcdf_fp_type<double>() { return ncDouble; }
    template<> NcType netcdf_fp_type<float>()  { return ncFloat; }

    template<typename TF> TF netcdf_fp_fillvalue();
    template<> double netcdf_fp_fillvalue<double>() { return NC_FILL_DOUBLE; }
    template<> float  netcdf_fp_fillvalue<float>()  { return NC_FILL_FLOAT; }
}


template<typename TF>
Column<TF>::Column(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin):
    master(masterin), grid(gridin), fields(fieldsin)
{
    swcolumn = inputin.get_item<bool>("column", "swcolumn", "", false);

    if (swcolumn)
    {
        sampletime = inputin.get_item<double>("column", "sampletime", "");
    }
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
void Column<TF>::create(int iotime, std::string sim_name)
{
    // do not create file if column is disabled
    if (!swcolumn)
        return;

    int nerror = 0;
    auto& gd = grid.get_grid_data();

    // create a NetCDF file for the statistics
    if (master.get_mpiid() == 0)
    {
        std::stringstream filename;
        filename << sim_name << "." << "column" << "." << std::setfill('0') << std::setw(7) << iotime << ".nc";

        try
        {
            data_file = new NcFile(filename.str(), NcFile::newFile);
        }
        catch(NcException& e)
        {
            master.print_error("NetCDF exception: %s\n",e.what());
            ++nerror;
        }
    }

    // Crash on all processes in case the file could not be written
    master.broadcast(&nerror, 1);
    if (nerror)
        throw 1;

    // create dimensions
    if (master.get_mpiid() == 0)
    {
        z_dim  = data_file->addDim("z" , gd.kmax);
        zh_dim = data_file->addDim("zh", gd.kmax+1);
        t_dim  = data_file->addDim("t");

        NcVar z_var;
        NcVar zh_var;

        // create variables belonging to dimensions
        iter_var = data_file->addVar("iter", ncInt, t_dim);
        iter_var.putAtt("units", "-");
        iter_var.putAtt("long_name", "Iteration number");

        t_var = data_file->addVar("time", ncDouble, t_dim);
        t_var.putAtt("units", "s");
        t_var.putAtt("long_name", "Time");

        z_var = data_file->addVar("z", ncDouble, z_dim);
        z_var.putAtt("units", "m");
        z_var.putAtt("long_name", "Full level height");

        zh_var = data_file->addVar("zh", ncDouble, zh_dim);
        zh_var.putAtt("units", "m");
        zh_var.putAtt("long_name", "Half level height");

        // save the grid variables
        z_var .putVar(&gd.z [gd.kstart]);
        zh_var.putVar(&gd.zh[gd.kstart]);

        // Synchronize the NetCDF file
        // BvS: only the last netCDF4-c++ includes the NcFile->sync()
        //      for now use sync() from the netCDF-C library to support older NetCDF4-c++ versions
        //data_file->sync();
        nc_sync(data_file->getId());
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

    auto& gd = grid.get_grid_data();
    // write message in case column is triggered
    master.print_message("Saving column for time %f\n", time);

    // put the data into the NetCDF file
    if (master.get_mpiid() == 0)
    {
        const std::vector<size_t> time_index = {static_cast<size_t>(statistics_counter)};

        t_var   .putVar(time_index, &time     );
        iter_var.putVar(time_index, &iteration);
        const std::vector<size_t> time_height_index = {static_cast<size_t>(statistics_counter), 0};
        std::vector<size_t> time_height_size  = {1, 0};
        for (auto& p : profs)
        {
            time_height_size[1] = p.second.ncvar.getDim(1).getSize();
            p.second.ncvar.putVar(time_height_index, time_height_size, &p.second.data.data()[gd.kstart]);
        }
        // Synchronize the NetCDF file
        // BvS: only the last netCDF4-c++ includes the NcFile->sync()
        //      for now use sync() from the netCDF-C library to support older NetCDF4-c++ versions
        //data_file->sync();
        nc_sync(data_file->getId());
    }

    ++statistics_counter;
}

template<typename TF>
void Column<TF>::add_prof(std::string name, std::string longname, std::string unit, std::string zloc)
{
    auto& gd = grid.get_grid_data();
    // create the NetCDF variable
    if (master.get_mpiid() == 0)
    {
        std::vector<NcDim> dim_vector = {t_dim};

        if (zloc == "z")
        {
            dim_vector.push_back(z_dim);
            profs[name].ncvar = data_file->addVar(name.c_str(), netcdf_fp_type<TF>(), dim_vector);
        }
        else if (zloc == "zh")
        {
            dim_vector.push_back(zh_dim);
            profs[name].ncvar = data_file->addVar(name.c_str(), netcdf_fp_type<TF>(), dim_vector);
        }
        profs.at(name).ncvar.putAtt("units", unit.c_str());
        profs.at(name).ncvar.putAtt("long_name", longname.c_str());
        profs.at(name).ncvar.putAtt("_FillValue", netcdf_fp_type<TF>(), netcdf_fp_fillvalue<TF>());

        nc_sync(data_file->getId());
    }
    profs[name].data.resize(gd.kcells);

}

template<typename TF>
void Column<TF>::calc_column(std::string profname, const TF* const restrict data,
                      const TF offset)
{
    auto& gd = grid.get_grid_data();
    const int jj = gd.icells;
    const int kk = gd.ijcells;

    for (int k=1; k<gd.kcells; k++)
    {
        const int ijk  = gd.istart + gd.jstart*jj + k*kk;
        profs.at(profname).data.data()[k] = (data[ijk] + offset);
    }

}

template class Column<double>;
template class Column<float>;
