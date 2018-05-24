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
#include "thermo_moist.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "model.h"
#include "diff_smag2.h"
#include "timeloop.h"

#include <netcdf>       // C++
#include <netcdf.h>     // C, for sync() using older netCDF-C++ versions
using namespace netCDF;
using namespace netCDF::exceptions;

Column::Column(Model* modelin, Input* inputin)
{
    model = modelin;
    master = model->master;

    // set the pointers to zero
    dataFile = 0;

    int nerror = 0;
    nerror += inputin->get_item(&swcolumn, "column", "swcolumn", "", "0");

    if (swcolumn == "1")
        nerror += inputin->get_item(&sampletime, "column", "sampletime", "");

    if (!(swcolumn == "0" || swcolumn == "1"))
    {
        ++nerror;
        master->print_error("\"%s\" is an illegal value for swcolumn\n", swcolumn.c_str());
    }

    if (nerror)
        throw 1;
}

Column::~Column()
{
    if (swcolumn == "0")
        return;

    delete dataFile;
    for (Column_map::const_iterator it=profs.begin(); it!=profs.end(); ++it)
        delete[] it->second.data;
}

void Column::init(double ifactor)
{
    // do not create file if column is disabled
    if (swcolumn == "0")
        return;

    // convenience pointers for short notation in class
    grid   = model->grid;
    fields = model->fields;
    master = model->master;

    isampletime = (unsigned long)(ifactor * sampletime);

    // set the number of column to zero
    ncolumn = 0;
}

void Column::create(int n)
{
    // do not create file if column is disabled
    if (swcolumn == "0")
        return;

    int nerror = 0;

        // create a NetCDF file for the statistics
        if (master->mpiid == 0)
        {
            std::stringstream filename;
            filename << master->simname << "." << "column" << "." << std::setfill('0') << std::setw(7) << n << ".nc";

            try
            {
                dataFile = new NcFile(filename.str(), NcFile::newFile);
            }
            catch(NcException& e)
            {
                master->print_error("NetCDF exception: %s\n",e.what());
                ++nerror;
            }
        }

        // Crash on all processes in case the file could not be written
        master->broadcast(&nerror, 1);
        if (nerror)
            throw 1;

        // create dimensions
        if (master->mpiid == 0)
        {
            z_dim  = dataFile->addDim("z" , grid->kmax);
            zh_dim = dataFile->addDim("zh", grid->kmax+1);
            t_dim  = dataFile->addDim("t");

            NcVar z_var;
            NcVar zh_var;

            // create variables belonging to dimensions
            iter_var = dataFile->addVar("iter", ncInt, t_dim);
            iter_var.putAtt("units", "-");
            iter_var.putAtt("long_name", "Iteration number");

            t_var = dataFile->addVar("t", ncDouble, t_dim);
            t_var.putAtt("units", "s");
            t_var.putAtt("long_name", "Time");

            z_var = dataFile->addVar("z", ncDouble, z_dim);
            z_var.putAtt("units", "m");
            z_var.putAtt("long_name", "Full level height");

            zh_var = dataFile->addVar("zh", ncDouble, zh_dim);
            zh_var.putAtt("units", "m");
            zh_var.putAtt("long_name", "Half level height");

            // save the grid variables
            z_var .putVar(&grid->z [grid->kstart]);
            zh_var.putVar(&grid->zh[grid->kstart]);

            // Synchronize the NetCDF file
            // BvS: only the last netCDF4-c++ includes the NcFile->sync()
            //      for now use sync() from the netCDF-C library to support older NetCDF4-c++ versions
            //dataFile->sync();
            nc_sync(dataFile->getId());
        }

}

unsigned long Column::get_time_limit(unsigned long itime)
{
    if (swcolumn == "0")
        return Constants::ulhuge;
    // If sampletime is negative, output column every timestep
    if (isampletime < 0)
        return Constants::ulhuge;

    unsigned long idtlim = isampletime - itime % isampletime;
    return idtlim;
}

bool Column::doColumn()
{
    // check if column are enabled
    if (swcolumn == "0")
        return false;

    // If sampletime is negative, output column every timestep
    if (isampletime < 0)
        return true;

    // check if time for execution
    if (model->timeloop->get_itime() % isampletime != 0)
        return false;

    // return true such that column are computed
    return true;
}

void Column::exec(int iteration, double time, unsigned long itime)
{
    // This function is only called when column are enabled no need for swcolumn check.

    // write message in case column is triggered
    master->print_message("Saving column for time %f\n", time);

    // put the data into the NetCDF file
    if (master->mpiid == 0)
    {
        const std::vector<size_t> time_index = {static_cast<size_t>(ncolumn)};

        t_var   .putVar(time_index, &time     );
        iter_var.putVar(time_index, &iteration);

        const std::vector<size_t> time_height_index = {static_cast<size_t>(ncolumn), 0};
        std::vector<size_t> time_height_size  = {1, 0};

        for (Column_map::const_iterator it=profs.begin(); it!=profs.end(); ++it)
        {
            time_height_size[1] = profs[it->first].ncvar.getDim(1).getSize();
            profs[it->first].ncvar.putVar(time_height_index, time_height_size, &profs[it->first].data[grid->kstart]);
        }

        // Synchronize the NetCDF file
        // BvS: only the last netCDF4-c++ includes the NcFile->sync()
        //      for now use sync() from the netCDF-C library to support older NetCDF4-c++ versions
        //dataFile->sync();
        nc_sync(dataFile->getId());
    }
    
    ++ncolumn;
}

std::string Column::get_switch()
{
    return swcolumn;
}


void Column::add_prof(std::string name, std::string longname, std::string unit, std::string zloc)
{
    // create the NetCDF variable
    if (master->mpiid == 0)
    {
        std::vector<NcDim> dim_vector = {t_dim};

        if (zloc == "z")
        {
            dim_vector.push_back(z_dim);
            profs[name].ncvar = dataFile->addVar(name, ncDouble, dim_vector);
            profs[name].data = NULL;
        }
        else if (zloc == "zh")
        {
            dim_vector.push_back(zh_dim);
            profs[name].ncvar = dataFile->addVar(name.c_str(), ncDouble, dim_vector);
            profs[name].data = NULL;
        }
        profs[name].ncvar.putAtt("units", unit.c_str());
        profs[name].ncvar.putAtt("long_name", longname.c_str());
        profs[name].ncvar.putAtt("_FillValue", ncDouble, NC_FILL_DOUBLE);
    }

    // and allocate the memory and initialize at zero
    profs[name].data = new double[grid->kcells];
    for (int k=0; k<grid->kcells; ++k)
        profs[name].data[k] = 0.;



}
void Column::calc_column(double* const restrict prof, const double* const restrict data,
                      const double offset)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;
   

    for (int k=1; k<grid->kcells; k++)
    {
        const int ijk  = grid->istart + grid->jstart*jj + k*kk;
        prof[k] = (data[ijk] + offset);
    }
}

