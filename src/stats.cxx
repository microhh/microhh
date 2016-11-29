/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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
#include "stats.h"
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

Stats::Stats(Model* modelin, Input* inputin)
{
    model = modelin;
    master = model->master;

    // set the pointers to zero
    nmask  = 0;
    nmaskh = 0;

    int nerror = 0;
    nerror += inputin->get_item(&swstats, "stats", "swstats", "", "0");

    if (swstats == "1")
        nerror += inputin->get_item(&sampletime, "stats", "sampletime", "");

    if (!(swstats == "0" || swstats == "1"))
    {
        ++nerror;
        master->print_error("\"%s\" is an illegal value for swstats\n", swstats.c_str());
    }

    if (nerror)
        throw 1;
}

Stats::~Stats()
{
    delete[] nmask;
    delete[] nmaskh;

    // delete the profiles
    for (Mask_map::iterator it=masks.begin(); it!=masks.end(); ++it)
    {
        delete it->second.dataFile;
        for (Prof_map::const_iterator it2=it->second.profs.begin(); it2!=it->second.profs.end(); ++it2)
            delete[] it2->second.data;
    }
}

void Stats::init(double ifactor)
{
    // convenience pointers for short notation in class
    grid   = model->grid;
    fields = model->fields;
    master = model->master;

    // add the default mask
    add_mask("default");

    isampletime = (unsigned long)(ifactor * sampletime);

    nmask  = new int[grid->kcells];
    nmaskh = new int[grid->kcells];

    // set the number of stats to zero
    nstats = 0;
}

void Stats::create(int n)
{
    // do not create file if stats is disabled
    if (swstats == "0")
        return;

    int nerror = 0;

    for (Mask_map::iterator it=masks.begin(); it!=masks.end(); ++it)
    {
        // shortcut
        Mask* m = &it->second;

        // create a NetCDF file for the statistics
        if (master->mpiid == 0)
        {
            std::stringstream filename;
            filename << master->simname << "." << m->name << "." << std::setfill('0') << std::setw(7) << n << ".nc";

            try
            {
                m->dataFile = new NcFile(filename.str(), NcFile::newFile);
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
            m->z_dim  = m->dataFile->addDim("z" , grid->kmax);
            m->zh_dim = m->dataFile->addDim("zh", grid->kmax+1);
            m->t_dim  = m->dataFile->addDim("t");

            NcVar z_var;
            NcVar zh_var;

            // create variables belonging to dimensions
            m->iter_var = m->dataFile->addVar("iter", ncInt, m->t_dim);
            m->iter_var.putAtt("units", "-");
            m->iter_var.putAtt("long_name", "Iteration number");

            m->t_var = m->dataFile->addVar("t", ncDouble, m->t_dim);
            m->t_var.putAtt("units", "s");
            m->t_var.putAtt("long_name", "Time");

            z_var = m->dataFile->addVar("z", ncDouble, m->z_dim);
            z_var.putAtt("units", "m");
            z_var.putAtt("long_name", "Full level height");

            zh_var = m->dataFile->addVar("zh", ncDouble, m->zh_dim);
            zh_var.putAtt("units", "m");
            zh_var.putAtt("long_name", "Half level height");

            // save the grid variables
            z_var .putVar(&grid->z [grid->kstart]);
            zh_var.putVar(&grid->zh[grid->kstart]);

            // Synchronize the NetCDF file
            // BvS: only the last netCDF4-c++ includes the NcFile->sync()
            //      for now use sync() from the netCDF-C library to support older NetCDF4-c++ versions
            //m->dataFile->sync();
            nc_sync(m->dataFile->getId());
        }
    }

    // for each mask add the area as a variable
    add_prof("area" , "Fractional area contained in mask", "-", "z" );
    add_prof("areah", "Fractional area contained in mask", "-", "zh");
}

unsigned long Stats::get_time_limit(unsigned long itime)
{
    if (swstats == "0")
        return Constants::ulhuge;

    unsigned long idtlim = isampletime - itime % isampletime;
    return idtlim;
}

bool Stats::doStats()
{
    // check if stats are enabled
    if (swstats == "0")
        return false;

    // check if time for execution
    if (model->timeloop->get_itime() % isampletime != 0)
        return false;

    // return true such that stats are computed
    return true;
}

void Stats::exec(int iteration, double time, unsigned long itime)
{
    // This function is only called when stats are enabled no need for swstats check.

    // check if time for execution
    if (itime % isampletime != 0)
        return;

    // write message in case stats is triggered
    master->print_message("Saving stats for time %f\n", model->timeloop->get_time());

    for (Mask_map::iterator it=masks.begin(); it!=masks.end(); ++it)
    {
        // shortcut
        Mask* m = &it->second;

        // put the data into the NetCDF file
        if (master->mpiid == 0)
        {
            const std::vector<size_t> time_index = {static_cast<size_t>(nstats)};

            m->t_var   .putVar(time_index, &time     );
            m->iter_var.putVar(time_index, &iteration);

            const std::vector<size_t> time_height_index = {static_cast<size_t>(nstats), 0};
            std::vector<size_t> time_height_size  = {1, 0};

            for (Prof_map::const_iterator it=m->profs.begin(); it!=m->profs.end(); ++it)
            {
                time_height_size[1] = m->profs[it->first].ncvar.getDim(1).getSize();
                m->profs[it->first].ncvar.putVar(time_height_index, time_height_size, &m->profs[it->first].data[grid->kstart]);
            }

            for (Time_series_map::const_iterator it=m->tseries.begin(); it!=m->tseries.end(); ++it)
                m->tseries[it->first].ncvar.putVar(time_index, &m->tseries[it->first].data);

            // Synchronize the NetCDF file
            // BvS: only the last netCDF4-c++ includes the NcFile->sync()
            //      for now use sync() from the netCDF-C library to support older NetCDF4-c++ versions
            //m->dataFile->sync();
            nc_sync(m->dataFile->getId());
        }
    }

    ++nstats;
}

std::string Stats::get_switch()
{
    return swstats;
}

void Stats::add_mask(const std::string maskname)
{
    masks[maskname].name = maskname;
    masks[maskname].dataFile = 0;
}

void Stats::add_prof(std::string name, std::string longname, std::string unit, std::string zloc)
{
    // add the profile to all files
    for (Mask_map::iterator it=masks.begin(); it!=masks.end(); ++it)
    {
        // shortcut
        Mask* m = &it->second;

        // create the NetCDF variable
        if (master->mpiid == 0)
        {
            std::vector<NcDim> dim_vector = {m->t_dim};

            if (zloc == "z")
            {
                dim_vector.push_back(m->z_dim);
                m->profs[name].ncvar = m->dataFile->addVar(name, ncDouble, dim_vector);
                m->profs[name].data = NULL;
            }
            else if (zloc == "zh")
            {
                dim_vector.push_back(m->zh_dim);
                m->profs[name].ncvar = m->dataFile->addVar(name.c_str(), ncDouble, dim_vector);
                m->profs[name].data = NULL;
            }
            m->profs[name].ncvar.putAtt("units", unit.c_str());
            m->profs[name].ncvar.putAtt("long_name", longname.c_str());
            m->profs[name].ncvar.putAtt("_FillValue", ncDouble, NC_FILL_DOUBLE);
        }

        // and allocate the memory and initialize at zero
        m->profs[name].data = new double[grid->kcells];
        for (int k=0; k<grid->kcells; ++k)
            m->profs[name].data[k] = 0.;
    }
}

void Stats::add_fixed_prof(std::string name, std::string longname, std::string unit, std::string zloc, double* restrict prof)
{
    // add the profile to all files
    for (Mask_map::iterator it=masks.begin(); it!=masks.end(); ++it)
    {
        // shortcut
        Mask* m = &it->second;

        // create the NetCDF variable
        if (master->mpiid == 0)
        {
            NcVar var;
            if (zloc == "z")
                var = m->dataFile->addVar(name.c_str(), ncDouble, m->z_dim);
            else if (zloc == "zh")
                var = m->dataFile->addVar(name.c_str(), ncDouble, m->zh_dim);
            var.putAtt("units", unit.c_str());
            var.putAtt("long_name", longname.c_str());
            var.putAtt("_FillValue", ncDouble, NC_FILL_DOUBLE);

            const std::vector<size_t> index = {0};
            if (zloc == "z")
            {
                const std::vector<size_t> size  = {static_cast<size_t>(grid->kmax)};
                var.putVar(index, size, &prof[grid->kstart]);
            }
            else if (zloc == "zh")
            {
                const std::vector<size_t> size  = {static_cast<size_t>(grid->kmax+1)};
                var.putVar(index, size, &prof[grid->kstart]);
            }
        }
    }
}

void Stats::add_time_series(std::string name, std::string longname, std::string unit)
{
    // add the series to all files
    for (Mask_map::iterator it=masks.begin(); it!=masks.end(); ++it)
    {
        // shortcut
        Mask* m = &it->second;

        // create the NetCDF variable
        if (master->mpiid == 0)
        {
            m->tseries[name].ncvar = m->dataFile->addVar(name.c_str(), ncDouble, m->t_dim);
            m->tseries[name].ncvar.putAtt("units", unit.c_str());
            m->tseries[name].ncvar.putAtt("long_name", longname.c_str());
            m->tseries[name].ncvar.putAtt("_FillValue", ncDouble, NC_FILL_DOUBLE);
        }

        // Initialize at zero
        m->tseries[name].data = 0.;
    }
}

void Stats::get_mask(Field3d* mfield, Field3d* mfieldh, Mask* m)
{
    calc_mask(mfield->data, mfieldh->data, mfieldh->databot,
              nmask, nmaskh, &nmaskbot);
}

// COMPUTATIONAL KERNELS BELOW
void Stats::calc_mask(double* restrict mask, double* restrict maskh, double* restrict maskbot,
                      int* restrict nmask, int* restrict nmaskh, int* restrict nmaskbot)
{
    int ijtot = grid->itot*grid->jtot;

    // set all the mask values to 1
    for (int n=0; n<grid->ncells; ++n)
        mask[n] = 1.;

    for (int n=0; n<grid->ncells; ++n)
        maskh[n] = 1.;

    for (int n=0; n<grid->ijcells; ++n)
        maskbot[n] = 1.;

    for (int k=0; k<grid->kcells; ++k)
    {
        nmask [k] = ijtot;
        nmaskh[k] = ijtot;
    }
    *nmaskbot = ijtot;
}

void Stats::calc_area(double* restrict area, const int loc[3], int* restrict nmask)
{
    const int ijtot = grid->itot*grid->jtot;

    for (int k=grid->kstart; k<grid->kend+loc[2]; k++)
    {
        if (nmask[k] > nthres)
            area[k] = (double)(nmask[k]) / (double)ijtot;
        else
            area[k] = 0.;
    }
}

void Stats::calc_mean(double* const restrict prof, const double* const restrict data,
                      const double offset, const int loc[3],
                      const double* const restrict mask, const int * const restrict nmask)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    for (int k=1; k<grid->kcells; k++)
    {
        prof[k] = 0.;
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk  = i + j*jj + k*kk;
                prof[k] += mask[ijk]*(data[ijk] + offset);
            }
    }

    master->sum(prof, grid->kcells);

    for (int k=1; k<grid->kcells; k++)
    {
        if (nmask[k] > nthres)
            prof[k] /= (double)(nmask[k]);
        else
            prof[k] = NC_FILL_DOUBLE;
    }
}

void Stats::calc_mean2d(double* const restrict mean, const double* const restrict data,
                        const double offset,
                        const double* const restrict mask, const int * const restrict nmask)
{
    const int jj = grid->icells;

    if (*nmask > nthres)
    {
        *mean = 0.;
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ij = i + j*jj;
                *mean += mask[ij]*(data[ij] + offset);
            }
        master->sum(mean,1);
        *mean /= (double)*nmask;
    }
    else
        *mean = NC_FILL_DOUBLE;
}

void Stats::calc_sorted_prof(double* restrict data, double* restrict bin, double* restrict prof)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;
    const int kend = grid->kend;

    double minval =  Constants::dhuge;
    double maxval = -Constants::dhuge;

    // first, get min and max
    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj + k*kk;
                if (data[ijk] < minval)
                    minval = data[ijk];
                if (data[ijk] > maxval)
                    maxval = data[ijk];
            }

    master->min(&minval, 1);
    master->max(&maxval, 1);

    // make sure that the max ends up in the last bin (introduce 1E-9 error)
    maxval *= (1.+Constants::dsmall);

    const double range = maxval-minval;

    // In case the field is entirely uniform, dbin becomes zero. In that case we set the profile to the minval.
    if (range < 1.e-16)
    {
        for (int k=grid->kstart; k<grid->kend; ++k)
            prof[k] = minval;
    }
    else
    {
        // create bins, equal to the number of grid cells per proc
        // make sure that bins is not larger than the memory of one 3d field
        const int bins = grid->nmax;

        // calculate bin width, subtract one to make the minimum and maximum
        // are in the middle of the bin range and add half a bin size on both sides
        // |----x----|----x----|----x----|
        const double dbin = range / (double)(bins-1);

        minval -= 0.5*dbin;
        maxval += 0.5*dbin;

        // set the bin array to zero
        for (int n=0; n<bins; ++n)
            bin[n] = 0;

        // calculate the division factor of one equivalent height unit
        // (the total volume saved is itot*jtot*zsize)
        const double nslice = (double)(grid->itot*grid->jtot);

        // check in which bin each value falls and increment the bin count
        for (int k=grid->kstart; k<grid->kend; ++k)
        {
            const double dzslice = grid->dz[k] / nslice;
            for (int j=grid->jstart; j<grid->jend; ++j)
                // do not add a ivdep pragma here, because multiple instances could write the same bin[index]
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const int index = (int)((data[ijk] - minval) / dbin);
                    bin[index] += dzslice;
                }
        }

        // get the bin count
        master->sum(bin, bins);

        // set the starting values of the loop
        int index = 0;
        double zbin = 0.5*bin[index];
        double profval = minval + 0.5*dbin;

        for (int k=grid->kstart; k<grid->kend; ++k)
        {
            // Integrate the profile up to the bin count.
            // Escape the while loop when the integrated profile
            // exceeds the next grid point.
            while (zbin < grid->z[k])
            {
                zbin += 0.5*(bin[index]+bin[index+1]);
                profval += dbin;
                ++index;
            }

            // In case the first bin is larger than the grid spacing, which can happen
            // in the inital phase of an MPI run, make sure that no out-of-bounds reads
            // happen.
            if (index == 0)
                prof[k] = profval;
            else
            {
                const double dzfrac = (zbin-grid->z[k]) / (0.5*(bin[index-1]+bin[index]));
                prof[k] = profval - dzfrac*dbin;
            }
        }
    }

    // now calculate the ghost cells
    // \TODO this might not be accurate enough, extrapolate properly
    double profbot = minval;
    double proftop = maxval;

    if (grid->swspatialorder == "2")
    {
        prof[kstart-1] = 2.*profbot - prof[kstart];
        prof[kend]     = 2.*proftop - prof[kend-1];
    }
    else if (grid->swspatialorder == "4")
    {
        prof[kstart-1] = (8./3.)*profbot - 2.*prof[kstart] + (1./3.)*prof[kstart+1];
        prof[kstart-2] = 8.*profbot      - 9.*prof[kstart] + 2.*prof[kstart+1];
        prof[kend]     = (8./3.)*proftop - 2.*prof[kend-1] + (1./3.)*prof[kend-2];
        prof[kend+1]   = 8.*proftop      - 9.*prof[kend-1] + 2.*prof[kend-2];
    }
}

// \TODO the count function assumes that the variable to count is at the mask location
void Stats::calc_count(double* restrict data, double* restrict prof, double threshold,
                       double* restrict mask, int* restrict nmask)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    for (int k=0; k<grid->kcells; ++k)
    {
        prof[k] = 0.;
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                if (data[ijk] > threshold)
                    prof[k] += mask[ijk]*1.;
            }
    }

    master->sum(prof, grid->kcells);

    for (int k=0; k<grid->kcells; k++)
    {
        if (nmask[k] > nthres)
            prof[k] /= (double)(nmask[k]);
        else
            prof[k] = NC_FILL_DOUBLE;
    }
}

void Stats::calc_moment(double* restrict data, double* restrict datamean, double* restrict prof, double power, const int loc[3],
                        double* restrict mask, int* restrict nmask)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    for (int k=grid->kstart; k<grid->kend+1; ++k)
    {
        prof[k] = 0.;
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                prof[k] += mask[ijk]*std::pow(data[ijk]-datamean[k], power);
            }
    }

    master->sum(prof, grid->kcells);

    for (int k=1; k<grid->kcells; k++)
    {
        if (nmask[k] > nthres)
            prof[k] /= (double)(nmask[k]);
        else
            prof[k] = NC_FILL_DOUBLE;
    }
}

void Stats::calc_flux_2nd(double* restrict data, double* restrict datamean, double* restrict w, double* restrict wmean,
                          double* restrict prof, double* restrict tmp1, const int loc[3],
                          double* restrict mask, int* restrict nmask)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    // set a pointer to the field that contains w, either interpolated or the original
    double* restrict calcw = w;

    // define the locations
    const int wloc [3] = {0,0,1};
    const int uwloc[3] = {1,0,1};
    const int vwloc[3] = {0,1,1};

    if (loc[0] == 1)
    {
        grid->interpolate_2nd(tmp1, w, wloc, uwloc);
        calcw = tmp1;
    }
    else if (loc[1] == 1)
    {
        grid->interpolate_2nd(tmp1, w, wloc, vwloc);
        calcw = tmp1;
    }

    for (int k=grid->kstart; k<grid->kend+1; ++k)
    {
        prof[k] = 0.;
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk  = i + j*jj + k*kk;
                prof[k] += mask[ijk]*(0.5*(data[ijk-kk]+data[ijk])-0.5*(datamean[k-1]+datamean[k]))*(calcw[ijk]-wmean[k]);
                // prof[k] += mask[ijk]*0.5*(data[ijk-kk]+data[ijk])*calcw[ijk];
            }
    }

    master->sum(prof, grid->kcells);

    for (int k=1; k<grid->kcells; ++k)
    {
        if (nmask[k] > nthres && datamean[k-1] != NC_FILL_DOUBLE && datamean[k] != NC_FILL_DOUBLE)
            prof[k] /= (double)(nmask[k]);
        else
            prof[k] = NC_FILL_DOUBLE;
    }
}

void Stats::calc_flux_4th(double* restrict data, double* restrict w, double* restrict prof, double* restrict tmp1, const int loc[3],
                          double* restrict mask, int* restrict nmask)
{
    using namespace Finite_difference::O4;

    const int jj  = 1*grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    // set a pointer to the field that contains w, either interpolated or the original
    double* restrict calcw = w;

    // define the locations
    const int wloc [3] = {0,0,1};
    const int uwloc[3] = {1,0,1};
    const int vwloc[3] = {0,1,1};

    if (loc[0] == 1)
    {
        grid->interpolate_4th(tmp1, w, wloc, uwloc);
        calcw = tmp1;
    }
    else if (loc[1] == 1)
    {
        grid->interpolate_4th(tmp1, w, wloc, vwloc);
        calcw = tmp1;
    }

    for (int k=grid->kstart; k<grid->kend+1; ++k)
    {
        prof[k] = 0.;
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk1;
                prof[k] += mask[ijk]*(ci0*data[ijk-kk2] + ci1*data[ijk-kk1] + ci2*data[ijk] + ci3*data[ijk+kk1])*calcw[ijk];
            }
    }

    master->sum(prof, grid->kcells);

    for (int k=1; k<grid->kcells; k++)
    {
        if (nmask[k] > nthres)
            prof[k] /= (double)(nmask[k]);
        else
            prof[k] = NC_FILL_DOUBLE;
    }
}

void Stats::calc_grad_2nd(double* restrict data, double* restrict prof, double* restrict dzhi, const int loc[3],
                          double* restrict mask, int* restrict nmask)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    for (int k=grid->kstart; k<grid->kend+1; ++k)
    {
        prof[k] = 0.;
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                prof[k] += mask[ijk]*(data[ijk]-data[ijk-kk])*dzhi[k];
            }
    }

    master->sum(prof, grid->kcells);

    for (int k=1; k<grid->kcells; k++)
    {
        if (nmask[k] > nthres)
            prof[k] /= (double)(nmask[k]);
        else
            prof[k] = NC_FILL_DOUBLE;
    }
}

void Stats::calc_grad_4th(double* restrict data, double* restrict prof, double* restrict dzhi4, const int loc[3],
                          double* restrict mask, int* restrict nmask)
{
    using namespace Finite_difference::O4;

    const int jj  = 1*grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    for (int k=grid->kstart; k<grid->kend+1; ++k)
    {
        prof[k] = 0.;
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk1;
                prof[k] += mask[ijk]*(cg0*data[ijk-kk2] + cg1*data[ijk-kk1] + cg2*data[ijk] + cg3*data[ijk+kk1])*dzhi4[k];
            }
    }

    master->sum(prof, grid->kcells);

    for (int k=1; k<grid->kcells; k++)
    {
        if (nmask[k] > nthres)
            prof[k] /= (double)(nmask[k]);
        else
            prof[k] = NC_FILL_DOUBLE;
    }
}

void Stats::calc_diff_4th(double* restrict data, double* restrict prof, double* restrict dzhi4, double visc, const int loc[3],
                          double* restrict mask, int* restrict nmask)
{
    using namespace Finite_difference::O4;

    const int jj  = 1*grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    for (int k=grid->kstart; k<grid->kend+1; ++k)
    {
        prof[k] = 0.;
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk1;
                prof[k] -= mask[ijk]*visc*(cg0*data[ijk-kk2] + cg1*data[ijk-kk1] + cg2*data[ijk] + cg3*data[ijk+kk1])*dzhi4[k];
            }
    }

    master->sum(prof, grid->kcells);

    for (int k=1; k<grid->kcells; k++)
    {
        if (nmask[k] > nthres)
            prof[k] /= (double)(nmask[k]);
        else
            prof[k] = NC_FILL_DOUBLE;
    }
}

void Stats::calc_diff_2nd(double* restrict data, double* restrict prof, double* restrict dzhi, double visc, const int loc[3],
                          double* restrict mask, int* restrict nmask)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    for (int k=grid->kstart; k<grid->kend+1; ++k)
    {
        prof[k] = 0.;
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                prof[k] -= mask[ijk]*visc*(data[ijk] - data[ijk-kk])*dzhi[k];
            }
    }

    master->sum(prof, grid->kcells);

    for (int k=1; k<grid->kcells; k++)
    {
        if (nmask[k] > nthres)
            prof[k] /= (double)(nmask[k]);
        else
            prof[k] = NC_FILL_DOUBLE;
    }
}


void Stats::calc_diff_2nd(double* restrict data, double* restrict w, double* restrict evisc,
                          double* restrict prof, double* restrict dzhi,
                          double* restrict fluxbot, double* restrict fluxtop, double tPr, const int loc[3],
                          double* restrict mask, int* restrict nmask)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;
    const int kend = grid->kend;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    // bottom boundary
    prof[kstart] = 0.;
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            prof[kstart] += mask[ijk]*fluxbot[ij];
        }

    // calculate the interior
    if (loc[0] == 1)
    {
        for (int k=grid->kstart+1; k<grid->kend; ++k)
        {
            prof[k] = 0.;
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk  = i + j*jj + k*kk;
                    // evisc * (du/dz + dw/dx)
                    const double eviscu = 0.25*(evisc[ijk-ii-kk]+evisc[ijk-ii]+evisc[ijk-kk]+evisc[ijk]);
                    prof[k] += -mask[ijk]*eviscu*( (data[ijk]-data[ijk-kk])*dzhi[k] + (w[ijk]-w[ijk-ii])*dxi );
                }
        }
    }
    else if (loc[1] == 1)
    {
        for (int k=grid->kstart+1; k<grid->kend; ++k)
        {
            prof[k] = 0.;
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    // evisc * (dv/dz + dw/dy)
                    const double eviscv = 0.25*(evisc[ijk-jj-kk]+evisc[ijk-jj]+evisc[ijk-kk]+evisc[ijk]);
                    prof[k] += -mask[ijk]*eviscv*( (data[ijk]-data[ijk-kk])*dzhi[k] + (w[ijk]-w[ijk-jj])*dyi );
                }
        }
    }
    else
    {
        for (int k=grid->kstart+1; k<grid->kend; ++k)
        {
            prof[k] = 0.;
            for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
                for (int i=grid->istart; i<grid->iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const double eviscs = 0.5*(evisc[ijk-kk]+evisc[ijk])/tPr;
                    prof[k] += -mask[ijk]*eviscs*(data[ijk]-data[ijk-kk])*dzhi[k];
                }
        }
    }

    // top boundary
    prof[kend] = 0.;
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kend*kk;
            prof[kend] += mask[ijk]*fluxtop[ij];
        }

    master->sum(prof, grid->kcells);

    for (int k=1; k<grid->kcells; k++)
    {
        if (nmask[k] > nthres)
            prof[k] /= (double)(nmask[k]);
        else
            prof[k] = NC_FILL_DOUBLE;
    }
}

void Stats::add_fluxes(double* restrict flux, double* restrict turb, double* restrict diff)
{
    for (int k=grid->kstart; k<grid->kend+1; ++k)
    {
        if (turb[k] == NC_FILL_DOUBLE || diff[k] == NC_FILL_DOUBLE)
            flux[k] = NC_FILL_DOUBLE;
        else
            flux[k] = turb[k] + diff[k];
    }
}

/**
 * This function calculates the total domain integrated path of variable data over maskbot
 */
void Stats::calc_path(double* restrict data, double* restrict maskbot, int* restrict nmaskbot, double* restrict path)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    *path = 0.;

    if (*nmaskbot > nthres)
    {
        // Integrate liquid water
        for (int j=grid->jstart; j<grid->jend; j++)
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ij = i + j*jj;
                if (maskbot[ij] == 1)
                    for (int k=kstart; k<grid->kend; k++)
                    {
                        const int ijk = i + j*jj + k*kk;
                        *path += fields->rhoref[k] * data[ijk] * grid->dz[k];
                    }
            }
        *path /= (double)*nmaskbot;
        master->sum(path, 1);
    }
    else
        *path = NC_FILL_DOUBLE;
}

/**
 * This function calculates the vertical projected cover of variable data over maskbot
 */
void Stats::calc_cover(double* restrict data, double* restrict maskbot, int* restrict nmaskbot, double* restrict cover, double threshold)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    *cover = 0.;

    if (*nmaskbot > nthres)
    {
        // Per column, check if cloud present
        for (int j=grid->jstart; j<grid->jend; j++)
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ij = i + j*jj;
                if (maskbot[ij] == 1)
                    for (int k=kstart; k<grid->kend; k++)
                    {
                        const int ijk = i + j*jj + k*kk;
                        if (data[ijk]>threshold)
                        {
                            *cover += 1.;
                            break;
                        }
                    }
            }
        *cover /= (double)*nmaskbot;
        grid->get_prof(cover,1);
    }
    else
        *cover = NC_FILL_DOUBLE;
}
