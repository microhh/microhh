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
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <netcdf>       // C++
#include <netcdf.h>     // C, for sync() using older netCDF-C++ versions
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
//#include "thermo_moist.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
//#include "diff_smag2.h"
#include "timeloop.h"

using namespace netCDF;
using namespace netCDF::exceptions;
using namespace Constants;

namespace
{
    // Help functions to switch between the different NetCDF data types
    template<typename TF> NcType netcdf_fp_type();
    template<> NcType netcdf_fp_type<double>() { return ncDouble; }
    template<> NcType netcdf_fp_type<float>()  { return ncFloat; }

    template<typename TF> TF netcdf_fp_fillvalue();
    template<> double netcdf_fp_fillvalue<double>() { return NC_FILL_DOUBLE; }
    template<> float  netcdf_fp_fillvalue<float>()  { return NC_FILL_FLOAT; }


    template<typename TF, Stats_mask_type mode>
    TF compare(const TF value, const TF threshold)
    {
        if (mode == Stats_mask_type::Plus)
            return (value > threshold);
        else if (mode == Stats_mask_type::Min)
            return (value <= threshold);
    }

    template<typename TF, Stats_mask_type mode>
    void calc_mask_thres(unsigned int* const restrict mfield, unsigned int* const restrict mfield_bot, const unsigned int flag, const unsigned int flagh,
                        const TF* const restrict fld, const TF* const restrict fldh,
                        const TF* const restrict fld_bot, const TF threshold,
                        const int istart, const int jstart, const int kstart,
                        const int iend,   const int jend,   const int kend,
                        const int icells, const int ijcells)
    {

        #pragma omp parallel for
        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    mfield[ijk] = (mfield[ijk] & flag ) * compare<TF, mode>(fld[ijk], threshold);
                    mfield[ijk] = (mfield[ijk] & flagh) * compare<TF, mode>(fldh[ijk], threshold);
                }

        // Set the mask for surface projected quantities
        #pragma omp parallel for
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij  = i + j*icells + kstart*ijcells;
                const int ij2 = i + j*icells + (kend+1)*ijcells;
                mfield_bot[ij] = (mfield_bot[ij] & flag)*compare<TF, mode>(fld_bot[ij], threshold);
                mfield[ij2] = (mfield[ij2] & flagh) * compare<TF, mode>(fldh[ij2], threshold);
            }
    }


    template<typename TF, Stats_mask_type mode>
    void calc_mask_thres_pert(unsigned int* const restrict mfield, unsigned int* const restrict mfield_bot, const unsigned int flag, const unsigned int flagh,
                        const TF* const restrict fld, const TF* const restrict fld_mean, const TF* const restrict fldh,
                        const TF* const restrict fldh_mean, const TF* const restrict fld_bot, const TF threshold,
                        const int istart, const int jstart, const int kstart,
                        const int iend,   const int jend,   const int kend,
                        const int icells, const int ijcells)
    {

        #pragma omp parallel for
        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    mfield[ijk] = (mfield[ijk] & flag ) * compare<TF, mode>(fld[ijk]-fld_mean[k], threshold);
                    mfield[ijk] = (mfield[ijk] & flagh) * compare<TF, mode>(fldh[ijk]-fldh_mean[k], threshold);
                }

        // Set the mask for surface projected quantities
        #pragma omp parallel for
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij  = i + j*icells + kstart*ijcells;
                const int ij2 = i + j*icells + (kend+1)*ijcells;
                mfield_bot[ij] = (mfield_bot[ij] & flag)*compare<TF, mode>(fld_bot[ij]-fld_mean[kstart], threshold);
                mfield[ij2] = (mfield[ij2] & flagh) * compare<TF, mode>(fldh[ij2]-fldh_mean[kend], threshold);
            }
    }

    template<typename TF>
    void calc_area(TF* const restrict area, const int loc[3], const int* const restrict nmask, const int kstart, const int kend, const int ijtot)
    {
        for (int k=kstart; k<kend+loc[2]; k++)
        {
            if (nmask[k])
                area[k] = static_cast<TF>(nmask[k]) / static_cast<TF>(ijtot);
            else
                area[k] = 0.;
        }
    }

    // Sets all the mask values to one (non-masked field)
    template<typename TF>
    void calc_nmask(int* restrict nmask_full, int* restrict nmask_half, int& nmask_bottom,
                    const unsigned int* const mfield, const unsigned int* const mfield_bot,const int flag,const int flagh,
                   const int istart, const int iend, const int jstart, const int jend,
                   const int kstart, const int kend, const int icells, const int ijcells, const int kcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
        {
            nmask_full[k] = 0;
            nmask_half[k] = 0;

            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk  = i + j*icells + k*ijcells;
                    nmask_full[k]+=(mfield[ijk] & flag);
                    nmask_half[k]+=(mfield[ijk] & flagh);
                }
        }

        nmask_bottom = 0;
        for (int j=jstart; j<jend; ++j)
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                nmask_bottom+=(mfield_bot[ij] & flag);
            }
    }

    template<typename TF>
    void calc_mean(TF* const restrict prof, const TF* const restrict fld, const TF offset,
                    const unsigned int* const mask, const unsigned int flag, const int* const nmask,
                    const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend+1; k++)
        {
            if (nmask[k])
            {
                prof[k] = 0.;
                for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk  = i + j*icells + k*ijcells;
                    prof[k] += (mask[ijk] & flag)*(fld[ijk] + offset);
                }
                prof[k] /= static_cast<TF>(nmask[k]);
            }
            else
                prof[k] = netcdf_fp_fillvalue<TF>();
        }
    }

    template<typename TF>
    void calc_moment(TF* const restrict prof, const TF* const restrict fld, const TF* const restrict fld_mean, const TF offset,
                    const unsigned int* const mask, const unsigned int flag, const int* const nmask, const int power,
                    const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend, const int icells, const int ijcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend+1; k++)
        {
            if (nmask[k])
            {
                prof[k] = 0.;
                for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk  = i + j*icells + k*ijcells;
                    prof[k] += static_cast<TF>(mask[ijk] & flag)*std::pow(fld[ijk] - fld_mean[k] + offset, power);
                }
                prof[k] /= static_cast<TF>(nmask[k]);
            }
            else
                prof[k] = netcdf_fp_fillvalue<TF>();
        }
    }

    template<typename TF>
    void add_fluxes(TF* restrict flux, TF* restrict turb, TF* restrict diff, const TF fillvalue, const int kstart, const int kend)
    {
        for (int k=kstart; k<kend+1; ++k)
        {
            if (turb[k] == fillvalue || diff[k] == fillvalue)
                flux[k] = fillvalue;
            else
                flux[k] = turb[k] + diff[k];
        }
    }

    bool has_only_digits(const std::string s)
    {
        return s.find_first_not_of( "23456789" ) == std::string::npos;
    }
}

template<typename TF>
Stats<TF>::Stats(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin):
    master(masterin), grid(gridin), fields(fieldsin),  boundary_cyclic(master, grid)

{
    swstats = inputin.get_item<bool>("stats", "swstats", "", false);

    if (swstats)
    {
        sampletime = inputin.get_item<double>("stats", "sampletime", "");
        masklist   = inputin.get_list<std::string>("stats", "masklist", "", std::vector<std::string>());
        masklist.push_back("default");  // Add the default mask, which calculates the domain mean without sampling.
    }
}

template<typename TF>
Stats<TF>::~Stats()
{
}

template<typename TF>
void Stats<TF>::init(double ifactor)
{
    if (!swstats)
        return;

    auto& gd = grid.get_grid_data();

    boundary_cyclic.init();

    isampletime = static_cast<unsigned long>(ifactor * sampletime);
    statistics_counter = 0;

    // Vectors which hold the amount of grid points sampled on each model level.
    mfield.resize(gd.ncells);
    mfield_bot.resize(gd.ijcells);
}

template<typename TF>
void Stats<TF>::create(int iotime, std::string sim_name)
{
    // Do not create statistics file if stats is disabled.
    if (!swstats)
        return;

    int nerror = 0;
    auto& gd = grid.get_grid_data();

    // Create a NetCDF file for each of the masks.
    for (auto& mask : masks)
    {
        Mask<TF>& m = mask.second;

        if (master.get_mpiid() == 0)
        {
            std::stringstream filename;
            filename << sim_name << "." << m.name << "." << std::setfill('0') << std::setw(7) << iotime << ".nc";

            // Create new NetCDF file, and catch any exceptions locally, to be able
            // to communicate them to the other processes.
            try
            {
                m.data_file = new NcFile(filename.str(), NcFile::newFile);
            }
            catch(NcException& e)
            {
                master.print_error("NetCDF exception: %s\n",e.what());
                ++nerror;
            }
        }

        // Crash on all processes in case the file could not be written.
        master.broadcast(&nerror, 1);
        if (nerror)
            throw 1;

        // Create dimensions.
        if (master.get_mpiid() == 0)
        {
            m.z_dim  = m.data_file->addDim("z" , gd.kmax);
            m.zh_dim = m.data_file->addDim("zh", gd.kmax+1);
            m.t_dim  = m.data_file->addDim("t");

            NcVar z_var;
            NcVar zh_var;

            // Create variables belonging to dimensions.
            m.iter_var = m.data_file->addVar("iter", ncInt, m.t_dim);
            m.iter_var.putAtt("units", "-");
            m.iter_var.putAtt("long_name", "Iteration number");

            m.t_var = m.data_file->addVar("t", ncDouble, m.t_dim);
            m.t_var.putAtt("units", "s");
            m.t_var.putAtt("long_name", "Time");

            z_var = m.data_file->addVar("z", netcdf_fp_type<TF>(), m.z_dim);
            z_var.putAtt("units", "m");
            z_var.putAtt("long_name", "Full level height");

            zh_var = m.data_file->addVar("zh", netcdf_fp_type<TF>(), m.zh_dim);
            zh_var.putAtt("units", "m");
            zh_var.putAtt("long_name", "Half level height");

            // Save the grid variables.
            z_var .putVar(&gd.z [gd.kstart]);
            zh_var.putVar(&gd.zh[gd.kstart]);

            // Synchronize the NetCDF file.
            // BvS: only the last netCDF4-c++ includes the NcFile->sync()
            //      for now use sync() from the netCDF-C library to support older NetCDF4-c++ versions
            //m.data_file->sync();
            nc_sync(m.data_file->getId());
        }

        m.nmask. resize(gd.kcells);
        m.nmaskh.resize(gd.kcells);
        std::cout << mask.first << m.nmask.size() << "\n";

    }

    // For each mask, add the area as a variable.
    add_prof("area" , "Fractional area contained in mask", "-", "z" );
    add_prof("areah", "Fractional area contained in mask", "-", "zh");

}

template<typename TF>
unsigned long Stats<TF>::get_time_limit(unsigned long itime)
{
    // If statistics is disabled, return large (... huge!) value.
    if (!swstats)
        return Constants::ulhuge;

    unsigned long idtlim = isampletime - itime % isampletime;
    return idtlim;
}

template<typename TF>
bool Stats<TF>::do_statistics(unsigned long itime)
{
    // Check if stats are enabled.
    if (!swstats)
        return false;

    // Check if time for execution.
    if (itime % isampletime != 0)
        return false;

    // Return true such that stats are computed.
    return true;
}

template<typename TF>
void Stats<TF>::exec(int iteration, double time, unsigned long itime)
{
    if (!swstats)
        return;

    auto& gd = grid.get_grid_data();

    // check if time for execution
    // BvS: why was this used? This function is only called after a stats->do_statistics(), which already checks the sampletime...
    //if (itime % isampletime != 0)
    //    return;

    // Write message in case stats is triggered
    master.print_message("Saving statistics for time %f\n", time);

    for (auto& mask : masks)
    {
        Mask<TF>& m = mask.second;

        // Put the data into the NetCDF file
        if (master.get_mpiid() == 0)
        {
            const std::vector<size_t> time_index = {static_cast<size_t>(statistics_counter)};

            // Write the time and iteration number
            m.t_var   .putVar(time_index, &time     );
            m.iter_var.putVar(time_index, &iteration);

            const std::vector<size_t> time_height_index = {static_cast<size_t>(statistics_counter), 0};
            std::vector<size_t> time_height_size  = {1, 0};

            //for (Prof_map::const_iterator it=m.profs.begin(); it!=m.profs.end(); ++it)
            for (auto& p : m.profs)
            {
                time_height_size[1] = m.profs[p.first].ncvar.getDim(1).getSize();
                m.profs[p.first].ncvar.putVar(time_height_index, time_height_size, &m.profs[p.first].data.data()[gd.kstart]);
            }

            for (auto& ts: m.tseries)
                m.tseries[ts.first].ncvar.putVar(time_index, &m.tseries[ts.first].data);

            // Synchronize the NetCDF file
            // BvS: only the last netCDF4-c++ includes the NcFile->sync()
            //      for now use sync() from the netCDF-C library to support older NetCDF4-c++ versions
            //m.dataFile->sync();
            nc_sync(m.data_file->getId());
        }
    }
    wmean_set = false;
    // Increment the statistics index
    ++statistics_counter;
}

// Retrieve the user input list of requested masks
template<typename TF>
const std::vector<std::string>& Stats<TF>::get_mask_list()
{
    return masklist;
}

// Add a new mask to the mask map
template<typename TF>
void Stats<TF>::add_mask(const std::string maskname)
{
    auto& gd = grid.get_grid_data();

    masks[maskname].name = maskname;
    masks[maskname].data_file = 0;
    int nmasks = masks.size();
    masks[maskname].flag = (1 << (2 * (nmasks - 1)));
    masks[maskname].flag = (1 << (2 * (nmasks-1) + 1));
}

// Add a new profile to each of the NetCDF files
template<typename TF>
void Stats<TF>::add_prof(std::string name, std::string longname, std::string unit, std::string zloc)
{
    auto& gd = grid.get_grid_data();

    // Add profile to all the NetCDF files
    for (auto& mask : masks)
    {
        Mask<TF>& m = mask.second;

        // Create the NetCDF variable
        if (master.get_mpiid() == 0)
        {
            std::vector<NcDim> dim_vector = {m.t_dim};

            if (zloc == "z")
            {
                dim_vector.push_back(m.z_dim);
                m.profs[name].ncvar = m.data_file->addVar(name, netcdf_fp_type<TF>(), dim_vector);
                //m.profs[name].data = NULL;
            }
            else if (zloc == "zh")
            {
                dim_vector.push_back(m.zh_dim);
                m.profs[name].ncvar = m.data_file->addVar(name.c_str(), netcdf_fp_type<TF>(), dim_vector);
                //m.profs[name].data = NULL;
            }
            m.profs[name].ncvar.putAtt("units", unit.c_str());
            m.profs[name].ncvar.putAtt("long_name", longname.c_str());
            m.profs[name].ncvar.putAtt("_FillValue", netcdf_fp_type<TF>(), netcdf_fp_fillvalue<TF>());

            nc_sync(m.data_file->getId());
        }

        // Resize the vector holding the data at all processes
        m.profs[name].data.resize(gd.kcells);
    }
}

template<typename TF>
void Stats<TF>::add_fixed_prof(std::string name, std::string longname, std::string unit, std::string zloc, TF* restrict prof)
{
    auto& gd = grid.get_grid_data();

    for (auto& mask : masks)
    {
        Mask<TF>& m = mask.second;

        // Create the NetCDF variable
        if (master.get_mpiid() == 0)
        {
           NcVar var;
           if (zloc == "z")
               var = m.data_file->addVar(name, netcdf_fp_type<TF>(), m.z_dim);
           else if (zloc == "zh")
               var = m.data_file->addVar(name, netcdf_fp_type<TF>(), m.zh_dim);
           var.putAtt("units", unit.c_str());
           var.putAtt("long_name", longname.c_str());
           var.putAtt("_FillValue", netcdf_fp_type<TF>(), netcdf_fp_fillvalue<TF>());

           const std::vector<size_t> index = {0};
           if (zloc == "z")
           {
               const std::vector<size_t> size  = {static_cast<size_t>(gd.kmax)};
               var.putVar(index, size, &prof[gd.kstart]);
           }
           else if (zloc == "zh")
           {
               const std::vector<size_t> size  = {static_cast<size_t>(gd.kmax+1)};
               var.putVar(index, size, &prof[gd.kstart]);
           }
       }
   }
}
//
template<typename TF>
void Stats<TF>::add_time_series(const std::string name, const std::string longname, const std::string unit)
{
    // add the series to all files
    for (auto& mask : masks)
    {
        // shortcut
        Mask<TF>& m = mask.second;

        // create the NetCDF variable
        if (master.get_mpiid() == 0)
        {
            m.tseries[name].ncvar = m.data_file->addVar(name.c_str(), netcdf_fp_type<TF>(), m.t_dim);
            m.tseries[name].ncvar.putAtt("units", unit.c_str());
            m.tseries[name].ncvar.putAtt("long_name", longname.c_str());
            m.tseries[name].ncvar.putAtt("_FillValue", netcdf_fp_type<TF>(), netcdf_fp_fillvalue<TF>());
        }

        // Initialize at zero
        m.tseries[name].data = 0.;
    }
}

template<typename TF>
void Stats<TF>::initialize_masks()
{
    auto& gd = grid.get_grid_data();
    for (int n=0; n<gd.ncells; ++n)
        mfield[n] = std::numeric_limits<int>::max();
    for (int n=0; n<gd.ijcells; ++n)
        mfield_bot[n] = std::numeric_limits<int>::max();
}


template<typename TF>
void Stats<TF>::finalize_masks()
{
    auto& gd = grid.get_grid_data();
    const int sloc[] = {0,0,0};
    const int wloc[] = {0,0,1};

    boundary_cyclic.exec(mfield.data());
    boundary_cyclic.exec_2d(mfield_bot.data());

    for (auto& it : masks)
    {
        calc_nmask<TF>(it.second.nmask.data(), it.second.nmaskh.data(), it.second.nmask_bot,
                       mfield.data(), mfield_bot.data(), it.second.flag, it.second.flagh,
                       gd.istart, gd.iend, gd.jstart, gd.jend, 0, gd.kcells,
                       gd.icells, gd.ijcells, gd.kcells);
        master.sum(it.second.nmask.data() , gd.kcells);
        master.sum(it.second.nmaskh.data(), gd.kcells);
        it.second.nmask_bot = it.second.nmaskh[gd.kstart];
        calc_area(it.second.profs["area" ].data.data(), sloc, it.second.nmask .data(), gd.kstart, gd.kend, gd.itot*gd.jtot);
        calc_area(it.second.profs["areah"].data.data(), wloc, it.second.nmaskh.data(), gd.kstart, gd.kend, gd.itot*gd.jtot);
    }
}

template<typename TF>
void Stats<TF>::set_mask_thres(std::string mask_name, Field3d<TF>& fld, Field3d<TF>& fldh, TF threshold, Stats_mask_type mode)
{
    auto& gd = grid.get_grid_data();
    unsigned int flag, flagh;
    bool found_mask = false;

    for(auto& it : masks)
    {
        if(it.second.name == mask_name)
        {
            found_mask = true;
            flag = it.second.flag;
            flagh = it.second.flagh;
        }
    }
    if(!found_mask)
        throw std::runtime_error("Invalid mask name in set_mask_thres()");

    if (mode == Stats_mask_type::Plus)
        calc_mask_thres<TF, Stats_mask_type::Plus>(mfield.data(), mfield_bot.data(), flag, flagh,
            fld.fld.data(), fldh.fld.data(), fldh.fld_bot.data(), threshold,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.ijcells);
    else if (mode == Stats_mask_type::Min)
        calc_mask_thres<TF, Stats_mask_type::Min>(mfield.data(), mfield_bot.data(), flag, flagh,
            fld.fld.data(), fldh.fld.data(), fldh.fld_bot.data(), threshold,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.ijcells);
    else
        throw std::runtime_error("Invalid mask type in set_mask_thres()");
}

template<typename TF>
void Stats<TF>::set_mask_thres_pert(std::string mask_name, Field3d<TF>& fld, Field3d<TF>& fldh, TF threshold, Stats_mask_type mode)
{
    auto& gd = grid.get_grid_data();
    unsigned int flag, flagh;
    bool found_mask = false;

    for(auto& it : masks)
    {
        if(it.second.name == mask_name)
        {
            found_mask = true;
            flag  = it.second.flag;
            flagh = it.second.flagh;
        }
    }
    if(!found_mask)
        throw std::runtime_error("Invalid mask name in set_mask_thres()");

    if (mode == Stats_mask_type::Plus)
        calc_mask_thres_pert<TF, Stats_mask_type::Plus>(mfield.data(), mfield_bot.data(), flag, flagh,
            fld.fld.data(), fld.fld_mean.data(), fldh.fld.data(), fldh.fld_mean.data(), fldh.fld_bot.data(), threshold,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.ijcells);
    else if (mode == Stats_mask_type::Min)
        calc_mask_thres_pert<TF, Stats_mask_type::Min>(mfield.data(), mfield_bot.data(), flag, flagh,
            fld.fld.data(), fld.fld_mean.data(), fldh.fld.data(), fldh.fld_mean.data(), fldh.fld_bot.data(), threshold,
            gd.istart, gd.jstart, gd.kstart,
            gd.iend,   gd.jend,   gd.kend,
            gd.icells, gd.ijcells);
    else
        throw std::runtime_error("Invalid mask type in set_mask_thres_pert()");
}


template<typename TF>
void Stats<TF>::set_prof(const std::string varname, const std::vector<TF> prof)
{
    for (auto& it : masks)
    {
        it.second.profs.at(varname).data = prof;
    }
}

template<typename TF>
void Stats<TF>::calc_stats(const std::string varname, const Field3d<TF>& fld, const int* loc, const TF offset, const TF threshold, std::vector<std::string> operations)
{
    auto& gd = grid.get_grid_data();


    unsigned int flag;
    std::string name;

    sanatize_operations_vector(operations);

    // Process mean first
    if (std::find(operations.begin(), operations.end(), "mean") != operations.end())
    {
        for (auto& m : masks)
        {
            if(loc[2]==0)
                flag = m.second.flag;
            else
                flag = m.second.flagh;

            calc_mean(m.second.profs.at(varname).data.data(), fld.fld.data(), offset, mfield.data(), flag, m.second.nmask.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
            master.sum(m.second.profs.at(varname).data.data(), gd.kcells);
        }
        if(varname == "w")
            wmean_set = true;
        operations.erase(it);
    }

    //Loop over all other operations.
    for(auto& it : operations)
    {
        name = varname+it;
        if(has_only_digits(it))
        {
            int power = std::stoi(it);
            for (auto& m : masks)
            {

                if(loc[2]==0)
                    flag = m.second.flag;
                else
                    flag = m.second.flagh;
                calc_moment(m.second.profs.at(name).data.data(), fld.fld.data(), m.second.profs.at(varname).data.data(), offset, mfield.data(), flag, m.second.nmask.data(),
                        power, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);

                master.sum(m.second.profs.at(name).data.data(), gd.kcells);
            }

        }
        else if (it == "w")
        {
            if(!wmean_set)
                throw std::runtime_error("W mean not calculated in stat - needed for flux");
        }
        else if (it == "diff")
        {

        }
        else if (it == "flux")
        {
/*            for (auto& m : masks)
            {
                add_fluxes(m.second.profs.at(name).data.data(), m.second.profs.at(varname+"w").data.data(), m.second.profs.at(varname+"diff").data.data(),
                        netcdf_fp_fillvalue<TF>(), gd.kstart, gd.kend);
                master.sum(m.second.profs.at(name).data.data(), gd.kcells);
            }
*/
        }
        else if (it == "grad")
        {

        }
        else if (it == "path")
        {

        }
        else if (it == "cover")
        {

        }
        else if (it == "frac")
        {

        }
        else
        {
            std::cout << varname << " " <<it<< "\n";
            throw std::runtime_error("Invalid operations in stat.");
        }
    }
}

template<typename TF>
void Stats<TF>::sanatize_operations_vector(std::vector<std::string> operations)
{
    // Sanatize the operations vector:
    //find instances that need a mean ({2,3,4,5}); if so, add it to the vector if necessary
    for(auto it : operations)
    {
        if(has_only_digits(it))
        {
            operations.push_back("mean");
        }
        if(it == "flux")
        {
            operations.push_back("diff");
            operations.push_back("w");
        }
    }
    // Check for duplicates
    std::sort( operations.begin(), operations.end() );
    operations.erase( unique( operations.begin(), operations.end() ), operations.end() );
    // Make sure that flux goes at the end
    for(auto& it : operations)
    {
        if(it == "flux" )
        {
            std::swap(it, operations.back());
            break;
        }
    }
}
template class Stats<double>;
template class Stats<float>;
