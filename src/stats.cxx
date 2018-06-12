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

    // Sets all the mask values to one (non-masked field)
    template<typename TF>
    void calc_mask(TF* restrict mask_full, TF* restrict mask_half, TF* restrict mask_bottom,
                   int* restrict nmask_full, int* restrict nmask_half, int& nmask_bottom,
                   const int itot, const int jtot, const int kcells, const int ijcells, const int ncells)
    {
        const int ijtot = itot*jtot;
        nmask_bottom = ijtot;

        for (int n=0; n<ncells; ++n)
        {
            mask_full[n] = 1.;
            mask_half[n] = 1.;
        }

        for (int n=0; n<ijcells; ++n)
            mask_bottom[n] = 1.;

        for (int k=0; k<kcells; ++k)
        {
            nmask_full[k] = ijtot;
            nmask_half[k] = ijtot;
        }
    }
    // Sets all the mask values to one (non-masked field)
    template<typename TF>
    void calc_mask_true(TF* restrict mask_full, TF* restrict mask_half, TF* restrict mask_bottom, const int ijcells, const int ncells)
    {
        for (int n=0; n<ncells; ++n)
        {
            mask_full[n] = 1.;
            mask_half[n] = 1.;
        }

        for (int n=0; n<ijcells; ++n)
            mask_bottom[n] = 1.;
    }

    template<typename TF>
    void calc_mask_thres(TF* const restrict mask, TF* const restrict maskh, TF* const restrict maskbot,
                         const TF* const restrict fld,const TF* const restrict fldh, const TF* const restrict fldbot, const TF threshold, const bool mode,
                         const int istart, const int jstart, const int kstart,
                         const int iend,   const int jend,   const int kend,
                         const int icells, const int ijcells)
    {
        if (mode)
        {
            for (int k=kstart; k<kend; k++)
                for (int j=jstart; j<jend; j++)
                    #pragma ivdep
                    for (int i=istart; i<iend; i++)
                    {
                        const int ijk = i + j*icells + k*ijcells;
                        mask[ijk] *= fld[ijk] > threshold;
                        maskh[ijk] *= fldh[ijk] > threshold;
                    }

            // Set the mask for surface projected quantities
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ij  = i + j*icells;
                    const int ijk = i + j*icells + kstart*ijcells;
                    maskbot[ijk] *= fldbot[ijk] > threshold;
                }
        }
        else
        {
            for (int k=kstart; k<kend; k++)
                for (int j=jstart; j<jend; j++)
                    #pragma ivdep
                    for (int i=istart; i<iend; i++)
                    {
                        const int ijk = i + j*icells + k*ijcells;
                        mask[ijk] *= fld[ijk] <= threshold;
                        maskh[ijk] *= fldh[ijk] <= threshold;
                    }

            // Set the mask for surface projected quantities
            // In this case: velocity at surface, so zero
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ij  = i + j*icells;
                    const int ijk = i + j*icells + kstart*ijcells;
                    maskbot[ijk] *= fldbot[ijk] <= threshold;
                }
        }
    }
    template<typename TF>
    void calc_mask_thres_pert(TF* const restrict mask, TF* const restrict maskh, TF* const restrict maskbot,
                     const TF* const restrict fld,const TF* const restrict fld_mean,const TF* const restrict fldh,const TF* const restrict fldh_mean, const TF* const restrict fldbot, const TF threshold, const bool mode,
                     const int istart, const int jstart, const int kstart,
                     const int iend,   const int jend,   const int kend,
                     const int icells, const int ijcells)
    {
        if (mode)
        {
            for (int k=kstart; k<kend; k++)
                for (int j=jstart; j<jend; j++)
                    #pragma ivdep
                    for (int i=istart; i<iend; i++)
                    {
                        const int ijk = i + j*icells + k*ijcells;
                        mask[ijk] *= (fld[ijk]-fld_mean[k]) > threshold;
                        maskh[ijk] *= (fldh[ijk]-fldh_mean[k]) > threshold;
                    }

            // Set the mask for surface projected quantities
            // In this case: velocity at surface, so zero
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ij  = i + j*icells;
                    const int ijk = i + j*icells + kstart*ijcells;
                    maskbot[ijk] *= (fldbot[ijk]-fld_mean[kstart]) > threshold;
                }
        }
        else
        {
            for (int k=kstart; k<kend; k++)
                for (int j=jstart; j<jend; j++)
                    #pragma ivdep
                    for (int i=istart; i<iend; i++)
                    {
                        const int ijk = i + j*icells + k*ijcells;
                        mask[ijk] *= (fld[ijk]-fld_mean[k]) <= threshold;
                        maskh[ijk] *= (fldh[ijk]-fldh_mean[k]) <= threshold;
                    }

            // Set the mask for surface projected quantities
            // In this case: velocity at surface, so zero
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ij  = i + j*icells;
                    const int ijk = i + j*icells + kstart*ijcells;
                    maskbot[ijk] *= (fldbot[ijk]-fld_mean[kstart]) <= threshold;
                }
        }
    }

    // Sets all the mask values to one (non-masked field)
    template<typename TF>
    void calc_nmask(TF* restrict mask_full, TF* restrict mask_half, TF* restrict mask_bottom,
                   int* restrict nmask_full, int* restrict nmask_half, int& nmask_bottom,
                   const int itot, const int jtot, const int ktot, const int ijcells)
    {
        for (int k=0; k<ktot; ++k)
        {
            nmask_full[k] = 0;
            nmask_half[k] = 0;

            for (int j=0; j<jtot; ++j)
                for (int i=0; i<itot; ++i)
                {
                    const int ijk  = i + j*itot + k*ijcells;
                    nmask_full[k]+=mask_full[ijk];
                    nmask_half[k]+=mask_half[ijk];
                }
        }

        nmask_bottom = 0;
        for (int j=0; j<jtot; ++j)
            for (int i=0; i<itot; ++i)
            {
                const int ijk  = i + j*itot;
                nmask_bottom+=mask_bottom[ijk];
            }
    }
}

template<typename TF>
Stats<TF>::Stats(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin):
    master(masterin), grid(gridin), fields(fieldsin)
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

    isampletime = static_cast<unsigned long>(ifactor * sampletime);
    statistics_counter = 0;

    // Vectors which hold the amount of grid points sampled on each model level.
    nmask. resize(gd.kcells);
    nmaskh.resize(gd.kcells);
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
    masks[maskname].name = maskname;
    masks[maskname].data_file = 0;
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
void Stats<TF>::get_mask(Field3d<TF>& mask_full, Field3d<TF>& mask_half)
{
    auto& gd = grid.get_grid_data();

    calc_mask<TF>(mask_full.fld.data(), mask_half.fld.data(), mask_half.fld_bot.data(),
                  nmask.data(), nmaskh.data(), nmaskbot,
                  gd.itot, gd.jtot, gd.kcells, gd.ijcells, gd.ncells);
}

template<typename TF>
void Stats<TF>::get_nmask(Field3d<TF>& mask_full, Field3d<TF>& mask_half)
{
    auto& gd = grid.get_grid_data();
    calc_nmask<TF>(mask_full.fld.data(), mask_half.fld.data(), mask_half.fld_bot.data(),
                  nmask.data(), nmaskh.data(), nmaskbot,
                  gd.itot, gd.jtot, gd.kcells, gd.ijcells);
}
template<typename TF>
void Stats<TF>::set_mask_true(Field3d<TF>& mask_full, Field3d<TF>& mask_half)
{
    auto& gd = grid.get_grid_data();
    calc_mask_true<TF>(mask_full.fld.data(), mask_half.fld.data(), mask_half.fld_bot.data(), gd.ijcells, gd.ncells);
}

template<typename TF>
void Stats<TF>::set_mask_thres(Field3d<TF>& mask_full, Field3d<TF>& mask_half, Field3d<TF>& fld, Field3d<TF>& fldh, TF threshold, Stats_mask_type mode)
{
    auto& gd = grid.get_grid_data();
    calc_mask_thres<TF>(mask_full.fld.data(), mask_half.fld.data(), mask_half.fld_bot.data(),
                  fld.fld.data(), fldh.fld.data(), fldh.fld_bot.data(), threshold, mode==Stats_mask_type::Plus,
                  gd.istart, gd.jstart, gd.kstart,
                  gd.iend,   gd.jend,   gd.kend,
                  gd.icells, gd.ijcells);
}

template<typename TF>
void Stats<TF>::set_mask_thres_pert(Field3d<TF>& mask_full, Field3d<TF>& mask_half, Field3d<TF>& fld, Field3d<TF>& fldh, TF threshold, Stats_mask_type mode)
{
    auto& gd = grid.get_grid_data();
    calc_mask_thres_pert<TF>(mask_full.fld.data(), mask_half.fld.data(), mask_half.fld_bot.data(),
                  fld.fld.data(), fld.fld_mean.data(), fldh.fld.data(), fldh.fld_mean.data(), fldh.fld_bot.data(), threshold, mode==Stats_mask_type::Plus,
                  gd.istart, gd.jstart, gd.kstart,
                  gd.iend,   gd.jend,   gd.kend,
                  gd.icells, gd.ijcells);
}


//// COMPUTATIONAL KERNELS BELOW
//void Stats::calc_mask(double* restrict mask, double* restrict maskh, double* restrict maskbot,
//                      int* restrict nmask, int* restrict nmaskh, int* restrict nmaskbot)
//{
//    int ijtot = grid->itot*grid->jtot;
//
//    // set all the mask values to 1
//    for (int n=0; n<grid->ncells; ++n)
//        mask[n] = 1.;
//
//    for (int n=0; n<grid->ncells; ++n)
//        maskh[n] = 1.;
//
//    for (int n=0; n<grid->ijcells; ++n)
//        maskbot[n] = 1.;
//
//    for (int k=0; k<grid->kcells; ++k)
//    {
//        nmask [k] = ijtot;
//        nmaskh[k] = ijtot;
//    }
//    *nmaskbot = ijtot;
//}


template<typename TF>
void Stats<TF>::calc_area(TF* restrict area, const int loc[3], int* restrict nmask)
{
    auto& gd = grid.get_grid_data();
    const int ijtot = gd.itot*gd.jtot;

    for (int k=gd.kstart; k<gd.kend+loc[2]; k++)
    {
        if (nmask[k] > nthres)
            area[k] = static_cast<TF>(nmask[k]) / static_cast<TF>(ijtot);
        else
            area[k] = 0.;
    }
}


template<typename TF>
void Stats<TF>::calc_mean(TF* const restrict prof, const TF* const restrict data,
                          const TF offset, const TF* const restrict mask, const int * const restrict nmask)
{
    auto& gd = grid.get_grid_data();

    for (int k=gd.kstart; k<gd.kend+1; k++)
    {
        prof[k] = 0.;
        for (int j=gd.jstart; j<gd.jend; j++)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; i++)
            {
                const int ijk  = i + j*gd.icells + k*gd.ijcells;
                prof[k] += mask[ijk]*(data[ijk] + offset);
            }
    }

    master.sum(prof, gd.kcells);

    for (int k=gd.kstart; k<gd.kend+1; k++)
    {
        if (nmask[k] > nthres)
            prof[k] /= static_cast<TF>(nmask[k]);
        else
            prof[k] = netcdf_fp_fillvalue<TF>();
    }
}

template<typename TF>
void Stats<TF>::calc_mean_2d(TF& mean, const TF* const restrict data,
                             const TF offset,
                             const TF* const restrict mask, const int nmask)
{
    auto& gd = grid.get_grid_data();
    const int jj = gd.icells;

    if (nmask > nthres)
    {
        mean = 0.;
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ij = i + j*jj;
                mean += mask[ij]*(data[ij] + offset);
            }
        master.sum(&mean,1);
        mean /= static_cast<TF>(nmask);
    }
    else
        mean = netcdf_fp_fillvalue<TF>();
}

template<typename TF>
void Stats<TF>::calc_max_2d(TF& max, const TF* const restrict data,
                            const TF offset, const TF* const restrict mask, const int nmask)
{
    auto& gd = grid.get_grid_data();
    const int jj = gd.icells;

    if (nmask > nthres)
    {
        max = -dbig;
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ij = i + j*jj;

                if (mask[ij] == 1)
                    max = std::max(max, (data[ij] + offset));
            }
        master.max(&max, 1);
    }
    else
        max = netcdf_fp_fillvalue<TF>();
}

//
//void Stats::calc_sorted_prof(double* restrict data, double* restrict bin, double* restrict prof)
//{
//    const int jj = grid->icells;
//    const int kk = grid->ijcells;
//    const int kstart = grid->kstart;
//    const int kend = grid->kend;
//
//    double minval =  Constants::dhuge;
//    double maxval = -Constants::dhuge;
//
//    // first, get min and max
//    for (int k=grid->kstart; k<grid->kend; ++k)
//        for (int j=grid->jstart; j<grid->jend; j++)
//#pragma ivdep
//            for (int i=grid->istart; i<grid->iend; i++)
//            {
//                const int ijk = i + j*jj + k*kk;
//                if (data[ijk] < minval)
//                    minval = data[ijk];
//                if (data[ijk] > maxval)
//                    maxval = data[ijk];
//            }
//
//    master->min(&minval, 1);
//    master->max(&maxval, 1);
//
//    // make sure that the max ends up in the last bin (introduce 1E-9 error)
//    maxval *= (1.+Constants::dsmall);
//
//    const double range = maxval-minval;
//
//    // In case the field is entirely uniform, dbin becomes zero. In that case we set the profile to the minval.
//    if (range < 1.e-16)
//    {
//        for (int k=grid->kstart; k<grid->kend; ++k)
//            prof[k] = minval;
//    }
//    else
//    {
//        // create bins, equal to the number of grid cells per proc
//        // make sure that bins is not larger than the memory of one 3d field
//        const int bins = grid->nmax;
//
//        // calculate bin width, subtract one to make the minimum and maximum
//        // are in the middle of the bin range and add half a bin size on both sides
//        // |----x----|----x----|----x----|
//        const double dbin = range / (double)(bins-1);
//
//        minval -= 0.5*dbin;
//        maxval += 0.5*dbin;
//
//        // set the bin array to zero
//        for (int n=0; n<bins; ++n)
//            bin[n] = 0;
//
//        // calculate the division factor of one equivalent height unit
//        // (the total volume saved is itot*jtot*zsize)
//        const double nslice = (double)(grid->itot*grid->jtot);
//
//        // check in which bin each value falls and increment the bin count
//        for (int k=grid->kstart; k<grid->kend; ++k)
//        {
//            const double dzslice = grid->dz[k] / nslice;
//            for (int j=grid->jstart; j<grid->jend; ++j)
//                // do not add a ivdep pragma here, because multiple instances could write the same bin[index]
//                for (int i=grid->istart; i<grid->iend; ++i)
//                {
//                    const int ijk = i + j*jj + k*kk;
//                    const int index = (int)((data[ijk] - minval) / dbin);
//                    bin[index] += dzslice;
//                }
//        }
//
//        // get the bin count
//        master->sum(bin, bins);
//
//        // set the starting values of the loop
//        int index = 0;
//        double zbin = 0.5*bin[index];
//        double profval = minval + 0.5*dbin;
//
//        for (int k=grid->kstart; k<grid->kend; ++k)
//        {
//            // Integrate the profile up to the bin count.
//            // Escape the while loop when the integrated profile
//            // exceeds the next grid point.
//            while (zbin < grid->z[k])
//            {
//                zbin += 0.5*(bin[index]+bin[index+1]);
//                profval += dbin;
//                ++index;
//            }
//
//            // In case the first bin is larger than the grid spacing, which can happen
//            // in the inital phase of an MPI run, make sure that no out-of-bounds reads
//            // happen.
//            if (index == 0)
//                prof[k] = profval;
//            else
//            {
//                const double dzfrac = (zbin-grid->z[k]) / (0.5*(bin[index-1]+bin[index]));
//                prof[k] = profval - dzfrac*dbin;
//            }
//        }
//    }
//
//    // now calculate the ghost cells
//    // \TODO this might not be accurate enough, extrapolate properly
//    double profbot = minval;
//    double proftop = maxval;
//
//    if (grid->swspatialorder == "2")
//    {
//        prof[kstart-1] = 2.*profbot - prof[kstart];
//        prof[kend]     = 2.*proftop - prof[kend-1];
//    }
//    else if (grid->swspatialorder == "4")
//    {
//        prof[kstart-1] = (8./3.)*profbot - 2.*prof[kstart] + (1./3.)*prof[kstart+1];
//        prof[kstart-2] = 8.*profbot      - 9.*prof[kstart] + 2.*prof[kstart+1];
//        prof[kend]     = (8./3.)*proftop - 2.*prof[kend-1] + (1./3.)*prof[kend-2];
//        prof[kend+1]   = 8.*proftop      - 9.*prof[kend-1] + 2.*prof[kend-2];
//    }
//}
//
//// \TODO the count function assumes that the variable to count is at the mask location
//void Stats::calc_count(double* restrict data, double* restrict prof, double threshold,
//                       double* restrict mask, int* restrict nmask)
//{
//    const int jj = grid->icells;
//    const int kk = grid->ijcells;
//
//    for (int k=0; k<grid->kcells; ++k)
//    {
//        prof[k] = 0.;
//        for (int j=grid->jstart; j<grid->jend; ++j)
//#pragma ivdep
//            for (int i=grid->istart; i<grid->iend; ++i)
//            {
//                const int ijk = i + j*jj + k*kk;
//                if (data[ijk] > threshold)
//                    prof[k] += mask[ijk]*1.;
//            }
//    }
//
//    master->sum(prof, grid->kcells);
//
//    for (int k=0; k<grid->kcells; k++)
//    {
//        if (nmask[k] > nthres)
//            prof[k] /= (double)(nmask[k]);
//        else
//            prof[k] = NC_FILL_DOUBLE;
//    }
//}

template<typename TF>
void Stats<TF>::calc_moment(TF* restrict data, TF* restrict fld_mean, TF* restrict prof, TF power,
                            TF* restrict mask, int* restrict nmask)
{
    auto& gd = grid.get_grid_data();

    for (int k=gd.kstart; k<gd.kend+1; ++k)
    {
        prof[k] = 0.;
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*gd.icells + k*gd.ijcells;
                prof[k] += mask[ijk]*std::pow(data[ijk]-fld_mean[k], power);
            }
    }

    master.sum(prof, gd.kcells);

    for (int k=gd.kstart; k<gd.kend+1; k++)
    {
        if (nmask[k] > nthres)
            prof[k] /= static_cast<TF>(nmask[k]);
        else
            prof[k] = netcdf_fp_fillvalue<TF>();
    }
}

template<typename TF>
void Stats<TF>::calc_flux_2nd(TF* restrict data, TF* restrict fld_mean, TF* restrict w, TF* restrict wmean,
                              TF* restrict prof, TF* restrict tmp1, const int loc[3],
                              TF* restrict mask, int* restrict nmask)
{
    auto& gd = grid.get_grid_data();
    const int kk = gd.ijcells;

    // set a pointer to the field that contains w, either interpolated or the original
    TF* restrict calcw = w;

    // define the locations
    const int wloc [3] = {0,0,1};
    const int uwloc[3] = {1,0,1};
    const int vwloc[3] = {0,1,1};

    if (loc[0] == 1)
    {
        grid.interpolate_2nd(tmp1, w, wloc, uwloc);
        calcw = tmp1;
    }
    else if (loc[1] == 1)
    {
        grid.interpolate_2nd(tmp1, w, wloc, vwloc);
        calcw = tmp1;
    }

    for (int k=gd.kstart; k<gd.kend+1; ++k)
    {
        prof[k] = 0.;
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk  = i + j*gd.icells + k*kk;
                prof[k] += mask[ijk]*(0.5*(data[ijk-kk]+data[ijk])-0.5*(fld_mean[k-1]+fld_mean[k]))*(calcw[ijk]-wmean[k]);
            }
    }

    master.sum(prof, gd.kcells);

    for (int k=1; k<gd.kcells; ++k)
    {
        if (nmask[k] > nthres && fld_mean[k-1] != netcdf_fp_fillvalue<TF>() && fld_mean[k] != netcdf_fp_fillvalue<TF>())
            prof[k] /= static_cast<TF>(nmask[k]);
        else
            prof[k] = netcdf_fp_fillvalue<TF>();
    }
}

//void Stats::calc_flux_4th(double* restrict data, double* restrict w, double* restrict prof, double* restrict tmp1, const int loc[3],
//                          double* restrict mask, int* restrict nmask)
//{
//    using namespace Finite_difference::O4;
//
//    const int jj  = 1*grid->icells;
//    const int kk1 = 1*grid->ijcells;
//    const int kk2 = 2*grid->ijcells;
//
//    // set a pointer to the field that contains w, either interpolated or the original
//    double* restrict calcw = w;
//
//    // define the locations
//    const int wloc [3] = {0,0,1};
//    const int uwloc[3] = {1,0,1};
//    const int vwloc[3] = {0,1,1};
//
//    if (loc[0] == 1)
//    {
//        grid->interpolate_4th(tmp1, w, wloc, uwloc);
//        calcw = tmp1;
//    }
//    else if (loc[1] == 1)
//    {
//        grid->interpolate_4th(tmp1, w, wloc, vwloc);
//        calcw = tmp1;
//    }
//
//    for (int k=grid->kstart; k<grid->kend+1; ++k)
//    {
//        prof[k] = 0.;
//        for (int j=grid->jstart; j<grid->jend; ++j)
//#pragma ivdep
//            for (int i=grid->istart; i<grid->iend; ++i)
//            {
//                const int ijk = i + j*jj + k*kk1;
//                prof[k] += mask[ijk]*(ci0*data[ijk-kk2] + ci1*data[ijk-kk1] + ci2*data[ijk] + ci3*data[ijk+kk1])*calcw[ijk];
//            }
//    }
//
//    master->sum(prof, grid->kcells);
//
//    for (int k=1; k<grid->kcells; k++)
//    {
//        if (nmask[k] > nthres)
//            prof[k] /= (double)(nmask[k]);
//        else
//            prof[k] = NC_FILL_DOUBLE;
//    }
//}

template<typename TF>
void Stats<TF>::calc_grad_2nd(TF* restrict data, TF* restrict prof, const TF* restrict dzhi,
                              TF* restrict mask, int* restrict nmask)
{
    auto& gd = grid.get_grid_data();

    for (int k=gd.kstart; k<gd.kend+1; ++k)
    {
        prof[k] = 0.;
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*gd.icells + k*gd.ijcells;
                prof[k] += mask[ijk]*(data[ijk]-data[ijk-gd.ijcells])*dzhi[k];
            }
    }

    master.sum(prof, gd.kcells);

    for (int k=gd.kstart; k<gd.kend+1; k++)
    {
        if (nmask[k] > nthres)
            prof[k] /= static_cast<TF>(nmask[k]);
        else
            prof[k] = netcdf_fp_fillvalue<TF>();
    }
}

//void Stats::calc_grad_4th(double* restrict data, double* restrict prof, double* restrict dzhi4, const int loc[3],
//                          double* restrict mask, int* restrict nmask)
//{
//    using namespace Finite_difference::O4;
//
//    const int jj  = 1*grid->icells;
//    const int kk1 = 1*grid->ijcells;
//    const int kk2 = 2*grid->ijcells;
//
//    for (int k=grid->kstart; k<grid->kend+1; ++k)
//    {
//        prof[k] = 0.;
//        for (int j=grid->jstart; j<grid->jend; ++j)
//#pragma ivdep
//            for (int i=grid->istart; i<grid->iend; ++i)
//            {
//                const int ijk = i + j*jj + k*kk1;
//                prof[k] += mask[ijk]*(cg0*data[ijk-kk2] + cg1*data[ijk-kk1] + cg2*data[ijk] + cg3*data[ijk+kk1])*dzhi4[k];
//            }
//    }
//
//    master->sum(prof, grid->kcells);
//
//    for (int k=1; k<grid->kcells; k++)
//    {
//        if (nmask[k] > nthres)
//            prof[k] /= (double)(nmask[k]);
//        else
//            prof[k] = NC_FILL_DOUBLE;
//    }
//}
//
//void Stats::calc_diff_4th(double* restrict data, double* restrict prof, double* restrict dzhi4, double visc, const int loc[3],
//                          double* restrict mask, int* restrict nmask)
//{
//    using namespace Finite_difference::O4;
//
//    const int jj  = 1*grid->icells;
//    const int kk1 = 1*grid->ijcells;
//    const int kk2 = 2*grid->ijcells;
//
//    for (int k=grid->kstart; k<grid->kend+1; ++k)
//    {
//        prof[k] = 0.;
//        for (int j=grid->jstart; j<grid->jend; ++j)
//#pragma ivdep
//            for (int i=grid->istart; i<grid->iend; ++i)
//            {
//                const int ijk = i + j*jj + k*kk1;
//                prof[k] -= mask[ijk]*visc*(cg0*data[ijk-kk2] + cg1*data[ijk-kk1] + cg2*data[ijk] + cg3*data[ijk+kk1])*dzhi4[k];
//            }
//    }
//
//    master->sum(prof, grid->kcells);
//
//    for (int k=1; k<grid->kcells; k++)
//    {
//        if (nmask[k] > nthres)
//            prof[k] /= (double)(nmask[k]);
//        else
//            prof[k] = NC_FILL_DOUBLE;
//    }
//}

template<typename TF>
void Stats<TF>::calc_diff_2nd(TF* restrict data, TF* restrict prof, const TF* restrict dzhi, TF visc, const int loc[3],
                              TF* restrict mask, int* restrict nmask)
{
    auto& gd = grid.get_grid_data();
    const int jj = gd.icells;
    const int kk = gd.ijcells;

    for (int k=gd.kstart; k<gd.kend+1; ++k)
    {
        prof[k] = 0.;
        for (int j=gd.jstart; j<gd.jend; ++j)
            #pragma ivdep
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                prof[k] -= mask[ijk]*visc*(data[ijk] - data[ijk-kk])*dzhi[k];
            }
    }

    master.sum(prof, gd.kcells);

    for (int k=gd.kstart; k<gd.kend+1; k++)
    {
        if (nmask[k] > nthres)
            prof[k] /= static_cast<TF>(nmask[k]);
        else
            prof[k] = netcdf_fp_fillvalue<TF>();
    }
}

template<typename TF>
void Stats<TF>::calc_diff_2nd(
        TF* restrict data, TF* restrict w, TF* restrict evisc,
        TF* restrict prof, const TF* restrict dzhi,
        TF* restrict fluxbot, TF* restrict fluxtop, const TF tPr, const int loc[3],
        TF* restrict mask, int* restrict nmask)
{
    auto& gd = grid.get_grid_data();
    const int ii = 1;
    const int jj = gd.icells;
    const int kk = gd.ijcells;
    const int kstart = gd.kstart;
    const int kend = gd.kend;

    const double dxi = 1./gd.dx;
    const double dyi = 1./gd.dy;

    // bottom boundary
    prof[kstart] = 0.;
    for (int j=gd.jstart; j<gd.jend; ++j)
        #pragma ivdep
        for (int i=gd.istart; i<gd.iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            prof[kstart] += mask[ijk]*fluxbot[ij];
        }

    // calculate the interior
    if (loc[0] == 1)
    {
        for (int k=gd.kstart+1; k<gd.kend; ++k)
        {
            prof[k] = 0.;
            for (int j=gd.jstart; j<gd.jend; ++j)
                #pragma ivdep
                for (int i=gd.istart; i<gd.iend; ++i)
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
        for (int k=gd.kstart+1; k<gd.kend; ++k)
        {
            prof[k] = 0.;
            for (int j=gd.jstart; j<gd.jend; ++j)
                #pragma ivdep
                for (int i=gd.istart; i<gd.iend; ++i)
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
        for (int k=gd.kstart+1; k<gd.kend; ++k)
        {
            prof[k] = 0.;
            for (int j=gd.jstart; j<gd.jend; ++j)
                #pragma ivdep
                for (int i=gd.istart; i<gd.iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    const double eviscs = 0.5*(evisc[ijk-kk]+evisc[ijk])/tPr;
                    prof[k] += -mask[ijk]*eviscs*(data[ijk]-data[ijk-kk])*dzhi[k];
                }
        }
    }

    // top boundary
    prof[kend] = 0.;
    for (int j=gd.jstart; j<gd.jend; ++j)
        #pragma ivdep
        for (int i=gd.istart; i<gd.iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kend*kk;
            prof[kend] += mask[ijk]*fluxtop[ij];
        }

    master.sum(prof, gd.kcells);

    for (int k=1; k<gd.kcells; k++)
    {
        if (nmask[k] > nthres)
            prof[k] /= (double)(nmask[k]);
        else
            prof[k] = NC_FILL_DOUBLE;
    }
}

template<typename TF>
void Stats<TF>::add_fluxes(TF* restrict flux, TF* restrict turb, TF* restrict diff)
{
    auto& gd = grid.get_grid_data();

    for (int k=gd.kstart; k<gd.kend+1; ++k)
    {
        if (turb[k] == netcdf_fp_fillvalue<TF>() || diff[k] == netcdf_fp_fillvalue<TF>())
            flux[k] = netcdf_fp_fillvalue<TF>();
        else
            flux[k] = turb[k] + diff[k];
    }
}

///**
// * This function calculates the total domain integrated path of variable data over maskbot
// */
//void Stats::calc_path(double* restrict data, double* restrict maskbot, int* restrict nmaskbot, double* restrict path)
//{
//    const int jj = grid->icells;
//    const int kk = grid->ijcells;
//    const int kstart = grid->kstart;
//
//    *path = 0.;
//
//    if (*nmaskbot > nthres)
//    {
//        // Integrate liquid water
//        for (int j=grid->jstart; j<grid->jend; j++)
//            for (int i=grid->istart; i<grid->iend; i++)
//            {
//                const int ij = i + j*jj;
//                if (maskbot[ij] == 1)
//                    for (int k=kstart; k<grid->kend; k++)
//                    {
//                        const int ijk = i + j*jj + k*kk;
//                        *path += fields->rhoref[k] * data[ijk] * grid->dz[k];
//                    }
//            }
//        *path /= (double)*nmaskbot;
//        master->sum(path, 1);
//    }
//    else
//        *path = NC_FILL_DOUBLE;
//}
//
///**
// * This function calculates the vertical projected cover of variable data over maskbot
// */
//void Stats::calc_cover(double* restrict data, double* restrict maskbot, int* restrict nmaskbot, double* restrict cover, double threshold)
//{
//    const int jj = grid->icells;
//    const int kk = grid->ijcells;
//    const int kstart = grid->kstart;
//
//    *cover = 0.;
//
//    if (*nmaskbot > nthres)
//    {
//        // Per column, check if cloud present
//        for (int j=grid->jstart; j<grid->jend; j++)
//            for (int i=grid->istart; i<grid->iend; i++)
//            {
//                const int ij = i + j*jj;
//                if (maskbot[ij] == 1)
//                    for (int k=kstart; k<grid->kend; k++)
//                    {
//                        const int ijk = i + j*jj + k*kk;
//                        if (data[ijk]>threshold)
//                        {
//                            *cover += 1.;
//                            break;
//                        }
//                    }
//            }
//        *cover /= (double)*nmaskbot;
//        master->sum(cover,1);
//    }
//    else
//        *cover = NC_FILL_DOUBLE;
//}

template class Stats<double>;
template class Stats<float>;
