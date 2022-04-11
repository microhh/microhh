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
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <vector>
#include <utility>
#include <netcdf.h>

#include "master.h"
#include "grid.h"
#include "soil_grid.h"
#include "fields.h"
#include "stats.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "timeloop.h"
#include "advec.h"
#include "diff.h"
#include "netcdf_interface.h"

namespace
{
    using namespace Constants;

    // Help function(s) to switch between the different NetCDF data types
    template<typename TF> TF netcdf_fp_fillvalue();
    template<> double netcdf_fp_fillvalue<double>() { return NC_FILL_DOUBLE; }
    template<> float  netcdf_fp_fillvalue<float>()  { return NC_FILL_FLOAT; }

    template<typename TF, Stats_mask_type mode>
    TF is_false(const TF value, const TF threshold)
    {
        if (mode == Stats_mask_type::Plus)
            return (value <= threshold);
        else if (mode == Stats_mask_type::Min)
            return (value > threshold);
    }

    template<typename TF>
    TF in_mask(const unsigned int mask, const unsigned int flag)
    {
        return static_cast<TF>( (mask & flag) != 0 );
    }

    template<typename TF>
    void set_flag(unsigned int& flag, const int*& restrict nmask, const Mask<TF>& m, const int loc)
    {
        if (loc == 0)
        {
            flag = m.flag;
            nmask = m.nmask.data();
        }
        else
        {
            flag = m.flagh;
            nmask = m.nmaskh.data();
        }
    }

    template<typename TF>
    void set_fillvalue_prof(TF* const restrict data, const int* const restrict nmask, const int kstart, const int kcells)
    {
        for (int k=0; k<kcells; ++k)
        {
            if (nmask[k] == 0)
                data[k] = netcdf_fp_fillvalue<TF>();
        }
    }

    template<typename TF, Stats_mask_type mode>
    void calc_mask_thres(
            unsigned int* const restrict mfield, unsigned int* const restrict mfield_bot,
            const unsigned int flag, const unsigned int flagh,
            const TF* const restrict fld, const TF* const restrict fldh,
            const TF* const restrict fld_bot, const TF threshold,
            const int istart, const int jstart, const int kstart,
            const int iend, const int jend, const int kend,
            const int icells, const int ijcells, const int kcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    mfield[ijk] -= (mfield[ijk] & flag ) * is_false<TF, mode>(fld [ijk], threshold);
                    mfield[ijk] -= (mfield[ijk] & flagh) * is_false<TF, mode>(fldh[ijk], threshold);
                }

        // Set the top value for the flux level.
        #pragma omp parallel for
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*icells + kend*ijcells;
                mfield[ijk] -= (mfield[ijk] & flagh) * is_false<TF, mode>(fldh[ijk], threshold);
            }

        // Set the ghost cells equal to the first model level.
        #pragma omp parallel for
        for (int k=0; k<kstart; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    const int ijk_ref = i + j*icells + kstart*ijcells;
                    mfield[ijk] -= (mfield[ijk] & flag ) * is_false<TF, mode>(fld [ijk_ref], threshold);
                    mfield[ijk] -= (mfield[ijk] & flagh) * is_false<TF, mode>(fldh[ijk_ref], threshold);
                }

        // Set the ghost cells for the full level equal to kend-1.
        #pragma omp parallel for
        for (int k=kend; k<kcells; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    const int ijk_ref = i + j*icells + (kend-1)*ijcells;
                    mfield[ijk] -= (mfield[ijk] & flag) * is_false<TF, mode>(fld [ijk_ref], threshold);
                }

        // Set the ghost cells for the flux level equal to kend.
        #pragma omp parallel for
        for (int k=kend+1; k<kcells; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    const int ijk_ref = i + j*icells + kend*ijcells;
                    mfield[ijk] -= (mfield[ijk] & flagh) * is_false<TF, mode>(fldh[ijk_ref], threshold);
                }

        // Set the mask for surface projected quantities
        #pragma omp parallel for
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij = i + j*icells;
                mfield_bot[ij] -= (mfield_bot[ij] & flag) * is_false<TF, mode>(fld_bot[ij], threshold);
            }
    }

    /*
    template<typename TF, Stats_mask_type mode>
    void calc_mask_thres_pert(
            unsigned int* const restrict mfield, unsigned int* const restrict mfield_bot,
            const unsigned int flag, const unsigned int flagh,
            const TF* const restrict fld, const TF* const restrict fld_mean, const TF* const restrict fldh,
            const TF* const restrict fldh_mean, const TF* const restrict fld_bot, const TF threshold,
            const int istart, const int jstart, const int kstart,
            const int iend,   const int jend,   const int kend,
            const int icells, const int ijcells)
    {

        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    mfield[ijk] -=  (mfield[ijk] & flag ) * is_false<TF, mode>(fld [ijk]-fld_mean [k], threshold);
                    mfield[ijk] -=  (mfield[ijk] & flagh) * is_false<TF, mode>(fldh[ijk]-fldh_mean[k], threshold);
                }

        // Set the mask for surface projected quantities
        #pragma omp parallel for
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                mfield_bot[ij] -= (mfield_bot[ij] & flag) * is_false<TF, mode>(fld_bot[ij]-fld_mean[kstart], threshold);
            }
    }
    */

    template<typename TF>
    void calc_area(
            TF* const restrict area, const int loc[3], const int* const restrict nmask,
            const int kstart, const int kend, const int ijtot)
    {
        for (int k=kstart; k<kend+loc[2]; ++k)
        {
            if (nmask[k])
                area[k] = static_cast<TF>(nmask[k]) / static_cast<TF>(ijtot);
            else
                area[k] = 0.;
        }
    }

    // Calculate the number of points contained in the mask.
    template<typename TF>
    void calc_nmask(
            int* restrict nmask_full, int* restrict nmask_half, int& nmask_bottom,
            const unsigned int* const mfield, const unsigned int* const mfield_bot,
            const unsigned int flag, const unsigned int flagh,
            const int istart, const int iend, const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells, const int kcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
        {
            nmask_full[k] = 0;
            nmask_half[k] = 0;

            // #pragma omp parallel for reduction (+:nmask_full[k], nmask_half[k]) collapse(2)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    nmask_full[k] += in_mask<int>(mfield[ijk], flag );
                    nmask_half[k] += in_mask<int>(mfield[ijk], flagh);
                }
        }

        nmask_bottom     = 0;
        // nmask_half[kend] = 0;
        // #pragma omp parallel for reduction (+:nmask_bottom) collapse(2)
        // #pragma omp parallel for reduction (+:nmask_bottom, nmask_half[kend]) collapse(2)
        for (int j=jstart; j<jend; ++j)
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                const int ijk = i + j*icells + kend*ijcells;
                nmask_bottom += in_mask<int>(mfield_bot[ij], flag);
                // nmask_half[kend] += in_mask<int>(mfield[ijk], flagh);
            }
    }

    template<typename TF>
    void calc_mean(
            TF* const restrict prof, const TF* const restrict fld,
            const unsigned int* const mask, const unsigned int flag, const int* const nmask,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend+1; ++k)
        {
            if (nmask[k])
            {
                double tmp = 0.;
                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk  = i + j*icells + k*ijcells;
                        tmp += in_mask<double>(mask[ijk], flag) * fld[ijk];
                    }

                prof[k] = tmp / nmask[k];
            }
        }
    }

    template<typename TF>
    void calc_mean_2d(
            TF& out, const TF* const restrict fld,
            const int istart, const int iend, const int jstart, const int jend,
            const int icells, const int itot, const int jtot)
    {
        double tmp = 0.;

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                tmp += fld[ij];
            }

        out = tmp / (itot*jtot);
    }


    template<typename TF>
    void calc_mean_projected_mask(
            TF* const restrict prof, const TF* const restrict fld,
            const unsigned int* const mask_bot, const int nmask_bot,
            const int flag,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        double tmp;

        if (nmask_bot == 0)
            return;

        for (int k=kstart; k<kend; ++k)
        {
            tmp = 0;
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i  + j*icells;
                    const int ijk = ij + k*ijcells;

                    tmp += in_mask<double>(mask_bot[ij], flag) * fld[ijk];
                }
            prof[k] = tmp / nmask_bot;
        }
    }


    template<typename TF>
    void calc_moment(
            TF* const restrict prof, const TF* const restrict fld, const TF* const restrict fld_mean, const TF offset,
            const unsigned int* const mask, const unsigned int flag, const int* const nmask, const int power,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend+1; ++k)
        {
            if (nmask[k])
            {
                double tmp = 0.;
                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk  = i + j*icells + k*ijcells;
                        tmp += in_mask<double>(mask[ijk], flag)*std::pow(fld[ijk] - fld_mean[k] + offset, power);
                    }

                prof[k] = tmp / nmask[k];
            }
        }
    }

    template<typename TF>
    void calc_cov(
            TF* const restrict prof, const TF* const restrict fld1, const TF* const restrict fld1_mean,
            const TF offset1, const int pow1,
            const TF* const restrict fld2, const TF* const restrict fld2_mean, const TF offset2, const int pow2,
            const unsigned int* const mask, const unsigned int flag, const int* const nmask,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend+1; ++k)
        {
            if (nmask[k])
            {
                double tmp = 0.;
                if ((fld1_mean[k] != netcdf_fp_fillvalue<TF>()) && (fld2_mean[k] != netcdf_fp_fillvalue<TF>()))
                {
                    for (int j=jstart; j<jend; ++j)
                        #pragma ivdep
                        for (int i=istart; i<iend; ++i)
                        {
                            const int ijk  = i + j*icells + k*ijcells;
                            tmp += in_mask<double>(mask[ijk], flag)
                                * std::pow(fld1[ijk] - fld1_mean[k] + offset1, pow1)
                                * std::pow(fld2[ijk] - fld2_mean[k] + offset2, pow2);
                        }

                    prof[k] = tmp / nmask[k];
                }
            }
        }
    }

    template<typename TF>
    void add_fluxes(
            TF* const restrict flux, const TF* const restrict turb, const TF* const restrict diff,
            const int kstart, const int kend)
    {
        for (int k=kstart; k<kend+1; ++k)
            flux[k] = turb[k] + diff[k];
    }


    template<typename TF>
    void calc_frac(
            TF* const restrict prof, const TF* const restrict fld, const TF offset, const TF threshold,
            const unsigned int* const mask, const unsigned int flag, const int* const nmask,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend+1; ++k)
        {
            if (nmask[k])
            {
                double tmp = 0.;
                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk  = i + j*icells + k*ijcells;
                        tmp += in_mask<double>(mask[ijk], flag)*((fld[ijk] + offset) > threshold);
                    }
                prof[k] = tmp / nmask[k];
            }
        }
    }

    template<typename TF>
    std::pair<TF, int> calc_path(
            const TF* const restrict data, const TF* const restrict dz, const TF* const restrict rho,
            const unsigned int* const mask, const unsigned int flag, const int* const nmask,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        int nmask_proj = 0;
        TF path = TF(0.);

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                for (int k=kstart; k<kend; ++k)
                {
                    const int ijk = i + j*jj + k*kk;
                    if (in_mask<bool>(mask[ijk], flag))
                    {
                        ++nmask_proj;
                        break;
                    }
                }

                if (nmask_proj > 0)
                {
                    for (int k=kstart; k<kend; ++k)
                    {
                        const int ijk = i + j*jj + k*kk;
                        if (in_mask<bool>(mask[ijk], flag))
                        {
                            path += data[ijk]*rho[k]*dz[k];
                        }
                    }
                }
            }

        return std::make_pair(path, nmask_proj);
    }

    template<typename TF>
    std::pair<int, int> calc_cover(
            const TF* const restrict fld, const TF offset, const TF threshold,
            const unsigned int* const mask, const unsigned int flag, const int* const nmask,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        int cover = 0.;
        int nmaskcover = 0;

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                int maskincolumn = 0;
                for (int k=kstart; k<kend; ++k)
                {
                    const int ijk  = i + j*icells + k*ijcells;
                    if (in_mask<bool>(mask[ijk], flag))
                    {
                        maskincolumn = 1;
                        if ((fld[ijk] + offset) > threshold)
                        {
                            ++cover;
                            break;
                        }
                    }
                }
                nmaskcover += maskincolumn;
            }

        return std::make_pair(cover, nmaskcover);
    }

    bool has_only_digits(const std::string& s)
    {
        return s.find_first_not_of( "23456789" ) == std::string::npos;
    }
}

template<typename TF>
Stats<TF>::Stats(
        Master& masterin, Grid<TF>& gridin, Soil_grid<TF>& soilgridin,
        Fields<TF>& fieldsin, Advec<TF>& advecin, Diff<TF>& diffin, Input& inputin):
    master(masterin), grid(gridin), soil_grid(soilgridin), fields(fieldsin), advec(advecin), diff(diffin),
    boundary_cyclic(master, grid)

{
    swstats = inputin.get_item<bool>("stats", "swstats", "", false);

    if (swstats)
    {
        sampletime = inputin.get_item<double>("stats", "sampletime", "");
        masklist   = inputin.get_list<std::string>("stats", "masklist", "", std::vector<std::string>());
        masklist.push_back("default"); // Add the default mask, which calculates the domain mean without sampling.

        // Add user XY masks
        std::vector<std::string> xymasklist = inputin.get_list<std::string>("stats", "xymasklist", "", std::vector<std::string>());
        masklist.insert(masklist.end(), xymasklist.begin(), xymasklist.end());

        swtendency = inputin.get_item<bool>("stats", "swtendency", "", false);
        std::vector<std::string> whitelistin = inputin.get_list<std::string>("stats", "whitelist", "", std::vector<std::string>());

        // Anything without an underscore is mean value, so should be on the whitelist
        // std::regex re("^[^_]*$");
        // whitelist.push_back(re);

        for (auto& it : whitelistin)
        {
            std::regex re(it);
            whitelist.push_back(re);
        }

        std::vector<std::string> blacklistin = inputin.get_list<std::string>("stats", "blacklist", "", std::vector<std::string>());

        for (auto& it : blacklistin)
        {
            std::regex re(it);
            blacklist.push_back(re);
        }
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
void Stats<TF>::create(const Timeloop<TF>& timeloop, std::string sim_name)
{
    // Do not create statistics file if stats is disabled.
    if (!swstats)
        return;

    int iotime = timeloop.get_iotime();

    auto& gd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    // Create a NetCDF file for each of the masks.
    for (auto& mask : masks)
    {
        Mask<TF>& m = mask.second;

        std::stringstream filename;
        filename << sim_name << "." << m.name << "." << std::setfill('0') << std::setw(7) << iotime << ".nc";

        // Create new NetCDF file
        m.data_file = std::make_unique<Netcdf_file>(master, filename.str(), Netcdf_mode::Create);

        // Create dimensions.
        m.data_file->add_dimension("z",  gd.kmax);
        m.data_file->add_dimension("zh", gd.kmax+1);
        m.data_file->add_dimension("time");

        if (sgd.is_enabled)
        {
            m.data_file->add_dimension("zs", sgd.kmax);
            m.data_file->add_dimension("zsh", sgd.kmax+1);
        }

        // Create variables belonging to dimensions.
        Netcdf_handle& iter_handle =
            m.data_file->group_exists("default") ? m.data_file->get_group("default") : m.data_file->add_group("default");

        m.iter_var = std::make_unique<Netcdf_variable<int>>(iter_handle.add_variable<int>("iter", {"time"}));
        m.iter_var->add_attribute("units", "-");
        m.iter_var->add_attribute("long_name", "Iteration number");

        m.time_var = std::make_unique<Netcdf_variable<TF>>(m.data_file->template add_variable<TF>("time", {"time"}));
        if (timeloop.has_utc_time())
            m.time_var->add_attribute("units", "seconds since " + timeloop.get_datetime_utc_start_string());
        else
            m.time_var->add_attribute("units", "seconds since start");
        m.time_var->add_attribute("long_name", "Time");

        // Add vertical grid variables (atmosphere)
        Netcdf_variable<TF> z_var = m.data_file->template add_variable<TF>("z", {"z"});
        z_var.add_attribute("units", "m");
        z_var.add_attribute("long_name", "Full level height");

        Netcdf_variable<TF> zh_var = m.data_file->template add_variable<TF>("zh", {"zh"});
        zh_var.add_attribute("units", "m");
        zh_var.add_attribute("long_name", "Half level height");

        // Save the grid variables.
        std::vector<TF> z_nogc (gd.z. begin() + gd.kstart, gd.z. begin() + gd.kend  );
        std::vector<TF> zh_nogc(gd.zh.begin() + gd.kstart, gd.zh.begin() + gd.kend+1);
        z_var .insert( z_nogc, {0});
        zh_var.insert(zh_nogc, {0});

        // Add vertical grid variables (soil)
        if (sgd.is_enabled)
        {
            Netcdf_variable<TF> zs_var = m.data_file->template add_variable<TF>("zs", {"zs"});
            zs_var.add_attribute("units", "m");
            zs_var.add_attribute("long_name", "Full level height soil");

            Netcdf_variable<TF> zsh_var = m.data_file->template add_variable<TF>("zsh", {"zsh"});
            zsh_var.add_attribute("units", "m");
            zsh_var.add_attribute("long_name", "Half level height soil");

            // Save the grid variables.
            std::vector<TF> zs_nogc (sgd.z. begin() + sgd.kstart, sgd.z. begin() + sgd.kend  );
            std::vector<TF> zsh_nogc(sgd.zh.begin() + sgd.kstart, sgd.zh.begin() + sgd.kend+1);
            zs_var .insert( zs_nogc, {0});
            zsh_var.insert(zsh_nogc, {0});
        }

        // Synchronize the NetCDF file.
        m.data_file->sync();

        m.nmask. resize(gd.kcells);
        m.nmaskh.resize(gd.kcells);
    }

    // For each mask, add the area as a variable.
    add_prof("area" , "Fractional area contained in mask", "-", "z" , "default");
    add_prof("areah", "Fractional area contained in mask", "-", "zh", "default");
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
void Stats<TF>::set_tendency(bool in)
{
    doing_tendency = in;
}

template<typename TF>
void Stats<TF>::exec(const int iteration, const double time, const unsigned long itime)
{
    if (!swstats)
        return;

    auto& agd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    // Write message in case stats is triggered.
    master.print_message("Saving statistics for time %f\n", time);

    // Finalize the total tendencies
    if (do_tendency())
    {
        for (auto& var: tendency_order)
        {
            //Subtract the previous tendency from the current one; skipping the last one (total tendency)
            for (auto tend = std::next(var.second.rbegin()); tend != std::prev(var.second.rend()); ++tend)
            {
                auto prev_tend = std::next(tend); // It is a reverse iterator....
                std::string name = var.first + "_" + *tend;
                std::string prev_name = var.first + "_" + *prev_tend;
                for (auto& mask : masks)
                {
                    Mask<TF>& m = mask.second;

                    std::transform (m.profs.at(name).data.begin(), m.profs.at(name).data.end(),
                                    m.profs.at(prev_name).data.begin(), m.profs.at(name).data.begin(), std::minus<TF>());

                    const int* nmask;
                    if (m.profs.at(name).level == Level_type::Full)
                        nmask = m.nmask.data();
                    else
                        nmask = m.nmaskh.data();
                    set_fillvalue_prof(m.profs.at(name).data.data(), nmask, agd.kstart, agd.kcells);
                }
            }

            var.second.clear();
        }
    }

    for (auto& mask : masks)
    {
        Mask<TF>& m = mask.second;

        // Put the data into the NetCDF file.
        const std::vector<int> time_index{statistics_counter};

        // Write the time and iteration number.
        m.time_var->insert(time     , time_index);
        m.iter_var->insert(iteration, time_index);

        const std::vector<int> time_height_index = {statistics_counter, 0};

        for (auto& p : m.profs)
        {
            const int ksize = p.second.ncvar.get_dim_sizes()[1];
            std::vector<int> time_height_size  = {1, ksize};

            std::vector<TF> prof_nogc(
                    p.second.data.begin() + agd.kstart,
                    p.second.data.begin() + agd.kstart + ksize);

            m.profs.at(p.first).ncvar.insert(prof_nogc, time_height_index, time_height_size);
        }

        for (auto& p : m.soil_profs)
        {
            const int ksize = p.second.ncvar.get_dim_sizes()[1];
            std::vector<int> time_height_size  = {1, ksize};

            std::vector<TF> prof_nogc(
                    p.second.data.begin() + sgd.kstart,
                    p.second.data.begin() + sgd.kstart + ksize);

            m.soil_profs.at(p.first).ncvar.insert(prof_nogc, time_height_index, time_height_size);
        }

        for (auto& ts : m.tseries)
            m.tseries.at(ts.first).ncvar.insert(m.tseries.at(ts.first).data, time_index);

        // Synchronize the NetCDF file.
        m.data_file->sync();
    }

    wmean_set = false;

    // Increment the statistics index.
    ++statistics_counter;
}

// Retrieve the user input list of requested masks.
template<typename TF>
const std::vector<std::string>& Stats<TF>::get_mask_list()
{
    return masklist;
}

// Add a new dimension to the stats file.
template<typename TF>
void Stats<TF>::add_dimension(
        const std::string& name, const int size)
{
    // Create a NetCDF file for each of the masks.
    for (auto& mask : masks)
        mask.second.data_file->add_dimension(name, size);
}

// Add a new mask to the mask map.
template<typename TF>
void Stats<TF>::add_mask(const std::string& maskname)
{
    masks.emplace(maskname, Mask<TF>{});
    masks.at(maskname).name = maskname;
    masks.at(maskname).data_file = 0;

    int nmasks = masks.size();
    masks.at(maskname).flag  = (1 << (2 * (nmasks-1)    ));
    masks.at(maskname).flagh = (1 << (2 * (nmasks-1) + 1));
}

// Add a new profile to each of the NetCDF files.
template<typename TF>
void Stats<TF>::add_profs(
        const Field3d<TF>& var, const std::string& zloc, std::vector<std::string> operations, const std::string& group_name)
{
    std::string zloc_alt;
    if (zloc == "z")
        zloc_alt = "zh";
    else if (zloc == "zh")
        zloc_alt = "z";
    else
        throw std::runtime_error(zloc + " is an invalid height in add_profs");

    sanitize_operations_vector(var.name, operations);

    for (auto& it : operations)
    {
        if (it == "mean")
        {
            add_prof(var.name, var.longname, var.unit, zloc, group_name, Stats_whitelist_type::White);
        }
        else if (has_only_digits(it))
        {
            add_prof(var.name + "_" + it, "Moment " + it + " of the " + var.longname,fields.simplify_unit(var.unit, "",std::stoi(it)), zloc, group_name);
        }
        else if (it == "w")
        {
            add_prof(var.name+"_w", "Turbulent flux of the " + var.longname, fields.simplify_unit(var.unit, "m s-1"), zloc_alt, group_name);
        }
        else if (it == "grad")
        {
            add_prof(var.name+"_grad", "Gradient of the " + var.longname, fields.simplify_unit(var.unit, "m-1"), zloc_alt, group_name);
        }
        else if (it == "flux")
        {
            add_prof(var.name+"_flux", "Total flux of the " + var.longname, fields.simplify_unit(var.unit, "m s-1"), zloc_alt, group_name);
        }
        else if (it == "diff")
        {
            add_prof(var.name+"_diff", "Diffusive flux of the " + var.longname, fields.simplify_unit(var.unit, "m s-1"), zloc_alt, group_name);
        }
        else if (it == "frac")
        {
            add_prof(var.name+"_frac", var.longname + " fraction", "-", zloc, group_name);
        }
        else if (it == "path")
        {
            add_time_series(var.name+"_path", var.longname + " path", fields.simplify_unit(var.unit, "kg m-2"), group_name);
        }
        else if (it == "cover")
        {
            add_time_series(var.name+"_cover", var.longname + " cover", "-", group_name);
        }
        else
        {
            throw std::runtime_error(it + "is an invalid operator name to add profs");
        }
    }
}

// Add a new profile to each of the NetCDF files.
template<typename TF>
void Stats<TF>::add_tendency(
        const Field3d<TF>& var, const std::string& zloc,
        const std::string& tend_name, const std::string& tend_longname, const std::string& group_name)
{
    if (swstats && swtendency)
    {
        add_prof(var.name + "_" + tend_name, tend_longname + " " + var.longname, var.unit, zloc, group_name);
        if (tendency_order.find(var.name) == tendency_order.end())
        {
            tendency_order[var.name];
            add_tendency(var, zloc, "total", "Total");
        }
    }
}

template<typename TF>
void Stats<TF>::add_covariance(const Field3d<TF>& var1, const Field3d<TF>& var2, const std::string& zloc, const std::string& group_name)
{
    for (int pow1 = 1; pow1<5; ++pow1)
    {
        for (int pow2 = 1; pow2<5; ++pow2)
        {
            std::string spow1 = std::to_string(pow1);
            std::string spow2 = std::to_string(pow2);

            std::string name = var1.name + "_" + spow1 + "_" + var2.name + "_" + spow2;
            std::string longname = "Covariance of " + var1.name + spow1 + " and " + var2.name + spow2;
            add_prof(name, longname, fields.simplify_unit(var1.unit, var2.unit, pow1, pow2), zloc, group_name, Stats_whitelist_type::Black);
        }
    }

}


template<typename TF>
void Stats<TF>::add_operation(
        std::vector<std::string>& operations, const std::string& varname, const std::string& op)
{
    operations.push_back(op);
    if (is_blacklisted(varname + "_" +  op))
    {
        std::regex re(varname + "_" + op);
        whitelist.push_back(re);
    }

}

template<typename TF>
void Stats<TF>::sanitize_operations_vector(
        const std::string& varname, std::vector<std::string>& operations)
{
    // Sanitize the operations vector:
    // find instances that need a mean ({2,3,4,5}); if so, add it to the vector if necessary
    for (auto it = operations.begin(); it != operations.end(); )
    {
        if (is_blacklisted(varname + "_" + *it))
        {
            it = operations.erase(it);
        }
        else
        {
            ++it;
        }
    }

    auto tmpvec = operations;
    for (auto& it : tmpvec)
    {
        if (it == "flux")
        {
            add_operation(operations, varname, "diff");
            add_operation(operations, varname, "w");
        }
        else if (has_only_digits(it))
        {
            add_operation(operations, varname, "mean");
        }

        else if (it == "w")
        {
            add_operation(operations, varname, "mean");
        }
    }

    // Check for duplicates
    std::sort( operations.begin(), operations.end() );
    operations.erase( std::unique( operations.begin(), operations.end() ), operations.end() );

    // Make sure that mean goes in the front
    for (auto& it : operations)
    {
        if (it == "mean" )
        {
            std::swap(it, operations.front());
            break;
        }
    }
    // Make sure that flux goes at the end
    for (auto& it : operations)
    {
        if (it == "flux")
        {
            std::swap(it, operations.back());
            break;
        }
    }
}

template<typename TF>
bool Stats<TF>::is_blacklisted(const std::string& name, Stats_whitelist_type wltype)
{
    if (wltype == Stats_whitelist_type::White)
        return false;

    for (const auto& it : whitelist)
    {
        if (std::regex_match(name, it))
            return false;
    }

    if (wltype == Stats_whitelist_type::Black)
        return true;

    for (const auto& it : blacklist)
    {
        if (std::regex_match(name, it))
            return true;
    }

    return false;
}

// Add a new profile to each of the NetCDF files
template<typename TF>
void Stats<TF>::add_prof(
        const std::string& name, const std::string& longname,
        const std::string& unit, const std::string& zloc, const std::string& group_name,
        Stats_whitelist_type wltype)
{
    auto& agd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    // Check whether variable is part of whitelist/blacklist and is not the area coverage of the mask.
    if (is_blacklisted(name, wltype))
        return;

    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
        return;

    Level_type level;
    if ((zloc == "z") || (zloc == "zs"))
        level = Level_type::Full;
    else
        level = Level_type::Half;

    // Add profile to all the NetCDF files.
    for (auto& mask : masks)
    {
        Mask<TF>& m = mask.second;

        Netcdf_handle& handle = (group_name == "") ? dynamic_cast<Netcdf_handle&>(*m.data_file) : dynamic_cast<Netcdf_handle&>
            (m.data_file->group_exists(group_name) ? m.data_file->get_group(group_name) : m.data_file->add_group(group_name));

        // Create the NetCDF variable.
        // Create the profile variable and the vector at the appropriate size.

        if ((zloc == "z") || (zloc == "zh"))
        {
            Prof_var<TF> tmp{handle.add_variable<TF>(name, {"time", zloc}), std::vector<TF>(agd.kcells), level};

            m.profs.emplace(std::piecewise_construct, std::forward_as_tuple(name), std::forward_as_tuple(std::move(tmp)));

            m.profs.at(name).ncvar.add_attribute("units", unit);
            m.profs.at(name).ncvar.add_attribute("long_name", longname);
        }
        else
        {
            Prof_var<TF> tmp{handle.add_variable<TF>(name, {"time", zloc}), std::vector<TF>(sgd.kcells), level};

            m.soil_profs.emplace(std::piecewise_construct, std::forward_as_tuple(name), std::forward_as_tuple(std::move(tmp)));

            m.soil_profs.at(name).ncvar.add_attribute("units", unit);
            m.soil_profs.at(name).ncvar.add_attribute("long_name", longname);
        }

        m.data_file->sync();
    }

    if ((zloc == "z") || (zloc == "zh"))
        varlist.push_back(name);
    else
        varlist_soil.push_back(name);
}

template<typename TF>
void Stats<TF>::add_fixed_prof(
        const std::string& name,
        const std::string& longname,
        const std::string& unit,
        const std::string& zloc,
        const std::string& group_name,
        const std::vector<TF>& prof)
{
    auto& gd = grid.get_grid_data();

    for (auto& mask : masks)
    {
        Mask<TF>& m = mask.second;

        Netcdf_handle& handle = (group_name == "") ? dynamic_cast<Netcdf_handle&>(*m.data_file) : dynamic_cast<Netcdf_handle&>
            (m.data_file->group_exists(group_name) ? m.data_file->get_group(group_name) : m.data_file->add_group(group_name));

        // Create the NetCDF variable.
        Netcdf_variable<TF> var = handle.add_variable<TF>(name, {zloc});

        var.add_attribute("units", unit.c_str());
        var.add_attribute("long_name", longname.c_str());
        // var.add_attribute("_FillValue", netcdf_fp_fillvalue<TF>());

        if (zloc == "z")
        {
            std::vector<TF> prof_nogc(prof.begin() + gd.kstart, prof.begin() + gd.kend);
            var.insert(prof_nogc, {0}, {gd.ktot});
        }
        else if (zloc == "zh")
        {
            std::vector<TF> prof_nogc(prof.begin() + gd.kstart, prof.begin() + gd.kend+1);
            var.insert(prof_nogc, {0}, {gd.ktot+1});
        }

        m.data_file->sync();
    }
}

template<typename TF>
void Stats<TF>::add_fixed_prof_raw(
        const std::string& name,
        const std::string& longname,
        const std::string& unit,
        const std::string& dim,
        const std::string& group_name,
        const std::vector<TF>& prof)
{
    for (auto& mask : masks)
    {
        Mask<TF>& m = mask.second;

        Netcdf_handle& handle = (group_name == "") ? dynamic_cast<Netcdf_handle&>(*m.data_file) : dynamic_cast<Netcdf_handle&>
            (m.data_file->group_exists(group_name) ? m.data_file->get_group(group_name) : m.data_file->add_group(group_name));

        // Create the NetCDF variable.
        Netcdf_variable<TF> var = handle.add_variable<TF>(name, {dim});

        var.add_attribute("units", unit.c_str());
        var.add_attribute("long_name", longname.c_str());

        const int size = prof.size();
        var.insert(prof, {0}, {size});

        m.data_file->sync();
    }
}

template<typename TF>
void Stats<TF>::add_time_series(
        const std::string& name, const std::string& longname,
        const std::string& unit, const std::string& group_name, Stats_whitelist_type wltype)
{
    // Check whether variable is part of whitelist/blacklist.
    if (is_blacklisted(name, wltype))
        return;

    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
        return;

    // Add the series to all files.
    for (auto& mask : masks)
    {
        // Shortcut
        Mask<TF>& m = mask.second;

        Netcdf_handle& handle = (group_name == "") ? dynamic_cast<Netcdf_handle&>(*m.data_file) : dynamic_cast<Netcdf_handle&>
            (m.data_file->group_exists(group_name) ? m.data_file->get_group(group_name) : m.data_file->add_group(group_name));

        // Create the NetCDF variable
        Time_series_var<TF> tmp{handle.add_variable<TF>(name, {"time"}), 0.};

        m.tseries.emplace(
                std::piecewise_construct, std::forward_as_tuple(name), std::forward_as_tuple(std::move(tmp)));

        m.tseries.at(name).ncvar.add_attribute("units", unit);
        m.tseries.at(name).ncvar.add_attribute("long_name", longname);
        // m.tseries.at(name).ncvar.add_attribute("_FillValue", netcdf_fp_fillvalue<TF>());
    }

    varlist.push_back(name);

}

template<typename TF>
void Stats<TF>::initialize_masks()
{
    auto& gd = grid.get_grid_data();
    unsigned int flagmax = 0;

    for (auto& it : masks)
        flagmax += it.second.flag + it.second.flagh;

    for (int n=0; n<gd.ncells; ++n)
        mfield[n] = flagmax;

    for (int n=0; n<gd.ijcells; ++n)
        mfield_bot[n] = flagmax;
}

template<typename TF>
void Stats<TF>::finalize_masks()
{
    auto& gd = grid.get_grid_data();

    boundary_cyclic.exec(mfield.data());
    boundary_cyclic.exec_2d(mfield_bot.data());

    for (auto& it : masks)
    {
        // CvH: compute the nmask over the entire depth. Masks need to provide the proper count for
        // the ghost cells in order to be able to calculate mean profile in ghost cells (needed for budgets).
        calc_nmask<TF>(
                it.second.nmask.data(), it.second.nmaskh.data(), it.second.nmask_bot,
                mfield.data(), mfield_bot.data(), it.second.flag, it.second.flagh,
                gd.istart, gd.iend, gd.jstart, gd.jend, 0, gd.kcells,
                gd.icells, gd.ijcells, gd.kcells);

        master.sum(it.second.nmask.data() , gd.kcells);
        master.sum(it.second.nmaskh.data(), gd.kcells);

        it.second.nmask_bot = it.second.nmaskh[gd.kstart];

        auto it1 = std::find(varlist.begin(), varlist.end(), "area");
        if (it1 != varlist.end())
            calc_area(it.second.profs.at("area").data.data(), gd.sloc.data(), it.second.nmask.data(),
                    gd.kstart, gd.kend, gd.itot*gd.jtot);

        it1 = std::find(varlist.begin(), varlist.end(), "areah");
        if (it1 != varlist.end())
            calc_area(it.second.profs.at("areah").data.data(), gd.wloc.data(), it.second.nmaskh.data(),
                    gd.kstart, gd.kend, gd.itot*gd.jtot);
    }
}

template<typename TF>
void Stats<TF>::set_mask_thres(
        std::string mask_name, Field3d<TF>& fld, Field3d<TF>& fldh, TF threshold, Stats_mask_type mode)
{
    auto& gd = grid.get_grid_data();
    unsigned int flag, flagh;
    bool found_mask = false;

    for (auto& it : masks)
    {
        if (it.second.name == mask_name)
        {
            found_mask = true;
            flag = it.second.flag;
            flagh = it.second.flagh;
        }
    }

    if (!found_mask)
        throw std::runtime_error("Invalid mask name in set_mask_thres()");

    if (mode == Stats_mask_type::Plus)
        calc_mask_thres<TF, Stats_mask_type::Plus>(
                mfield.data(), mfield_bot.data(), flag, flagh,
                fld.fld.data(), fldh.fld.data(), fldh.fld_bot.data(), threshold,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells, gd.kcells);
    else if (mode == Stats_mask_type::Min)
        calc_mask_thres<TF, Stats_mask_type::Min>(
                mfield.data(), mfield_bot.data(), flag, flagh,
                fld.fld.data(), fldh.fld.data(), fldh.fld_bot.data(), threshold,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.ijcells, gd.kcells);
    else
        throw std::runtime_error("Invalid mask type in set_mask_thres()");
}

template<typename TF>
void Stats<TF>::set_prof(const std::string& varname, const std::vector<TF>& prof)
{
    auto it = std::find(varlist.begin(), varlist.end(), varname);
    if (it == varlist.end())
        throw std::runtime_error("Set_prof: Variable " + varname + " does not exist");
    else
    {
        for (auto& it : masks)
            it.second.profs.at(varname).data = prof;
    }
}

template<typename TF>
void Stats<TF>::set_time_series(const std::string& varname, const TF val)
{
    auto it = std::find(varlist.begin(), varlist.end(), varname);
    if (it != varlist.end())
    {
        for (auto& it : masks)
            it.second.tseries.at(varname).data = val;
    }
}

template<typename TF>
void Stats<TF>::calc_mask_mean_profile(
        std::vector<TF>& prof,
        const std::pair<const std::string, Mask<TF>>& m,
        const Field3d<TF>& fld)
{
    auto& gd = grid.get_grid_data();

    unsigned int flag;
    const int* nmask;

    set_flag(flag, nmask, m.second, fld.loc[2]);

    // CvH. Do the mean over the entire depth. The calc_mean function always add 1 to the specified
    // kend, so I send kcells-1 as the limit. This is not elegant, yet it works.
    calc_mean(
            prof.data(), fld.fld.data(), mfield.data(), flag, nmask,
            gd.istart, gd.iend, gd.jstart, gd.jend, 0, gd.kcells-1, gd.icells, gd.ijcells);

    master.sum(prof.data(), gd.kcells);
}

template<typename TF>
void Stats<TF>::calc_mask_stats(
        std::pair<const std::string, Mask<TF>>& m,
        const std::string& varname, const Field3d<TF>& fld, const TF offset, const TF threshold)
{
    auto& gd = grid.get_grid_data();

    unsigned int flag;
    const int* nmask;
    std::string name;

    // Calc mean
    if (std::find(varlist.begin(), varlist.end(), varname) != varlist.end())
    {
        set_flag(flag, nmask, m.second, fld.loc[2]);
        calc_mean(m.second.profs.at(varname).data.data(), fld.fld.data(), mfield.data(), flag, nmask,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
        master.sum(m.second.profs.at(varname).data.data(), gd.kcells);

        // Add the offset.
        for (auto& value : m.second.profs.at(varname).data)
            value += offset;

        set_fillvalue_prof(m.second.profs.at(varname).data.data(), nmask, gd.kstart, gd.kcells);
    }

    // Calc moments
    for (int power=2; power<=4; power++)
    {
        name = varname + "_" + std::to_string(power);
        if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
        {
            set_flag(flag, nmask, m.second, fld.loc[2]);
            calc_moment(
                    m.second.profs.at(name).data.data(), fld.fld.data(),
                    m.second.profs.at(varname).data.data(), offset, mfield.data(), flag, nmask,
                    power, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

            master.sum(m.second.profs.at(name).data.data(), gd.kcells);
            set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);
        }
    }

    // Calc Resolved Flux
    name = varname + "_w";
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        auto advec_flux = fields.get_tmp();
        advec.get_advec_flux(*advec_flux, fld);

        set_flag(flag, nmask, m.second, !fld.loc[2]);
        calc_mean(
                m.second.profs.at(name).data.data(), advec_flux->fld.data(), mfield.data(), flag, nmask,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        master.sum(m.second.profs.at(name).data.data(), gd.kcells);
        set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);

        fields.release_tmp(advec_flux);
    }

    // Calc Diffusive Flux
    name = varname + "_diff";
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        auto diff_flux = fields.get_tmp();
        diff.diff_flux(*diff_flux, fld);

        set_flag(flag, nmask, m.second, !fld.loc[2]);
        calc_mean(
                m.second.profs.at(name).data.data(), diff_flux->fld.data(), mfield.data(), flag, nmask,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        master.sum(m.second.profs.at(name).data.data(), gd.kcells);
        set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);

        fields.release_tmp(diff_flux);
    }

    // Calc Total Flux
    name = varname + "_flux";
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        // No sum is required in this routine as values all.
        set_flag(flag, nmask, m.second, !fld.loc[2]);
        add_fluxes(
                m.second.profs.at(name).data.data(), m.second.profs.at(varname+"_w").data.data(), m.second.profs.at(varname+"_diff").data.data(),
                gd.kstart, gd.kend);
        set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);
    }

    // Calc Gradient
    name = varname + "_grad";
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        set_flag(flag, nmask, m.second, !fld.loc[2]);

        if (grid.get_spatial_order() == Grid_order::Second)
        {
            calc_grad_2nd(
                    m.second.profs.at(name).data.data(), fld.fld.data(), gd.dzhi.data(), mfield.data(), flag, nmask,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
        }
        else if (grid.get_spatial_order() == Grid_order::Fourth)
        {
            calc_grad_4th(
                    m.second.profs.at(name).data.data(), fld.fld.data(), gd.dzhi4.data(), mfield.data(), flag, nmask,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
        }

        master.sum(m.second.profs.at(name).data.data(), gd.kcells);
        set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);
    }

    // Calc Integrated Path
    name = varname + "_path";
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        set_flag(flag, nmask, m.second, fld.loc[2]);

        std::pair<TF, int> path = calc_path(
                fld.fld.data(), gd.dz.data(), fields.rhoref.data(),
                mfield.data(), flag, nmask,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        master.sum(&path.first, 1);
        master.sum(&path.second, 1);

        m.second.tseries.at(name).data = path.first / path.second;
    }

    // Calc Cover
    name = varname + "_cover";
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        set_flag(flag, nmask, m.second, fld.loc[2]);

        // Function returns number of poinst covered (cover.first) and number of points in mask (cover.second).
        std::pair<int, int> cover = calc_cover(
                fld.fld.data(), offset, threshold, mfield.data(), flag, nmask,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        master.sum(&cover.first, 1);
        master.sum(&cover.second, 1);

        // Only assign if number of points in mask is positive.
        m.second.tseries.at(name).data = (cover.second > 0) ? TF(cover.first)/TF(cover.second) : 0.;
    }

    // Calc Fraction
    name = varname + "_frac";
    auto it1 = std::find(varlist.begin(), varlist.end(), name);
    if (it1 != varlist.end())
    {
        set_flag(flag, nmask, m.second, fld.loc[2]);

        calc_frac(
                m.second.profs.at(name).data.data(), fld.fld.data(),
                offset, threshold, mfield.data(), flag, nmask,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        master.sum(m.second.profs.at(name).data.data(), gd.kcells);
        set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);
    }
}

template<typename TF>
void Stats<TF>::calc_stats(
        const std::string& varname, const Field3d<TF>& fld, const TF offset, const TF threshold)
{
    auto& gd = grid.get_grid_data();

    unsigned int flag;
    const int* nmask;
    std::string name;

    // Calc mean of atmospheric variables
    if (std::find(varlist.begin(), varlist.end(), varname) != varlist.end())
    {
        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, fld.loc[2]);
            calc_mean(m.second.profs.at(varname).data.data(), fld.fld.data(), mfield.data(), flag, nmask,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
            master.sum(m.second.profs.at(varname).data.data(), gd.kcells);

            // Add the offset.
            for (auto& value : m.second.profs.at(varname).data)
                value += offset;

            set_fillvalue_prof(m.second.profs.at(varname).data.data(), nmask, gd.kstart, gd.kcells);
        }
    }

    // Calc moments
    for (int power=2; power<=4; power++)
    {
        name = varname + "_" + std::to_string(power);
        if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
        {
            for (auto& m : masks)
            {
                set_flag(flag, nmask, m.second, fld.loc[2]);
                calc_moment(
                        m.second.profs.at(name).data.data(), fld.fld.data(),
                        m.second.profs.at(varname).data.data(), offset, mfield.data(), flag, nmask,
                        power, gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                        gd.icells, gd.ijcells);

                master.sum(m.second.profs.at(name).data.data(), gd.kcells);
                set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);
            }
        }
    }

    // Calc Resolved Flux
    name = varname + "_w";
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        auto advec_flux = fields.get_tmp();
        advec.get_advec_flux(*advec_flux, fld);

        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, !fld.loc[2]);
            calc_mean(
                    m.second.profs.at(name).data.data(), advec_flux->fld.data(), mfield.data(), flag, nmask,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

            master.sum(m.second.profs.at(name).data.data(), gd.kcells);
            set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);
        }
        fields.release_tmp(advec_flux);
    }

    // Calc Diffusive Flux
    name = varname + "_diff";
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        auto diff_flux = fields.get_tmp();
        diff.diff_flux(*diff_flux, fld);

        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, !fld.loc[2]);
            calc_mean(
                    m.second.profs.at(name).data.data(), diff_flux->fld.data(), mfield.data(), flag, nmask,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

            master.sum(m.second.profs.at(name).data.data(), gd.kcells);
            set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);
        }

        fields.release_tmp(diff_flux);
    }

    // Calc Total Flux
    name = varname + "_flux";
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        for (auto& m : masks)
        {
            // No sum is required in this routine as values all.
            set_flag(flag, nmask, m.second, !fld.loc[2]);
            add_fluxes(
                    m.second.profs.at(name).data.data(), m.second.profs.at(varname+"_w").data.data(),
                    m.second.profs.at(varname+"_diff").data.data(),
                    gd.kstart, gd.kend);
            set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);
        }
    }

    // Calc Gradient
    name = varname + "_grad";
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, !fld.loc[2]);

            if (grid.get_spatial_order() == Grid_order::Second)
            {
                calc_grad_2nd(
                        m.second.profs.at(name).data.data(), fld.fld.data(), gd.dzhi.data(), mfield.data(), flag, nmask,
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                        gd.icells, gd.ijcells);
            }
            else if (grid.get_spatial_order() == Grid_order::Fourth)
            {
                calc_grad_4th(
                        m.second.profs.at(name).data.data(), fld.fld.data(), gd.dzhi4.data(), mfield.data(), flag, nmask,
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                        gd.icells, gd.ijcells);
            }

            master.sum(m.second.profs.at(name).data.data(), gd.kcells);
            set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);
        }
    }

    // Calc Integrated Path
    name = varname + "_path";
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, fld.loc[2]);

            std::pair<TF, int> path = calc_path(
                    fld.fld.data(), gd.dz.data(), fields.rhoref.data(),
                    mfield.data(), flag, nmask,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

            master.sum(&path.first, 1);
            master.sum(&path.second, 1);

            m.second.tseries.at(name).data = path.first / path.second;
        }
    }

    // Calc Cover
    name = varname + "_cover";
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, fld.loc[2]);

            // Function returns number of poinst covered (cover.first) and number of points in mask (cover.second).
            std::pair<int, int> cover = calc_cover(
                    fld.fld.data(), offset, threshold, mfield.data(), flag, nmask,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

            master.sum(&cover.first, 1);
            master.sum(&cover.second, 1);

            // Only assign if number of points in mask is positive.
            m.second.tseries.at(name).data = (cover.second > 0) ? TF(cover.first)/TF(cover.second) : 0.;
        }
    }

    // Calc Fraction
    name = varname + "_frac";
    auto it1 = std::find(varlist.begin(), varlist.end(), name);
    if (it1 != varlist.end())
    {
        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, fld.loc[2]);

            calc_frac(
                    m.second.profs.at(name).data.data(), fld.fld.data(), offset, threshold, mfield.data(), flag, nmask,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);

            master.sum(m.second.profs.at(name).data.data(), gd.kcells);
            set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);
        }
    }
}

template<typename TF>
void Stats<TF>::calc_tend(Field3d<TF>& fld, const std::string& tend_name)
{
    if (!doing_tendency)
        return;

    auto& gd = grid.get_grid_data();
    unsigned int flag;
    const int* nmask;

    std::string name = fld.name + "_" + tend_name;
    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        tendency_order.at(fld.name).push_back(tend_name);
        #ifdef USECUDA
        fields.backward_field_device_3d(fld.fld.data(), fld.fld_g);
        #endif
        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, fld.loc[2]);
            calc_mean(m.second.profs.at(name).data.data(), fld.fld.data(), mfield.data(), flag, nmask,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend, gd.icells, gd.ijcells);
            master.sum(m.second.profs.at(name).data.data(), gd.kcells);

            set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);
        }
    }
}

template<typename TF>
void Stats<TF>::calc_stats_2d(
        const std::string& varname, const std::vector<TF>& fld, const TF offset)
{
    auto& gd = grid.get_grid_data();

    if (std::find(varlist.begin(), varlist.end(), varname) != varlist.end())
    {
        for (auto& m : masks)
        {
            calc_mean_2d(m.second.tseries.at(varname).data, fld.data(),
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.icells, gd.itot, gd.jtot);
            master.sum(&m.second.tseries.at(varname).data, 1);
            m.second.tseries.at(varname).data += offset;
        }
    }
}

template<typename TF>
void Stats<TF>::calc_stats_soil(
        const std::string varname, const std::vector<TF>& fld, const TF offset)
{
    /*
       Calculate soil statistics, using the surface mask
       projected on the entire soil column.
    */
    auto& agd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    if (std::find(varlist_soil.begin(), varlist_soil.end(), varname) != varlist_soil.end())
    {
        for (auto& m : masks)
        {
            if (m.second.nmask_bot > 0)
            {
                calc_mean_projected_mask(
                        m.second.soil_profs.at(varname).data.data(), fld.data(),
                        mfield_bot.data(), m.second.nmask_bot,
                        m.second.flag,
                        agd.istart, agd.iend,
                        agd.jstart, agd.jend,
                        sgd.kstart, sgd.kend,
                        agd.icells, agd.ijcells);

                // Add offset
                for (auto& value : m.second.soil_profs.at(varname).data)
                    value += offset;

                master.sum(m.second.soil_profs.at(varname).data.data(), sgd.kmax);
            }
            else
            {
                for (auto& value : m.second.soil_profs.at(varname).data)
                    value = netcdf_fp_fillvalue<TF>();
            }
        }
    }
}

template<typename TF>
void Stats<TF>::calc_covariance(
        const std::string& varname1, const Field3d<TF>& fld1, const TF offset1, const TF threshold1, const int power1,
        const std::string& varname2, const Field3d<TF>& fld2, const TF offset2, const TF threshold2, const int power2)
{
    auto& gd = grid.get_grid_data();

    std::string name = varname1 + "_" + std::to_string(power1) + "_" + varname2 + "_" + std::to_string(power2);
    unsigned int flag;

    int* nmask;

    if (std::find(varlist.begin(), varlist.end(), name) != varlist.end())
    {
        if (fld1.loc == fld2.loc)
        {
            TF fld1_mean[gd.kcells];
            for (auto& m : masks)
            {
                if (fld2.loc[2] == 0)
                {
                    flag = m.second.flag;
                    for (int k = gd.kstart; k<gd.kend+1; ++k)
                    {
                        fld1_mean[k] = m.second.profs.at(varname1).data[k];
                    }
                    nmask = m.second.nmask.data();
                }
                else
                {
                    flag = m.second.flagh;
                    for (int k = gd.kstart; k<gd.kend+1; ++k)
                    {
                        if (fld1_mean[k-1] != netcdf_fp_fillvalue<TF>() && fld1_mean[k] != netcdf_fp_fillvalue<TF>())
                            fld1_mean[k] = 0.5*(m.second.profs.at(varname1).data[k]+m.second.profs.at(varname1).data[k-1]);
                        else
                            fld1_mean[k] = netcdf_fp_fillvalue<TF>();

                    }
                    nmask = m.second.nmaskh.data();
                }

                calc_cov(
                        m.second.profs.at(name).data.data(), fld1.fld.data(), fld1_mean, offset1, power1,
                        fld2.fld.data(), m.second.profs.at(varname2).data.data(), offset2, power2,
                        mfield.data(), flag, nmask,
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                        gd.icells, gd.ijcells);

                master.sum(m.second.profs.at(name).data.data(), gd.kcells);
                set_fillvalue_prof(m.second.profs.at(name).data.data(), nmask, gd.kstart, gd.kcells);
            }
        }
        else
        {
            auto tmp = fields.get_tmp();
            if (grid.get_spatial_order() == Grid_order::Second)
                grid.interpolate_2nd(tmp->fld.data(), fld1.fld.data(), fld1.loc.data(), fld2.loc.data());
            else if (grid.get_spatial_order() == Grid_order::Fourth)
                grid.interpolate_4th(tmp->fld.data(), fld1.fld.data(), fld1.loc.data(), fld2.loc.data());

            for (auto& m : masks)
            {
                if (fld2.loc[2] == 0)
                {
                    flag = m.second.flag;
                    nmask = m.second.nmask.data();
                }
                else
                {
                    flag = m.second.flagh;
                    nmask = m.second.nmaskh.data();
                }

                calc_cov(
                        m.second.profs.at(name).data.data(), tmp->fld.data(),
                        m.second.profs.at(varname1).data.data(), offset1, power1,
                        fld2.fld.data(), m.second.profs.at(varname2).data.data(), offset2, power2,
                        mfield.data(), flag, nmask,
                        gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                        gd.icells, gd.ijcells);

                master.sum(m.second.profs.at(name).data.data(), gd.kcells);
            }

            fields.release_tmp(tmp);
        }
    }
}

template<typename TF>
void Stats<TF>::calc_flux_2nd(
        TF* const restrict prof, const TF* const restrict data, const TF* const restrict fld_mean, const TF offset,
        TF* const restrict w, const TF* const restrict wmean,
        TF* restrict tmp1, const int loc[3], const unsigned int* const mask, const unsigned int flag, const int* const nmask,
        const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
        const int icells, const int ijcells)
{
    auto& gd = grid.get_grid_data();

    // set a pointer to the field that contains w, either interpolated or the original
    TF* restrict calcw = w;

    if (loc[0] == 1)
    {
        grid.interpolate_2nd(tmp1, w, gd.wloc.data(), gd.uwloc.data());
        calcw = tmp1;
    }
    else if (loc[1] == 1)
    {
        grid.interpolate_2nd(tmp1, w, gd.wloc.data(), gd.vwloc.data());
        calcw = tmp1;
    }

    #pragma omp parallel for
    for (int k=kstart; k<kend+1; ++k)
    {
        if (nmask[k])
        {
            prof[k] = 0.; // This assignment is crucial to avoid unitialized values in case if below is false.

            // Check whether mean is contained in the mask as well. It cannot do this check at the top point, which exceptionally could lead to problems.
            if ((fld_mean[k-1] != netcdf_fp_fillvalue<TF>()) && (fld_mean[k] != netcdf_fp_fillvalue<TF>()) && (k != kend))
            {
                double tmp = 0;

                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk = i + j*icells + k*ijcells;
                        const TF data_prime = 0.5*(data[ijk-ijcells]+data[ijk]) + offset - 0.5*(fld_mean[k-1]+fld_mean[k]);
                        const TF w_prime = calcw[ijk] - wmean[k];
                        tmp += in_mask<double>(mask[ijk], flag) * data_prime * w_prime;
                    }

                prof[k] = tmp / nmask[k];
            }
        }
    }
}

template<typename TF>
void Stats<TF>::calc_flux_4th(
        TF* const restrict prof, const TF* const restrict data, const TF* const restrict fld_mean, TF* const restrict w, const TF* const restrict wmean,
        TF* restrict tmp1, const int loc[3], const unsigned int* const mask, const unsigned int flag, const int* const nmask,
        const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
        const int icells, const int ijcells)
{
    using namespace Finite_difference::O4;

    auto& gd = grid.get_grid_data();

    const int kk1 = 1*ijcells;
    const int kk2 = 2*ijcells;

    // set a pointer to the field that contains w, either interpolated or the original
    TF* restrict calcw = w;

    if (loc[0] == 1)
    {
        grid.interpolate_4th(tmp1, w, gd.wloc.data(), gd.uwloc.data());
        calcw = tmp1;
    }
    else if (loc[1] == 1)
    {
        grid.interpolate_4th(tmp1, w, gd.wloc.data(), gd.vwloc.data());
        calcw = tmp1;
    }

    #pragma omp parallel for
    for (int k=kstart; k<kend+1; ++k)
    {
        if (nmask[k] && fld_mean[k-1] != netcdf_fp_fillvalue<TF>() && fld_mean[k] != netcdf_fp_fillvalue<TF>())
        {
            double tmp = 0.;
            if ((fld_mean[k-1] != netcdf_fp_fillvalue<TF>()) && (fld_mean[k] != netcdf_fp_fillvalue<TF>()))
            {
                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk  = i + j*icells + k*ijcells;
                        tmp = in_mask<double>(mask[ijk], flag)*(ci0<double>*data[ijk-kk2] + ci1<double>*data[ijk-kk1] + ci2<double>*data[ijk] + ci3<double>*data[ijk+kk1])*calcw[ijk];
                    }
                prof[k] = tmp / nmask[k];
            }
        }
    }
}

template<typename TF>
void Stats<TF>::calc_grad_2nd(
        TF* const restrict prof, const TF* const restrict data, const TF* const restrict dzhi,
        const unsigned int* const mask, const unsigned int flag, const int* const nmask,
        const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
        const int icells, const int ijcells)
{
    #pragma omp parallel for
    for (int k=kstart; k<kend+1; ++k)
    {
        if (nmask[k])
        {
            double tmp = 0.;
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    tmp += in_mask<double>(mask[ijk], flag)*(data[ijk]-data[ijk-ijcells])*dzhi[k];
                }

            prof[k] = tmp / nmask[k];
        }
    }
}

template<typename TF>
void Stats<TF>::calc_grad_4th(
        TF* const restrict prof, const TF* const restrict data, const TF* const restrict dzhi4,
        const unsigned int* const mask, const unsigned int flag, const int* const nmask,
        const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
        const int icells, const int ijcells)
{
    using namespace Finite_difference::O4;

    const int jj  = 1*icells;
    const int kk1 = 1*ijcells;
    const int kk2 = 2*ijcells;

    #pragma omp parallel for
    for (int k=kstart; k<kend+1; ++k)
    {
        if (nmask[k])
        {
            double tmp = 0.;
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk1;
                    tmp += in_mask<double>(mask[ijk], flag)*(cg0<double>*data[ijk-kk2] + cg1<double>*data[ijk-kk1] + cg2<double>*data[ijk] + cg3<double>*data[ijk+kk1])*dzhi4[k];
                }

            prof[k] = tmp / nmask[k];
        }
    }
}

template class Stats<double>;
template class Stats<float>;
