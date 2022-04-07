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

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <boost/algorithm/string.hpp>

#include "master.h"
#include "grid.h"
#include "fields.h"
#include "field3d.h"
#include "soil_field3d.h"
#include "input.h"
#include "netcdf_interface.h"
#include "defines.h"
#include "finite_difference.h"
#include "stats.h"
#include "column.h"
#include "cross.h"
#include "dump.h"
#include "diff.h"

namespace
{
    enum class Mask_type {Wplus, Wmin};

    template<typename TF, Mask_type mode>
    void calc_mask_w(TF* const restrict mask, TF* const restrict maskh, TF* const restrict maskbot,
                     int* const restrict nmask, int* const restrict nmaskh, const TF* const restrict w,
                     const int istart, const int jstart, const int kstart,
                     const int iend,   const int jend,   const int kend,
                     const int icells, const int ijcells)
    {
        int ntmp = 0;

        // Calculate mask at full levels
        for (int k=kstart; k<kend; k++)
        {
            nmask[k] = 0;
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;

                    if (mode == Mask_type::Wplus)
                        ntmp = (w[ijk] + w[ijk+ijcells]) >  0.;
                    else
                        ntmp = (w[ijk] + w[ijk+ijcells]) <= 0.;

                    nmask[k] += ntmp;
                    mask[ijk] = static_cast<TF>(ntmp);
                }
        }

        // Calculate mask at half levels
        for (int k=kstart; k<kend+1; k++)
        {
            nmaskh[k] = 0;
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ijk = i + j*icells + k*ijcells;

                    if (mode == Mask_type::Wplus)
                        ntmp = w[ijk] >  0.;
                    else
                        ntmp = w[ijk] <= 0.;

                    nmaskh[k] += ntmp;
                    maskh[ijk] = static_cast<TF>(ntmp);
                }
        }

        // Set the mask for surface projected quantities
        // In this case: velocity at surface, so zero
        for (int j=jstart; j<jend; j++)
            #pragma ivdep
            for (int i=istart; i<iend; i++)
            {
                const int ij  = i + j*icells;
                const int ijk = i + j*icells + kstart*ijcells;

                maskbot[ij] = maskh[ijk];
            }
    }

    template<typename TF>
    void set_xy_mask(
            TF* const restrict mask,
            TF* const restrict maskh,
            const TF* const restrict xymask,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        for (int k=kstart; k<kend; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ij = i + j*icells;
                    const int ijk = ij + k*ijcells;
                    mask[ijk] = xymask[ij] > TF(0.5) ? TF(1) : TF(0);
                }

        for (int k=kstart; k<kend+1; k++)
            for (int j=jstart; j<jend; j++)
                #pragma ivdep
                for (int i=istart; i<iend; i++)
                {
                    const int ij = i + j*icells;
                    const int ijk = ij + k*ijcells;
                    maskh[ijk] = xymask[ij] > TF(0.5) ? TF(1) : TF(0);
                }
    }

    template<typename TF>
    TF calc_momentum_2nd(
            const TF* restrict u, const TF* restrict v, const TF* restrict w,
            const TF* restrict dz, const TF itot_jtot_zsize,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk,
            Master& master)
    {
        using Finite_difference::O2::interp2;

        const int ii = 1;

        TF momentum = 0;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    momentum += (interp2(u[ijk], u[ijk+ii]) + interp2(v[ijk], v[ijk+jj]) + interp2(w[ijk], w[ijk+kk]))*dz[k];
                }

        master.sum(&momentum, 1);
        momentum /= itot_jtot_zsize;

        return momentum;
    }

    template<typename TF>
    TF calc_tke_2nd(
            const TF* restrict u, const TF* restrict v, const TF* restrict w,
            const TF* restrict dz, const TF itot_jtot_zsize,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk,
            Master& master)
    {
        using Finite_difference::O2::interp2;

        const int ii = 1;

        TF tke = 0;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    tke += ( interp2(u[ijk]*u[ijk], u[ijk+ii]*u[ijk+ii])
                           + interp2(v[ijk]*v[ijk], v[ijk+jj]*v[ijk+jj])
                           + interp2(w[ijk]*w[ijk], w[ijk+kk]*w[ijk+kk]))*dz[k];
                }

        master.sum(&tke, 1);
        tke /= itot_jtot_zsize;
        tke *= 0.5;

        return tke;
    }

    template<typename TF>
    TF calc_mass(
            const TF* restrict s,
            const TF* restrict dz, const TF itot_jtot_zsize,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jj, const int kk,
            Master& master)
    {
        TF mass = 0;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    mass += s[ijk]*dz[k];
                }

        master.sum(&mass, 1);
        mass /= itot_jtot_zsize;

        return mass;
    }

    std::pair<std::string, int> split_unit(const std::string s, const int pow)
    {
        std::string unit;
        int power;

        size_t delim = s.find_first_of("-123456789");
        if (delim == std::string::npos)
        {
            unit  = s;
            power = pow;
        }
        else
        {
            unit  = s.substr(0,delim);
            power = pow * std::stoi(s.substr(delim));
        }
        return std::make_pair(unit, power);
    }

    std::vector<std::pair<std::string, int>> get_units_vector(const std::string str, const int pow)
    {
        std::vector<std::string> fields;
        std::vector<std::pair<std::string, int>> unit;

        boost::split(fields, str, boost::is_any_of( " " ), boost::token_compress_on );
        for (auto& field : fields)
        {
            if (field != "-" && field != "")
                unit.push_back(split_unit(field, pow));
        }

        return unit;
    }

    template<typename TF>
    void reset_field(TF* const restrict fld, const int ncells)
    {
        for (int i=0; i<ncells; ++i)
            fld[i] = TF(0);
    }
}

template<typename TF>
Fields<TF>::Fields(Master& masterin, Grid<TF>& gridin, Soil_grid<TF>& soilgridin, Input& input) :
    master(masterin),
    grid(gridin),
    soil_grid(soilgridin),
    field3d_io(master, grid),
    field3d_operators(master, grid, *this),
    boundary_cyclic(master, grid)
{
    auto& gd = grid.get_grid_data();
    calc_mean_profs = false;

    // Initialize GPU pointers
    // rhoref_g  = 0;
    // rhorefh_g = 0;

    // obligatory parameters
    visc = input.get_item<TF>("fields", "visc", "");

    const std::string group_name = "default";

    // Initialize the passive scalars
    std::vector<std::string> slist = input.get_list<std::string>(
            "fields", "slist", "", std::vector<std::string>());
    for (auto& s : slist)
    {
        init_prognostic_field(s, s, "-", group_name, gd.sloc);
        sp.at(s)->visc = input.get_item<TF>("fields", "svisc", s);
    }

    // Initialize the basic set of fields.
    init_momentum_field("u", "U velocity", "m s-1", group_name, gd.uloc);
    init_momentum_field("v", "V velocity", "m s-1", group_name, gd.vloc);
    init_momentum_field("w", "Vertical velocity", "m s-1", group_name, gd.wloc);

    mp.at("u")->visc = visc;
    mp.at("v")->visc = visc;
    mp.at("w")->visc = visc;

    init_diagnostic_field("p", "Pressure", "Pa", group_name, gd.sloc);

    // Set a default of 4 temporary fields. Other classes can increase this number
    // before the init phase, where they are initialized in Fields::init()
    n_tmp_fields = 4;

    // Specify the masks that fields can provide / calculate
    available_masks.insert(available_masks.end(), {"default", "wplus", "wmin"});

    // Add user specified XY masks as available masks
    xymasklist = input.get_list<std::string>("stats", "xymasklist", "", std::vector<std::string>());
    available_masks.insert(available_masks.end(), xymasklist.begin(), xymasklist.end());
}

template<typename TF>
Fields<TF>::~Fields()
{
}

template<typename TF>
void Fields<TF>::init(Input& input, Dump<TF>& dump, Cross<TF>& cross, const Sim_mode sim_mode)
{
    auto& gd = grid.get_grid_data();

    boundary_cyclic.init();

    int nerror = 0;
    // ALLOCATE ALL THE FIELDS
    // allocate the prognostic velocity fields
    for (auto& it : mp)
        nerror += it.second->init();

    // allocate the velocity tendency fields
    for (auto& it : mt)
        nerror += it.second->init();

    // allocate the prognostic scalar fields
    for (auto& it : sp)
        nerror += it.second->init();

    // allocate the scalar tendency fields
    for (auto& it : st)
        nerror += it.second->init();

    // allocate the diagnostic scalars
    for (auto& it : sd)
        nerror += it.second->init();

    // allocate the prognostic soil fields
    for (auto& it : sps)
        nerror += it.second->init();
    for (auto& it : sts)
        nerror += it.second->init();

    // Allocate the prognostic 2d fields
    for (auto& it : ap2d)
        it.second->fld.resize(gd.ijcells);
    for (auto& it : at2d)
        it.second->fld.resize(gd.ijcells);

    // now that all classes have been able to set the minimum number of tmp fields, initialize them
    for (int i=0; i<n_tmp_fields; ++i)
        init_tmp_field();

    // allocate the tmp fields
    for (auto& tmp : atmp)
        nerror += tmp->init();

    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("Error allocating fields");

    rhoref .resize(gd.kcells);
    rhorefh.resize(gd.kcells);

    // \TODO Define a reference density. Needs to be replaced once anelastic is there
    // BvS: Always init rhoref at 1 for situation with e.g. thermo=0? For anelastic, overwrite it.
    std::fill(rhoref .begin(), rhoref .end(), 1.);
    std::fill(rhorefh.begin(), rhorefh.end(), 1.);

    // Create help arrays for statistics.
    umodel.resize(gd.kcells);
    vmodel.resize(gd.kcells);

    // Allocate user XY masks
    for (auto& mask : xymasklist)
        xymasks.emplace(mask, std::vector<TF>(gd.ijcells));

    // Set up output classes
    create_dump(dump);
    create_cross(cross);

    // Flag the data from the input that is not used outside of Init mode.
    if (sim_mode != Sim_mode::Init)
    {
        input.flag_as_used("fields", "rndamp", "");
        input.flag_as_used("fields", "rndexp", "");
        input.flag_as_used("fields", "rndz"  , "");

        // Also, flag the subspecified items.
        for (const auto& a : ap)
        {
            input.flag_as_used("fields", "rndamp", a.first);
            input.flag_as_used("fields", "rndexp", a.first);
            input.flag_as_used("fields", "rndz"  , a.first);
        }

        input.flag_as_used("fields", "rndseed", "");

        input.flag_as_used("fields", "vortexnpair", "");
        input.flag_as_used("fields", "vortexamp", "");
        input.flag_as_used("fields", "vortexaxis", "" );
    }
}

#ifndef USECUDA
template<typename TF>
void Fields<TF>::exec()
{
    // calculate the means for the prognostic scalars
    if (calc_mean_profs)
    {
        for (auto& it : ap)
            field3d_operators.calc_mean_profile(it.second->fld_mean.data(), it.second->fld.data());
    }
}
#endif


template<typename TF>
void Fields<TF>::create_dump(Dump<TF>& dump)
{
    // add the profiles to the columns
    if (dump.get_switch())
    {
        // Get global dump-list from dump.cxx
        std::vector<std::string>& dumplist_global = dump.get_dumplist();

        // Check if fields in dumplist are diagnostic fields, if not delete them and print warning
        std::vector<std::string>::iterator dumpvar = dumplist_global.begin();
        while (dumpvar != dumplist_global.end())
        {
            if (a.count(*dumpvar))
            {
                // Remove variable from global list, put in local list
                dumplist.push_back(*dumpvar);
                dumplist_global.erase(dumpvar); // erase() returns iterator of next element..
            }
            else
                ++dumpvar;
        }
    }
}

template<typename TF>
void Fields<TF>::create_cross(Cross<TF>& cross)
{

    if (cross.get_switch())
    {

        // Get global cross-list from cross.cxx
        std::vector<std::string>& crosslist_global = cross.get_crosslist();

        // Check different type of crosses and put them in their respective lists
        for (auto& it : ap)
        {
            check_added_cross(it.first, "",        crosslist_global, cross_simple);
            check_added_cross(it.first, "lngrad",  crosslist_global, cross_lngrad);
            check_added_cross(it.first, "bot",     crosslist_global, cross_bot);
            check_added_cross(it.first, "top",     crosslist_global, cross_top);
            check_added_cross(it.first, "fluxbot", crosslist_global, cross_fluxbot);
            check_added_cross(it.first, "fluxtop", crosslist_global, cross_fluxtop);
            check_added_cross(it.first, "path",    crosslist_global, cross_path);
        }

        for (auto& it : sd)
        {
            check_added_cross(it.first, "",        crosslist_global, cross_simple);
            check_added_cross(it.first, "lngrad",  crosslist_global, cross_lngrad);
        }
    }
}

template<typename TF>
void Fields<TF>::check_added_cross(
        const std::string& var,
        const std::string& type,
        std::vector<std::string>& crosslist,
        std::vector<std::string>& typelist)
{
    std::vector<std::string>::iterator position;

    position = std::find(crosslist.begin(), crosslist.end(), var + type);
    if (position != crosslist.end())
    {
        // Don't allow lngrad in 2nd order mode.
        if (!(type == "lngrad" && grid.get_spatial_order() == Grid_order::Second))
        {
            typelist.push_back(var);
            crosslist.erase(position);
        }
    }
}

template<typename TF>
std::shared_ptr<Field3d<TF>> Fields<TF>::get_tmp()
{
    std::shared_ptr<Field3d<TF>> tmp;

    #pragma omp critical
    {
        // In case of insufficient tmp fields, allocate a new one.
        if (atmp.empty())
        {
            init_tmp_field();
            tmp = atmp.back();
            tmp->init();
        }
        else
            tmp = atmp.back();

        atmp.pop_back();
    }
    return tmp;
}

template<typename TF>
void Fields<TF>::release_tmp(std::shared_ptr<Field3d<TF>>& tmp)
{
    if (tmp == nullptr)
        throw std::runtime_error("Cannot release a tmp field with value nullptr");

    atmp.push_back(std::move(tmp));
}

template<typename TF>
std::shared_ptr<std::vector<TF>> Fields<TF>::get_tmp_xy()
{
    auto& gd = grid.get_grid_data();
    std::shared_ptr<std::vector<TF>> tmp;

    #pragma omp critical
    {
        // In case of insufficient tmp fields, allocate a new one.
        if (atmp_xy.empty())
        {
            static int ntmp_xy = 0;
            ++ntmp_xy;
            std::string fldname = "tmp_xy" + std::to_string(ntmp_xy);
            std::string message = "Allocating temporary XY field: " + fldname;
            master.print_message(message);

            atmp_xy.push_back(std::make_shared<std::vector<TF>>(gd.ijcells));
            tmp = atmp_xy.back();
        }
        else
            tmp = atmp_xy.back();

        atmp_xy.pop_back();
    }

    return tmp;
}

template<typename TF>
void Fields<TF>::release_tmp_xy(std::shared_ptr<std::vector<TF>>& tmp)
{
    if (tmp == nullptr)
        throw std::runtime_error("Cannot release a tmp field with value nullptr");

    atmp_xy.push_back(std::move(tmp));
}

template<typename TF>
void Fields<TF>::get_mask(Stats<TF>& stats, std::string mask_name)
{
    // We don't have to do anything for the default mask
    if (mask_name == "default")
        return;

    auto& gd = grid.get_grid_data();

    // User XY masks
    if (xymasks.find(mask_name) != xymasks.end())
    {
        auto mask  = get_tmp();
        auto maskh = get_tmp();

        std::fill(mask->fld.begin(), mask->fld.end(), TF(0));
        std::fill(maskh->fld.begin(), maskh->fld.end(), TF(0));

        set_xy_mask(
                mask->fld.data(),
                maskh->fld.data(),
                xymasks.at(mask_name).data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        const TF threshold = TF(0.5);
        stats.set_mask_thres(mask_name, *mask, *maskh, threshold, Stats_mask_type::Plus);

        release_tmp(mask);
        release_tmp(maskh);
    }
    else
    {
        // Interpolate w to full level:
        auto wf = get_tmp();
        grid.interpolate_2nd(wf->fld.data(), mp.at("w")->fld.data(), gd.wloc.data(), gd.sloc.data());

        // Calculate masks
        const TF threshold = 0;
        if (mask_name == "wplus")
            stats.set_mask_thres(mask_name, *mp.at("w"), *wf, threshold, Stats_mask_type::Plus);
        else if (mask_name == "wmin")
            stats.set_mask_thres(mask_name, *mp.at("w"), *wf, threshold, Stats_mask_type::Min);

        release_tmp(wf);
    }
}

template<typename TF>
void Fields<TF>::exec_stats(Stats<TF>& stats)
{
    const TF no_offset = 0.;
    const TF no_threshold = 0.;

    stats.calc_stats("w", *mp.at("w"), no_offset, no_threshold);
    stats.calc_stats("u", *mp.at("u"), grid.utrans, no_threshold);
    stats.calc_stats("v", *mp.at("v"), grid.vtrans, no_threshold);

    for (auto& it : sp)
        stats.calc_stats(it.first, *it.second, no_offset, no_threshold);

    for (auto& it : sp)
        stats.calc_stats_2d(it.first + "_bot", it.second->fld_bot, no_offset);

    stats.calc_stats("p", *sd.at("p"), no_offset, no_threshold);

    // Calculate covariances
    for (auto& it1 : ap)
    {
        for (auto& it2 : ap)
        {
            for (int pow1 = 1; pow1<5; ++pow1)
            {
                for (int pow2 = 1; pow2<5; ++pow2)
                {
                    stats.calc_covariance(it1.first, *it1.second, no_offset, no_threshold, pow1,
                                          it2.first, *it2.second, no_offset, no_threshold, pow2);
                }
            }
        }
    }
}

template<typename TF>
void Fields<TF>::set_calc_mean_profs(bool sw)
{
    calc_mean_profs = sw;
}

template<typename TF>
void Fields<TF>::init_momentum_field(
        const std::string& fldname, const std::string& longname,
        const std::string& unit, const std::string& groupname,
        const std::array<int,3>& loc)
{
    if (mp.find(fldname) != mp.end())
    {
        std::string msg = fldname + " already exists";
        throw std::runtime_error(msg);
    }

    // Add a new prognostic momentum variable.
    mp[fldname] = std::make_shared<Field3d<TF>>(master, grid, fldname, longname, unit, groupname, loc);

    // Add a new tendency for momentum variable.
    std::string fldtname  = fldname + "t";
    std::string tunit     = simplify_unit(unit, "s-1");
    std::string tlongname = "Tendency of " + longname;
    mt[fldname] = std::make_shared<Field3d<TF>>(master, grid, fldtname, tlongname, tunit, groupname, loc);

    // Add the prognostic variable and its tendency to the collection
    // of all fields and tendencies.
    a [fldname] = mp[fldname];
    ap[fldname] = mp[fldname];
    at[fldname] = mt[fldname];
}

template<typename TF>
void Fields<TF>::init_prognostic_field(
        const std::string& fldname, const std::string& longname,
        const std::string& unit, const std::string& groupname,
        const std::array<int,3>& loc)
{
    if (sp.find(fldname)!=sp.end())
    {
        std::string msg = fldname + " already exists";
        throw std::runtime_error(msg);
    }

    // add a new scalar variable
    sp[fldname] = std::make_shared<Field3d<TF>>(master, grid, fldname, longname, unit, groupname, loc);

    // add a new tendency for scalar variable
    std::string fldtname  = fldname + "t";
    std::string tlongname = "Tendency of " + longname;
    std::string tunit     = simplify_unit(unit, "s-1");
    st[fldname] = std::make_shared<Field3d<TF>>(master, grid, fldtname, tlongname, tunit, groupname, loc);

    // add the prognostic variable and its tendency to the collection
    // of all fields and tendencies
    a [fldname] = sp[fldname];
    ap[fldname] = sp[fldname];
    at[fldname] = st[fldname];
}

template<typename TF>
void Fields<TF>::init_prognostic_soil_field(
        const std::string& fldname, const std::string& longname, const std::string& unit)
{
    if (sps.find(fldname)!=sps.end())
    {
        std::string msg = fldname + " already exists";
        throw std::runtime_error(msg);
    }

    // add a new scalar variable
    sps[fldname] = std::make_shared<Soil_field3d<TF>>(master, grid, soil_grid, fldname, longname, unit);

    // add a new tendency for scalar variable
    std::string fldtname  = fldname + "t";
    std::string tlongname = "Tendency of " + longname;
    std::string tunit     = simplify_unit(unit, "s-1");
    sts[fldname] = std::make_shared<Soil_field3d<TF>>(master, grid, soil_grid, fldtname, tlongname, tunit);
}

template<typename TF>
void Fields<TF>::init_prognostic_2d_field(const std::string& fldname)
{
    if (ap2d.find(fldname)!=ap2d.end())
    {
        std::string msg = fldname + " already exists";
        throw std::runtime_error(msg);
    }

    // Add a new scalar variable
    ap2d[fldname] = std::make_shared<Field2d<TF>>();
    at2d[fldname] = std::make_shared<Field2d<TF>>();
}

template<typename TF>
void Fields<TF>::init_diagnostic_field(
        const std::string& fldname, const std::string& longname,
        const std::string& unit, const std::string& groupname,
        const std::array<int,3>& loc)
{
    if (sd.find(fldname)!=sd.end())
    {
        std::string msg = fldname + " already exists";
        throw std::runtime_error(msg);
    }

    sd[fldname] = std::make_shared<Field3d<TF>>(master, grid, fldname, longname, unit, groupname, loc);
    a [fldname] = sd[fldname];
}

template<typename TF>
void Fields<TF>::init_tmp_field()
{
    static int ntmp = 0;
    ++ntmp;
    std::string fldname = "tmp" + std::to_string(ntmp);
    std::string longname = "";
    std::string unit = "";
    std::string group = "tmp_group";
    std::array<int,3> loc = {0,0,0};

    std::string message = "Allocating temporary field: " + fldname;
    master.print_message(message);
    atmp.push_back(std::make_shared<Field3d<TF>>(master, grid, fldname, longname, unit, group, loc));
}

#ifdef USECUDA
template<typename TF>
void Fields<TF>::init_tmp_field_g()
{
    static int ntmp = 0;
    ++ntmp;
    std::string fldname = "tmp_gpu" + std::to_string(ntmp);
    std::string longname = "";
    std::string unit = "";
    std::string group = "tmp_group";
    std::array<int,3> loc = {0,0,0};

    std::string message = "Allocating temporary field: " + fldname;
    master.print_message(message);
    atmp_g.push_back(std::make_shared<Field3d<TF>>(master, grid, fldname, longname, unit, group, loc));
}
#endif


template<typename TF>
void Fields<TF>::create(Input& input, Netcdf_file& input_nc)
{
    // Randomize the momentum
    randomize(input, "u", mp.at("u")->fld.data());
    randomize(input, "w", mp.at("w")->fld.data());

    // Only add perturbation to v in case of a 3d run.
    const Grid_data<TF>& gd = grid.get_grid_data();
    if (gd.jtot > 1)
        randomize(input, "v", mp.at("v")->fld.data());

    // Randomize the scalars
    for (auto& it : sp)
        randomize(input, it.first, it.second->fld.data());

    // Add Vortices
    add_vortex_pair(input);

    // Add the mean profiles to the fields
    add_mean_profs(input_nc);

    /*
    nerror += add_mean_prof(inputin, "u", mp.at("u")->data, grid.utrans);
    nerror += add_mean_prof(inputin, "v", mp.at("v")->data, grid.vtrans);

    for (auto& it : sp)
        nerror += add_mean_prof(inputin, it.first, it.second->data, 0.);
    */

    // Make sure the boundaries of w are zero. Non-zero initialization is not recoverable.
    int lbot = gd.kstart*gd.ijcells;
    int ltop = gd.kend  *gd.ijcells;

    for (int l=0; l<gd.ijcells; ++l)
    {
        mp.at("w")->fld[lbot+l] = 0.;
        mp.at("w")->fld[ltop+l] = 0.;
    }
}

template<typename TF>
void Fields<TF>::randomize(Input& input, std::string fld, TF* const restrict data)
{
    // Set mpiid as random seed to avoid having the same field at all procs
    int static seed = 0;

    if (!seed)
    {
        seed = input.get_item<int>("fields", "rndseed", "", 0);
        seed += master.get_mpiid() + 2;
        std::srand(seed);
    }

    const Grid_data<TF>& gd = grid.get_grid_data();

    const int jj = gd.icells;
    const int kk = gd.ijcells;

    // Look up the specific randomizer variables.
    rndamp = input.get_item<TF>("fields", "rndamp", fld, 0.);
    rndz   = input.get_item<TF>("fields", "rndz"  , fld, 0.);
    rndexp = input.get_item<TF>("fields", "rndexp", fld, 0.);

    if (rndz > gd.zsize)
    {
        std::string msg = "randomizer height rndz (" + std::to_string(rndz) + ") higher than domain top (" + std::to_string(gd.zsize) +")";
        throw std::runtime_error(msg);
    }

    // Find the location of the randomizer height.
    int kendrnd = gd.kstart;
    while (gd.z[kendrnd] < rndz)
        ++kendrnd;

    // Issue a warning if the randomization depth is larger than zero, but less than the first model level.
    if (kendrnd == gd.kstart && rndz > 0.)
        master.print_warning("randomization depth is less than the height of the first model level\n");

    for (int k=gd.kstart; k<kendrnd; ++k)
    {
        const TF rndfac = std::pow((rndz-gd.z[k])/rndz, rndexp);
        for (int j=gd.jstart; j<gd.jend; ++j)
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                data[ijk] = rndfac * rndamp * ((TF) std::rand() / (TF) RAND_MAX - 0.5);
            }
    }
}

namespace
{
    template<typename TF>
    void add_mean_prof_to_field(TF* restrict const data,
                                const TF* restrict const dataprof,
                                const TF offset,
                                const int istart, const int iend,
                                const int jstart, const int jend,
                                const int kstart, const int kend,
                                const int jj, const int kk)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    data[ijk] += dataprof[k-kstart] - offset;
                }
    }
}

template<typename TF>
void Fields<TF>::add_mean_profs(Netcdf_handle& input_nc)
{
    const Grid_data<TF>& gd = grid.get_grid_data();
    std::vector<TF> prof(gd.ktot);

    const std::vector<int> start = {0};
    const std::vector<int> count = {gd.ktot};

    Netcdf_group& group_nc = input_nc.get_group("init");
    group_nc.get_variable(prof, "u", start, count);

    add_mean_prof_to_field<TF>(mp.at("u")->fld.data(), prof.data(), grid.utrans,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    group_nc.get_variable(prof, "v", start, count);
    add_mean_prof_to_field<TF>(mp.at("v")->fld.data(), prof.data(), grid.vtrans,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells);

    for (auto& f : sp)
    {
        group_nc.get_variable(prof, f.first, start, count);
        add_mean_prof_to_field<TF>(f.second->fld.data(), prof.data(), 0.,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells);

        // For cold start, initialise surface values with first model level values
        std::fill(f.second->fld_bot.begin(), f.second->fld_bot.end(), prof[0]);
    }
}

template<typename TF>
void Fields<TF>::add_vortex_pair(Input& inputin)
{
    auto& gd = grid.get_grid_data();

    // Optional parameters.
    vortexnpair = inputin.get_item<int>        ("fields", "vortexnpair", "", 0    );
    vortexamp   = inputin.get_item<TF>         ("fields", "vortexamp"  , "", 1.e-3);
    vortexaxis  = inputin.get_item<std::string>("fields", "vortexaxis" , "", "y"  );

    // Add a double vortex to the initial conditions.
    const double pi = std::acos((double)-1.);

    if (vortexnpair > 0)
    {
        if (vortexaxis == "y")
            for (int k=gd.kstart; k<gd.kend; ++k)
                for (int j=gd.jstart; j<gd.jend; ++j)
                    for (int i=gd.istart; i<gd.iend; ++i)
                    {
                        const int ijk = i + j*gd.icells + k*gd.ijcells;
                        mp.at("u")->fld[ijk] +=  vortexamp*std::sin(vortexnpair*2.*pi*(gd.xh[i])/gd.xsize)*std::cos(pi*gd.z [k]/gd.zsize);
                        mp.at("w")->fld[ijk] += -vortexamp*std::cos(vortexnpair*2.*pi*(gd.x [i])/gd.xsize)*std::sin(pi*gd.zh[k]/gd.zsize);
                    }
        else if (vortexaxis == "x")
            for (int k=gd.kstart; k<gd.kend; ++k)
                for (int j=gd.jstart; j<gd.jend; ++j)
                    for (int i=gd.istart; i<gd.iend; ++i)
                    {
                        const int ijk = i + j*gd.icells + k*gd.ijcells;
                        mp.at("v")->fld[ijk] +=  vortexamp*std::sin(vortexnpair*2.*pi*(gd.yh[j])/gd.ysize)*std::cos(pi*gd.z [k]/gd.zsize);
                        mp.at("w")->fld[ijk] += -vortexamp*std::cos(vortexnpair*2.*pi*(gd.y [j])/gd.ysize)*std::sin(pi*gd.zh[k]/gd.zsize);
                    }
    }
}

template <typename TF>
void Fields<TF>::create_stats(Stats<TF>& stats)
{
    const std::string group_name = "default";

    const std::vector<std::string> stat_op_def = {"mean", "2", "3", "4", "w", "grad", "diff", "flux", "path"};
    const std::vector<std::string> stat_op_w = {"mean", "2", "3", "4"};
    const std::vector<std::string> stat_op_p = {"mean", "2", "w", "grad"};

    // Add the profiles to te statistics
    if (stats.get_switch())
    {
        for (auto& it : ap)
        {
            if (it.first == "w")
                stats.add_profs(*it.second, "zh", stat_op_w, it.second->group);
            else
                stats.add_profs(*it.second, "z", stat_op_def, it.second->group);
        }
        stats.add_profs(*sd.at("p"), "z", stat_op_p, group_name);

        // Covariances
        for (auto& it1 : ap)
        {
            for (auto& it2 : ap)
            {
                std::string locstring;
                if (it2.first == "w")
                    locstring = "zh";
                else
                    locstring = "z";

                stats.add_covariance(*it1.second, *it2.second, locstring);
            }
        }
    }

    // Add time series of scalar surface values
    for (auto& it : sp)
        stats.add_time_series(it.first + "_bot", "Surface " + it.second->longname, it.second->unit, it.second->group);
}

template<typename TF>
void Fields<TF>::create_column(Column<TF>& column)
{
    // add the profiles to the columns
    if (column.get_switch())
    {
        // add variables to the statistics
        column.add_prof(mp.at("u")->name, mp.at("u")->longname, mp.at("u")->unit, "z" );
        column.add_prof(mp.at("v")->name, mp.at("v")->longname, mp.at("v")->unit, "z" );
        column.add_prof(mp.at("w")->name, mp.at("w")->longname, mp.at("w")->unit, "zh");

        for (auto& it : sp)
            column.add_prof(it.first, it.second->longname, it.second->unit, "z");

        for (auto& it : sp)
            column.add_time_series(it.first + "_bot", "Surface " + it.second->longname, it.second->unit);

        column.add_prof(sd.at("p")->name, sd.at("p")->longname, sd.at("p")->unit, "z");
    }
}

template<typename TF>
void Fields<TF>::save(int n)
{
    auto& gd = grid.get_grid_data();
    const TF no_offset = 0.;

    auto tmp1 = get_tmp();
    auto tmp2 = get_tmp();

    int nerror = 0;
    for (auto& f : ap)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", f.second->name.c_str(), n);
        master.print_message("Saving \"%s\" ... ", filename);

        // The offset is kept at zero, because otherwise bitwise identical restarts are not possible.
        if (field3d_io.save_field3d(
                    f.second->fld.data(),
                    tmp1->fld.data(), tmp2->fld.data(),
                    filename, no_offset,
                    gd.kstart, gd.kend))
        {
            master.print_message("FAILED\n");
            ++nerror;
        }
        else
        {
            master.print_message("OK\n");
        }
    }

    release_tmp(tmp1);
    release_tmp(tmp2);

    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("Error saving 3D fields");
}

template<typename TF>
void Fields<TF>::load(int n)
{
    auto& gd = grid.get_grid_data();
    const TF no_offset = 0.;

    auto tmp1 = get_tmp();
    auto tmp2 = get_tmp();

    int nerror = 0;

    for (auto& f : ap)
    {
        // The offset is kept at zero, otherwise bitwise identical restarts is not possible.
        char filename[256];
        std::sprintf(filename, "%s.%07d", f.second->name.c_str(), n);
        master.print_message("Loading \"%s\" ... ", filename);

        if (field3d_io.load_field3d(
                    f.second->fld.data(),
                    tmp1->fld.data(), tmp2->fld.data(),
                    filename, no_offset,
                    gd.kstart, gd.kend))
        {
            master.print_message("FAILED\n");
            ++nerror;
        }
        else
        {
            master.print_message("OK\n");
        }
    }

    // Load surface (XY) masks
    for (auto& mask : xymasks)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", mask.first.c_str(), n);
        master.print_message("Loading \"%s\" ... ", filename);

        if (field3d_io.load_xy_slice(
                mask.second.data(), tmp1->fld.data(), filename))
        {
            master.print_message("FAILED\n");
            ++nerror;
        }
        else
            master.print_message("OK\n");
    }

    release_tmp(tmp1);
    release_tmp(tmp2);

    master.sum(&nerror, 1);

    if (nerror)
        throw std::runtime_error("Error loading fields");
}

#ifndef USECUDA
template<typename TF>
TF Fields<TF>::check_momentum()
{
    auto& gd = grid.get_grid_data();
    return calc_momentum_2nd(
            mp.at("u")->fld.data(), mp.at("v")->fld.data(), mp.at("w")->fld.data(),
            gd.dz.data(), gd.itot*gd.jtot*gd.zsize,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells,
            master);
}
#endif

#ifndef USECUDA
template<typename TF>
TF Fields<TF>::check_tke()
{
    auto& gd = grid.get_grid_data();
    return calc_tke_2nd(
            mp.at("u")->fld.data(), mp.at("v")->fld.data(), mp.at("w")->fld.data(),
            gd.dz.data(), gd.itot*gd.jtot*gd.zsize,
            gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
            gd.icells, gd.ijcells,
            master);
}
#endif

#ifndef USECUDA
template<typename TF>
TF Fields<TF>::check_mass()
{
    auto& gd = grid.get_grid_data();

    auto it = sp.begin();
    if (sp.begin() != sp.end())
        return calc_mass(
                it->second->fld.data(),
                gd.dz.data(), gd.itot*gd.jtot*gd.zsize,
                gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                gd.icells, gd.ijcells,
                master);
    else
        return 0.;
}
#endif

template<typename TF>
void Fields<TF>::exec_cross(Cross<TF>& cross, unsigned long iotime)
{
    for (auto& it : cross_simple)
        cross.cross_simple(a.at(it)->fld.data(), a.at(it)->name, iotime, a.at(it)->loc);

    for (auto& it : cross_lngrad)
        cross.cross_lngrad(a.at(it)->fld.data(), a.at(it)->name+"lngrad", iotime);

    for (auto& it : cross_fluxbot)
        cross.cross_plane(a.at(it)->flux_bot.data(), a.at(it)->name+"fluxbot", iotime);

    for (auto& it : cross_fluxtop)
        cross.cross_plane(a.at(it)->flux_top.data(), a.at(it)->name+"fluxtop", iotime);

    for (auto& it : cross_bot)
        cross.cross_plane(a.at(it)->fld_bot.data(), a.at(it)->name+"bot", iotime);

    for (auto& it : cross_top)
        cross.cross_plane(a.at(it)->fld_top.data(), a.at(it)->name+"top", iotime);

    for (auto& it : cross_path)
        cross.cross_path(a.at(it)->fld.data(), a.at(it)->name+"path", iotime);
}

template<typename TF>
void Fields<TF>::exec_dump(Dump<TF>& dump, unsigned long iotime)
{
    for (auto& it : dumplist)
        dump.save_dump(a.at(it)->fld.data(), a.at(it)->name, iotime);
}

#ifndef USECUDA
template<typename TF>
void Fields<TF>::exec_column(Column<TF>& column)
{
    const TF no_offset = 0.;

    column.calc_column("u", mp.at("u")->fld.data(), grid.utrans);
    column.calc_column("v", mp.at("v")->fld.data(), grid.vtrans);
    column.calc_column("w", mp.at("w")->fld.data(), no_offset);

    for (auto& it : sp)
        column.calc_column(it.first, it.second->fld.data(), no_offset);

    for (auto& it : sp)
        column.calc_time_series(it.first + "_bot", it.second->fld_bot.data(), no_offset);

    column.calc_column("p", sd.at("p")->fld.data(), no_offset);
}
#endif

template<typename TF>
bool Fields<TF>::has_mask(std::string mask_name)
{
    if (std::find(available_masks.begin(), available_masks.end(), mask_name) != available_masks.end())
        return true;
    else
        return false;
}

template<typename TF>
std::string Fields<TF>::simplify_unit(const std::string str1, const std::string str2, const int pow1, const int pow2)
{
    std::vector<std::pair<std::string, int>> unit1, unit2;

    //Split each string in separate unit strings; split those in pairs of unit and power
    unit1 = get_units_vector(str1, pow1);
    unit2 = get_units_vector(str2, pow2);

    //Loop through units to find matches; in which case add the powers
    int unit1_size = unit1.size();
    for (auto& u2 : unit2)
    {
        int i;
        for (i = 0 ; i < unit1_size; i++)
        {
            if (u2.first == unit1[i].first)
            {
                if (u2.first == "kg") //Special case: there could be a kg/kg here to simplify
                {
                    int j;
                    for (j = i++ ; j < unit1_size; j++)
                    {
                        if (u2.first == unit1[j].first)
                            break;
                    }
                    if (j == unit1_size)
                        unit1[i].second += u2.second;
                    else if (u2.second * unit1[j].second < 0)
                        unit1[j].second += u2.second;
                    else
                        unit1[i].second += u2.second;
                    break;
                }
                else
                {
                    unit1[i].second += u2.second;
                    break;
                }
            }
        }
        if (i == unit1_size)
            unit1.push_back(u2);
    }

    // Remove the entries with zero power
    for (auto u1 = unit1.begin(); u1 != unit1.end(); )
    {
        if ((*u1).second  == 0)
            u1 = unit1.erase(u1);
        else
            ++u1;
    }

    //Convert pairs back into strings
    std::string output;
    if (unit1.size() == 0)
    {
        output = "-";
    }
    else
    {
        std::ostringstream ostream;
        for (auto& u1 : unit1)
        {
            if (u1.second == 1)
            {
                ostream << u1.first << " ";
            }
            else
            {
                ostream << u1.first << u1.second << " ";
            }
        }
        output = ostream.str();
        output.erase(output.end()-1); //remove final space
    }

    return output;
}

template<typename TF>
void Fields<TF>::reset_tendencies()
{
    auto& gd = grid.get_grid_data();

    for (auto& fld3d : at)
        reset_field(fld3d.second->fld.data(), gd.ncells);
}


template class Fields<double>;
template class Fields<float>;
