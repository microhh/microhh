/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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

//#include <cstdio>
#include <algorithm>
#include <iostream>
#include <math.h>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "thermo.h"
#include "stats.h"
#include "decay.h"

namespace
{
    template<typename TF>
    void enforce_exponential_decay(
            TF* restrict tend, const TF* var, const TF decaytime, const TF dt,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        const TF rate = 1./(std::max(decaytime, dt));

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    tend[ijk] -= rate * var[ijk];
                }
    }

    template<typename TF>
    void enforce_reset_ql(
            TF* restrict fld, const TF* ql,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;
                    if (ql[ijk] > 0.)
                        fld[ijk] = 1.;
                }
    }
}

template<typename TF>
Decay<TF>::Decay(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin, Thermo<TF>& thermoin):
    master(masterin), grid(gridin), fields(fieldsin), thermo(thermoin)
{
}

template <typename TF>
Decay<TF>::~Decay()
{
}

template <typename TF>
void Decay<TF>::init(Input& inputin)
{
    for (auto& it : fields.st)
    {
        const std::string type = inputin.get_item<std::string>("decay", "swdecay", it.first, "0");
        if (type == "0")
        {
            // Cycle to avoid reading unneeded namelist options.
            continue;
        }
        else if (type == "exponential")
        {
            dmap[it.first].type = Decay_type::exponential;
            dmap[it.first].timescale = inputin.get_item<TF>("decay", "timescale", it.first);
        }
        else
            throw std::runtime_error("Invalid option for \"decay type\"");
    }

    // Read the setting only if map is not empty.
    if (!dmap.empty())
        nstd_couvreux = inputin.get_item<TF>("decay", "nstd_couvreux", "", 1.);

    for (auto& it : fields.st)
    {
        const std::string type = inputin.get_item<std::string>("decay", "swreset", it.first, "0");
        if (type == "0")
        {
            // Cycle to avoid reading unneeded namelist options.
            continue;
        }
        else if (type == "ql")
        {
            dmap[it.first].reset_type= Reset_type::ql;
        }
        else
            throw std::runtime_error("Invalid option for \"reset type\"");
    }


}

template <typename TF>
void Decay<TF>::create(Input& inputin, Stats<TF>& stats)
{
}

#ifndef USECUDA
template <typename TF>
void Decay<TF>::exec(double dt, Stats<TF>& stats, Thermo<TF>& thermo)
{
    auto& gd = grid.get_grid_data();
    for (auto& it : dmap)
    {
        if (it.second.type == Decay_type::exponential)
        {
            enforce_exponential_decay<TF>(
                    fields.st.at(it.first)->fld.data(), fields.sp.at(it.first)->fld.data(), it.second.timescale, dt,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
            stats.calc_tend(*fields.st.at(it.first), tend_name);
        }
    }
    for (auto& it : dmap)
    {
        if (it.second.reset_type == Reset_type::ql)
        {
            auto ql = fields.get_tmp();
            thermo.get_thermo_field(*ql, "ql", true, false);
            enforce_reset_ql<TF>(
                    fields.sp.at(it.first)->fld.data(), ql,
                    gd.istart, gd.iend, gd.jstart, gd.jend, gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
            // stats.calc_tend(*fields.st.at(it.first), tend_name);
        }
    }
}
#endif

template<typename TF>
bool Decay<TF>::has_mask(std::string name)
{
    if (std::find(available_masks.begin(), available_masks.end(), name) != available_masks.end())
        return true;
    else
        return false;
}

template<typename TF>
void Decay<TF>::get_mask(Stats<TF>& stats, std::string mask_name)
{
    auto& gd = grid.get_grid_data();

    if (mask_name == "couvreux")
    {
        if (fields.sp.find("couvreux") == fields.sp.end())
        {
            std::string message = "Couvreux mask not available without couvreux scalar";
            throw std::runtime_error(message);
        }

        TF ijtot = static_cast<TF>(gd.itot*gd.jtot);

        auto couvreux = fields.get_tmp();
        auto couvreuxh = fields.get_tmp();

        // Calculate mean and variance
        for (int k=gd.kstart; k<gd.kend; ++k)
        {
            TF mean = 0.;
            TF var = 0.;
            for (int j=gd.jstart; j<gd.jend; ++j)
            {
                #pragma ivdep
                for (int i=gd.istart; i<gd.iend; ++i)
                {
                    const int ijk = i + j*gd.icells + k*gd.ijcells;
                    mean += fields.sp.at("couvreux")->fld[ijk];
                    var  += fields.sp.at("couvreux")->fld[ijk]*fields.sp.at("couvreux")->fld[ijk];
                }
            }
            master.sum(&mean,1);
            master.sum(&var,1);
            mean /= ijtot;
            var /= ijtot;
            var -= mean*mean;
            TF std = sqrt(var);

            for (int j=gd.jstart; j<gd.jend; ++j)
            {
                #pragma ivdep
                for (int i=gd.istart; i<gd.iend; ++i)
                {
                    const int ijk = i + j*gd.icells + k*gd.ijcells;
                    couvreux->fld[ijk] = fields.sp.at("couvreux")->fld[ijk] - mean - nstd_couvreux * std;
                }
            }
        }
        grid.interpolate_2nd(couvreuxh->fld.data(), couvreux->fld.data(), gd.sloc.data(), gd.wloc.data());

        // Calculate masks
        TF threshold = 0.;
        stats.set_mask_thres(mask_name, *couvreux, *couvreuxh, threshold, Stats_mask_type::Plus);

        fields.release_tmp(couvreux);
        fields.release_tmp(couvreuxh);
    }
    else if (mask_name == "cldshell" || mask_name == "cldshelldown")
    {
        auto shellh = fields.get_tmp();
        grid.interpolate_2nd(shellh->fld.data(), fields.sp.at("shell")->fld.data(), gd.sloc.data(), gd.wloc.data());
        stats.set_mask_thres(mask_name, *fields.sp.at("shell"), *shellh, 0.5, Stats_mask_type::Plus);
        fields.release_tmp(shellh);

        auto ql = fields.get_tmp();
        auto qlh = fields.get_tmp();
        thermo.get_thermo_field(*ql, "ql", true, false);
        thermo.get_thermo_field(*qlh, "ql_h", true, false);

        stats.set_mask_thres(mask_name, *ql, *qlh, 1e-10, Stats_mask_type::Min);

        fields.release_tmp(ql);
        fields.release_tmp(qlh);

        if (mask_name == "cldshelldown")
        {
            auto wf = fields.get_tmp();
            grid.interpolate_2nd(wf->fld.data(), fields.mp.at("w")->fld.data(), gd.wloc.data(), gd.sloc.data());
            stats.set_mask_thres(mask_name, *wf, *fields.mp.at("w"), 0., Stats_mask_type::Min);
            fields.release_tmp(wf);
        }
    }

    else
    {
        std::string message = "Decay cannot provide mask: \"" + mask_name +"\"";
        throw std::runtime_error(message);
    }
}


#ifdef FLOAT_SINGLE
template class Decay<float>;
#else
template class Decay<double>;
#endif
