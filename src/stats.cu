/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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
// #include "background_profs.h"
#include "fields.h"
#include "stats.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "timeloop.h"
#include "advec.h"
#include "diff.h"
#include "netcdf_interface.h"
#include "stats_functions.h"
#include "tools.h"

using namespace Stats_functions;

#ifdef USECUDA
template<typename TF>
void Stats<TF>::calc_stats_mean(
        const std::string& varname, const Field3d<TF>& fld, const TF offset)
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;

    // For 2D field excluding ghost cells
    int gridi = gd.itot/blocki + (gd.itot%blocki > 0);
    int gridj = gd.jtot/blockj + (gd.jtot%blockj > 0);
    dim3 gridGPU (gridi,  gridj,  1);
    dim3 blockGPU(blocki, blockj, 1);

    unsigned int flag;
    const int* nmask;
    std::string name;
    auto slice = fields.get_tmp_xy_g();
    // Calc mean of atmospheric variables
    if (std::find(varlist.begin(), varlist.end(), varname) != varlist.end())
    {
        for (auto& m : masks)
        {
            set_flag(flag, nmask, m.second, fld.loc[2]);
            for (int k = gd.kstart; k < gd.kend; ++k)
            {
                if (nmask[k] > 0)
                {
                    int kk = k * gd.ijcells;
                    apply_mask_g<<<gridGPU, blockGPU>>>(slice->data(), &(fld.fld_g[kk]), &(mfield_g[kk]), flag, gd.icells, gd.jcells);
                    m.second.profs.at(varname).data[k] = ( gd.itot * gd.jtot ) / nmask[k] * field3d_operators.calc_mean_2d_g(slice->data());  //Assuming that nmask is per process
                }              
            }
            master.sum(m.second.profs.at(varname).data.data(), gd.kcells);

            // Add the offset.
            for (auto& value : m.second.profs.at(varname).data)
                value += offset;

            set_fillvalue_prof(m.second.profs.at(varname).data.data(), nmask, gd.kstart, gd.kcells);
        }
    }
    fields.release_tmp_xy_g(slice);
}
#endif

#ifdef USECUDA
template<typename TF>
void Stats<TF>::initialize_masks()
{
    auto& gd = grid.get_grid_data();
    unsigned int flagmax = 0;
    using namespace Tools_g;

    for (auto& it : masks){
        flagmax += it.second.flag + it.second.flagh;
        std::cout << "flagmax: " << flagmax << std::endl;
    }
    const int nblock = 256;
    int ngrid  = gd.ncells/nblock + (gd.ncells%nblock > 0);
    set_to_val<<<ngrid, nblock>>>(mfield_g.data(), gd.ncells, flagmax);

    ngrid  = gd.ijcells/nblock + (gd.ncells%nblock > 0);
    set_to_val<<<ngrid, nblock>>>(mfield_bot_g.data(), gd.ijcells, flagmax);
}

template<typename TF>
void Stats<TF>::finalize_masks()
{
    using namespace Tools_g;
    auto& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 gridGPU(gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);

    boundary_cyclic.exec_g(mfield_g.data());
    boundary_cyclic.exec_2d_g(mfield_bot_g.data());

    auto slice = fields.get_tmp_xy_g();
    auto ones = fields.get_tmp_xy_g();
    
    set_to_val<<<gridGPU, blockGPU>>>(ones->data(), gd.ijcells, TF(1.0));
    for (auto& it : masks)
    {
        unsigned int flag = it.second.flag;
        unsigned int flagh = it.second.flagh;

        for (int k = gd.kstart; k < gd.kend; ++k)
        {
            int kk = k * gd.ijcells;
            apply_mask_g<<<gridGPU, blockGPU>>>(slice->data(), ones->data(), &(mfield_g[kk]), flag, gd.icells, gd.jcells);
            it.second.nmask[k] = field3d_operators.calc_mean_2d_g(slice->data());  
            apply_mask_g<<<gridGPU, blockGPU>>>(slice->data(), ones->data(), &(mfield_g[kk]), flagh, gd.icells, gd.jcells);
            it.second.nmaskh[k] = ( gd.itot * gd.jtot ) * field3d_operators.calc_mean_2d_g(slice->data());  
        }
        int ijk = gd.istart + gd.jstart * gd.icells + gd.kstart * gd.ijcells;
        std::cout << "nmask: "<< it.first << it.second.nmask[gd.kstart]<< " " << slice->data()[ijk]<< std::endl;
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
    fields.release_tmp_xy_g(slice);
    fields.release_tmp_xy_g(ones);
}

template<typename TF>
void Stats<TF>::set_mask_thres(
        std::string mask_name, Field3d<TF>& fld, Field3d<TF>& fldh, TF threshold, Stats_mask_type mode)
{
    auto& gd = grid.get_grid_data();
    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.icells/blocki + (gd.icells%blocki > 0);
    const int gridj  = gd.jcells/blockj + (gd.jcells%blockj > 0);

    dim3 gridGPU(gridi, gridj, gd.kcells);
    dim3 blockGPU(blocki, blockj, 1);


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
        calc_mask_thres_g<TF, Stats_mask_type::Plus><<<gridGPU, blockGPU>>>(
                mfield_g.data(), mfield_bot_g.data(), flag, flagh,
                fld.fld_g.data(), fldh.fld_g.data(), fldh.fld_bot_g.data(), threshold,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.jcells, gd.ijcells, gd.kcells);
    else if (mode == Stats_mask_type::Min)
        calc_mask_thres_g<TF, Stats_mask_type::Min><<<gridGPU, blockGPU>>>(
                mfield_g.data(), mfield_bot_g.data(), flag, flagh,
                fld.fld_g.data(), fldh.fld_g.data(), fldh.fld_bot_g.data(), threshold,
                gd.istart, gd.jstart, gd.kstart,
                gd.iend,   gd.jend,   gd.kend,
                gd.icells, gd.jcells, gd.ijcells, gd.kcells);
    else
        throw std::runtime_error("Invalid mask type in set_mask_thres()");
}
#endif


#ifdef FLOAT_SINGLE
template class Stats<float>;
#else
template class Stats<double>;
#endif