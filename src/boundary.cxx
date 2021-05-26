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
#include <algorithm>
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "diff.h"
#include "defines.h"
#include "timeloop.h"
#include "timedep.h"
#include "finite_difference.h"

#include "boundary_cyclic.h"

// Boundary schemes.
#include "boundary.h"
#include "boundary_surface.h"
#include "boundary_surface_bulk.h"
#include "boundary_surface_lsm.h"

namespace
{
    template<typename TF>
    void set_bc(TF* const restrict a, TF* const restrict agrad, TF* const restrict aflux,
                const Boundary_type sw, const TF aval, const TF visc, const TF offset,
                const int icells, const int jcells)
    {
        const int jj = icells;

        if (sw == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    a[ij] = aval - offset;
                }
        }
        else if (sw == Boundary_type::Neumann_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    agrad[ij] = aval;
                    aflux[ij] = -aval*visc;
                }
        }
        else if (sw == Boundary_type::Flux_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij = i + j*jj;
                    aflux[ij] = aval;
                    agrad[ij] = -aval/visc;
                }
        }
    }
}

template<typename TF>
Boundary<TF>::Boundary(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin),
    grid(gridin),
    fields(fieldsin),
    boundary_cyclic(master, grid),
    field3d_io(master, grid)
{
    swboundary = "default";
}

template<typename TF>
Boundary<TF>::~Boundary()
{
    // empty the map
    // CvH: is this necessary?
    sbc.clear();

    // clean up time dependent data
    // for (auto& i : timedepdata)
    //     delete[] i.second;
}

template<typename TF>
std::string Boundary<TF>::get_switch()
{
    return swboundary;
}


template<typename TF>
void Boundary<TF>::process_bcs(Input& input)
{
    std::string swbot = input.get_item<std::string>("boundary", "mbcbot", "");
    std::string swtop = input.get_item<std::string>("boundary", "mbctop", "");

    ubot = input.get_item<TF>("boundary", "ubot", "", 0.);
    utop = input.get_item<TF>("boundary", "utop", "", 0.);
    vbot = input.get_item<TF>("boundary", "vbot", "", 0.);
    vtop = input.get_item<TF>("boundary", "vtop", "", 0.);

    // set the bottom bc
    if (swbot == "noslip")
        mbcbot = Boundary_type::Dirichlet_type;
    else if (swbot == "freeslip" )
        mbcbot = Boundary_type::Neumann_type;
    else if (swbot == "neumann" )
        mbcbot = Boundary_type::Neumann_type;
    else if (swbot == "ustar")
        mbcbot = Boundary_type::Ustar_type;
    else
    {
        std::string msg = swbot + " is an illegal value for mbcbot";
        throw std::runtime_error(msg);
    }

    // set the top bc
    if (swtop == "noslip")
        mbctop = Boundary_type::Dirichlet_type;
    else if (swtop == "freeslip")
        mbctop = Boundary_type::Neumann_type;
    else if (swtop == "neumann")
        mbctop = Boundary_type::Neumann_type;
    else if (swtop == "freeslip")
        mbctop = Boundary_type::Neumann_type;
    else if (swtop == "off")
        mbctop = Boundary_type::Off_type;
    else if (swtop == "ustar")
        mbctop = Boundary_type::Ustar_type;
    else
    {
        std::string msg = swtop + " is an illegal value for mbctop";
        throw std::runtime_error(msg);
    }

    // read the boundaries per field
    for (auto& it : fields.sp)
    {
        sbc.emplace(it.first, Field3dBc<TF>());
        swbot = input.get_item<std::string>("boundary", "sbcbot", it.first);
        swtop = input.get_item<std::string>("boundary", "sbctop", it.first);
        sbc.at(it.first).bot = input.get_item<TF>("boundary", "sbot", it.first);
        sbc.at(it.first).top = input.get_item<TF>("boundary", "stop", it.first);

        // set the bottom bc
        if (swbot == "dirichlet")
            sbc.at(it.first).bcbot = Boundary_type::Dirichlet_type;
        else if (swbot == "flux")
            sbc.at(it.first).bcbot = Boundary_type::Flux_type;
        else if (swbot == "neumann")
            sbc.at(it.first).bcbot = Boundary_type::Neumann_type;
        else
        {
            std::string msg = swbot + " is an illegal value for sbcbot";
            throw std::runtime_error(msg);
        }

        // set the top bc
        if (swtop == "dirichlet")
            sbc.at(it.first).bctop = Boundary_type::Dirichlet_type;
        else if (swtop == "neumann")
            sbc.at(it.first).bctop = Boundary_type::Neumann_type;
        else if (swtop == "flux")
            sbc.at(it.first).bctop = Boundary_type::Flux_type;
        else if (swtop == "off")
            sbc.at(it.first).bctop = Boundary_type::Off_type;
        else
        {
            std::string msg = swbot + " is an illegal value for sbctop";
            throw std::runtime_error(msg);
        }
    }

    sbot_2d_list = input.get_list<std::string>("boundary", "sbot_2d_list", "", std::vector<std::string>());
}

template<typename TF>
void Boundary<TF>::init(Input& input, Thermo<TF>& thermo)
{
    // Read the boundary information from the ini files, it throws at error.
    process_bcs(input);

    // there is no option (yet) for prescribing ustar without surface model
    if (mbcbot == Boundary_type::Ustar_type || mbctop == Boundary_type::Ustar_type)
    {
        throw std::runtime_error("Cannot use ustar bc for default boundary");
    }

    // Initialize the boundary cyclic.
    boundary_cyclic.init();
}

template<typename TF>
void Boundary<TF>::create(
        Input& input, Netcdf_handle& input_nc,
        Stats<TF>& stats, Column<TF>& column, Cross<TF>& cross)
{
    process_time_dependent(input, input_nc);
}


template<typename TF>
void Boundary<TF>::process_time_dependent(Input& input, Netcdf_handle& input_nc)
{
    // get the list of time varying variables
    bool swtimedep = input.get_item<bool>("boundary", "swtimedep"  , "", false);
    std::vector<std::string> timedeplist = input.get_list<std::string>("boundary", "timedeplist", "", std::vector<std::string>());

    if (swtimedep)
    {
        if (!sbot_2d_list.empty())
            master.print_warning("Provided 2D sbot fields are potentially overwritten by timedep");

        // Create temporary list to check which entries are used.
        std::vector<std::string> tmplist = timedeplist;

        // See if there is data available for the surface boundary conditions.
        for (auto& it : fields.sp)
        {
            std::string name = it.first+"_sbot";
            if (std::find(timedeplist.begin(), timedeplist.end(), name) != timedeplist.end())
            {
                // Process the time dependent data.
                tdep_bc.emplace(it.first, new Timedep<TF>(master, grid, name, true));
                tdep_bc.at(it.first)->create_timedep(input_nc);

                // Remove the item from the tmplist.
                std::vector<std::string>::iterator ittmp = std::find(tmplist.begin(), tmplist.end(), name);
                if (ittmp != tmplist.end())
                    tmplist.erase(ittmp);
            }
        }

        // Display a warning for the non-supported.
        for (std::vector<std::string>::const_iterator ittmp=tmplist.begin(); ittmp!=tmplist.end(); ++ittmp)
            master.print_warning("%s is not supported (yet) as a time dependent parameter\n", ittmp->c_str());
    }

}
#ifndef USECUDA
template <typename TF>
void Boundary<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    const TF no_offset = 0.;

    for (auto& it : tdep_bc)
    {
        it.second->update_time_dependent(sbc.at(it.first).bot, timeloop);
        set_bc<TF>(fields.sp.at(it.first)->fld_bot.data(), fields.sp.at(it.first)->grad_bot.data(), fields.sp.at(it.first)->flux_bot.data(),
                sbc.at(it.first).bcbot, sbc.at(it.first).bot, fields.sp.at(it.first)->visc, no_offset, gd.icells, gd.jcells);
    }
}
#endif

template<typename TF>
void Boundary<TF>::set_values()
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    set_bc<TF>(fields.mp.at("u")->fld_bot.data(), fields.mp.at("u")->grad_bot.data(), fields.mp.at("u")->flux_bot.data(),
           mbcbot, ubot, fields.visc, grid.utrans,
           gd.icells, gd.jcells);
    set_bc<TF>(fields.mp.at("v")->fld_bot.data(), fields.mp.at("v")->grad_bot.data(), fields.mp.at("v")->flux_bot.data(),
           mbcbot, vbot, fields.visc, grid.vtrans,
           gd.icells, gd.jcells);

    set_bc<TF>(fields.mp.at("u")->fld_top.data(), fields.mp.at("u")->grad_top.data(), fields.mp.at("u")->flux_top.data(),
           mbctop, utop, fields.visc, grid.utrans,
           gd.icells, gd.jcells);
    set_bc<TF>(fields.mp.at("v")->fld_top.data(), fields.mp.at("v")->grad_top.data(), fields.mp.at("v")->flux_top.data(),
           mbctop, vtop, fields.visc, grid.vtrans,
           gd.icells, gd.jcells);

    const TF no_offset = 0.;

    for (auto& it : fields.sp)
    {
        // Load 2D fields for bottom boundary from disk.
        if (std::find(sbot_2d_list.begin(), sbot_2d_list.end(), it.first) != sbot_2d_list.end())
        {
            std::string filename = it.first + "_bot.0000000";
            master.print_message("Loading \"%s\" ... ", filename.c_str());

            auto tmp = fields.get_tmp();
            TF* fld_2d_ptr = nullptr;
            if (sbc.at(it.first).bcbot == Boundary_type::Dirichlet_type)
                fld_2d_ptr = it.second->fld_bot.data();
            else if (sbc.at(it.first).bcbot == Boundary_type::Neumann_type)
                fld_2d_ptr = it.second->grad_bot.data();
            else if (sbc.at(it.first).bcbot == Boundary_type::Flux_type)
                fld_2d_ptr = it.second->flux_bot.data();

            if (field3d_io.load_xy_slice(fld_2d_ptr, tmp->fld.data(), filename.c_str()))
            {
                master.print_message("FAILED\n");
                throw std::runtime_error("Error loading 2D field of bottom boundary");
            }
            else
                master.print_message("OK\n");

            fields.release_tmp(tmp);
        }
        else
        {
            set_bc<TF>(it.second->fld_bot.data(), it.second->grad_bot.data(), it.second->flux_bot.data(),
                   sbc.at(it.first).bcbot, sbc.at(it.first).bot, it.second->visc, no_offset,
                   gd.icells, gd.jcells);
            set_bc<TF>(it.second->fld_top.data(), it.second->grad_top.data(), it.second->flux_top.data(),
                   sbc.at(it.first).bctop, sbc.at(it.first).top, it.second->visc, no_offset,
                   gd.icells, gd.jcells);
        }
    }
}

namespace
{
    template<typename TF>
    void calc_ghost_cells_bot_2nd(TF* const restrict a, const TF* const restrict dzh, Boundary_type boundary_type,
                                  TF* const restrict abot, TF* const restrict agradbot,
                                  const int kstart, const int icells, const int jcells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        if (boundary_type == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    a[ijk-kk] = TF(2.)*abot[ij] - a[ijk];
                }
        }
        else if (boundary_type == Boundary_type::Neumann_type || boundary_type == Boundary_type::Flux_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk;
                    a[ijk-kk] = -agradbot[ij]*dzh[kstart] + a[ijk];
                }
        }
    }

    template<typename TF>
    void calc_ghost_cells_top_2nd(TF* const restrict a, const TF* const restrict dzh, Boundary_type boundary_type,
                                  TF* const restrict atop, TF* const restrict agradtop,
                                  const int kend, const int icells, const int jcells, const int ijcells)
    {
        const int jj = icells;
        const int kk = ijcells;

        if (boundary_type == Boundary_type::Dirichlet_type || boundary_type == Boundary_type::Off_type)
        {
            if (boundary_type == Boundary_type::Off_type)
            {
                for (int j=0; j<jcells; ++j)
                    #pragma ivdep
                    for (int i=0; i<icells; ++i)
                    {
                        const int ij  = i + j*jj;
                        const int ijk = i + j*jj + (kend-1)*kk;
                        atop[ij] = TF(3./2.)*a[ijk] - TF(1./2.)*a[ijk-kk];
                    }
            }

            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + (kend-1)*kk;
                    a[ijk+kk] = TF(2.)*atop[ij] - a[ijk];
                }
        }
        else if (boundary_type == Boundary_type::Neumann_type || boundary_type == Boundary_type::Flux_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + (kend-1)*kk;
                    a[ijk+kk] = agradtop[ij]*dzh[kend] + a[ijk];
                }
        }
    }

    template<typename TF>
    void calc_ghost_cells_bot_4th(TF* const restrict a, const TF* const restrict z, Boundary_type boundary_type,
                                  TF* const restrict abot, TF* restrict agradbot,
                                  const int kstart, const int icells, const int jcells, const int ijcells)
    {
        const int jj  = icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        if (boundary_type == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk1;
                    a[ijk-kk1] = TF(8./3.)*abot[ij] - TF(2.)*a[ijk] + TF(1./3.)*a[ijk+kk1];
                    a[ijk-kk2] = TF(8.)*abot[ij] - TF(9.)*a[ijk] + TF(2.)*a[ijk+kk1];
                }
        }
        else if (boundary_type == Boundary_type::Neumann_type || boundary_type == Boundary_type::Flux_type)
        {
            using Finite_difference::O4::grad4;

            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk1;
                    a[ijk-kk1] = TF(-1.)*grad4(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk    ];
                    a[ijk-kk2] = TF(-3.)*grad4(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk+kk1];
                }
        }
    }

    template<typename TF>
    void calc_ghost_cells_top_4th(TF* const restrict a, const TF* const restrict z, Boundary_type boundary_type,
                                  TF* const restrict atop, TF* const restrict agradtop,
                                  const int kend, const int icells, const int jcells, const int ijcells)
    {
        const int jj  = icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        if (boundary_type == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + (kend-1)*kk1;
                    a[ijk+kk1] = TF(8./3.)*atop[ij] - TF(2.)*a[ijk] + TF(1./3.)*a[ijk-kk1];
                    a[ijk+kk2] = TF(8.)*atop[ij] - TF(9.)*a[ijk] + TF(2.)*a[ijk-kk1];
                }
        }
        else if (boundary_type == Boundary_type::Neumann_type || boundary_type == Boundary_type::Flux_type)
        {
            using Finite_difference::O4::grad4;

            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + (kend-1)*kk1;
                    a[ijk+kk1] = TF(1.)*grad4(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk    ];
                    a[ijk+kk2] = TF(3.)*grad4(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk-kk1];
                }
        }
    }

    // BOUNDARY CONDITIONS FOR THE VERTICAL VELOCITY (NO PENETRATION)
    template<typename TF>
    void calc_ghost_cells_botw_cons_4th(TF* const restrict w,
            const int kstart, const int icells, const int jcells, const int ijcells)
    {
        const int jj  = icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ijk = i + j*jj + kstart*kk1;
                w[ijk-kk1] = -w[ijk+kk1];
                w[ijk-kk2] = -w[ijk+kk2];
            }
    }

    template<typename TF>
    void calc_ghost_cells_topw_cons_4th(TF* const restrict w,
            const int kend, const int icells, const int jcells, const int ijcells)
    {
        const int jj  = icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ijk = i + j*jj + kend*kk1;
                w[ijk+kk1] = -w[ijk-kk1];
                w[ijk+kk2] = -w[ijk-kk2];
            }
    }

    template<typename TF>
    void calc_ghost_cells_botw_4th(TF* restrict w,
            const int kstart, const int icells, const int jcells, const int ijcells)
    {
        const int jj  = icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;
        const int kk3 = 3*ijcells;

        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ijk = i + j*jj + kstart*kk1;
                w[ijk-kk1] = TF(-6.)*w[ijk+kk1] + TF(4.)*w[ijk+kk2] - w[ijk+kk3];
            }
    }

    template<typename TF>
    void calc_ghost_cells_topw_4th(TF* restrict w,
            const int kend, const int icells, const int jcells, const int ijcells)
    {
        const int jj  = icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;
        const int kk3 = 3*ijcells;

        for (int j=0; j<jcells; ++j)
            #pragma ivdep
            for (int i=0; i<icells; ++i)
            {
                const int ijk = i + j*jj + kend*kk1;
                w[ijk+kk1] = TF(-6.)*w[ijk-kk1] + TF(4.)*w[ijk-kk2] - w[ijk-kk3];
            }
    }
}

template<typename TF>
void Boundary<TF>::exec(Thermo<TF>& thermo)
{
}

#ifndef USECUDA
template<typename TF>
void Boundary<TF>::set_ghost_cells()
{
    auto& gd = grid.get_grid_data();

    if (grid.get_spatial_order() == Grid_order::Second)
    {
        calc_ghost_cells_bot_2nd<TF>(fields.mp.at("u")->fld.data(), gd.dzh.data(), mbcbot,
                fields.mp.at("u")->fld_bot.data(), fields.mp.at("u")->grad_bot.data(),
                gd.kstart, gd.icells, gd.jcells, gd.ijcells);
        calc_ghost_cells_top_2nd<TF>(fields.mp.at("u")->fld.data(), gd.dzh.data(), mbctop,
                fields.mp.at("u")->fld_top.data(), fields.mp.at("u")->grad_top.data(),
                gd.kend, gd.icells, gd.jcells, gd.ijcells);

        calc_ghost_cells_bot_2nd<TF>(fields.mp.at("v")->fld.data(), gd.dzh.data(), mbcbot,
                fields.mp.at("v")->fld_bot.data(), fields.mp.at("v")->grad_bot.data(),
                gd.kstart, gd.icells, gd.jcells, gd.ijcells);
        calc_ghost_cells_top_2nd<TF>(fields.mp.at("v")->fld.data(), gd.dzh.data(), mbctop,
                fields.mp.at("v")->fld_top.data(), fields.mp.at("v")->grad_top.data(),
                gd.kend, gd.icells, gd.jcells, gd.ijcells);

        for (auto& it : fields.sp)
        {
            calc_ghost_cells_bot_2nd<TF>(it.second->fld.data(), gd.dzh.data(),
                    sbc.at(it.first).bcbot, it.second->fld_bot.data(), it.second->grad_bot.data(),
                    gd.kstart, gd.icells, gd.jcells, gd.ijcells);
            calc_ghost_cells_top_2nd<TF>(it.second->fld.data(), gd.dzh.data(),
                    sbc.at(it.first).bctop, it.second->fld_top.data(), it.second->grad_top.data(),
                    gd.kend, gd.icells, gd.jcells, gd.ijcells);
        }
    }
    else if (grid.get_spatial_order() == Grid_order::Fourth)
    {
        calc_ghost_cells_bot_4th<TF>(fields.mp.at("u")->fld.data(), gd.z.data(), mbcbot,
                fields.mp.at("u")->fld_bot.data(), fields.mp.at("u")->grad_bot.data(),
                gd.kstart, gd.icells, gd.jcells, gd.ijcells);
        calc_ghost_cells_top_4th<TF>(fields.mp.at("u")->fld.data(), gd.z.data(), mbctop,
                fields.mp.at("u")->fld_top.data(), fields.mp.at("u")->grad_top.data(),
                gd.kend, gd.icells, gd.jcells, gd.ijcells);

        calc_ghost_cells_bot_4th<TF>(fields.mp.at("v")->fld.data(), gd.z.data(), mbcbot,
                fields.mp.at("v")->fld_bot.data(), fields.mp.at("v")->grad_bot.data(),
                gd.kstart, gd.icells, gd.jcells, gd.ijcells);
        calc_ghost_cells_top_4th<TF>(fields.mp.at("v")->fld.data(), gd.z.data(), mbctop,
                fields.mp.at("v")->fld_top.data(), fields.mp.at("v")->grad_top.data(),
                gd.kend, gd.icells, gd.jcells, gd.ijcells);

        calc_ghost_cells_botw_4th<TF>(fields.mp.at("w")->fld.data(),
                gd.kstart, gd.icells, gd.jcells, gd.ijcells);
        calc_ghost_cells_topw_4th<TF>(fields.mp.at("w")->fld.data(),
                gd.kend, gd.icells, gd.jcells, gd.ijcells);

        for (auto& it : fields.sp)
        {
            calc_ghost_cells_bot_4th<TF>(it.second->fld.data(), gd.z.data(), sbc.at(it.first).bcbot,
                    it.second->fld_bot.data(), it.second->grad_bot.data(),
                    gd.kstart, gd.icells, gd.jcells, gd.ijcells);
            calc_ghost_cells_top_4th<TF>(it.second->fld.data(), gd.z.data(), sbc.at(it.first).bctop,
                    it.second->fld_top.data(), it.second->grad_top.data(),
                    gd.kend, gd.icells, gd.jcells, gd.ijcells);
        }
    }

    // Update the boundary fields that are a slave of the boundary condition.
    update_slave_bcs();
}

template<typename TF>
void Boundary<TF>::set_ghost_cells_w(const Boundary_w_type boundary_w_type)
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    if (grid.get_spatial_order() == Grid_order::Fourth)
    {
        if (boundary_w_type == Boundary_w_type::Normal_type)
        {
            calc_ghost_cells_botw_4th<TF>(fields.mp.at("w")->fld.data(),
                    gd.kstart, gd.icells, gd.jcells, gd.ijcells);
            calc_ghost_cells_topw_4th<TF>(fields.mp.at("w")->fld.data(),
                    gd.kend, gd.icells, gd.jcells, gd.ijcells);
        }
        else if (boundary_w_type == Boundary_w_type::Conservation_type)
        {
            calc_ghost_cells_botw_cons_4th<TF>(fields.mp.at("w")->fld.data(),
                    gd.kstart, gd.icells, gd.jcells, gd.ijcells);
            calc_ghost_cells_topw_cons_4th<TF>(fields.mp.at("w")->fld.data(),
                    gd.kend, gd.icells, gd.jcells, gd.ijcells);
        }
    }
}
#endif

/*
template<typename TF>
void Boundary<TF>::exec_cross()
{
}
*/

template<typename TF>
void Boundary<TF>::exec_stats(Stats<TF>&)
{
}

template<typename TF>
void Boundary<TF>::exec_column(Column<TF>&)
{
}

// Computational kernel for boundary calculation.
namespace
{
    template<typename TF, int spatial_order>
    void calc_slave_bc_bot(TF* const restrict abot, TF* const restrict agradbot, TF* const restrict afluxbot,
                           const TF* const restrict a, const TF* const restrict dzhi,
                           const Boundary_type boundary_type, const TF visc,
                           const int kstart, const int icells, const int jcells, const int ijcells)
    {
        const int jj = icells;
        const int kk1 = 1*ijcells;
        const int kk2 = 2*ijcells;

        using namespace Finite_difference;

        // Variable dzhi in this case is dzhi for 2nd order and dzhi4 for 4th order.
        if (boundary_type == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk1;
                    if (spatial_order == 2)
                        agradbot[ij] = O2::grad2(a[ijk-kk1], a[ijk]) * dzhi[kstart];
                    else if (spatial_order == 4)
                        agradbot[ij] = O4::grad4(a[ijk-kk2], a[ijk-kk1], a[ijk], a[ijk+kk1]) * dzhi[kstart];
                    afluxbot[ij] = -visc*agradbot[ij];
                }
        }
        else if (boundary_type == Boundary_type::Neumann_type || boundary_type == Boundary_type::Flux_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk1;
                    if (spatial_order == 2)
                        abot[ij] = O2::interp2(a[ijk-kk1], a[ijk]);
                    else if (spatial_order == 4)
                        abot[ij] = O4::interp4c(a[ijk-kk2], a[ijk-kk1], a[ijk], a[ijk+kk1]);
                }
        }
    }
}

template<typename TF>
void Boundary<TF>::update_slave_bcs()
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    if (grid.get_spatial_order() == Grid_order::Second)
    {
        calc_slave_bc_bot<TF,2>(fields.mp.at("u")->fld_bot.data(), fields.mp.at("u")->grad_bot.data(), fields.mp.at("u")->flux_bot.data(),
                                fields.mp.at("u")->fld.data(), gd.dzhi.data(),
                                mbcbot, fields.mp.at("u")->visc,
                                gd.kstart, gd.icells, gd.jcells, gd.ijcells);

        calc_slave_bc_bot<TF,2>(fields.mp.at("v")->fld_bot.data(), fields.mp.at("v")->grad_bot.data(), fields.mp.at("v")->flux_bot.data(),
                                fields.mp.at("v")->fld.data(), gd.dzhi.data(),
                                mbcbot, fields.mp.at("v")->visc,
                                gd.kstart, gd.icells, gd.jcells, gd.ijcells);

        for (auto& it : fields.sp)
            calc_slave_bc_bot<TF,2>(it.second->fld_bot.data(), it.second->grad_bot.data(), it.second->flux_bot.data(),
                                    it.second->fld.data(), gd.dzhi.data(),
                                    sbc.at(it.first).bcbot, it.second->visc,
                                    gd.kstart, gd.icells, gd.jcells, gd.ijcells);
    }
    else if (grid.get_spatial_order() == Grid_order::Fourth)
    {
        calc_slave_bc_bot<TF,4>(fields.mp.at("u")->fld_bot.data(), fields.mp.at("u")->grad_bot.data(), fields.mp.at("u")->flux_bot.data(),
                                fields.mp.at("u")->fld.data(), gd.dzhi4.data(),
                                mbcbot, fields.mp.at("u")->visc,
                                gd.kstart, gd.icells, gd.jcells, gd.ijcells);

        calc_slave_bc_bot<TF,4>(fields.mp.at("v")->fld_bot.data(), fields.mp.at("v")->grad_bot.data(), fields.mp.at("v")->flux_bot.data(),
                                fields.mp.at("v")->fld.data(), gd.dzhi4.data(),
                                mbcbot, fields.mp.at("v")->visc,
                                gd.kstart, gd.icells, gd.jcells, gd.ijcells);

        for (auto& it : fields.sp)
            calc_slave_bc_bot<TF,4>(it.second->fld_bot.data(), it.second->grad_bot.data(), it.second->flux_bot.data(),
                                    it.second->fld.data(), gd.dzhi4.data(),
                                    sbc.at(it.first).bcbot, it.second->visc,
                                    gd.kstart, gd.icells, gd.jcells, gd.ijcells);
    }
}

template<typename TF>
const std::vector<TF>& Boundary<TF>::get_z0m() const
{
    throw std::runtime_error("Function get_z0m() not implemented in base boundary.");
}

template<typename TF>
const std::vector<TF>& Boundary<TF>::get_dudz() const
{
    throw std::runtime_error("Function get_dudz() not implemented in base boundary.");
}

template<typename TF>
const std::vector<TF>& Boundary<TF>::get_dvdz() const
{
    throw std::runtime_error("Function get_dvdz() not implemented in base boundary.");
}

template<typename TF>
const std::vector<TF>& Boundary<TF>::get_dbdz() const
{
    throw std::runtime_error("Function get_dbdz() not implemented in base boundary.");
}

template<typename TF>
std::shared_ptr<Boundary<TF>> Boundary<TF>::factory(
        Master& master, Grid<TF>& grid, Fields<TF>& fields, Input& input)
{
    std::string swboundary;
    swboundary = input.get_item<std::string>("boundary", "swboundary", "", "default");

    if (swboundary == "default")
        return std::make_shared<Boundary<TF>>(master, grid, fields, input);
    else if (swboundary == "surface")
        return std::make_shared<Boundary_surface<TF>>(master, grid, fields, input);
    else if (swboundary == "surface_bulk")
        return std::make_shared<Boundary_surface_bulk<TF>>(master, grid, fields, input);
    else if (swboundary == "surface_lsm")
        return std::make_shared<Boundary_surface_lsm<TF>>(master, grid, fields, input);
    else
    {
        std::string msg = swboundary + " is an illegal value for swboundary";
        throw std::runtime_error(msg);
    }
}
/*
template<typename TF>
void Boundary<TF>::get_mask(Field3d* field, Field3d* fieldh, Mask* m)
{
    // Set surface mask
    for (int i=0; i<grid.ijcells; ++i)
        fieldh->fld_bot[i] = 1;

    // Set atmospheric mask
    for (int i=0; i<grid.ncells; ++i)
    {
        field ->data[i] = 1;
        fieldh->data[i] = 1;
    }
}

template<typename TF>
void Boundary<TF>::get_surface_mask(Field3d* field)
{
    // Set surface mask
    for (int i=0; i<grid.ijcells; ++i)
        field->fld_bot[i] = 1;
}

*/
template<typename TF>
void Boundary<TF>::prepare_device()
{
}

template<typename TF>
void Boundary<TF>::forward_device()
{
}

template<typename TF>
void Boundary<TF>::backward_device()
{
}

template class Boundary<double>;
template class Boundary<float>;
