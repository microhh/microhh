/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
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
#include "defines.h"
#include "timeloop.h"
#include "finite_difference.h"

#include "boundary_cyclic.h"

// Boundary schemes.
#include "boundary.h"
// #include "boundary_surface.h"
// #include "boundary_surface_bulk.h"
// #include "boundary_surface_patch.h"
// #include "boundary_patch.h"

template<typename TF>
Boundary<TF>::Boundary(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin),
    grid(gridin),
    fields(fieldsin),
    boundary_cyclic(master, grid)
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
    int nerror = 0;

    std::string swbot = input.get_item<std::string>("boundary", "mbcbot", "");
    std::string swtop = input.get_item<std::string>("boundary", "mbctop", "");

    ubot = input.get_item<TF>("boundary", "ubot", "", 0.);
    utop = input.get_item<TF>("boundary", "utop", "", 0.);
    vbot = input.get_item<TF>("boundary", "vbot", "", 0.);
    vtop = input.get_item<TF>("boundary", "vtop", "", 0.);

    // set the bottom bc
    if (swbot == "noslip")
        mbcbot = Boundary_type::Dirichlet_type;
    else if (swbot == "freeslip")
        mbcbot = Boundary_type::Neumann_type;
    else if (swbot == "ustar")
        mbcbot = Boundary_type::Ustar_type;
    else
    {
        master.print_error("%s is illegal value for mbcbot\n", swbot.c_str());
        nerror++;
    }

    // set the top bc
    if (swtop == "noslip")
        mbctop = Boundary_type::Dirichlet_type;
    else if (swtop == "freeslip")
        mbctop = Boundary_type::Neumann_type;
    else if (swtop == "ustar")
        mbctop = Boundary_type::Ustar_type;
    else
    {
        master.print_error("%s is illegal value for mbctop\n", swtop.c_str());
        nerror++;
    }

    // read the boundaries per field
    for (auto& it : fields.sp)
    {
        sbc.emplace(it.first, Field3dBc());
        swbot = input.get_item<std::string>("boundary", "sbcbot", it.first);
        swtop = input.get_item<std::string>("boundary", "sbctop", it.first);
        sbc.at(it.first).bot = input.get_item<double>("boundary", "sbot", it.first);
        sbc.at(it.first).top = input.get_item<double>("boundary", "stop", it.first);

        // set the bottom bc
        if (swbot == "dirichlet")
            sbc.at(it.first).bcbot = Boundary_type::Dirichlet_type;
        else if (swbot == "neumann")
            sbc.at(it.first).bcbot = Boundary_type::Neumann_type;
        else if (swbot == "flux")
            sbc.at(it.first).bcbot = Boundary_type::Flux_type;
        else
        {
            master.print_error("%s is illegal value for sbcbot\n", swbot.c_str());
            nerror++;
        }

        // set the top bc
        if (swtop == "dirichlet")
            sbc.at(it.first).bctop = Boundary_type::Dirichlet_type;
        else if (swtop == "neumann")
            sbc.at(it.first).bctop = Boundary_type::Neumann_type;
        else if (swtop == "flux")
            sbc.at(it.first).bctop = Boundary_type::Flux_type;
        else
        {
            master.print_error("%s is illegal value for sbctop\n", swtop.c_str());
            nerror++;
        }
    }

    // get the list of time varying variables
    swtimedep   = input.get_item<std::string>("boundary", "swtimedep"  , "", "0");
    timedeplist = input.get_list<std::string>("boundary", "timedeplist", "", std::vector<std::string>());

    if (nerror)
        throw 1;
}

template<typename TF>
void Boundary<TF>::init(Input& input)
{
    // Read the boundary information from the ini files, it throws at error.
    process_bcs(input);

    int nerror = 0;

    // there is no option (yet) for prescribing ustar without surface model
    if (mbcbot == Boundary_type::Ustar_type || mbctop == Boundary_type::Ustar_type)
    {
        master.print_error("ustar bc is not supported for default boundary\n");
        ++nerror;
    }

    if (nerror)
        throw std::runtime_error("Cannot use ustar bc for default boundary");

    // Initialize the boundary cyclic.
    boundary_cyclic.init();
}

template<typename TF>
void Boundary<TF>::create(Input& input)
{
    // process_time_dependent(input);
}

/*
template<typename TF>
void Boundary<TF>::process_time_dependent(Input* inputin)
{
    int nerror = 0;

    if (swtimedep == "1")
    {
        // create temporary list to check which entries are used
        std::vector<std::string> tmplist = timedeplist;

        // see if there is data available for the surface boundary conditions
        for (auto& it : fields.sp)
        {
            std::string name = "sbot[" + it.first + "]";
            if (std::find(timedeplist.begin(), timedeplist.end(), name) != timedeplist.end())
            {
                nerror += inputin->get_time(&timedepdata[name], &timedeptime, name);

                // remove the item from the tmplist
                std::vector<std::string>::iterator ittmp = std::find(tmplist.begin(), tmplist.end(), name);
                if (ittmp != tmplist.end())
                    tmplist.erase(ittmp);
            }
        }

        // display a warning for the non-supported
        for (std::vector<std::string>::const_iterator ittmp=tmplist.begin(); ittmp!=tmplist.end(); ++ittmp)
            master.print_warning("%s is not supported (yet) as a time dependent parameter\n", ittmp->c_str());
    }

    if (nerror)
        throw 1;
}

template<typename TF>
void Boundary<TF>::update_time_dependent()
{
    if (swtimedep == "0")
        return;

    // first find the index for the time entries
    unsigned int index0 = 0;
    unsigned int index1 = 0;
    for (std::vector<double>::const_iterator it=timedeptime.begin(); it!=timedeptime.end(); ++it)
    {
        if (model->timeloop->get_time() < *it)
            break;
        else
            ++index1;
    }

    // second, calculate the weighting factor
    double fac0, fac1;

    // correct for out of range situations where the simulation is longer than the time range in input
    if (index1 == 0)
    {
        fac0 = 0.;
        fac1 = 1.;
        index0 = 0;
    }
    else if (index1 == timedeptime.size())
    {
        fac0 = 1.;
        fac1 = 0.;
        index0 = index1-1;
        index1 = index0;
    }
    else
    {
        index0 = index1-1;
        double timestep;
        timestep = timedeptime[index1] - timedeptime[index0];
        fac0 = (timedeptime[index1] - model->timeloop->get_time()) / timestep;
        fac1 = (model->timeloop->get_time() - timedeptime[index0]) / timestep;
    }

    // process time dependent bcs for the surface fluxes
    for (FieldMap::const_iterator it1=fields.sp.begin(); it1!=fields.sp.end(); ++it1)
    {
        std::string name = "sbot[" + it1->first + "]";
        std::map<std::string, double *>::const_iterator it2 = timedepdata.find(name);
        if (it2 != timedepdata.end())
        {
            sbc[it1->first]->bot = fac0*it2->second[index0] + fac1*it2->second[index1];

            // BvS: for now branched here; seems a bit wasteful to copy the entire settimedep to boundary.cu?
            const double noOffset = 0.;

#ifndef USECUDA
            set_bc(it1->second->fld_bot, it1->second->grad_bot, it1->second->flux_bot, sbc[it1->first]->bcbot, sbc[it1->first]->bot, it1->second->visc, noOffset);
#else
            set_bc_g(it1->second->fld_bot_g, it1->second->grad_bot_g, it1->second->flux_bot_g, sbc[it1->first]->bcbot, sbc[it1->first]->bot, it1->second->visc, noOffset);
#endif
        }
    }
}
*/

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

    const double no_offset = 0.;

    for (auto& it : fields.sp)
    {
        set_bc<TF>(it.second->fld_bot.data(), it.second->grad_bot.data(), it.second->flux_bot.data(),
               sbc.at(it.first).bcbot, sbc.at(it.first).bot, it.second->visc, no_offset,
               gd.icells, gd.jcells);
        set_bc<TF>(it.second->fld_top.data(), it.second->grad_top.data(), it.second->flux_top.data(),
               sbc.at(it.first).bctop, sbc.at(it.first).top, it.second->visc, no_offset,
               gd.icells, gd.jcells);
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
                    a[ijk-kk] = 2.*abot[ij] - a[ijk];
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

        if (boundary_type == Boundary_type::Dirichlet_type)
        {
            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + (kend-1)*kk;
                    a[ijk+kk] = 2.*atop[ij] - a[ijk];
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
                    a[ijk-kk1] = (8./3.)*abot[ij] - 2.*a[ijk] + (1./3.)*a[ijk+kk1];
                    a[ijk-kk2] = 8.*abot[ij] - 9.*a[ijk] + 2.*a[ijk+kk1];
                }
        }
        else if (boundary_type == Boundary_type::Neumann_type || boundary_type == Boundary_type::Flux_type)
        {
            using Finite_difference::O4::grad4x;

            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + kstart*kk1;
                    a[ijk-kk1] = -(1./24.)*grad4x(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk    ];
                    a[ijk-kk2] = -(1./ 8.)*grad4x(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk+kk1];
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
                    a[ijk+kk1] = (8./3.)*atop[ij] - 2.*a[ijk] + (1./3.)*a[ijk-kk1];
                    a[ijk+kk2] = 8.*atop[ij] - 9.*a[ijk] + 2.*a[ijk-kk1];
                }
        }
        else if (boundary_type == Boundary_type::Neumann_type || boundary_type == Boundary_type::Flux_type)
        {
            using Finite_difference::O4::grad4x;

            for (int j=0; j<jcells; ++j)
                #pragma ivdep
                for (int i=0; i<icells; ++i)
                {
                    const int ij  = i + j*jj;
                    const int ijk = i + j*jj + (kend-1)*kk1;
                    a[ijk+kk1] = (1./24.)*grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk    ];
                    a[ijk+kk2] = (1./ 8.)*grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk-kk1];
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
                w[ijk-kk1] = -6.*w[ijk+kk1] + 4.*w[ijk+kk2] - w[ijk+kk3];
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
                w[ijk+kk1] = -6.*w[ijk-kk1] + 4.*w[ijk-kk2] - w[ijk-kk3];
            }
    }
}

#ifndef USECUDA
template<typename TF>
void Boundary<TF>::exec()
{
    /* CvH TEMPORARY DISABLED FOR TESTING.
    // Cyclic boundary conditions, do this before the bottom BC's
    grid.boundary_cyclic(fields.mp.at("u")->fld.data());
    grid.boundary_cyclic(fields.mp.at("v")->fld.data());
    grid.boundary_cyclic(fields.mp.at("w")->fld.data());

    for (auto& it : fields.sp)
        grid.boundary_cyclic(it.second->fld.data());
    CvH END DISABLED */

    boundary_cyclic.exec(fields.mp.at("u")->fld.data());
    boundary_cyclic.exec(fields.mp.at("v")->fld.data());
    boundary_cyclic.exec(fields.mp.at("w")->fld.data());

    for (auto& it : fields.sp)
        boundary_cyclic.exec(it.second->fld.data());

    // Update the boundary values.
    update_bcs();

    const Grid_data<TF>& gd = grid.get_grid_data();

    if (grid.swspatialorder == "2")
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
    else if (grid.swspatialorder == "4")
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

    if (grid.swspatialorder == "4")
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

template<typename TF>
void Boundary<TF>::exec_stats(Mask* m)
{
}
*/

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
                        agradbot[ij] = O2::grad2x(a[ijk-kk1], a[ijk]) * dzhi[kstart];
                    else if (spatial_order == 4)
                        agradbot[ij] = O4::grad4x(a[ijk-kk2], a[ijk-kk1], a[ijk], a[ijk+kk1]) * dzhi[kstart];
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
                        abot[ij] = O4::interp4(a[ijk-kk2], a[ijk-kk1], a[ijk], a[ijk+kk1]);
                }
        }
    }
}

template<typename TF>
void Boundary<TF>::update_slave_bcs()
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    if (grid.swspatialorder == "2")
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
    else if (grid.swspatialorder == "4")
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
void Boundary<TF>::update_bcs()
{
}

template<typename TF>
std::shared_ptr<Boundary<TF>> Boundary<TF>::factory(Master& master, Grid<TF>& grid, Fields<TF>& fields, Input& input)
{
    std::string swboundary;
    swboundary = input.get_item<std::string>("boundary", "swboundary", "", "default");

    // if (swboundary == "surface")
    //     return new Boundary_surface(modelin, inputin);
    // if (swboundary == "surface_bulk")
    //     return new Boundary_surface_bulk(modelin, inputin);
    // else if (swboundary == "surface_patch")
    //     return new Boundary_surface_patch(modelin, inputin);
    // else if (swboundary == "patch")
    //     return new Boundary_patch(modelin, inputin);
    // else if (swboundary == "default")
    if (swboundary == "default")
        return std::make_shared<Boundary<TF>>(master, grid, fields, input);
    else
    {
        master.print_error("\"%s\" is an illegal value for swboundary\n", swboundary.c_str());
        throw std::runtime_error("Illegal value for swboundary");
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
*/

template class Boundary<double>;
template class Boundary<float>;
