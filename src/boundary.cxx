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
#include <algorithm>
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "model.h"
#include "timeloop.h"
#include "finite_difference.h"

// Boundary schemes.
#include "boundary.h"
#include "boundary_surface.h"
#include "boundary_user.h"

Boundary::Boundary(Model* modelin, Input* inputin)
{
    swboundary = "default";

    model  = modelin;
    grid   = model->grid;
    fields = model->fields;
    master = model->master;
}

Boundary::~Boundary()
{
    for (BcMap::const_iterator it=sbc.begin(); it!=sbc.end(); ++it)
        delete it->second;

    // empty the map
    sbc.clear();

    // clean up time dependent data
    for (std::map<std::string, double *>::const_iterator it=timedepdata.begin(); it!=timedepdata.end(); ++it)
        delete[] it->second;
}

std::string Boundary::get_switch()
{
    return swboundary;
}

void Boundary::process_bcs(Input* inputin)
{
    int nerror = 0;

    std::string swbot, swtop;

    nerror += inputin->get_item(&swbot, "boundary", "mbcbot", "");
    nerror += inputin->get_item(&swtop, "boundary", "mbctop", "");

    nerror += inputin->get_item(&ubot, "boundary", "ubot", "", 0.);
    nerror += inputin->get_item(&utop, "boundary", "utop", "", 0.);
    nerror += inputin->get_item(&vbot, "boundary", "vbot", "", 0.);
    nerror += inputin->get_item(&vtop, "boundary", "vtop", "", 0.);

    // set the bottom bc
    if (swbot == "noslip")
        mbcbot = Dirichlet_type;
    else if (swbot == "freeslip")
        mbcbot = Neumann_type;
    else if (swbot == "ustar")
        mbcbot = Ustar_type;
    else
    {
        master->print_error("%s is illegal value for mbcbot\n", swbot.c_str());
        nerror++;
    }

    // set the top bc
    if (swtop == "noslip")
        mbctop = Dirichlet_type;
    else if (swtop == "freeslip")
        mbctop = Neumann_type;
    else if (swtop == "ustar")
        mbctop = Ustar_type;
    else
    {
        master->print_error("%s is illegal value for mbctop\n", swtop.c_str());
        nerror++;
    }

    // read the boundaries per field
    for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
        sbc[it->first] = new Field3dBc;
        nerror += inputin->get_item(&swbot, "boundary", "sbcbot", it->first);
        nerror += inputin->get_item(&swtop, "boundary", "sbctop", it->first);
        nerror += inputin->get_item(&sbc[it->first]->bot, "boundary", "sbot", it->first);
        nerror += inputin->get_item(&sbc[it->first]->top, "boundary", "stop", it->first);

        // set the bottom bc
        if (swbot == "dirichlet")
            sbc[it->first]->bcbot = Dirichlet_type;
        else if (swbot == "neumann")
            sbc[it->first]->bcbot = Neumann_type;
        else if (swbot == "flux")
            sbc[it->first]->bcbot = Flux_type;
        else
        {
            master->print_error("%s is illegal value for sbcbot\n", swbot.c_str());
            nerror++;
        }

        // set the top bc
        if (swtop == "dirichlet")
            sbc[it->first]->bctop = Dirichlet_type;
        else if (swtop == "neumann")
            sbc[it->first]->bctop = Neumann_type;
        else if (swtop == "flux")
            sbc[it->first]->bctop = Flux_type;
        else
        {
            master->print_error("%s is illegal value for sbctop\n", swtop.c_str());
            nerror++;
        }
    }

    // get the list of time varying variables
    nerror += inputin->get_item(&swtimedep  , "boundary", "swtimedep"  , "", "0");
    nerror += inputin->get_list(&timedeplist, "boundary", "timedeplist", "");

    if (nerror)
        throw 1;
}

void Boundary::init(Input *inputin)
{
    // Read the boundary information from the ini files, it throws at error.
    process_bcs(inputin);

    int nerror = 0;

    // there is no option (yet) for prescribing ustar without surface model
    if (mbcbot == Ustar_type || mbctop == Ustar_type)
    {
        master->print_error("ustar bc is not supported for default boundary\n");
        ++nerror;
    }

    if (nerror)
        throw 1;
}

void Boundary::create(Input* inputin)
{
    process_time_dependent(inputin);
}

void Boundary::process_time_dependent(Input* inputin)
{
    int nerror = 0;

    if (swtimedep == "1")
    {
        // create temporary list to check which entries are used
        std::vector<std::string> tmplist = timedeplist;

        // see if there is data available for the surface boundary conditions
        for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
        {
            std::string name = "sbot[" + it->first + "]";
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
            master->print_warning("%s is not supported (yet) as a time dependent parameter\n", ittmp->c_str());
    }

    if (nerror)
        throw 1;
}

void Boundary::update_time_dependent()
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
    for (FieldMap::const_iterator it1=fields->sp.begin(); it1!=fields->sp.end(); ++it1)
    {
        std::string name = "sbot[" + it1->first + "]";
        std::map<std::string, double *>::const_iterator it2 = timedepdata.find(name);
        if (it2 != timedepdata.end())
        {
            sbc[it1->first]->bot = fac0*it2->second[index0] + fac1*it2->second[index1];

            // BvS: for now branched here; seems a bit wasteful to copy the entire settimedep to boundary.cu?
            const double noOffset = 0.;

#ifndef USECUDA
            set_bc(it1->second->databot, it1->second->datagradbot, it1->second->datafluxbot, sbc[it1->first]->bcbot, sbc[it1->first]->bot, it1->second->visc, noOffset);
#else
            set_bc_g(it1->second->databot_g, it1->second->datagradbot_g, it1->second->datafluxbot_g, sbc[it1->first]->bcbot, sbc[it1->first]->bot, it1->second->visc, noOffset);
#endif
        }
    }
}

void Boundary::set_values()
{
    const double noOffset = 0.;

    set_bc(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot, mbcbot, ubot, fields->visc, grid->utrans);
    set_bc(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot, mbcbot, vbot, fields->visc, grid->vtrans);

    set_bc(fields->u->datatop, fields->u->datagradtop, fields->u->datafluxtop, mbctop, utop, fields->visc, grid->utrans);
    set_bc(fields->v->datatop, fields->v->datagradtop, fields->v->datafluxtop, mbctop, vtop, fields->visc, grid->vtrans);

    for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
        set_bc(it->second->databot, it->second->datagradbot, it->second->datafluxbot, sbc[it->first]->bcbot, sbc[it->first]->bot, it->second->visc, noOffset);
        set_bc(it->second->datatop, it->second->datagradtop, it->second->datafluxtop, sbc[it->first]->bctop, sbc[it->first]->top, it->second->visc, noOffset);
    }
}

#ifndef USECUDA
void Boundary::exec()
{
    // Cyclic boundary conditions, do this before the bottom BC's
    grid->boundary_cyclic(fields->u->data);
    grid->boundary_cyclic(fields->v->data);
    grid->boundary_cyclic(fields->w->data);

    for (FieldMap::const_iterator it = fields->sp.begin(); it!=fields->sp.end(); ++it)
        grid->boundary_cyclic(it->second->data);

    // Update the boundary values.
    update_bcs();

    if (grid->swspatialorder == "2")
    {
        calc_ghost_cells_bot_2nd(fields->u->data, grid->dzh, mbcbot, fields->u->databot, fields->u->datagradbot);
        calc_ghost_cells_top_2nd(fields->u->data, grid->dzh, mbctop, fields->u->datatop, fields->u->datagradtop);

        calc_ghost_cells_bot_2nd(fields->v->data, grid->dzh, mbcbot, fields->v->databot, fields->v->datagradbot);
        calc_ghost_cells_top_2nd(fields->v->data, grid->dzh, mbctop, fields->v->datatop, fields->v->datagradtop);

        for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
        {
            calc_ghost_cells_bot_2nd(it->second->data, grid->dzh, sbc[it->first]->bcbot, it->second->databot, it->second->datagradbot);
            calc_ghost_cells_top_2nd(it->second->data, grid->dzh, sbc[it->first]->bctop, it->second->datatop, it->second->datagradtop);
        }
    }
    else if (grid->swspatialorder == "4")
    {
        calc_ghost_cells_bot_4th(fields->u->data, grid->z, mbcbot, fields->u->databot, fields->u->datagradbot);
        calc_ghost_cells_top_4th(fields->u->data, grid->z, mbctop, fields->u->datatop, fields->u->datagradtop);

        calc_ghost_cells_bot_4th(fields->v->data, grid->z, mbcbot, fields->v->databot, fields->v->datagradbot);
        calc_ghost_cells_top_4th(fields->v->data, grid->z, mbctop, fields->v->datatop, fields->v->datagradtop);

        calc_ghost_cells_botw_4th(fields->w->data);
        calc_ghost_cells_topw_4th(fields->w->data);

        for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
        {
            calc_ghost_cells_bot_4th(it->second->data, grid->z, sbc[it->first]->bcbot, it->second->databot, it->second->datagradbot);
            calc_ghost_cells_top_4th(it->second->data, grid->z, sbc[it->first]->bctop, it->second->datatop, it->second->datagradtop);
        }
    }

    // Update the boundary fields that are a slave of the boundary condition.
    update_slave_bcs();
}

void Boundary::set_ghost_cells_w(const Boundary_w_type boundary_w_type)
{
    if (grid->swspatialorder == "4")
    {
        if (boundary_w_type == Normal_type)
        {
            calc_ghost_cells_botw_4th(fields->w->data);
            calc_ghost_cells_topw_4th(fields->w->data);
        }
        else if (boundary_w_type == Conservation_type)
        {
            calc_ghost_cells_botw_cons_4th(fields->w->data);
            calc_ghost_cells_topw_cons_4th(fields->w->data);
        }
    }
}
#endif

void Boundary::exec_cross()
{
}

void Boundary::exec_stats(Mask* m)
{
}

// Computational kernel for boundary calculation.
namespace
{
    template<int spatial_order>
    void calc_slave_bc_bot(double* const restrict abot, double* const restrict agradbot, double* const restrict afluxbot,
                           const double* const restrict a,
                           const Grid* const grid, const double* const restrict dzhi,
                           const Boundary::Boundary_type boundary_type, const double visc)
    {
        const int jj = grid->icells;
        const int kk1 = 1*grid->ijcells;
        const int kk2 = 2*grid->ijcells;

        const int kstart = grid->kstart;

        using namespace Finite_difference;

        // Variable dzhi in this case is dzhi for 2nd order and dzhi4 for 4th order.
        if (boundary_type == Boundary::Dirichlet_type)
        {
            for (int j=0; j<grid->jcells; ++j)
                #pragma ivdep
                for (int i=0; i<grid->icells; ++i)
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
        else if (boundary_type == Boundary::Neumann_type || boundary_type == Boundary::Flux_type)
        {
            for (int j=0; j<grid->jcells; ++j)
                #pragma ivdep
                for (int i=0; i<grid->icells; ++i)
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

void Boundary::update_slave_bcs()
{
    if (grid->swspatialorder == "2")
    {
        calc_slave_bc_bot<2>(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot,
                             fields->u->data,
                             grid, grid->dzhi,
                             mbcbot, fields->u->visc);

        calc_slave_bc_bot<2>(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot,
                             fields->v->data,
                             grid, grid->dzhi,
                             mbcbot, fields->v->visc);

        for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
            calc_slave_bc_bot<2>(it->second->databot, it->second->datagradbot, it->second->datafluxbot,
                                 it->second->data,
                                 grid, grid->dzhi,
                                 sbc[it->first]->bcbot, it->second->visc);
    }
    else if (grid->swspatialorder == "4")
    {
        calc_slave_bc_bot<4>(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot,
                             fields->u->data,
                             grid, grid->dzhi4,
                             mbcbot, fields->u->visc);

        calc_slave_bc_bot<4>(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot,
                             fields->v->data,
                             grid, grid->dzhi4,
                             mbcbot, fields->v->visc);

        for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
            calc_slave_bc_bot<4>(it->second->databot, it->second->datagradbot, it->second->datafluxbot,
                                 it->second->data,
                                 grid, grid->dzhi4,
                                 sbc[it->first]->bcbot, it->second->visc);
    }
}

void Boundary::update_bcs()
{
}

Boundary* Boundary::factory(Master* masterin, Input* inputin, Model* modelin)
{
    std::string swboundary;
    if (inputin->get_item(&swboundary, "boundary", "swboundary", "", "default"))
        return 0;

    if (swboundary == "surface")
        return new Boundary_surface(modelin, inputin);
    else if (swboundary == "user")
        return new Boundary_user(modelin, inputin);
    else if (swboundary == "default")
        return new Boundary(modelin, inputin);
    else
    {
        masterin->print_error("\"%s\" is an illegal value for swboundary\n", swboundary.c_str());
        throw 1;
    }
}

void Boundary::set_bc(double* restrict a, double* restrict agrad, double* restrict aflux, Boundary_type sw, double aval, double visc, double offset)
{
    int ij,jj;
    jj = grid->icells;

    if (sw == Dirichlet_type)
    {
        for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                ij = i + j*jj;
                a[ij] = aval - offset;
            }
    }
    else if (sw == Neumann_type)
    {
        for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                ij = i + j*jj;
                agrad[ij] = aval;
                aflux[ij] = -aval*visc;
            }
    }
    else if (sw == Flux_type)
    {
        for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                ij = i + j*jj;
                aflux[ij] = aval;
                agrad[ij] = -aval/visc;
            }
    }
}

// BOUNDARY CONDITIONS THAT CONTAIN A 2D PATTERN
void Boundary::calc_ghost_cells_bot_2nd(double* restrict a, double* restrict dzh, Boundary_type boundary_type,
                                        double* restrict abot, double* restrict agradbot)
{
    int ij,ijk,jj,kk,kstart;

    jj = grid->icells;
    kk = grid->ijcells;

    kstart = grid->kstart;

    if (boundary_type == Dirichlet_type)
    {
        for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                ij  = i + j*jj;
                ijk = i + j*jj + kstart*kk;
                a[ijk-kk] = 2.*abot[ij] - a[ijk];
            }
    }
    else if (boundary_type == Neumann_type || boundary_type == Flux_type)
    {
        for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                ij  = i + j*jj;
                ijk = i + j*jj + kstart*kk;
                a[ijk-kk] = -agradbot[ij]*dzh[kstart] + a[ijk];
            }
    }
}

void Boundary::calc_ghost_cells_top_2nd(double* restrict a, double* restrict dzh, Boundary_type boundary_type,
                                        double* restrict atop, double* restrict agradtop)
{
    int ij,ijk,jj,kk,kend;

    kend = grid->kend;

    jj = grid->icells;
    kk = grid->ijcells;

    if (boundary_type == Dirichlet_type)
    {
        for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                ij  = i + j*jj;
                ijk = i + j*jj + (kend-1)*kk;
                a[ijk+kk] = 2.*atop[ij] - a[ijk];
            }
    }
    else if (boundary_type == Neumann_type || boundary_type == Flux_type)
    {
        for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                ij  = i + j*jj;
                ijk = i + j*jj + (kend-1)*kk;
                a[ijk+kk] = agradtop[ij]*dzh[kend] + a[ijk];
            }
    }
}

void Boundary::calc_ghost_cells_bot_4th(double* restrict a, double* restrict z, Boundary_type boundary_type,
                                        double* restrict abot, double* restrict agradbot)
{
    int ij,ijk,jj,kk1,kk2,kstart;

    jj  = grid->icells;
    kk1 = 1*grid->ijcells;
    kk2 = 2*grid->ijcells;

    kstart = grid->kstart;

    if (boundary_type == Dirichlet_type)
    {
        for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                ij  = i + j*jj;
                ijk = i + j*jj + kstart*kk1;
                a[ijk-kk1] = (8./3.)*abot[ij] - 2.*a[ijk] + (1./3.)*a[ijk+kk1];
                a[ijk-kk2] = 8.*abot[ij] - 9.*a[ijk] + 2.*a[ijk+kk1];
            }
    }
    else if (boundary_type == Neumann_type || boundary_type == Flux_type)
    {
        using Finite_difference::O4::grad4x;

        for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                ij  = i + j*jj;
                ijk = i + j*jj + kstart*kk1;
                a[ijk-kk1] = -(1./24.)*grad4x(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk    ];
                a[ijk-kk2] = -(1./ 8.)*grad4x(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk+kk1];
            }
    }
}

void Boundary::calc_ghost_cells_top_4th(double* restrict a, double* restrict z, Boundary_type boundary_type,
                                        double* restrict atop, double* restrict agradtop)
{
    const int kend = grid->kend;

    const int jj  = grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    if (boundary_type == Dirichlet_type)
    {
        for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + (kend-1)*kk1;
                a[ijk+kk1] = (8./3.)*atop[ij] - 2.*a[ijk] + (1./3.)*a[ijk-kk1];
                a[ijk+kk2] = 8.*atop[ij] - 9.*a[ijk] + 2.*a[ijk-kk1];
            }
    }
    else if (boundary_type == Neumann_type || boundary_type == Flux_type)
    {
        using Finite_difference::O4::grad4x;

        for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + (kend-1)*kk1;
                a[ijk+kk1] = (1./24.)*grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk    ];
                a[ijk+kk2] = (1./ 8.)*grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk-kk1];
            }
    }
}

// BOUNDARY CONDITIONS FOR THE VERTICAL VELOCITY (NO PENETRATION)
void Boundary::calc_ghost_cells_botw_cons_4th(double* restrict w)
{
    const int jj  = grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    const int kstart = grid->kstart;

    for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
        for (int i=0; i<grid->icells; ++i)
        {
            const int ijk = i + j*jj + kstart*kk1;
            w[ijk-kk1] = -w[ijk+kk1];
            w[ijk-kk2] = -w[ijk+kk2];
        }
}

void Boundary::calc_ghost_cells_topw_cons_4th(double* restrict w)
{
    const int jj  = grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    const int kend = grid->kend;

    for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
        for (int i=0; i<grid->icells; ++i)
        {
            const int ijk = i + j*jj + kend*kk1;
            w[ijk+kk1] = -w[ijk-kk1];
            w[ijk+kk2] = -w[ijk-kk2];
        }
}

void Boundary::calc_ghost_cells_botw_4th(double* restrict w)
{
    const int jj  = grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;
    const int kk3 = 3*grid->ijcells;

    const int kstart = grid->kstart;

    for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
        for (int i=0; i<grid->icells; ++i)
        {
            const int ijk = i + j*jj + kstart*kk1;
            w[ijk-kk1] = -6.*w[ijk+kk1] + 4.*w[ijk+kk2] - w[ijk+kk3];
        }
}

void Boundary::calc_ghost_cells_topw_4th(double* restrict w)
{
    const int jj  = grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;
    const int kk3 = 3*grid->ijcells;

    const int kend = grid->kend;

    for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
        for (int i=0; i<grid->icells; ++i)
        {
            const int ijk = i + j*jj + kend*kk1;
            w[ijk+kk1] = -6.*w[ijk-kk1] + 4.*w[ijk-kk2] - w[ijk-kk3];
        }
}

void Boundary::prepare_device()
{
}

void Boundary::forward_device()
{
}

void Boundary::backward_device()
{
}
