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
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "boundary_surface_bulk.h"
#include "defines.h"
#include "constants.h"
#include "thermo.h"
#include "model.h"
#include "master.h"
#include "cross.h"

Boundary_surface_bulk::Boundary_surface_bulk(Model* modelin, Input* inputin) : Boundary_surface(modelin, inputin)
{
    #ifdef USECUDA
    master->print_error("swboundary=\"surface_bulk\" not (yet) implemented in CUDA\n");
    throw 1;
    #endif
}

Boundary_surface_bulk::~Boundary_surface_bulk()
{
}

void Boundary_surface_bulk::init(Input *inputin)
{    
    stats = model->stats;

    // 1. Process the boundary conditions now all fields are registered
    process_bcs(inputin);

    int nerror = 0;
    nerror += inputin->get_item(&z0m, "boundary", "z0m", "");
    nerror += inputin->get_item(&z0h, "boundary", "z0h", "");

    // Read list of cross sections
    nerror += inputin->get_list(&crosslist , "boundary", "crosslist" , "");

    // Allow only mbcbot = noslip with sbcbot = dirichlet
    if(mbcbot != Dirichlet_type)
    {
        master->print_error("Only \"noslip\" is allowed as mbcbot with swboundary=\"bulk\"\n");
        ++nerror;
    }

    // Read bulk transfer coefficient momentum
    nerror += inputin->get_item(&bulk_cm, "boundary", "bulk_cm", "");

    // process the scalars
    for (BcMap::const_iterator it=sbc.begin(); it!=sbc.end(); ++it)
    {
        // crash in case fixed gradient is prescribed
        if (it->second->bcbot != Dirichlet_type)
        {
            master->print_error("Only \"dirichlet\" is allowed as sbcbot with swboundary=\"bulk\"\n");
            ++nerror;
        }

        // Read bulk transfer coefficient
        nerror += inputin->get_item(&bulk_cs[it->first], "boundary", "bulk_cs", it->first);
    }

    if (nerror)
        throw 1;

    // 2. Allocate the fields
    obuk  = new double[grid->ijcells];
    ustar = new double[grid->ijcells];

    // Cross sections
    allowedcrossvars.push_back("ustar");
    allowedcrossvars.push_back("obuk");

    // Get global cross-list from cross.cxx
    std::vector<std::string> *crosslist_global = model->cross->get_crosslist(); 

    // Check input list of cross variables (crosslist)
    std::vector<std::string>::iterator it2=crosslist_global->begin();
    while (it2 != crosslist_global->end())
    {
        if (std::count(allowedcrossvars.begin(),allowedcrossvars.end(),*it2))
        {
            // Remove variable from global list, put in local list
            crosslist.push_back(*it2);
            crosslist_global->erase(it2); // erase() returns iterator of next element..
        }
        else
            ++it2;
    }
}

void Boundary_surface_bulk::set_values()
{
    const double no_velocity = 0.;
    const double no_offset = 0.;

    // grid transformation is properly taken into account by setting the databot and top values
    set_bc(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot, mbcbot, no_velocity, fields->visc, grid->utrans);
    set_bc(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot, mbcbot, no_velocity, fields->visc, grid->vtrans);

    set_bc(fields->u->datatop, fields->u->datagradtop, fields->u->datafluxtop, mbctop, no_velocity, fields->visc, grid->utrans);
    set_bc(fields->v->datatop, fields->v->datagradtop, fields->v->datafluxtop, mbctop, no_velocity, fields->visc, grid->vtrans);

    for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
        set_bc(it->second->databot, it->second->datagradbot, it->second->datafluxbot, sbc[it->first]->bcbot, sbc[it->first]->bot, it->second->visc, no_offset);
        set_bc(it->second->datatop, it->second->datagradtop, it->second->datafluxtop, sbc[it->first]->bctop, sbc[it->first]->top, it->second->visc, no_offset);
    }
}

//#ifndef USECUDA
void Boundary_surface_bulk::calculate_du(double* restrict dutot, double* restrict u, double* restrict v, double* restrict ubot, double* restrict vbot)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    // calculate total wind
    double du2;
    const double minval = 1.e-1;
    for (int j=grid->jstart; j<grid->jend; ++j)
        #pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            du2 = std::pow(0.5*(u[ijk] + u[ijk+ii]) - 0.5*(ubot[ij] + ubot[ij+ii]), 2)
                + std::pow(0.5*(v[ijk] + v[ijk+jj]) - 0.5*(vbot[ij] + vbot[ij+jj]), 2);
            dutot[ij] = std::max(std::pow(du2, 0.5), minval);
        }

    grid->boundary_cyclic_2d(dutot);
}

void Boundary_surface_bulk::momentum_fluxgrad(double* restrict ufluxbot, double* restrict vfluxbot, 
                                      double* restrict ugradbot, double* restrict vgradbot,
                                      double* restrict u, double* restrict v, 
                                      double* restrict ubot, double* restrict vbot, 
                                      double* restrict dutot, const double Cm, const double zsl)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    for (int j=grid->jstart; j<grid->jend; ++j)
        #pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;

            ufluxbot[ij] = -Cm * dutot[ij] * (u[ijk]-ubot[ij]);
            vfluxbot[ij] = -Cm * dutot[ij] * (v[ijk]-vbot[ij]);
        }

    grid->boundary_cyclic_2d(ufluxbot);
    grid->boundary_cyclic_2d(vfluxbot);

    for (int j=0; j<grid->jcells; ++j)
        #pragma ivdep
        for (int i=0; i<grid->icells; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            // use the linearly interpolated grad, rather than the MO grad,
            // to prevent giving unresolvable gradients to advection schemes
            ugradbot[ij] = (u[ijk]-ubot[ij])/zsl;
            vgradbot[ij] = (v[ijk]-vbot[ij])/zsl;
        }
}

void Boundary_surface_bulk::scalar_fluxgrad(double* restrict sfluxbot, double* restrict sgradbot, double* restrict s, double* restrict sbot,
                                    double* restrict dutot, const double Cs, const double zsl) 
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    for (int j=grid->jstart; j<grid->jend; ++j)
        #pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;

            sfluxbot[ij] = -Cs * dutot[ij] * (s[ijk]-sbot[ij]);
        }

    grid->boundary_cyclic_2d(sfluxbot);

    for (int j=0; j<grid->jcells; ++j)
        #pragma ivdep
        for (int i=0; i<grid->icells; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            sgradbot[ij] = (s[ijk]-sbot[ij])/zsl;
        }
}

void Boundary_surface_bulk::surface_scaling(double* restrict ustar, double* restrict obuk, double* restrict dutot, double* restrict bfluxbot, const double Cm)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const double sqrt_Cm = pow(Cm, 0.5);

    for (int j=grid->jstart; j<grid->jend; ++j)
        #pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij  = i + j*jj;

            ustar[ij] = sqrt_Cm * dutot[ij];
            obuk[ij] = -pow(ustar[ij], 3) / (Constants::kappa * bfluxbot[ij]);
        }
}

void Boundary_surface_bulk::update_bcs()
{
    const double zsl = grid->z[grid->kstart];

    // Calculate total wind speed difference with surface
    calculate_du(fields->atmp["tmp1"]->data, fields->u->data, fields->v->data, fields->u->databot, fields->v->databot);

    // Calculate surface momentum fluxes and gradients
    momentum_fluxgrad(fields->u->datafluxbot, fields->v->datafluxbot, fields->u->datagradbot, fields->v->datagradbot,
                      fields->u->data, fields->v->data, fields->u->databot, fields->v->databot, fields->atmp["tmp1"]->data, bulk_cm, zsl);

    // Calculate surface scalar fluxes and gradients
    for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
        scalar_fluxgrad(it->second->datafluxbot, it->second->datagradbot, it->second->data, it->second->databot, fields->atmp["tmp1"]->data, bulk_cs[it->first], zsl);
    
    // Calculate Obukhov length and ustar
    model->thermo->get_buoyancy_fluxbot(fields->atmp["tmp2"]);
    surface_scaling(ustar, obuk, fields->atmp["tmp1"]->data, fields->atmp["tmp2"]->datafluxbot, bulk_cm); 
}
//#endif
