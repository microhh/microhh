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
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "boundary_bulk.h"
#include "defines.h"
#include "constants.h"
#include "thermo.h"
#include "model.h"
#include "master.h"
#include "cross.h"

Boundary_bulk::Boundary_bulk(Model* modelin, Input* inputin) : Boundary_surface(modelin, inputin)
{
}

Boundary_bulk::~Boundary_bulk()
{
}

void Boundary_bulk::init(Input *inputin)
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
        master->print_error("Only \"noslip\" is allowed as mbcbot with swboundary=bulk\n");
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
            master->print_error("Only \"dirichlet\" is allowed as sbcbot with swboundary=bulk\n");
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

void Boundary_bulk::set_values()
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

#ifndef USECUDA
void Boundary_bulk::update_bcs()
{
}
#endif
