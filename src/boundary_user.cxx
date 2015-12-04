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
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "boundary_user.h"
#include "defines.h"
#include "model.h"

Boundary_user::Boundary_user(Model* modelin, Input* inputin) : Boundary(modelin, inputin)
{
}

void Boundary_user::init(Input* inputin)
{
    int nerror = 0;

    process_bcs(inputin);

    // Patch type.
    nerror += inputin->get_item(&patch_dim,  "boundary", "patch_dim" , "", 2 );
    nerror += inputin->get_item(&patch_xh,   "boundary", "patch_xh"  , "", 1.);
    nerror += inputin->get_item(&patch_xr,   "boundary", "patch_xr"  , "", 1.);
    nerror += inputin->get_item(&patch_xi,   "boundary", "patch_xi"  , "", 0.);
    nerror += inputin->get_item(&patch_facr, "boundary", "patch_facr", "", 1.);
    nerror += inputin->get_item(&patch_facl, "boundary", "patch_facl", "", 0.);

    if (nerror)
        throw 1;
}

void Boundary_user::set_values()
{
    const double no_offset = 0.;

    set_bc(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot, mbcbot, ubot, fields->visc, grid->utrans);
    set_bc(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot, mbcbot, vbot, fields->visc, grid->vtrans);

    set_bc(fields->u->datatop, fields->u->datagradtop, fields->u->datafluxtop, mbctop, utop, fields->visc, grid->utrans);
    set_bc(fields->v->datatop, fields->v->datagradtop, fields->v->datafluxtop, mbctop, vtop, fields->visc, grid->vtrans);

    for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
        set_bc_patch(it->second->databot, it->second->datagradbot, it->second->datafluxbot,
                     sbc[it->first]->bcbot, sbc[it->first]->bot, it->second->visc, no_offset, fields->atmp["tmp1"]->data, patch_facl, patch_facr);
        set_bc      (it->second->datatop, it->second->datagradtop, it->second->datafluxtop,
                     sbc[it->first]->bctop, sbc[it->first]->top, it->second->visc, no_offset);
    }
}

void Boundary_user::set_bc_patch(double* restrict a, double* restrict agrad, double* restrict aflux, int sw, double aval, double visc, double offset,
                                 double* restrict tmp, double facl, double facr)
{
    const int jj = grid->icells;

    const double avall = facl*aval;
    const double avalr = facr*aval;

    double errvalx, errvaly;

    // save the pattern
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij = i + j*jj;
            const double xmod = fmod(grid->x[i], patch_xh);
            const double ymod = fmod(grid->y[j], patch_xh);

            errvalx = 0.5 - 0.5*erf(2.*(std::abs(2.*xmod - patch_xh) - patch_xr) / patch_xi);

            if (patch_dim == 2)
                errvaly = 0.5 - 0.5*erf(2.*(std::abs(2.*ymod - patch_xh) - patch_xr) / patch_xi);
            else
                errvaly = 1.;

            tmp[ij] = errvalx*errvaly;
        }

    if (sw == Dirichlet_type)
    {
        for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                const int ij = i + j*jj;
                a[ij] = avall + (avalr-avall)*tmp[ij] - offset;
            }
    }
    else if (sw == Neumann_type)
    {
        for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                const int ij = i + j*jj;
                agrad[ij] = avall + (avalr-avall)*tmp[ij];
                aflux[ij] = -agrad[ij]*visc;
            }
    }
    else if (sw == Flux_type)
    {
        for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                const int ij = i + j*jj;
                aflux[ij] = avall + (avalr-avall)*tmp[ij];
                agrad[ij] = -aflux[ij]/visc;
            }
    }
}
