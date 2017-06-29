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
#include <cmath>
#include "input.h"
#include "grid.h"
#include "master.h"
#include "fields.h"
#include "boundary_patch.h"
#include "boundary.h"
#include "defines.h"
#include "model.h"
#include "stats.h"

Boundary_patch::Boundary_patch(Model* modelin, Input* inputin) : Boundary(modelin, inputin)
{
}

void Boundary_patch::init(Input* inputin)
{
    int nerror = 0;

    // 1. Process the boundary conditions now all fields are registered
    process_bcs(inputin);

    // Patch type.
    nerror += inputin->get_item(&patch_dim,   "boundary", "patch_dim"  , "", 2 );
    nerror += inputin->get_item(&patch_xh,    "boundary", "patch_xh"   , "", 1.);
    nerror += inputin->get_item(&patch_xr,    "boundary", "patch_xr"   , "", 1.);
    nerror += inputin->get_item(&patch_xi,    "boundary", "patch_xi"   , "", 0.);
    nerror += inputin->get_item(&patch_xoffs, "boundary", "patch_xoffs", "", 0.);
    nerror += inputin->get_item(&patch_yh,    "boundary", "patch_yh"   , "", 1.);
    nerror += inputin->get_item(&patch_yr,    "boundary", "patch_yr"   , "", 1.);
    nerror += inputin->get_item(&patch_yi,    "boundary", "patch_yi"   , "", 0.);
    nerror += inputin->get_item(&patch_yoffs, "boundary", "patch_yoffs", "", 0.);

    // Process the patch factors for all scalars
    for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
        nerror += inputin->get_item(&patch_facr_map[it->first], "boundary", "patch_facr", it->first, 1.);
        nerror += inputin->get_item(&patch_facl_map[it->first], "boundary", "patch_facl", it->first, 0.);
    }

    if (nerror)
        throw 1;
}

void Boundary_patch::set_values()
{
    const double no_offset = 0.;

    set_bc(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot, mbcbot, ubot, fields->visc, grid->utrans);
    set_bc(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot, mbcbot, vbot, fields->visc, grid->vtrans);

    set_bc(fields->u->datatop, fields->u->datagradtop, fields->u->datafluxtop, mbctop, utop, fields->visc, grid->utrans);
    set_bc(fields->v->datatop, fields->v->datagradtop, fields->v->datafluxtop, mbctop, vtop, fields->visc, grid->vtrans);

    calc_patch(fields->atmp["tmp1"]->databot, grid->x, grid->y, patch_dim, 
               patch_xh, patch_xr, patch_xi, 
               patch_yh, patch_yr, patch_yi, 
               patch_xoffs, patch_yoffs);

    for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
        set_bc_patch(it->second->databot, it->second->datagradbot, it->second->datafluxbot, 
                     fields->atmp["tmp1"]->databot, patch_facl_map[it->first], patch_facr_map[it->first],
                     sbc[it->first]->bcbot, sbc[it->first]->bot, it->second->visc, no_offset);
        set_bc      (it->second->datatop, it->second->datagradtop, it->second->datafluxtop,
                     sbc[it->first]->bctop, sbc[it->first]->top, it->second->visc, no_offset);
    }
}

void Boundary_patch::get_mask(Field3d* field, Field3d* fieldh, Mask* m)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    // Switch between patch - no patch
    int sw;
    m->name == "patch_high" ? sw = 1 : sw = 0;

    // Calculate surface pattern
    calc_patch(fields->atmp["tmp1"]->databot, grid->x, grid->y, patch_dim, 
               patch_xh, patch_xr, patch_xi, 
               patch_yh, patch_yr, patch_yi, 
               patch_xoffs, patch_yoffs);

    // Set the values ranging between 0....1 to 0 or 1
    model->stats->nmaskh[grid->kstart] = 0;
    for (int j=grid->jstart; j<grid->jend; ++j)
        #pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij = i + j*jj;

            if (fields->atmp["tmp1"]->databot[ij] >= 0.5)
                fieldh->databot[ij] = sw;
            else
                fieldh->databot[ij] = 1-sw;

            if (fieldh->databot[ij] == 1)
                model->stats->nmaskh[grid->kstart] += 1;
        }

    // Sum the area, and propagate info
    master->sum(&model->stats->nmaskh[grid->kstart], 1);

    for (int k=grid->kstart+1; k<grid->kend; ++k)
        model->stats->nmaskh[k] = model->stats->nmaskh[grid->kstart];
    for (int k=grid->kstart; k<grid->kend; ++k)
        model->stats->nmask[k]  = model->stats->nmaskh[grid->kstart];
    model->stats->nmaskbot = model->stats->nmaskh[grid->kstart]; 

    // Set the atmospheric values
    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            #pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ij  = i + j*jj;
                const int ijk = i + j*jj + k*kk;

                field ->data[ijk] = fieldh->databot[ij];
                fieldh->data[ijk] = fieldh->databot[ij];
            }

    grid->boundary_cyclic   (field->data    );
    grid->boundary_cyclic   (fieldh->data   );
    grid->boundary_cyclic_2d(fieldh->databot);
}

void Boundary_patch::get_surface_mask(Field3d* field)
{
    const int jj = grid->icells;

    calc_patch(fields->atmp["tmp1"]->databot, grid->x, grid->y, patch_dim, 
               patch_xh, patch_xr, patch_xi, 
               patch_yh, patch_yr, patch_yi, 
               patch_xoffs, patch_yoffs);

    // Set the values ranging between 0....1 to 0 or 1
    for (int j=grid->jstart; j<grid->jend; ++j)
        #pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij = i + j*jj;

            if (fields->atmp["tmp1"]->databot[ij] >= 0.5)
                field->databot[ij] = 1;
            else
                field->databot[ij] = 0;
        }
}

void Boundary_patch::calc_patch(double* const restrict patch, const double* const restrict x, const double* const restrict y,
                                        const int patch_dim, 
                                        const double patch_xh, const double patch_xr, const double patch_xi,
                                        const double patch_yh, const double patch_yr, const double patch_yi,
                                        const double patch_xoffs, const double patch_yoffs) 
{
    const int jj = grid->icells;
    double errvalx, errvaly;

    for (int j=grid->jstart; j<grid->jend; ++j)
        #pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij = i + j*jj;
            const double xmod = fmod(x[i]-patch_xoffs, patch_xh);
            const double ymod = fmod(y[j]-patch_yoffs, patch_yh);

            errvalx = 0.5 - 0.5*erf(2.*(std::abs(2.*xmod - patch_xh) - patch_xr) / patch_xi);

            if (patch_dim == 2)
                errvaly = 0.5 - 0.5*erf(2.*(std::abs(2.*ymod - patch_yh) - patch_yr) / patch_yi);
            else
                errvaly = 1.;

            patch[ij] = errvalx*errvaly;
        }
}

void Boundary_patch::set_bc_patch(double* restrict a, double* restrict agrad, double* restrict aflux, 
                                          double* restrict patch, const double patch_facl, const double patch_facr,
                                          const int sw, const double aval, const double visc, const double offset)
{
    const int jj = grid->icells;

    const double avall = patch_facl*aval;
    const double avalr = patch_facr*aval;

    if (sw == Dirichlet_type)
    {
        for (int j=0; j<grid->jcells; ++j)
            #pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                const int ij = i + j*jj;
                a[ij] = avall + (avalr-avall)*patch[ij] - offset;
            }
    }
    else if (sw == Neumann_type)
    {
        for (int j=0; j<grid->jcells; ++j)
            #pragma ivdep
            for (int i=0; i<grid->icells; ++i)
            {
                const int ij = i + j*jj;
                agrad[ij] = avall + (avalr-avall)*patch[ij];
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
                aflux[ij] = avall + (avalr-avall)*patch[ij];
                agrad[ij] = -aflux[ij]/visc;
            }
    }
}
