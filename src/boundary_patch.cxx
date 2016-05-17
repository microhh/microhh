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
#include "master.h"
#include "fields.h"
#include "boundary_patch.h"
#include "defines.h"
#include "model.h"

Boundary_patch::Boundary_patch(Model* modelin, Input* inputin) : Boundary_surface(modelin, inputin)
{
}

void Boundary_patch::init(Input* inputin)
{
    int nerror = 0;

    // 1. Process the boundary conditions now all fields are registered
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

    // 2. Read and check the boundary_surface specific settings
    process_input(inputin);

    // 3. Allocate and initialize the 2D surface fields
    init_surface();
}

void Boundary_patch::set_values()
{
    const double no_offset = 0.;

    set_bc(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot, mbcbot, ubot, fields->visc, grid->utrans);
    set_bc(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot, mbcbot, vbot, fields->visc, grid->vtrans);

    set_bc(fields->u->datatop, fields->u->datagradtop, fields->u->datafluxtop, mbctop, utop, fields->visc, grid->utrans);
    set_bc(fields->v->datatop, fields->v->datagradtop, fields->v->datafluxtop, mbctop, vtop, fields->visc, grid->vtrans);

    calc_patch(fields->atmp["tmp1"]->databot, grid->x, grid->y, patch_dim, patch_xh, patch_xr, patch_xi, patch_facr, patch_facl);

    for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
        set_bc_patch(it->second->databot, it->second->datagradbot, it->second->datafluxbot, fields->atmp["tmp1"]->databot,
                     sbc[it->first]->bcbot, sbc[it->first]->bot, it->second->visc, no_offset);
        set_bc      (it->second->datatop, it->second->datagradtop, it->second->datafluxtop,
                     sbc[it->first]->bctop, sbc[it->first]->top, it->second->visc, no_offset);
    }

    // in case the momentum has a fixed ustar, set the value to that of the input
    if (mbcbot == Ustar_type)
        set_ustar();

    // Prepare the lookup table for the surface solver
    init_solver();
}

void Boundary_patch::get_mask(Field3d* field, Field3d* fieldh, Mask* m)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    // Switch between patch - no patch
    int sw;
    m->name == "patch_high" ? sw = 1 : sw = 0;

    // Calculate surface pattern
    calc_patch(fieldh->databot, grid->x, grid->y, patch_dim, patch_xh, patch_xr, patch_xi, patch_facr, patch_facl);

    // Set the values ranging between 0....1 to 0 or 1
    const double threshold = 0.5 * (patch_facr + patch_facl);
    model->stats->nmaskh[grid->kstart] = 0;
    for (int j=grid->jstart; j<grid->jend; ++j)
        #pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij = i + j*jj;

            if (fieldh->databot[ij] >= threshold)
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

    calc_patch(field->databot, grid->x, grid->y, patch_dim, patch_xh, patch_xr, patch_xi, patch_facr, patch_facl);

    // Set the values ranging between 0....1 to 0 or 1
    const double threshold = 0.5 * (patch_facr + patch_facl);

    for (int j=grid->jstart; j<grid->jend; ++j)
        #pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij = i + j*jj;

            if (field->databot[ij] >= threshold)
                field->databot[ij] = 1;
            else
                field->databot[ij] = 0;
        }
}

void Boundary_patch::calc_patch(double* const restrict patch, const double* const restrict x, const double* const restrict y,
                                const int patch_dim, const double patch_xh, const double patch_xr, const double patch_xi,
                                const double patch_facr, const double patch_facl) 
{
    const int jj = grid->icells;
    double errvalx, errvaly;

    for (int j=grid->jstart; j<grid->jend; ++j)
        #pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij = i + j*jj;
            const double xmod = fmod(x[i], patch_xh);
            const double ymod = fmod(y[j], patch_xh);

            errvalx = 0.5 - 0.5*erf(2.*(std::abs(2.*xmod - patch_xh) - patch_xr) / patch_xi);

            if (patch_dim == 2)
                errvaly = 0.5 - 0.5*erf(2.*(std::abs(2.*ymod - patch_xh) - patch_xr) / patch_xi);
            else
                errvaly = 1.;

            patch[ij] = errvalx*errvaly;
        }
}

void Boundary_patch::set_bc_patch(double* restrict a, double* restrict agrad, double* restrict aflux, double* restrict patch, 
                                  int sw, double aval, double visc, double offset)
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
