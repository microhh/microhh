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
#include "grid.h"
#include "fields.h"
#include "thermo_buoy.h"
#include "defines.h"
#include "fd.h"

using fd::o2::interp2;
using fd::o4::interp4;

Thermo_buoy::Thermo_buoy(Model* modelin, Input* inputin) : Thermo(modelin, inputin)
{
    swthermo = "buoy";

    fields->init_prognostic_field("b", "Buoyancy", "m s-2");

    int nerror = 0;
    nerror += inputin->get_item(&fields->sp["b"]->visc, "fields", "svisc", "b");

    if (nerror)
        throw 1;
}

Thermo_buoy::~Thermo_buoy()
{
}

#ifndef USECUDA
void Thermo_buoy::exec()
{
    // extend later for gravity vector not normal to surface
    if (grid->swspatialorder== "2")
        calc_buoyancy_tend_2nd(fields->wt->data, fields->sp["b"]->data);
    else if (grid->swspatialorder == "4")
        calc_buoyancy_tend_4th(fields->wt->data, fields->sp["b"]->data);
}
#endif

void Thermo_buoy::get_buoyancy(Field3d* bfield, Field3d* tmp)
{
    calc_buoyancy(bfield->data, fields->sp["b"]->data);
}

void Thermo_buoy::get_buoyancy_fluxbot(Field3d* bfield)
{
    calc_buoyancy_fluxbot(bfield->datafluxbot, fields->sp["b"]->datafluxbot);
}

void Thermo_buoy::get_buoyancy_surf(Field3d* bfield)
{
    calc_buoyancy_bot(bfield->data, bfield->databot,
                      fields->sp["b"]->data, fields->sp["b"]->databot);
    calc_buoyancy_fluxbot(bfield->datafluxbot, fields->sp["b"]->datafluxbot);
}

bool Thermo_buoy::check_thermo_field(std::string name)
{
    if (name == "b")
        return false;
    else
        return true;
}

void Thermo_buoy::get_thermo_field(Field3d *field, Field3d *tmp, std::string name)
{
    calc_buoyancy(field->data, fields->sp["b"]->data);
}

void Thermo_buoy::get_prog_vars(std::vector<std::string> *list)
{
    list->push_back("b");
}

void Thermo_buoy::calc_buoyancy(double* restrict b, double* restrict bin)
{
    for (int n=0; n<grid->ncells; ++n)
        b[n] = bin[n];
}

void Thermo_buoy::calc_buoyancy_bot(double* restrict b  , double* restrict bbot,
                                    double* restrict bin, double* restrict binbot)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
        for (int i=0; i<grid->icells; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            bbot[ij] = binbot[ij];
            b[ijk]   = bin[ijk];
        }
}

void Thermo_buoy::calc_buoyancy_fluxbot(double* restrict bfluxbot, double* restrict binfluxbot)
{
    const int jj = grid->icells;

    for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
        for (int i=0; i<grid->icells; ++i)
        {
            const int ij  = i + j*jj;
            bfluxbot[ij] = binfluxbot[ij];
        }
}

void Thermo_buoy::calc_buoyancy_tend_2nd(double* restrict wt, double* restrict b)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    for (int k=grid->kstart+1; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                wt[ijk] += interp2(b[ijk-kk], b[ijk]);
            }
}

void Thermo_buoy::calc_buoyancy_tend_4th(double* restrict wt, double* restrict b)
{
    const int jj  = grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    for (int k=grid->kstart+1; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk1;
                wt[ijk] += interp4(b[ijk-kk2], b[ijk-kk1], b[ijk], b[ijk+kk1]);
            }
}
