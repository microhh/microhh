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
#include <iostream>
#include "grid.h"
#include "fields.h"
#include "thermo_buoy.h"
#include "defines.h"
#include "finite_difference.h"
#include "model.h"
#include "master.h"

using Finite_difference::O2::interp2;
using Finite_difference::O4::interp4;

Thermo_buoy::Thermo_buoy(Model *modelin, Input *inputin) : Thermo(modelin, inputin)
{
    swthermo = "buoy";
    master = modelin->master;

    fields->init_prognostic_field("b", "Buoyancy", "m s-2");

    int nerror = 0;
    nerror += inputin->get_item(&alpha, "thermo", "alpha", "");
    nerror += inputin->get_item(&n2   , "thermo", "N2"   , "");
    nerror += inputin->get_item(&fields->sp["b"]->visc, "fields", "svisc", "b");
	
	slope = alpha > 0;

    if (nerror)
        throw 1;
}

Thermo_buoy::~Thermo_buoy()
{
}

#ifndef USECUDA
void Thermo_buoy::exec()
{
    if (grid->swspatialorder == "2") 
    {
	    if (slope) 
	    {
            calc_buoyancy_tend_u_2nd(fields->ut->data, fields->sp["b"]->data);
            calc_buoyancy_tend_w_2nd(fields->wt->data, fields->sp["b"]->data);
            calc_buoyancy_tend_b_2nd(fields->st["b"]->data, fields->u->data, fields->w->data);
        } 
        else 
        {
	        calc_buoyancy_tend_2nd(fields->wt->data, fields->sp["b"]->data);
        }
    } 
    else if (grid->swspatialorder == "4") 
    {    
	    if (slope) 
	    {
		    calc_buoyancy_tend_u_4th(fields->ut->data, fields->sp["b"]->data);
            calc_buoyancy_tend_w_4th(fields->wt->data, fields->sp["b"]->data);
            calc_buoyancy_tend_b_4th(fields->st["b"]->data, fields->u->data, fields->w->data);
	    }
	    else 
	    {
		    calc_buoyancy_tend_4th(fields->wt->data, fields->sp["b"]->data);
	    }
    }
}
#endif

void Thermo_buoy::get_thermo_field(Field3d* field, Field3d* tmp, const std::string name)
{
    calc_buoyancy(field->data, fields->sp["b"]->data);
}

void Thermo_buoy::get_prog_vars(std::vector<std::string>* list)
{
    list->push_back("b");
}

void Thermo_buoy::get_buoyancy_fluxbot(Field3d* bfield)
{
    calc_buoyancy_fluxbot(bfield->datafluxbot, fields->sp["b"]->datafluxbot);
}

void Thermo_buoy::get_buoyancy_surf(Field3d *bfield)
{
    calc_buoyancy_bot(bfield->data         , bfield->databot,
                      fields->sp["b"]->data, fields->sp["b"]->databot);
    calc_buoyancy_fluxbot(bfield->datafluxbot, fields->sp["b"]->datafluxbot);
}

bool Thermo_buoy::check_field_exists(std::string name)
{
    if (name == "b")
        return true;
    else
        return false;
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

void Thermo_buoy::calc_buoyancy_tend_u_2nd(double* restrict ut, double* restrict b)
{
    const int ii1 = 1;
    const int jj  = grid->icells;
    const int kk  = grid->ijcells;

    const double sinalpha = std::sin(this->alpha);
    
    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                ut[ijk] += sinalpha * interp2(b[ijk-ii1], b[ijk]);
            }
}

void Thermo_buoy::calc_buoyancy_tend_w_2nd(double* restrict wt, double* restrict b)
{
    const int jj  = grid->icells;
    const int kk1 = 1*grid->ijcells;

    const double cosalpha = std::cos(this->alpha);
    
    for (int k=grid->kstart+1; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk1;
                wt[ijk] += cosalpha * interp2(b[ijk-kk1], b[ijk]);
            }
}

void Thermo_buoy::calc_buoyancy_tend_b_2nd(double* restrict bt, double* restrict u, double* restrict w)
{
    const int ii1 = 1;
    const int jj  = 1*grid->icells;
    const int kk1 = 1*grid->ijcells;

    const double sinalpha = std::sin(this->alpha);
    const double cosalpha = std::cos(this->alpha);
    const double n2 = this->n2;
    const double utrans = grid->utrans;
    
    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk1;
                bt[ijk] -= n2 * ( sinalpha * (interp2(u[ijk], u[ijk+ii1]) + utrans)
                                + cosalpha *  interp2(w[ijk], w[ijk+kk1]) );
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

void Thermo_buoy::calc_buoyancy_tend_u_4th(double* restrict ut, double* restrict b)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int jj  = grid->icells;
    const int kk  = grid->ijcells;

    const double sinalpha = std::sin(this->alpha);
    
    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                ut[ijk] += sinalpha * interp4(b[ijk-ii2], b[ijk-ii1], b[ijk], b[ijk+ii1]);
            }
}

void Thermo_buoy::calc_buoyancy_tend_w_4th(double* restrict wt, double* restrict b)
{
    const int jj  = grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    const double cosalpha = std::cos(this->alpha);
    
    for (int k=grid->kstart+1; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk1;
                wt[ijk] += cosalpha * interp4(b[ijk-kk2], b[ijk-kk1], b[ijk], b[ijk+kk1]);
            }
}

void Thermo_buoy::calc_buoyancy_tend_b_4th(double* restrict bt, double* restrict u, double* restrict w)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int jj  = 1*grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    const double sinalpha = std::sin(this->alpha);
    const double cosalpha = std::cos(this->alpha);
    const double n2 = this->n2;
    const double utrans = grid->utrans;
    
    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk1;
                bt[ijk] -= n2 * ( sinalpha * (interp4(u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2]) + utrans)
                                + cosalpha *  interp4(w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2]) );
            }
}