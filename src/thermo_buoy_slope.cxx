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
#include "grid.h"
#include "fields.h"
#include "thermo_buoy_slope.h"
#include "defines.h"
#include "fd.h"
#include "model.h"
#include "master.h"

using fd::o4::interp4;

Thermo_buoy_slope::Thermo_buoy_slope(Model *modelin, Input *inputin) : Thermo(modelin, inputin)
{
  swthermo = "slope";
  master = modelin->master;

  fields->init_prognostic_field("b", "Buoyancy", "m s-2");

  int nerror = 0;
  nerror += inputin->getItem(&alpha, "thermo", "alpha", "");
  nerror += inputin->getItem(&n2   , "thermo", "n2"   , "");
  nerror += inputin->getItem(&fields->sp["b"]->visc, "fields", "svisc", "b");

  if (grid->swspatialorder == "2")
  {
    master->print_error("swthermo = buoy_slope is incompatible with swspatialorder = 2\n");
    ++nerror;
  }

  if (nerror)
    throw 1;
}

Thermo_buoy_slope::~Thermo_buoy_slope()
{
}

#ifndef USECUDA
void Thermo_buoy_slope::exec()
{
  if (grid->swspatialorder == "2")
  {
    master->print_error("Second order not implemented for slope flow thermodynamics\n");
    throw 1;
  }
  else if (grid->swspatialorder == "4")
  {
    calcBuoyancyTend_u_4th(fields->ut->data, fields->sp["b"]->data);
    calcBuoyancyTend_w_4th(fields->wt->data, fields->sp["b"]->data);
    calcBuoyancyTend_b_4th(fields->st["b"]->data, fields->u->data, fields->w->data);
  }
}
#endif

bool Thermo_buoy_slope::checkThermoField(std::string name)
{
  if (name == "b")
    return false;
  else
    return true;
}

void Thermo_buoy_slope::getThermoField(Field3d *field, Field3d *tmp, std::string name)
{
  calcBuoyancy(field->data, fields->sp["b"]->data);
}

void Thermo_buoy_slope::getBuoyancyFluxbot(Field3d *bfield)
{
  calcBuoyancyFluxbot(bfield->datafluxbot, fields->sp["b"]->datafluxbot);
}

void Thermo_buoy_slope::getBuoyancySurf(Field3d *bfield)
{
  calcBuoyancyBot(bfield->data        , bfield->databot,
                  fields->sp["b"]->data, fields->sp["b"]->databot);
  calcBuoyancyFluxbot(bfield->datafluxbot, fields->sp["b"]->datafluxbot);
}

void Thermo_buoy_slope::getProgVars(std::vector<std::string> *list)
{
  list->push_back("b");
}

void Thermo_buoy_slope::calcBuoyancy(double * restrict b, double * restrict s)
{
  const int jj = grid->icells;
  const int kk = grid->ijcells;

  for (int k=0; k<grid->kcells; ++k)
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for (int i=grid->istart; i<grid->iend; ++i)
      {
        const int ijk = i + j*jj + k*kk;
        b[ijk] = s[ijk];
      }
}

void Thermo_buoy_slope::calcBuoyancyBot(double * restrict b, double * restrict bbot,
                                      double * restrict s, double * restrict sbot)
{
  int ij,ijk,jj,kk,kstart;
  jj = grid->icells;
  kk = grid->ijcells;
  kstart = grid->kstart;

  for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for (int i=0; i<grid->icells; ++i)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      bbot[ij] = sbot[ij];
      b[ijk]   = s[ijk];
    }
}

void Thermo_buoy_slope::calcBuoyancyFluxbot(double * restrict bfluxbot, double * restrict sfluxbot)
{
  const int jj = grid->icells;

  for (int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for (int i=0; i<grid->icells; ++i)
    {
      const int ij  = i + j*jj;
      bfluxbot[ij] = sfluxbot[ij];
    }
}

void Thermo_buoy_slope::calcBuoyancyTend_u_4th(double * restrict ut, double * restrict b)
{
  int ijk,ii1,ii2,jj,kk;

  ii1 = 1;
  ii2 = 2;
  jj  = grid->icells;
  kk  = grid->ijcells;

  const double sinalpha = std::sin(this->alpha);

  for (int k=grid->kstart; k<grid->kend; ++k)
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for (int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        ut[ijk] += sinalpha * interp4(b[ijk-ii2], b[ijk-ii1], b[ijk], b[ijk+ii1]);
      }
}

void Thermo_buoy_slope::calcBuoyancyTend_w_4th(double * restrict wt, double * restrict b)
{
  int ijk,jj;
  int kk1,kk2;

  jj  = grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  const double cosalpha = std::cos(this->alpha);

  for (int k=grid->kstart+1; k<grid->kend; ++k)
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for (int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk1;
        wt[ijk] += cosalpha * interp4(b[ijk-kk2], b[ijk-kk1], b[ijk], b[ijk+kk1]);
      }
}

void Thermo_buoy_slope::calcBuoyancyTend_b_4th(double * restrict bt, double * restrict u, double * restrict w)
{
  int ijk,ii1,ii2,jj,kk1,kk2;

  ii1 = 1;
  ii2 = 2;
  jj  = 1*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  const double sinalpha = std::sin(this->alpha);
  const double cosalpha = std::cos(this->alpha);
  const double n2 = this->n2;
  const double utrans = grid->utrans;

  for (int k=grid->kstart; k<grid->kend; ++k)
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for (int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk1;
        bt[ijk] -= n2 * ( sinalpha * (interp4(u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2]) + utrans)
                        + cosalpha *  interp4(w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2]) );
      }
}
