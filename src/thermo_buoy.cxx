/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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

cthermo_buoy::cthermo_buoy(cmodel *modelin, cinput *inputin) : cthermo(modelin, inputin)
{
  swthermo = "buoy";

  int nerror = 0;
  nerror += fields->initpfld("b", "Buoyancy", "m s-2");
  nerror += inputin->getItem(&fields->sp["b"]->visc, "fields", "svisc", "b");

  if(nerror)
    throw 1;
}

cthermo_buoy::~cthermo_buoy()
{
}

#ifndef USECUDA
int cthermo_buoy::exec()
{
  // extend later for gravity vector not normal to surface
  if(grid->swspatialorder== "2")
    calcbuoyancytend_2nd(fields->wt->data, fields->s["b"]->data);
  else if(grid->swspatialorder == "4")
    calcbuoyancytend_4th(fields->wt->data, fields->s["b"]->data);

  return 0;
}
#endif

int cthermo_buoy::getbuoyancy(cfield3d *bfield, cfield3d *tmp)
{
  calcbuoyancy(bfield->data, fields->s["b"]->data);
  return 0;
}

int cthermo_buoy::getbuoyancyfluxbot(cfield3d *bfield)
{
  calcbuoyancyfluxbot(bfield->datafluxbot, fields->s["b"]->datafluxbot);
  return 0;
}

int cthermo_buoy::getbuoyancysurf(cfield3d *bfield)
{
  calcbuoyancybot(bfield->data, bfield->databot,
                  fields->s["b"]->data, fields->s["b"]->databot);
  calcbuoyancyfluxbot(bfield->datafluxbot, fields->s["b"]->datafluxbot);
  return 0;
}

int cthermo_buoy::checkthermofield(std::string name)
{
  if(name == "b")
    return 0;
  else
    return 1;
}

int cthermo_buoy::getthermofield(cfield3d *field, cfield3d *tmp, std::string name)
{
  calcbuoyancy(field->data, fields->s["b"]->data);
  return 0;
}

int cthermo_buoy::getprogvars(std::vector<std::string> *list)
{
  list->push_back("b");
  return 0;
}

int cthermo_buoy::calcbuoyancy(double * restrict b, double * restrict bin)
{
  int ijk,jj,kk;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(int n=0; n<grid->ncells; ++n)
    b[n] = bin[n];

  return 0;
}

int cthermo_buoy::calcbuoyancybot(double * restrict b  , double * restrict bbot,
                                  double * restrict bin, double * restrict binbot)
{
  int ij,ijk,jj,kk,kstart;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  kstart = grid->kstart;

  for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for(int i=0; i<grid->icells; ++i)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      bbot[ij] = binbot[ij];
      b[ijk]   = bin[ijk];
    }

  return 0;
}

int cthermo_buoy::calcbuoyancyfluxbot(double * restrict bfluxbot, double * restrict binfluxbot)
{
  int ij,jj;
  jj = grid->icells;

  for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for(int i=0; i<grid->icells; ++i)
    {
      ij  = i + j*jj;
      bfluxbot[ij] = binfluxbot[ij];
    }

  return 0;
}

int cthermo_buoy::calcbuoyancytend_2nd(double * restrict wt, double * restrict b)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(int k=grid->kstart+1; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        wt[ijk] += interp2(b[ijk-kk], b[ijk]);
      }

  return 0;
}

int cthermo_buoy::calcbuoyancytend_4th(double * restrict wt, double * restrict b)
{
  int ijk,jj;
  int kk1,kk2;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  for(int k=grid->kstart+1; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk1;
        wt[ijk] += interp4(b[ijk-kk2], b[ijk-kk1], b[ijk], b[ijk+kk1]);
      }

  return 0;
}
