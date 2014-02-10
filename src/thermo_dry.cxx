/*
 * MicroHH
 * Copyright (c) 2011-2013 Chiel van Heerwaarden
 * Copyright (c) 2011-2013 Thijs Heus
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
#include "thermo_dry.h"
#include "defines.h"

#define grav 9.81
#define Rd 287.04
#define cp 1005.

// make input options later
#define thvref 300.
#define psurf 100000.

cthermo_dry::cthermo_dry(cgrid *gridin, cfields *fieldsin, cmaster *masterin) : cthermo(gridin, fieldsin, masterin)
{
  swthermo = "dry";
}

cthermo_dry::~cthermo_dry()
{
}

int cthermo_dry::readinifile(cinput *inputin)
{
  int nerror = 0;
  nerror += inputin->getItem(&thref, "thermo", "thref", "");

  nerror += fields->initpfld("th");
  nerror += inputin->getItem(&fields->sp["th"]->visc, "fields", "svisc", "th");

  return nerror;
}

int cthermo_dry::init()
{
  // fields for anelastic solver
  pref    = new double[grid->kcells];
  exner   = new double[grid->kcells];
  rhoref  = new double[grid->kcells];

  prefh   = new double[grid->kcells];
  exnerh  = new double[grid->kcells];
  rhorefh = new double[grid->kcells];

  return 0;
}

int cthermo_dry::create()
{
  // ANELASTIC
  // calculate the base state pressure and density
  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    pref   [k] = psurf*std::exp(-grav/(Rd*thvref)*grid->z[k]);
    exner  [k] = std::pow(pref[k]/psurf, Rd/cp);
    rhoref [k] = pref[k] / (Rd*exner[k]*thvref);
  }

  for(int k=grid->kstart; k<grid->kend+1; ++k)
  {
    prefh  [k] = psurf*std::exp(-grav/(Rd*thvref)*grid->zh[k]);
    exnerh [k] = std::pow(prefh[k]/psurf, Rd/cp);
    rhorefh[k] = prefh[k] / (Rd*exnerh[k]*thvref);
  }

  return 0;
}

int cthermo_dry::exec()
{
  if(grid->swspatialorder== "2")
    calcbuoyancytend_2nd(fields->wt->data, fields->s["th"]->data);
  else if(grid->swspatialorder == "4")
    calcbuoyancytend_4th(fields->wt->data, fields->s["th"]->data);

  return 0;
}

int cthermo_dry::getbuoyancy(cfield3d *bfield, cfield3d *tmp)
{
  calcbuoyancy(bfield->data, fields->s["th"]->data);
  return 0;
}

int cthermo_dry::getbuoyancyfluxbot(cfield3d *bfield)
{
  calcbuoyancyfluxbot(bfield->datafluxbot, fields->s["th"]->datafluxbot);
  return 0;
}

int cthermo_dry::getbuoyancysurf(cfield3d *bfield)
{
  calcbuoyancybot(bfield->data, bfield->databot,
                  fields->s["th"]->data, fields->s["th"]->databot);
  calcbuoyancyfluxbot(bfield->datafluxbot, fields->s["th"]->datafluxbot);
  return 0;
}

int cthermo_dry::calcbuoyancy(double * restrict b, double * restrict th)
{
  int ijk,jj,kk;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  double thref  = this->thref;
  double gthref = gravity/this->thref;

  for(int k=0; k<grid->kcells; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        b[ijk] = gthref * (th[ijk] - thref);
      }

  return 0;
}

int cthermo_dry::calcbuoyancybot(double * restrict b , double * restrict bbot,
                                 double * restrict th, double * restrict thbot)
{
  int ij,ijk,jj,kk,kstart;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  kstart = grid->kstart;

  double thref  = this->thref;
  double gthref = gravity/this->thref;

  for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for(int i=0; i<grid->icells; ++i)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      bbot[ij] = gthref * (thbot[ij] - thref);
      b[ijk]   = gthref * (th[ijk] - thref);
    }

  return 0;
}

int cthermo_dry::calcbuoyancyfluxbot(double * restrict bfluxbot, double * restrict thfluxbot)
{
  int ij,jj;
  jj = grid->icells;

  double gthref = gravity/this->thref;

  for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for(int i=0; i<grid->icells; ++i)
    {
      ij  = i + j*jj;
      bfluxbot[ij] = gthref*thfluxbot[ij];
    }

  return 0;
}

int cthermo_dry::calcbuoyancytend_2nd(double * restrict wt, double * restrict th)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  double thref  = this->thref;
  double gthref = gravity/this->thref;

  for(int k=grid->kstart+1; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        wt[ijk] += gthref * (interp2(th[ijk-kk], th[ijk]) - thref);
      }

  return 0;
}

int cthermo_dry::calcbuoyancytend_4th(double * restrict wt, double * restrict th)
{
  int ijk,jj;
  int kk1,kk2;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  double thref  = this->thref;
  double gthref = gravity/this->thref;

  for(int k=grid->kstart+1; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk1;
        wt[ijk] += gthref * (interp4(th[ijk-kk2], th[ijk-kk1], th[ijk], th[ijk+kk1]) - thref);
      }

  return 0;
}

inline double cthermo_dry::interp2(const double a, const double b)
{
  return 0.5*(a + b);
}

inline double cthermo_dry::interp4(const double a, const double b, const double c, const double d)
{
  return (-a + 9.*b + 9.*c - d) / 16.;
}

