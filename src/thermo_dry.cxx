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

#define gravity 9.81
#define Rd 287.04
#define cp 1005.

// make input options later
#define psurf 100000.
#define gammath 0.003
#define thstart 300.

cthermo_dry::cthermo_dry(cmodel *modelin) : cthermo(modelin)
{
  swthermo = "dry";
}

cthermo_dry::~cthermo_dry()
{
}

int cthermo_dry::readinifile(cinput *inputin)
{
  int nerror = 0;
  // nerror += inputin->getItem(&thref0, "thermo", "thref0", "");

  nerror += fields->initpfld("th");
  nerror += inputin->getItem(&fields->sp["th"]->visc, "fields", "svisc", "th");

  return nerror;
}

int cthermo_dry::init()
{
  // fields for anelastic solver
  thref   = new double[grid->kcells];
  pref    = new double[grid->kcells];
  exner   = new double[grid->kcells];
  // rhoref  = new double[grid->kcells];

  threfh  = new double[grid->kcells];
  prefh   = new double[grid->kcells];
  exnerh  = new double[grid->kcells];
  // rhorefh = new double[grid->kcells];

  return 0;
}

int cthermo_dry::create(cinput *inputin)
{
  // take the initial profile as the reference
  if(inputin->getProf(&thref[grid->kstart], "th", grid->kmax))
    return 1;

  int kstart = grid->kstart;
  int kend   = grid->kend;

  // extrapolate the profile to get the bottom value
  threfh[kstart] = thref[kstart] - grid->z[kstart]*(thref[kstart+1]-thref[kstart])*grid->dzhi[kstart+1];

  // extrapolate the profile to get the top value
  threfh[kend] = thref[kend-1] + (grid->zh[kend]-grid->z[kend-1])*(thref[kend-1]-thref[kend-2])*grid->dzhi[kend-1];

  // set the ghost cells for the reference temperature
  thref[kstart-1] = 2.*threfh[kstart] - thref[kstart];
  thref[kend]     = 2.*threfh[kend]   - thref[kend-1];

  // interpolate the reference temperature profile
  for(int k=grid->kstart+1; k<grid->kend; ++k)
    threfh[k] = 0.5*(thref[k-1] + thref[k]);

  // ANELASTIC
  // calculate the base state pressure and density
  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    pref [k] = psurf*std::exp(-gravity/(Rd*thref[k])*grid->z[k]);
    exner[k] = std::pow(pref[k]/psurf, Rd/cp);

    // set the base density for the entire model
    fields->rhoref[k] = pref[k] / (Rd*exner[k]*thref[k]);
  }

  for(int k=grid->kstart; k<grid->kend+1; ++k)
  {
    prefh [k] = psurf*std::exp(-gravity/(Rd*threfh[k])*grid->zh[k]);
    exnerh[k] = std::pow(prefh[k]/psurf, Rd/cp);

    // set the base density for the entire model
    fields->rhorefh[k] = prefh[k] / (Rd*exnerh[k]*threfh[k]);
  }

  // set the ghost cells for the reference variables
  // CvH for now in 2nd order
  pref [kstart-1] = 2.*prefh [kstart] - pref [kstart];
  exner[kstart-1] = 2.*exnerh[kstart] - exner[kstart];
  fields->rhoref[kstart-1] = 2.*fields->rhorefh[kstart] - fields->rhoref[kstart];

  pref [kend] = 2.*prefh [kend] - pref [kend-1];
  exner[kend] = 2.*exnerh[kend] - exner[kend-1];
  fields->rhoref[kend] = 2.*fields->rhorefh[kend] - fields->rhoref[kend-1];

  return 0;
}

int cthermo_dry::exec()
{
  if(grid->swspatialorder== "2")
    calcbuoyancytend_2nd(fields->wt->data, fields->s["th"]->data, threfh);
  else if(grid->swspatialorder == "4")
    calcbuoyancytend_4th(fields->wt->data, fields->s["th"]->data, threfh);

  return 0;
}

int cthermo_dry::getbuoyancy(cfield3d *bfield, cfield3d *tmp)
{
  calcbuoyancy(bfield->data, fields->s["th"]->data, thref);
  return 0;
}

int cthermo_dry::getN2(cfield3d *bfield, cfield3d *tmp)
{
  calcN2(bfield->data, fields->s["th"]->data, grid->dzi, thref);
  return 0;
}

int cthermo_dry::getbuoyancyfluxbot(cfield3d *bfield)
{
  calcbuoyancyfluxbot(bfield->datafluxbot, fields->s["th"]->datafluxbot, threfh);
  return 0;
}

int cthermo_dry::getbuoyancysurf(cfield3d *bfield)
{
  calcbuoyancybot(bfield->data, bfield->databot,
                  fields->s["th"]->data, fields->s["th"]->databot, thref, threfh);
  calcbuoyancyfluxbot(bfield->datafluxbot, fields->s["th"]->datafluxbot, threfh);
  return 0;
}

int cthermo_dry::calcbuoyancy(double * restrict b, double * restrict th, double * restrict thref)
{
  int ijk,jj,kk;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(int k=0; k<grid->kcells; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        b[ijk] = gravity/thref[k] * (th[ijk] - thref[k]);
      }

  return 0;
}

int cthermo_dry::calcN2(double * restrict N2, double * restrict th, double * restrict dzi, double * restrict thref)
{
  int ijk,jj,kk;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(int k=0; k<grid->kcells; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        N2[ijk] = gravity/thref[k]*0.5*(th[ijk+kk] - th[ijk-kk])*dzi[k];
      }

  return 0;
}

int cthermo_dry::calcbuoyancybot(double * restrict b , double * restrict bbot,
                                 double * restrict th, double * restrict thbot,
                                 double * restrict thref, double * restrict threfh)
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
      bbot[ij] = gravity/threfh[kstart] * (thbot[ij] - threfh[kstart]);
      b[ijk]   = gravity/thref [kstart] * (th[ijk]   - thref [kstart]);
    }

  return 0;
}

int cthermo_dry::calcbuoyancyfluxbot(double * restrict bfluxbot, double * restrict thfluxbot, double * restrict threfh)
{
  int ij,jj;
  jj = grid->icells;

  int kstart = grid->kstart;

  for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for(int i=0; i<grid->icells; ++i)
    {
      ij  = i + j*jj;
      bfluxbot[ij] = gravity/threfh[kstart]*thfluxbot[ij];
    }

  return 0;
}

int cthermo_dry::calcbuoyancytend_2nd(double * restrict wt, double * restrict th, double * restrict threfh)
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
        wt[ijk] += gravity/threfh[k] * (interp2(th[ijk-kk], th[ijk]) - threfh[k]);
      }

  return 0;
}

int cthermo_dry::calcbuoyancytend_4th(double * restrict wt, double * restrict th, double * restrict threfh)
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
        wt[ijk] += gravity/threfh[k] * (interp4(th[ijk-kk2], th[ijk-kk1], th[ijk], th[ijk+kk1]) - threfh[k]);
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

