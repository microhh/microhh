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
#include "grid.h"
#include "fields.h"
#include "thermo_dry.h"
#include "defines.h"

cthermo_dry::cthermo_dry(cgrid *gridin, cfields *fieldsin, cmpi *mpiin) : cthermo(gridin, fieldsin, mpiin)
{
  swthermo = "dry";
}

cthermo_dry::~cthermo_dry()
{
}

int cthermo_dry::readinifile(cinput *inputin)
{
  int nerror = 0;
  nerror += inputin->getItem(&gravitybeta, "thermo", "gravitybeta", "");

  nerror += fields->initpfld("s");
  nerror += inputin->getItem(&fields->sp["s"]->visc, "fields", "svisc", "s");

  return nerror;
}

int cthermo_dry::exec()
{
  // extend later for gravity vector not normal to surface
  if(grid->swspatialorder== "2")
    calcbuoyancytend_2nd(fields->wt->data, fields->s["s"]->data);
  else if(grid->swspatialorder == "4")
    calcbuoyancytend_4th(fields->wt->data, fields->s["s"]->data);

  return 0;
}

int cthermo_dry::getbuoyancy(cfield3d *bfield, cfield3d *tmp)
{
  calcbuoyancy(bfield->data, fields->s["s"]->data);
  return 0;
}

int cthermo_dry::getbuoyancyfluxbot(cfield3d *bfield)
{
  calcbuoyancyfluxbot(bfield->datafluxbot, fields->s["s"]->datafluxbot);
  return 0;
}

int cthermo_dry::getbuoyancysurf(cfield3d *bfield)
{
  calcbuoyancybot(bfield->data        , bfield->databot,
                  fields->s["s"]->data, fields->s["s"]->databot);
  calcbuoyancyfluxbot(bfield->datafluxbot, fields->s["s"]->datafluxbot);
  return 0;
}

int cthermo_dry::calcbuoyancy(double * restrict b, double * restrict s)
{
  int ijk,jj,kk;
  double ql;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  double gravitybeta = this->gravitybeta;

  for(int k=0; k<grid->kcells; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        b[ijk] = gravitybeta*s[ijk];
      }

  return 0;
}

int cthermo_dry::calcbuoyancybot(double * restrict b , double * restrict bbot,
                                 double * restrict s , double * restrict sbot)
{
  int ij,ijk,jj,kk,kstart;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  kstart = grid->kstart;

  double gravitybeta = this->gravitybeta;

  for(int j=0; j<grid->jcells; j++)
#pragma ivdep
    for(int i=0; i<grid->icells; i++)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      bbot[ij] = gravitybeta*sbot[ij];
      b[ijk]   = gravitybeta*s[ijk];
    }

  return 0;
}

int cthermo_dry::calcbuoyancyfluxbot(double * restrict bfluxbot, double * restrict sfluxbot)
{
  int ij,jj,kk;
  jj = grid->icells;

  double gravitybeta = this->gravitybeta;

  for(int j=0; j<grid->jcells; j++)
#pragma ivdep
    for(int i=0; i<grid->icells; i++)
    {
      ij  = i + j*jj;
      bfluxbot[ij] = gravitybeta*sfluxbot[ij];
    }

  return 0;
}

int cthermo_dry::calcbuoyancytend_2nd(double * restrict wt, double * restrict s)
{
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  double gravitybeta = this->gravitybeta;

  // CvH check the usage of the gravity term here, in case of scaled DNS we use one. But thermal expansion coeff??
  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        wt[ijk] += gravitybeta * interp2(s[ijk-kk], s[ijk]);
      }

  return 0;
}

int cthermo_dry::calcbuoyancytend_4th(double * restrict wt, double * restrict s)
{
  int ijk,jj;
  int kk1,kk2;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  double gravitybeta = this->gravitybeta;

  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk1;
        wt[ijk] += gravitybeta * interp4(s[ijk-kk2], s[ijk-kk1], s[ijk], s[ijk+kk1]);
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

