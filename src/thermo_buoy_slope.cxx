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
#include "thermo_buoy_slope.h"
#include "defines.h"
#include "model.h"
#include "master.h"

cthermo_buoy_slope::cthermo_buoy_slope(cmodel *modelin, cinput *inputin) : cthermo(modelin, inputin)
{
  swthermo = "slope";
  master = modelin->master;

  int nerror = 0;

  nerror += inputin->getItem(&alpha, "thermo", "alpha", "");
  nerror += inputin->getItem(&n2   , "thermo", "n2"   , "");

  nerror += fields->initpfld("b", "Buoyancy", "m s-2");
  nerror += inputin->getItem(&fields->sp["b"]->visc, "fields", "svisc", "b");

  if(grid->swspatialorder == "2")
  {
    if(master->mpiid == 0) std::printf("ERROR swthermo = buoy_slope is incompatible with swspatialorder = 2\n");
    ++nerror;
  }

  if(nerror)
    throw 1;
}

cthermo_buoy_slope::~cthermo_buoy_slope()
{
}

int cthermo_buoy_slope::exec()
{
  calcbuoyancytendu_4th(fields->ut->data, fields->s["b"]->data);
  calcbuoyancytendw_4th(fields->wt->data, fields->s["b"]->data);
  calcbuoyancytendb_4th(fields->st["b"]->data, fields->u->data, fields->w->data);

  return 0;
}

int cthermo_buoy_slope::execcross()
{
  return 0;
}

int cthermo_buoy_slope::checkthermofield(std::string name)
{
  if(name == "b")
    return 0;
  else
    return 1;
}

int cthermo_buoy_slope::getthermofield(cfield3d *field, cfield3d *tmp, std::string name)
{
  calcbuoyancy(field->data, fields->s["b"]->data);
  return 0;
}

int cthermo_buoy_slope::getbuoyancyfluxbot(cfield3d *bfield)
{
  calcbuoyancyfluxbot(bfield->datafluxbot, fields->s["b"]->datafluxbot);
  return 0;
}

int cthermo_buoy_slope::getbuoyancysurf(cfield3d *bfield)
{
  calcbuoyancybot(bfield->data        , bfield->databot,
                  fields->s["b"]->data, fields->s["b"]->databot);
  calcbuoyancyfluxbot(bfield->datafluxbot, fields->s["b"]->datafluxbot);
  return 0;
}

int cthermo_buoy_slope::getprogvars(std::vector<std::string> *list)
{
  list->push_back("b");
  return 0;
}

int cthermo_buoy_slope::calcbuoyancy(double * restrict b, double * restrict s)
{
  int ijk,jj,kk;
  double ql;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  for(int k=0; k<grid->kcells; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        b[ijk] = s[ijk];
      }

  return 0;
}

int cthermo_buoy_slope::calcbuoyancybot(double * restrict b , double * restrict bbot,
                                 double * restrict s , double * restrict sbot)
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
      bbot[ij] = sbot[ij];
      b[ijk]   = s[ijk];
    }

  return 0;
}

int cthermo_buoy_slope::calcbuoyancyfluxbot(double * restrict bfluxbot, double * restrict sfluxbot)
{
  int ij,jj,kk;
  jj = grid->icells;

  for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for(int i=0; i<grid->icells; ++i)
    {
      ij  = i + j*jj;
      bfluxbot[ij] = sfluxbot[ij];
    }

  return 0;
}

int cthermo_buoy_slope::calcbuoyancytendw_4th(double * restrict wt, double * restrict s)
{
  int ijk,jj;
  int kk1,kk2;

  jj  = grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  double alpha = this->alpha;

  for(int k=grid->kstart+1; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk1;
        wt[ijk] += std::cos(alpha) * interp4(s[ijk-kk2], s[ijk-kk1], s[ijk], s[ijk+kk1]);
      }

  return 0;
}

int cthermo_buoy_slope::calcbuoyancytendu_4th(double * restrict ut, double * restrict s)
{
  int ijk,ii1,ii2,jj,kk;

  ii1 = 1;
  ii2 = 2;
  jj  = grid->icells;
  kk  = grid->ijcells;

  double alpha = this->alpha;

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        ut[ijk] += std::sin(alpha) * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk], s[ijk+ii1]);
      }

  return 0;
}

int cthermo_buoy_slope::calcbuoyancytendb_4th(double * restrict st, double * restrict u, double * restrict w)
{
  int ijk,ii1,ii2,jj,kk1,kk2;

  ii1 = 1;
  ii2 = 2;
  jj  = 1*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  double alpha = this->alpha;
  double n2    = this->n2;

  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk1;
        st[ijk] -= ( n2*std::sin(alpha)*interp4(u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2])
                   + n2*std::cos(alpha)*interp4(w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2]) );
      }

  return 0;
}

inline double cthermo_buoy_slope::interp4(const double a, const double b, const double c, const double d)
{
  return (-a + 9.*b + 9.*c - d) / 16.;
}

