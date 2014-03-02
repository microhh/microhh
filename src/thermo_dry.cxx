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
#include "model.h"
#include "stats.h"
#include "diff_les2s.h"

#define gravity 9.81
#define NO_OFFSET 0.

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
  nerror += inputin->getItem(&thref, "thermo", "thref", "");

  nerror += fields->initpfld("th", "Potential Temperature", "K");
  nerror += inputin->getItem(&fields->sp["th"]->visc, "fields", "svisc", "th");

  return nerror;
}

int cthermo_dry::init()
{
  stats = model->stats;
  return 0;
}

int cthermo_dry::create()
{
  stats->addprof("b", "Buoyancy", "m s-2", "z");
  for(int n=2; n<5; ++n)
  {
    std::stringstream ss;
    ss << n;
    std::string sn = ss.str();
    stats->addprof("b"+sn, "Moment " +sn+" of the buoyancy", "(m s-2)"+sn,"z");
  }

  stats->addprof("bgrad", "Gradient of the buoyancy", "m s-3", "zh");
  stats->addprof("bw"   , "Turbulent flux of the buoyancy", "m2 s-3", "zh");
  stats->addprof("bdiff", "Diffusive flux of the buoyancy", "m2 s-3", "zh");
  stats->addprof("bflux", "Total flux of the buoyancy", "m2 s-3", "zh");

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

int cthermo_dry::statsexec()
{
  // calculate the buoyancy and its surface flux for the profiles
  calcbuoyancy(fields->s["tmp1"]->data, fields->s["th"]->data);
  calcbuoyancyfluxbot(fields->s["tmp1"]->datafluxbot, fields->s["th"]->datafluxbot);

  // calculate the mean
  stats->calcmean(fields->s["tmp1"]->data, stats->profs["b"].data, NO_OFFSET);

  // calculate the moments
  for(int n=2; n<5; ++n)
  {
    std::stringstream ss;
    ss << n;
    std::string sn = ss.str();
    stats->calcmoment(fields->s["tmp1"]->data, stats->profs["b"].data, stats->profs["b"+sn].data, n, 0);
  }

  // calculate the gradients
  if(grid->swspatialorder == "2")
    stats->calcgrad_2nd(fields->s["tmp1"]->data, stats->profs["bgrad"].data, grid->dzhi);
  if(grid->swspatialorder == "4")
    stats->calcgrad_4th(fields->s["tmp1"]->data, stats->profs["bgrad"].data, grid->dzhi4);

  // calculate turbulent fluxes
  if(grid->swspatialorder == "2")
    stats->calcflux_2nd(fields->s["tmp1"]->data, fields->w->data, stats->profs["bw"].data, fields->s["tmp2"]->data, 0, 0);
  if(grid->swspatialorder == "4")
    stats->calcflux_4th(fields->s["tmp1"]->data, fields->w->data, stats->profs["bw"].data, fields->s["tmp2"]->data, 0, 0);

  // calculate diffusive fluxes
  if(model->diff->getname() == "les2s")
  {
    cdiff_les2s *diffptr = static_cast<cdiff_les2s *>(model->diff);
    stats->calcdiff_2nd(fields->s["tmp1"]->data, fields->s["evisc"]->data, stats->profs["bdiff"].data, grid->dzhi, fields->s["tmp1"]->datafluxbot, fields->s["tmp1"]->datafluxtop, diffptr->tPr);
  }
  else
    stats->calcdiff_4th(fields->s["tmp1"]->data, stats->profs["bdiff"].data, grid->dzhi4, fields->s["th"]->visc);

  // calculate the total fluxes
  stats->addfluxes(stats->profs["bflux"].data, stats->profs["bw"].data, stats->profs["bdiff"].data);

  return 0;
}

int cthermo_dry::checkthermofield(std::string name)
{
  if(name == "b")
    return 0;
  else
    return 1;
}

int cthermo_dry::getthermofield(cfield3d *field, cfield3d *tmp, std::string name)
{
  calcbuoyancy(field->data, fields->s["th"]->data);
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

int cthermo_dry::getprogvars(std::vector<std::string> *list)
{
  list->push_back("th");
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

