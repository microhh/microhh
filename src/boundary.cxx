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
#include <cmath>
#include <algorithm>
#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "model.h"
#include "timeloop.h"

// boundary schemes
#include "boundary.h"
#include "boundary_surface.h"
#include "boundary_user.h"

#define NO_OFFSET 0.
#define NO_VELOCITY 0.

#define BC_DIRICHLET 0
#define BC_NEUMANN 1
#define BC_FLUX 2
#define BC_USTAR 3

cboundary::cboundary(cmodel *modelin, cinput *inputin)
{
  model  = modelin;
  grid   = model->grid;
  fields = model->fields;
  master = model->master;
}

cboundary::~cboundary()
{
  for(bcmap::const_iterator it=sbc.begin(); it!=sbc.end(); ++it)
    delete it->second;

  // empty the map
  sbc.clear();

  // clean up time dependent data
  for(std::map<std::string, double *>::const_iterator it=timedepdata.begin(); it!=timedepdata.end(); ++it)
    delete[] it->second;
}

int cboundary::processbcs(cinput *inputin)
{
  int nerror = 0;

  std::string swbot, swtop;

  nerror += inputin->getItem(&swbot, "boundary", "mbcbot", "");
  nerror += inputin->getItem(&swtop, "boundary", "mbctop", "");

  // set the bottom bc
  if(swbot == "noslip")
    mbcbot = BC_DIRICHLET;
  else if(swbot == "freeslip")
    mbcbot = BC_NEUMANN;
  else if(swbot == "ustar")
    mbcbot = BC_USTAR;
  else
  {
    if(master->mpiid == 0) std::printf("ERROR %s is illegal value for mbcbot\n", swbot.c_str());
    nerror++;
  }

  // set the top bc
  if(swtop == "noslip")
    mbctop = BC_DIRICHLET;
  else if(swtop == "freeslip")
    mbctop = BC_NEUMANN;
  else if(swtop == "ustar")
    mbctop = BC_USTAR;
  else
  {
    if(master->mpiid == 0) std::printf("ERROR %s is illegal value for mbctop\n", swtop.c_str());
    nerror++;
  }

  // read the boundaries per field
  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  {
    sbc[it->first] = new field3dbc;
    nerror += inputin->getItem(&swbot, "boundary", "sbcbot", it->first);
    nerror += inputin->getItem(&swtop, "boundary", "sbctop", it->first);
    nerror += inputin->getItem(&sbc[it->first]->bot, "boundary", "sbot", it->first);
    nerror += inputin->getItem(&sbc[it->first]->top, "boundary", "stop", it->first);

    // set the bottom bc
    if(swbot == "dirichlet")
      sbc[it->first]->bcbot = BC_DIRICHLET;
    else if(swbot == "neumann")
      sbc[it->first]->bcbot = BC_NEUMANN;
    else if(swbot == "flux")
      sbc[it->first]->bcbot = BC_FLUX;
    else
    {
      if(master->mpiid == 0) std::printf("ERROR %s is illegal value for sbcbot\n", swbot.c_str());
      nerror++;
    }

    // set the top bc
    if(swtop == "dirichlet")
      sbc[it->first]->bctop = BC_DIRICHLET;
    else if(swtop == "neumann")
      sbc[it->first]->bctop = BC_NEUMANN;
    else if(swtop == "flux")
      sbc[it->first]->bctop = BC_FLUX;
    else
    {
      if(master->mpiid == 0) std::printf("ERROR %s is illegal value for sbctop\n", swtop.c_str());
      nerror++;
    }
  }

  // get the list of time varying variables
  nerror += inputin->getItem(&swtimedep  , "boundary", "swtimedep"  , "", "0");
  nerror += inputin->getList(&timedeplist, "boundary", "timedeplist", "");

  return nerror;
}

void cboundary::init(cinput *inputin)
{
  int nerror = 0;
  nerror += processbcs(inputin);

  // there is no option (yet) for prescribing ustar without surface model
  if(mbcbot == BC_USTAR || mbctop == BC_USTAR)
  {
    if(master->mpiid == 0) std::printf("ERROR ustar bc is not supported for default boundary\n");
    ++nerror;
  }

  if(nerror)
    throw 1;
}

void cboundary::create(cinput *inputin)
{
  int nerror = 0;
  nerror += processtimedep(inputin);
  if(nerror)
    throw 1;
}

int cboundary::processtimedep(cinput *inputin)
{
  int nerror = 0;

  if(swtimedep == "1")
  {
    // create temporary list to check which entries are used
    std::vector<std::string> tmplist = timedeplist;

    // see if there is data available for the surface boundary conditions
    for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
      std::string name = "sbot[" + it->first + "]";
      if(std::find(timedeplist.begin(), timedeplist.end(), name) != timedeplist.end()) 
      {
        nerror += inputin->getTime(&timedepdata[name], &timedeptime, name);

        // remove the item from the tmplist
        std::vector<std::string>::iterator ittmp = std::find(tmplist.begin(), tmplist.end(), name);
        if(ittmp != tmplist.end())
          tmplist.erase(ittmp);
      }
    }

    // display a warning for the non-supported 
    for(std::vector<std::string>::const_iterator ittmp=tmplist.begin(); ittmp!=tmplist.end(); ++ittmp)
      if(master->mpiid == 0) std::printf("WARNING %s is not supported (yet) as a time dependent parameter\n", ittmp->c_str());
  }

  return nerror;
}

int cboundary::settimedep()
{
  if(swtimedep == "0")
    return 0;

  // first find the index for the time entries
  unsigned int index0 = 0;
  unsigned int index1 = 0;
  for(std::vector<double>::const_iterator it=timedeptime.begin(); it!=timedeptime.end(); ++it)
  {
    if(model->timeloop->time < *it)
      break;
    else
      ++index1;
  }

  // second, calculate the weighting factor
  double fac0, fac1;

  // correct for out of range situations where the simulation is longer than the time range in input
  if(index1 == 0)
  {
    fac0 = 0.;
    fac1 = 1.;
    index0 = 0;
  }
  else if(index1 == timedeptime.size())
  {
    fac0 = 1.;
    fac1 = 0.;
    index0 = index1-1;
    index1 = index0;
  }
  else
  {
    index0 = index1-1;
    double timestep;
    timestep = timedeptime[index1] - timedeptime[index0];
    fac0 = (timedeptime[index1] - model->timeloop->time) / timestep;
    fac1 = (model->timeloop->time - timedeptime[index0]) / timestep;
  }

  // process time dependent bcs for the surface fluxes
  for(fieldmap::const_iterator it1=fields->sp.begin(); it1!=fields->sp.end(); ++it1)
  {
    std::string name = "sbot[" + it1->first + "]";
    std::map<std::string, double *>::const_iterator it2 = timedepdata.find(name);
    if(it2 != timedepdata.end())
    {
      sbc[it1->first]->bot = fac0*it2->second[index0] + fac1*it2->second[index1];

// BvS: for now branched here; seems a bit wasteful to copy the entire settimedep to boundary.cu?
#ifndef USECUDA
      setbc(it1->second->databot, it1->second->datagradbot, it1->second->datafluxbot, sbc[it1->first]->bcbot, sbc[it1->first]->bot, it1->second->visc, NO_OFFSET);
#else
      setbc_g(it1->second->databot_g, it1->second->datagradbot_g, it1->second->datafluxbot_g, sbc[it1->first]->bcbot, sbc[it1->first]->bot, it1->second->visc, NO_OFFSET);
#endif
    }
  }

  return 0;
}

void cboundary::save(int iotime)
{
}

void cboundary::load(int iotime)
{
}

void cboundary::setvalues()
{
  setbc(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot, mbcbot, NO_VELOCITY, fields->visc, grid->utrans);
  setbc(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot, mbcbot, NO_VELOCITY, fields->visc, grid->vtrans);

  setbc(fields->u->datatop, fields->u->datagradtop, fields->u->datafluxtop, mbctop, NO_VELOCITY, fields->visc, grid->utrans);
  setbc(fields->v->datatop, fields->v->datagradtop, fields->v->datafluxtop, mbctop, NO_VELOCITY, fields->visc, grid->vtrans);

  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  {
    setbc(it->second->databot, it->second->datagradbot, it->second->datafluxbot, sbc[it->first]->bcbot, sbc[it->first]->bot, it->second->visc, NO_OFFSET);
    setbc(it->second->datatop, it->second->datagradtop, it->second->datafluxtop, sbc[it->first]->bctop, sbc[it->first]->top, it->second->visc, NO_OFFSET);
  }
}

#ifndef USECUDA
int cboundary::exec()
{
  // cyclic boundary conditions, do this before the bottom BC's
  grid->boundary_cyclic(fields->u->data);
  grid->boundary_cyclic(fields->v->data);
  grid->boundary_cyclic(fields->w->data);

  for(fieldmap::const_iterator it = fields->sp.begin(); it!=fields->sp.end(); ++it)
    grid->boundary_cyclic(it->second->data);

  // calculate boundary values
  bcvalues();

  if(grid->swspatialorder == "2")
  {
    setgcbot_2nd(fields->u->data, grid->dzh, mbcbot, fields->u->databot, fields->u->datagradbot);
    setgctop_2nd(fields->u->data, grid->dzh, mbctop, fields->u->datatop, fields->u->datagradtop);

    setgcbot_2nd(fields->v->data, grid->dzh, mbcbot, fields->v->databot, fields->v->datagradbot);
    setgctop_2nd(fields->v->data, grid->dzh, mbctop, fields->v->datatop, fields->v->datagradtop);

    for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
      setgcbot_2nd(it->second->data, grid->dzh, sbc[it->first]->bcbot, it->second->databot, it->second->datagradbot);
      setgctop_2nd(it->second->data, grid->dzh, sbc[it->first]->bctop, it->second->datatop, it->second->datagradtop);
    }
  }
  else if(grid->swspatialorder == "4")
  {
    setgcbot_4th(fields->u->data, grid->z, mbcbot, fields->u->databot, fields->u->datagradbot);
    setgctop_4th(fields->u->data, grid->z, mbctop, fields->u->datatop, fields->u->datagradtop);

    setgcbot_4th(fields->v->data, grid->z, mbcbot, fields->v->databot, fields->v->datagradbot);
    setgctop_4th(fields->v->data, grid->z, mbctop, fields->v->datatop, fields->v->datagradtop);

    setgcbotw_4th(fields->w->data);
    setgctopw_4th(fields->w->data);

    for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
      setgcbot_4th(it->second->data, grid->z, sbc[it->first]->bcbot, it->second->databot, it->second->datagradbot);
      setgctop_4th(it->second->data, grid->z, sbc[it->first]->bctop, it->second->datatop, it->second->datagradtop);
    }
  }

  return 0;
}
#endif

void cboundary::execcross()
{
}

int cboundary::execstats(mask *m)
{
  return 0;
}

int cboundary::bcvalues()
{
  return 0;
}

cboundary* cboundary::factory(cmaster *masterin, cinput *inputin, cmodel *modelin)
{
  std::string swboundary;
  if(inputin->getItem(&swboundary, "boundary", "swboundary", "", "default"))
    return 0;

  if(swboundary == "surface")
    return new cboundary_surface(modelin, inputin);
  else if(swboundary == "user")
    return new cboundary_user(modelin, inputin);
  else if(swboundary == "default")
    return new cboundary(modelin, inputin);
  else
  {
    masterin->printError("\"%s\" is an illegal value for swboundary\n", swboundary.c_str());
    return 0;
  }
}

int cboundary::setbc(double * restrict a, double * restrict agrad, double * restrict aflux, int sw, double aval, double visc, double offset)
{
  int ij,jj;
  jj = grid->icells;

  if(sw == BC_DIRICHLET)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij = i + j*jj;
        a[ij] = aval - offset;
      }
  }
  else if(sw == BC_NEUMANN)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij = i + j*jj;
        agrad[ij] = aval;
        aflux[ij] = -aval*visc;
      }
  }
  else if(sw == BC_FLUX)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij = i + j*jj;
        aflux[ij] = aval;
        agrad[ij] = -aval/visc;
      }
  }

  return 0;
}

// BOUNDARY CONDITIONS THAT CONTAIN A 2D PATTERN
int cboundary::setgcbot_2nd(double * restrict a, double * restrict dzh, int sw, double * restrict abot, double * restrict agradbot)
{
  int ij,ijk,jj,kk,kstart;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  kstart = grid->kstart;

  if(sw == BC_DIRICHLET)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        a[ijk-kk] = 2.*abot[ij] - a[ijk];
      }
  }
  else if(sw == BC_NEUMANN || sw == BC_FLUX)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk;
        a[ijk-kk] = -agradbot[ij]*dzh[kstart] + a[ijk];
      }
  }

  return 0;
}

int cboundary::setgctop_2nd(double * restrict a, double * restrict dzh, int sw, double * restrict atop, double * restrict agradtop)
{
  int ij,ijk,jj,kk,kend;

  kend = grid->kend;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  if(sw == BC_DIRICHLET)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + (kend-1)*kk;
        a[ijk+kk] = 2.*atop[ij] - a[ijk];
      }
  }
  else if(sw == BC_NEUMANN || sw == BC_FLUX)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + (kend-1)*kk;
        a[ijk+kk] = agradtop[ij]*dzh[kend] + a[ijk];
      }
  }

  return 0;
}

int cboundary::setgcbot_4th(double * restrict a, double * restrict z, int sw, double * restrict abot, double * restrict agradbot)
{
  int ij,ijk,jj,kk1,kk2,kstart;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  kstart = grid->kstart;

  if(sw == BC_DIRICHLET)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk1;
        a[ijk-kk1] = (8./3.)*abot[ij] - 2.*a[ijk] + (1./3.)*a[ijk+kk1];
        a[ijk-kk2] = 8.*abot[ij] - 9.*a[ijk] + 2.*a[ijk+kk1];
      }
  }
  else if(sw == BC_NEUMANN || sw == BC_FLUX)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + kstart*kk1;
        a[ijk-kk1] = -(1./24.)*grad4x(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk    ];
        a[ijk-kk2] = -(1./ 8.)*grad4x(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk+kk1];
      }
  }

  return 0;
}

int cboundary::setgctop_4th(double * restrict a, double * restrict z, int sw, double * restrict atop, double * restrict agradtop)
{
  int ij,ijk,jj,kend,kk1,kk2;

  kend = grid->kend;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  if(sw == BC_DIRICHLET)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + (kend-1)*kk1;
        a[ijk+kk1] = (8./3.)*atop[ij] - 2.*a[ijk] + (1./3.)*a[ijk-kk1];
        a[ijk+kk2] = 8.*atop[ij] - 9.*a[ijk] + 2.*a[ijk-kk1];
      }
  }
  else if(sw == BC_NEUMANN || sw == BC_FLUX)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + (kend-1)*kk1;
        a[ijk+kk1] = (1./24.)*grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk    ];
        a[ijk+kk2] = (1./ 8.)*grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk-kk1];
      }
  }

  return 0;
}

// BOUNDARY CONDITIONS FOR THE VERTICAL VELOCITY (NO PENETRATION)
int cboundary::setgcbotw_4th(double * restrict w)
{
  int ijk,jj,kk1,kk2,kstart;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  kstart = grid->kstart;

  for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for(int i=0; i<grid->icells; ++i)
    {
      ijk = i + j*jj + kstart*kk1;
      w[ijk-kk1] = -w[ijk+kk1];
      w[ijk-kk2] = -w[ijk+kk2];
    }

  return 0;
}

int cboundary::setgctopw_4th(double * restrict w)
{
  int ijk,jj,kk1,kk2,kend;

  jj  = grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  kend = grid->kend;

  for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
    for(int i=0; i<grid->icells; ++i)
    {
      ijk = i + j*jj + kend*kk1;
      w[ijk+kk1] = -w[ijk-kk1];
      w[ijk+kk2] = -w[ijk-kk2];
    }

  return 0;
}

int cboundary::prepareDevice()
{
  return 0;
}

int cboundary::forwardDevice()
{
  return 0;
}

int cboundary::backwardDevice()
{
  return 0;
}


inline double cboundary::grad4x(const double a, const double b, const double c, const double d)
{
  return (-(d-a) + 27.*(c-b));
}
