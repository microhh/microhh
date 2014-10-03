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
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "boundary_user.h"
#include "defines.h"
#include "model.h"

Boundary_user::Boundary_user(Model *modelin, Input *inputin) : Boundary(modelin, inputin)
{
}

void Boundary_user::init(Input *inputin)
{
  int nerror = 0;

  processbcs(inputin);

  // patch type
  nerror += inputin->getItem(&patch_dim,  "boundary", "patch_dim" , "", 2 );
  nerror += inputin->getItem(&patch_xh,   "boundary", "patch_xh"  , "", 1.);
  nerror += inputin->getItem(&patch_xr,   "boundary", "patch_xr"  , "", 1.);
  nerror += inputin->getItem(&patch_xi,   "boundary", "patch_xi"  , "", 0.);
  nerror += inputin->getItem(&patch_facr, "boundary", "patch_facr", "", 1.);
  nerror += inputin->getItem(&patch_facl, "boundary", "patch_facl", "", 0.);
 
  if(nerror)
    throw 1;
}

void Boundary_user::setvalues()
{
  setbc(fields->u->databot, fields->u->datagradbot, fields->u->datafluxbot, mbcbot, noVelocity, fields->visc, grid->utrans);
  setbc(fields->v->databot, fields->v->datagradbot, fields->v->datafluxbot, mbcbot, noVelocity, fields->visc, grid->vtrans);

  setbc(fields->u->datatop, fields->u->datagradtop, fields->u->datafluxtop, mbctop, noVelocity, fields->visc, grid->utrans);
  setbc(fields->v->datatop, fields->v->datagradtop, fields->v->datafluxtop, mbctop, noVelocity, fields->visc, grid->vtrans);

  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  {
    setbc_patch(it->second->databot, it->second->datagradbot, it->second->datafluxbot,
                sbc[it->first]->bcbot, sbc[it->first]->bot, it->second->visc, noOffset, fields->s["tmp1"]->data, patch_facl, patch_facr);
    setbc      (it->second->datatop, it->second->datagradtop, it->second->datafluxtop,
                sbc[it->first]->bctop, sbc[it->first]->top, it->second->visc, noOffset);
  }
}

int Boundary_user::setbc_patch(double * restrict a, double * restrict agrad, double * restrict aflux, int sw, double aval, double visc, double offset,
                                double * restrict tmp, double facl, double facr)
{
  double avall, avalr;
  double xmod, ymod;
  double errvalx, errvaly;

  int ij,jj;
  jj = grid->icells;

  avall = facl*aval;
  avalr = facr*aval;

  // save the pattern
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij = i + j*jj;
      xmod = fmod(grid->x[i], patch_xh);
      ymod = fmod(grid->y[j], patch_xh);

      errvalx = 0.5 - 0.5*erf(2.*(std::abs(2.*xmod - patch_xh) - patch_xr) / patch_xi);

      if(patch_dim == 2)
        errvaly = 0.5 - 0.5*erf(2.*(std::abs(2.*ymod - patch_xh) - patch_xr) / patch_xi);
      else
        errvaly = 1.;

      tmp[ij] = errvalx*errvaly;
    }

  if(sw == DirichletType)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij = i + j*jj;
        a[ij] = avall + (avalr-avall)*tmp[ij] - offset;
      }
  }
  else if(sw == NeumannType)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij = i + j*jj;
        agrad[ij] = avall + (avalr-avall)*tmp[ij];
        aflux[ij] = -agrad[ij]*visc;
      }
  }
  else if(sw == FluxType)
  {
    for(int j=0; j<grid->jcells; ++j)
#pragma ivdep
      for(int i=0; i<grid->icells; ++i)
      {
        ij = i + j*jj;
        aflux[ij] = avall + (avalr-avall)*tmp[ij];
        agrad[ij] = -aflux[ij]/visc;
      }
  }

  return 0;
}
