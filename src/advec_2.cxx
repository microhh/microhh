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
#include "grid.h"
#include "fields.h"
#include "advec_2.h"
#include "defines.h"
#include "constants.h"
#include "fd.h"
#include "model.h"

using namespace fd::o2;

cadvec_2::cadvec_2(cmodel *modelin, cinput *inputin) : cadvec(modelin, inputin)
{
}

cadvec_2::~cadvec_2()
{
}

double cadvec_2::getcfl(double dt)
{
  double cfl;

  cfl = calccfl(fields->u->data, fields->v->data, fields->w->data, grid->dzi, dt);

  return cfl;
}

unsigned long cadvec_2::gettimelim(unsigned long idt, double dt)
{
  unsigned long idtlim;
  double cfl;

  cfl = calccfl(fields->u->data, fields->v->data, fields->w->data, grid->dzi, dt);
  // avoid zero divisons
  cfl = std::max(constants::dsmall, cfl);

  idtlim = idt * cflmax / cfl;

  return idtlim;
}

#ifndef USECUDA
void cadvec_2::exec()
{
  advecu(fields->ut->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi,
         fields->rhoref, fields->rhorefh);
  advecv(fields->vt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi,
         fields->rhoref, fields->rhorefh);
  advecw(fields->wt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzhi,
         fields->rhoref, fields->rhorefh);

  for(fieldmap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    advecs(it->second->data, fields->s[it->first]->data, fields->u->data, fields->v->data, fields->w->data,
           grid->dzi, fields->rhoref, fields->rhorefh);
}
#endif

#ifndef USECUDA
double cadvec_2::calccfl(double * restrict u, double * restrict v, double * restrict w, double * restrict dzi, double dt)
{
  int    ijk,ii,jj,kk;
  double dxi,dyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  double cfl = 0;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk;
        cfl = std::max(cfl, std::abs(interp2(u[ijk], u[ijk+ii]))*dxi + std::abs(interp2(v[ijk], v[ijk+jj]))*dyi + std::abs(interp2(w[ijk], w[ijk+kk]))*dzi[k]);
      }

  grid->getmax(&cfl);

  cfl = cfl*dt;

  return cfl;
}
#endif

void cadvec_2::advecu(double * restrict ut, double * restrict u, double * restrict v, double * restrict w,
                      double * restrict dzi, double * restrict rhoref, double * restrict rhorefh)
{
  int    ijk,ii,jj,kk;
  double dxi,dyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ut[ijk] +=
              - (  interp2(u[ijk   ], u[ijk+ii]) * interp2(u[ijk   ], u[ijk+ii])
                 - interp2(u[ijk-ii], u[ijk   ]) * interp2(u[ijk-ii], u[ijk   ]) ) * dxi

              - (  interp2(v[ijk-ii+jj], v[ijk+jj]) * interp2(u[ijk   ], u[ijk+jj])
                 - interp2(v[ijk-ii   ], v[ijk   ]) * interp2(u[ijk-jj], u[ijk   ]) ) * dyi

              - (  rhorefh[k+1] * interp2(w[ijk-ii+kk], w[ijk+kk]) * interp2(u[ijk   ], u[ijk+kk])
                 - rhorefh[k  ] * interp2(w[ijk-ii   ], w[ijk   ]) * interp2(u[ijk-kk], u[ijk   ]) ) / rhoref[k] * dzi[k];
      }
}

void cadvec_2::advecv(double * restrict vt, double * restrict u, double * restrict v, double * restrict w,
                      double * restrict dzi, double * restrict rhoref, double * restrict rhorefh)
{
  int    ijk,ii,jj,kk;
  double dxi,dyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        vt[ijk] +=
              - (  interp2(u[ijk+ii-jj], u[ijk+ii]) * interp2(v[ijk   ], v[ijk+ii])
                 - interp2(u[ijk   -jj], u[ijk   ]) * interp2(v[ijk-ii], v[ijk   ]) ) * dxi

              - (  interp2(v[ijk   ], v[ijk+jj]) * interp2(v[ijk   ], v[ijk+jj])
                 - interp2(v[ijk-jj], v[ijk   ]) * interp2(v[ijk-jj], v[ijk   ]) ) * dyi

              - (  rhorefh[k+1] * interp2(w[ijk-jj+kk], w[ijk+kk]) * interp2(v[ijk   ], v[ijk+kk])
                 - rhorefh[k  ] * interp2(w[ijk-jj   ], w[ijk   ]) * interp2(v[ijk-kk], v[ijk   ]) ) / rhoref[k] * dzi[k];
      }
}

void cadvec_2::advecw(double * restrict wt, double * restrict u, double * restrict v, double * restrict w,
                      double * restrict dzhi, double * restrict rhoref, double * restrict rhorefh)
{
  int    ijk,ii,jj,kk;
  double dxi,dyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        wt[ijk] +=
              - (  interp2(u[ijk+ii-kk], u[ijk+ii]) * interp2(w[ijk   ], w[ijk+ii])
                 - interp2(u[ijk   -kk], u[ijk   ]) * interp2(w[ijk-ii], w[ijk   ]) ) * dxi

              - (  interp2(v[ijk+jj-kk], v[ijk+jj]) * interp2(w[ijk   ], w[ijk+jj])
                 - interp2(v[ijk   -kk], v[ijk   ]) * interp2(w[ijk-jj], w[ijk   ]) ) * dyi

              - (  rhoref[k  ] * interp2(w[ijk   ], w[ijk+kk]) * interp2(w[ijk   ], w[ijk+kk])
                 - rhoref[k-1] * interp2(w[ijk-kk], w[ijk   ]) * interp2(w[ijk-kk], w[ijk   ]) ) / rhorefh[k] * dzhi[k];
      }
}

void cadvec_2::advecs(double * restrict st, double * restrict s, double * restrict u, double * restrict v, double * restrict w,
                      double * restrict dzi, double * restrict rhoref, double * restrict rhorefh)
{
  int    ijk,ii,jj,kk;
  double dxi,dyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        st[ijk] +=
              - (  u[ijk+ii] * interp2(s[ijk   ], s[ijk+ii])
                 - u[ijk   ] * interp2(s[ijk-ii], s[ijk   ]) ) * dxi

              - (  v[ijk+jj] * interp2(s[ijk   ], s[ijk+jj])
                 - v[ijk   ] * interp2(s[ijk-jj], s[ijk   ]) ) * dyi

              - (  rhorefh[k+1] * w[ijk+kk] * interp2(s[ijk   ], s[ijk+kk])
                 - rhorefh[k  ] * w[ijk   ] * interp2(s[ijk-kk], s[ijk   ]) ) / rhoref[k] * dzi[k];
      }
}
