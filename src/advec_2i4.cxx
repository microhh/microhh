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
#include "advec_2i4.h"
#include "defines.h"
#include "constants.h"
#include "fd.h"
#include "model.h"

using namespace fd::o4;
using namespace fd::o2;

Advec2i4::Advec2i4(Model *modelin, Input *inputin) : Advec(modelin, inputin)
{
  // BvS Where to put this?
  grid->igc = std::max(grid->igc,2);
  grid->jgc = std::max(grid->jgc,2);
}

Advec2i4::~Advec2i4()
{
}

unsigned long Advec2i4::getTimeLimit(unsigned long idt, double dt)
{
  unsigned long idtlim;
  double cfl;

  cfl = calc_cfl(fields->u->data, fields->v->data, fields->w->data, grid->dzi, dt);
  // avoid zero divisons
  cfl = std::max(constants::dsmall, cfl);
  idtlim = idt * cflmax / cfl;

  return idtlim;
}

double Advec2i4::get_cfl(double dt)
{
  double cfl;

  cfl = calc_cfl(fields->u->data, fields->v->data, fields->w->data, grid->dzi, dt);

  return cfl;
}

void Advec2i4::exec()
{
  advecu(fields->ut->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi );
  advecv(fields->vt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi );
  advecw(fields->wt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzhi);

  for(FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    advecs(it->second->data, fields->sp[it->first]->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi);
}

double Advec2i4::calc_cfl(double * restrict u, double * restrict v, double * restrict w, double * restrict dzi, double dt)
{
  int    ijk,ii1,jj1,kk1,ii2,jj2,kk2,kk3,k;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;
  kk3 = 2*grid->ijcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  double cfl = 0;

  int kstart = grid->kstart;
  int kend   = grid->kend;

  k = kstart;
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk  = i + j*jj1 + k*kk1;
      cfl = std::max(cfl, std::abs(interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]))*dxi 
                        + std::abs(interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]))*dyi 
                        + std::abs(interp2(w[ijk    ], w[ijk+kk1]))*dzi[k]);
    }

  for(k=grid->kstart+1; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj1 + k*kk1;
        cfl = std::max(cfl, std::abs(interp4(u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2]))*dxi 
                          + std::abs(interp4(v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2]))*dyi 
                          + std::abs(interp4(w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2]))*dzi[k]);
      }

  k = kend-1;
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk  = i + j*jj1 + k*kk1;
      cfl = std::max(cfl, std::abs(interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]))*dxi 
                        + std::abs(interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]))*dyi 
                        + std::abs(interp2(w[ijk    ], w[ijk+kk1]))*dzi[k]);
    }

  grid->getMax(&cfl);

  cfl = cfl*dt;

  return cfl;
}

void Advec2i4::advecu(double * restrict ut, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi)
{
  int    ijk,ii1,jj1,kk1,ii2,jj2,kk2,k;
  double dxi,dyi;
  int    kstart, kend;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  kstart = grid->kstart;
  kend   = grid->kend;

  k = kstart; 
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + k*kk1;
      ut[ijk] += 
            - (  interp4(u[ijk-ii1],     u[ijk    ],     u[ijk+ii1], u[ijk+ii2])     * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
               - interp4(u[ijk-ii2],     u[ijk-ii1],     u[ijk    ], u[ijk+ii1])     * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

            - (  interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
               - interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 

            - (  interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp2(u[ijk], u[ijk+kk1]) ) * dzi[k];
    }

  k = kstart + 1; 
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + k*kk1;
      ut[ijk] += 
            - (  interp4(u[ijk-ii1],     u[ijk    ],     u[ijk+ii1], u[ijk+ii2])     * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
               - interp4(u[ijk-ii2],     u[ijk-ii1],     u[ijk    ], u[ijk+ii1])     * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

            - (  interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
               - interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 

            - (  interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp4(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2])
               - interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp2(u[ijk-kk1], u[ijk    ]) ) * dzi[k];
    }


  for(k=grid->kstart+2; k<grid->kend-2; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        ut[ijk] += 
              - (  interp4(u[ijk-ii1],     u[ijk    ],     u[ijk+ii1], u[ijk+ii2])     * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                 - interp4(u[ijk-ii2],     u[ijk-ii1],     u[ijk    ], u[ijk+ii1])     * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

              - (  interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                 - interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 

              - (  interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp4(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2])
                 - interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp4(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) * dzi[k];
      }

  k = kend - 2; 
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + k*kk1;
      ut[ijk] += 
            - (  interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]) * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
               - interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

            - (  interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
               - interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 

            - (  interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp2(u[ijk    ], u[ijk+kk1])
               - interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp4(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) * dzi[k];
    }

  k = kend - 1; 
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + k*kk1;
      ut[ijk] += 
            - (  interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]) * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
               - interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

            - (  interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
               - interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 

            - (- interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp2(u[ijk-kk1], u[ijk    ]) ) * dzi[k];
    }
}

void Advec2i4::advecv(double * restrict vt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi)
{
  int    ijk,ii1,jj1,kk1,ii2,jj2,kk2,k;
  double dxi,dyi;
  int    kstart, kend;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  kstart = grid->kstart;
  kend   = grid->kend;

  k = kstart;
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + k*kk1;
      vt[ijk] += 
            - (  interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
               - interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

            - (  interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
               - interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

            - (  interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp2(v[ijk    ], v[ijk+kk1]) ) * dzi[k];
    }

  k = kstart+1;
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + k*kk1;
      vt[ijk] += 
            - (  interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
               - interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

            - (  interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
               - interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

            - (  interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp4(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])
               - interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp2(v[ijk-kk1], v[ijk    ]) ) * dzi[k];
    }

  for(k=grid->kstart+2; k<grid->kend-2; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        vt[ijk] += 
              - (  interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                 - interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

              - (  interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                 - interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

              - (  interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp4(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])
                 - interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp4(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) * dzi[k];
      }

  k = kend-2;
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + k*kk1;
      vt[ijk] += 
            - (  interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
               - interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

            - (  interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
               - interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

            - (  interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp2(v[ijk    ], v[ijk+kk1])
               - interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp4(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) * dzi[k];
    }

  k = kend-1;
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + k*kk1;
      vt[ijk] += 
            - (  interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
               - interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

            - (  interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
               - interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

            - (- interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp2(v[ijk-kk1], v[ijk    ]) ) * dzi[k];
    }
}

void Advec2i4::advecw(double * restrict wt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzhi)
{
  int    ijk,ii1,jj1,kk1,ii2,jj2,kk2,k;
  double dxi,dyi;
  int    kstart, kend;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  kstart = grid->kstart;
  kend   = grid->kend;

  k = kstart+1;
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + k*kk1;
      wt[ijk] += 
            - (  interp2(u[ijk+ii1-kk1], u[ijk+ii1]) * interp4(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
               - interp2(u[ijk    -kk1], u[ijk    ]) * interp4(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi

            - (  interp2(v[ijk+jj1-kk1], v[ijk+jj1]) * interp4(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
               - interp2(v[ijk    -kk1], v[ijk    ]) * interp4(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi

            - (  interp4(w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2]) * interp4   (w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2])
               - interp2(w[ijk], w[ijk-kk1]) * interp2(w[ijk], w[ijk-kk1]) ) * dzhi[k];
    }

  for(k=grid->kstart+2; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        wt[ijk] += 
              - (  interp4(u[ijk+ii1-kk2], u[ijk+ii1-kk1], u[ijk+ii1], u[ijk+ii1+kk1]) * interp4(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
                 - interp4(u[ijk    -kk2], u[ijk    -kk1], u[ijk    ], u[ijk    +kk1]) * interp4(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi

              - (  interp4(v[ijk+jj1-kk2], v[ijk+jj1-kk1], v[ijk+jj1], v[ijk+jj1+kk1]) * interp4(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
                 - interp4(v[ijk    -kk2], v[ijk    -kk1], v[ijk    ], v[ijk    +kk1]) * interp4(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi

              - (  interp4(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) * interp4(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2])
                 - interp4(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) * interp4(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) * dzhi[k];
      }

  k = kend-1;
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + k*kk1;
      wt[ijk] += 
            - (  interp2(u[ijk+ii1-kk1], u[ijk+ii1]) * interp4(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
               - interp2(u[ijk    -kk1], u[ijk    ]) * interp4(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi

            - (  interp2(v[ijk+jj1-kk1], v[ijk+jj1]) * interp4(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
               - interp2(v[ijk    -kk1], v[ijk    ]) * interp4(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi

            - (  interp2(w[ijk+kk1], w[ijk]) * interp2(w[ijk+kk1], w[ijk])
               - interp4(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) * interp4   (w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) * dzhi[k];
    }
}

void Advec2i4::advecs(double * restrict st, double * restrict s, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi)
{
  int    ijk,ii1,jj1,kk1,ii2,jj2,kk2,k;
  double dxi,dyi;
  int    kstart, kend;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  kstart = grid->kstart;
  kend   = grid->kend;
 
  // assume that w at the boundary equals zero...
  k = kstart;
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + k*kk1;
      st[ijk] += 
            - (  u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
               - u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

            - (  v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
               - v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 

            - (  w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1]) ) * dzi[k];
    }

  k = kstart+1;
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + k*kk1;
      st[ijk] += 
            - (  u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
               - u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

            - (  v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
               - v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 

            - (  w[ijk+kk1] * interp4(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
               - w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]) ) * dzi[k];
    }

  for(k=grid->kstart+2; k<grid->kend-2; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        st[ijk] += 
              - (  u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                 - u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

              - (  v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                 - v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 

              - (  w[ijk+kk1] * interp4(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                 - w[ijk    ] * interp4(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) * dzi[k];
      }

  k = kend-2;
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + k*kk1;
      st[ijk] += 
            - (  u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
               - u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

            - (  v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
               - v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 

            - (  w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1])
               - w[ijk    ] * interp4(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) * dzi[k];
    }

  // assume that w at the boundary equals zero...
  k = kend-1;
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + k*kk1;
      st[ijk] += 
            - (  u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
               - u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

            - (  v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
               - v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 

            - (- w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]) ) * dzi[k];
    }
}
