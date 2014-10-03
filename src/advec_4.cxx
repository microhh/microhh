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
#include "advec_4.h"
#include "defines.h"
#include "constants.h"
#include "fd.h"
#include "model.h"

using namespace fd::o4;

Advec_4::Advec_4(Model *modelin, Input *inputin) : Advec(modelin, inputin)
{
}

Advec_4::~Advec_4()
{
}

unsigned long Advec_4::gettimelim(unsigned long idt, double dt)
{
  unsigned long idtlim;
  double cfl;

  cfl = calccfl(fields->u->data, fields->v->data, fields->w->data, grid->dzi, dt);
  // avoid zero divisons
  cfl = std::max(constants::dsmall, cfl);
  idtlim = idt * cflmax / cfl;

  return idtlim;
}

double Advec_4::getcfl(double dt)
{
  double cfl;

  cfl = calccfl((*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi, dt);

  return cfl;
}

#ifndef USECUDA
void Advec_4::exec()
{
  advecu(fields->ut->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi4 );
  advecv(fields->vt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi4 );
  advecw(fields->wt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzhi4);

  for(fieldmap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    advecs(it->second->data, fields->sp[it->first]->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi4);
}
#endif

#ifndef USECUDA
double Advec_4::calccfl(double * restrict u, double * restrict v, double * restrict w, double * restrict dzi, double dt)
{
  int    ijk;
  int    ii1,ii2,jj1,jj2,kk1,kk2;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  double cfl = 0;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj1 + k*kk1;
        cfl = std::max(cfl, std::abs(ci0*u[ijk-ii1] + ci1*u[ijk] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2])*dxi 
                          + std::abs(ci0*v[ijk-jj1] + ci1*v[ijk] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2])*dyi 
                          + std::abs(ci0*w[ijk-kk1] + ci1*w[ijk] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])*dzi[k] );
      }

  grid->getmax(&cfl);

  cfl = cfl*dt;

  return cfl;
}
#endif

void Advec_4::advecu(double * restrict ut, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi4)
{
  int    ijk,kstart,kend;
  int    ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

  kstart = grid->kstart;
  kend   = grid->kend;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + kstart*kk1;
      ut[ijk] -= ( cg0*((ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]) * (ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]))
                 + cg1*((ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]) * (ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]))
                 + cg2*((ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]) * (ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]))
                 + cg3*((ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3]) * (ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3])) ) * cgi*dxi;

      ut[ijk] -= ( cg0*((ci0*v[ijk-ii2-jj1] + ci1*v[ijk-ii1-jj1] + ci2*v[ijk-jj1] + ci3*v[ijk+ii1-jj1]) * (ci0*u[ijk-jj3] + ci1*u[ijk-jj2] + ci2*u[ijk-jj1] + ci3*u[ijk    ]))
                 + cg1*((ci0*v[ijk-ii2    ] + ci1*v[ijk-ii1    ] + ci2*v[ijk    ] + ci3*v[ijk+ii1    ]) * (ci0*u[ijk-jj2] + ci1*u[ijk-jj1] + ci2*u[ijk    ] + ci3*u[ijk+jj1]))
                 + cg2*((ci0*v[ijk-ii2+jj1] + ci1*v[ijk-ii1+jj1] + ci2*v[ijk+jj1] + ci3*v[ijk+ii1+jj1]) * (ci0*u[ijk-jj1] + ci1*u[ijk    ] + ci2*u[ijk+jj1] + ci3*u[ijk+jj2]))
                 + cg3*((ci0*v[ijk-ii2+jj2] + ci1*v[ijk-ii1+jj2] + ci2*v[ijk+jj2] + ci3*v[ijk+ii1+jj2]) * (ci0*u[ijk    ] + ci1*u[ijk+jj1] + ci2*u[ijk+jj2] + ci3*u[ijk+jj3])) ) * cgi*dyi;

      ut[ijk] -= ( cg0*((ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk-kk1] + ci3*w[ijk+ii1-kk1]) * (bi0*u[ijk-kk2] + bi1*u[ijk-kk1] + bi2*u[ijk    ] + bi3*u[ijk+kk1]))
                 + cg1*((ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk    ] + ci3*w[ijk+ii1    ]) * (ci0*u[ijk-kk2] + ci1*u[ijk-kk1] + ci2*u[ijk    ] + ci3*u[ijk+kk1]))
                 + cg2*((ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk+kk1] + ci3*w[ijk+ii1+kk1]) * (ci0*u[ijk-kk1] + ci1*u[ijk    ] + ci2*u[ijk+kk1] + ci3*u[ijk+kk2]))
                 + cg3*((ci0*w[ijk-ii2+kk2] + ci1*w[ijk-ii1+kk2] + ci2*w[ijk+kk2] + ci3*w[ijk+ii1+kk2]) * (ci0*u[ijk    ] + ci1*u[ijk+kk1] + ci2*u[ijk+kk2] + ci3*u[ijk+kk3])) )
                 * dzi4[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        ut[ijk] -= ( cg0*((ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]) * (ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]))
                   + cg1*((ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]) * (ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]))
                   + cg2*((ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]) * (ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]))
                   + cg3*((ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3]) * (ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3])) ) * cgi*dxi;

        ut[ijk] -= ( cg0*((ci0*v[ijk-ii2-jj1] + ci1*v[ijk-ii1-jj1] + ci2*v[ijk-jj1] + ci3*v[ijk+ii1-jj1]) * (ci0*u[ijk-jj3] + ci1*u[ijk-jj2] + ci2*u[ijk-jj1] + ci3*u[ijk    ]))
                   + cg1*((ci0*v[ijk-ii2    ] + ci1*v[ijk-ii1    ] + ci2*v[ijk    ] + ci3*v[ijk+ii1    ]) * (ci0*u[ijk-jj2] + ci1*u[ijk-jj1] + ci2*u[ijk    ] + ci3*u[ijk+jj1]))
                   + cg2*((ci0*v[ijk-ii2+jj1] + ci1*v[ijk-ii1+jj1] + ci2*v[ijk+jj1] + ci3*v[ijk+ii1+jj1]) * (ci0*u[ijk-jj1] + ci1*u[ijk    ] + ci2*u[ijk+jj1] + ci3*u[ijk+jj2]))
                   + cg3*((ci0*v[ijk-ii2+jj2] + ci1*v[ijk-ii1+jj2] + ci2*v[ijk+jj2] + ci3*v[ijk+ii1+jj2]) * (ci0*u[ijk    ] + ci1*u[ijk+jj1] + ci2*u[ijk+jj2] + ci3*u[ijk+jj3])) ) * cgi*dyi;

        ut[ijk] -= ( cg0*((ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk-kk1] + ci3*w[ijk+ii1-kk1]) * (ci0*u[ijk-kk3] + ci1*u[ijk-kk2] + ci2*u[ijk-kk1] + ci3*u[ijk    ]))
                   + cg1*((ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk    ] + ci3*w[ijk+ii1    ]) * (ci0*u[ijk-kk2] + ci1*u[ijk-kk1] + ci2*u[ijk    ] + ci3*u[ijk+kk1]))
                   + cg2*((ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk+kk1] + ci3*w[ijk+ii1+kk1]) * (ci0*u[ijk-kk1] + ci1*u[ijk    ] + ci2*u[ijk+kk1] + ci3*u[ijk+kk2]))
                   + cg3*((ci0*w[ijk-ii2+kk2] + ci1*w[ijk-ii1+kk2] + ci2*w[ijk+kk2] + ci3*w[ijk+ii1+kk2]) * (ci0*u[ijk    ] + ci1*u[ijk+kk1] + ci2*u[ijk+kk2] + ci3*u[ijk+kk3])) )
                   * dzi4[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      ut[ijk] -= ( cg0*((ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]) * (ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]))
                 + cg1*((ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]) * (ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]))
                 + cg2*((ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]) * (ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]))
                 + cg3*((ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3]) * (ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3])) ) * cgi*dxi;

      ut[ijk] -= ( cg0*((ci0*v[ijk-ii2-jj1] + ci1*v[ijk-ii1-jj1] + ci2*v[ijk-jj1] + ci3*v[ijk+ii1-jj1]) * (ci0*u[ijk-jj3] + ci1*u[ijk-jj2] + ci2*u[ijk-jj1] + ci3*u[ijk    ]))
                 + cg1*((ci0*v[ijk-ii2    ] + ci1*v[ijk-ii1    ] + ci2*v[ijk    ] + ci3*v[ijk+ii1    ]) * (ci0*u[ijk-jj2] + ci1*u[ijk-jj1] + ci2*u[ijk    ] + ci3*u[ijk+jj1]))
                 + cg2*((ci0*v[ijk-ii2+jj1] + ci1*v[ijk-ii1+jj1] + ci2*v[ijk+jj1] + ci3*v[ijk+ii1+jj1]) * (ci0*u[ijk-jj1] + ci1*u[ijk    ] + ci2*u[ijk+jj1] + ci3*u[ijk+jj2]))
                 + cg3*((ci0*v[ijk-ii2+jj2] + ci1*v[ijk-ii1+jj2] + ci2*v[ijk+jj2] + ci3*v[ijk+ii1+jj2]) * (ci0*u[ijk    ] + ci1*u[ijk+jj1] + ci2*u[ijk+jj2] + ci3*u[ijk+jj3])) ) * cgi*dyi;

      ut[ijk] -= ( cg0*((ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk-kk1] + ci3*w[ijk+ii1-kk1]) * (ci0*u[ijk-kk3] + ci1*u[ijk-kk2] + ci2*u[ijk-kk1] + ci3*u[ijk    ]))
                 + cg1*((ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk    ] + ci3*w[ijk+ii1    ]) * (ci0*u[ijk-kk2] + ci1*u[ijk-kk1] + ci2*u[ijk    ] + ci3*u[ijk+kk1]))
                 + cg2*((ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk+kk1] + ci3*w[ijk+ii1+kk1]) * (ci0*u[ijk-kk1] + ci1*u[ijk    ] + ci2*u[ijk+kk1] + ci3*u[ijk+kk2]))
                 + cg3*((ci0*w[ijk-ii2+kk2] + ci1*w[ijk-ii1+kk2] + ci2*w[ijk+kk2] + ci3*w[ijk+ii1+kk2]) * (ti0*u[ijk-kk1] + ti1*u[ijk    ] + ti2*u[ijk+kk1] + ti3*u[ijk+kk2])) )
                 * dzi4[kend-1];
    }
}

void Advec_4::advecv(double * restrict vt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi4)
{
  int    ijk,kstart,kend;
  int    ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

  kstart = grid->kstart;
  kend   = grid->kend;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + kstart*kk1;
      vt[ijk] -= ( cg0*((ci0*u[ijk-ii1-jj2] + ci1*u[ijk-ii1-jj1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+jj1]) * (ci0*v[ijk-ii3] + ci1*v[ijk-ii2] + ci2*v[ijk-ii1] + ci3*v[ijk    ]))
                 + cg1*((ci0*u[ijk    -jj2] + ci1*u[ijk    -jj1] + ci2*u[ijk    ] + ci3*u[ijk    +jj1]) * (ci0*v[ijk-ii2] + ci1*v[ijk-ii1] + ci2*v[ijk    ] + ci3*v[ijk+ii1]))
                 + cg2*((ci0*u[ijk+ii1-jj2] + ci1*u[ijk+ii1-jj1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+jj1]) * (ci0*v[ijk-ii1] + ci1*v[ijk    ] + ci2*v[ijk+ii1] + ci3*v[ijk+ii2]))
                 + cg3*((ci0*u[ijk+ii2-jj2] + ci1*u[ijk+ii2-jj1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+jj1]) * (ci0*v[ijk    ] + ci1*v[ijk+ii1] + ci2*v[ijk+ii2] + ci3*v[ijk+ii3])) ) * cgi*dxi;

      vt[ijk] -= ( cg0*((ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]) * (ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]))
                 + cg1*((ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]) * (ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]))
                 + cg2*((ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]) * (ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]))
                 + cg3*((ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3]) * (ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3])) ) * cgi*dyi;

      vt[ijk] -= ( cg0*((ci0*w[ijk-jj2-kk1] + ci1*w[ijk-jj1-kk1] + ci2*w[ijk-kk1] + ci3*w[ijk+jj1-kk1]) * (bi0*v[ijk-kk2] + bi1*v[ijk-kk1] + bi2*v[ijk    ] + bi3*v[ijk+kk1]))
                 + cg1*((ci0*w[ijk-jj2    ] + ci1*w[ijk-jj1    ] + ci2*w[ijk    ] + ci3*w[ijk+jj1    ]) * (ci0*v[ijk-kk2] + ci1*v[ijk-kk1] + ci2*v[ijk    ] + ci3*v[ijk+kk1]))
                 + cg2*((ci0*w[ijk-jj2+kk1] + ci1*w[ijk-jj1+kk1] + ci2*w[ijk+kk1] + ci3*w[ijk+jj1+kk1]) * (ci0*v[ijk-kk1] + ci1*v[ijk    ] + ci2*v[ijk+kk1] + ci3*v[ijk+kk2]))
                 + cg3*((ci0*w[ijk-jj2+kk2] + ci1*w[ijk-jj1+kk2] + ci2*w[ijk+kk2] + ci3*w[ijk+jj1+kk2]) * (ci0*v[ijk    ] + ci1*v[ijk+kk1] + ci2*v[ijk+kk2] + ci3*v[ijk+kk3])) )
                 * dzi4[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        vt[ijk] -= ( cg0*((ci0*u[ijk-ii1-jj2] + ci1*u[ijk-ii1-jj1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+jj1]) * (ci0*v[ijk-ii3] + ci1*v[ijk-ii2] + ci2*v[ijk-ii1] + ci3*v[ijk    ]))
                   + cg1*((ci0*u[ijk    -jj2] + ci1*u[ijk    -jj1] + ci2*u[ijk    ] + ci3*u[ijk    +jj1]) * (ci0*v[ijk-ii2] + ci1*v[ijk-ii1] + ci2*v[ijk    ] + ci3*v[ijk+ii1]))
                   + cg2*((ci0*u[ijk+ii1-jj2] + ci1*u[ijk+ii1-jj1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+jj1]) * (ci0*v[ijk-ii1] + ci1*v[ijk    ] + ci2*v[ijk+ii1] + ci3*v[ijk+ii2]))
                   + cg3*((ci0*u[ijk+ii2-jj2] + ci1*u[ijk+ii2-jj1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+jj1]) * (ci0*v[ijk    ] + ci1*v[ijk+ii1] + ci2*v[ijk+ii2] + ci3*v[ijk+ii3])) ) * cgi*dxi;

        vt[ijk] -= ( cg0*((ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]) * (ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]))
                   + cg1*((ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]) * (ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]))
                   + cg2*((ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]) * (ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]))
                   + cg3*((ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3]) * (ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3])) ) * cgi*dyi;

        vt[ijk] -= ( cg0*((ci0*w[ijk-jj2-kk1] + ci1*w[ijk-jj1-kk1] + ci2*w[ijk-kk1] + ci3*w[ijk+jj1-kk1]) * (ci0*v[ijk-kk3] + ci1*v[ijk-kk2] + ci2*v[ijk-kk1] + ci3*v[ijk    ]))
                   + cg1*((ci0*w[ijk-jj2    ] + ci1*w[ijk-jj1    ] + ci2*w[ijk    ] + ci3*w[ijk+jj1    ]) * (ci0*v[ijk-kk2] + ci1*v[ijk-kk1] + ci2*v[ijk    ] + ci3*v[ijk+kk1]))
                   + cg2*((ci0*w[ijk-jj2+kk1] + ci1*w[ijk-jj1+kk1] + ci2*w[ijk+kk1] + ci3*w[ijk+jj1+kk1]) * (ci0*v[ijk-kk1] + ci1*v[ijk    ] + ci2*v[ijk+kk1] + ci3*v[ijk+kk2]))
                   + cg3*((ci0*w[ijk-jj2+kk2] + ci1*w[ijk-jj1+kk2] + ci2*w[ijk+kk2] + ci3*w[ijk+jj1+kk2]) * (ci0*v[ijk    ] + ci1*v[ijk+kk1] + ci2*v[ijk+kk2] + ci3*v[ijk+kk3])) )
                   * dzi4[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      vt[ijk] -= ( cg0*((ci0*u[ijk-ii1-jj2] + ci1*u[ijk-ii1-jj1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+jj1]) * (ci0*v[ijk-ii3] + ci1*v[ijk-ii2] + ci2*v[ijk-ii1] + ci3*v[ijk    ]))
                 + cg1*((ci0*u[ijk    -jj2] + ci1*u[ijk    -jj1] + ci2*u[ijk    ] + ci3*u[ijk    +jj1]) * (ci0*v[ijk-ii2] + ci1*v[ijk-ii1] + ci2*v[ijk    ] + ci3*v[ijk+ii1]))
                 + cg2*((ci0*u[ijk+ii1-jj2] + ci1*u[ijk+ii1-jj1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+jj1]) * (ci0*v[ijk-ii1] + ci1*v[ijk    ] + ci2*v[ijk+ii1] + ci3*v[ijk+ii2]))
                 + cg3*((ci0*u[ijk+ii2-jj2] + ci1*u[ijk+ii2-jj1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+jj1]) * (ci0*v[ijk    ] + ci1*v[ijk+ii1] + ci2*v[ijk+ii2] + ci3*v[ijk+ii3])) ) * cgi*dxi;

      vt[ijk] -= ( cg0*((ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]) * (ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]))
                 + cg1*((ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]) * (ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]))
                 + cg2*((ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]) * (ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]))
                 + cg3*((ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3]) * (ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3])) ) * cgi*dyi;

      vt[ijk] -= ( cg0*((ci0*w[ijk-jj2-kk1] + ci1*w[ijk-jj1-kk1] + ci2*w[ijk-kk1] + ci3*w[ijk+jj1-kk1]) * (ci0*v[ijk-kk3] + ci1*v[ijk-kk2] + ci2*v[ijk-kk1] + ci3*v[ijk    ]))
                 + cg1*((ci0*w[ijk-jj2    ] + ci1*w[ijk-jj1    ] + ci2*w[ijk    ] + ci3*w[ijk+jj1    ]) * (ci0*v[ijk-kk2] + ci1*v[ijk-kk1] + ci2*v[ijk    ] + ci3*v[ijk+kk1]))
                 + cg2*((ci0*w[ijk-jj2+kk1] + ci1*w[ijk-jj1+kk1] + ci2*w[ijk+kk1] + ci3*w[ijk+jj1+kk1]) * (ci0*v[ijk-kk1] + ci1*v[ijk    ] + ci2*v[ijk+kk1] + ci3*v[ijk+kk2]))
                 + cg3*((ci0*w[ijk-jj2+kk2] + ci1*w[ijk-jj1+kk2] + ci2*w[ijk+kk2] + ci3*w[ijk+jj1+kk2]) * (ti0*v[ijk-kk1] + ti1*v[ijk    ] + ti2*v[ijk+kk1] + ti3*v[ijk+kk2])) )
                 * dzi4[kend-1];
    }
}

void Advec_4::advecw(double * restrict wt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzhi4)
{
  int    ijk,kstart,kend;
  int    ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

  kstart = grid->kstart;
  kend   = grid->kend;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kstart+1)*kk1;
      wt[ijk] -= ( cg0*((ci0*u[ijk-ii1-kk2] + ci1*u[ijk-ii1-kk1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+kk1]) * (ci0*w[ijk-ii3] + ci1*w[ijk-ii2] + ci2*w[ijk-ii1] + ci3*w[ijk    ]))
                 + cg1*((ci0*u[ijk    -kk2] + ci1*u[ijk    -kk1] + ci2*u[ijk    ] + ci3*u[ijk    +kk1]) * (ci0*w[ijk-ii2] + ci1*w[ijk-ii1] + ci2*w[ijk    ] + ci3*w[ijk+ii1]))
                 + cg2*((ci0*u[ijk+ii1-kk2] + ci1*u[ijk+ii1-kk1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+kk1]) * (ci0*w[ijk-ii1] + ci1*w[ijk    ] + ci2*w[ijk+ii1] + ci3*w[ijk+ii2]))
                 + cg3*((ci0*u[ijk+ii2-kk2] + ci1*u[ijk+ii2-kk1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+ii1] + ci2*w[ijk+ii2] + ci3*w[ijk+ii3])) ) * cgi*dxi;

      wt[ijk] -= ( cg0*((ci0*v[ijk-jj1-kk2] + ci1*v[ijk-jj1-kk1] + ci2*v[ijk-jj1] + ci3*v[ijk-jj1+kk1]) * (ci0*w[ijk-jj3] + ci1*w[ijk-jj2] + ci2*w[ijk-jj1] + ci3*w[ijk    ]))
                 + cg1*((ci0*v[ijk    -kk2] + ci1*v[ijk    -kk1] + ci2*v[ijk    ] + ci3*v[ijk    +kk1]) * (ci0*w[ijk-jj2] + ci1*w[ijk-jj1] + ci2*w[ijk    ] + ci3*w[ijk+jj1]))
                 + cg2*((ci0*v[ijk+jj1-kk2] + ci1*v[ijk+jj1-kk1] + ci2*v[ijk+jj1] + ci3*v[ijk+jj1+kk1]) * (ci0*w[ijk-jj1] + ci1*w[ijk    ] + ci2*w[ijk+jj1] + ci3*w[ijk+jj2]))
                 + cg3*((ci0*v[ijk+jj2-kk2] + ci1*v[ijk+jj2-kk1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+jj1] + ci2*w[ijk+jj2] + ci3*w[ijk+jj3])) ) * cgi*dyi;

      wt[ijk] -= ( cg0*((bi0*w[ijk-kk2] + bi1*w[ijk-kk1] + bi2*w[ijk    ] + bi3*w[ijk+kk1]) * (bi0*w[ijk-kk2] + bi1*w[ijk-kk1] + bi2*w[ijk    ] + bi3*w[ijk+kk1]))
                 + cg1*((ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1]) * (ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1]))
                 + cg2*((ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2]) * (ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2]))
                 + cg3*((ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3]) * (ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3])) )
                 * dzhi4[kstart+1];
    }

  for(int k=grid->kstart+2; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        wt[ijk] -= ( cg0*((ci0*u[ijk-ii1-kk2] + ci1*u[ijk-ii1-kk1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+kk1]) * (ci0*w[ijk-ii3] + ci1*w[ijk-ii2] + ci2*w[ijk-ii1] + ci3*w[ijk    ]))
                   + cg1*((ci0*u[ijk    -kk2] + ci1*u[ijk    -kk1] + ci2*u[ijk    ] + ci3*u[ijk    +kk1]) * (ci0*w[ijk-ii2] + ci1*w[ijk-ii1] + ci2*w[ijk    ] + ci3*w[ijk+ii1]))
                   + cg2*((ci0*u[ijk+ii1-kk2] + ci1*u[ijk+ii1-kk1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+kk1]) * (ci0*w[ijk-ii1] + ci1*w[ijk    ] + ci2*w[ijk+ii1] + ci3*w[ijk+ii2]))
                   + cg3*((ci0*u[ijk+ii2-kk2] + ci1*u[ijk+ii2-kk1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+ii1] + ci2*w[ijk+ii2] + ci3*w[ijk+ii3])) ) * cgi*dxi;

        wt[ijk] -= ( cg0*((ci0*v[ijk-jj1-kk2] + ci1*v[ijk-jj1-kk1] + ci2*v[ijk-jj1] + ci3*v[ijk-jj1+kk1]) * (ci0*w[ijk-jj3] + ci1*w[ijk-jj2] + ci2*w[ijk-jj1] + ci3*w[ijk    ]))
                   + cg1*((ci0*v[ijk    -kk2] + ci1*v[ijk    -kk1] + ci2*v[ijk    ] + ci3*v[ijk    +kk1]) * (ci0*w[ijk-jj2] + ci1*w[ijk-jj1] + ci2*w[ijk    ] + ci3*w[ijk+jj1]))
                   + cg2*((ci0*v[ijk+jj1-kk2] + ci1*v[ijk+jj1-kk1] + ci2*v[ijk+jj1] + ci3*v[ijk+jj1+kk1]) * (ci0*w[ijk-jj1] + ci1*w[ijk    ] + ci2*w[ijk+jj1] + ci3*w[ijk+jj2]))
                   + cg3*((ci0*v[ijk+jj2-kk2] + ci1*v[ijk+jj2-kk1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+jj1] + ci2*w[ijk+jj2] + ci3*w[ijk+jj3])) ) * cgi*dyi;

        wt[ijk] -= ( cg0*((ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ]) * (ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ]))
                   + cg1*((ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1]) * (ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1]))
                   + cg2*((ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2]) * (ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2]))
                   + cg3*((ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3]) * (ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3])) )
                   * dzhi4[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      wt[ijk] -= ( cg0*((ci0*u[ijk-ii1-kk2] + ci1*u[ijk-ii1-kk1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+kk1]) * (ci0*w[ijk-ii3] + ci1*w[ijk-ii2] + ci2*w[ijk-ii1] + ci3*w[ijk    ]))
                 + cg1*((ci0*u[ijk    -kk2] + ci1*u[ijk    -kk1] + ci2*u[ijk    ] + ci3*u[ijk    +kk1]) * (ci0*w[ijk-ii2] + ci1*w[ijk-ii1] + ci2*w[ijk    ] + ci3*w[ijk+ii1]))
                 + cg2*((ci0*u[ijk+ii1-kk2] + ci1*u[ijk+ii1-kk1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+kk1]) * (ci0*w[ijk-ii1] + ci1*w[ijk    ] + ci2*w[ijk+ii1] + ci3*w[ijk+ii2]))
                 + cg3*((ci0*u[ijk+ii2-kk2] + ci1*u[ijk+ii2-kk1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+ii1] + ci2*w[ijk+ii2] + ci3*w[ijk+ii3])) ) * cgi*dxi;

      wt[ijk] -= ( cg0*((ci0*v[ijk-jj1-kk2] + ci1*v[ijk-jj1-kk1] + ci2*v[ijk-jj1] + ci3*v[ijk-jj1+kk1]) * (ci0*w[ijk-jj3] + ci1*w[ijk-jj2] + ci2*w[ijk-jj1] + ci3*w[ijk    ]))
                 + cg1*((ci0*v[ijk    -kk2] + ci1*v[ijk    -kk1] + ci2*v[ijk    ] + ci3*v[ijk    +kk1]) * (ci0*w[ijk-jj2] + ci1*w[ijk-jj1] + ci2*w[ijk    ] + ci3*w[ijk+jj1]))
                 + cg2*((ci0*v[ijk+jj1-kk2] + ci1*v[ijk+jj1-kk1] + ci2*v[ijk+jj1] + ci3*v[ijk+jj1+kk1]) * (ci0*w[ijk-jj1] + ci1*w[ijk    ] + ci2*w[ijk+jj1] + ci3*w[ijk+jj2]))
                 + cg3*((ci0*v[ijk+jj2-kk2] + ci1*v[ijk+jj2-kk1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+jj1] + ci2*w[ijk+jj2] + ci3*w[ijk+jj3])) ) * cgi*dyi;

      wt[ijk] -= ( cg0*((ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ]) * (ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ]))
                 + cg1*((ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1]) * (ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1]))
                 + cg2*((ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2]) * (ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2]))
                 + cg3*((ti0*w[ijk-kk1] + ti1*w[ijk    ] + ti2*w[ijk+kk1] + ti3*w[ijk+kk2]) * (ti0*w[ijk-kk1] + ti1*w[ijk    ] + ti2*w[ijk+kk1] + ti3*w[ijk+kk2])) )
                 * dzhi4[kend-1];
    }
}

void Advec_4::advecs(double * restrict st, double * restrict s, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi4)
{
  int    ijk,kstart,kend;
  int    ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  kstart = grid->kstart;
  kend   = grid->kend;

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + kstart*kk1;
      st[ijk] -= ( cg0*(u[ijk-ii1] * (ci0*s[ijk-ii3] + ci1*s[ijk-ii2] + ci2*s[ijk-ii1] + ci3*s[ijk    ]))
                 + cg1*(u[ijk    ] * (ci0*s[ijk-ii2] + ci1*s[ijk-ii1] + ci2*s[ijk    ] + ci3*s[ijk+ii1]))
                 + cg2*(u[ijk+ii1] * (ci0*s[ijk-ii1] + ci1*s[ijk    ] + ci2*s[ijk+ii1] + ci3*s[ijk+ii2]))
                 + cg3*(u[ijk+ii2] * (ci0*s[ijk    ] + ci1*s[ijk+ii1] + ci2*s[ijk+ii2] + ci3*s[ijk+ii3])) ) * cgi*dxi;

      st[ijk] -= ( cg0*(v[ijk-jj1] * (ci0*s[ijk-jj3] + ci1*s[ijk-jj2] + ci2*s[ijk-jj1] + ci3*s[ijk    ]))
                 + cg1*(v[ijk    ] * (ci0*s[ijk-jj2] + ci1*s[ijk-jj1] + ci2*s[ijk    ] + ci3*s[ijk+jj1]))
                 + cg2*(v[ijk+jj1] * (ci0*s[ijk-jj1] + ci1*s[ijk    ] + ci2*s[ijk+jj1] + ci3*s[ijk+jj2]))
                 + cg3*(v[ijk+jj2] * (ci0*s[ijk    ] + ci1*s[ijk+jj1] + ci2*s[ijk+jj2] + ci3*s[ijk+jj3])) ) * cgi*dyi;

      st[ijk] -= ( cg0*(w[ijk-kk1] * (bi0*s[ijk-kk2] + bi1*s[ijk-kk1] + bi2*s[ijk    ] + bi3*s[ijk+kk1]))
                 + cg1*(w[ijk    ] * (ci0*s[ijk-kk2] + ci1*s[ijk-kk1] + ci2*s[ijk    ] + ci3*s[ijk+kk1]))
                 + cg2*(w[ijk+kk1] * (ci0*s[ijk-kk1] + ci1*s[ijk    ] + ci2*s[ijk+kk1] + ci3*s[ijk+kk2]))
                 + cg3*(w[ijk+kk2] * (ci0*s[ijk    ] + ci1*s[ijk+kk1] + ci2*s[ijk+kk2] + ci3*s[ijk+kk3])) )
                 * dzi4[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        st[ijk] -= ( cg0*(u[ijk-ii1] * (ci0*s[ijk-ii3] + ci1*s[ijk-ii2] + ci2*s[ijk-ii1] + ci3*s[ijk    ]))
                   + cg1*(u[ijk    ] * (ci0*s[ijk-ii2] + ci1*s[ijk-ii1] + ci2*s[ijk    ] + ci3*s[ijk+ii1]))
                   + cg2*(u[ijk+ii1] * (ci0*s[ijk-ii1] + ci1*s[ijk    ] + ci2*s[ijk+ii1] + ci3*s[ijk+ii2]))
                   + cg3*(u[ijk+ii2] * (ci0*s[ijk    ] + ci1*s[ijk+ii1] + ci2*s[ijk+ii2] + ci3*s[ijk+ii3])) ) * cgi*dxi;

        st[ijk] -= ( cg0*(v[ijk-jj1] * (ci0*s[ijk-jj3] + ci1*s[ijk-jj2] + ci2*s[ijk-jj1] + ci3*s[ijk    ]))
                   + cg1*(v[ijk    ] * (ci0*s[ijk-jj2] + ci1*s[ijk-jj1] + ci2*s[ijk    ] + ci3*s[ijk+jj1]))
                   + cg2*(v[ijk+jj1] * (ci0*s[ijk-jj1] + ci1*s[ijk    ] + ci2*s[ijk+jj1] + ci3*s[ijk+jj2]))
                   + cg3*(v[ijk+jj2] * (ci0*s[ijk    ] + ci1*s[ijk+jj1] + ci2*s[ijk+jj2] + ci3*s[ijk+jj3])) ) * cgi*dyi;

        st[ijk] -= ( cg0*(w[ijk-kk1] * (ci0*s[ijk-kk3] + ci1*s[ijk-kk2] + ci2*s[ijk-kk1] + ci3*s[ijk    ]))
                   + cg1*(w[ijk    ] * (ci0*s[ijk-kk2] + ci1*s[ijk-kk1] + ci2*s[ijk    ] + ci3*s[ijk+kk1]))
                   + cg2*(w[ijk+kk1] * (ci0*s[ijk-kk1] + ci1*s[ijk    ] + ci2*s[ijk+kk1] + ci3*s[ijk+kk2]))
                   + cg3*(w[ijk+kk2] * (ci0*s[ijk    ] + ci1*s[ijk+kk1] + ci2*s[ijk+kk2] + ci3*s[ijk+kk3])) )
                   * dzi4[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      st[ijk] -= ( cg0*(u[ijk-ii1] * (ci0*s[ijk-ii3] + ci1*s[ijk-ii2] + ci2*s[ijk-ii1] + ci3*s[ijk    ]))
                 + cg1*(u[ijk    ] * (ci0*s[ijk-ii2] + ci1*s[ijk-ii1] + ci2*s[ijk    ] + ci3*s[ijk+ii1]))
                 + cg2*(u[ijk+ii1] * (ci0*s[ijk-ii1] + ci1*s[ijk    ] + ci2*s[ijk+ii1] + ci3*s[ijk+ii2]))
                 + cg3*(u[ijk+ii2] * (ci0*s[ijk    ] + ci1*s[ijk+ii1] + ci2*s[ijk+ii2] + ci3*s[ijk+ii3])) ) * cgi*dxi;

      st[ijk] -= ( cg0*(v[ijk-jj1] * (ci0*s[ijk-jj3] + ci1*s[ijk-jj2] + ci2*s[ijk-jj1] + ci3*s[ijk    ]))
                 + cg1*(v[ijk    ] * (ci0*s[ijk-jj2] + ci1*s[ijk-jj1] + ci2*s[ijk    ] + ci3*s[ijk+jj1]))
                 + cg2*(v[ijk+jj1] * (ci0*s[ijk-jj1] + ci1*s[ijk    ] + ci2*s[ijk+jj1] + ci3*s[ijk+jj2]))
                 + cg3*(v[ijk+jj2] * (ci0*s[ijk    ] + ci1*s[ijk+jj1] + ci2*s[ijk+jj2] + ci3*s[ijk+jj3])) ) * cgi*dyi;

      st[ijk] -= ( cg0*(w[ijk-kk1] * (ci0*s[ijk-kk3] + ci1*s[ijk-kk2] + ci2*s[ijk-kk1] + ci3*s[ijk    ]))
                 + cg1*(w[ijk    ] * (ci0*s[ijk-kk2] + ci1*s[ijk-kk1] + ci2*s[ijk    ] + ci3*s[ijk+kk1]))
                 + cg2*(w[ijk+kk1] * (ci0*s[ijk-kk1] + ci1*s[ijk    ] + ci2*s[ijk+kk1] + ci3*s[ijk+kk2]))
                 + cg3*(w[ijk+kk2] * (ti0*s[ijk-kk1] + ti1*s[ijk    ] + ti2*s[ijk+kk1] + ti3*s[ijk+kk2])) )
                 * dzi4[kend-1];
    }
}

