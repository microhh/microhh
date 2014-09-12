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
#include "master.h"
#include "diff_4.h"
#include "defines.h"
#include "model.h"

cdiff_4::cdiff_4(cmodel *modelin, cinput *inputin) : cdiff(modelin, inputin)
{
  swdiff = "4";
}

cdiff_4::~cdiff_4()
{
}

int cdiff_4::setvalues()
{
  // get the maximum time step for diffusion
  double viscmax = fields->visc;
  for(fieldmap::iterator it = fields->sp.begin(); it!=fields->sp.end(); it++)
    viscmax = std::max(it->second->visc, viscmax);

  dnmul = 0;
  for(int k=grid->kstart; k<grid->kend; k++)
    dnmul = std::max(dnmul, std::abs(viscmax * (1./(grid->dx*grid->dx) + 1./(grid->dy*grid->dy) + 1./(grid->dz[k]*grid->dz[k]))));

  return 0;
}

unsigned long cdiff_4::gettimelim(unsigned long idt, double dt)
{
  unsigned long idtlim;

  idtlim = idt * dnmax / (dt * dnmul);

  return idtlim;
}

double cdiff_4::getdn(double dt)
{
  double dn;

  dn = dnmul*dt;

  return dn;
}

#ifndef USECUDA
int cdiff_4::exec()
{
  diffc(fields->ut->data, fields->u->data, grid->dzi4, grid->dzhi4, fields->visc);
  diffc(fields->vt->data, fields->v->data, grid->dzi4, grid->dzhi4, fields->visc);
  diffw(fields->wt->data, fields->w->data, grid->dzi4, grid->dzhi4, fields->visc);

  for(fieldmap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    diffc(it->second->data, fields->s[it->first]->data, grid->dzi4, grid->dzhi4, fields->s[it->first]->visc);

  return 0;
}
#endif

int cdiff_4::diffc(double * restrict at, double * restrict a, double * restrict dzi4, double * restrict dzhi4, double visc)
{
  int    ijk,kstart,kend;
  int    ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  double dxidxi,dyidyi;

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

  dxidxi = 1./(grid->dx * grid->dx);
  dyidyi = 1./(grid->dy * grid->dy);

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + kstart*kk1;
      at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
      at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
      at[ijk] += visc * ( cg0*(bg0*a[ijk-kk2] + bg1*a[ijk-kk1] + bg2*a[ijk    ] + bg3*a[ijk+kk1]) * dzhi4[kstart-1]
                        + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzhi4[kstart  ]
                        + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[kstart+1]
                        + cg3*(cg0*a[ijk    ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzhi4[kstart+2] )
                        * dzi4[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
        at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
        at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*a[ijk    ]) * dzhi4[k-1]
                          + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzhi4[k  ]
                          + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[k+1]
                          + cg3*(cg0*a[ijk    ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzhi4[k+2] )
                          * dzi4[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
      at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
      at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*a[ijk    ]) * dzhi4[kend-2]
                        + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzhi4[kend-1]
                        + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[kend  ]
                        + cg3*(tg0*a[ijk-kk1] + tg1*a[ijk    ] + tg2*a[ijk+kk1] + tg3*a[ijk+kk2]) * dzhi4[kend+1] )
                        * dzi4[kend-1];
    }

  return 0;
}

int cdiff_4::diffw(double * restrict at, double * restrict a, double * restrict dzi4, double * restrict dzhi4, double visc)
{
  int    ijk,kstart,kend;
  int    ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  double dxidxi,dyidyi;

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

  dxidxi = 1./(grid->dx * grid->dx);
  dyidyi = 1./(grid->dy * grid->dy);

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kstart+1)*kk1;
      at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
      at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
      at[ijk] += visc * ( cg0*(bg0*a[ijk-kk2] + bg1*a[ijk-kk1] + bg2*a[ijk    ] + bg3*a[ijk+kk1]) * dzi4[kstart-1]
                        + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzi4[kstart  ]
                        + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzi4[kstart+1]
                        + cg3*(cg0*a[ijk    ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzi4[kstart+2] )
                        * dzhi4[kstart+1];
    }

  for(int k=grid->kstart+2; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
        at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
        at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*a[ijk    ]) * dzi4[k-2]
                          + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzi4[k-1]
                          + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzi4[k  ]
                          + cg3*(cg0*a[ijk    ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzi4[k+1] )
                          * dzhi4[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
      at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
      at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*a[ijk    ]) * dzi4[kend-3]
                        + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzi4[kend-2]
                        + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzi4[kend-1]
                        + cg3*(tg0*a[ijk-kk1] + tg1*a[ijk    ] + tg2*a[ijk+kk1] + tg3*a[ijk+kk2]) * dzi4[kend  ] )
                        * dzhi4[kend-1];
    }

  return 0;
}
