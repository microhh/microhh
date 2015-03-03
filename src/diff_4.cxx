/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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
#include "fd.h"
#include "model.h"

using namespace fd::o4;

Diff_4::Diff_4(Model *modelin, Input *inputin) : Diff(modelin, inputin)
{
  swdiff = "4";
}

Diff_4::~Diff_4()
{
}

void Diff_4::set_values()
{
  // get the maximum time step for diffusion
  double viscmax = fields->visc;
  for (FieldMap::iterator it = fields->sp.begin(); it!=fields->sp.end(); it++)
    viscmax = std::max(it->second->visc, viscmax);

  dnmul = 0;
  for (int k=grid->kstart; k<grid->kend; k++)
    dnmul = std::max(dnmul, std::abs(viscmax * (1./(grid->dx*grid->dx) + 1./(grid->dy*grid->dy) + 1./(grid->dz[k]*grid->dz[k]))));
}

unsigned long Diff_4::get_time_limit(unsigned long idt, double dt)
{
  unsigned long idtlim;

  idtlim = idt * dnmax / (dt * dnmul);

  return idtlim;
}

double Diff_4::get_dn(double dt)
{
  double dn;

  dn = dnmul*dt;

  return dn;
}

#ifndef USECUDA
void Diff_4::exec()
{
  // In case of a two-dimensional run, strip v component out of all kernels and do 
  // not calculate v-diffusion tendency.
  if (grid->jtot == 1)
  {
    diff_c<false>(fields->ut->data, fields->u->data, grid->dzi4, grid->dzhi4, fields->visc);
    diff_w<false>(fields->wt->data, fields->w->data, grid->dzi4, grid->dzhi4, fields->visc);

    for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
      diff_c<false>(it->second->data, fields->sp[it->first]->data, grid->dzi4, grid->dzhi4, fields->sp[it->first]->visc);
  }
  else
  {
    diff_c<true>(fields->ut->data, fields->u->data, grid->dzi4, grid->dzhi4, fields->visc);
    diff_c<true>(fields->vt->data, fields->v->data, grid->dzi4, grid->dzhi4, fields->visc);
    diff_w<true>(fields->wt->data, fields->w->data, grid->dzi4, grid->dzhi4, fields->visc);

    for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
      diff_c<true>(it->second->data, fields->sp[it->first]->data, grid->dzi4, grid->dzhi4, fields->sp[it->first]->visc);
  }
}
#endif

template<bool dim3>
void Diff_4::diff_c(double* restrict at, double* restrict a, double* restrict dzi4, double* restrict dzhi4, const double visc)
{
  const int ii1 = 1;
  const int ii2 = 2;
  const int ii3 = 3;
  const int jj1 = 1*grid->icells;
  const int jj2 = 2*grid->icells;
  const int jj3 = 3*grid->icells;
  const int kk1 = 1*grid->ijcells;
  const int kk2 = 2*grid->ijcells;
  const int kk3 = 3*grid->ijcells;

  const int kstart = grid->kstart;
  const int kend   = grid->kend;

  const double dxidxi = 1./(grid->dx * grid->dx);
  const double dyidyi = 1./(grid->dy * grid->dy);

  int ijk;

  // bottom boundary
  for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for (int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + kstart*kk1;
      at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
      if (dim3)
        at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
      at[ijk] += visc * ( cg0*(bg0*a[ijk-kk2] + bg1*a[ijk-kk1] + bg2*a[ijk    ] + bg3*a[ijk+kk1]) * dzhi4[kstart-1]
                        + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzhi4[kstart  ]
                        + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[kstart+1]
                        + cg3*(cg0*a[ijk    ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzhi4[kstart+2] )
                        * dzi4[kstart];
    }

  for (int k=grid->kstart+1; k<grid->kend-1; k++)
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for (int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
        if (dim3)
          at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
        at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*a[ijk    ]) * dzhi4[k-1]
                          + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzhi4[k  ]
                          + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[k+1]
                          + cg3*(cg0*a[ijk    ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzhi4[k+2] )
                          * dzi4[k];
      }

  // top boundary
  for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for (int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
      if (dim3)
        at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
      at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*a[ijk    ]) * dzhi4[kend-2]
                        + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzhi4[kend-1]
                        + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzhi4[kend  ]
                        + cg3*(tg0*a[ijk-kk1] + tg1*a[ijk    ] + tg2*a[ijk+kk1] + tg3*a[ijk+kk2]) * dzhi4[kend+1] )
                        * dzi4[kend-1];
    }
}

template<bool dim3>
void Diff_4::diff_w(double* restrict at, double* restrict a, double* restrict dzi4, double* restrict dzhi4, double visc)
{
  const int ii1 = 1;
  const int ii2 = 2;
  const int ii3 = 3;
  const int jj1 = 1*grid->icells;
  const int jj2 = 2*grid->icells;
  const int jj3 = 3*grid->icells;
  const int kk1 = 1*grid->ijcells;
  const int kk2 = 2*grid->ijcells;
  const int kk3 = 3*grid->ijcells;

  const int kstart = grid->kstart;
  const int kend   = grid->kend;

  const double dxidxi = 1./(grid->dx * grid->dx);
  const double dyidyi = 1./(grid->dy * grid->dy);

  int ijk;

  // bottom boundary
  for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for (int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kstart+1)*kk1;
      at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
      if (dim3)
        at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
      at[ijk] += visc * ( cg0*(bg0*a[ijk-kk2] + bg1*a[ijk-kk1] + bg2*a[ijk    ] + bg3*a[ijk+kk1]) * dzi4[kstart-1]
                        + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzi4[kstart  ]
                        + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzi4[kstart+1]
                        + cg3*(cg0*a[ijk    ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzi4[kstart+2] )
                        * dzhi4[kstart+1];
    }

  for (int k=grid->kstart+2; k<grid->kend-1; k++)
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for (int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
        if (dim3)
          at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
        at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*a[ijk    ]) * dzi4[k-2]
                          + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzi4[k-1]
                          + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzi4[k  ]
                          + cg3*(cg0*a[ijk    ] + cg1*a[ijk+kk1] + cg2*a[ijk+kk2] + cg3*a[ijk+kk3]) * dzi4[k+1] )
                          * dzhi4[k];
      }

  // top boundary
  for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for (int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      at[ijk] += visc * (cdg3*a[ijk-ii3] + cdg2*a[ijk-ii2] + cdg1*a[ijk-ii1] + cdg0*a[ijk] + cdg1*a[ijk+ii1] + cdg2*a[ijk+ii2] + cdg3*a[ijk+ii3])*dxidxi;
      if (dim3)
        at[ijk] += visc * (cdg3*a[ijk-jj3] + cdg2*a[ijk-jj2] + cdg1*a[ijk-jj1] + cdg0*a[ijk] + cdg1*a[ijk+jj1] + cdg2*a[ijk+jj2] + cdg3*a[ijk+jj3])*dyidyi;
      at[ijk] += visc * ( cg0*(cg0*a[ijk-kk3] + cg1*a[ijk-kk2] + cg2*a[ijk-kk1] + cg3*a[ijk    ]) * dzi4[kend-3]
                        + cg1*(cg0*a[ijk-kk2] + cg1*a[ijk-kk1] + cg2*a[ijk    ] + cg3*a[ijk+kk1]) * dzi4[kend-2]
                        + cg2*(cg0*a[ijk-kk1] + cg1*a[ijk    ] + cg2*a[ijk+kk1] + cg3*a[ijk+kk2]) * dzi4[kend-1]
                        + cg3*(tg0*a[ijk-kk1] + tg1*a[ijk    ] + tg2*a[ijk+kk1] + tg3*a[ijk+kk2]) * dzi4[kend  ] )
                        * dzhi4[kend-1];
    }
}
