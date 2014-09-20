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
#include "diff_2.h"
#include "defines.h"
#include "model.h"

cdiff_2::cdiff_2(cmodel *modelin, cinput *inputin) : cdiff(modelin, inputin)
{
  swdiff = "2";
}

cdiff_2::~cdiff_2()
{
}

void cdiff_2::setvalues()
{
  // get the maximum time step for diffusion
  double viscmax = fields->visc;
  for(fieldmap::iterator it = fields->sp.begin(); it!=fields->sp.end(); it++)
    viscmax = std::max(it->second->visc, viscmax);

  dnmul = 0;
  for(int k=grid->kstart; k<grid->kend; k++)
    dnmul = std::max(dnmul, std::abs(viscmax * (1./(grid->dx*grid->dx) + 1./(grid->dy*grid->dy) + 1./(grid->dz[k]*grid->dz[k]))));
}

unsigned long cdiff_2::gettimelim(unsigned long idt, double dt)
{
  unsigned long idtlim;

  idtlim = idt * dnmax / (dt * dnmul);

  return idtlim;
}

double cdiff_2::getdn(double dt)
{
  double dn;

  dn = dnmul*dt;

  return dn;
}

#ifndef USECUDA
int cdiff_2::exec()
{
  diffc(fields->ut->data, fields->u->data, grid->dzi, grid->dzhi, fields->visc);
  diffc(fields->vt->data, fields->v->data, grid->dzi, grid->dzhi, fields->visc);
  diffw(fields->wt->data, fields->w->data, grid->dzi, grid->dzhi, fields->visc);

  for(fieldmap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    diffc(it->second->data, fields->s[it->first]->data, grid->dzi, grid->dzhi, fields->s[it->first]->visc);
  
  return 0;
}
#endif

int cdiff_2::diffc(double * restrict at, double * restrict a, double * restrict dzi, double * restrict dzhi, double visc)
{
  int    ijk,ii,jj,kk;
  double dxidxi,dyidyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dxidxi = 1./(grid->dx * grid->dx);
  dyidyi = 1./(grid->dy * grid->dy);

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        at[ijk] += visc * (
              + (  (a[ijk+ii] - a[ijk   ]) 
                 - (a[ijk   ] - a[ijk-ii]) ) * dxidxi 
              + (  (a[ijk+jj] - a[ijk   ]) 
                 - (a[ijk   ] - a[ijk-jj]) ) * dyidyi
              + (  (a[ijk+kk] - a[ijk   ]) * dzhi[k+1]
                 - (a[ijk   ] - a[ijk-kk]) * dzhi[k]   ) * dzi[k] );
      }

  return 0;
}

int cdiff_2::diffw(double * restrict wt, double * restrict w, double * restrict dzi, double * restrict dzhi, double visc)
{
  int    ijk,ii,jj,kk;
  double dxidxi,dyidyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dxidxi = 1./(grid->dx*grid->dx);
  dyidyi = 1./(grid->dy*grid->dy);

  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        wt[ijk] += visc * (
              + (  (w[ijk+ii] - w[ijk   ]) 
                 - (w[ijk   ] - w[ijk-ii]) ) * dxidxi 
              + (  (w[ijk+jj] - w[ijk   ]) 
                 - (w[ijk   ] - w[ijk-jj]) ) * dyidyi
              + (  (w[ijk+kk] - w[ijk   ]) * dzi[k]
                 - (w[ijk   ] - w[ijk-kk]) * dzi[k-1] ) * dzhi[k] );
      }

  return 0;
}

