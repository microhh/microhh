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
#include "boundary.h"
#include "defines.h"
#include "model.h"
#include "timeloop.h"
#include "fd.h"

#define NO_OFFSET 0.
#define NO_VELOCITY 0.

#define BC_DIRICHLET 0
#define BC_NEUMANN 1
#define BC_FLUX 2
#define BC_USTAR 3

using namespace fd::o4;

__global__ void boundary_setgcbot_2nd(double * __restrict__ a, double * __restrict__ dzh, int sw, 
                                      double * __restrict__ abot, double * __restrict__ agradbot,
                                      const int icells, const int icellsp,
                                      const int jcells, const int kstart)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;

  int kk  = icellsp*jcells;
  int ij  = i + j*icellsp;
  int ijk = i + j*icellsp + kstart*kk;

  if(i < icells && j < jcells)
  {
    if(sw == BC_DIRICHLET)
      a[ijk-kk] = 2.*abot[ij] - a[ijk];

    else if(sw == BC_NEUMANN || sw == BC_FLUX)
      a[ijk-kk] = -agradbot[ij]*dzh[kstart] + a[ijk];
  }
} 

__global__ void boundary_setgctop_2nd(double * __restrict__ a, double * __restrict__ dzh, int sw,
                                      double * __restrict__ atop, double * __restrict__ agradtop,
                                      const int icells, const int icellsp,
                                      const int jcells, const int kend)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x;
  int j = blockIdx.y*blockDim.y + threadIdx.y;

  int kk  = icellsp*jcells;
  int ij  = i + j*icellsp;
  int ijk = i + j*icellsp + (kend-1)*kk;

  if(i < icells && j < jcells)
  {
    if(sw == BC_DIRICHLET)
      a[ijk+kk] = 2.*atop[ij] - a[ijk];

    else if(sw == BC_NEUMANN || sw == BC_FLUX)
      a[ijk+kk] = agradtop[ij]*dzh[kend] + a[ijk];
  }
}

__device__ double fd_grad4x(const double a, const double b, const double c, const double d)
{
  return (-(d-a) + 27.*(c-b));
}

__global__ void boundary_setgcbot_4th(double * __restrict__ a, const int sw,
                                      double * __restrict__ abot, double * __restrict__ agradbot,
                                      double * __restrict z,
                                      const int icells, const int icellsp,
                                      const int jcells, const int kstart)
{
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  const int j = blockIdx.y*blockDim.y + threadIdx.y;

  const int kk1 = 1*icellsp*jcells;
  const int kk2 = 2*icellsp*jcells;

  const int ij  = i + j*icellsp;
  const int ijk = i + j*icellsp + kstart*kk1;

  if(i < icells && j < jcells)
  {
    if(sw == BC_DIRICHLET)
    {
      a[ijk-kk1] = (8./3.)*abot[ij] - 2.*a[ijk] + (1./3.)*a[ijk+kk1];
      a[ijk-kk2] = 8.*abot[ij] - 9.*a[ijk] + 2.*a[ijk+kk1];
    }

    else if(sw == BC_NEUMANN || sw == BC_FLUX)
    {
      a[ijk-kk1] = -(1./24.)*fd_grad4x(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk    ];
      a[ijk-kk2] = -(1./ 8.)*fd_grad4x(z[kstart-2], z[kstart-1], z[kstart], z[kstart+1])*agradbot[ij] + a[ijk+kk1];
    }
  }
} 

__global__ void boundary_setgctop_4th(double * __restrict__ a, const int sw,
                                      double * __restrict__ atop, double * __restrict__ agradtop,
                                      double * __restrict__ z,
                                      const int icells, const int icellsp,
                                      const int jcells, const int kend)
{
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  const int j = blockIdx.y*blockDim.y + threadIdx.y;

  const int kk1 = 1*icellsp*jcells;
  const int kk2 = 2*icellsp*jcells;

  const int ij  = i + j*icellsp;
  const int ijk = i + j*icellsp + (kend-1)*kk1;

  if(i < icells && j < jcells)
  {
    if(sw == BC_DIRICHLET)
    {
      a[ijk+kk1] = (8./3.)*atop[ij] - 2.*a[ijk] + (1./3.)*a[ijk-kk1];
      a[ijk+kk2] = 8.*atop[ij] - 9.*a[ijk] + 2.*a[ijk-kk1];
    }

    else if(sw == BC_NEUMANN || sw == BC_FLUX)
    {
      a[ijk+kk1] = (1./24.)*fd_grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk    ];
      a[ijk+kk2] = (1./ 8.)*fd_grad4x(z[kend-2], z[kend-1], z[kend], z[kend+1])*agradtop[ij] + a[ijk-kk1];
    }
  }
} 

__global__ void boundary_setgcbotw_4th(double * __restrict__ w,
                                       const int icells, const int icellsp,
                                       const int jcells, const int kstart)
{
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  const int j = blockIdx.y*blockDim.y + threadIdx.y;

  const int kk1 = 1*icellsp*jcells;
  const int kk2 = 2*icellsp*jcells;

  const int ijk = i + j*icellsp + kstart*kk1;

  if(i < icells && j < jcells)
  {
    w[ijk-kk1] = -w[ijk+kk1];
    w[ijk-kk2] = -w[ijk+kk2];
  }
}

__global__ void boundary_setgctopw_4th(double * __restrict__ w,
                                       const int icells, const int icellsp,
                                       const int jcells, const int kend)
{
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  const int j = blockIdx.y*blockDim.y + threadIdx.y;

  const int kk1 = 1*icellsp*jcells;
  const int kk2 = 2*icellsp*jcells;

  const int ijk = i + j*icellsp + (kend-1)*kk1;

  if(i < icells && j < jcells)
  {
    w[ijk+kk1] = -w[ijk-kk1];
    w[ijk+kk2] = -w[ijk-kk2];
  }
}

#ifdef USECUDA
int cboundary::exec()
{
  fields->forwardGPU();

  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->icells/blocki + (grid->icells%blocki > 0);
  const int gridj  = grid->jcells/blockj + (grid->jcells%blockj > 0);

  dim3 grid2dGPU (gridi, gridj);
  dim3 block2dGPU(blocki, blockj);

  const int offs = grid->memoffset;

  // Cyclic boundary conditions, do this before the bottom BC's
  grid->boundary_cyclic_g(&fields->u->data_g[offs]);
  grid->boundary_cyclic_g(&fields->v->data_g[offs]);
  grid->boundary_cyclic_g(&fields->w->data_g[offs]);

  for(fieldmap::const_iterator it = fields->sp.begin(); it!=fields->sp.end(); ++it)
    grid->boundary_cyclic_g(&it->second->data_g[offs]);

  // Calculate the boundary values
  bcvalues();

  fields->backwardGPU();

  if(grid->swspatialorder == "2")
  {
    boundary_setgcbot_2nd<<<grid2dGPU, block2dGPU>>>(&fields->u->data_g[offs], grid->dzh_g, mbcbot, 
                                                     &fields->u->databot_g[offs], &fields->u->datagradbot_g[offs],
                                                     grid->icells, grid->icellsp,
                                                     grid->jcells, grid->kstart);
    boundary_setgctop_2nd<<<grid2dGPU, block2dGPU>>>(&fields->u->data_g[offs], grid->dzh_g, mbctop, 
                                                     &fields->u->datatop_g[offs], &fields->u->datagradtop_g[offs],
                                                     grid->icells, grid->icellsp,
                                                     grid->jcells, grid->kend);

    boundary_setgcbot_2nd<<<grid2dGPU, block2dGPU>>>(&fields->v->data_g[offs], grid->dzh_g, mbcbot, 
                                                     &fields->v->databot_g[offs], &fields->v->datagradbot_g[offs],
                                                     grid->icells, grid->icellsp,
                                                     grid->jcells, grid->kstart);
    boundary_setgctop_2nd<<<grid2dGPU, block2dGPU>>>(&fields->v->data_g[offs], grid->dzh_g, mbctop, 
                                                     &fields->v->datatop_g[offs], &fields->v->datagradtop_g[offs],
                                                     grid->icells, grid->icellsp,
                                                     grid->jcells, grid->kend);

    for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
      boundary_setgcbot_2nd<<<grid2dGPU, block2dGPU>>>(&it->second->data_g[offs], grid->dzh_g, sbc[it->first]->bcbot, 
                                                       &it->second->databot_g[offs], &it->second->datagradbot_g[offs],
                                                       grid->icells, grid->icellsp,
                                                       grid->jcells, grid->kstart);
      boundary_setgctop_2nd<<<grid2dGPU, block2dGPU>>>(&it->second->data_g[offs], grid->dzh_g, sbc[it->first]->bctop, 
                                                       &it->second->datatop_g[offs], &it->second->datagradtop_g[offs],
                                                       grid->icells, grid->icellsp,
                                                       grid->jcells, grid->kend);
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
    /*
    boundary_setgcbot_4th<<<grid2dGPU, block2dGPU>>>(&fields->u->data_g[offs], mbcbot, 
                                                     &fields->u->databot_g[offs], &fields->u->datagradbot_g[offs],
                                                     grid->z_g,
                                                     grid->icells, grid->icellsp,
                                                     grid->jcells, grid->kstart);
    boundary_setgctop_4th<<<grid2dGPU, block2dGPU>>>(&fields->u->data_g[offs], mbctop, 
                                                     &fields->u->datatop_g[offs], &fields->u->datagradtop_g[offs],
                                                     grid->z_g,
                                                     grid->icells, grid->icellsp,
                                                     grid->jcells, grid->kend);

    boundary_setgcbot_4th<<<grid2dGPU, block2dGPU>>>(&fields->v->data_g[offs], mbcbot, 
                                                     &fields->v->databot_g[offs], &fields->v->datagradbot_g[offs],
                                                     grid->z_g,
                                                     grid->icells, grid->icellsp,
                                                     grid->jcells, grid->kstart);
    boundary_setgctop_4th<<<grid2dGPU, block2dGPU>>>(&fields->v->data_g[offs], mbctop, 
                                                     &fields->v->datatop_g[offs], &fields->v->datagradtop_g[offs],
                                                     grid->z_g,
                                                     grid->icells, grid->icellsp,
                                                     grid->jcells, grid->kend);

    boundary_setgcbotw_4th<<<grid2dGPU, block2dGPU>>>(&fields->w->data_g[offs],
                                                      grid->icells, grid->icellsp,
                                                      grid->jcells, grid->kstart);
    boundary_setgctopw_4th<<<grid2dGPU, block2dGPU>>>(&fields->w->data_g[offs],
                                                      grid->icells, grid->icellsp,
                                                      grid->jcells, grid->kend);

    for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
      boundary_setgcbot_4th<<<grid2dGPU, block2dGPU>>>(&it->second->data_g[offs], sbc[it->first]->bcbot,
                                                       &it->second->databot_g[offs], &it->second->datagradbot_g[offs],
                                                       grid->z_g,
                                                       grid->icells, grid->icellsp,
                                                       grid->jcells, grid->kstart);
      boundary_setgctop_4th<<<grid2dGPU, block2dGPU>>>(&it->second->data_g[offs], sbc[it->first]->bctop, 
                                                       &it->second->datatop_g[offs], &it->second->datagradtop_g[offs],
                                                       grid->z_g,
                                                       grid->icells, grid->icellsp,
                                                       grid->jcells, grid->kend);
    }*/
  }

  return 0;
}
#endif
