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

#define NO_OFFSET 0.
#define NO_VELOCITY 0.

#define BC_DIRICHLET 0
#define BC_NEUMANN 1
#define BC_FLUX 2
#define BC_USTAR 3

__global__ void boundary_setgcbot_2nd(double * __restrict__ a, double * __restrict__ dzh, int sw, 
                                      double * __restrict__ abot, double * __restrict__ agradbot,
                                      const unsigned int icells, const unsigned int icellsp, const unsigned int jcells, const unsigned int kstart)
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
                                      const unsigned int icells, const unsigned int icellsp, const unsigned int jcells, const unsigned int kend)
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

  // cyclic boundary conditions, do this before the bottom BC's
  //grid->boundary_cyclic(fields->u->data_g);
  //grid->boundary_cyclic(fields->v->data_g);
  //grid->boundary_cyclic(fields->w->data_g);

  const int offs = grid->memoffset;

  grid->boundary_cyclic_gpu(&fields->u->data_g[offs]);
  grid->boundary_cyclic_gpu(&fields->v->data_g[offs]);
  grid->boundary_cyclic_gpu(&fields->w->data_g[offs]);

  for(fieldmap::const_iterator it = fields->sp.begin(); it!=fields->sp.end(); ++it)
    grid->boundary_cyclic_gpu(&it->second->data_g[offs]);

  // calculate boundary values
  bcvalues();

  if(grid->swspatialorder == "2")
  {
    boundary_setgcbot_2nd<<<grid2dGPU, block2dGPU>>>(&fields->u->data_g[offs], grid->dzh_g, mbcbot, 
                                                     &fields->u->databot_g[offs], &fields->u->datagradbot_g[offs],
                                                     grid->icells, grid->icellsp, grid->jcells, grid->kstart);
    boundary_setgctop_2nd<<<grid2dGPU, block2dGPU>>>(&fields->u->data_g[offs], grid->dzh_g, mbctop, 
                                                     &fields->u->datatop_g[offs], &fields->u->datagradtop_g[offs],
                                                     grid->icells, grid->icellsp, grid->jcells, grid->kend);

    boundary_setgcbot_2nd<<<grid2dGPU, block2dGPU>>>(&fields->v->data_g[offs], grid->dzh_g, mbcbot, 
                                                     &fields->v->databot_g[offs], &fields->v->datagradbot_g[offs],
                                                     grid->icells, grid->icellsp, grid->jcells, grid->kstart);
    boundary_setgctop_2nd<<<grid2dGPU, block2dGPU>>>(&fields->v->data_g[offs], grid->dzh_g, mbctop, 
                                                     &fields->v->datatop_g[offs], &fields->v->datagradtop_g[offs],
                                                     grid->icells, grid->icellsp, grid->jcells, grid->kend);

    for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
      boundary_setgcbot_2nd<<<grid2dGPU, block2dGPU>>>(&it->second->data_g[offs], grid->dzh_g, sbc[it->first]->bcbot, 
                                                       &it->second->databot_g[offs], &it->second->datagradbot_g[offs],
                                                       grid->icells, grid->icellsp, grid->jcells, grid->kstart);
      boundary_setgctop_2nd<<<grid2dGPU, block2dGPU>>>(&it->second->data_g[offs], grid->dzh_g, sbc[it->first]->bctop, 
                                                       &it->second->datatop_g[offs], &it->second->datagradtop_g[offs],
                                                       grid->icells, grid->icellsp, grid->jcells, grid->kend);
    }
  }
  //else if(grid->swspatialorder == "4")
  //{
  //  setgcbot_4th(fields->u->data, grid->z, mbcbot, fields->u->databot, fields->u->datagradbot);
  //  setgctop_4th(fields->u->data, grid->z, mbctop, fields->u->datatop, fields->u->datagradtop);

  //  setgcbot_4th(fields->v->data, grid->z, mbcbot, fields->v->databot, fields->v->datagradbot);
  //  setgctop_4th(fields->v->data, grid->z, mbctop, fields->v->datatop, fields->v->datagradtop);

  //  setgcbotw_4th(fields->w->data);
  //  setgctopw_4th(fields->w->data);

  //  for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
  //  {
  //    setgcbot_4th(it->second->data, grid->z, sbc[it->first]->bcbot, it->second->databot, it->second->datagradbot);
  //    setgctop_4th(it->second->data, grid->z, sbc[it->first]->bctop, it->second->datatop, it->second->datagradtop);
  //  }
  //}


  fields->backwardGPU();

  return 0;
}
#endif





