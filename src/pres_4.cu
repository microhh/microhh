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

#include "grid.h"
#include "fields.h"
#include "pres_4.h"
#include "fd.h"

using namespace fd::o4;

__global__ void pres_4_calcdivergence(double * __restrict__ div,
                                      double * __restrict__ u, double * __restrict__ v, double * __restrict__ w,
                                      double * __restrict__ dzi4,
                                      double dxi, double dyi,
                                      int jj, int kk,
                                      int istart, int jstart, int kstart,
                                      int iend, int jend, int kend)
{
  const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
  const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
  const int k = blockIdx.z + kstart;

  const int ii1 = 1;
  const int ii2 = 2;
  const int jj1 = 1*jj;
  const int jj2 = 2*jj;
  const int kk1 = 1*kk;
  const int kk2 = 2*kk;

  if(i < iend && j < jend && k < kend)
  {
    const int ijk = i + j*jj + k*kk;
    div[ijk] = (cg0*u[ijk-ii1] + cg1*u[ijk] + cg2*u[ijk+ii1] + cg3*u[ijk+ii2]) * cgi*dxi
             + (cg0*v[ijk-jj1] + cg1*v[ijk] + cg2*v[ijk+jj1] + cg3*v[ijk+jj2]) * cgi*dyi
             + (cg0*w[ijk-kk1] + cg1*w[ijk] + cg2*w[ijk+kk1] + cg3*w[ijk+kk2]) * dzi4[k];
  }
}

#ifdef USECUDA
void cpres_4::exec(double dt)
{
  // 1. Create the input for the pressure solver.
  pres_in(fields->sd["p"]->data,
          fields->u ->data, fields->v ->data, fields->w ->data,
          fields->ut->data, fields->vt->data, fields->wt->data, 
          grid->dzi4, dt);

  // 2. Solve the Poisson equation using FFTs and a heptadiagonal solver
  // Take slices out of a temporary field to save memory. The temp arrays
  // are always big enough, this cannot fail.
  double *tmp2 = fields->sd["tmp2"]->data;
  double *tmp3 = fields->sd["tmp3"]->data;
  const int ns = grid->iblock*(grid->kmax+4);
  pres_solve(fields->sd["p"]->data, fields->sd["tmp1"]->data, grid->dz,
             m1, m2, m3, m4,
             m5, m6, m7,
             &tmp2[0*ns], &tmp2[1*ns], &tmp2[2*ns], &tmp2[3*ns], 
             &tmp3[0*ns], &tmp3[1*ns], &tmp3[2*ns], &tmp3[3*ns], 
             bmati, bmatj);

  // 3. Get the pressure tendencies from the pressure field.
  pres_out(fields->ut->data, fields->vt->data, fields->wt->data, 
           fields->sd["p"]->data, grid->dzhi4);
}

double cpres_4::check()
{
  fields->forwardGPU();

  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  const int offs = grid->memoffset;

  pres_4_calcdivergence<<<gridGPU, blockGPU>>>(&fields->a["tmp1"]->data_g[offs],
                                               &fields->u->data_g[offs], &fields->v->data_g[offs], &fields->w->data_g[offs],
                                               grid->dzi4_g,
                                               grid->dxi, grid->dyi,
                                               grid->icellsp, grid->ijcellsp,
                                               grid->istart,  grid->jstart, grid->kstart,
                                               grid->iend,    grid->jend,   grid->kend);

  double divmax = grid->getmax_g(&fields->a["tmp1"]->data_g[offs], fields->a["tmp2"]->data_g);
  grid->getmax(&divmax);

  fields->backwardGPU();

  return divmax;
}
#endif
