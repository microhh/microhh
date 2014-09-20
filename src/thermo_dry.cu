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
#include "grid.h"
#include "fields.h"
#include "thermo_dry.h"
#include "defines.h"
#include "constants.h"
#include "master.h"

__global__ void thermo_dry_calcbuoyancytend_2nd(double * __restrict__ wt, double * __restrict__ th,
                                                double * __restrict__ threfh, double grav, 
                                                int istart, int jstart, int kstart,
                                                int iend,   int jend,   int kend,
                                                int jj, int kk)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
  int k = blockIdx.z + kstart; 

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;

    wt[ijk] += grav/threfh[k] * (0.5*(th[ijk-kk]+th[ijk]) - threfh[k]);
  }
}


__global__ void thermo_dry_calcbuoyancy(double * __restrict__ b, double * __restrict__ th,
                                        double * __restrict__ thref, 
                                        int istart, int jstart,
                                        int iend,   int jend,   int kcells,
                                        int jj, int kk)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
  int k = blockIdx.z; 

  if(i < iend && j < jend && k < kcells)
  {
    int ijk = i + j*jj + k*kk;
    b[ijk] = constants::grav/thref[k] * (th[ijk] - thref[k]);
  }
}

__global__ void thermo_dry_calcbuoyancybot(double * __restrict__ b,     double * __restrict__ bbot,
                                           double * __restrict__ th,    double * __restrict__ thbot, 
                                           double * __restrict__ thref, double * __restrict__ threfh,
                                           double grav, int kstart, int icells, int jcells,  
                                           int jj, int kk)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x; 
  int j = blockIdx.y*blockDim.y + threadIdx.y; 

  if(i < icells && j < jcells)
  {
    const int ij  = i + j*jj;
    const int ijk = i + j*jj + kstart*kk;

    bbot[ij] = grav/threfh[kstart] * (thbot[ij] - threfh[kstart]);
    b[ijk]   = grav/thref [kstart] * (th[ijk]   - thref [kstart]);
  }
}

__global__ void thermo_dry_calcbuoyancyfluxbot(double * __restrict__ bfluxbot, double * __restrict__ thfluxbot,
                                               double * __restrict__ threfh, 
                                               double grav, int kstart, int icells, int jcells,  
                                               int jj, int kk)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x; 
  int j = blockIdx.y*blockDim.y + threadIdx.y; 

  if(i < icells && j < jcells)
  {
    const int ij  = i + j*jj;
    bfluxbot[ij] = grav/threfh[kstart]*thfluxbot[ij];
  }
}

__global__ void thermo_dry_calcN2(double * __restrict__ N2, double * __restrict__ th,
                                  double * __restrict__ thref, double * __restrict__ dzi, 
                                  int istart, int jstart,
                                  int iend,   int jend,   int kcells,
                                  int jj, int kk)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
  int k = blockIdx.z; 

  if(i < iend && j < jend && k < kcells)
  {
    int ijk = i + j*jj + k*kk;
    N2[ijk] = constants::grav/thref[k]*0.5*(th[ijk+kk] - th[ijk-kk])*dzi[k];
  }
}

int cthermo_dry::prepareDevice()
{
  const int nmemsize = grid->kcells*sizeof(double);

  // Allocate fields for Boussinesq and anelastic solver
  cudaMalloc(&thref_g,  nmemsize);
  cudaMalloc(&threfh_g, nmemsize);
  cudaMalloc(&pref_g,   nmemsize);
  cudaMalloc(&prefh_g,  nmemsize);
  cudaMalloc(&exner_g,  nmemsize);
  cudaMalloc(&exnerh_g, nmemsize);

  // Copy fields to device
  cudaMemcpy(thref_g,  thref,  nmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(threfh_g, threfh, nmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(pref_g,   pref,   nmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(prefh_g,  prefh,  nmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(exner_g,  exner,  nmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(exnerh_g, exnerh, nmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(thref_g,  thref,  nmemsize, cudaMemcpyHostToDevice);

  return 0;
}

#ifdef USECUDA
int cthermo_dry::exec()
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);
  
  const int offs = grid->memoffset;

  if(grid->swspatialorder== "2")
    thermo_dry_calcbuoyancytend_2nd<<<gridGPU, blockGPU>>>(&fields->wt->data_g[offs], &fields->s["th"]->data_g[offs], threfh_g, constants::grav, 
                                                           grid->istart, grid->jstart, grid->kstart+1,
                                                           grid->iend,   grid->jend, grid->kend,
                                                           grid->icellsp, grid->ijcellsp);

  //else if(grid->swspatialorder == "4")
  //  calcbuoyancytend_4th(fields->wt->data, fields->s["th"]->data, threfh);

  return 0;
}
#endif

#ifdef USECUDA
int cthermo_dry::getthermofield(cfield3d *fld, cfield3d *tmp, std::string name)
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kcells);
  dim3 blockGPU(blocki, blockj, 1);
  
  const int offs = grid->memoffset;

  if(name == "b")
    thermo_dry_calcbuoyancy<<<gridGPU, blockGPU>>>(&fld->data_g[offs], &fields->s["th"]->data_g[offs], 
                                                   thref_g, grid->istart, grid->jstart, grid->iend, grid->jend, grid->kcells,
                                                   grid->icellsp, grid->ijcellsp);
  else if(name == "N2")
    thermo_dry_calcN2<<<gridGPU, blockGPU>>>(&fld->data_g[offs], &fields->s["th"]->data_g[offs], 
                                             thref_g, grid->dzi_g, grid->istart, grid->jstart, grid->iend, grid->jend, grid->kcells,
                                             grid->icellsp, grid->ijcellsp);
  else
    return 1;

  return 0;
}
#endif

#ifdef USECUDA
int cthermo_dry::getbuoyancyfluxbot(cfield3d *bfield)
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->icells/blocki + (grid->icells%blocki > 0);
  const int gridj  = grid->jcells/blockj + (grid->jcells%blockj > 0);

  dim3 gridGPU (gridi, gridj, 1);
  dim3 blockGPU(blocki, blockj, 1);
  
  const int offs = grid->memoffset;

  thermo_dry_calcbuoyancyfluxbot<<<gridGPU, blockGPU>>>(&bfield->datafluxbot_g[offs], &fields->s["th"]->datafluxbot_g[offs], 
                                                        threfh_g, constants::grav, grid->kstart, grid->icells, grid->jcells, 
                                                        grid->icellsp, grid->ijcellsp);

  return 0;
}
#endif

#ifdef USECUDA
int cthermo_dry::getbuoyancysurf(cfield3d *bfield)
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->icells/blocki + (grid->icells%blocki > 0);
  const int gridj  = grid->jcells/blockj + (grid->jcells%blockj > 0);

  dim3 gridGPU (gridi, gridj, 1);
  dim3 blockGPU(blocki, blockj, 1);
  
  const int offs = grid->memoffset;

  thermo_dry_calcbuoyancybot<<<gridGPU, blockGPU>>>(&bfield->data_g[offs], &bfield->databot_g[offs], 
                                                    &fields->s["th"]->data_g[offs], &fields->s["th"]->databot_g[offs],
                                                    thref_g, threfh_g, constants::grav, grid->kstart, grid->icells, grid->jcells, 
                                                    grid->icellsp, grid->ijcellsp);

  thermo_dry_calcbuoyancyfluxbot<<<gridGPU, blockGPU>>>(&bfield->datafluxbot_g[offs], &fields->s["th"]->datafluxbot_g[offs], 
                                                        threfh_g, constants::grav, grid->kstart, grid->icells, grid->jcells, 
                                                        grid->icellsp, grid->ijcellsp);

  return 0;
}
#endif
