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

#include "fields.h"
#include "grid.h"
#include "master.h"
#include "boundary.h" // TMP BVS
#include "constants.h"
#include "tools.h"

// TODO use interp2 functions instead of manual interpolation
__global__ void fields_calcmom_2nd(double * __restrict__ u, double * __restrict__ v, double * __restrict__ w, 
                                   double * __restrict__ mom, double * __restrict__ dz,
                                   int istart, int jstart, int kstart,
                                   int iend,   int jend,   int kend,
                                   int jj, int kk)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
  int k = blockIdx.z + kstart; 
  int ii = 1;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    mom[ijk] = (0.5*(u[ijk]+u[ijk+ii]) + 0.5*(v[ijk]+v[ijk+jj]) + 0.5*(w[ijk]+w[ijk+kk]))*dz[k];
  }
}

__global__ void fields_calctke_2nd(double * __restrict__ u, double * __restrict__ v, double * __restrict__ w, 
                                   double * __restrict__ tke, double * __restrict__ dz,
                                   int istart, int jstart, int kstart,
                                   int iend,   int jend,   int kend,
                                   int jj, int kk)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
  int k = blockIdx.z + kstart; 
  int ii = 1;

  if(i < iend && j < jend && k < kend)
  {
    int ijk = i + j*jj + k*kk;
    tke[ijk] = ( 0.5*(pow(u[ijk],2)+pow(u[ijk+ii],2)) 
               + 0.5*(pow(v[ijk],2)+pow(v[ijk+jj],2)) 
               + 0.5*(pow(w[ijk],2)+pow(w[ijk+kk],2)))*dz[k];
  }
}

__global__ void fields_calcmass_2nd(double * __restrict__ s, double * __restrict__ mass, double * __restrict__ dz,
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
    mass[ijk] = s[ijk]*dz[k];
  }
}

#ifdef USECUDA
int Fields::exec()
{
  // calculate the means for the prognostic scalars
  if(calcprofs)
  {
    for(fieldmap::iterator it=sp.begin(); it!=sp.end(); ++it)
      grid->calcmean_g(it->second->datamean_g, &it->second->data_g[grid->memoffset], s["tmp1"]->data_g);
  }

  return 0;
}
#endif

#ifdef USECUDA
double Fields::checkmom()
{
  const int blocki = cuda::blockSizeI;
  const int blockj = cuda::blockSizeJ;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kcells);
  dim3 blockGPU(blocki, blockj, 1);

  const int offs = grid->memoffset;

  fields_calcmom_2nd<<<gridGPU, blockGPU>>>(&u->data_g[offs], &v->data_g[offs], &w->data_g[offs], 
                                            &a["tmp1"]->data_g[offs], grid->dz_g,
                                            grid->istart,  grid->jstart, grid->kstart,
                                            grid->iend,    grid->jend,   grid->kend,
                                            grid->icellsp, grid->ijcellsp);
  cudaCheckError();

  double mom = grid->getsum_g(&a["tmp1"]->data_g[offs], a["tmp2"]->data_g); 
  grid->getsum(&mom);
  mom /= (grid->itot*grid->jtot*grid->zsize);

  return mom;
}
#endif

#ifdef USECUDA
double Fields::checktke()
{
  const int blocki = cuda::blockSizeI;
  const int blockj = cuda::blockSizeJ;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kcells);
  dim3 blockGPU(blocki, blockj, 1);

  const int offs = grid->memoffset;

  fields_calctke_2nd<<<gridGPU, blockGPU>>>(&u->data_g[offs], &v->data_g[offs], &w->data_g[offs], 
                                            &a["tmp1"]->data_g[offs], grid->dz_g,
                                            grid->istart,  grid->jstart, grid->kstart,
                                            grid->iend,    grid->jend,   grid->kend,
                                            grid->icellsp, grid->ijcellsp);
  cudaCheckError();

  double tke = grid->getsum_g(&a["tmp1"]->data_g[offs], a["tmp2"]->data_g); 

  grid->getsum(&tke);
  tke /= (grid->itot*grid->jtot*grid->zsize);
  tke *= 0.5;

  return tke;
}
#endif

#ifdef USECUDA
double Fields::checkmass()
{
  const int blocki = cuda::blockSizeI;
  const int blockj = cuda::blockSizeJ;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kcells);
  dim3 blockGPU(blocki, blockj, 1);

  const int offs = grid->memoffset;
  double mass;

  // CvH for now, do the mass check on the first scalar... Do we want to change this?
  fieldmap::iterator itProg=sp.begin();
  if(sp.begin() != sp.end())
  {
    fields_calcmass_2nd<<<gridGPU, blockGPU>>>(&itProg->second->data_g[offs], &a["tmp1"]->data_g[offs], grid->dz_g,
                                               grid->istart,  grid->jstart, grid->kstart,
                                               grid->iend,    grid->jend,   grid->kend,
                                               grid->icellsp, grid->ijcellsp);
    cudaCheckError();

    mass = grid->getsum_g(&a["tmp1"]->data_g[offs], a["tmp2"]->data_g); 
    grid->getsum(&mass);
    mass /= (grid->itot*grid->jtot*grid->zsize);
  }
  else
    mass = 0; 

  return mass;
}
#endif

int Fields::prepareDevice()
{
  const int nmemsize   = grid->ncellsp*sizeof(double);
  const int nmemsize1d = grid->kcells*sizeof(double);
  const int nmemsize2d = (grid->ijcellsp+grid->memoffset)*sizeof(double);

  for(fieldmap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
  {
    cudaSafeCall(cudaMalloc(&it->second->data_g,        nmemsize  ));
    cudaSafeCall(cudaMalloc(&it->second->databot_g,     nmemsize2d));
    cudaSafeCall(cudaMalloc(&it->second->datatop_g,     nmemsize2d));
    cudaSafeCall(cudaMalloc(&it->second->datagradbot_g, nmemsize2d));
    cudaSafeCall(cudaMalloc(&it->second->datagradtop_g, nmemsize2d));
    cudaSafeCall(cudaMalloc(&it->second->datafluxbot_g, nmemsize2d));
    cudaSafeCall(cudaMalloc(&it->second->datafluxtop_g, nmemsize2d));
    cudaSafeCall(cudaMalloc(&it->second->datamean_g,    nmemsize1d));
  }

  // BvS: slightly wasteful, but make sure we have all fields...
  for(fieldmap::const_iterator it=sd.begin(); it!=sd.end(); ++it)
  {
    cudaSafeCall(cudaMalloc(&it->second->data_g,        nmemsize  ));
    cudaSafeCall(cudaMalloc(&it->second->databot_g,     nmemsize2d));
    cudaSafeCall(cudaMalloc(&it->second->datatop_g,     nmemsize2d));
    cudaSafeCall(cudaMalloc(&it->second->datagradbot_g, nmemsize2d));
    cudaSafeCall(cudaMalloc(&it->second->datagradtop_g, nmemsize2d));
    cudaSafeCall(cudaMalloc(&it->second->datafluxbot_g, nmemsize2d));
    cudaSafeCall(cudaMalloc(&it->second->datafluxtop_g, nmemsize2d));
    cudaSafeCall(cudaMalloc(&it->second->datamean_g,    nmemsize1d));
  }

  for(fieldmap::const_iterator it=at.begin(); it!=at.end(); ++it)
  {
    cudaSafeCall(cudaMalloc(&it->second->data_g, nmemsize));
  }

  cudaSafeCall(cudaMalloc(&rhoref_g,  nmemsize1d));
  cudaSafeCall(cudaMalloc(&rhorefh_g, nmemsize1d));

  // copy all the data to the GPU
  forwardDevice();

  return 0;
}

int Fields::forwardDevice()
{
  const int jcells     = grid->jcells;
  const int jkcells    = grid->jcells * grid->kcells;
  const int nmemsize1d = grid->kcells*sizeof(double);
  const int imemsizep  = grid->icellsp * sizeof(double);
  const int imemsize   = grid->icells  * sizeof(double);

  for(fieldmap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
  {
    cudaSafeCall(cudaMemcpy2D(&it->second->data_g[grid->memoffset],        imemsizep,  it->second->data,        imemsize, imemsize, jkcells, cudaMemcpyHostToDevice));
    cudaSafeCall(cudaMemcpy2D(&it->second->databot_g[grid->memoffset],     imemsizep,  it->second->databot,     imemsize, imemsize, jcells,  cudaMemcpyHostToDevice));
    cudaSafeCall(cudaMemcpy2D(&it->second->datatop_g[grid->memoffset],     imemsizep,  it->second->datatop,     imemsize, imemsize, jcells,  cudaMemcpyHostToDevice));
    cudaSafeCall(cudaMemcpy2D(&it->second->datagradbot_g[grid->memoffset], imemsizep,  it->second->datagradbot, imemsize, imemsize, jcells,  cudaMemcpyHostToDevice));
    cudaSafeCall(cudaMemcpy2D(&it->second->datagradtop_g[grid->memoffset], imemsizep,  it->second->datagradtop, imemsize, imemsize, jcells,  cudaMemcpyHostToDevice));
    cudaSafeCall(cudaMemcpy2D(&it->second->datafluxbot_g[grid->memoffset], imemsizep,  it->second->datafluxbot, imemsize, imemsize, jcells,  cudaMemcpyHostToDevice));
    cudaSafeCall(cudaMemcpy2D(&it->second->datafluxtop_g[grid->memoffset], imemsizep,  it->second->datafluxtop, imemsize, imemsize, jcells,  cudaMemcpyHostToDevice));
    cudaSafeCall(cudaMemcpy(it->second->datamean_g, it->second->datamean, nmemsize1d, cudaMemcpyHostToDevice));
  }

  for(fieldmap::const_iterator it=sd.begin(); it!=sd.end(); ++it)
  {
    cudaSafeCall(cudaMemcpy2D(&it->second->data_g[grid->memoffset],        imemsizep,  it->second->data,        imemsize, imemsize, jkcells, cudaMemcpyHostToDevice));
    cudaSafeCall(cudaMemcpy2D(&it->second->databot_g[grid->memoffset],     imemsizep,  it->second->databot,     imemsize, imemsize, jcells,  cudaMemcpyHostToDevice));
    cudaSafeCall(cudaMemcpy2D(&it->second->datatop_g[grid->memoffset],     imemsizep,  it->second->datatop,     imemsize, imemsize, jcells,  cudaMemcpyHostToDevice));
    cudaSafeCall(cudaMemcpy2D(&it->second->datagradbot_g[grid->memoffset], imemsizep,  it->second->datagradbot, imemsize, imemsize, jcells,  cudaMemcpyHostToDevice));
    cudaSafeCall(cudaMemcpy2D(&it->second->datagradtop_g[grid->memoffset], imemsizep,  it->second->datagradtop, imemsize, imemsize, jcells,  cudaMemcpyHostToDevice));
    cudaSafeCall(cudaMemcpy2D(&it->second->datafluxbot_g[grid->memoffset], imemsizep,  it->second->datafluxbot, imemsize, imemsize, jcells,  cudaMemcpyHostToDevice));
    cudaSafeCall(cudaMemcpy2D(&it->second->datafluxtop_g[grid->memoffset], imemsizep,  it->second->datafluxtop, imemsize, imemsize, jcells,  cudaMemcpyHostToDevice));
    cudaSafeCall(cudaMemcpy(it->second->datamean_g, it->second->datamean, nmemsize1d, cudaMemcpyHostToDevice));
  }

  for(fieldmap::const_iterator it=at.begin(); it!=at.end(); ++it)
    cudaSafeCall(cudaMemcpy2D(&it->second->data_g[grid->memoffset],        imemsizep, it->second->data,        imemsize, imemsize, jkcells, cudaMemcpyHostToDevice));

  cudaSafeCall(cudaMemcpy(rhoref_g,  rhoref,  nmemsize1d, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(rhorefh_g, rhorefh, nmemsize1d, cudaMemcpyHostToDevice));

  //master->printMessage("Synchronized GPU with CPU (forward)\n");

  return 0;
}

int Fields::backwardDevice()
{
  const int jcells     = grid->jcells;
  const int jkcells    = grid->jcells * grid->kcells;
  const int nmemsize1d = grid->kcells*sizeof(double);
  const int imemsizep  = grid->icellsp * sizeof(double);
  const int imemsize   = grid->icells  * sizeof(double);

  for(fieldmap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
  {
    cudaSafeCall(cudaMemcpy2D(it->second->data,        imemsize, &it->second->data_g[grid->memoffset],        imemsizep, imemsize, jkcells, cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaMemcpy2D(it->second->databot,     imemsize, &it->second->databot_g[grid->memoffset],     imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaMemcpy2D(it->second->datatop,     imemsize, &it->second->datatop_g[grid->memoffset],     imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaMemcpy2D(it->second->datagradbot, imemsize, &it->second->datagradbot_g[grid->memoffset], imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaMemcpy2D(it->second->datagradtop, imemsize, &it->second->datagradtop_g[grid->memoffset], imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaMemcpy2D(it->second->datafluxbot, imemsize, &it->second->datafluxbot_g[grid->memoffset], imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaMemcpy2D(it->second->datafluxtop, imemsize, &it->second->datafluxtop_g[grid->memoffset], imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaMemcpy(it->second->datamean, it->second->datamean_g, nmemsize1d, cudaMemcpyDeviceToHost));
  }

  for(fieldmap::const_iterator it=sd.begin(); it!=sd.end(); ++it)
  {
    cudaSafeCall(cudaMemcpy2D(it->second->data,        imemsize, &it->second->data_g[grid->memoffset],        imemsizep, imemsize, jkcells, cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaMemcpy2D(it->second->databot,     imemsize, &it->second->databot_g[grid->memoffset],     imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaMemcpy2D(it->second->datatop,     imemsize, &it->second->datatop_g[grid->memoffset],     imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaMemcpy2D(it->second->datagradbot, imemsize, &it->second->datagradbot_g[grid->memoffset], imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaMemcpy2D(it->second->datagradtop, imemsize, &it->second->datagradtop_g[grid->memoffset], imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaMemcpy2D(it->second->datafluxbot, imemsize, &it->second->datafluxbot_g[grid->memoffset], imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaMemcpy2D(it->second->datafluxtop, imemsize, &it->second->datafluxtop_g[grid->memoffset], imemsizep, imemsize, jcells,  cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaMemcpy(it->second->datamean, it->second->datamean_g, nmemsize1d, cudaMemcpyDeviceToHost));
  }

  for(fieldmap::const_iterator it=at.begin(); it!=at.end(); ++it)
    cudaSafeCall(cudaMemcpy2D(it->second->data,        imemsize, &it->second->data_g[grid->memoffset],        imemsizep, imemsize, jkcells, cudaMemcpyDeviceToHost));

  cudaSafeCall(cudaMemcpy(rhoref,  rhoref_g,  nmemsize1d, cudaMemcpyDeviceToHost));
  cudaSafeCall(cudaMemcpy(rhorefh, rhorefh_g, nmemsize1d, cudaMemcpyDeviceToHost));

  //master->printMessage("Synchronized CPU with GPU (backward)\n");

  return 0;
}

int Fields::clearDevice()
{
  for(fieldmap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
  {
    cudaSafeCall(cudaFree(it->second->data_g));
    cudaSafeCall(cudaFree(it->second->databot_g));
    cudaSafeCall(cudaFree(it->second->datatop_g));
    cudaSafeCall(cudaFree(it->second->datagradbot_g));
    cudaSafeCall(cudaFree(it->second->datagradtop_g));
    cudaSafeCall(cudaFree(it->second->datafluxbot_g));
    cudaSafeCall(cudaFree(it->second->datafluxtop_g));
    cudaSafeCall(cudaFree(it->second->datamean_g));
  }

  for(fieldmap::const_iterator it=sd.begin(); it!=sd.end(); ++it)
  {
    cudaSafeCall(cudaFree(it->second->data_g));
    cudaSafeCall(cudaFree(it->second->databot_g));
    cudaSafeCall(cudaFree(it->second->datatop_g));
    cudaSafeCall(cudaFree(it->second->datagradbot_g));
    cudaSafeCall(cudaFree(it->second->datagradtop_g));
    cudaSafeCall(cudaFree(it->second->datafluxbot_g));
    cudaSafeCall(cudaFree(it->second->datafluxtop_g));
    cudaSafeCall(cudaFree(it->second->datamean_g));
  }

  for(fieldmap::const_iterator it=at.begin(); it!=at.end(); ++it)
  {
    cudaSafeCall(cudaFree(it->second->data_g));
  }

  cudaSafeCall(cudaFree(rhoref_g));
  cudaSafeCall(cudaFree(rhorefh_g));

  return 0;
}

