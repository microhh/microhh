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

namespace Fields_g
{
  // TODO use interp2 functions instead of manual interpolation
  __global__ void calcmom_2nd(double * __restrict__ u, double * __restrict__ v, double * __restrict__ w, 
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
  
  __global__ void calctke_2nd(double * __restrict__ u, double * __restrict__ v, double * __restrict__ w, 
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
  
  __global__ void calcmass_2nd(double * __restrict__ s, double * __restrict__ mass, double * __restrict__ dz,
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
}

#ifdef USECUDA
int Fields::exec()
{
  // calculate the means for the prognostic scalars
  if(calcMeanProfs)
  {
    for(FieldMap::iterator it=sp.begin(); it!=sp.end(); ++it)
      grid->calcMean_g(it->second->datamean_g, &it->second->data_g[grid->memoffset], atmp["tmp1"]->data_g);
  }

  return 0;
}
#endif

#ifdef USECUDA
double Fields::checkMomentum()
{
  const int blocki = cuda::blockSizeI;
  const int blockj = cuda::blockSizeJ;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kcells);
  dim3 blockGPU(blocki, blockj, 1);

  const int offs = grid->memoffset;

  Fields_g::calcmom_2nd<<<gridGPU, blockGPU>>>(&u->data_g[offs], &v->data_g[offs], &w->data_g[offs], 
                                               &atmp["tmp1"]->data_g[offs], grid->dz_g,
                                               grid->istart,  grid->jstart, grid->kstart,
                                               grid->iend,    grid->jend,   grid->kend,
                                               grid->icellsp, grid->ijcellsp);
  cudaCheckError();

  double mom = grid->getSum_g(&atmp["tmp1"]->data_g[offs], atmp["tmp2"]->data_g); 
  grid->getSum(&mom);
  mom /= (grid->itot*grid->jtot*grid->zsize);

  return mom;
}
#endif

#ifdef USECUDA
double Fields::checkTke()
{
  const int blocki = cuda::blockSizeI;
  const int blockj = cuda::blockSizeJ;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kcells);
  dim3 blockGPU(blocki, blockj, 1);

  const int offs = grid->memoffset;

  Fields_g::calctke_2nd<<<gridGPU, blockGPU>>>(&u->data_g[offs], &v->data_g[offs], &w->data_g[offs], 
                                               &atmp["tmp1"]->data_g[offs], grid->dz_g,
                                               grid->istart,  grid->jstart, grid->kstart,
                                               grid->iend,    grid->jend,   grid->kend,
                                               grid->icellsp, grid->ijcellsp);
  cudaCheckError();

  double tke = grid->getSum_g(&atmp["tmp1"]->data_g[offs], atmp["tmp2"]->data_g); 

  grid->getSum(&tke);
  tke /= (grid->itot*grid->jtot*grid->zsize);
  tke *= 0.5;

  return tke;
}
#endif

#ifdef USECUDA
double Fields::checkMass()
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
  FieldMap::iterator itProg=sp.begin();
  if(sp.begin() != sp.end())
  {
    Fields_g::calcmass_2nd<<<gridGPU, blockGPU>>>(&itProg->second->data_g[offs], &atmp["tmp1"]->data_g[offs], grid->dz_g,
                                                  grid->istart,  grid->jstart, grid->kstart,
                                                  grid->iend,    grid->jend,   grid->kend,
                                                  grid->icellsp, grid->ijcellsp);
    cudaCheckError();

    mass = grid->getSum_g(&atmp["tmp1"]->data_g[offs], atmp["tmp2"]->data_g); 
    grid->getSum(&mass);
    mass /= (grid->itot*grid->jtot*grid->zsize);
  }
  else
    mass = 0; 

  return mass;
}
#endif

/**
 * This function allocates all field3d instances and fields at device 
 */
void Fields::prepareDevice()
{
  const int nmemsize   = grid->ncellsp*sizeof(double);
  const int nmemsize1d = grid->kcells*sizeof(double);

  // Prognostic fields
  for(FieldMap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
    it->second->initDevice();
 
  // Diagnostic fields 
  for(FieldMap::const_iterator it=sd.begin(); it!=sd.end(); ++it)
    it->second->initDevice();

  // Tendencies
  for(FieldMap::const_iterator it=at.begin(); it!=at.end(); ++it)
    cudaSafeCall(cudaMalloc(&it->second->data_g, nmemsize));

  // Temporary fields
  atmp["tmp1"]->initDevice();
  atmp["tmp2"]->initDevice();

  // Reference profiles
  cudaSafeCall(cudaMalloc(&rhoref_g,  nmemsize1d));
  cudaSafeCall(cudaMalloc(&rhorefh_g, nmemsize1d));

  // copy all the data to the GPU
  forwardDevice();
}

/**
 * This function deallocates all field3d instances and fields at device 
 */
void Fields::clearDevice()
{
  for(FieldMap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
    it->second->clearDevice();

  for(FieldMap::const_iterator it=sd.begin(); it!=sd.end(); ++it)
    it->second->clearDevice();

  for(FieldMap::const_iterator it=at.begin(); it!=at.end(); ++it)
    cudaSafeCall(cudaFree(it->second->data_g));

  atmp["tmp1"]->clearDevice();
  atmp["tmp2"]->clearDevice();
  
  cudaSafeCall(cudaFree(rhoref_g));
  cudaSafeCall(cudaFree(rhorefh_g));
}

/**
 * This function copies all fields from host to device 
 */
void Fields::forwardDevice()
{
  for(FieldMap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
    forwardField3dDevice(it->second);

  for(FieldMap::const_iterator it=sd.begin(); it!=sd.end(); ++it)
    forwardField3dDevice(it->second);

  for(FieldMap::const_iterator it=at.begin(); it!=at.end(); ++it)
    forward3DFieldDevice(it->second->data_g, it->second->data, Offset);

  forwardField3dDevice(atmp["tmp1"]);
  forwardField3dDevice(atmp["tmp2"]);

  forward1DFieldDevice(rhoref_g,  rhoref , grid->kcells);
  forward1DFieldDevice(rhorefh_g, rhorefh, grid->kcells);

  //master->printMessage("Synchronized GPU with CPU (forward)\n");
}

/**
 * This function copies all fields required for statistics and output from device to host 
 */
void Fields::backwardDevice()
{
  for(FieldMap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
    backwardField3dDevice(it->second);

  for(FieldMap::const_iterator it=sd.begin(); it!=sd.end(); ++it)
    backwardField3dDevice(it->second);

  //for(FieldMap::const_iterator it=at.begin(); it!=at.end(); ++it)
  //  backward3DFieldDevice(it->second->data,        it->second->data_g,        Offset);

  //backwardField3dDevice(atmp["tmp1"]);
  //backwardField3dDevice(atmp["tmp2"]);

  //backward1DFieldDevice(rhoref,  rhoref_g,  grid->kcells);
  //backward1DFieldDevice(rhorefh, rhorefh_g, grid->kcells);

  //master->printMessage("Synchronized CPU with GPU (backward)\n");
}

/* BvS: it would make more sense to put this routine in field3d.cu, but how to solve this with the calls to fields.cu? */
/**
 * This function copies a field3d instance from host to device 
 * @param fld Pointer to field3d instance  
 */
void Fields::forwardField3dDevice(Field3d *fld)
{
  forward3DFieldDevice(fld->data_g,        fld->data,        Offset);
  forward2DFieldDevice(fld->databot_g,     fld->databot,     Offset);
  forward2DFieldDevice(fld->datatop_g,     fld->datatop,     Offset);
  forward2DFieldDevice(fld->datagradbot_g, fld->datagradbot, Offset);
  forward2DFieldDevice(fld->datagradtop_g, fld->datagradtop, Offset);
  forward2DFieldDevice(fld->datafluxbot_g, fld->datafluxbot, Offset);
  forward2DFieldDevice(fld->datafluxtop_g, fld->datafluxtop, Offset);
  forward1DFieldDevice(fld->datamean_g,    fld->datamean, grid->kcells);
}

/* BvS: it would make more sense to put this routine in field3d.cu, but how to solve this with the calls to fields.cu? */
/**
 * This function copies a field3d instance from device to host 
 * @param fld Pointer to field3d instance  
 */
void Fields::backwardField3dDevice(Field3d *fld)
{
  backward3DFieldDevice(fld->data,        fld->data_g,        Offset);
  backward2DFieldDevice(fld->databot,     fld->databot_g,     Offset);
  backward2DFieldDevice(fld->datatop,     fld->datatop_g,     Offset);
  backward2DFieldDevice(fld->datagradbot, fld->datagradbot_g, Offset);
  backward2DFieldDevice(fld->datagradtop, fld->datagradtop_g, Offset);
  backward2DFieldDevice(fld->datafluxbot, fld->datafluxbot_g, Offset);
  backward2DFieldDevice(fld->datafluxtop, fld->datafluxtop_g, Offset);
  backward1DFieldDevice(fld->datamean,    fld->datamean_g, grid->kcells);
}

/**
 * This function copies a single 3d field from host to device 
 * @param field_g Pointer to 3d field at device  
 * @param field Pointer to 3d field at host 
 * @param sw Switch to align the host field to device memory 
 */
void Fields::forward3DFieldDevice(double * field_g, double * field, OffsetType sw)
{
  const int imemsizep  = grid->icellsp * sizeof(double);
  const int imemsize   = grid->icells  * sizeof(double);

  if(sw == Offset)
    cudaSafeCall(cudaMemcpy2D(&field_g[grid->memoffset], imemsizep,  field, imemsize, imemsize, grid->jcells*grid->kcells, cudaMemcpyHostToDevice));
  else if(sw == NoOffset)
    cudaSafeCall(cudaMemcpy(field_g, field, grid->ncells*sizeof(double), cudaMemcpyHostToDevice));
}

/**
 * This function copies a single 2d field from host to device 
 * @param field_g Pointer to 2d field at device  
 * @param field Pointer to 2d field at host 
 * @param sw Switch to align the host field to device memory 
 */
void Fields::forward2DFieldDevice(double * field_g, double * field, OffsetType sw)
{
  const int imemsizep  = grid->icellsp * sizeof(double);
  const int imemsize   = grid->icells  * sizeof(double);
 
  if(sw == Offset)
    cudaSafeCall(cudaMemcpy2D(&field_g[grid->memoffset], imemsizep,  field, imemsize, imemsize, grid->jcells,  cudaMemcpyHostToDevice));
  else if(sw == NoOffset)
    cudaSafeCall(cudaMemcpy(field_g, field, grid->ijcells*sizeof(double), cudaMemcpyHostToDevice));
}

/**
 * This function copies an array from host to device 
 * @param field_g Pointer array at device  
 * @param field Pointer to array at host 
 * @param ncells Number of (double precision) values to copy 
 */
void Fields::forward1DFieldDevice(double * field_g, double * field, int ncells)
{
  cudaSafeCall(cudaMemcpy(field_g, field, ncells*sizeof(double), cudaMemcpyHostToDevice));
}

/**
 * This function copies a single 3d field from device to host 
 * @param field Pointer to 3d field at host 
 * @param field_g Pointer to 3d field at device  
 * @param sw Switch to align the host field to device memory 
 */
void Fields::backward3DFieldDevice(double * field, double * field_g, OffsetType sw)
{
  const int imemsizep  = grid->icellsp * sizeof(double);
  const int imemsize   = grid->icells  * sizeof(double);

  if(sw == Offset)
    cudaSafeCall(cudaMemcpy2D(field, imemsize, &field_g[grid->memoffset], imemsizep, imemsize, grid->jcells*grid->kcells, cudaMemcpyDeviceToHost));
  else if(sw == NoOffset)
    cudaSafeCall(cudaMemcpy(field, field_g, grid->ncells*sizeof(double), cudaMemcpyDeviceToHost));
}

/**
 * This function copies a single 2d field from device to host 
 * @param field Pointer to 2d field at host 
 * @param field_g Pointer to 2d field at device  
 * @param sw Switch to align the host field to device memory 
 */
void Fields::backward2DFieldDevice(double * field, double * field_g, OffsetType sw)
{
  const int imemsizep  = grid->icellsp * sizeof(double);
  const int imemsize   = grid->icells  * sizeof(double);
 
  if(sw == Offset)
    cudaSafeCall(cudaMemcpy2D(field, imemsize, &field_g[grid->memoffset], imemsizep, imemsize, grid->jcells, cudaMemcpyDeviceToHost));
  else if(sw == NoOffset)
    cudaSafeCall(cudaMemcpy(field, field_g, grid->ijcells*sizeof(double), cudaMemcpyDeviceToHost));
}

/**
 * This function copies an array from device to host 
 * @param field Pointer to array at host 
 * @param field_g Pointer array at device  
 * @param ncells Number of (double precision) values to copy 
 */
void Fields::backward1DFieldDevice(double * field, double * field_g, int ncells)
{
  cudaSafeCall(cudaMemcpy(field, field_g, ncells*sizeof(double), cudaMemcpyDeviceToHost));
}
