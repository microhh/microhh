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

#include "fields.h"
#include "grid.h"
#include "master.h"
#include "boundary.h" // TMP BVS
#include "constants.h"
#include "tools.h"

namespace Fields_g
{
    // TODO use interp2 functions instead of manual interpolation
    __global__ 
    void calc_mom_2nd(double* __restrict__ u, double* __restrict__ v, double* __restrict__ w, 
                     double* __restrict__ mom, double* __restrict__ dz,
                     int istart, int jstart, int kstart,
                     int iend,   int jend,   int kend,
                     int jj,     int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
        const int k = blockIdx.z + kstart; 
        const int ii = 1;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            mom[ijk] = (0.5*(u[ijk]+u[ijk+ii]) + 0.5*(v[ijk]+v[ijk+jj]) + 0.5*(w[ijk]+w[ijk+kk]))*dz[k];
        }
    }

    __global__ 
    void calc_tke_2nd(double* __restrict__ u, double* __restrict__ v, double* __restrict__ w, 
                      double* __restrict__ tke, double* __restrict__ dz,
                      int istart, int jstart, int kstart,
                      int iend,   int jend,   int kend,
                      int jj,     int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
        const int k = blockIdx.z + kstart; 
        const int ii = 1;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            tke[ijk] = ( 0.5*(pow(u[ijk],2)+pow(u[ijk+ii],2)) 
                       + 0.5*(pow(v[ijk],2)+pow(v[ijk+jj],2)) 
                       + 0.5*(pow(w[ijk],2)+pow(w[ijk+kk],2)))*dz[k];
        }
    }

    __global__ 
    void calc_mass_2nd(double* __restrict__ s, double* __restrict__ mass, double* __restrict__ dz,
                       int istart, int jstart, int kstart,
                       int iend,   int jend,   int kend,
                       int jj,     int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
        const int k = blockIdx.z + kstart; 

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            mass[ijk] = s[ijk]*dz[k];
        }
    }
}

#ifdef USECUDA
void Fields::exec()
{
    // calculate the means for the prognostic scalars
    if (calc_mean_profs)
    {
        for (FieldMap::iterator it=sp.begin(); it!=sp.end(); ++it)
            grid->calcMean_g(it->second->datamean_g, &it->second->data_g[grid->memoffset], atmp["tmp1"]->data_g);
    }
}
#endif

#ifdef USECUDA
double Fields::check_momentum()
{
    const int blocki = grid->iThreadBlock;
    const int blockj = grid->jThreadBlock;
    const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
    const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, grid->kcells);
    dim3 blockGPU(blocki, blockj, 1);

    const int offs = grid->memoffset;

    Fields_g::calc_mom_2nd<<<gridGPU, blockGPU>>>(
        &u->data_g[offs], &v->data_g[offs], &w->data_g[offs], 
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
double Fields::check_tke()
{
    const int blocki = grid->iThreadBlock;
    const int blockj = grid->jThreadBlock;
    const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
    const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, grid->kcells);
    dim3 blockGPU(blocki, blockj, 1);

    const int offs = grid->memoffset;

    Fields_g::calc_tke_2nd<<<gridGPU, blockGPU>>>(
        &u->data_g[offs], &v->data_g[offs], &w->data_g[offs], 
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
double Fields::check_mass()
{
    const int blocki = grid->iThreadBlock;
    const int blockj = grid->jThreadBlock;
    const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
    const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, grid->kcells);
    dim3 blockGPU(blocki, blockj, 1);

    const int offs = grid->memoffset;
    double mass;

    // CvH for now, do the mass check on the first scalar... Do we want to change this?
    FieldMap::iterator itProg=sp.begin();
    if (sp.begin() != sp.end())
    {
        Fields_g::calc_mass_2nd<<<gridGPU, blockGPU>>>(
            &itProg->second->data_g[offs], &atmp["tmp1"]->data_g[offs], grid->dz_g,
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
void Fields::prepare_device()
{
    const int nmemsize   = grid->ncellsp*sizeof(double);
    const int nmemsize1d = grid->kcells*sizeof(double);

    // Prognostic fields
    for (FieldMap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
        it->second->init_device();

    // Diagnostic fields 
    for (FieldMap::const_iterator it=sd.begin(); it!=sd.end(); ++it)
        it->second->init_device();

    // Tendencies
    for (FieldMap::const_iterator it=at.begin(); it!=at.end(); ++it)
        cudaSafeCall(cudaMalloc(&it->second->data_g, nmemsize));

    // Temporary fields
    atmp["tmp1"]->init_device();
    atmp["tmp2"]->init_device();

    // Reference profiles
    cudaSafeCall(cudaMalloc(&rhoref_g,  nmemsize1d));
    cudaSafeCall(cudaMalloc(&rhorefh_g, nmemsize1d));

    // copy all the data to the GPU
    forward_device();
}

/**
 * This function deallocates all field3d instances and fields at device 
 */
void Fields::clear_device()
{
    for (FieldMap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
        it->second->clear_device();

    for (FieldMap::const_iterator it=sd.begin(); it!=sd.end(); ++it)
        it->second->clear_device();

    for (FieldMap::const_iterator it=at.begin(); it!=at.end(); ++it)
        cudaSafeCall(cudaFree(it->second->data_g));

    atmp["tmp1"]->clear_device();
    atmp["tmp2"]->clear_device();

    cudaSafeCall(cudaFree(rhoref_g));
    cudaSafeCall(cudaFree(rhorefh_g));
}

/**
 * This function copies all fields from host to device 
 */
void Fields::forward_device()
{
    for (FieldMap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
        forward_field3d_device(it->second);

    for (FieldMap::const_iterator it=sd.begin(); it!=sd.end(); ++it)
        forward_field3d_device(it->second);

    for (FieldMap::const_iterator it=at.begin(); it!=at.end(); ++it)
        forward_field_device_3d(it->second->data_g, it->second->data, Offset);

    forward_field3d_device(atmp["tmp1"]);
    forward_field3d_device(atmp["tmp2"]);

    forward_field_device_1d(rhoref_g,  rhoref , grid->kcells);
    forward_field_device_1d(rhorefh_g, rhorefh, grid->kcells);
}

/**
 * This function copies all fields required for statistics and output from device to host 
 */
void Fields::backward_device()
{
    for (FieldMap::const_iterator it=ap.begin(); it!=ap.end(); ++it)
        backward_field3d_device(it->second);

    for (FieldMap::const_iterator it=sd.begin(); it!=sd.end(); ++it)
        backward_field3d_device(it->second);

    //master->printMessage("Synchronized CPU with GPU (backward)\n");
}

/* BvS: it would make more sense to put this routine in field3d.cu, but how to solve this with the calls to fields.cu? */
/**
 * This function copies a field3d instance from host to device 
 * @param fld Pointer to field3d instance  
 */
void Fields::forward_field3d_device(Field3d* fld)
{
    forward_field_device_3d(fld->data_g,        fld->data,        Offset);
    forward_field_device_2d(fld->databot_g,     fld->databot,     Offset);
    forward_field_device_2d(fld->datatop_g,     fld->datatop,     Offset);
    forward_field_device_2d(fld->datagradbot_g, fld->datagradbot, Offset);
    forward_field_device_2d(fld->datagradtop_g, fld->datagradtop, Offset);
    forward_field_device_2d(fld->datafluxbot_g, fld->datafluxbot, Offset);
    forward_field_device_2d(fld->datafluxtop_g, fld->datafluxtop, Offset);
    forward_field_device_1d(fld->datamean_g,    fld->datamean,    grid->kcells);
}

/* BvS: it would make more sense to put this routine in field3d.cu, but how to solve this with the calls to fields.cu? */
/**
 * This function copies a field3d instance from device to host 
 * @param fld Pointer to field3d instance  
 */
void Fields::backward_field3d_device(Field3d* fld)
{
    backward_field_device_3d(fld->data,        fld->data_g,        Offset);
    backward_field_device_2d(fld->databot,     fld->databot_g,     Offset);
    backward_field_device_2d(fld->datatop,     fld->datatop_g,     Offset);
    backward_field_device_2d(fld->datagradbot, fld->datagradbot_g, Offset);
    backward_field_device_2d(fld->datagradtop, fld->datagradtop_g, Offset);
    backward_field_device_2d(fld->datafluxbot, fld->datafluxbot_g, Offset);
    backward_field_device_2d(fld->datafluxtop, fld->datafluxtop_g, Offset);
    backward_field_device_1d(fld->datamean,    fld->datamean_g,    grid->kcells);
}

/**
 * This function copies a single 3d field from host to device 
 * @param field_g Pointer to 3d field at device  
 * @param field Pointer to 3d field at host 
 * @param sw Switch to align the host field to device memory 
 */
void Fields::forward_field_device_3d(double* field_g, double* field, Offset_type sw)
{
    const int imemsizep  = grid->icellsp * sizeof(double);
    const int imemsize   = grid->icells  * sizeof(double);

    if (sw == Offset)
        cudaSafeCall(cudaMemcpy2D(&field_g[grid->memoffset], imemsizep,  field, imemsize, imemsize, grid->jcells*grid->kcells, cudaMemcpyHostToDevice));
    else if (sw == No_offset)
        cudaSafeCall(cudaMemcpy(field_g, field, grid->ncells*sizeof(double), cudaMemcpyHostToDevice));
}

/**
 * This function copies a single 2d field from host to device 
 * @param field_g Pointer to 2d field at device  
 * @param field Pointer to 2d field at host 
 * @param sw Switch to align the host field to device memory 
 */
void Fields::forward_field_device_2d(double* field_g, double* field, Offset_type sw)
{
    const int imemsizep  = grid->icellsp * sizeof(double);
    const int imemsize   = grid->icells  * sizeof(double);

    if (sw == Offset)
        cudaSafeCall(cudaMemcpy2D(&field_g[grid->memoffset], imemsizep,  field, imemsize, imemsize, grid->jcells,  cudaMemcpyHostToDevice));
    else if (sw == No_offset)
        cudaSafeCall(cudaMemcpy(field_g, field, grid->ijcells*sizeof(double), cudaMemcpyHostToDevice));
}

/**
 * This function copies an array from host to device 
 * @param field_g Pointer array at device  
 * @param field Pointer to array at host 
 * @param ncells Number of (double precision) values to copy 
 */
void Fields::forward_field_device_1d(double* field_g, double* field, int ncells)
{
    cudaSafeCall(cudaMemcpy(field_g, field, ncells*sizeof(double), cudaMemcpyHostToDevice));
}

/**
 * This function copies a single 3d field from device to host 
 * @param field Pointer to 3d field at host 
 * @param field_g Pointer to 3d field at device  
 * @param sw Switch to align the host field to device memory 
 */
void Fields::backward_field_device_3d(double* field, double* field_g, Offset_type sw)
{
    const int imemsizep  = grid->icellsp * sizeof(double);
    const int imemsize   = grid->icells  * sizeof(double);

    if (sw == Offset)
        cudaSafeCall(cudaMemcpy2D(field, imemsize, &field_g[grid->memoffset], imemsizep, imemsize, grid->jcells*grid->kcells, cudaMemcpyDeviceToHost));
    else if (sw == No_offset)
        cudaSafeCall(cudaMemcpy(field, field_g, grid->ncells*sizeof(double), cudaMemcpyDeviceToHost));
}

/**
 * This function copies a single 2d field from device to host 
 * @param field Pointer to 2d field at host 
 * @param field_g Pointer to 2d field at device  
 * @param sw Switch to align the host field to device memory 
 */
void Fields::backward_field_device_2d(double* field, double* field_g, Offset_type sw)
{
    const int imemsizep  = grid->icellsp * sizeof(double);
    const int imemsize   = grid->icells  * sizeof(double);

    if (sw == Offset)
        cudaSafeCall(cudaMemcpy2D(field, imemsize, &field_g[grid->memoffset], imemsizep, imemsize, grid->jcells, cudaMemcpyDeviceToHost));
    else if (sw == No_offset)
        cudaSafeCall(cudaMemcpy(field, field_g, grid->ijcells*sizeof(double), cudaMemcpyDeviceToHost));
}

/**
 * This function copies an array from device to host 
 * @param field Pointer to array at host 
 * @param field_g Pointer array at device  
 * @param ncells Number of (double precision) values to copy 
 */
void Fields::backward_field_device_1d(double* field, double* field_g, int ncells)
{
    cudaSafeCall(cudaMemcpy(field, field_g, ncells*sizeof(double), cudaMemcpyDeviceToHost));
}
