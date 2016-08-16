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
#include "grid.h"
#include "fields.h"
#include "thermo_dry.h"
#include "defines.h"
#include "constants.h"
#include "master.h"
#include "tools.h"

namespace
{
    __global__ 
    void calc_buoyancy_tend_2nd_g(double* __restrict__ wt, 
                                  double* __restrict__ th, double* __restrict__ threfh, 
                                  int istart, int jstart, int kstart,
                                  int iend,   int jend,   int kend,
                                  int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
        const int k = blockIdx.z + kstart; 

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            wt[ijk] += Constants::grav/threfh[k] * (0.5*(th[ijk-kk]+th[ijk]) - threfh[k]);
        }
    }


    __global__ 
    void calc_buoyancy_g(double* __restrict__ b,
                         double* __restrict__ th, double* __restrict__ thref, 
                         int istart, int jstart,
                         int iend,   int jend,   int kcells,
                         int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
        const int k = blockIdx.z; 

        if (i < iend && j < jend && k < kcells)
        {
            const int ijk = i + j*jj + k*kk;
            b[ijk] = Constants::grav/thref[k] * (th[ijk] - thref[k]);
        }
    }

    __global__ 
    void calc_buoyancy_bot_g(double* __restrict__ b,     double* __restrict__ bbot,
                             double* __restrict__ th,    double* __restrict__ thbot, 
                             double* __restrict__ thref, double* __restrict__ threfh,
                             double grav, int kstart, int icells, int jcells,  
                             int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y; 

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;

            bbot[ij] = grav/threfh[kstart] * (thbot[ij] - threfh[kstart]);
            b[ijk]   = grav/thref [kstart] * (th[ijk]   - thref [kstart]);
        }
    }

    __global__ 
    void calc_buoyancy_flux_bot_g(double* __restrict__ bfluxbot, double* __restrict__ thfluxbot,
                                  double* __restrict__ threfh, 
                                  double grav, int kstart, int icells, int jcells,  
                                  int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y; 

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;
            bfluxbot[ij] = grav/threfh[kstart]*thfluxbot[ij];
        }
    }

    __global__ 
    void calc_N2_g(double* __restrict__ N2,    double* __restrict__ th,
                   double* __restrict__ thref, double* __restrict__ dzi, 
                   int istart, int jstart, int kstart,
                   int iend,   int jend,   int kend,
                   int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
        const int k = blockIdx.z + kstart; 

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            N2[ijk] = Constants::grav/thref[k]*0.5*(th[ijk+kk] - th[ijk-kk])*dzi[k];
        }
    }
} // end namespace

void Thermo_dry::prepare_device()
{
    const int nmemsize = grid->kcells*sizeof(double);

    // Allocate fields for Boussinesq and anelastic solver
    cuda_safe_call(cudaMalloc(&thref_g,   nmemsize));
    cuda_safe_call(cudaMalloc(&threfh_g,  nmemsize));
    cuda_safe_call(cudaMalloc(&pref_g,    nmemsize));
    cuda_safe_call(cudaMalloc(&prefh_g,   nmemsize));
    cuda_safe_call(cudaMalloc(&exnref_g,  nmemsize));
    cuda_safe_call(cudaMalloc(&exnrefh_g, nmemsize));

    // Copy fields to device
    cuda_safe_call(cudaMemcpy(thref_g,   thref,   nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(threfh_g,  threfh,  nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(pref_g,    pref,    nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(prefh_g,   prefh,   nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(exnref_g,  exnref,  nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(exnrefh_g, exnrefh, nmemsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(thref_g,   thref,   nmemsize, cudaMemcpyHostToDevice));
}

void Thermo_dry::clear_device()
{
    cuda_safe_call(cudaFree(thref_g ));
    cuda_safe_call(cudaFree(threfh_g));
    cuda_safe_call(cudaFree(pref_g  ));
    cuda_safe_call(cudaFree(prefh_g ));
    cuda_safe_call(cudaFree(exnref_g ));
    cuda_safe_call(cudaFree(exnrefh_g));
}

#ifdef USECUDA
void Thermo_dry::exec()
{
    const int blocki = grid->ithread_block;
    const int blockj = grid->jthread_block;
    const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
    const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, grid->kmax);
    dim3 blockGPU(blocki, blockj, 1);

    const int offs = grid->memoffset;

    if (grid->swspatialorder== "2")
    {
        calc_buoyancy_tend_2nd_g<<<gridGPU, blockGPU>>>(
            &fields->wt->data_g[offs], &fields->sp["th"]->data_g[offs], threfh_g, 
            grid->istart,  grid->jstart, grid->kstart+1,
            grid->iend,    grid->jend,   grid->kend,
            grid->icellsp, grid->ijcellsp);

        cuda_check_error();
    }
    else if (grid->swspatialorder == "4")
    {
        master->print_message("4th order thermo_dry not (yet) implemented\n");  
        throw 1;
    }
}
#endif

#ifdef USECUDA
void Thermo_dry::get_thermo_field(Field3d *fld, Field3d *tmp, std::string name, bool cyclic)
{
    const int blocki = grid->ithread_block;
    const int blockj = grid->jthread_block;
    const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
    const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, grid->kcells);
    dim3 blockGPU(blocki, blockj, 1);

    dim3 gridGPU2 (gridi, gridj, grid->kmax);
    dim3 blockGPU2(blocki, blockj, 1);

    const int offs = grid->memoffset;

    if (name == "b")
    {
        calc_buoyancy_g<<<gridGPU, blockGPU>>>(
            &fld->data_g[offs], &fields->sp["th"]->data_g[offs], thref_g, 
            grid->istart, grid->jstart, 
            grid->iend, grid->jend, grid->kcells,
            grid->icellsp, grid->ijcellsp);
        cuda_check_error();
    }
    else if (name == "N2")
    {
        calc_N2_g<<<gridGPU2, blockGPU2>>>(
            &fld->data_g[offs], &fields->sp["th"]->data_g[offs], thref_g, grid->dzi_g, 
            grid->istart,  grid->jstart, grid->kstart, 
            grid->iend,    grid->jend,   grid->kend,
            grid->icellsp, grid->ijcellsp);
        cuda_check_error();
    }
    else
    {
        master->print_error("get_thermo_field \"%s\" not supported\n",name.c_str());
        throw 1;
    }

    if (cyclic)
        grid->boundary_cyclic_g(&fld->data_g[offs]);
}
#endif

#ifdef USECUDA
void Thermo_dry::get_buoyancy_fluxbot(Field3d *bfield)
{
    const int blocki = grid->ithread_block;
    const int blockj = grid->jthread_block;
    const int gridi  = grid->icells/blocki + (grid->icells%blocki > 0);
    const int gridj  = grid->jcells/blockj + (grid->jcells%blockj > 0);

    dim3 gridGPU (gridi, gridj, 1);
    dim3 blockGPU(blocki, blockj, 1);

    const int offs = grid->memoffset;

    calc_buoyancy_flux_bot_g<<<gridGPU, blockGPU>>>(
        &bfield->datafluxbot_g[offs], &fields->sp["th"]->datafluxbot_g[offs], 
        threfh_g, Constants::grav, grid->kstart, grid->icells, grid->jcells, 
        grid->icellsp, grid->ijcellsp);
    cuda_check_error();
}
#endif

#ifdef USECUDA
void Thermo_dry::get_buoyancy_surf(Field3d *bfield)
{
    const int blocki = grid->ithread_block;
    const int blockj = grid->jthread_block;
    const int gridi  = grid->icells/blocki + (grid->icells%blocki > 0);
    const int gridj  = grid->jcells/blockj + (grid->jcells%blockj > 0);

    dim3 gridGPU (gridi, gridj, 1);
    dim3 blockGPU(blocki, blockj, 1);

    const int offs = grid->memoffset;

    calc_buoyancy_bot_g<<<gridGPU, blockGPU>>>(
        &bfield->data_g[offs], &bfield->databot_g[offs], 
        &fields->sp["th"]->data_g[offs], &fields->sp["th"]->databot_g[offs],
        thref_g, threfh_g, Constants::grav, grid->kstart, grid->icells, grid->jcells, 
        grid->icellsp, grid->ijcellsp);
    cuda_check_error();

    calc_buoyancy_flux_bot_g<<<gridGPU, blockGPU>>>(
        &bfield->datafluxbot_g[offs], &fields->sp["th"]->datafluxbot_g[offs], 
        threfh_g, Constants::grav, grid->kstart, grid->icells, grid->jcells, 
        grid->icellsp, grid->ijcellsp);
    cuda_check_error();
}
#endif
