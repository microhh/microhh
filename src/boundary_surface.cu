/*
 * MicroHH
 * Copyright (c) 2011-2018 Chiel van Heerwaarden
 * Copyright (c) 2011-2018 Thijs Heus
 * Copyright (c) 2014-2018 Bart van Stratum
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
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "boundary_surface.h"
#include "defines.h"
#include "constants.h"
#include "thermo.h"
#include "model.h"
#include "master.h"
#include "tools.h"
#include "monin_obukhov.h"

namespace
{
    namespace most = Monin_obukhov;
    const int nzL = 10000; // Size of the lookup table for MO iterations.

    template<typename TF>
    __device__ 
    double find_Obuk_g(const float* const __restrict__ zL, const float* const __restrict__ f,
                       int &n, const double Ri, const double zsl)
    {
        // Determine search direction.
        if ((f[n]-Ri) > 0)
            while ( (f[n-1]-Ri) > 0 && n > 0) { --n; }
        else
            while ( (f[n]-Ri) < 0 && n < (nzL-1) ) { ++n; }

        const double zL0 = (n == 0 || n == nzL-1) ? zL[n] : zL[n-1] + (Ri-f[n-1]) / (f[n]-f[n-1]) * (zL[n]-zL[n-1]);

        return zsl/zL0;
    }


    template<typename TF>
    __device__ 
    double calc_Obuk_noslip_flux_g(float* __restrict__ zL, float* __restrict__ f, int& n, double du, double bfluxbot, double zsl)
    {
        // Calculate the appropriate Richardson number.
        const double Ri = -Constants::kappa * bfluxbot * zsl / pow(du, 3);
        return find_Obuk_g(zL, f, n, Ri, zsl);
    }

    template<typename TF>
    __device__ 
    double calc_Obuk_noslip_dirichlet_g(float* __restrict__ zL, float* __restrict__ f, int& n, double du, double db, double zsl)
    {
        // Calculate the appropriate Richardson number.
        const double Ri = Constants::kappa * db * zsl / pow(du, 2);
        return find_Obuk_g(zL, f, n, Ri, zsl);
    }

    /* Calculate absolute wind speed */
    template<typename TF>
    __global__ 
    void du_tot_g(double* __restrict__ dutot, 
                  double* __restrict__ u,    double* __restrict__ v,
                  double* __restrict__ ubot, double* __restrict__ vbot, 
                  int istart, int jstart, int kstart,
                  int iend,   int jend, int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 

        if (i < iend && j < jend)
        {
            const int ii  = 1;
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            const double minval = 1.e-1;

            const double du2 = pow(0.5*(u[ijk] + u[ijk+ii]) - 0.5*(ubot[ij] + ubot[ij+ii]), 2)
                             + pow(0.5*(v[ijk] + v[ijk+jj]) - 0.5*(vbot[ij] + vbot[ij+jj]), 2);
            dutot[ij] = fmax(pow(du2, 0.5), minval);
        }
    }

    template<typename TF>
    __global__ 
    void stability_g(double* __restrict__ ustar, double* __restrict__ obuk,
                     double* __restrict__ b, double* __restrict__ bbot, double* __restrict__ bfluxbot,
                     double* __restrict__ dutot, float* __restrict__ zL_sl_g, float* __restrict__ f_sl_g, 
                     int* __restrict__ nobuk_g,
                     double z0m, double z0h, double zsl,
                     int icells, int jcells, int kstart, int jj, int kk, 
                     Boundary::Boundary_type mbcbot, int thermobc)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y; 

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;

            // case 1: fixed buoyancy flux and fixed ustar
            if (mbcbot == Boundary::Ustar_type && thermobc == Boundary::Flux_type)
            {
                obuk[ij] = -pow(ustar[ij], 3) / (Constants::kappa*bfluxbot[ij]);
            }
            // case 2: fixed buoyancy flux and free ustar
            else if (mbcbot == Boundary::Dirichlet_type && thermobc == Boundary::Flux_type)
            {
                obuk [ij] = calc_Obuk_noslip_flux_g(zL_sl_g, f_sl_g, nobuk_g[ij], dutot[ij], bfluxbot[ij], zsl);
                ustar[ij] = dutot[ij] * most::fm(zsl, z0m, obuk[ij]);
            }
            // case 3: fixed buoyancy surface value and free ustar
            else if (mbcbot == Boundary::Dirichlet_type && thermobc == Boundary::Dirichlet_type)
            {
                double db = b[ijk] - bbot[ij];
                obuk [ij] = calc_Obuk_noslip_dirichlet_g(zL_sl_g, f_sl_g, nobuk_g[ij], dutot[ij], db, zsl);
                ustar[ij] = dutot[ij] * most::fm(zsl, z0m, obuk[ij]);
            }
        }
    }

    template<typename TF>
    __global__ 
    void stability_neutral_g(double* __restrict__ ustar, double* __restrict__ obuk,
                             double* __restrict__ dutot, double z0m, double z0h, double zsl,
                             int icells, int jcells, int kstart, int jj, int kk,
                             Boundary::Boundary_type mbcbot, int thermobc)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y; 

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;

            // case 1: fixed buoyancy flux and fixed ustar
            if (mbcbot == Boundary::Ustar_type && thermobc == Boundary::Flux_type)
            {
                obuk[ij] = -Constants::dbig;
            }
            // case 2: fixed buoyancy flux and free ustar
            else if (mbcbot == Boundary::Dirichlet_type && thermobc == Boundary::Flux_type)
            {
                obuk [ij] = -Constants::dbig;
                ustar[ij] = dutot[ij] * most::fm(zsl, z0m, obuk[ij]);
            }
            // case 3: fixed buoyancy surface value and free ustar
            else if (mbcbot == Boundary::Dirichlet_type && thermobc == Boundary::Dirichlet_type)
            {
                obuk [ij] = -Constants::dbig;
                ustar[ij] = dutot[ij] * most::fm(zsl, z0m, obuk[ij]);
            }
        }
    }

    template<typename TF>
    __global__ 
    void surfm_flux_g(double* __restrict__ ufluxbot, double* __restrict__ vfluxbot,
                      double* __restrict__ u,        double* __restrict__ v,
                      double* __restrict__ ubot,     double* __restrict__ vbot, 
                      double* __restrict__ ustar,    double* __restrict__ obuk, 
                      double zsl, double z0m,
                      int istart, int jstart, int kstart,
                      int iend,   int jend, int jj, int kk,
                      Boundary::Boundary_type bcbot)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 

        if (i < iend && j < jend)
        {
            const int ii  = 1;
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;

            if (bcbot == Boundary::Dirichlet_type)
            {
                // interpolate the whole stability function rather than ustar or obuk
                ufluxbot[ij] = -(u[ijk]-ubot[ij])*0.5*(ustar[ij-ii]*most::fm(zsl, z0m, obuk[ij-ii]) + ustar[ij]*most::fm(zsl, z0m, obuk[ij]));
                vfluxbot[ij] = -(v[ijk]-vbot[ij])*0.5*(ustar[ij-jj]*most::fm(zsl, z0m, obuk[ij-jj]) + ustar[ij]*most::fm(zsl, z0m, obuk[ij]));
            }
            else if (bcbot == Boundary::Ustar_type)
            {
                const double minval = 1.e-2;

                // minimize the wind at 0.01, thus the wind speed squared at 0.0001
                const double vonu2 = fmax(minval, 0.25*( pow(v[ijk-ii]-vbot[ij-ii], 2) + pow(v[ijk-ii+jj]-vbot[ij-ii+jj], 2)
                                                       + pow(v[ijk   ]-vbot[ij   ], 2) + pow(v[ijk   +jj]-vbot[ij   +jj], 2)) );
                const double uonv2 = fmax(minval, 0.25*( pow(u[ijk-jj]-ubot[ij-jj], 2) + pow(u[ijk+ii-jj]-ubot[ij+ii-jj], 2)
                                                       + pow(u[ijk   ]-ubot[ij   ], 2) + pow(u[ijk+ii   ]-ubot[ij+ii   ], 2)) );

                const double u2 = fmax(minval, pow(u[ijk]-ubot[ij], 2));
                const double v2 = fmax(minval, pow(v[ijk]-vbot[ij], 2));

                const double ustaronu4 = 0.5*(pow(ustar[ij-ii], 4) + pow(ustar[ij], 4));
                const double ustaronv4 = 0.5*(pow(ustar[ij-jj], 4) + pow(ustar[ij], 4));

                ufluxbot[ij] = -copysign(1., u[ijk]-ubot[ij]) * pow(ustaronu4 / (1. + vonu2 / u2), 0.5);
                vfluxbot[ij] = -copysign(1., v[ijk]-vbot[ij]) * pow(ustaronv4 / (1. + uonv2 / v2), 0.5);
            }
        }
    }

    template<typename TF>
    __global__ 
    void surfm_grad_g(double* __restrict__ ugradbot, double* __restrict__ vgradbot,
                      double* __restrict__ u,        double* __restrict__ v, 
                      double* __restrict__ ubot,     double* __restrict__ vbot, double zsl, 
                      int icells, int jcells, int kstart, int jj, int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y; 

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;

            ugradbot[ij] = (u[ijk]-ubot[ij])/zsl;
            vgradbot[ij] = (v[ijk]-vbot[ij])/zsl;
        }
    }

    template<typename TF>
    __global__ 
    void surfs_g(double* __restrict__ varfluxbot, double* __restrict__ vargradbot, 
                 double* __restrict__ varbot,     double* __restrict__ var, 
                 double* __restrict__ ustar,      double* __restrict__ obuk, double zsl, double z0h,
                 int icells, int jcells, int kstart,
                 int jj, int kk,
                 Boundary::Boundary_type bcbot)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y; 

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;

            if (bcbot == Boundary::Dirichlet_type)
            {
                varfluxbot[ij] = -(var[ijk]-varbot[ij])*ustar[ij]*most::fh(zsl, z0h, obuk[ij]);
                vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
            }
            else if (bcbot == Boundary::Flux_type)
            {
                varbot[ij]     = varfluxbot[ij] / (ustar[ij]*most::fh(zsl, z0h, obuk[ij])) + var[ijk];
                vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
            }
        }
    }
}

template<typename TF>
void Boundary_surface::prepare_device()
{
    const int dmemsize2d  = (grid->ijcellsp+grid->memoffset)*sizeof(double);
    const int imemsize2d  = (grid->ijcellsp+grid->memoffset)*sizeof(int);
    const int dimemsizep  = grid->icellsp * sizeof(double);
    const int dimemsize   = grid->icells  * sizeof(double);
    const int iimemsizep  = grid->icellsp * sizeof(int);
    const int iimemsize   = grid->icells  * sizeof(int);

    cuda_safe_call(cudaMalloc(&obuk_g,  dmemsize2d));
    cuda_safe_call(cudaMalloc(&ustar_g, dmemsize2d));
    cuda_safe_call(cudaMalloc(&nobuk_g, imemsize2d));

    cuda_safe_call(cudaMalloc(&zL_sl_g, nzL*sizeof(float)));
    cuda_safe_call(cudaMalloc(&f_sl_g,  nzL*sizeof(float)));

    cuda_safe_call(cudaMemcpy2D(&obuk_g[grid->memoffset],  dimemsizep, obuk,  dimemsize, dimemsize, grid->jcells, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy2D(&ustar_g[grid->memoffset], dimemsizep, ustar, dimemsize, dimemsize, grid->jcells, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy2D(&nobuk_g[grid->memoffset], iimemsizep, nobuk, iimemsize, iimemsize, grid->jcells, cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy(zL_sl_g, zL_sl, nzL*sizeof(float), cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(f_sl_g,  f_sl,  nzL*sizeof(float), cudaMemcpyHostToDevice));
}

// TMP BVS
template<typename TF>
void Boundary_surface::forward_device()
{
    const int dimemsizep  = grid->icellsp * sizeof(double);
    const int dimemsize   = grid->icells  * sizeof(double);
    const int iimemsizep  = grid->icellsp * sizeof(int);
    const int iimemsize   = grid->icells  * sizeof(int);

    cuda_safe_call(cudaMemcpy2D(&obuk_g[grid->memoffset],  dimemsizep, obuk,  dimemsize, dimemsize, grid->jcells,  cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy2D(&ustar_g[grid->memoffset], dimemsizep, ustar, dimemsize, dimemsize, grid->jcells,  cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy2D(&nobuk_g[grid->memoffset], iimemsizep, nobuk, iimemsize, iimemsize, grid->jcells,  cudaMemcpyHostToDevice));
}

// TMP BVS
template<typename TF>
void Boundary_surface::backward_device()
{
    const int dimemsizep  = grid->icellsp * sizeof(double);
    const int dimemsize   = grid->icells  * sizeof(double);
    const int iimemsizep  = grid->icellsp * sizeof(int);
    const int iimemsize   = grid->icells  * sizeof(int);

    cuda_safe_call(cudaMemcpy2D(obuk,  dimemsize, &obuk_g[grid->memoffset],  dimemsizep, dimemsize, grid->jcells,  cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy2D(ustar, dimemsize, &ustar_g[grid->memoffset], dimemsizep, dimemsize, grid->jcells,  cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy2D(nobuk, iimemsize, &nobuk_g[grid->memoffset], iimemsizep, iimemsize, grid->jcells,  cudaMemcpyDeviceToHost));
}

template<typename TF>
void Boundary_surface::clear_device()
{
    cuda_safe_call(cudaFree(obuk_g ));
    cuda_safe_call(cudaFree(ustar_g));
    cuda_safe_call(cudaFree(nobuk_g));
    cuda_safe_call(cudaFree(zL_sl_g));
    cuda_safe_call(cudaFree(f_sl_g ));
}

#ifdef USECUDA
template<typename TF>
void Boundary_surface::update_bcs()
{
    const int blocki = grid->ithread_block;
    const int blockj = grid->jthread_block;

    // For 2D field excluding ghost cells
    int gridi = grid->imax/blocki + (grid->imax%blocki > 0);
    int gridj = grid->jmax/blockj + (grid->jmax%blockj > 0);
    dim3 gridGPU (gridi,  gridj,  1);
    dim3 blockGPU(blocki, blockj, 1);

    // For 2D field including ghost cells
    gridi = grid->icells/blocki + (grid->icells%blocki > 0);
    gridj = grid->jcells/blockj + (grid->jcells%blockj > 0);
    dim3 gridGPU2 (gridi,  gridj,  1);
    dim3 blockGPU2(blocki, blockj, 1);

    const int offs = grid->memoffset;

    // Calculate dutot in tmp2
    du_tot_g<<<gridGPU, blockGPU>>>(
        &fields->atmp["tmp2"]->data_g[offs], 
        &fields->u->data_g[offs],    &fields->v->data_g[offs],
        &fields->u->databot_g[offs], &fields->v->databot_g[offs],
        grid->istart, grid->jstart, grid->kstart,
        grid->iend,   grid->jend,   grid->icellsp, grid->ijcellsp);
    cuda_check_error();

    // 2D cyclic boundaries on dutot  
    grid->boundary_cyclic2d_g(&fields->atmp["tmp2"]->data_g[offs]);

    // start with retrieving the stability information
    if (model->thermo->get_switch() == "0")
    {
        // Calculate ustar and Obukhov length, including ghost cells
        stability_neutral_g<<<gridGPU2, blockGPU2>>>(
            &ustar_g[offs], &obuk_g[offs], 
            &fields->atmp["tmp2"]->data_g[offs], z0m, z0h, grid->z[grid->kstart],
            grid->icells, grid->jcells, grid->kstart, grid->icellsp, grid->ijcellsp, mbcbot, thermobc); 
        cuda_check_error();
    }
    else
    {
        // store the buoyancy in tmp1
        model->thermo->get_buoyancy_surf(fields->atmp["tmp1"]);

        // Calculate ustar and Obukhov length, including ghost cells
        stability_g<<<gridGPU2, blockGPU2>>>(
            &ustar_g[offs], &obuk_g[offs], 
            &fields->atmp["tmp1"]->data_g[offs], &fields->atmp["tmp1"]->databot_g[offs], &fields->atmp["tmp1"]->datafluxbot_g[offs],
            &fields->atmp["tmp2"]->data_g[offs], 
            zL_sl_g, f_sl_g, &nobuk_g[offs],
            z0m, z0h, grid->z[grid->kstart],
            grid->icells, grid->jcells, grid->kstart, grid->icellsp, grid->ijcellsp, mbcbot, thermobc); 
        cuda_check_error();
    }

    // Calculate surface momentum fluxes, excluding ghost cells
    surfm_flux_g<<<gridGPU, blockGPU>>>(
        &fields->u->datafluxbot_g[offs], &fields->v->datafluxbot_g[offs],
        &fields->u->data_g[offs],        &fields->v->data_g[offs],
        &fields->u->databot_g[offs],     &fields->v->databot_g[offs], 
        &ustar_g[offs], &obuk_g[offs], grid->z[grid->kstart], z0m,
        grid->istart, grid->jstart, grid->kstart,
        grid->iend,   grid->jend,   grid->icellsp, grid->ijcellsp, mbcbot);
    cuda_check_error();

    // 2D cyclic boundaries on the surface fluxes  
    grid->boundary_cyclic2d_g(&fields->u->datafluxbot_g[offs]);
    grid->boundary_cyclic2d_g(&fields->v->datafluxbot_g[offs]);

    // Calculate surface gradients, including ghost cells
    surfm_grad_g<<<gridGPU2, blockGPU2>>>(
        &fields->u->datagradbot_g[offs], &fields->v->datagradbot_g[offs],
        &fields->u->data_g[offs],        &fields->v->data_g[offs],
        &fields->u->databot_g[offs],     &fields->v->databot_g[offs],
        grid->z[grid->kstart], grid->icells, grid->jcells, grid->kstart, grid->icellsp, grid->ijcellsp);  
    cuda_check_error();

    // Calculate scalar fluxes, gradients and/or values, including ghost cells
    for (FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
        surfs_g<<<gridGPU2, blockGPU2>>>(
            &it->second->datafluxbot_g[offs], &it->second->datagradbot_g[offs],
            &it->second->databot_g[offs],     &it->second->data_g[offs],
            &ustar_g[offs], &obuk_g[offs], grid->z[grid->kstart], z0h,            
            grid->icells,  grid->jcells, grid->kstart,
            grid->icellsp, grid->ijcellsp, sbc[it->first]->bcbot);
    cuda_check_error();
}
#endif
