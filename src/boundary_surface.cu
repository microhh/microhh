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

//texture<float, 1, cudaReadModeElementType> zl_tex;
//texture<float, 1, cudaReadModeElementType> f_tex;

namespace BoundarySurface_g
{
  // a sign function
  __device__ double sign(double n) 
  { 
    //return n > 0 ? 1 : (n < 0 ? -1 : 0); // Doesn't work?
    if(n > 0)
      return 1.;
    else if(n < 0)
      return -1.;
    else
      return 0.;
  }
  
  __device__ double psim(double zeta)
  {
    double psim;
    double x;
    if(zeta <= 0.)
    {
      x    = pow(1. + pow(3.6 * fabs(zeta),2./3.), -0.5);
      psim = 3.*log( (1. + 1./x) / 2.);
    }
    else
    {
      psim = -2./3.*(zeta - 5./0.35) * exp(-0.35 * zeta) - zeta - (10./3.) / 0.35;
    }
    return psim;
  }
  
  __device__ double psih(double zeta)
  {
    double psih;
    double x;
    if(zeta <= 0.)
    {
      x    = pow(1. + pow(7.9*fabs(zeta), (2./3.)), -0.5);
      psih = 3. * log( (1. + 1. / x) / 2.);
    }
    else
    {
      psih  = (-2./3.) * (zeta-5./0.35) * exp(-0.35*zeta) - pow(1. + (2./3.) * zeta, 1.5) - (10./3.) / 0.35 + 1.;
    }
    return psih;
  }
  
  __device__ double phim(double zeta)
  {
    double phim;
    if(zeta <= 0.)
    {
      phim = pow(1. + 3.6*pow(fabs(zeta), 2./3.), -1./2.);
    }
    else
      phim = 1. + 5.*zeta;
  
    return phim;
  }
  
  __device__ double phih(double zeta)
  {
    double phih;
    if(zeta <= 0.)
    {
      phih = pow(1. + 7.9*pow(fabs(zeta), 2./3.), -1./2.);
    }
    else
      phih = 1. + 5.*zeta;
  
    return phih;
  }
  
  __device__ double fm(double zsl, double z0m, double L) 
  { 
    return constants::kappa / (log(zsl/z0m) - psim(zsl/L) + psim(z0m/L)); 
  }
  
  __device__ double fh(double zsl, double z0h, double L) 
  { 
    return constants::kappa / (log(zsl/z0h) - psih(zsl/L) + psih(z0h/L)); 
  }

  __device__ double findObuk(const float* const __restrict__ zL, const float* const __restrict__ f,
                  int &n, const double Ri, const double zsl)
  {
    // Determine search direction.
    if ( (f[n]-Ri) > 0 )
      while ( (f[n-1]-Ri) > 0 && n > 0) { --n; }
    else
      while ( (f[n]-Ri) < 0 && n < (nzL-1) ) { ++n; }

    const double zL0 = (n == 0 || n == nzL-1) ? zL[n] : zL[n-1] + (Ri-f[n-1]) / (f[n]-f[n-1]) * (zL[n]-zL[n-1]);

    return zsl/zL0;
  }
 
  __device__ double calcobuk_noslip_flux(float * __restrict__ zL, float * __restrict__ f, int& n, double du, double bfluxbot, double zsl)
  {
    // Calculate the appropriate Richardson number.
    const double Ri = -constants::kappa * bfluxbot * zsl / pow(du, 3);
    return findObuk(zL, f, n, Ri, zsl);
  }
  
  __device__ double calcobuk_noslip_dirichlet(float * __restrict__ zL, float * __restrict__ f, int& n, double du, double db, double zsl)
  {
    // Calculate the appropriate Richardson number.
    const double Ri = constants::kappa * db * zsl / pow(du, 2);
    return findObuk(zL, f, n, Ri, zsl);
  }
  
  /* Calculate absolute wind speed */
  __global__ void dutot(double * __restrict__ dutot, 
                        double * __restrict__ u,    double * __restrict__ v,
                        double * __restrict__ ubot, double * __restrict__ vbot, 
                        int istart, int jstart, int kstart,
                        int iend,   int jend, int jj, int kk)
  {
    int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
    int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
  
    if(i < iend && j < jend)
    {
      int ii  = 1;
      int ij  = i + j*jj;
      int ijk = i + j*jj + kstart*kk;
      const double minval = 1.e-1;
  
      double du2 = pow(0.5*(u[ijk] + u[ijk+ii]) - 0.5*(ubot[ij] + ubot[ij+ii]), 2)
                 + pow(0.5*(v[ijk] + v[ijk+jj]) - 0.5*(vbot[ij] + vbot[ij+jj]), 2);
      dutot[ij] = fmax(pow(du2, 0.5), minval);
    }
  }
  
  //template <int mbcbot, int thermobc>  // BvS for now normal parameter. Make template again...
  __global__ void stability(double * __restrict__ ustar, double * __restrict__ obuk,
                            double * __restrict__ b, double * __restrict__ bbot, double * __restrict__ bfluxbot,
                            double * __restrict__ dutot, float * __restrict__ zL_sl_g, float * __restrict__ f_sl_g, 
                            int * __restrict__ nobuk_g,
                            double z0m, double z0h, double zsl,
                            int icells, int jcells, int kstart, int jj, int kk, 
                            Boundary::BoundaryType mbcbot, int thermobc)
  {
    int i = blockIdx.x*blockDim.x + threadIdx.x; 
    int j = blockIdx.y*blockDim.y + threadIdx.y; 
  
    if(i < icells && j < jcells)
    {
      int ij  = i + j*jj;
      int ijk = i + j*jj + kstart*kk;
  
      // case 1: fixed buoyancy flux and fixed ustar
      if(mbcbot == Boundary::UstarType && thermobc == Boundary::FluxType)
      {
        obuk[ij] = -pow(ustar[ij], 3) / (constants::kappa*bfluxbot[ij]);
      }
      // case 2: fixed buoyancy flux and free ustar
      else if(mbcbot == Boundary::DirichletType && thermobc == Boundary::FluxType)
      {
        obuk [ij] = calcobuk_noslip_flux(zL_sl_g, f_sl_g, nobuk_g[ij], dutot[ij], bfluxbot[ij], zsl);
        ustar[ij] = dutot[ij] * fm(zsl, z0m, obuk[ij]);
      }
      // case 3: fixed buoyancy surface value and free ustar
      else if(mbcbot == Boundary::DirichletType && thermobc == Boundary::DirichletType)
      {
        double db = b[ijk] - bbot[ij];
        obuk [ij] = calcobuk_noslip_dirichlet(zL_sl_g, f_sl_g, nobuk_g[ij], dutot[ij], db, zsl);
        ustar[ij] = dutot[ij] * fm(zsl, z0m, obuk[ij]);
      }
    }
  }
  
  //template <int mbcbot, int thermobc>  // BvS for now normal parameter. Make template again...
  __global__ void stability_neutral(double * __restrict__ ustar, double * __restrict__ obuk,
                                    double * __restrict__ dutot, double z0m, double z0h, double zsl,
                                    int icells, int jcells, int kstart, int jj, int kk,
                                    Boundary::BoundaryType mbcbot, int thermobc)
  {
    int i = blockIdx.x*blockDim.x + threadIdx.x; 
    int j = blockIdx.y*blockDim.y + threadIdx.y; 
  
    if(i < icells && j < jcells)
    {
      int ij  = i + j*jj;
  
      // case 1: fixed buoyancy flux and fixed ustar
      if(mbcbot == Boundary::UstarType && thermobc == Boundary::FluxType)
      {
        obuk[ij] = -constants::dbig;
      }
      // case 2: fixed buoyancy flux and free ustar
      else if(mbcbot == Boundary::DirichletType && thermobc == Boundary::FluxType)
      {
        obuk [ij] = -constants::dbig;
        ustar[ij] = dutot[ij] * fm(zsl, z0m, obuk[ij]);
      }
      // case 3: fixed buoyancy surface value and free ustar
      else if(mbcbot == Boundary::DirichletType && thermobc == Boundary::DirichletType)
      {
        obuk [ij] = -constants::dbig;
        ustar[ij] = dutot[ij] * fm(zsl, z0m, obuk[ij]);
      }
    }
  }
  
  __global__ void surfm_flux(double * __restrict__ ufluxbot, double * __restrict__ vfluxbot,
                             double * __restrict__ u,        double * __restrict__ v,
                             double * __restrict__ ubot,     double * __restrict__ vbot, 
                             double * __restrict__ ustar,    double * __restrict__ obuk, 
                             double zsl, double z0m,
                             int istart, int jstart, int kstart,
                             int iend,   int jend, int jj, int kk,
                             Boundary::BoundaryType bcbot)
  {
    int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
    int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
  
    if(i < iend && j < jend)
    {
      int ii  = 1;
      int ij  = i + j*jj;
      int ijk = i + j*jj + kstart*kk;
  
      if(bcbot == Boundary::DirichletType)
      {
        // interpolate the whole stability function rather than ustar or obuk
        ufluxbot[ij] = -(u[ijk]-ubot[ij])*0.5*(ustar[ij-ii]*fm(zsl, z0m, obuk[ij-ii]) + ustar[ij]*fm(zsl, z0m, obuk[ij]));
        vfluxbot[ij] = -(v[ijk]-vbot[ij])*0.5*(ustar[ij-jj]*fm(zsl, z0m, obuk[ij-jj]) + ustar[ij]*fm(zsl, z0m, obuk[ij]));
      }
      else if(bcbot == Boundary::UstarType)
      {
        double u2,v2,vonu2,uonv2,ustaronu4,ustaronv4;
        const double minval = 1.e-2;
  
        // minimize the wind at 0.01, thus the wind speed squared at 0.0001
        vonu2 = fmax(minval, 0.25*( pow(v[ijk-ii]-vbot[ij-ii], 2) + pow(v[ijk-ii+jj]-vbot[ij-ii+jj], 2)
                                  + pow(v[ijk   ]-vbot[ij   ], 2) + pow(v[ijk   +jj]-vbot[ij   +jj], 2)) );
        uonv2 = fmax(minval, 0.25*( pow(u[ijk-jj]-ubot[ij-jj], 2) + pow(u[ijk+ii-jj]-ubot[ij+ii-jj], 2)
                                  + pow(u[ijk   ]-ubot[ij   ], 2) + pow(u[ijk+ii   ]-ubot[ij+ii   ], 2)) );
  
        u2 = fmax(minval, pow(u[ijk]-ubot[ij], 2));
        v2 = fmax(minval, pow(v[ijk]-vbot[ij], 2));
  
        ustaronu4 = 0.5*(pow(ustar[ij-ii], 4) + pow(ustar[ij], 4));
        ustaronv4 = 0.5*(pow(ustar[ij-jj], 4) + pow(ustar[ij], 4));
  
        ufluxbot[ij] = -sign(u[ijk]-ubot[ij]) * pow(ustaronu4 / (1. + vonu2 / u2), 0.5);
        vfluxbot[ij] = -sign(v[ijk]-vbot[ij]) * pow(ustaronv4 / (1. + uonv2 / v2), 0.5);
      }
    }
  }
  
  __global__ void surfm_grad(double * __restrict__ ugradbot, double * __restrict__ vgradbot,
                             double * __restrict__ u,        double * __restrict__ v, 
                             double * __restrict__ ubot,     double * __restrict__ vbot, double zsl, 
                             int icells, int jcells, int kstart, int jj, int kk)
  {
    int i = blockIdx.x*blockDim.x + threadIdx.x; 
    int j = blockIdx.y*blockDim.y + threadIdx.y; 
  
    if(i < icells && j < jcells)
    {
      int ij  = i + j*jj;
      int ijk = i + j*jj + kstart*kk;
  
      ugradbot[ij] = (u[ijk]-ubot[ij])/zsl;
      vgradbot[ij] = (v[ijk]-vbot[ij])/zsl;
    }
  }
  
  __global__ void surfs(double * __restrict__ varfluxbot, double * __restrict__ vargradbot, 
                        double * __restrict__ varbot, double * __restrict__ var, 
                        double * __restrict__ ustar, double * __restrict__ obuk, double zsl, double z0h,
                        int icells, int jcells, int kstart,
                        int jj, int kk,
                        Boundary::BoundaryType bcbot)
  {
    int i = blockIdx.x*blockDim.x + threadIdx.x; 
    int j = blockIdx.y*blockDim.y + threadIdx.y; 
  
    if(i < icells && j < jcells)
    {
      int ij  = i + j*jj;
      int ijk = i + j*jj + kstart*kk;
  
      if(bcbot == Boundary::DirichletType)
      {
        varfluxbot[ij] = -(var[ijk]-varbot[ij])*ustar[ij]*fh(zsl, z0h, obuk[ij]);
        vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
      }
      else if(bcbot == Boundary::FluxType)
      {
        varbot[ij]     = varfluxbot[ij] / (ustar[ij]*fh(zsl, z0h, obuk[ij])) + var[ijk];
        vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
      }
    }
  }
}

void BoundarySurface::prepareDevice()
{
  const int dmemsize2d  = (grid->ijcellsp+grid->memoffset)*sizeof(double);
  const int imemsize2d  = (grid->ijcellsp+grid->memoffset)*sizeof(int);
  const int dimemsizep  = grid->icellsp * sizeof(double);
  const int dimemsize   = grid->icells  * sizeof(double);
  const int iimemsizep  = grid->icellsp * sizeof(int);
  const int iimemsize   = grid->icells  * sizeof(int);

  cudaSafeCall(cudaMalloc(&obuk_g,  dmemsize2d));
  cudaSafeCall(cudaMalloc(&ustar_g, dmemsize2d));
  cudaSafeCall(cudaMalloc(&nobuk_g, imemsize2d));

  cudaSafeCall(cudaMalloc(&zL_sl_g, nzL*sizeof(float)));
  cudaSafeCall(cudaMalloc(&f_sl_g,  nzL*sizeof(float)));

  cudaSafeCall(cudaMemcpy2D(&obuk_g[grid->memoffset],  dimemsizep, obuk,  dimemsize, dimemsize, grid->jcells, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy2D(&ustar_g[grid->memoffset], dimemsizep, ustar, dimemsize, dimemsize, grid->jcells, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy2D(&nobuk_g[grid->memoffset], iimemsizep, nobuk, iimemsize, iimemsize, grid->jcells, cudaMemcpyHostToDevice));

  cudaSafeCall(cudaMemcpy(zL_sl_g, zL_sl, nzL*sizeof(float), cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(f_sl_g,  f_sl,  nzL*sizeof(float), cudaMemcpyHostToDevice));

  //cudaBindTexture(0, zL_tex, zL_sl_g, nzL*sizeof(float)); 
  //cudaBindTexture(0, f_tex,  f_sl_g,  nzL*sizeof(float)); 
}

// TMP BVS
void BoundarySurface::forwardDevice()
{
  const int dimemsizep  = grid->icellsp * sizeof(double);
  const int dimemsize   = grid->icells  * sizeof(double);
  const int iimemsizep  = grid->icellsp * sizeof(int);
  const int iimemsize   = grid->icells  * sizeof(int);

  cudaSafeCall(cudaMemcpy2D(&obuk_g[grid->memoffset],  dimemsizep, obuk,  dimemsize, dimemsize, grid->jcells,  cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy2D(&ustar_g[grid->memoffset], dimemsizep, ustar, dimemsize, dimemsize, grid->jcells,  cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy2D(&nobuk_g[grid->memoffset], iimemsizep, nobuk, iimemsize, iimemsize, grid->jcells,  cudaMemcpyHostToDevice));
}

// TMP BVS
void BoundarySurface::backwardDevice()
{
  const int dimemsizep  = grid->icellsp * sizeof(double);
  const int dimemsize   = grid->icells  * sizeof(double);
  const int iimemsizep  = grid->icellsp * sizeof(int);
  const int iimemsize   = grid->icells  * sizeof(int);

  cudaSafeCall(cudaMemcpy2D(obuk,  dimemsize, &obuk_g[grid->memoffset],  dimemsizep, dimemsize, grid->jcells,  cudaMemcpyDeviceToHost));
  cudaSafeCall(cudaMemcpy2D(ustar, dimemsize, &ustar_g[grid->memoffset], dimemsizep, dimemsize, grid->jcells,  cudaMemcpyDeviceToHost));
  cudaSafeCall(cudaMemcpy2D(nobuk, iimemsize, &nobuk_g[grid->memoffset], iimemsizep, iimemsize, grid->jcells,  cudaMemcpyDeviceToHost));
}

void BoundarySurface::clearDevice()
{
  cudaSafeCall(cudaFree(obuk_g ));
  cudaSafeCall(cudaFree(ustar_g));
  cudaSafeCall(cudaFree(nobuk_g));
  cudaSafeCall(cudaFree(zL_sl_g));
  cudaSafeCall(cudaFree(f_sl_g ));
}

#ifdef USECUDA
void BoundarySurface::updateBcs()
{
  int gridi, gridj;
  const int blocki = grid->iThreadBlock;
  const int blockj = grid->jThreadBlock;

  // For 2D field excluding ghost cells
  gridi = grid->imax/blocki + (grid->imax%blocki > 0);
  gridj = grid->jmax/blockj + (grid->jmax%blockj > 0);
  dim3 gridGPU (gridi,  gridj,  1);
  dim3 blockGPU(blocki, blockj, 1);
 
  // For 2D field including ghost cells
  gridi = grid->icells/blocki + (grid->icells%blocki > 0);
  gridj = grid->jcells/blockj + (grid->jcells%blockj > 0);
  dim3 gridGPU2 (gridi,  gridj,  1);
  dim3 blockGPU2(blocki, blockj, 1);
 
  const int offs = grid->memoffset;

  // Calculate dutot in tmp2
  BoundarySurface_g::dutot<<<gridGPU, blockGPU>>>(&fields->atmp["tmp2"]->data_g[offs], 
                                                  &fields->u->data_g[offs],    &fields->v->data_g[offs],
                                                  &fields->u->databot_g[offs], &fields->v->databot_g[offs],
                                                  grid->istart, grid->jstart, grid->kstart,
                                                  grid->iend,   grid->jend,   grid->icellsp, grid->ijcellsp);
  cudaCheckError();

  // 2D cyclic boundaries on dutot  
  grid->boundaryCyclic2d_g(&fields->atmp["tmp2"]->data_g[offs]);

  // start with retrieving the stability information
  if(model->thermo->getSwitch() == "0")
  {
    // Calculate ustar and Obukhov length, including ghost cells
    BoundarySurface_g::stability_neutral<<<gridGPU2, blockGPU2>>>(&ustar_g[offs], &obuk_g[offs], 
                                                  &fields->atmp["tmp2"]->data_g[offs], z0m, z0h, grid->z[grid->kstart],
                                                  grid->icells, grid->jcells, grid->kstart, grid->icellsp, grid->ijcellsp, mbcbot, thermobc); 
    cudaCheckError();
  }
  else
  {
    // store the buoyancy in tmp1
    model->thermo->getBuoyancySurf(fields->atmp["tmp1"]);

    // Calculate ustar and Obukhov length, including ghost cells
    BoundarySurface_g::stability<<<gridGPU2, blockGPU2>>>(&ustar_g[offs], &obuk_g[offs], 
                                                  &fields->atmp["tmp1"]->data_g[offs], &fields->atmp["tmp1"]->databot_g[offs], &fields->atmp["tmp1"]->datafluxbot_g[offs],
                                                  &fields->atmp["tmp2"]->data_g[offs], 
                                                  zL_sl_g, f_sl_g, &nobuk_g[offs],
                                                  z0m, z0h, grid->z[grid->kstart],
                                                  grid->icells, grid->jcells, grid->kstart, grid->icellsp, grid->ijcellsp, mbcbot, thermobc); 

    cudaCheckError();
  }

  // Calculate surface momentum fluxes, excluding ghost cells
  BoundarySurface_g::surfm_flux<<<gridGPU, blockGPU>>>(&fields->u->datafluxbot_g[offs], &fields->v->datafluxbot_g[offs],
                                                       &fields->u->data_g[offs],        &fields->v->data_g[offs],
                                                       &fields->u->databot_g[offs],     &fields->v->databot_g[offs], 
                                                       &ustar_g[offs], &obuk_g[offs], grid->z[grid->kstart], z0m,
                                                       grid->istart, grid->jstart, grid->kstart,
                                                       grid->iend,   grid->jend,   grid->icellsp, grid->ijcellsp, mbcbot);
  cudaCheckError();

  // 2D cyclic boundaries on the surface fluxes  
  grid->boundaryCyclic2d_g(&fields->u->datafluxbot_g[offs]);
  grid->boundaryCyclic2d_g(&fields->v->datafluxbot_g[offs]);

  // Calculate surface gradients, including ghost cells
  BoundarySurface_g::surfm_grad<<<gridGPU2, blockGPU2>>>(&fields->u->datagradbot_g[offs], &fields->v->datagradbot_g[offs],
                                                         &fields->u->data_g[offs],        &fields->v->data_g[offs],
                                                         &fields->u->databot_g[offs],     &fields->v->databot_g[offs],
                                                         grid->z[grid->kstart], grid->icells, grid->jcells, grid->kstart, grid->icellsp, grid->ijcellsp);  
  cudaCheckError();

  // Calculate scalar fluxes, gradients and/or values, including ghost cells
  for(FieldMap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    BoundarySurface_g::surfs<<<gridGPU2, blockGPU2>>>(&it->second->datafluxbot_g[offs], &it->second->datagradbot_g[offs],
                                                      &it->second->databot_g[offs],     &it->second->data_g[offs],
                                                      &ustar_g[offs], &obuk_g[offs], grid->z[grid->kstart], z0h,            
                                                      grid->icells, grid->jcells, grid->kstart,
                                                      grid->icellsp, grid->ijcellsp, sbc[it->first]->bcbot);
  cudaCheckError();
}
#endif
