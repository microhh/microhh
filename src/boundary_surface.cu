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

#define NO_VELOCITY 0.
#define NO_OFFSET 0.

#define BC_DIRICHLET 0
#define BC_NEUMANN 1
#define BC_FLUX 2
#define BC_USTAR 3

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

__device__ double boundary_surface_psim(double zeta)
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

__device__ double boundary_surface_psih(double zeta)
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

__device__ double boundary_surface_phim(double zeta)
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

__device__ double boundary_surface_phih(double zeta)
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

__device__ double boundary_surface_fm(double zsl, double z0m, double L) 
{ 
  return constants::kappa / (log(zsl/z0m) - boundary_surface_psim(zsl/L) + boundary_surface_psim(z0m/L)); 
}

__device__ double boundary_surface_fh(double zsl, double z0h, double L) 
{ 
  return constants::kappa / (log(zsl/z0h) - boundary_surface_psih(zsl/L) + boundary_surface_psih(z0h/L)); 
}

__device__ double boundary_surface_calcobuk_noslip_flux(double L, double du, double bfluxbot, double zsl, double z0m)
{
  double L0;
  double Lstart, Lend;
  double fx, fxdif;

  int m = 0;
  int nlim = 10;

  const double Lmax = 1.e20;

  // avoid bfluxbot to be zero
  if(bfluxbot >= 0.)
    bfluxbot = fmax(constants::dsmall, bfluxbot);
  else
    bfluxbot = fmin(-constants::dsmall, bfluxbot);

  // allow for one restart
  while(m <= 1)
  {
    // if L and bfluxbot are of the same sign, or the last calculation did not converge,
    // the stability has changed and the procedure needs to be reset
    if(L*bfluxbot >= 0.)
    {
      nlim = 200;
      if(bfluxbot >= 0.)
        L = -constants::dsmall;
      else
        L = constants::dsmall;
    }

    if(bfluxbot >= 0.)
      L0 = -constants::dhuge;
    else
      L0 = constants::dhuge;

    int n = 0;

    // exit on convergence or on iteration count
    while(fabs((L - L0)/L0) > 0.001 && n < nlim && fabs(L) < Lmax)
    {
      L0     = L;
      // fx     = Rib - zsl/L * (std::log(zsl/z0h) - psih(zsl/L) + psih(z0h/L)) / std::pow(std::log(zsl/z0m) - psim(zsl/L) + psim(z0m/L), 2.);
      fx     = zsl/L + constants::kappa*zsl*bfluxbot / pow(du * boundary_surface_fm(zsl, z0m, L), 3);
      Lstart = L - 0.001*L;
      Lend   = L + 0.001*L;
      fxdif  = ( (zsl/Lend   + constants::kappa*zsl*bfluxbot / pow(du * boundary_surface_fm(zsl, z0m, Lend),   3))
               - (zsl/Lstart + constants::kappa*zsl*bfluxbot / pow(du * boundary_surface_fm(zsl, z0m, Lstart), 3)) )
             / (Lend - Lstart);
      L      = L - fx/fxdif;
      ++n;
    }

    // convergence has been reached
    if(n < nlim && fabs(L) < Lmax)
      break;
    // convergence has not been reached, procedure restarted once
    else
    {
      L = constants::dsmall;
      ++m;
      nlim = 200;
    }
  }

  if(m > 1)
    printf("ERROR convergence has not been reached in Obukhov length calculation\n");

  return L;
}

__device__ double boundary_surface_calcobuk_noslip_dirichlet(double L, double du, double db, double zsl, double z0m, double z0h)
{
  double L0;
  double Lstart, Lend;
  double fx, fxdif;

  int m = 0;
  int nlim = 10;

  const double Lmax = 1.e20;

  // avoid db to be zero
  if(db >= 0.)
    db = fmax(constants::dsmall, db);
  else
    db = fmin(-constants::dsmall, db);

  // allow for one restart
  while(m <= 1)
  {
    // if L and db are of different sign, or the last calculation did not converge,
    // the stability has changed and the procedure needs to be reset
    if(L*db <= 0.)
    {
      nlim = 200;
      if(db >= 0.)
        L = constants::dsmall;
      else
        L = -constants::dsmall;
    }

    if(db >= 0.)
      L0 = constants::dhuge;
    else
      L0 = -constants::dhuge;

    int n = 0;

    // exit on convergence or on iteration count
    while(fabs((L - L0)/L0) > 0.001 && n < nlim && fabs(L) < Lmax)
    {
      L0     = L;
      // fx     = Rib - zsl/L * (std::log(zsl/z0h) - psih(zsl/L) + psih(z0h/L)) / std::pow(std::log(zsl/z0m) - psim(zsl/L) + psim(z0m/L), 2.);
      fx     = zsl/L - constants::kappa*zsl*db*boundary_surface_fh(zsl, z0h, L) / pow(du * boundary_surface_fm(zsl, z0m, L), 2);
      Lstart = L - 0.001*L;
      Lend   = L + 0.001*L;
      fxdif  = ( (zsl/Lend   - constants::kappa*zsl*db*boundary_surface_fh(zsl, z0h, Lend)   / pow(du * boundary_surface_fm(zsl, z0m, Lend),   2))
               - (zsl/Lstart - constants::kappa*zsl*db*boundary_surface_fh(zsl, z0h, Lstart) / pow(du * boundary_surface_fm(zsl, z0m, Lstart), 2)) )
             / (Lend - Lstart);
      L      = L - fx/fxdif;
      ++n;
    }

    // convergence has been reached
    if(n < nlim && fabs(L) < Lmax)
      break;
    // convergence has not been reached, procedure restarted once
    else
    {
      L = constants::dsmall;
      ++m;
      nlim = 200;
    }
  }

  if(m > 1)
    printf("ERROR convergence has not been reached in Obukhov length calculation\n");

  return L;
}

/* Calculate absolute wind speed */
__global__ void boundary_surface_dutot(double * __restrict__ dutot, 
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
__global__ void boundary_surface_stability(double * __restrict__ ustar, double * __restrict__ obuk,
                                           double * __restrict__ b, double * __restrict__ bbot, double * __restrict__ bfluxbot,
                                           double * __restrict__ dutot, double z0m, double z0h, double zsl,
                                           int icells, int jcells, int kstart, int jj, int kk, int mbcbot, int thermobc)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x; 
  int j = blockIdx.y*blockDim.y + threadIdx.y; 

  if(i < icells && j < jcells)
  {
    int ij  = i + j*jj;
    int ijk = i + j*jj + kstart*kk;

    // case 1: fixed buoyancy flux and fixed ustar
    if(mbcbot == BC_USTAR && thermobc == BC_FLUX)
    {
      obuk[ij] = -pow(ustar[ij], 3) / (constants::kappa*bfluxbot[ij]);
    }
    // case 2: fixed buoyancy flux and free ustar
    else if(mbcbot == BC_DIRICHLET && thermobc == BC_FLUX)
    {
      obuk [ij] = boundary_surface_calcobuk_noslip_flux(obuk[ij], dutot[ij], bfluxbot[ij], zsl, z0m);
      ustar[ij] = dutot[ij] * boundary_surface_fm(zsl, z0m, obuk[ij]);
    }
    // case 3: fixed buoyancy surface value and free ustar
    else if(mbcbot == BC_DIRICHLET && thermobc == BC_DIRICHLET)
    {
      double db = b[ijk] - bbot[ij];
      obuk [ij] = boundary_surface_calcobuk_noslip_dirichlet(obuk[ij], dutot[ij], db, zsl, z0m, z0h);
      ustar[ij] = dutot[ij] * boundary_surface_fm(zsl, z0m, obuk[ij]);
    }
  }
}

__global__ void boundary_surface_surfm_flux(double * __restrict__ ufluxbot, double * __restrict__ vfluxbot,
                                            double * __restrict__ u,        double * __restrict__ v,
                                            double * __restrict__ ubot,     double * __restrict__ vbot, 
                                            double * __restrict__ ustar,    double * __restrict__ obuk, 
                                            double zsl, double z0m,
                                            int istart, int jstart, int kstart,
                                            int iend,   int jend, int jj, int kk, int bcbot)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
  int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 

  if(i < iend && j < jend)
  {
    int ii  = 1;
    int ij  = i + j*jj;
    int ijk = i + j*jj + kstart*kk;

    if(bcbot == BC_DIRICHLET)
    {
      // interpolate the whole stability function rather than ustar or obuk
      ufluxbot[ij] = -(u[ijk]-ubot[ij])*0.5*(ustar[ij-ii]*boundary_surface_fm(zsl, z0m, obuk[ij-ii]) + ustar[ij]*boundary_surface_fm(zsl, z0m, obuk[ij]));
      vfluxbot[ij] = -(v[ijk]-vbot[ij])*0.5*(ustar[ij-jj]*boundary_surface_fm(zsl, z0m, obuk[ij-jj]) + ustar[ij]*boundary_surface_fm(zsl, z0m, obuk[ij]));
    }
    else if(bcbot == BC_USTAR)
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

__global__ void boundary_surface_surfm_grad(double * __restrict__ ugradbot, double * __restrict__ vgradbot,
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

__global__ void boundary_surface_surfs(double * __restrict__ varfluxbot, double * __restrict__ vargradbot, 
                                       double * __restrict__ varbot, double * __restrict__ var, 
                                       double * __restrict__ ustar, double * __restrict__ obuk, double zsl, double z0h,
                                       int icells, int jcells, int kstart,
                                       int jj, int kk, int bcbot)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x; 
  int j = blockIdx.y*blockDim.y + threadIdx.y; 

  if(i < icells && j < jcells)
  {
    int ij  = i + j*jj;
    int ijk = i + j*jj + kstart*kk;

    if(bcbot == BC_DIRICHLET)
    {
      varfluxbot[ij] = -(var[ijk]-varbot[ij])*ustar[ij]*boundary_surface_fh(zsl, z0h, obuk[ij]);
      vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
    }
    else if(bcbot == BC_FLUX)
    {
      varbot[ij]     = varfluxbot[ij] / (ustar[ij]*boundary_surface_fh(zsl, z0h, obuk[ij])) + var[ijk];
      vargradbot[ij] = (var[ijk]-varbot[ij])/zsl;
    }
  }
}

int cboundary_surface::prepareDevice()
{
  const int nmemsize2d = (grid->ijcellsp+grid->memoffset)*sizeof(double);
  const int imemsizep  = grid->icellsp * sizeof(double);
  const int imemsize   = grid->icells  * sizeof(double);

  cudaMalloc(&obuk_g,  nmemsize2d);
  cudaMalloc(&ustar_g, nmemsize2d);

  cudaMemcpy2D(&obuk_g[grid->memoffset],  imemsizep, obuk, imemsize, imemsize, grid->jcells,  cudaMemcpyHostToDevice);
  cudaMemcpy2D(&ustar_g[grid->memoffset], imemsizep, ustar, imemsize, imemsize, grid->jcells,  cudaMemcpyHostToDevice);

  return 0;
}

// TMP BVS
int cboundary_surface::forwardDevice()
{
  const int imemsizep  = grid->icellsp * sizeof(double);
  const int imemsize   = grid->icells  * sizeof(double);

  cudaMemcpy2D(&obuk_g[grid->memoffset],  imemsizep, obuk,  imemsize, imemsize, grid->jcells,  cudaMemcpyHostToDevice);
  cudaMemcpy2D(&ustar_g[grid->memoffset], imemsizep, ustar, imemsize, imemsize, grid->jcells,  cudaMemcpyHostToDevice);

  return 0;
}

// TMP BVS
int cboundary_surface::backwardDevice()
{
  const int imemsizep  = grid->icellsp * sizeof(double);
  const int imemsize   = grid->icells  * sizeof(double);

  cudaMemcpy2D(obuk,  imemsize, &obuk_g[grid->memoffset],  imemsizep, imemsize, grid->jcells,  cudaMemcpyDeviceToHost);
  cudaMemcpy2D(ustar, imemsize, &ustar_g[grid->memoffset], imemsizep, imemsize, grid->jcells,  cudaMemcpyDeviceToHost);

  return 0;
}

#ifdef USECUDA
int cboundary_surface::bcvalues()
{
  //fields->forwardDevice();

  int gridi, gridj;
  const int blocki = 128;
  const int blockj = 2;

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

  // start with retrieving the stability information
  if(model->thermo->getsw() == "0")
  {
    master->printMessage("Neutral boundary_surface not yet implemented on GPU\n");
    //stability_neutral(ustar, obuk,
    //                  fields->u->data, fields->v->data,
    //                  fields->u->databot, fields->v->databot,
    //                  fields->sd["tmp1"]->data, grid->z);
  }
  else
  {
    // store the buoyancy in tmp1
    model->thermo->getbuoyancysurf(fields->sd["tmp1"]);

    // Calculate dutot in tmp2
    boundary_surface_dutot<<<gridGPU, blockGPU>>>(&fields->sd["tmp2"]->data_g[offs], 
                                                  &fields->u->data_g[offs],    &fields->v->data_g[offs],
                                                  &fields->u->databot_g[offs], &fields->v->databot_g[offs],
                                                  grid->istart, grid->jstart, grid->kstart,
                                                  grid->iend,   grid->jend,   grid->icellsp, grid->ijcellsp);

    // 2D cyclic boundaries on dutot  
    grid->boundary_cyclic2d_g(&fields->sd["tmp2"]->data_g[offs]);

    // Calculate ustar and Obukhov length, including ghost cells
    boundary_surface_stability<<<gridGPU2, blockGPU2>>>(&ustar_g[offs], &obuk_g[offs], 
                                                  &fields->sd["tmp1"]->data_g[offs], &fields->sd["tmp1"]->databot_g[offs], &fields->sd["tmp1"]->datafluxbot_g[offs],
                                                  &fields->sd["tmp2"]->data_g[offs], z0m, z0h, grid->z[grid->kstart],
                                                  grid->icells, grid->jcells, grid->kstart, grid->icellsp, grid->ijcellsp, mbcbot, thermobc); 

    // Calculate surface momentum fluxes, excluding ghost cells
    boundary_surface_surfm_flux<<<gridGPU, blockGPU>>>(&fields->u->datafluxbot_g[offs], &fields->v->datafluxbot_g[offs],
                                                       &fields->u->data_g[offs],        &fields->v->data_g[offs],
                                                       &fields->u->databot_g[offs],     &fields->v->databot_g[offs], 
                                                       &ustar_g[offs], &obuk_g[offs], grid->z[grid->kstart], z0m,
                                                       grid->istart, grid->jstart, grid->kstart,
                                                       grid->iend,   grid->jend,   grid->icellsp, grid->ijcellsp, mbcbot);

    // 2D cyclic boundaries on the surface fluxes  
    grid->boundary_cyclic2d_g(&fields->u->datafluxbot_g[offs]);
    grid->boundary_cyclic2d_g(&fields->v->datafluxbot_g[offs]);

    // Calculate surface gradients, including ghost cells
    boundary_surface_surfm_grad<<<gridGPU2, blockGPU2>>>(&fields->u->datagradbot_g[offs], &fields->v->datagradbot_g[offs],
                                                         &fields->u->data_g[offs],        &fields->v->data_g[offs],
                                                         &fields->u->databot_g[offs],     &fields->v->databot_g[offs],
                                                         grid->z[grid->kstart], grid->icells, grid->jcells, grid->kstart, grid->icellsp, grid->ijcellsp);  

    // Calculate scalar fluxes, gradients and/or values, including ghost cells
    for(fieldmap::const_iterator it=fields->sp.begin(); it!=fields->sp.end(); ++it)
    {
      boundary_surface_surfs<<<gridGPU2, blockGPU2>>>(&it->second->datafluxbot_g[offs], &it->second->datagradbot_g[offs],
                                                      &it->second->databot_g[offs],     &it->second->data_g[offs],
                                                      &ustar_g[offs], &obuk_g[offs], grid->z[grid->kstart], z0h,            
                                                      grid->icells, grid->jcells, grid->kstart,
                                                      grid->icellsp, grid->ijcellsp, sbc[it->first]->bcbot);
    }
  }

  
  //fields->backwardDevice();
  //backwardDevice();

  return 0;
}
#endif
