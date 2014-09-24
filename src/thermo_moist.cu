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
#include "thermo_moist.h"
#include "defines.h"
#include "constants.h"
#include "master.h"

using namespace constants;

__device__ double thermo_moist_exn(double p)
{
  return pow((p/p0),(Rd/cp));
}

__device__ double thermo_moist_bu(double p, double s, double qt, double ql, double thvref)
{
  return grav * ((s + Lv*ql/(cp*thermo_moist_exn(p))) * (1. - (1. - Rv/Rd)*qt - Rv/Rd*ql) - thvref) / thvref;
}

__device__ double thermo_moist_bunoql(double s, double qt, double thvref)
{
  return grav * (s * (1. - (1. - Rv/Rd)*qt) - thvref) / thvref;
}

__device__ double thermo_moist_bufluxnoql(double s, double sflux, double qt, double qtflux, double thvref)
{
  return grav/thvref * (sflux * (1. - (1.-Rv/Rd)*qt) - (1.-Rv/Rd)*s*qtflux);
}

__device__ double thermo_moist_esl(double T)
{
  const double x=fmax(-80.,T-T0);
  return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))));
}

__device__ double thermo_moist_rslf(double p, double T)
{
  double esl = thermo_moist_esl(T);
  return ep*esl/(p-(1-ep)*esl);
}

__device__ double thermo_moist_calcql(double s, double qt, double p, double exn)
{
  int niter = 0; //, nitermax = 5;
  double ql, tl, tnr_old = 1.e9, tnr, qs;
  tl = s * exn;
  tnr = tl;
  while (fabs(tnr-tnr_old)/tnr_old> 1e-5)// && niter < nitermax)
  {
    ++niter;
    tnr_old = tnr;
    qs = thermo_moist_rslf(p,tnr);
    tnr = tnr - (tnr+(Lv/cp)*qs-tl-(Lv/cp)*qt)/(1+(pow(Lv,2)*qs)/ (Rv*cp*pow(tnr,2)));
  }
  ql = fmax(0.,qt-qs);
  return ql;
}


__global__ void thermo_moist_calcbuoyancytend_2nd(double * __restrict__ wt, double * __restrict__ th, double * __restrict__ qt,
                                                  double * __restrict__ thvrefh, double * __restrict__ exnh, double * __restrict__ ph,  
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

    // Half level temperature and moisture content
    double thh = 0.5 * (th[ijk-kk] + th[ijk]);        // Half level liq. water pot. temp.
    double qth = 0.5 * (qt[ijk-kk] + qt[ijk]);        // Half level specific hum.
    double tl  = thh * exnh[k];                       // Half level liq. water temp.
    double ql  = qth - thermo_moist_rslf(ph[k], tl);

    // If ql(Tl)>0, saturation adjustment routine needed. 
    if(ql > 0)
      ql = thermo_moist_calcql(thh, qth, ph[k], exnh[k]);
    else
      ql = 0.;

    // Calculate tendency
    wt[ijk] += thermo_moist_bu(ph[k], thh, qth, ql, thvrefh[k]);
  }
}

__global__ void thermo_moist_calcbuoyancy(double * __restrict__ b,  double * __restrict__ th, 
                                          double * __restrict__ qt, double * __restrict__ thvref,
                                          double * __restrict__ p,  double * __restrict__ exn, 
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
    double ql = thermo_moist_calcql(th[ijk], qt[ijk], p[k], exn[k]);
    b[ijk] = thermo_moist_bu(p[k], th[ijk], qt[ijk], ql, thvref[k]);
  }
}

__global__ void thermo_moist_calcbuoyancybot(double * __restrict__ b,      double * __restrict__ bbot,
                                             double * __restrict__ th,     double * __restrict__ thbot, 
                                             double * __restrict__ qt,     double * __restrict__ qtbot,
                                             double * __restrict__ thvref, double * __restrict__ thvrefh,
                                             int kstart, int icells, int jcells,  
                                             int jj, int kk)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x; 
  int j = blockIdx.y*blockDim.y + threadIdx.y; 

  if(i < icells && j < jcells)
  {
    const int ij  = i + j*jj;
    const int ijk = i + j*jj + kstart*kk;

    bbot[ij ] = thermo_moist_bunoql(thbot[ij], qtbot[ij], thvrefh[kstart]);
    b   [ijk] = thermo_moist_bunoql(th[ijk],   qt[ijk],   thvref[kstart]);
  }
}

__global__ void thermo_moist_calcbuoyancyfluxbot(double * __restrict__ bfluxbot,
                                                 double * __restrict__ thbot, double * __restrict__ thfluxbot, 
                                                 double * __restrict__ qtbot, double * __restrict__ qtfluxbot,
                                                 double * __restrict__ thvrefh, 
                                                 int kstart, int icells, int jcells,  
                                                 int jj, int kk)
{
  int i = blockIdx.x*blockDim.x + threadIdx.x; 
  int j = blockIdx.y*blockDim.y + threadIdx.y; 

  if(i < icells && j < jcells)
  {
    const int ij  = i + j*jj;
    bfluxbot[ij] = thermo_moist_bufluxnoql(thbot[ij], thfluxbot[ij], qtbot[ij], qtfluxbot[ij], thvrefh[kstart]);
  }
}

__global__ void thermo_moist_calcN2(double * __restrict__ N2, double * __restrict__ th,
                                    double * __restrict__ thvref, double * __restrict__ dzi, 
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
    N2[ijk] = grav/thvref[k]*0.5*(th[ijk+kk] - th[ijk-kk])*dzi[k];
  }
}

__global__ void thermo_moist_calcqlfield(double * __restrict__ ql, double * __restrict__ th, double * __restrict__ qt,
                                         double * __restrict__ exn, double * __restrict__ p,  
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
    ql[ijk] = thermo_moist_calcql(th[ijk], qt[ijk], p[k], exn[k]);
  }
}

int cthermo_moist::prepareDevice()
{
  const int nmemsize = grid->kcells*sizeof(double);

  // Allocate fields for Boussinesq and anelastic solver
  cudaMalloc(&thvref_g,  nmemsize);
  cudaMalloc(&thvrefh_g, nmemsize);
  cudaMalloc(&pref_g,    nmemsize);
  cudaMalloc(&prefh_g,   nmemsize);
  cudaMalloc(&exnref_g,  nmemsize);
  cudaMalloc(&exnrefh_g, nmemsize);

  // Copy fields to device
  cudaMemcpy(thvref_g,  thvref,  nmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(thvrefh_g, thvrefh, nmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(pref_g,    pref,    nmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(prefh_g,   prefh,   nmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(exnref_g,  exnref,  nmemsize, cudaMemcpyHostToDevice);
  cudaMemcpy(exnrefh_g, exnrefh, nmemsize, cudaMemcpyHostToDevice);

  return 0;
}

int cthermo_moist::clearDevice()
{
  cudaFree(thvref_g);
  cudaFree(thvrefh_g);
  cudaFree(pref_g);
  cudaFree(prefh_g);
  cudaFree(exnref_g);
  cudaFree(exnrefh_g);

  return 0;
}

#ifdef USECUDA
int cthermo_moist::exec()
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);
  
  const int offs = grid->memoffset;

  // Re-calculate hydrostatic pressure and exner, pass dummy as rhoref,thvref to prevent overwriting base state 
  //double * restrict tmp2 = fields->s["tmp2"]->data;
  //if(swupdatebasestate)
  //  calcbasestate(pref, prefh, &tmp2[0*kcells], &tmp2[1*kcells], &tmp2[2*kcells], &tmp2[3*kcells], exnref, exnrefh, 
  //                fields->s["s"]->datamean, fields->s["qt"]->datamean);

  if(grid->swspatialorder== "2")
    thermo_moist_calcbuoyancytend_2nd<<<gridGPU, blockGPU>>>(&fields->wt->data_g[offs], &fields->s["s"]->data_g[offs], 
                                                             &fields->s["qt"]->data_g[offs], thvrefh_g, exnrefh_g, prefh_g,  
                                                             grid->istart,  grid->jstart, grid->kstart+1,
                                                             grid->iend,    grid->jend,   grid->kend,
                                                             grid->icellsp, grid->ijcellsp);
  //else if(grid->swspatialorder == "4")
  //  calcbuoyancytend_4th(fields->wt->data, fields->s["th"]->data, threfh);

  return 0;
}
#endif

#ifdef USECUDA
int cthermo_moist::getthermofield(cfield3d *fld, cfield3d *tmp, std::string name)
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kcells);
  dim3 blockGPU(blocki, blockj, 1);

  dim3 gridGPU2 (gridi, gridj, grid->kmax);
  dim3 blockGPU2(blocki, blockj, 1);
  
  const int offs = grid->memoffset;

  if(name == "b")
    thermo_moist_calcbuoyancy<<<gridGPU, blockGPU>>>(&fld->data_g[offs], &fields->s["s"]->data_g[offs], &fields->s["qt"]->data_g[offs],
                                                     thvref_g, pref_g, exnref_g,
                                                     grid->istart, grid->jstart, grid->iend, grid->jend, grid->kcells,
                                                     grid->icellsp, grid->ijcellsp);
  else if(name == "ql")
    thermo_moist_calcqlfield<<<gridGPU2, blockGPU2>>>(&fld->data_g[offs], &fields->s["s"]->data_g[offs], &fields->s["qt"]->data_g[offs], exnref_g, pref_g, 
                                                      grid->istart,  grid->jstart,  grid->kstart, 
                                                      grid->iend,    grid->jend,    grid->kend,
                                                      grid->icellsp, grid->ijcellsp);
  else if(name == "N2")
    thermo_moist_calcN2<<<gridGPU2, blockGPU2>>>(&fld->data_g[offs], &fields->s["s"]->data_g[offs], thvref_g, grid->dzi_g, 
                                                 grid->istart, grid->jstart, grid->kstart, 
                                                 grid->iend,   grid->jend,   grid->kend,
                                                 grid->icellsp, grid->ijcellsp);
  else
    return 1;

  return 0;
}
#endif

#ifdef USECUDA
int cthermo_moist::getbuoyancyfluxbot(cfield3d *bfield)
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->icells/blocki + (grid->icells%blocki > 0);
  const int gridj  = grid->jcells/blockj + (grid->jcells%blockj > 0);

  dim3 gridGPU (gridi, gridj, 1);
  dim3 blockGPU(blocki, blockj, 1);
  
  const int offs = grid->memoffset;

  thermo_moist_calcbuoyancyfluxbot<<<gridGPU, blockGPU>>>(&bfield->datafluxbot_g[offs], 
                                                          &fields->s["s"] ->databot_g[offs], &fields->s["s"] ->datafluxbot_g[offs], 
                                                          &fields->s["qt"]->databot_g[offs], &fields->s["qt"]->datafluxbot_g[offs], 
                                                          thvrefh_g, grid->kstart, grid->icells, grid->jcells, 
                                                          grid->icellsp, grid->ijcellsp);
  return 0;
}
#endif

#ifdef USECUDA
int cthermo_moist::getbuoyancysurf(cfield3d *bfield)
{
  const int blocki = 128;
  const int blockj = 2;
  const int gridi  = grid->icells/blocki + (grid->icells%blocki > 0);
  const int gridj  = grid->jcells/blockj + (grid->jcells%blockj > 0);

  dim3 gridGPU (gridi, gridj, 1);
  dim3 blockGPU(blocki, blockj, 1);
  
  const int offs = grid->memoffset;

  thermo_moist_calcbuoyancybot<<<gridGPU, blockGPU>>>(&bfield->data_g[offs], &bfield->databot_g[offs], 
                                                      &fields->s["s"] ->data_g[offs], &fields->s["s"] ->databot_g[offs],
                                                      &fields->s["qt"]->data_g[offs], &fields->s["qt"]->databot_g[offs],
                                                      thvref_g, thvrefh_g, grid->kstart, grid->icells, grid->jcells, 
                                                      grid->icellsp, grid->ijcellsp);

  thermo_moist_calcbuoyancyfluxbot<<<gridGPU, blockGPU>>>(&bfield->datafluxbot_g[offs], 
                                                          &fields->s["s"] ->databot_g[offs], &fields->s["s"] ->datafluxbot_g[offs], 
                                                          &fields->s["qt"]->databot_g[offs], &fields->s["qt"]->datafluxbot_g[offs], 
                                                          thvrefh_g, grid->kstart, grid->icells, grid->jcells, 
                                                          grid->icellsp, grid->ijcellsp);
  return 0;
}
#endif
