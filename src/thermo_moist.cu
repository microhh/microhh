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
#include "fd.h"
#include "master.h"
#include "tools.h"

using namespace constants;
using namespace fd::o4;

__device__ double thermo_moist_interp2(double a, double b)
{
  return 0.5*(a + b);
}

__device__ double thermo_moist_interp4(double a, double b, double c, double d) 
{
  return ci0*a + ci1*b + ci2*c + ci3*d;
}

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

/*
s = liquid water potential temperature [K]
qt = moisture mixing ratio [kg kg-1]
p = pressure [pa]
exn = exner [-]
*/
__device__ double thermo_moist_calcql(double s, double qt, double p, double exn)
{
  double tl = s * exn;  // Liquid water temperature

  // Calculate if q-qs(Tl) <= 0. If so, return 0. Else continue with saturation adjustment
  if(qt-thermo_moist_rslf(p, tl) <= 0)
    return 0.;

  int niter = 0; //, nitermax = 5;
  double ql, tnr_old = 1.e9, tnr, qs;
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

/* This routine needs to be solved level by level, so one thread does everything */
template <int swspatialorder>
__global__ void thermo_moist_calcbasestate(double * __restrict__ pref,     double * __restrict__ prefh,
                                           double * __restrict__ rho,      double * __restrict__ rhoh,
                                           double * __restrict__ thv,      double * __restrict__ thvh,
                                           double * __restrict__ ex,       double * __restrict__ exh,
                                           double * __restrict__ thlmean,  double * __restrict__ qtmean,
                                           double * __restrict__ z,        double * __restrict__ dz,
                                           double * __restrict__ dzh,
                                           double pbot, int kstart, int kend)
{
  double ssurf, qtsurf, ql, si, qti, qli;
  double rdcp = Rd/cp;

  if(swspatialorder == 2)
  {
    ssurf  = thermo_moist_interp2(thlmean[kstart-1], thlmean[kstart]);
    qtsurf = thermo_moist_interp2(qtmean[kstart-1],  qtmean[kstart]);
  }
  else if(swspatialorder == 4)
  {
    ssurf  = thermo_moist_interp4(thlmean[kstart-2], thlmean[kstart-1], thlmean[kstart], thlmean[kstart+1]);
    qtsurf = thermo_moist_interp4(qtmean[kstart-2],  qtmean[kstart-1],  qtmean[kstart],  qtmean[kstart+1]);
  }

  // Calculate surface (half=kstart) values
  exh[kstart]   = thermo_moist_exn(pbot);
  ql            = thermo_moist_calcql(ssurf,qtsurf,pbot,exh[kstart]); 
  thvh[kstart]  = (ssurf + Lv*ql/(cp*exh[kstart])) * (1. - (1. - Rv/Rd)*qtsurf - Rv/Rd*ql);
  prefh[kstart] = pbot;
  rhoh[kstart]  = pbot / (Rd * exh[kstart] * thvh[kstart]);

  // First full grid level pressure
  pref[kstart] = pow((pow(pbot,rdcp) - grav * pow(p0,rdcp) * z[kstart] / (cp * thvh[kstart])),(1./rdcp)); 

  for(int k=kstart+1; k<kend+1; k++)
  {
    // 1. Calculate values at full level below zh[k] 
    ex[k-1]  = thermo_moist_exn(pref[k-1]);
    ql       = thermo_moist_calcql(thlmean[k-1],qtmean[k-1],pref[k-1],ex[k-1]); 
    thv[k-1] = (thlmean[k-1] + Lv*ql/(cp*ex[k-1])) * (1. - (1. - Rv/Rd)*qtmean[k-1] - Rv/Rd*ql); 
    rho[k-1] = pref[k-1] / (Rd * ex[k-1] * thv[k-1]);
 
    // 2. Calculate half level pressure at zh[k] using values at z[k-1]
    prefh[k] = pow((pow(prefh[k-1],rdcp) - grav * pow(p0,rdcp) * dz[k-1] / (cp * thv[k-1])),(1./rdcp));

    // 3. Interpolate conserved variables to zh[k] and calculate virtual temp and ql
    if(swspatialorder == 2)
    {
      si     = thermo_moist_interp2(thlmean[k-1],thlmean[k]);
      qti    = thermo_moist_interp2(qtmean[k-1],qtmean[k]);
    }
    else if(swspatialorder == 4)
    {
      si     = thermo_moist_interp4(thlmean[k-2],thlmean[k-1],thlmean[k],thlmean[k+1]);
      qti    = thermo_moist_interp4(qtmean[k-2],qtmean[k-1],qtmean[k],qtmean[k+1]);
    }

    exh[k]   = thermo_moist_exn(prefh[k]);
    qli      = thermo_moist_calcql(si,qti,prefh[k],exh[k]);
    thvh[k]  = (si + Lv*qli/(cp*exh[k])) * (1. - (1. - Rv/Rd)*qti - Rv/Rd*qli); 
    rhoh[k]  = prefh[k] / (Rd * exh[k] * thvh[k]); 

    // 4. Calculate full level pressure at z[k]
    pref[k]  = pow((pow(pref[k-1],rdcp) - grav * pow(p0,rdcp) * dzh[k] / (cp * thvh[k])),(1./rdcp)); 
  }

  // Fill bottom and top full level ghost cells 
  if(swspatialorder == 2)
  {
    pref[kstart-1] = 2.*prefh[kstart] - pref[kstart];
    pref[kend]     = 2.*prefh[kend]   - pref[kend-1];
  }
  else if(swspatialorder == 4)
  {
    pref[kstart-1] = (8./3.)*prefh[kstart] - 2.*pref[kstart] + (1./3.)*pref[kstart+1];
    pref[kstart-2] = 8.*prefh[kstart]      - 9.*pref[kstart] + 2.*pref[kstart+1];
    pref[kend]     = (8./3.)*prefh[kend]   - 2.*pref[kend-1] + (1./3.)*pref[kend-2];
    pref[kend+1]   = 8.*prefh[kend]        - 9.*pref[kend-1] + 2.*pref[kend-2];
  }
}

/* This routine needs to be solved level by level, so one thread does everything */
template <int swspatialorder>
__global__ void thermo_moist_calchydropres(double * __restrict__ pref,     double * __restrict__ prefh,
                                           double * __restrict__ ex,       double * __restrict__ exh,
                                           double * __restrict__ thlmean,  double * __restrict__ qtmean,
                                           double * __restrict__ z,        double * __restrict__ dz,
                                           double * __restrict__ dzh,
                                           double pbot, int kstart, int kend)
{
  double ssurf, qtsurf, ql, si, qti, qli, thvh, thv;
  double rdcp = Rd/cp;

  if(swspatialorder == 2)
  {
    ssurf  = thermo_moist_interp2(thlmean[kstart-1], thlmean[kstart]);
    qtsurf = thermo_moist_interp2(qtmean[kstart-1],  qtmean[kstart]);
  }
  else if(swspatialorder == 4)
  {
    ssurf  = thermo_moist_interp4(thlmean[kstart-2], thlmean[kstart-1], thlmean[kstart], thlmean[kstart+1]);
    qtsurf = thermo_moist_interp4(qtmean[kstart-2],  qtmean[kstart-1],  qtmean[kstart],  qtmean[kstart+1]);
  }

  // Calculate surface (half=kstart) values
  ql            = thermo_moist_calcql(ssurf,qtsurf,pbot,exh[kstart]); 
  thvh          = (ssurf + Lv*ql/(cp*exh[kstart])) * (1. - (1. - Rv/Rd)*qtsurf - Rv/Rd*ql);
  prefh[kstart] = pbot;

  // First full grid level pressure
  pref[kstart] = pow((pow(pbot,rdcp) - grav * pow(p0,rdcp) * z[kstart] / (cp * thvh)),(1./rdcp)); 
  for(int k=kstart+1; k<kend+1; k++)
  {
    // 1. Calculate values at full level below zh[k] 
    ex[k-1]  = thermo_moist_exn(pref[k-1]);
    ql       = thermo_moist_calcql(thlmean[k-1],qtmean[k-1],pref[k-1],ex[k-1]); 
    thv      = (thlmean[k-1] + Lv*ql/(cp*ex[k-1])) * (1. - (1. - Rv/Rd)*qtmean[k-1] - Rv/Rd*ql); 
 
    // 2. Calculate half level pressure at zh[k] using values at z[k-1]
    prefh[k] = pow((pow(prefh[k-1],rdcp) - grav * pow(p0,rdcp) * dz[k-1] / (cp * thv)),(1./rdcp));

    // 3. Interpolate conserved variables to zh[k] and calculate virtual temp and ql
    if(swspatialorder == 2)
    {
      si     = thermo_moist_interp2(thlmean[k-1],thlmean[k]);
      qti    = thermo_moist_interp2(qtmean[k-1],qtmean[k]);
    }
    else if(swspatialorder == 4)
    {
      si     = thermo_moist_interp4(thlmean[k-2],thlmean[k-1],thlmean[k],thlmean[k+1]);
      qti    = thermo_moist_interp4(qtmean[k-2],qtmean[k-1],qtmean[k],qtmean[k+1]);
    }

    exh[k]   = thermo_moist_exn(prefh[k]);
    qli      = thermo_moist_calcql(si,qti,prefh[k],exh[k]);
    thvh     = (si + Lv*qli/(cp*exh[k])) * (1. - (1. - Rv/Rd)*qti - Rv/Rd*qli); 

    // 4. Calculate full level pressure at z[k]
    pref[k]  = pow((pow(pref[k-1],rdcp) - grav * pow(p0,rdcp) * dzh[k] / (cp * thvh)),(1./rdcp)); 
  }

  // Fill bottom and top full level ghost cells 
  if(swspatialorder == 2)
  {
    pref[kstart-1] = 2.*prefh[kstart] - pref[kstart];
    pref[kend]     = 2.*prefh[kend]   - pref[kend-1];
  }
  else if(swspatialorder == 4)
  {
    pref[kstart-1] = (8./3.)*prefh[kstart] - 2.*pref[kstart] + (1./3.)*pref[kstart+1];
    pref[kstart-2] = 8.*prefh[kstart]      - 9.*pref[kstart] + 2.*pref[kstart+1];
    pref[kend]     = (8./3.)*prefh[kend]   - 2.*pref[kend-1] + (1./3.)*pref[kend-2];
    pref[kend+1]   = 8.*prefh[kend]        - 9.*pref[kend-1] + 2.*pref[kend-2];
  }
}

void ThermoMoist::prepareDevice()
{
  const int nmemsize = grid->kcells*sizeof(double);

  // Allocate fields for Boussinesq and anelastic solver
  cudaSafeCall(cudaMalloc(&thvref_g,  nmemsize));
  cudaSafeCall(cudaMalloc(&thvrefh_g, nmemsize));
  cudaSafeCall(cudaMalloc(&pref_g,    nmemsize));
  cudaSafeCall(cudaMalloc(&prefh_g,   nmemsize));
  cudaSafeCall(cudaMalloc(&exnref_g,  nmemsize));
  cudaSafeCall(cudaMalloc(&exnrefh_g, nmemsize));

  // Copy fields to device
  cudaSafeCall(cudaMemcpy(thvref_g,  thvref,  nmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(thvrefh_g, thvrefh, nmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(pref_g,    pref,    nmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(prefh_g,   prefh,   nmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(exnref_g,  exnref,  nmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(exnrefh_g, exnrefh, nmemsize, cudaMemcpyHostToDevice));
}

void ThermoMoist::clearDevice()
{
  cudaSafeCall(cudaFree(thvref_g ));
  cudaSafeCall(cudaFree(thvrefh_g));
  cudaSafeCall(cudaFree(pref_g   ));
  cudaSafeCall(cudaFree(prefh_g  ));
  cudaSafeCall(cudaFree(exnref_g ));
  cudaSafeCall(cudaFree(exnrefh_g));
}

#ifdef USECUDA
void ThermoMoist::exec()
{
  const int blocki = cuda::blockSizeI;
  const int blockj = cuda::blockSizeJ;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);
  
  const int offs = grid->memoffset;
  int kk = grid->kcells;

  // Re-calculate hydrostatic pressure and exner, pass dummy as rhoref,thvref to prevent overwriting base state
  double * restrict tmp2 = fields->s["tmp2"]->data_g;
  if(swupdatebasestate)
  {
    if(grid->swspatialorder == "2")
      thermo_moist_calchydropres<2><<<1, 1>>>(pref_g, prefh_g, exnref_g, exnrefh_g, 
                                              fields->s["s"]->datamean_g, fields->s["qt"]->datamean_g, 
                                              grid->z_g, grid->dz_g, grid->dzh_g, pbot, grid->kstart, grid->kend);
    else if(grid->swspatialorder == "4")
      thermo_moist_calchydropres<4><<<1, 1>>>(pref_g, prefh_g, exnref_g, exnrefh_g, 
                                              fields->s["s"]->datamean_g, fields->s["qt"]->datamean_g, 
                                              grid->z_g, grid->dz_g, grid->dzh_g, pbot, grid->kstart, grid->kend);
    cudaCheckError();
  }

  if(grid->swspatialorder== "2")
  {
    thermo_moist_calcbuoyancytend_2nd<<<gridGPU, blockGPU>>>(&fields->wt->data_g[offs], &fields->s["s"]->data_g[offs], 
                                                             &fields->s["qt"]->data_g[offs], thvrefh_g, exnrefh_g, prefh_g,  
                                                             grid->istart,  grid->jstart, grid->kstart+1,
                                                             grid->iend,    grid->jend,   grid->kend,
                                                             grid->icellsp, grid->ijcellsp);
    cudaCheckError();
  }
  else if(grid->swspatialorder == "4")
  {
    master->printMessage("4th order thermo_moist not (yet) implemented\n");  
  //  calcbuoyancytend_4th(fields->wt->data, fields->s["th"]->data, threfh);
    throw 1;
  }
}
#endif

#ifdef USECUDA
void ThermoMoist::getThermoField(Field3d *fld, Field3d *tmp, std::string name)
{
  const int blocki = cuda::blockSizeI;
  const int blockj = cuda::blockSizeJ;
  const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
  const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

  dim3 gridGPU (gridi, gridj, grid->kcells);
  dim3 blockGPU(blocki, blockj, 1);

  dim3 gridGPU2 (gridi, gridj, grid->kmax);
  dim3 blockGPU2(blocki, blockj, 1);
  
  const int offs = grid->memoffset;
  int kk = grid->kcells;

  // BvS: getthermofield() is called from subgrid-model, before thermo(), so re-calculate the hydrostatic pressure
  // Pass dummy as rhoref,thvref to prevent overwriting base state 
  double * restrict tmp2 = fields->s["tmp2"]->data_g;
  if(swupdatebasestate && (name == "b" || name == "ql"))
  {
    if(grid->swspatialorder == "2")
      thermo_moist_calchydropres<2><<<1, 1>>>(pref_g, prefh_g, exnref_g, exnrefh_g, 
                                              fields->s["s"]->datamean_g, fields->s["qt"]->datamean_g, 
                                              grid->z_g, grid->dz_g, grid->dzh_g, pbot, grid->kstart, grid->kend);
    else if(grid->swspatialorder == "4")
      thermo_moist_calchydropres<4><<<1, 1>>>(pref_g, prefh_g, exnref_g, exnrefh_g, 
                                              fields->s["s"]->datamean_g, fields->s["qt"]->datamean_g, 
                                              grid->z_g, grid->dz_g, grid->dzh_g, pbot, grid->kstart, grid->kend);
    cudaCheckError();
  }

  if(name == "b")
  {
    thermo_moist_calcbuoyancy<<<gridGPU, blockGPU>>>(&fld->data_g[offs], &fields->s["s"]->data_g[offs], &fields->s["qt"]->data_g[offs],
                                                     thvref_g, pref_g, exnref_g,
                                                     grid->istart, grid->jstart, grid->iend, grid->jend, grid->kcells,
                                                     grid->icellsp, grid->ijcellsp);
    cudaCheckError();
  }
  else if(name == "ql")
  {
    thermo_moist_calcqlfield<<<gridGPU2, blockGPU2>>>(&fld->data_g[offs], &fields->s["s"]->data_g[offs], &fields->s["qt"]->data_g[offs], exnref_g, pref_g, 
                                                      grid->istart,  grid->jstart,  grid->kstart, 
                                                      grid->iend,    grid->jend,    grid->kend,
                                                      grid->icellsp, grid->ijcellsp);
    cudaCheckError();
  }
  else if(name == "N2")
  {
    thermo_moist_calcN2<<<gridGPU2, blockGPU2>>>(&fld->data_g[offs], &fields->s["s"]->data_g[offs], thvref_g, grid->dzi_g, 
                                                 grid->istart, grid->jstart, grid->kstart, 
                                                 grid->iend,   grid->jend,   grid->kend,
                                                 grid->icellsp, grid->ijcellsp);
    cudaCheckError();
  }
  else
    throw 1;
}
#endif

#ifdef USECUDA
void ThermoMoist::getBuoyancyFluxbot(Field3d *bfield)
{
  const int blocki = cuda::blockSizeI;
  const int blockj = cuda::blockSizeJ;
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
  cudaCheckError();
}
#endif

#ifdef USECUDA
void ThermoMoist::getBuoyancySurf(Field3d *bfield)
{
  const int blocki = cuda::blockSizeI;
  const int blockj = cuda::blockSizeJ;
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
  cudaCheckError();

  thermo_moist_calcbuoyancyfluxbot<<<gridGPU, blockGPU>>>(&bfield->datafluxbot_g[offs], 
                                                          &fields->s["s"] ->databot_g[offs], &fields->s["s"] ->datafluxbot_g[offs], 
                                                          &fields->s["qt"]->databot_g[offs], &fields->s["qt"]->datafluxbot_g[offs], 
                                                          thvrefh_g, grid->kstart, grid->icells, grid->jcells, 
                                                          grid->icellsp, grid->ijcellsp);
  cudaCheckError();
}
#endif
