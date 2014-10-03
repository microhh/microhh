/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
 *
 * The heptadiagonal matrix solver is
 * Copyright (c) 2014 Juan Pedro Mellado
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
#include <cmath>
#include <algorithm>
#include <fftw3.h>
#include "master.h"
#include "grid.h"
#include "fields.h"
#include "pres_4.h"
#include "defines.h"
#include "fd.h"
#include "model.h"

using namespace fd::o4;

Pres_4::Pres_4(Model *modelin, Input *inputin) : Pres(modelin, inputin)
{
  m1 = 0;
  m2 = 0;
  m3 = 0;
  m4 = 0;
  m5 = 0;
  m6 = 0;
  m7 = 0;

  bmati = 0;
  bmatj = 0;
}

Pres_4::~Pres_4()
{
  delete[] m1;
  delete[] m2;
  delete[] m3;
  delete[] m4;
  delete[] m5;
  delete[] m6;
  delete[] m7;

  delete[] bmati;
  delete[] bmatj;

  #ifdef USECUDA
  clearDevice();
  #endif
}

#ifndef USECUDA
void Pres_4::exec(double dt)
{
  // 1. Create the input for the pressure solver.
  pres_in(fields->sd["p"]->data,
          fields->u ->data, fields->v ->data, fields->w ->data,
          fields->ut->data, fields->vt->data, fields->wt->data, 
          grid->dzi4, dt);

  // 2. Solve the Poisson equation using FFTs and a heptadiagonal solver

  /* Find the thickness of a vectorizable slice. There is a need for 8 slices for the pressure
     solver and we use two three dimensional temp fields, so there are 4 slices per field.
     The thickness is therefore jblock/4. Since there are always three ghost cells, even in a 2D
     run the fields are large enough. */
  // const int jslice = std::max(grid->jblock/4, 1);

  /* The CPU version gives the best performance in case jslice = 1, due to cache misses.
     In case this value will be set to larger than 1, checks need to be build in for out of bounds
     reads in case jblock does not divide by 4. */
  const int jslice = 1;

  double *tmp2 = fields->atmp["tmp2"]->data;
  double *tmp3 = fields->atmp["tmp3"]->data;

  const int ns = grid->iblock*jslice*(grid->kmax+4);

  pres_solve(fields->sd["p"]->data, fields->atmp["tmp1"]->data, grid->dz,
             m1, m2, m3, m4,
             m5, m6, m7,
             &tmp2[0*ns], &tmp2[1*ns], &tmp2[2*ns], &tmp2[3*ns], 
             &tmp3[0*ns], &tmp3[1*ns], &tmp3[2*ns], &tmp3[3*ns], 
             bmati, bmatj,
             jslice);

  // 3. Get the pressure tendencies from the pressure field.
  pres_out(fields->ut->data, fields->vt->data, fields->wt->data, 
           fields->sd["p"]->data, grid->dzhi4);
}

double Pres_4::check()
{
  double divmax = 0.;

  divmax = calcdivergence(fields->u->data, fields->v->data, fields->w->data, grid->dzi4);

  return divmax;
}
#endif

void Pres_4::init()
{
  int imax, jmax, kmax;
  int itot, jtot, kstart;

  itot   = grid->itot;
  jtot   = grid->jtot;
  imax   = grid->imax;
  jmax   = grid->jmax;
  kmax   = grid->kmax;
  kstart = grid->kstart;

  bmati = new double[itot];
  bmatj = new double[jtot];

  m1 = new double[kmax];
  m2 = new double[kmax];
  m3 = new double[kmax];
  m4 = new double[kmax];
  m5 = new double[kmax];
  m6 = new double[kmax];
  m7 = new double[kmax];
}

void Pres_4::setvalues()
{
  int imax, jmax, kmax;
  int itot, jtot, kstart;

  itot   = grid->itot;
  jtot   = grid->jtot;
  imax   = grid->imax;
  jmax   = grid->jmax;
  kmax   = grid->kmax;
  kstart = grid->kstart;

  // compute the modified wave numbers of the 4th order scheme
  double dxidxi = 1./(grid->dx*grid->dx);
  double dyidyi = 1./(grid->dy*grid->dy);

  const double pi = std::acos(-1.);

  for(int j=0; j<jtot/2+1; j++)
    bmatj[j] = ( 2.* (1./576.)    * std::cos(6.*pi*(double)j/(double)jtot)
               - 2.* (54./576.)   * std::cos(4.*pi*(double)j/(double)jtot)
               + 2.* (783./576.)  * std::cos(2.*pi*(double)j/(double)jtot)
               -     (1460./576.) ) * dyidyi;

  for(int j=jtot/2+1; j<jtot; j++)
    bmatj[j] = bmatj[jtot-j];

  for(int i=0; i<itot/2+1; i++)
    bmati[i] = ( 2.* (1./576.)    * std::cos(6.*pi*(double)i/(double)itot)
               - 2.* (54./576.)   * std::cos(4.*pi*(double)i/(double)itot)
               + 2.* (783./576.)  * std::cos(2.*pi*(double)i/(double)itot)
               -     (1460./576.) ) * dxidxi;

  for(int i=itot/2+1; i<itot; i++)
    bmati[i] = bmati[itot-i];

  double *dzi4, *dzhi4;
  dzi4  = grid->dzi4;
  dzhi4 = grid->dzhi4;

  int k,kc;
  // create vectors that go into the matrix solver
  // bottom boundary, taking into account that w is mirrored over the wall to conserve global momentum
  k  = 0;
  kc = kstart+k;
  m1[k] = 0.;
  m2[k] = (                 -  27.*dzhi4[kc]                                      ) * dzi4[kc];
  m3[k] = ( -1.*dzhi4[kc+1] + 729.*dzhi4[kc] +  27.*dzhi4[kc+1]                   ) * dzi4[kc];
  m4[k] = ( 27.*dzhi4[kc+1] - 729.*dzhi4[kc] - 729.*dzhi4[kc+1] -  1.*dzhi4[kc+2] ) * dzi4[kc];
  m5[k] = (-27.*dzhi4[kc+1] +  27.*dzhi4[kc] + 729.*dzhi4[kc+1] + 27.*dzhi4[kc+2] ) * dzi4[kc];
  m6[k] = (  1.*dzhi4[kc+1]                  -  27.*dzhi4[kc+1] - 27.*dzhi4[kc+2] ) * dzi4[kc];
  m7[k] = (                                                     +  1.*dzhi4[kc+2] ) * dzi4[kc];
  
  for(int k=1; k<kmax-1; k++)
  {
    kc = kstart+k;
    m1[k] = (   1.*dzhi4[kc-1]                                                       ) * dzi4[kc];
    m2[k] = ( -27.*dzhi4[kc-1] -  27.*dzhi4[kc]                                      ) * dzi4[kc];
    m3[k] = (  27.*dzhi4[kc-1] + 729.*dzhi4[kc] +  27.*dzhi4[kc+1]                   ) * dzi4[kc];
    m4[k] = (  -1.*dzhi4[kc-1] - 729.*dzhi4[kc] - 729.*dzhi4[kc+1] -  1.*dzhi4[kc+2] ) * dzi4[kc];
    m5[k] = (                  +  27.*dzhi4[kc] + 729.*dzhi4[kc+1] + 27.*dzhi4[kc+2] ) * dzi4[kc];
    m6[k] = (                                   -  27.*dzhi4[kc+1] - 27.*dzhi4[kc+2] ) * dzi4[kc];
    m7[k] = (                                                      +  1.*dzhi4[kc+2] ) * dzi4[kc];
  }                                                                                                                                       

  // top boundary, taking into account that w is mirrored over the wall to conserve global momentum
  k  = kmax-1;
  kc = kstart+k;
  m1[k] = (   1.*dzhi4[kc-1]                                                     ) * dzi4[kc];
  m2[k] = ( -27.*dzhi4[kc-1] -  27.*dzhi4[kc]                    +  1.*dzhi4[kc] ) * dzi4[kc];
  m3[k] = (  27.*dzhi4[kc-1] + 729.*dzhi4[kc] +  27.*dzhi4[kc+1] - 27.*dzhi4[kc] ) * dzi4[kc];
  m4[k] = (  -1.*dzhi4[kc-1] - 729.*dzhi4[kc] - 729.*dzhi4[kc+1] + 27.*dzhi4[kc] ) * dzi4[kc];
  m5[k] = (                  +  27.*dzhi4[kc] + 729.*dzhi4[kc+1] -  1.*dzhi4[kc] ) * dzi4[kc];
  m6[k] = (                                   -  27.*dzhi4[kc+1]                 ) * dzi4[kc];
  m7[k] = 0.;
}

void Pres_4::pres_in(double * restrict p, 
                      double * restrict u , double * restrict v , double * restrict w ,
                      double * restrict ut, double * restrict vt, double * restrict wt,
                      double * restrict dzi4, double dt)
{
  int    ijk,ijkp,jjp,kkp;
  int    ii1,ii2,jj1,jj2,kk1,kk2;
  int    igc,jgc,kgc,kmax;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  jjp = grid->imax;
  kkp = grid->imax*grid->jmax;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  igc = grid->igc;
  jgc = grid->jgc;
  kgc = grid->kgc;

  kmax = grid->kmax;

  // Set the cyclic boundary conditions for the tendencies.
  grid->boundary_cyclic(ut);
  grid->boundary_cyclic(vt);
  grid->boundary_cyclic(wt);

  // Set the bc. 
  for(int j=0; j<grid->jmax; j++)
#pragma ivdep
    for(int i=0; i<grid->imax; i++)
    {
      ijk  = i+igc + (j+jgc)*jj1 + kgc*kk1;
      wt[ijk-kk1] = -wt[ijk+kk1];
    }
  for(int j=0; j<grid->jmax; j++)
#pragma ivdep
    for(int i=0; i<grid->imax; i++)
    {
      ijk  = i+igc + (j+jgc)*jj1 + (kmax+kgc)*kk1;
      wt[ijk+kk1] = -wt[ijk-kk1];
    }

  for(int k=0; k<grid->kmax; k++)
    for(int j=0; j<grid->jmax; j++)
#pragma ivdep
      for(int i=0; i<grid->imax; i++)
      {
        ijkp = i + j*jjp + k*kkp;
        ijk  = i+igc + (j+jgc)*jj1 + (k+kgc)*kk1;
        p[ijkp]  = (cg0*(ut[ijk-ii1] + u[ijk-ii1]/dt) + cg1*(ut[ijk] + u[ijk]/dt) + cg2*(ut[ijk+ii1] + u[ijk+ii1]/dt) + cg3*(ut[ijk+ii2] + u[ijk+ii2]/dt)) * cgi*dxi;
        p[ijkp] += (cg0*(vt[ijk-jj1] + v[ijk-jj1]/dt) + cg1*(vt[ijk] + v[ijk]/dt) + cg2*(vt[ijk+jj1] + v[ijk+jj1]/dt) + cg3*(vt[ijk+jj2] + v[ijk+jj2]/dt)) * cgi*dyi;
        p[ijkp] += (cg0*(wt[ijk-kk1] + w[ijk-kk1]/dt) + cg1*(wt[ijk] + w[ijk]/dt) + cg2*(wt[ijk+kk1] + w[ijk+kk1]/dt) + cg3*(wt[ijk+kk2] + w[ijk+kk2]/dt)) * dzi4[k+kgc];
      }
}

void Pres_4::pres_solve(double * restrict p, double * restrict work3d, double * restrict dz,
                         double * restrict m1, double * restrict m2, double * restrict m3, double * restrict m4,
                         double * restrict m5, double * restrict m6, double * restrict m7,
                         double * restrict m1temp, double * restrict m2temp, double * restrict m3temp, double * restrict m4temp,
                         double * restrict m5temp, double * restrict m6temp, double * restrict m7temp, double * restrict ptemp,
                         double * restrict bmati, double * restrict bmatj,
                         const int jslice)
{
  int jj,kk,ijk;
  int imax,jmax,kmax;
  int itot,jtot;
  int iblock,jblock,kblock;
  int igc,jgc,kgc;
  int iindex,jindex;

  imax   = grid->imax;
  jmax   = grid->jmax;
  kmax   = grid->kmax;
  itot   = grid->itot;
  jtot   = grid->jtot;
  iblock = grid->iblock;
  jblock = grid->jblock;
  kblock = grid->kblock;
  igc    = grid->igc;
  jgc    = grid->jgc;
  kgc    = grid->kgc;

  grid->fftforward(p, work3d, grid->fftini, grid->fftouti, grid->fftinj, grid->fftoutj);

  jj = iblock;
  kk = iblock*jblock;

  int ik;

  const int mpicoordx = master->mpicoordx;
  const int mpicoordy = master->mpicoordy;

  // Calculate the step size.
  const int nj = jblock/jslice;

  const int kki1 = 1*iblock*jslice;
  const int kki2 = 2*iblock*jslice;
  const int kki3 = 3*iblock*jslice;

  for(int n=0; n<nj; ++n)
  {
    for(int j=0; j<jslice; ++j)
#pragma ivdep
      for(int i=0; i<iblock; ++i)
      {
        // Set a zero gradient bc at the bottom.
        ik = i + j*jj;
        m1temp[ik] =  0.;
        m2temp[ik] =  0.;
        m3temp[ik] =  0.;
        m4temp[ik] =  1.;
        m5temp[ik] =  0.;
        m6temp[ik] =  0.;
        m7temp[ik] = -1.;
        ptemp [ik] =  0.;
      }

    for(int j=0; j<jslice; ++j)
#pragma ivdep
      for(int i=0; i<iblock; ++i)
      {
        ik = i + j*jj;
        m1temp[ik+kki1] =  0.;
        m2temp[ik+kki1] =  0.;
        m3temp[ik+kki1] =  0.;
        m4temp[ik+kki1] =  1.;
        m5temp[ik+kki1] = -1.;
        m6temp[ik+kki1] =  0.;
        m7temp[ik+kki1] =  0.;
        ptemp [ik+kki1] =  0.;
      }

    for(int k=0; k<kmax; ++k)
      for(int j=0; j<jslice; ++j)
      {
        jindex = mpicoordx*jblock + n*jslice + j;
#pragma ivdep
        for(int i=0; i<iblock; ++i)
        {
          // Swap the mpicoords, because domain is turned 90 degrees to avoid two mpi transposes.
          iindex = mpicoordy*iblock + i;

          ijk = i + (j + n*jslice)*jj + k*kk;
          ik  = i + j*jj + k*kki1;
          m1temp[ik+kki2] = m1[k];
          m2temp[ik+kki2] = m2[k];
          m3temp[ik+kki2] = m3[k];
          m4temp[ik+kki2] = m4[k] + bmati[iindex] + bmatj[jindex];
          m5temp[ik+kki2] = m5[k];
          m6temp[ik+kki2] = m6[k];
          m7temp[ik+kki2] = m7[k];
          ptemp [ik+kki2] = p[ijk];
        }
      }

    for(int j=0; j<jslice; ++j)
    {
      jindex = mpicoordx*jblock + n*jslice + j;
#pragma ivdep
      for(int i=0; i<iblock; ++i)
      {
        // Swap the mpicoords, because domain is turned 90 degrees to avoid two mpi transposes.
        iindex = mpicoordy*iblock + i;

        // Set the top boundary.
        ik = i + j*jj + kmax*kki1;
        if(iindex == 0 && jindex == 0)
        {
          m1temp[ik+kki2] =    0.;
          m2temp[ik+kki2] = -1/3.;
          m3temp[ik+kki2] =    2.;
          m4temp[ik+kki2] =    1.;

          m1temp[ik+kki3] =   -2.;
          m2temp[ik+kki3] =    9.;
          m3temp[ik+kki3] =    0.;
          m4temp[ik+kki3] =    1.;
        }
        // Set dp/dz at top to zero.
        else
        {
          m1temp[ik+kki2] =  0.;
          m2temp[ik+kki2] =  0.;
          m3temp[ik+kki2] = -1.;
          m4temp[ik+kki2] =  1.;

          m1temp[ik+kki3] = -1.;
          m2temp[ik+kki3] =  0.;
          m3temp[ik+kki3] =  0.;
          m4temp[ik+kki3] =  1.;
        }
      }
    }

    for(int j=0; j<jslice; ++j)
#pragma ivdep
      for(int i=0; i<iblock; ++i)
      {
        // Set the top boundary.
        ik = i + j*jj + kmax*kki1;
        m5temp[ik+kki2] = 0.;
        m6temp[ik+kki2] = 0.;
        m7temp[ik+kki2] = 0.;
        ptemp [ik+kki2] = 0.;

        m5temp[ik+kki3] = 0.;
        m6temp[ik+kki3] = 0.;
        m7temp[ik+kki3] = 0.;
        ptemp [ik+kki3] = 0.;
      }

    hdma(m1temp, m2temp, m3temp, m4temp, m5temp, m6temp, m7temp, ptemp, jslice);

    // Put back the solution.
    for(int k=0; k<kmax; ++k)
      for(int j=0; j<jslice; ++j)
#pragma ivdep
        for(int i=0; i<iblock; ++i)
        {
          ik  = i + j*jj + k*kki1;
          ijk = i + (j + n*jslice)*jj + k*kk;
          p[ijk] = ptemp[ik+kki2];
        }
  }

  grid->fftbackward(p, work3d, grid->fftini, grid->fftouti, grid->fftinj, grid->fftoutj);

  // Put the pressure back onto the original grid including ghost cells.
  jj = imax;
  kk = imax*jmax;

  int ijkp,jjp,kkp1,kkp2;
  jjp = grid->icells;
  kkp1 = 1*grid->ijcells;
  kkp2 = 2*grid->ijcells;

  for(int k=0; k<grid->kmax; k++)
    for(int j=0; j<grid->jmax; j++)
#pragma ivdep
      for(int i=0; i<grid->imax; i++)
      {
        ijkp = i+igc + (j+jgc)*jjp + (k+kgc)*kkp1;
        ijk  = i + j*jj + k*kk;
        p[ijkp] = work3d[ijk];
      }

  // Set a zero gradient boundary at the bottom.
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jjp + grid->kstart*kkp1;
      p[ijk-kkp1] = p[ijk     ];
      p[ijk-kkp2] = p[ijk+kkp1];
    }

  // Set a zero gradient boundary at the top.
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jjp + (grid->kend-1)*kkp1;
      p[ijk+kkp1] = p[ijk     ];
      p[ijk+kkp2] = p[ijk-kkp1];
    }

  // Set the cyclic boundary conditions.
  grid->boundary_cyclic(p);
}

void Pres_4::pres_out(double * restrict ut, double * restrict vt, double * restrict wt, 
                       double * restrict p , double * restrict dzhi4)
{
  int    ijk,ii1,ii2,jj1,jj2,kk1,kk2;
  int    kstart;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  kstart = grid->kstart;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + kstart*kk1;
      ut[ijk] -= (cg0*p[ijk-ii2] + cg1*p[ijk-ii1] + cg2*p[ijk] + cg3*p[ijk+ii1]) * cgi*dxi;
      vt[ijk] -= (cg0*p[ijk-jj2] + cg1*p[ijk-jj1] + cg2*p[ijk] + cg3*p[ijk+jj1]) * cgi*dyi;
    }

  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        ut[ijk] -= (cg0*p[ijk-ii2] + cg1*p[ijk-ii1] + cg2*p[ijk] + cg3*p[ijk+ii1]) * cgi*dxi;
        vt[ijk] -= (cg0*p[ijk-jj2] + cg1*p[ijk-jj1] + cg2*p[ijk] + cg3*p[ijk+jj1]) * cgi*dyi;
        wt[ijk] -= (cg0*p[ijk-kk2] + cg1*p[ijk-kk1] + cg2*p[ijk] + cg3*p[ijk+kk1]) * dzhi4[k];
      }
}

void Pres_4::hdma(double * restrict m1, double * restrict m2, double * restrict m3, double * restrict m4,
                   double * restrict m5, double * restrict m6, double * restrict m7, double * restrict p,
                   const int jslice)
{
  const int kmax   = grid->kmax;
  const int iblock = grid->iblock;

  const int jj = grid->iblock;

  const int kk1 = 1*grid->iblock*jslice;
  const int kk2 = 2*grid->iblock*jslice;
  const int kk3 = 3*grid->iblock*jslice;

  int k,ik;

  // Use LU factorization.
  k = 0;
  for(int j=0; j<jslice; ++j)
#pragma ivdep
    for(int i=0; i<iblock; ++i)
    {
      ik = i + j*jj;
      m1[ik] = 1.;
      m2[ik] = 1.;
      m3[ik] = 1.            / m4[ik];
      m4[ik] = 1.;
      m5[ik] = m5[ik]*m3[ik];
      m6[ik] = m6[ik]*m3[ik];
      m7[ik] = m7[ik]*m3[ik];
    }

  k = 1;
  for(int j=0; j<jslice; ++j)
#pragma ivdep
    for(int i=0; i<iblock; ++i)
    {
      ik = i + j*jj + k*kk1;
      m1[ik] = 1.;
      m2[ik] = 1.;
      m3[ik] = m3[ik]                     / m4[ik-kk1];
      m4[ik] = m4[ik] - m3[ik]*m5[ik-kk1];
      m5[ik] = m5[ik] - m3[ik]*m6[ik-kk1];
      m6[ik] = m6[ik] - m3[ik]*m7[ik-kk1];
    }

  k = 2;
  for(int j=0; j<jslice; ++j)
#pragma ivdep
    for(int i=0; i<iblock; ++i)
    {
      ik = i + j*jj + k*kk1;
      m1[ik] = 1.;
      m2[ik] =   m2[ik]                                           / m4[ik-kk2];
      m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] ) / m4[ik-kk1];
      m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2];
      m5[ik] =   m5[ik] - m3[ik]*m6[ik-kk1] - m2[ik]*m7[ik-kk2];
      m6[ik] =   m6[ik] - m3[ik]*m7[ik-kk1];
    }

  for(k=3; k<kmax+2; ++k)
    for(int j=0; j<jslice; ++j)
#pragma ivdep
      for(int i=0; i<iblock; ++i)
      {
        ik = i + j*jj + k*kk1;
        m1[ik] = ( m1[ik]                                                            ) / m4[ik-kk3];
        m2[ik] = ( m2[ik]                                         - m1[ik]*m5[ik-kk3]) / m4[ik-kk2];
        m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] - m1[ik]*m6[ik-kk3]) / m4[ik-kk1];
        m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2] - m1[ik]*m7[ik-kk3];
        m5[ik] =   m5[ik] - m3[ik]*m6[ik-kk1] - m2[ik]*m7[ik-kk2];
        m6[ik] =   m6[ik] - m3[ik]*m7[ik-kk1];
      }

  k = kmax+1;
  for(int j=0; j<jslice; ++j)
#pragma ivdep
    for(int i=0; i<iblock; ++i)
    {
      ik = i + j*jj + k*kk1;
      m7[ik] = 1.;
    }

  k = kmax+2;
  for(int j=0; j<jslice; ++j)
#pragma ivdep
    for(int i=0; i<iblock; ++i)
    {
      ik = i + j*jj + k*kk1;
      m1[ik] = ( m1[ik]                                                            ) / m4[ik-kk3];
      m2[ik] = ( m2[ik]                                         - m1[ik]*m5[ik-kk3]) / m4[ik-kk2];
      m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] - m1[ik]*m6[ik-kk3]) / m4[ik-kk1];
      m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2] - m1[ik]*m7[ik-kk3];
      m5[ik] =   m5[ik] - m3[ik]*m6[ik-kk1] - m2[ik]*m7[ik-kk2];
      m6[ik] = 1.;
      m7[ik] = 1.;
    }

  k = kmax+3;
  for(int j=0; j<jslice; ++j)
#pragma ivdep
    for(int i=0; i<iblock; ++i)
    {
      ik = i + j*jj + k*kk1;
      m1[ik] = ( m1[ik]                                                            ) / m4[ik-kk3];
      m2[ik] = ( m2[ik]                                         - m1[ik]*m5[ik-kk3]) / m4[ik-kk2];
      m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] - m1[ik]*m6[ik-kk3]) / m4[ik-kk1];
      m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2] - m1[ik]*m7[ik-kk3];
      m5[ik] = 1.;
      m6[ik] = 1.;
      m7[ik] = 1.;
    }

  // Do the backward substitution.
  // First, solve Ly = p, forward.
  for(int j=0; j<jslice; ++j)
#pragma ivdep
    for(int i=0; i<iblock; ++i)
    {
      ik = i + j*jj;
      p[ik    ] =             p[ik    ]*m3[ik    ];
      p[ik+kk1] = p[ik+kk1] - p[ik    ]*m3[ik+kk1];
      p[ik+kk2] = p[ik+kk2] - p[ik+kk1]*m3[ik+kk2] - p[ik]*m2[ik+kk2];
    }

  for(k=3; k<kmax+4; ++k)
    for(int j=0; j<jslice; ++j)
#pragma ivdep
      for(int i=0; i<iblock; ++i)
      {
        ik = i + j*jj + k*kk1;
        p[ik] = p[ik] - p[ik-kk1]*m3[ik] - p[ik-kk2]*m2[ik] - p[ik-kk3]*m1[ik];
      }

  // Second, solve Ux=y, backward.
  k = kmax+3;
  for(int j=0; j<jslice; ++j)
#pragma ivdep
    for(int i=0; i<iblock; ++i)
    {
      ik = i + j*jj + k*kk1;
      p[ik    ] =   p[ik    ]                                             / m4[ik    ];
      p[ik-kk1] = ( p[ik-kk1] - p[ik    ]*m5[ik-kk1] )                    / m4[ik-kk1];
      p[ik-kk2] = ( p[ik-kk2] - p[ik-kk1]*m5[ik-kk2] - p[ik]*m6[ik-kk2] ) / m4[ik-kk2];
    }

  for(k=kmax; k>=0; --k)
    for(int j=0; j<jslice; ++j)
#pragma ivdep
      for(int i=0; i<iblock; ++i)
      {
        ik = i + j*jj + k*kk1;
        p[ik] = ( p[ik] - p[ik+kk1]*m5[ik] - p[ik+kk2]*m6[ik] - p[ik+kk3]*m7[ik] ) / m4[ik];
      }
}

double Pres_4::calcdivergence(double * restrict u, double * restrict v, double * restrict w, double * restrict dzi4)
{
  int    ijk,ii1,ii2,jj1,jj2,kk1,kk2;
  int    kstart,kend;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk1 = 1*grid->ijcells;
  kk2 = 2*grid->ijcells;

  kstart = grid->kstart;
  kend   = grid->kend;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  double div, divmax;
  divmax = 0;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        div = (cg0*u[ijk-ii1] + cg1*u[ijk] + cg2*u[ijk+ii1] + cg3*u[ijk+ii2]) * cgi*dxi
            + (cg0*v[ijk-jj1] + cg1*v[ijk] + cg2*v[ijk+jj1] + cg3*v[ijk+jj2]) * cgi*dyi
            + (cg0*w[ijk-kk1] + cg1*w[ijk] + cg2*w[ijk+kk1] + cg3*w[ijk+kk2]) * dzi4[k];

        divmax = std::max(divmax, std::abs(div));
      }

  grid->getmax(&divmax);

  return divmax;
}
