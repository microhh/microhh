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

cpres_4::cpres_4(cmodel *modelin, cinput *inputin) : cpres(modelin, inputin)
{
  m1 = m1temp = 0;
  m2 = m2temp = 0;
  m3 = m3temp = 0;
  m4 = m4temp = 0;
  m5 = m5temp = 0;
  m6 = m6temp = 0;
  m7 = m7temp = 0;

  bmati = 0;
  bmatj = 0;

  ptemp = 0;
}

cpres_4::~cpres_4()
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

  // CvH temporary, remove later...
  delete[] m1temp;
  delete[] m2temp;
  delete[] m3temp;
  delete[] m4temp;
  delete[] m5temp;
  delete[] m6temp;
  delete[] m7temp;
  delete[] ptemp;
}

int cpres_4::exec(double dt)
{
  // create the input for the pressure solver
  pres_in(fields->sd["p"]->data,
          fields->u ->data, fields->v ->data, fields->w ->data,
          fields->ut->data, fields->vt->data, fields->wt->data, 
          grid->dzi4, dt);

  // solve the system
  pres_solve(fields->sd["p"]->data, fields->sd["tmp1"]->data, grid->dz,
             m1, m2, m3, m4,
             m5, m6, m7,
             m1temp, m2temp, m3temp, m4temp,
             m5temp, m6temp, m7temp, ptemp,
             bmati, bmatj);

  // get the pressure tendencies from the pressure field
  pres_out(fields->ut->data, fields->vt->data, fields->wt->data, 
           fields->sd["p"]->data, grid->dzhi4);

  return 0;
}

double cpres_4::check()
{
  double divmax = 0.;

  divmax = calcdivergence((*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi4);

  return divmax;
}

void cpres_4::init()
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

  // allocate help variables for the matrix solver
  m1 = new double[kmax];
  m2 = new double[kmax];
  m3 = new double[kmax];
  m4 = new double[kmax];
  m5 = new double[kmax];
  m6 = new double[kmax];
  m7 = new double[kmax];

  // CvH temporary, remove later...
  m1temp = new double[grid->iblock*(grid->kmax+4)];
  m2temp = new double[grid->iblock*(grid->kmax+4)];
  m3temp = new double[grid->iblock*(grid->kmax+4)];
  m4temp = new double[grid->iblock*(grid->kmax+4)];
  m5temp = new double[grid->iblock*(grid->kmax+4)];
  m6temp = new double[grid->iblock*(grid->kmax+4)];
  m7temp = new double[grid->iblock*(grid->kmax+4)];
  ptemp  = new double[grid->iblock*(grid->kmax+4)];
}

void cpres_4::setvalues()
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

int cpres_4::pres_in(double * restrict p, 
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

  // set the cyclic boundary conditions for the tendencies
  grid->boundary_cyclic(ut);
  grid->boundary_cyclic(vt);
  grid->boundary_cyclic(wt);

  // set the bc 
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

  return 0;
}

int cpres_4::pres_solve(double * restrict p, double * restrict work3d, double * restrict dz,
                        double * restrict m1, double * restrict m2, double * restrict m3, double * restrict m4,
                        double * restrict m5, double * restrict m6, double * restrict m7,
                        double * restrict m1temp, double * restrict m2temp, double * restrict m3temp, double * restrict m4temp,
                        double * restrict m5temp, double * restrict m6temp, double * restrict m7temp, double * restrict ptemp,
                        double * restrict bmati, double * restrict bmatj)
{
  int i,j,k,jj,kk,ijk;
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

  int ik,kki1,kki2,kki3;
  kki1 = 1*iblock;
  kki2 = 2*iblock;
  kki3 = 3*iblock;

  int mpicoordx = master->mpicoordx;
  int mpicoordy = master->mpicoordy;

  for(j=0; j<jblock; ++j)
  {
    jindex = mpicoordx * jblock + j;
 
#pragma ivdep
    for(i=0; i<iblock; ++i)
    {
      // set a zero gradient bc at the bottom
      ik = i;
      m1temp[ik] =  0.;
      m2temp[ik] =  0.;
      m3temp[ik] =  0.;
      m4temp[ik] =  1.;
      m5temp[ik] =  0.;
      m6temp[ik] =  0.;
      m7temp[ik] = -1.;
      ptemp [ik] =  0.;
    }

#pragma ivdep
    for(i=0; i<iblock; ++i)
    {
      ik = i;
      m1temp[ik+kki1] =  0.;
      m2temp[ik+kki1] =  0.;
      m3temp[ik+kki1] =  0.;
      m4temp[ik+kki1] =  1.;
      m5temp[ik+kki1] = -1.;
      m6temp[ik+kki1] =  0.;
      m7temp[ik+kki1] =  0.;
      ptemp [ik+kki1] =  0.;
    }

    for(k=0; k<kmax; ++k)
    {
#pragma ivdep
      for(i=0; i<iblock; ++i)
      {
        // swap the mpicoords, because domain is turned 90 degrees to avoid two mpi transposes
        iindex = mpicoordy * iblock + i;

        ijk = i + j*jj + k*kk;
        ik  = i + k*kki1;
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
        
#pragma ivdep
    for(i=0; i<iblock; ++i)
    {
      // set the top boundary
      ik = i + kmax*kki1;
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
      // set dp/dz at top to zero
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

#pragma ivdep
    for(i=0; i<iblock; ++i)
    {
      // set the top boundary
      ik = i + kmax*kki1;
      m5temp[ik+kki2] = 0.;
      m6temp[ik+kki2] = 0.;
      m7temp[ik+kki2] = 0.;
      ptemp [ik+kki2] = 0.;

      m5temp[ik+kki3] = 0.;
      m6temp[ik+kki3] = 0.;
      m7temp[ik+kki3] = 0.;
      ptemp [ik+kki3] = 0.;
    }

    hdma(m1temp, m2temp, m3temp, m4temp, m5temp, m6temp, m7temp, ptemp);

    // put back the solution
    for(k=0; k<kmax; ++k)
      for(int i=0; i<iblock; ++i)
      {
        ik  = i + k*kki1;
        ijk = i + j*jj + k*kk;
        p[ijk] = ptemp[ik+kki2];
      }
  }

  grid->fftbackward(p, work3d, grid->fftini, grid->fftouti, grid->fftinj, grid->fftoutj);

  // put the pressure back onto the original grid including ghost cells
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

  // set the boundary conditions
  // set a zero gradient boundary at the bottom
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jjp + grid->kstart*kkp1;
      p[ijk-kkp1] = p[ijk     ];
      p[ijk-kkp2] = p[ijk+kkp1];
    }

  // set a zero gradient boundary at the top
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jjp + (grid->kend-1)*kkp1;
      p[ijk+kkp1] = p[ijk     ];
      p[ijk+kkp2] = p[ijk-kkp1];
    }

  // set the cyclic boundary conditions
  grid->boundary_cyclic(p);

  return 0;
}

int cpres_4::pres_out(double * restrict ut, double * restrict vt, double * restrict wt, 
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

  return 0;
}

void cpres_4::hdma(double * restrict m1, double * restrict m2, double * restrict m3, double * restrict m4,
                   double * restrict m5, double * restrict m6, double * restrict m7,
                   double * restrict p)
{
  int k,ik;
  int kmax   = grid->kmax;
  int iblock = grid->iblock;
  int kk1 = 1*grid->iblock;
  int kk2 = 2*grid->iblock;
  int kk3 = 3*grid->iblock;

  // LU factorization
  k = 0;
#pragma ivdep
  for(int i=0; i<iblock; ++i)
  {
    ik = i;
    m1[ik] = 1.;
    m2[ik] = 1.;
    m3[ik] = 1.            / m4[ik];
    m4[ik] = 1.;
    m5[ik] = m5[ik]*m3[ik];
    m6[ik] = m6[ik]*m3[ik];
    m7[ik] = m7[ik]*m3[ik];
  }

  k = 1;
#pragma ivdep
  for(int i=0; i<iblock; ++i)
  {
    ik = i + k*kk1;
    m1[ik] = 1.;
    m2[ik] = 1.;
    m3[ik] = m3[ik]                     / m4[ik-kk1];
    m4[ik] = m4[ik] - m3[ik]*m5[ik-kk1];
    m5[ik] = m5[ik] - m3[ik]*m6[ik-kk1];
    m6[ik] = m6[ik] - m3[ik]*m7[ik-kk1];
  }

  k = 2;
#pragma ivdep
  for(int i=0; i<iblock; ++i)
  {
    ik = i + k*kk1;
    m1[ik] = 1.;
    m2[ik] =   m2[ik]                                           / m4[ik-kk2];
    m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] ) / m4[ik-kk1];
    m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2];
    m5[ik] =   m5[ik] - m3[ik]*m6[ik-kk1] - m2[ik]*m7[ik-kk2];
    m6[ik] =   m6[ik] - m3[ik]*m7[ik-kk1];
  }

  for(k=3; k<kmax+2; ++k)
#pragma ivdep
    for(int i=0; i<iblock; ++i)
    {
      ik = i + k*kk1;
      m1[ik] = ( m1[ik]                                                            ) / m4[ik-kk3];
      m2[ik] = ( m2[ik]                                         - m1[ik]*m5[ik-kk3]) / m4[ik-kk2];
      m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] - m1[ik]*m6[ik-kk3]) / m4[ik-kk1];
      m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2] - m1[ik]*m7[ik-kk3];
      m5[ik] =   m5[ik] - m3[ik]*m6[ik-kk1] - m2[ik]*m7[ik-kk2];
      m6[ik] =   m6[ik] - m3[ik]*m7[ik-kk1];
    }

  k = kmax+1;
#pragma ivdep
  for(int i=0; i<iblock; ++i)
  {
    ik = i + k*kk1;
    m7[ik] = 1.;
  }

  k = kmax+2;
#pragma ivdep
  for(int i=0; i<iblock; ++i)
  {
    ik = i + k*kk1;
    m1[ik] = ( m1[ik]                                                            ) / m4[ik-kk3];
    m2[ik] = ( m2[ik]                                         - m1[ik]*m5[ik-kk3]) / m4[ik-kk2];
    m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] - m1[ik]*m6[ik-kk3]) / m4[ik-kk1];
    m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2] - m1[ik]*m7[ik-kk3];
    m5[ik] =   m5[ik] - m3[ik]*m6[ik-kk1] - m2[ik]*m7[ik-kk2];
    m6[ik] = 1.;
    m7[ik] = 1.;
  }

  k = kmax+3;
#pragma ivdep
  for(int i=0; i<iblock; ++i)
  {
    ik = i + k*kk1;
    m1[ik] = ( m1[ik]                                                            ) / m4[ik-kk3];
    m2[ik] = ( m2[ik]                                         - m1[ik]*m5[ik-kk3]) / m4[ik-kk2];
    m3[ik] = ( m3[ik]                     - m2[ik]*m5[ik-kk2] - m1[ik]*m6[ik-kk3]) / m4[ik-kk1];
    m4[ik] =   m4[ik] - m3[ik]*m5[ik-kk1] - m2[ik]*m6[ik-kk2] - m1[ik]*m7[ik-kk3];
    m5[ik] = 1.;
    m6[ik] = 1.;
    m7[ik] = 1.;
  }

  // Backward substitution 
  // Solve Ly=p, forward
#pragma ivdep
  for(int i=0; i<iblock; ++i)
  {
    ik = i;
    p[ik    ] =             p[ik    ]*m3[ik    ];
    p[ik+kk1] = p[ik+kk1] - p[ik    ]*m3[ik+kk1];
    p[ik+kk2] = p[ik+kk2] - p[ik+kk1]*m3[ik+kk2] - p[ik]*m2[ik+kk2];
  }

  for(k=3; k<kmax+4; ++k)
#pragma ivdep
    for(int i=0; i<iblock; ++i)
    {
      ik = i + k*kk1;
      p[ik] = p[ik] - p[ik-kk1]*m3[ik] - p[ik-kk2]*m2[ik] - p[ik-kk3]*m1[ik];
    }

  // Solve Ux=y, backward
  k = kmax+3;
#pragma ivdep
  for(int i=0; i<iblock; ++i)
  {
    ik = i + k*kk1;
    p[ik    ] =   p[ik    ]                                             / m4[ik    ];
    p[ik-kk1] = ( p[ik-kk1] - p[ik    ]*m5[ik-kk1] )                    / m4[ik-kk1];
    p[ik-kk2] = ( p[ik-kk2] - p[ik-kk1]*m5[ik-kk2] - p[ik]*m6[ik-kk2] ) / m4[ik-kk2];
  }

  for(k=kmax; k>=0; --k)
#pragma ivdep
    for(int i=0; i<iblock; ++i)
    {
      ik = i + k*kk1;
      p[ik] = ( p[ik] - p[ik+kk1]*m5[ik] - p[ik+kk2]*m6[ik] - p[ik+kk3]*m7[ik] ) / m4[ik];
    }
}

double cpres_4::calcdivergence(double * restrict u, double * restrict v, double * restrict w, double * restrict dzi4)
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
