/*
 * MicroHH
 * Copyright (c) 2011-2013 Chiel van Heerwaarden
 * Copyright (c) 2011-2013 Thijs Heus
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
#include "grid.h"
#include "input.h"
#include "defines.h"

// build the grid
cgrid::cgrid(cmpi *mpiin)
{
  mpi = mpiin;

  allocated = false;
  mpitypes  = false;
  fftwplan  = false;
}

cgrid::~cgrid()
{
  if(fftwplan)
  {
    fftw_destroy_plan(iplanf);
    fftw_destroy_plan(iplanb);
    fftw_destroy_plan(jplanf);
    fftw_destroy_plan(jplanb);
  }

  if(allocated)
  { 
    delete[] x;
    delete[] xh;
    delete[] y;
    delete[] yh;
    delete[] z;
    delete[] zh;
    delete[] dz;
    delete[] dzh;
    delete[] dzi;
    delete[] dzhi;
    delete[] dzi4;
    delete[] dzhi4;

    fftw_free(fftini);
    fftw_free(fftouti);
    fftw_free(fftinj);
    fftw_free(fftoutj);

    fftw_cleanup();
  }

  exitmpi();
}

int cgrid::readinifile(cinput *inputin)
{
  int n = 0;

  n += inputin->getItem(&xsize, "grid", "xsize", "");
  n += inputin->getItem(&ysize, "grid", "ysize", "");
  n += inputin->getItem(&zsize, "grid", "zsize", "");

  n += inputin->getItem(&itot, "grid", "itot", "");
  n += inputin->getItem(&jtot, "grid", "jtot", "");
  n += inputin->getItem(&ktot, "grid", "ktot", "");

  // velocity of the grid for gaelian transformation
  n += inputin->getItem(&u, "grid", "u", "", 0.);
  n += inputin->getItem(&v, "grid", "v", "", 0.);

  n += inputin->getItem(&swspatialorder, "grid", "swspatialorder", "");

  if(n > 0)
    return 1;

  if(!(swspatialorder == "2" || swspatialorder == "4"))
  {
    if(mpi->mpiid == 0) std::printf("ERROR \"%s\" is an illegal value for swspatialorder\n", swspatialorder.c_str());
    return 1;
  }
 
  // 2nd order scheme requires only 1 ghost cell, a 4th order requires 2
  if(swspatialorder == "2")
  {
    igc = 1;
    jgc = 1;
    kgc = 1;
  }
  else if(swspatialorder == "4")
  {
    igc = 3;
    jgc = 3;
    kgc = 3;
  }

  return 0;
}

int cgrid::init()
{
  // check whether the grid fits the processor configuration
  if(itot % mpi->npx != 0)
  {
    if(mpi->mpiid == 0) std::printf("ERROR itot = %d is not a multiple of npx = %d\n", itot, mpi->npx);
    return 1;
  }
  if(itot % mpi->npy != 0)
  {
    if(mpi->mpiid == 0) std::printf("ERROR itot = %d is not a multiple of npy = %d\n", itot, mpi->npy);
    return 1;
  }
  // check this one only when npy > 1, since the transpose in that direction only happens then
  //if(jtot % mpi->npx != 0 && mpi->npy > 1)
  if(jtot % mpi->npx != 0)
  {
    if(mpi->mpiid == 0) std::printf("ERROR jtot = %d is not a multiple of npx = %d\n", jtot, mpi->npx);
    return 1;
  }
  if(jtot % mpi->npy != 0)
  {
    if(mpi->mpiid == 0) std::printf("ERROR jtot = %d is not a multiple of npy = %d\n", jtot, mpi->npy);
    return 1;
  }
  if(ktot % mpi->npx != 0)
  {
    if(mpi->mpiid == 0) std::printf("ERROR ktot = %d is not a multiple of npx = %d\n", ktot, mpi->npx);
    return 1;
  }

  imax   = itot / mpi->npx;
  jmax   = jtot / mpi->npy;
  kmax   = ktot;

  iblock = itot / mpi->npy;
  jblock = jtot / mpi->npx;
  kblock = ktot / mpi->npx;

  icells = (imax+2*igc);
  jcells = (jmax+2*jgc);
  kcells = (kmax+2*kgc);
  ncells = (imax+2*igc)*(jmax+2*jgc)*(kmax+2*kgc);

  istart = igc;
  jstart = jgc;
  kstart = kgc;

  iend   = imax + igc;
  jend   = jmax + jgc;
  kend   = kmax + kgc;

  x     = new double[imax+2*igc];
  xh    = new double[imax+2*igc];
  y     = new double[jmax+2*jgc];
  yh    = new double[jmax+2*jgc];
  z     = new double[kmax+2*kgc];
  zh    = new double[kmax+2*kgc];
  dz    = new double[kmax+2*kgc];
  dzh   = new double[kmax+2*kgc];
  dzi   = new double[kmax+2*kgc];
  dzhi  = new double[kmax+2*kgc];
  dzi4  = new double[kmax+2*kgc];
  dzhi4 = new double[kmax+2*kgc];

  // allocate the data for the fourier transforms
  fftini  = fftw_alloc_real(itot);
  fftouti = fftw_alloc_real(itot);
  fftinj  = fftw_alloc_real(jtot);
  fftoutj = fftw_alloc_real(jtot);

  allocated = true;

  // initialize the communication functions
  initmpi();

  return 0;
}

int cgrid::create(cinput *inputin)
{
  if(inputin->getProf(&z[kstart], "z", kmax))
    return 1;

  calculate();

  return 0;
}

int cgrid::calculate()
{
  int i,j,k;

  dx = xsize / itot;
  dy = ysize / jtot;

  double xoff = mpi->mpicoordx * xsize / mpi->npx;
  double yoff = mpi->mpicoordy * ysize / mpi->npy;

  // calculate the x and y coordinates
  for(i=0; i<icells; i++)
  {
    x [i] = 0.5*dx + (i-igc)*dx + xoff;
    xh[i] = (i-igc)*dx + xoff;
  }

  for(j=0; j<jcells; j++)
  {
    y [j] = 0.5*dy + (j-jgc)*dy + yoff;
    yh[j] = (j-jgc)*dy + yoff;
  }

  // the calculation of ghost cells and flux levels has to go according to numerical scheme
  if(swspatialorder == "2")
  {
    z[kstart-1] = -z[kstart];
    z[kend]     = 2.*zsize - z[kend-1];

    for(k=kstart+1; k<kend; k++)
      zh[k] = 0.5*(z[k-1]+z[k]);
    zh[kstart] = 0.;
    zh[kend]   = zsize;

    // calculate the half levels according to the numerical scheme
    // compute the height of the grid cells
    for(k=1; k<kcells; k++)
    {
      dzh [k] = z[k] - z[k-1];
      dzhi[k] = 1./dzh[k];
    }
    dzh [kstart-1] = dzh [kstart+1];
    dzhi[kstart-1] = dzhi[kstart+1];

    // compute the height of the grid cells
    for(k=1; k<kcells-1; k++)
    {
      dz [k] = zh[k+1] - zh[k];
      dzi[k] = 1./dz[k];
    }
    dz [kstart-1] = dz [kstart];
    dzi[kstart-1] = dzi[kstart];
    dz [kend]     = dz [kend-1];
    dzi[kend]     = dzi[kend-1];

    // do not calculate 4th order gradients for 2nd order
  }

  if(swspatialorder == "4")
  {
    // calculate the height of the ghost cell
    z[kstart-1] = -2.*z[kstart] + (1./3.)*z[kstart+1];
    z[kstart-2] = -9.*z[kstart] +      2.*z[kstart+1];

    z[kend  ] = (8./3.)*zsize - 2.*z[kend-1] + (1./3.)*z[kend-2];
    z[kend+1] =      8.*zsize - 9.*z[kend-1] +      2.*z[kend-2];

    zh[kstart  ] = 0.;
    for(k=kstart+1; k<kend; k++)
      zh[k] = ci0*z[k-2] + ci1*z[k-1] + ci2*z[k] + ci3*z[k+1];
    zh[kend] = zsize;

    zh[kstart-1] = bi0*z[kstart-2] + bi1*z[kstart-1] + bi2*z[kstart] + bi3*z[kstart+1];
    zh[kend+1]   = ti0*z[kend-2  ] + ti1*z[kend-1  ] + ti2*z[kend  ] + ti3*z[kend+1  ];

    // calculate the half levels according to the numerical scheme
    // compute the height of the grid cells
    for(k=1; k<kcells; k++)
    {
      dzh [k] = z[k] - z[k-1];
      dzhi[k] = 1./dzh[k];
    }
    dzh [kstart-3] = dzh [kstart+3];
    dzhi[kstart-3] = dzhi[kstart+3];

    // compute the height of the grid cells
    for(k=1; k<kcells-1; k++)
    {
      dz [k] = zh[k+1] - zh[k];
      dzi[k] = 1./dz[k];
    }
    dz [kstart-3] = dz [kstart+2];
    dzi[kstart-3] = dzi[kstart+2];
    dz [kend+2] = dz [kend-3];
    dzi[kend+2] = dzi[kend-3];

    // calculate the fourth order gradients
    for(k=kstart; k<kend; k++)
    {
      dzi4 [k] = 1./(cg0*zh[k-1] + cg1*zh[k  ] + cg2*zh[k+1] + cg3*zh[k+2]);
      dzhi4[k] = 1./(cg0*z [k-2] + cg1*z [k-1] + cg2*z [k  ] + cg3*z [k+1]);
    }
    dzhi4[kend  ] = 1./(cg0*z[kend-2] + cg1*z[kend-1] + cg2*z[kend] + cg3*z[kend+1]);

    // bc's
    dzi4 [kstart-1] = 1./(bg0*zh[kstart-1] + bg1*zh[kstart  ] + bg2*zh[kstart+1] + bg3*zh[kstart+2]);
    dzhi4[kstart-1] = 1./(bg0*z [kstart-2] + bg1*z [kstart-1] + bg2*z [kstart  ] + bg3*z [kstart+1]);

    dzi4 [kend  ] = 1./(tg0*zh[kend-2] + tg1*zh[kend-1] + tg2*zh[kend] + tg3*zh[kend+1]);
    dzhi4[kend+1] = 1./(tg0*z [kend-2] + tg1*z [kend-1] + tg2*z [kend] + tg3*z [kend+1]);
  }

  return 0;
}

// interpolation functions
// CvH merge interpolate functions later to something more consise but still vectorizable
int cgrid::interpolatex_2nd(double * restrict out, double * restrict in, int locx)
{
  // interpolation function, locx = 1 indicates that the reference is at the half level
  int ijk,ii,jj,kk,ihlf;

  ii = 1;
  jj = icells;
  kk = icells*jcells;

  ihlf = locx*ii;

  // interpolate in x
  for(int k=0; k<kcells; k++)
    for(int j=jstart; j<jend; j++)
#pragma ivdep
      for(int i=istart; i<iend; i++)
      {
        ijk = i + j*jj + k*kk;
        out[ijk] = 0.5*(in[ijk-ii+ihlf] + in[ijk+ihlf]);
      }

  return 0;
}

int cgrid::interpolatey_2nd(double * restrict out, double * restrict in, int locy)
{
  // interpolation function, locy = 1 indicates that the reference is at the half level
  int ijk,ii,jj,kk,jhlf;

  ii = 1;
  jj = icells;
  kk = icells*jcells;

  jhlf = locy*jj;

  // interpolate in y
  for(int k=0; k<kcells; k++)
    for(int j=jstart; j<jend; j++)
#pragma ivdep
      for(int i=istart; i<iend; i++)
      {
        ijk = i + j*jj + k*kk;
        out[ijk] = 0.5*(in[ijk-jj+jhlf] + in[ijk+jhlf]);
      }

  return 0;
}

// interpolation functions
// CvH merge interpolate functions later to something more consise but still vectorizable
int cgrid::interpolatex_4th(double * restrict out, double * restrict in, int locx)
{
  // interpolation function, locx = 1 indicates that the reference is at the half level
  int ijk,ii1,ii2,jj1,jj2,kk1,kk2,ihlf;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*icells;
  jj2 = 2*icells;
  kk1 = 1*icells*jcells;
  kk2 = 2*icells*jcells;

  ihlf = locx*ii1;

  // interpolate in x
  for(int k=0; k<kcells; k++)
    for(int j=jstart; j<jend; j++)
#pragma ivdep
      for(int i=istart; i<iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        out[ijk] = ci0*in[ijk-ii2+ihlf] + ci1*in[ijk-ii1+ihlf] + ci2*in[ijk+ihlf] + ci3*in[ijk+ii1+ihlf];
      }

  return 0;
}

int cgrid::interpolatey_4th(double * restrict out, double * restrict in, int locy)
{
  // interpolation function, locy = 1 indicates that the reference is at the half level
  int ijk,ii1,ii2,jj1,jj2,kk1,kk2,jhlf;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*icells;
  jj2 = 2*icells;
  kk1 = 1*icells*jcells;
  kk2 = 2*icells*jcells;

  jhlf = locy*jj1;

  // interpolate in y
  for(int k=0; k<kcells; k++)
    for(int j=jstart; j<jend; j++)
#pragma ivdep
      for(int i=istart; i<iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        out[ijk] = ci0*in[ijk-jj2+jhlf] + ci1*in[ijk-jj1+jhlf] + ci2*in[ijk+jhlf] + ci3*in[ijk+jj1+jhlf];
      }

  return 0;
}

int cgrid::calcmean(double * restrict prof, const double * restrict data, const int krange)
{
  int ijk,ii,jj,kk;

  ii = 1;
  jj = icells;
  kk = icells*jcells;
  
  for(int k=0; k<krange; k++)
  {
    prof[k] = 0.;
    for(int j=jstart; j<jend; j++)
#pragma ivdep
      for(int i=istart; i<iend; i++)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] += data[ijk];
      }
  }

  double n = imax*jmax;

  for(int k=0; k<krange; k++)
    prof[k] /= n;

  getprof(prof, krange);

  return 0;
}

/*
int cgrid::interpolatez_4th(double * restrict out, double * restrict in, int locz)
{
  // interpolation function, locz = 1 indicates that the reference is at the half level
  int ijk,ii1,ii2,jj1,jj2,kk1,kk2,khlf;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*icells;
  jj2 = 2*icells;
  kk1 = 1*icells*jcells;
  kk2 = 2*icells*jcells;

  khlf = locz*kk1;

  // interpolate in x, add aditional level when interpolation goes to half levels to have values at wall
  for(int k=kstart; k<kend+(1-locz); k++)
    for(int j=jstart; j<jend; j++)
#pragma ivdep
      for(int i=istart; i<iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        out[ijk] = ci0*in[ijk-kk2+khlf] + ci1*in[ijk-kk1+khlf] + ci2*in[ijk+khlf] + ci3*in[ijk+kk1+khlf];
      }

  return 0;
}
*/
// end of interpolation functions

