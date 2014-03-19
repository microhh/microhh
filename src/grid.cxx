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
#include "master.h"
#include "grid.h"
#include "input.h"
#include "defines.h"
#include "model.h"

/**
 * This function constructs the grid class.
 * @param modelin Pointer to the model class.
 */
cgrid::cgrid(cmodel *modelin)
{
  master = modelin->master;

  allocated = false;
  mpitypes  = false;
  fftwplan  = false;
}

/**
 * This function destructs the grid class.
 */
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

/**
 * This function processes the input data and stores them in the class
 * variables.
 * @param inputin Pointer to the input class.
 * @return Returns 1 on error, 0 otherwise.
 */
int cgrid::readinifile(cinput *inputin)
{
  int nerror = 0;

  nerror += inputin->getItem(&xsize, "grid", "xsize", "");
  nerror += inputin->getItem(&ysize, "grid", "ysize", "");
  nerror += inputin->getItem(&zsize, "grid", "zsize", "");

  nerror += inputin->getItem(&itot, "grid", "itot", "");
  nerror += inputin->getItem(&jtot, "grid", "jtot", "");
  nerror += inputin->getItem(&ktot, "grid", "ktot", "");

  // velocity of the grid for gaelian transformation
  nerror += inputin->getItem(&utrans, "grid", "utrans", "", 0.);
  nerror += inputin->getItem(&vtrans, "grid", "vtrans", "", 0.);

  nerror += inputin->getItem(&swspatialorder, "grid", "swspatialorder", "");

  if(nerror > 0)
    return nerror;

  if(!(swspatialorder == "2" || swspatialorder == "4"))
  {
    if(master->mpiid == 0) std::printf("ERROR \"%s\" is an illegal value for swspatialorder\n", swspatialorder.c_str());
    return 1;
  }
 
  // 2nd order scheme requires only 1 ghost cell
  if(swspatialorder == "2")
  {
    igc = 1;
    jgc = 1;
    kgc = 1;
  }
  // 4th order scheme requires 3 ghost cells
  else if(swspatialorder == "4")
  {
    igc = 3;
    jgc = 3;
    kgc = 3;
  }

  return 0;
}

/**
 * This function allocates the dynamic arrays in the field class
 * variables and calculates the derived grid indices and dimensions.
 * @return Returns 1 on error, 0 otherwise.
 */
int cgrid::init()
{
  // check whether the grid fits the processor configuration
  if(itot % master->npx != 0)
  {
    if(master->mpiid == 0) std::printf("ERROR itot = %d is not a multiple of npx = %d\n", itot, master->npx);
    return 1;
  }
  if(itot % master->npy != 0)
  {
    if(master->mpiid == 0) std::printf("ERROR itot = %d is not a multiple of npy = %d\n", itot, master->npy);
    return 1;
  }
  // check this one only when npy > 1, since the transpose in that direction only happens then
  if(jtot % master->npx != 0)
  {
    if(master->mpiid == 0) std::printf("ERROR jtot = %d is not a multiple of npx = %d\n", jtot, master->npx);
    return 1;
  }
  if(jtot % master->npy != 0)
  {
    if(master->mpiid == 0) std::printf("ERROR jtot = %d is not a multiple of npy = %d\n", jtot, master->npy);
    return 1;
  }
  if(ktot % master->npx != 0)
  {
    if(master->mpiid == 0) std::printf("ERROR ktot = %d is not a multiple of npx = %d\n", ktot, master->npx);
    return 1;
  }

  // calculate the grid dimensions per process
  imax   = itot / master->npx;
  jmax   = jtot / master->npy;
  kmax   = ktot;

  // calculate the block sizes for the transposes
  iblock = itot / master->npy;
  jblock = jtot / master->npx;
  kblock = ktot / master->npx;

  // calculate the grid dimensions including ghost cells
  icells  = (imax+2*igc);
  jcells  = (jmax+2*jgc);
  ijcells = icells*jcells;
  kcells  = (kmax+2*kgc);
  ncells  = (imax+2*igc)*(jmax+2*jgc)*(kmax+2*kgc);

  // calculate the starting and ending points for loops over the grid
  istart = igc;
  jstart = jgc;
  kstart = kgc;

  iend   = imax + igc;
  jend   = jmax + jgc;
  kend   = kmax + kgc;

  // allocate all arrays
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
  fftini  = fftw_alloc_real(itot*jmax);
  fftouti = fftw_alloc_real(itot*jmax);
  fftinj  = fftw_alloc_real(jtot*iblock);
  fftoutj = fftw_alloc_real(jtot*iblock);

  allocated = true;

  // initialize the communication functions
  initmpi();

  return 0;
}

/**
 * This function initializes the fields containing the grid dimensions based
 * on the profiles in the input file.
 * @param inputin Pointer to the input class.
 * @return Returns 1 on error, 0 otherwise.
 */
int cgrid::create(cinput *inputin)
{
  // get the grid coordinates from the input
  if(inputin->getProf(&z[kstart], "z", kmax))
    return 1;

  // calculate the grid
  calculate();

  return 0;
}

/**
 * This function calculates the scalars and arrays that contain the information
 * on the grid spacing.
 * @return Returns 0.
 */
int cgrid::calculate()
{
  int i,j,k;

  // calculate the grid spacing
  dx = xsize / itot;
  dy = ysize / jtot;

  // calculate the offset per process to get the true x- and y-coordinate
  double xoff = master->mpicoordx * xsize / master->npx;
  double yoff = master->mpicoordy * ysize / master->npy;

  // calculate the x and y coordinates
  for(i=0; i<icells; ++i)
  {
    x [i] = 0.5*dx + (i-igc)*dx + xoff;
    xh[i] = (i-igc)*dx + xoff;
  }

  for(j=0; j<jcells; ++j)
  {
    y [j] = 0.5*dy + (j-jgc)*dy + yoff;
    yh[j] = (j-jgc)*dy + yoff;
  }

  // the calculation of ghost cells and flux levels has to go according to numerical scheme
  if(swspatialorder == "2")
  {
    z[kstart-1] = -z[kstart];
    z[kend]     = 2.*zsize - z[kend-1];

    for(k=kstart+1; k<kend; ++k)
      zh[k] = 0.5*(z[k-1]+z[k]);
    zh[kstart] = 0.;
    zh[kend]   = zsize;

    // calculate the half levels according to the numerical scheme
    // compute the height of the grid cells
    for(k=1; k<kcells; ++k)
    {
      dzh [k] = z[k] - z[k-1];
      dzhi[k] = 1./dzh[k];
    }
    dzh [kstart-1] = dzh [kstart+1];
    dzhi[kstart-1] = dzhi[kstart+1];

    // compute the height of the grid cells
    for(k=1; k<kcells-1; ++k)
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
    for(k=kstart+1; k<kend; ++k)
      zh[k] = ci0*z[k-2] + ci1*z[k-1] + ci2*z[k] + ci3*z[k+1];
    zh[kend] = zsize;

    zh[kstart-1] = bi0*z[kstart-2] + bi1*z[kstart-1] + bi2*z[kstart] + bi3*z[kstart+1];
    zh[kend+1]   = ti0*z[kend-2  ] + ti1*z[kend-1  ] + ti2*z[kend  ] + ti3*z[kend+1  ];

    // calculate the half levels according to the numerical scheme
    // compute the height of the grid cells
    for(k=1; k<kcells; ++k)
    {
      dzh [k] = z[k] - z[k-1];
      dzhi[k] = 1./dzh[k];
    }
    dzh [kstart-3] = dzh [kstart+3];
    dzhi[kstart-3] = dzhi[kstart+3];

    // compute the height of the grid cells
    for(k=1; k<kcells-1; ++k)
    {
      dz [k] = zh[k+1] - zh[k];
      dzi[k] = 1./dz[k];
    }
    dz [kstart-3] = dz [kstart+2];
    dzi[kstart-3] = dzi[kstart+2];
    dz [kend+2] = dz [kend-3];
    dzi[kend+2] = dzi[kend-3];

    // calculate the fourth order gradients
    for(k=kstart; k<kend; ++k)
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

/**
 * This function does a second order horizontal interpolation in the x-direction
 * to the selected location on the grid.
 * @param out Pointer to the output field.
 * @param in Pointer to the input field.
 * @param locx Integer containing the location of the input field,
 * where a value of 1 refers to the flux level.
 * @return Returns 0.
 */
int cgrid::interpolate_2nd(double * restrict out, double * restrict in, const int locin[3], const int locout[3])
{
  int ijk,ii,jj,kk,iih,jjh,kkh;

  ii = 1;
  jj = icells;
  kk = ijcells;

  iih = (locin[0]-locout[0])*ii;
  jjh = (locin[1]-locout[1])*jj;
  kkh = (locin[2]-locout[2])*kk;

  // interpolate the field
  // \TODO add the vertical component
  for(int k=kstart; k<kend+locout[2]; ++k)
    for(int j=jstart; j<jend; ++j)
#pragma ivdep
      for(int i=istart; i<iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        out[ijk] = 0.5*(0.5*in[ijk    ] + 0.5*in[ijk+iih    ])
                 + 0.5*(0.5*in[ijk+jjh] + 0.5*in[ijk+iih+jjh]);
      }
}

/**
 * This function does a fourth order horizontal interpolation in the x-direction
 * to the selected location on the grid.
 * @param out Pointer to the output field.
 * @param in Pointer to the input field.
 * @param locx Integer containing the location of the input field,
 * where a value of 1 refers to the flux level.
 * @return Returns 0.
 */
int cgrid::interpolate_4th(double * restrict out, double * restrict in, const int locin[3], const int locout[3])
{
  // interpolation function, locx = 1 indicates that the reference is at the half level
  int ijk,ii,jj,kk,iih1,jjh1,kkh1,iih2,jjh2;

  ii = 1;
  jj = icells;
  kk = ijcells;

  // a shift to the left gives minus 1 a shift to the right +1
  iih1 = 1*(locin[0]-locout[0])*ii;
  iih2 = 2*(locin[0]-locout[0])*ii;
  jjh1 = 1*(locin[1]-locout[1])*jj;
  jjh2 = 2*(locin[1]-locout[1])*jj;
  kkh1 = 1*(locin[2]-locout[2])*kk;

  // \TODO add the vertical component
  for(int k=kstart; k<kend+locout[2]; ++k)
    for(int j=jstart; j<jend; ++j)
#pragma ivdep
      for(int i=istart; i<iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        out[ijk] = ci0*(ci0*in[ijk-iih1-jjh1] + ci1*in[ijk-jjh1] + ci2*in[ijk+iih1-jjh1] + ci3*in[ijk+iih2-jjh1])
                 + ci1*(ci0*in[ijk-iih1     ] + ci1*in[ijk     ] + ci2*in[ijk+iih1     ] + ci3*in[ijk+iih2     ])
                 + ci2*(ci0*in[ijk-iih1+jjh1] + ci1*in[ijk+jjh1] + ci2*in[ijk+iih1+jjh1] + ci3*in[ijk+iih2+jjh1])
                 + ci3*(ci0*in[ijk-iih1+jjh2] + ci1*in[ijk+jjh2] + ci2*in[ijk+iih1+jjh2] + ci3*in[ijk+iih2+jjh2]);
      }

  return 0;
}

/**
 * This function does a fourth order horizontal interpolation in the y-direction
 * to the selected location on the grid.
 * @param prof Pointer to the output field.
 * @param data Pointer to the input field.
 * @param krange Number of vertical levels over which the profile is to be calculated.
 * @return Returns 0.
 */
int cgrid::calcmean(double * restrict prof, const double * restrict data, const int krange)
{
  int ijk,ii,jj,kk;

  ii = 1;
  jj = icells;
  kk = ijcells;
  
  for(int k=0; k<krange; ++k)
  {
    prof[k] = 0.;
    for(int j=jstart; j<jend; ++j)
#pragma ivdep
      for(int i=istart; i<iend; ++i)
      {
        ijk  = i + j*jj + k*kk;
        prof[k] += data[ijk];
      }
  }

  double n = imax*jmax;

  for(int k=0; k<krange; ++k)
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
  for(int k=kstart; k<kend+(1-locz); ++k)
    for(int j=jstart; j<jend; ++j)
#pragma ivdep
      for(int i=istart; i<iend; ++i)
      {
        ijk = i + j*jj1 + k*kk1;
        out[ijk] = ci0*in[ijk-kk2+khlf] + ci1*in[ijk-kk1+khlf] + ci2*in[ijk+khlf] + ci3*in[ijk+kk1+khlf];
      }

  return 0;
}
*/
// end of interpolation functions

