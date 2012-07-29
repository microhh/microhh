#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fftw3.h>
#include "grid.h"
#include "fields.h"
#include "pres_g4.h"
#include "defines.h"

cpres_g4::cpres_g4(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  std::printf("Creating instance of object pres_g4\n");
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;

  allocated = false;
}

cpres_g4::~cpres_g4()
{
  if(allocated)
  {
    fftw_destroy_plan(iplanf);
    fftw_destroy_plan(iplanb);
    fftw_destroy_plan(jplanf);
    fftw_destroy_plan(jplanb);

    fftw_free(fftini);
    fftw_free(fftouti);
    fftw_free(fftinj);
    fftw_free(fftoutj);

    delete[] m0;
    delete[] m1;
    delete[] m2;
    delete[] m3;
    delete[] m4;
    delete[] m5;
    delete[] m6;
    delete[] m7;
    delete[] m8;

    delete[] work2d;

    delete[] bmati;
    delete[] bmatj;

    // CvH temporary, remove later...
    delete[] m0temp;
    delete[] m1temp;
    delete[] m2temp;
    delete[] m3temp;
    delete[] m4temp;
    delete[] m5temp;
    delete[] m6temp;
    delete[] m7temp;
    delete[] m8temp;
    delete[] ptemp;
  }

  std::printf("Destroying instance of object pres_g4\n");
}

int cpres_g4::init()
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

  // allocate help variables for the matrix solver
  m0 = new double[kmax];
  m1 = new double[kmax];
  m2 = new double[kmax];
  m3 = new double[kmax];
  m4 = new double[kmax];
  m5 = new double[kmax];
  m6 = new double[kmax];
  m7 = new double[kmax];
  m8 = new double[kmax];

  // CvH temporary, remove later...
  m0temp = new double[kmax+2];
  m1temp = new double[kmax+2];
  m2temp = new double[kmax+2];
  m3temp = new double[kmax+2];
  m4temp = new double[kmax+2];
  m5temp = new double[kmax+2];
  m6temp = new double[kmax+2];
  m7temp = new double[kmax+2];
  m8temp = new double[kmax+2];
  ptemp  = new double[kmax+2];

  work2d = new double[imax*jmax];

  double *zh, *z;
  z  = grid->z;
  zh = grid->zh;

  int k,kc;
  // create vectors that go into the tridiagonal matrix solver
  
  // bottom boundary
  k  = 0;
  kc = kstart+k;
  m0[k] = 0.;
  m1[k] = 0.;
  m2[k] = 0.;
  m3[k] = (  529.*grad4xbiasbot(z[kc-1], z[kc], z[kc+1], z[kc+2]) +  21.*grad4x(z[kc-1], z[kc], z[kc+1], z[kc+2])                                                                                                 ) / grad4xbiasbot(zh[kc], zh[kc+1], zh[kc+2], zh[kc+3]);
  m4[k] = ( -483.*grad4xbiasbot(z[kc-1], z[kc], z[kc+1], z[kc+2]) - 567.*grad4x(z[kc-1], z[kc], z[kc+1], z[kc+2]) +  3.*grad4x(z[kc], z[kc+1], z[kc+2], z[kc+3])                                                  ) / grad4xbiasbot(zh[kc], zh[kc+1], zh[kc+2], zh[kc+3]);
  m5[k] = (  -69.*grad4xbiasbot(z[kc-1], z[kc], z[kc+1], z[kc+2]) + 567.*grad4x(z[kc-1], z[kc], z[kc+1], z[kc+2]) - 81.*grad4x(z[kc], z[kc+1], z[kc+2], z[kc+3]) -  1.*grad4x(z[kc+1], z[kc+2], z[kc+3], z[kc+4]) ) / grad4xbiasbot(zh[kc], zh[kc+1], zh[kc+2], zh[kc+3]);
  m6[k] = (   23.*grad4xbiasbot(z[kc-1], z[kc], z[kc+1], z[kc+2]) -  21.*grad4x(z[kc-1], z[kc], z[kc+1], z[kc+2]) + 81.*grad4x(z[kc], z[kc+1], z[kc+2], z[kc+3]) + 27.*grad4x(z[kc+1], z[kc+2], z[kc+3], z[kc+4]) ) / grad4xbiasbot(zh[kc], zh[kc+1], zh[kc+2], zh[kc+3]);
  m7[k] = (                                                                                                       -  3.*grad4x(z[kc], z[kc+1], z[kc+2], z[kc+3]) - 27.*grad4x(z[kc+1], z[kc+2], z[kc+3], z[kc+4]) ) / grad4xbiasbot(zh[kc], zh[kc+1], zh[kc+2], zh[kc+3]);
  m8[k] = (                                                                                                                                                      +  1.*grad4x(z[kc+1], z[kc+2], z[kc+3], z[kc+4]) ) / grad4xbiasbot(zh[kc], zh[kc+1], zh[kc+2], zh[kc+3]);

  // bottom boundary + 1
  k  = 1;
  kc = kstart+k;
  m0[k] = 0.;
  m1[k] = 0.;
  m2[k] = (-23.*grad4xbiasbot(z[kc-2], z[kc-1], z[kc], z[kc+1]) -  27.*grad4x(z[kc-2], z[kc-1], z[kc], z[kc+1])                                                                                                ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
  m3[k] = ( 21.*grad4xbiasbot(z[kc-2], z[kc-1], z[kc], z[kc+1]) + 729.*grad4x(z[kc-2], z[kc-1], z[kc], z[kc+1]) +  27.*grad4x(z[kc-1], z[kc], z[kc+1], z[kc+2])                                                ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
  m4[k] = (  3.*grad4xbiasbot(z[kc-2], z[kc-1], z[kc], z[kc+1]) - 729.*grad4x(z[kc-2], z[kc-1], z[kc], z[kc+1]) - 729.*grad4x(z[kc-1], z[kc], z[kc+1], z[kc+2]) -  1.*grad4x(z[kc], z[kc+1], z[kc+2], z[kc+3]) ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
  m5[k] = ( -1.*grad4xbiasbot(z[kc-2], z[kc-1], z[kc], z[kc+1]) +  27.*grad4x(z[kc-2], z[kc-1], z[kc], z[kc+1]) + 729.*grad4x(z[kc-1], z[kc], z[kc+1], z[kc+2]) + 27.*grad4x(z[kc], z[kc+1], z[kc+2], z[kc+3]) ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
  m6[k] = (                                                                                                     -  27.*grad4x(z[kc-1], z[kc], z[kc+1], z[kc+2]) - 27.*grad4x(z[kc], z[kc+1], z[kc+2], z[kc+3]) ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
  m7[k] = (                                                                                                                                                     +  1.*grad4x(z[kc], z[kc+1], z[kc+2], z[kc+3]) ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
  m8[k] = 0.;
  
  for(int k=2; k<kmax-2; k++)
  {
    kc = kstart+k;
    m0[k] = 0.;
    m1[k] = (   1.*grad4x(z[kc-3], z[kc-2], z[kc-1], z[kc])                                                                                                                                                ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
    m2[k] = ( -27.*grad4x(z[kc-3], z[kc-2], z[kc-1], z[kc]) -  27.*grad4x(z[kc-2], z[kc-1], z[kc], z[kc+1])                                                                                                ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
    m3[k] = (  27.*grad4x(z[kc-3], z[kc-2], z[kc-1], z[kc]) + 729.*grad4x(z[kc-2], z[kc-1], z[kc], z[kc+1]) +  27.*grad4x(z[kc-1], z[kc], z[kc+1], z[kc+2])                                                ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
    m4[k] = (  -1.*grad4x(z[kc-3], z[kc-2], z[kc-1], z[kc]) - 729.*grad4x(z[kc-2], z[kc-1], z[kc], z[kc+1]) - 729.*grad4x(z[kc-1], z[kc], z[kc+1], z[kc+2]) -  1.*grad4x(z[kc], z[kc+1], z[kc+2], z[kc+3]) ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
    m5[k] = (                                               +  27.*grad4x(z[kc-2], z[kc-1], z[kc], z[kc+1]) + 729.*grad4x(z[kc-1], z[kc], z[kc+1], z[kc+2]) + 27.*grad4x(z[kc], z[kc+1], z[kc+2], z[kc+3]) ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
    m6[k] = (                                                                                               -  27.*grad4x(z[kc-1], z[kc], z[kc+1], z[kc+2]) - 27.*grad4x(z[kc], z[kc+1], z[kc+2], z[kc+3]) ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
    m7[k] = (                                                                                                                                               +  1.*grad4x(z[kc], z[kc+1], z[kc+2], z[kc+3]) ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
    m8[k] = 0.;
  }                                                                                                                                       

  // top boundary - 1
  k  = kmax-2;
  kc = kstart+k;
  m0[k] = 0.;
  m1[k] = (   1.*grad4x(z[kc-3], z[kc-2], z[kc-1], z[kc])                                                                                                                                                       ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
  m2[k] = ( -27.*grad4x(z[kc-3], z[kc-2], z[kc-1], z[kc]) -  27.*grad4x(z[kc-2], z[kc-1], z[kc], z[kc+1])                                                                                                       ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
  m3[k] = (  27.*grad4x(z[kc-3], z[kc-2], z[kc-1], z[kc]) + 729.*grad4x(z[kc-2], z[kc-1], z[kc], z[kc+1]) +  27.*grad4x(z[kc-1], z[kc], z[kc+1], z[kc+2]) -  1.*grad4xbiastop(z[kc-1], z[kc], z[kc+1], z[kc+2]) ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
  m4[k] = (  -1.*grad4x(z[kc-3], z[kc-2], z[kc-1], z[kc]) - 729.*grad4x(z[kc-2], z[kc-1], z[kc], z[kc+1]) - 729.*grad4x(z[kc-1], z[kc], z[kc+1], z[kc+2]) +  3.*grad4xbiastop(z[kc-1], z[kc], z[kc+1], z[kc+2]) ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
  m5[k] = (                                               +  27.*grad4x(z[kc-2], z[kc-1], z[kc], z[kc+1]) + 729.*grad4x(z[kc-1], z[kc], z[kc+1], z[kc+2]) + 21.*grad4xbiastop(z[kc-1], z[kc], z[kc+1], z[kc+2]) ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
  m6[k] = (                                                                                               -  27.*grad4x(z[kc-1], z[kc], z[kc+1], z[kc+2]) - 23.*grad4xbiastop(z[kc-1], z[kc], z[kc+1], z[kc+2]) ) / grad4x(zh[kc-1], zh[kc], zh[kc+1], zh[kc+2]);
  m7[k] = 0.;
  m8[k] = 0.;
  
  // top boundary
  k  = kmax-1;
  kc = kstart+k;
  m0[k] = (  1.*grad4x(z[kc-4], z[kc-3], z[kc-2], z[kc-1])                                                                                                                                                       ) / grad4xbiastop(zh[kc-2], zh[kc-1], zh[kc], zh[kc+1]);
  m1[k] = (-27.*grad4x(z[kc-4], z[kc-3], z[kc-2], z[kc-1]) -  3.*grad4x(z[kc-3], z[kc-2], z[kc-1], z[kc])                                                                                                        ) / grad4xbiastop(zh[kc-2], zh[kc-1], zh[kc], zh[kc+1]);
  m2[k] = ( 27.*grad4x(z[kc-4], z[kc-3], z[kc-2], z[kc-1]) + 81.*grad4x(z[kc-3], z[kc-2], z[kc-1], z[kc]) -  21.*grad4x(z[kc-2], z[kc-1], z[kc], z[kc+1]) +  23.*grad4xbiastop(z[kc-2], z[kc-1], z[kc], z[kc+1]) ) / grad4xbiastop(zh[kc-2], zh[kc-1], zh[kc], zh[kc+1]);
  m3[k] = ( -1.*grad4x(z[kc-4], z[kc-3], z[kc-2], z[kc-1]) - 81.*grad4x(z[kc-3], z[kc-2], z[kc-1], z[kc]) + 567.*grad4x(z[kc-2], z[kc-1], z[kc], z[kc+1]) -  69.*grad4xbiastop(z[kc-2], z[kc-1], z[kc], z[kc+1]) ) / grad4xbiastop(zh[kc-2], zh[kc-1], zh[kc], zh[kc+1]);
  m4[k] = (                                                +  3.*grad4x(z[kc-3], z[kc-2], z[kc-1], z[kc]) - 567.*grad4x(z[kc-2], z[kc-1], z[kc], z[kc+1]) - 483.*grad4xbiastop(z[kc-2], z[kc-1], z[kc], z[kc+1]) ) / grad4xbiastop(zh[kc-2], zh[kc-1], zh[kc], zh[kc+1]);
  m5[k] = (                                                                                               +  21.*grad4x(z[kc-2], z[kc-1], z[kc], z[kc+1]) + 529.*grad4xbiastop(z[kc-2], z[kc-1], z[kc], z[kc+1]) ) / grad4xbiastop(zh[kc-2], zh[kc-1], zh[kc], zh[kc+1]);
  m6[k] = 0.;
  m7[k] = 0.;
  m8[k] = 0.;
  
  fftini  = fftw_alloc_real(itot);
  fftouti = fftw_alloc_real(itot);
  fftinj  = fftw_alloc_real(jtot);
  fftoutj = fftw_alloc_real(jtot);

  allocated = true;

  return 0;
}

int cpres_g4::load()
{ 
  int itot, jtot;

  itot = grid->itot;
  jtot = grid->jtot;

  char filename[256];
  std::sprintf(filename, "%s.%07d", "fftwplan", 0);

  if(mpi->mpiid == 0)
    std::printf("Loading \"%s\"\n", filename);

  int n = fftw_import_wisdom_from_filename(filename);
  if(n == 0)
  {
    if(mpi->mpiid == 0)
      std::printf("ERROR \"%s\" does not exist\n", filename);
    return 1;
  }

  iplanf = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_R2HC, FFTW_EXHAUSTIVE);
  iplanb = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_HC2R, FFTW_EXHAUSTIVE);
  jplanf = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_R2HC, FFTW_EXHAUSTIVE);
  jplanb = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_HC2R, FFTW_EXHAUSTIVE);

  fftw_forget_wisdom();

  return 0;
}

int cpres_g4::save()
{
  int itot, jtot;

  itot = grid->itot;
  jtot = grid->jtot;

  iplanf = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_R2HC, FFTW_EXHAUSTIVE);
  iplanb = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_HC2R, FFTW_EXHAUSTIVE);
  jplanf = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_R2HC, FFTW_EXHAUSTIVE);
  jplanb = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_HC2R, FFTW_EXHAUSTIVE);

  if(mpi->mpiid == 0)
  {
    char filename[256];
    std::sprintf(filename, "%s.%07d", "fftwplan", 0);

    std::printf("Saving \"%s\"\n", filename);

    int n = fftw_export_wisdom_to_filename(filename);
    if(n == 0)
    {
      std::printf("ERROR \"%s\" cannot be saved\n", filename);
      return 1;
    }
  }

  return 0;
}

int cpres_g4::pres_in(double * restrict p, 
                      double * restrict u , double * restrict v , double * restrict w , 
                      double * restrict ut, double * restrict vt, double * restrict wt, 
                      double * restrict zh, double dt)
{
  int    ijk,ijkp,jjp,kkp;
  int    ii1,ii2,jj1,jj2,kk1,kk2,kk3;
  int    igc,jgc,kgc,kmax;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

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
  mpi->waitall();

  // write pressure as a 3d array without ghost cells
  // bottom boundary
  for(int j=0; j<grid->jmax; j++)
#pragma ivdep
    for(int i=0; i<grid->imax; i++)
    {
      ijkp = i + j*jjp;
      ijk  = i+igc + (j+jgc)*jj1 + kgc*kk1;
      p[ijkp] = grad4        ( ut[ijk-ii1] + u[ijk-ii1] / dt, ut[ijk    ] + u[ijk    ] / dt, ut[ijk+ii1] + u[ijk+ii1] / dt, ut[ijk+ii2] + u[ijk+ii2] / dt, dxi)
              + grad4        ( vt[ijk-jj1] + v[ijk-jj1] / dt, vt[ijk    ] + v[ijk    ] / dt, vt[ijk+jj1] + v[ijk+jj1] / dt, vt[ijk+jj2] + v[ijk+jj2] / dt, dyi)
              + grad4xbiasbot( wt[ijk    ] + w[ijk    ] / dt, wt[ijk+kk1] + w[ijk+kk1] / dt, wt[ijk+kk2] + w[ijk+kk2] / dt, wt[ijk+kk3] + w[ijk+kk3] / dt)
                / grad4xbiasbot( zh[kgc], zh[kgc+1], zh[kgc+2], zh[kgc+3]);
    }

  for(int k=1; k<grid->kmax-1; k++)
    for(int j=0; j<grid->jmax; j++)
#pragma ivdep
      for(int i=0; i<grid->imax; i++)
      {
        ijkp = i + j*jjp + k*kkp;
        ijk  = i+igc + (j+jgc)*jj1 + (k+kgc)*kk1;
        p[ijkp] = grad4 ( ut[ijk-ii1] + u[ijk-ii1] / dt, ut[ijk] + u[ijk] / dt, ut[ijk+ii1] + u[ijk+ii1] / dt, ut[ijk+ii2] + u[ijk+ii2] / dt, dxi)
                + grad4 ( vt[ijk-jj1] + v[ijk-jj1] / dt, vt[ijk] + v[ijk] / dt, vt[ijk+jj1] + v[ijk+jj1] / dt, vt[ijk+jj2] + v[ijk+jj2] / dt, dyi)
                + grad4x( wt[ijk-kk1] + w[ijk-kk1] / dt, wt[ijk] + w[ijk] / dt, wt[ijk+kk1] + w[ijk+kk1] / dt, wt[ijk+kk2] + w[ijk+kk2] / dt)
                  / grad4x( zh[(k-1)+kgc], zh[k+kgc], zh[(k+1)+kgc], zh[(k+2)+kgc]);
      }

  // top boundary
  for(int j=0; j<grid->jmax; j++)
#pragma ivdep
    for(int i=0; i<grid->imax; i++)
    {
      ijkp = i + j*jjp + (kmax-1)*kkp;
      ijk  = i+igc + (j+jgc)*jj1 + (kmax+kgc-1)*kk1;
      p[ijkp] = grad4        ( ut[ijk-ii1] + u[ijk-ii1] / dt, ut[ijk    ] + u[ijk    ] / dt, ut[ijk+ii1] + u[ijk+ii1] / dt, ut[ijk+ii2] + u[ijk+ii2] / dt, dxi)
              + grad4        ( vt[ijk-jj1] + v[ijk-jj1] / dt, vt[ijk    ] + v[ijk    ] / dt, vt[ijk+jj1] + v[ijk+jj1] / dt, vt[ijk+jj2] + v[ijk+jj2] / dt, dyi)
              + grad4xbiastop( wt[ijk-kk2] + w[ijk-kk2] / dt, wt[ijk-kk1] + w[ijk-kk1] / dt, wt[ijk    ] + w[ijk    ] / dt, wt[ijk+kk1] + w[ijk+kk1] / dt)
                / grad4xbiastop( zh[kmax+kgc-3], zh[kmax+kgc-2], zh[kmax+kgc-1], zh[kmax+kgc]);
    }

  return 0;
}

int cpres_g4::pres_solve(double * restrict p, double * restrict work3d, double * restrict m5calc, double * restrict dz,
                         double * restrict fftini, double * restrict fftouti, 
                         double * restrict fftinj, double * restrict fftoutj)

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

  // transpose the pressure field
  grid->transposezx(work3d,p);
  mpi->waitall();

  jj = itot;
  kk = itot*jmax;

  // do the first fourier transform
  for(int k=0; k<kblock; k++)
    for(int j=0; j<jmax; j++)
    {
#pragma ivdep
      for(int i=0; i<itot; i++)
      { 
        ijk = i + j*jj + k*kk;
        fftini[i] = work3d[ijk];
      }

      fftw_execute(iplanf);

#pragma ivdep
      for(int i=0; i<itot; i++)
      {
        ijk = i + j*jj + k*kk;
        work3d[ijk] = fftouti[i];
      }
    }

  // transpose again
  grid->transposexy(p,work3d);
  mpi->waitall();

  jj = iblock;
  kk = iblock*jtot;

  // do the second fourier transform
  for(int k=0; k<kblock; k++)
    for(int i=0; i<iblock; i++)
    {
      for(int j=0; j<jtot; j++)
      { 
        ijk = i + j*jj + k*kk;
        fftinj[j] = p[ijk];
      }

      fftw_execute(jplanf);

      for(int j=0; j<jtot; j++)
      {
        ijk = i + j*jj + k*kk;
        // shift to use p in pressure solver
        work3d[ijk] = fftoutj[j];
      }
    }

  // transpose back to original orientation
  grid->transposeyz(p,work3d);
  mpi->waitall();

  jj = iblock;
  kk = iblock*jblock;

  // solve the nonadiagonal system

  // CvH reenable later
  /*
  // create vectors that go into the solver
  //
  for(k=1; k<kmax+1; k++)
    for(j=0; j<jblock; j++)
#pragma ivdep
      for(i=0; i<iblock; i++)
      {
        // swap the mpicoords, because domain is turned 90 degrees to avoid two mpi transposes
        iindex = mpi->mpicoordy * iblock + i;
        jindex = mpi->mpicoordx * jblock + j;

        ijk  = i + j*jj + k*kk;
        m5calc[ijk] = bmati[iindex] + bmatj[jindex] + m5[k];
      }
      */

  for(j=0; j<jblock; j++)
#pragma ivdep
    for(i=0; i<iblock; i++)
    {
      // swap the mpicoords, because domain is turned 90 degrees to avoid two mpi transposes
      iindex = mpi->mpicoordy * iblock + i;
      jindex = mpi->mpicoordx * jblock + j;

      /*
      // set the BC's
      ijk = i + j*jj;
      m4[0] =       1.;
      m5[0] = -21./23.;
      m6[0] =  -3./23.;
      m7[0] =   1./23.;

      // for wave number 0, which contains average, set pressure at top to zero
      ijk  = i + j*jj + (kmax-1)*kk;
      if(iindex == 0 && jindex == 0)
      {
        m1[kmax+1] =  1./5.; 
        m2[kmax+1] = -5./5.;
        m3[kmax+1] = 15./5.;
        m4[kmax+1] =     1.;
      }
      // set dp/dz at top to zero
      else
      {
        m1[kmax+1] =   1./23.;
        m2[kmax+1] =  -3./23.;
        m3[kmax+1] = -21./23.;
        m4[kmax+1] =       1.;
      }
      */

      // set a zero gradient bc at the bottom
      m0temp[0] =      0.;
      m1temp[0] =      0.;
      m2temp[0] =      0.;
      m3temp[0] =      0.;
      m4temp[0] =      1.;
      m5temp[0] =-21./23.;
      m6temp[0] = -3./23.;
      m7temp[0] =  1./23.;
      m8temp[0] =      0.;
      ptemp [0] =      0.;

      // fill the matrix
      for(k=0; k<kmax; k++)
      {
        ijk  = i + j*jj + k*kk;
        m0temp[k+1] = m0[k];
        m1temp[k+1] = m1[k];
        m2temp[k+1] = m2[k];
        m3temp[k+1] = m3[k];
        m4temp[k+1] = m4[k] + bmati[iindex] + bmatj[jindex];
        m5temp[k+1] = m5[k];
        m6temp[k+1] = m6[k];
        m7temp[k+1] = m7[k];
        m8temp[k+1] = m8[k];
        ptemp [k+1] = p[ijk];
      }

      // set the top boundary
      m0temp[kmax+1] = 0.;
      if(iindex == 0 && jindex == 0)
      {
        m1temp[kmax+1] =  1./5.; 
        m2temp[kmax+1] = -5./5.;
        m3temp[kmax+1] = 15./5.;
        m4temp[kmax+1] =     1.;
      }
      // set dp/dz at top to zero
      else
      {
        m1temp[kmax+1] =   1./23.;
        m2temp[kmax+1] =  -3./23.;
        m3temp[kmax+1] = -21./23.;
        m4temp[kmax+1] =       1.;
      }
      m5temp[kmax+1] = 0.;
      m6temp[kmax+1] = 0.;
      m7temp[kmax+1] = 0.;
      m8temp[kmax+1] = 0.;
      ptemp [kmax+1] = 0.;

      // for now, call the solver here
      // ndma(m0temp, m1temp, m2temp, m3temp, m4temp, m5temp, m6temp, m7temp, m8temp, ptemp);

      // put back the solution
      for(k=0; k<kmax; k++)
      {
        ijk  = i + j*jj + k*kk;
        // p[ijk] = ptemp[k+1];
      }
    }

  // call ndma solver
  // ndma(m0, m1, m2, m3, m4, m5, m6, m7, m8, p);
        
  // transpose back to y
  grid->transposezy(work3d, p);
  mpi->waitall();
  
  jj = iblock;
  kk = iblock*jtot;

  // transform the second transform back
  for(int k=0; k<kblock; k++)
    for(int i=0; i<iblock; i++)
    {
      for(int j=0; j<jtot; j++)
      { 
        ijk = i + j*jj + k*kk;
        fftinj[j] = work3d[ijk];
      }

      fftw_execute(jplanb);

      for(int j=0; j<jtot; j++)
      {
        ijk = i + j*jj + k*kk;
        p[ijk] = fftoutj[j] / jtot;
      }
    }

  // transpose back to x
  grid->transposeyx(work3d, p);
  mpi->waitall();
    
  jj = itot;
  kk = itot*jmax;

  // transform the first transform back
  for(int k=0; k<kblock; k++)
    for(int j=0; j<jmax; j++)
    {
#pragma ivdep
      for(int i=0; i<itot; i++)
      { 
        ijk = i + j*jj + k*kk;
        fftini[i] = work3d[ijk];
      }

      fftw_execute(iplanb);

#pragma ivdep
      for(int i=0; i<itot; i++)
      {
        ijk = i + j*jj + k*kk;
        // swap array here to avoid unncessary 3d loop
        p[ijk] = fftouti[i] / itot;
      }
    }

  // and transpose back...
  grid->transposexz(work3d, p);
  mpi->waitall();

  jj = imax;
  kk = imax*jmax;

  int ijkp,jjp,kkp1,kkp2,kkp3;
  jjp  = grid->icells;
  kkp1 = 1*grid->icells*grid->jcells;
  kkp2 = 2*grid->icells*grid->jcells;
  kkp3 = 3*grid->icells*grid->jcells;

  // put the pressure back onto the original grid including ghost cells
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
      ijk = i + j*jjp + (grid->kstart-1)*kkp1;
      p[ijk] = (21./23.)*p[ijk+kkp1] + (3./23.)*p[ijk+kkp2] - (1./23.)*p[ijk+kkp3];
    }

  // set a zero gradient boundary at the top
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jjp + (grid->kend)*kkp1;
      p[ijk] = (21./23.)*p[ijk-kkp1] + (3./23.)*p[ijk-kkp2] - (1./23.)*p[ijk-kkp3];
    }

  // set the cyclic boundary conditions
  grid->boundary_cyclic(p);
  mpi->waitall();

  return 0;
}

int cpres_g4::pres_out(double * restrict ut, double * restrict vt, double * restrict wt, 
                       double * restrict p , double * restrict z)
{
  int    ijk,ii1,ii2,jj1,jj2,kk1,kk2;
  int    kstart;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  kstart = grid->kstart;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + kstart*kk1;
      ut[ijk] -= grad4 ( p[ijk-ii2],  p[ijk-ii1],  p[ijk],  p[ijk+ii1], dxi);
      vt[ijk] -= grad4 ( p[ijk-jj2],  p[ijk-jj1],  p[ijk],  p[ijk+jj1], dyi);
    }

  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        ut[ijk] -= grad4 ( p[ijk-ii2],  p[ijk-ii1],  p[ijk],  p[ijk+ii1], dxi);
        vt[ijk] -= grad4 ( p[ijk-jj2],  p[ijk-jj1],  p[ijk],  p[ijk+jj1], dyi);
        wt[ijk] -= grad4x( p[ijk-kk2],  p[ijk-kk1],  p[ijk],  p[ijk+kk1])
                   / grad4x( z[k-2], z[k-1], z[k], z[k+1]);
      }

  return 0;
}

int cpres_g4::ndma(double * restrict m0, double * restrict m1, double * restrict m2, double * restrict m3, double * restrict m4,
                   double * restrict m5, double * restrict m6, double * restrict m7, double * restrict m8, double * restrict p)
{
  int kmax;
  int inorm = 1;
  int k;

  kmax = grid->kmax;

  // LU factorization
  k = 0;
  m0[k] = 1.;
  m1[k] = 1.;
  m2[k] = 1.;
  if(inorm == 1)
  {
     m3[k] = 1.          / m4[k]; // padding, and used in nonadss to normalize 1st eqn.
     m4[k] = 1.;
     m5[k] = m5[k]*m3[k];
     m6[k] = m6[k]*m3[k];
     m7[k] = m7[k]*m3[k];
     m8[k] = m8[k]*m3[k];
  }
  else
     m3[k] = 1.; // padding

  k = 1;
  m0[k] = 1.; // padding
  m1[k] = 1.; // padding
  m2[k] = 1.; // padding
  m3[k] = m3[k]                 / m4[k-1];
  m4[k] = m4[k] - m3[k]*m5[k-1];
  m5[k] = m5[k] - m3[k]*m6[k-1];
  m6[k] = m6[k] - m3[k]*m7[k-1];
  m7[k] = m7[k] - m3[k]*m8[k-1];

  k = 2;
  m0[k] = 1.; // padding
  m1[k] = 1.; // padding
  m2[k] =   m2[k]                                   / m4[k-2];
  m3[k] = ( m3[k]                 - m2[k]*m5[k-2] ) / m4[k-1];
  m4[k] =   m4[k] - m3[k]*m5[k-1] - m2[k]*m6[k-2];
  m5[k] =   m5[k] - m3[k]*m6[k-1] - m2[k]*m7[k-2];
  m6[k] =   m6[k] - m3[k]*m7[k-1] - m2[k]*m8[k-2];
  m7[k] =   m7[k] - m3[k]*m8[k-1];

  k = 3;
  m0[k] = 1.; // padding
  m1[k] =   m1[k]                                                   / m4[k-3];
  m2[k] = ( m2[k]                                 - m1[k]*m5[k-3] ) / m4[k-2];
  m3[k] = ( m3[k]                 - m2[k]*m5[k-2] - m1[k]*m6[k-3] ) / m4[k-1];
  m4[k] =   m4[k] - m3[k]*m5[k-1] - m2[k]*m6[k-2] - m1[k]*m7[k-3]; 
  m5[k] =   m5[k] - m3[k]*m6[k-1] - m2[k]*m7[k-2] - m1[k]*m8[k-3];
  m6[k] =   m6[k] - m3[k]*m7[k-1] - m2[k]*m8[k-2];
  m7[k] =   m7[k] - m3[k]*m8[k-1];

  for(k=4; k<kmax-1; k++)
  {
    m0[k] =   m0[k]                                                                  / m4[k-4];
    m1[k] = ( m1[k]                                                 - m0[k]*m5[k-4]) / m4[k-3];
    m2[k] = ( m2[k]                                 - m1[k]*m5[k-3] - m0[k]*m6[k-4]) / m4[k-2];
    m3[k] = ( m3[k]                 - m2[k]*m5[k-2] - m1[k]*m6[k-3] - m0[k]*m7[k-4]) / m4[k-1];
    m4[k] =   m4[k] - m3[k]*m5[k-1] - m2[k]*m6[k-2] - m1[k]*m7[k-3] - m0[k]*m8[k-4];
    m5[k] =   m5[k] - m3[k]*m6[k-1] - m2[k]*m7[k-2] - m1[k]*m8[k-3];
    m6[k] =   m6[k] - m3[k]*m7[k-1] - m2[k]*m8[k-2];
    m7[k] =   m7[k] - m3[k]*m8[k-1];
  }
  m8[k-1] = 1.; // padding

  k = kmax-1;
  m0[k] =   m0[k]                                                                  / m4[k-4];
  m1[k] = ( m1[k]                                                 - m0[k]*m5[k-4]) / m4[k-3];
  m2[k] = ( m2[k]                                 - m1[k]*m5[k-3] - m0[k]*m6[k-4]) / m4[k-2];
  m3[k] = ( m3[k]                 - m2[k]*m5[k-2] - m1[k]*m6[k-3] - m0[k]*m7[k-4]) / m4[k-1];
  m4[k] =   m4[k] - m3[k]*m5[k-1] - m2[k]*m6[k-2] - m1[k]*m7[k-3] - m0[k]*m8[k-4];
  m5[k] =   m5[k] - m3[k]*m6[k-1] - m2[k]*m7[k-2] - m1[k]*m8[k-3];
  m6[k] =   m6[k] - m3[k]*m7[k-1] - m2[k]*m8[k-2];
  m7[k] =  1.; // padding
  m8[k] =  1.; // padding

  k = kmax;
  m0[k] =   m0[k]                                                                  / m4[k-4];
  m1[k] = ( m1[k]                                                 - m0[k]*m5[k-4]) / m4[k-3];
  m2[k] = ( m2[k]                                 - m1[k]*m5[k-3] - m0[k]*m6[k-4]) / m4[k-2];
  m3[k] = ( m3[k]                 - m2[k]*m5[k-2] - m1[k]*m6[k-3] - m0[k]*m7[k-4]) / m4[k-1];
  m4[k] =   m4[k] - m3[k]*m5[k-1] - m2[k]*m6[k-2] - m1[k]*m7[k-3] - m0[k]*m8[k-4];
  m5[k] =   m5[k] - m3[k]*m6[k-1] - m2[k]*m7[k-2] - m1[k]*m8[k-3];
  m6[k] = 1.; // padding
  m7[k] = 1.; // padding
  m8[k] = 1.; // padding

  k = kmax+1;
  m0[k] =   m0[k]                                                                  / m4[k-4];
  m1[k] = ( m1[k]                                                 - m0[k]*m5[k-4]) / m4[k-3];
  m2[k] = ( m2[k]                                 - m1[k]*m5[k-3] - m0[k]*m6[k-4]) / m4[k-2];
  m3[k] = ( m3[k]                 - m2[k]*m5[k-2] - m1[k]*m6[k-3] - m0[k]*m7[k-4]) / m4[k-1];
  m4[k] =   m4[k] - m3[k]*m5[k-1] - m2[k]*m6[k-2] - m1[k]*m7[k-3] - m0[k]*m8[k-4];
  m5[k] = 1.; // padding
  m6[k] = 1.; // padding
  m7[k] = 1.; // padding
  m8[k] = 1.; // padding

  // Backward substitution 
  int i,j,jj,ijk,ij;
  int kk1,kk2,kk3,kk4;
  int iblock,jblock;

  iblock = 1; // grid->iblock; CvH vectorize later
  jblock = 1; // grid->jblock;

  jj  = iblock;
  kk1 = 1*iblock*jblock;
  kk2 = 2*iblock*jblock;
  kk3 = 3*iblock*jblock;
  kk4 = 4*iblock*jblock;

  // Solve Ly=p, forward
  for(j=0;j<jblock;j++)
#pragma ivdep
    for(i=0;i<iblock;i++)
    {
      ij = i + j*jj;
      p[ij    ] =             p[ij    ]*m3[0]; // Normalize first eqn. See NONADFS
      p[ij+kk1] = p[ij+kk1] - p[ij    ]*m3[1];
      p[ij+kk2] = p[ij+kk2] - p[ij+kk1]*m3[2] - p[ij    ]*m2[2];
      p[ij+kk3] = p[ij+kk3] - p[ij+kk2]*m3[3] - p[ij+kk1]*m2[3] -p[ij]*m1[3];
    }

  for(k=4; k<kmax+2; k++)
    for(j=0;j<jblock;j++)
#pragma ivdep
      for(i=0;i<iblock;i++)
      {
        ijk = i + j*jj + k*kk1;
        p[ijk] = p[ijk] - p[ijk-kk1]*m3[k] - p[ijk-kk2]*m2[k] - p[ijk-kk3]*m1[k] - p[ijk-kk4]*m0[k];
      }

  // Solve Ux=y, backward
  for(j=0;j<jblock;j++)
#pragma ivdep
    for(i=0;i<iblock;i++)
    {
      ijk = i + j*jj + k*(kmax+1);
      p[ijk    ] =   p[ijk    ]                                                              / m4[k  ];
      p[ijk-kk1] = ( p[ijk-kk1] - p[ijk    ]*m5[k-1] )                                       / m4[k-1];
      p[ijk-kk2] = ( p[ijk-kk2] - p[ijk-kk1]*m5[k-2] - p[ijk    ]*m6[k-2] )                  / m4[k-2];
      p[ijk-kk3] = ( p[ijk-kk3] - p[ijk-kk2]*m5[k-3] - p[ijk-kk1]*m6[k-3] - p[ijk]*m7[k-3] ) / m4[k-3];
    }

  for(k=kmax-3; k>=0; k--)
    for(j=0;j<jblock;j++)
#pragma ivdep
      for(i=0;i<iblock;i++)
      {
        ijk = i + j*jj + k*kk1;
        p[ijk] = ( p[ijk] - p[ijk+kk1]*m5[k] - p[ijk+kk2]*m6[k] - p[ijk+kk3]*m7[k] - p[ijk+kk4]*m8[k] ) / m4[k];
      }

  return 0;
}

/*
// tridiagonal matrix solver, taken from Numerical Recipes, Press
int cpres_g4::tdma(double * restrict a, double * restrict b, double * restrict c, 
                double * restrict p, double * restrict work2d, double * restrict work3d)
                
{
  int i,j,k,jj,kk,ijk,ij;
  int iblock,jblock,kmax;

  iblock = grid->iblock;
  jblock = grid->jblock;
  kmax = grid->kmax;

  jj = iblock;
  kk = iblock*jblock;

  for(j=0;j<jblock;j++)
#pragma ivdep
    for(i=0;i<iblock;i++)
    {
      ij = i + j*jj;
      work2d[ij] = b[ij];
    }

  for(j=0;j<jblock;j++)
#pragma ivdep
    for(i=0;i<iblock;i++)
    {
      ij = i + j*jj;
      p[ij] /= work2d[ij];
    }

  for(k=1; k<kmax; k++)
  {
    for(j=0;j<jblock;j++)
#pragma ivdep
      for(i=0;i<iblock;i++)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + k*kk;
        work3d[ijk] = c[k-1] / work2d[ij];
      }
    for(j=0;j<jblock;j++)
#pragma ivdep
      for(i=0;i<iblock;i++)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + k*kk;
        work2d[ij] = b[ijk] - a[k]*work3d[ijk];
      }
    for(j=0;j<jblock;j++)
#pragma ivdep
      for(i=0;i<iblock;i++)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + k*kk;
        p[ijk] -= a[k]*p[ijk-kk];
        p[ijk] /= work2d[ij];
      }
  }

  for(k=kmax-2; k>=0; k--)
    for(j=0;j<jblock;j++)
#pragma ivdep
      for(i=0;i<iblock;i++)
      {
        ijk = i + j*jj + k*kk;
        p[ijk] -= work3d[ijk+kk]*p[ijk+kk];
      }

  return 0;
}
*/

double cpres_g4::calcdivergence(double * restrict u, double * restrict v, double * restrict w, double * restrict zh)
{
  int    ijk,ii1,ii2,jj1,jj2,kk1,kk2,kk3;
  int    kstart,kend;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

  kstart = grid->kstart;
  kend   = grid->kend;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  double div, divmax;
  divmax = 0;

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + kstart*kk1;
      div = grad4        ( u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2], dxi)
          + grad4        ( v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2], dyi)
          + grad4xbiasbot( w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3])
            / grad4xbiasbot( zh[kstart], zh[kstart+1], zh[kstart+2], zh[kstart+3]);

      divmax = std::max(divmax, std::abs(div));
    }

  for(int k=grid->kstart+1; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        div = grad4 ( u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2], dxi)
            + grad4 ( v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2], dyi)
            + grad4x( w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2])
              / grad4x( zh[k-1], zh[k], zh[k+1], zh[k+2]);

        divmax = std::max(divmax, std::abs(div));
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      div = grad4        ( u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2], dxi)
          + grad4        ( v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2], dyi)
          + grad4xbiastop( w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1])
            / grad4xbiastop( zh[kend-3], zh[kend-2], zh[kend-1], zh[kend]);

      divmax = std::max(divmax, std::abs(div));
    }


  grid->getmax(&divmax);

  return divmax;
}

inline double cpres_g4::grad4(const double a, const double b, const double c, const double d, const double dxi)
{
  return ( -(1./24.)*(d-a) + (27./24.)*(c-b) ) * dxi;
}

inline double cpres_g4::grad4x(const double a, const double b, const double c, const double d)
{
  return (-(d-a) + 27.*(c-b)); 
}

inline double cpres_g4::grad4xbiasbot(const double a, const double b, const double c, const double d)
{
  return (-23.*a + 21.*b + 3.*c - d);
}

inline double cpres_g4::grad4xbiastop(const double a, const double b, const double c, const double d)
{
  return ( 23.*d - 21.*c - 3.*b + a);
}

