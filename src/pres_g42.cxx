#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fftw3.h>
#include "grid.h"
#include "fields.h"
#include "pres_g42.h"
#include "defines.h"

cpres_g42::cpres_g42(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  // std::printf("Creating instance of object pres_g42\n");
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;

  allocated = false;
}

cpres_g42::~cpres_g42()
{
  if(allocated)
  {
    delete[] a;
    delete[] c;
    delete[] work2d;

    delete[] bmati;
    delete[] bmatj;
  }

  // std::printf("Destroying instance of object pres_g42\n");
}

int cpres_g42::init()
{
  int imax, jmax, kmax;
  int itot, jtot;

  itot = grid->itot;
  jtot = grid->jtot;
  imax = grid->imax;
  jmax = grid->jmax;
  kmax = grid->kmax;

  bmati = new double[itot];
  bmatj = new double[jtot];
  
  // allocate help variables for the matrix solver
  a = new double[kmax];
  c = new double[kmax];
  work2d = new double[imax*jmax];

  allocated = true;

  return 0;
}


int cpres_g42::setvalues()
{
  int imax, jmax, kmax;
  int itot, jtot, kgc;

  itot = grid->itot;
  jtot = grid->jtot;
  imax = grid->imax;
  jmax = grid->jmax;
  kmax = grid->kmax;
  kgc  = grid->kgc;

  // compute the modified wave numbers of the 42 order scheme
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

  // create vectors that go into the tridiagonal matrix solver
  for(int k=0; k<kmax; k++)
  {
    a[k] = grid->dz[k+kgc] * grid->dzhi[k+kgc  ];
    c[k] = grid->dz[k+kgc] * grid->dzhi[k+kgc+1];
  }

  return 0;
}

int cpres_g42::pres_in(double * restrict p, 
                      double * restrict u , double * restrict v , double * restrict w , 
                      double * restrict ut, double * restrict vt, double * restrict wt, 
                      double * restrict dzi,
                      double dt)
{
  int    ijk,kk,ijkp,jjp,kkp;
  int    ii1,ii2,jj1,jj2;
  int    igc,jgc,kgc;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk  = grid->icells*grid->jcells;

  jjp = grid->imax;
  kkp = grid->imax*grid->jmax;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  igc = grid->igc;
  jgc = grid->jgc;
  kgc = grid->kgc;

  // set the cyclic boundary conditions for the tendencies
  grid->boundary_cyclic(ut);
  grid->boundary_cyclic(vt);
  grid->boundary_cyclic(wt);

  // write pressure as a 3d array without ghost cells
  for(int k=0; k<grid->kmax; k++)
    for(int j=0; j<grid->jmax; j++)
#pragma ivdep
      for(int i=0; i<grid->imax; i++)
      {
        ijkp = i + j*jjp + k*kkp;
        ijk  = i+igc + (j+jgc)*jj1 + (k+kgc)*kk;
        p[ijkp] = grad4( ut[ijk-ii1] + u[ijk-ii1] / dt, ut[ijk] + u[ijk] / dt, ut[ijk+ii1] + u[ijk+ii1] / dt, ut[ijk+ii2] + u[ijk+ii2] / dt, dxi)
                + grad4( vt[ijk-jj1] + v[ijk-jj1] / dt, vt[ijk] + v[ijk] / dt, vt[ijk+jj1] + v[ijk+jj1] / dt, vt[ijk+jj2] + v[ijk+jj2] / dt, dyi)
                + ( (wt[ijk+kk] + w[ijk+kk] / dt) - (wt[ijk] + w[ijk] / dt) ) * dzi[k+kgc];
      }

  return 0;
}

int cpres_g42::pres_solve(double * restrict p, double * restrict work3d, double * restrict b, double * restrict dz,
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

  grid->fftforward(p, work3d, fftini, fftouti, fftinj, fftoutj);

  jj = iblock;
  kk = iblock*jblock;

  // solve the tridiagonal system
  // create vectors that go into the tridiagonal matrix solver
  for(k=0; k<kmax; k++)
    for(j=0; j<jblock; j++)
#pragma ivdep
      for(i=0; i<iblock; i++)
      {
        // swap the mpicoords, because domain is turned 90 degrees to avoid two mpi transposes
        iindex = mpi->mpicoordy * iblock + i;
        jindex = mpi->mpicoordx * jblock + j;

        ijk  = i + j*jj + k*kk;
        b[ijk] = dz[k+kgc]*dz[k+kgc] * (bmati[iindex]+bmatj[jindex]) - (a[k]+c[k]);
        p[ijk] = dz[k+kgc]*dz[k+kgc] * p[ijk];
      }

  for(j=0; j<jblock; j++)
#pragma ivdep
    for(i=0; i<iblock; i++)
    {
      iindex = mpi->mpicoordy * iblock + i;
      jindex = mpi->mpicoordx * jblock + j;

      // substitute BC's
      ijk = i + j*jj;
      b[ijk] += a[0];

      // for wave number 0, which contains average, set pressure at top to zero
      ijk  = i + j*jj + (kmax-1)*kk;
      if(iindex == 0 && jindex == 0)
        b[ijk] -= c[kmax-1];
      // set dp/dz at top to zero
      else
        b[ijk] += c[kmax-1];
    }

  // call tdma solver
  tdma(a, b, c, p, work2d, work3d);

  grid->fftbackward(p, work3d, fftini, fftouti, fftinj, fftoutj);
        
  jj = imax;
  kk = imax*jmax;

  int ijkp,jjp,kkp;
  jjp = grid->icells;
  kkp = grid->icells*grid->jcells;

  // put the pressure back onto the original grid including ghost cells
  for(int k=0; k<grid->kmax; k++)
    for(int j=0; j<grid->jmax; j++)
#pragma ivdep
      for(int i=0; i<grid->imax; i++)
      {
        ijkp = i+igc + (j+jgc)*jjp + (k+kgc)*kkp;
        ijk  = i + j*jj + k*kk;
        p[ijkp] = work3d[ijk];
      }

  // set the boundary conditions
  // set a zero gradient boundary at the bottom
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jjp + grid->kstart*kkp;
      p[ijk-kkp] = p[ijk];
    }

  // set the cyclic boundary conditions
  grid->boundary_cyclic(p);

  return 0;
}

int cpres_g42::pres_out(double * restrict ut, double * restrict vt, double * restrict wt, 
                       double * restrict p , double * restrict dzhi)
{
  int    ijk,ii1,ii2,jj1,jj2,kk;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk  = grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk;
        ut[ijk] -= grad4( p[ijk-ii2],  p[ijk-ii1],  p[ijk],  p[ijk+ii1], dxi);
        vt[ijk] -= grad4( p[ijk-jj2],  p[ijk-jj1],  p[ijk],  p[ijk+jj1], dyi);
        wt[ijk] -= (p[ijk] - p[ijk-kk]) * dzhi[k];
      }

  return 0;
}

// tridiagonal matrix solver, taken from Numerical Recipes, Press
int cpres_g42::tdma(double * restrict a, double * restrict b, double * restrict c, 
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

double cpres_g42::calcdivergence(double * restrict u, double * restrict v, double * restrict w, double * restrict dzi)
{
  int    ijk,ii1,ii2,jj1,jj2,kk;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk  = grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  double div, divmax;
  div    = 0;
  divmax = 0;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk;
        div = grad4( u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2], dxi)
            + grad4( v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2], dyi)
            + (w[ijk+kk]-w[ijk])*dzi[k];

        divmax = std::max(divmax, std::abs(div));
      }

  grid->getmax(&divmax);

  return divmax;
}

inline double cpres_g42::grad4(const double a, const double b, const double c, const double d, const double dxi)
{
  return ( -(1./24.)*(d-a) + (27./24.)*(c-b) ) * dxi;
}

