#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fftw3.h>
#include "grid.h"
#include "fields.h"
#include "pres.h"

cpres::cpres(cgrid *gridin, cfields *fieldsin)
{
  std::printf("Creating instance of object pres\n");
  grid   = gridin;
  fields = fieldsin;
}

cpres::~cpres()
{
  std::printf("Destroying instance of object pres\n");
}

int cpres::exec(double dt)
{
  // cyclic boundaries for tendencies 
  (*fields->ut).boundary_cyclic();
  (*fields->vt).boundary_cyclic();
  (*fields->wt).boundary_cyclic();

  // create the input for the pressure solver
  pres_2nd_in((*fields->p ).data,
              (*fields->u ).data, (*fields->v ).data, (*fields->w ).data,
              (*fields->ut).data, (*fields->vt).data, (*fields->wt).data, 
              grid->dzi, dt);

  // solve the system
  pres_2nd_solve((*fields->p).data, grid->dz);

  // set the boundary conditions
  (*fields->p).boundary_cyclic();
  (*fields->p).boundary_bottop(1);
  
  // get the pressure tendencies from the pressure field
  pres_2nd_out((*fields->ut).data, (*fields->vt).data, (*fields->wt).data, 
               (*fields->p ).data, grid->dzhi);

  return 0;
}

int cpres::init()
{
  pres_2nd_init();

  return 0;
}

int cpres::divergence()
{
  double divmax;
  divmax = calcdivergence((*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi);

  std::printf("divmax = %24.14f\n", divmax);

  return 0;
}

int cpres::pres_2nd_init()
{
  int itot, jtot, ktot, kgc;

  itot = grid->itot;
  jtot = grid->jtot;
  ktot = grid->ktot;
  kgc  = grid->kgc;

  fftini  = new double[itot];
  fftouti = new double[itot];
  fftinj  = new double[jtot];
  fftoutj = new double[jtot];

  iplanf = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_R2HC, FFTW_PATIENT);
  iplanb = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_HC2R, FFTW_PATIENT);
  jplanf = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_R2HC, FFTW_PATIENT);
  jplanb = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_HC2R, FFTW_PATIENT);

  bmati = new double[itot];
  bmatj = new double[jtot];
  
  // compute the modified wave numbers of the 2nd order scheme
  double dxidxi = 1./(grid->dx*grid->dx);
  double dyidyi = 1./(grid->dy*grid->dy);

  const double pi = std::acos(-1.);
  std::printf("%24.20f\n", pi);

  for(int j=0; j<jtot/2+1; j++)
    bmatj[j] = 2. * (std::cos(2.*pi*(double)j/(double)jtot)-1.) * dyidyi;

  for(int j=jtot/2+1; j<jtot; j++)
    bmatj[j] = bmatj[jtot-j];

  for(int i=0; i<itot/2+1; i++)
    bmati[i] = 2. * (std::cos(2.*pi*(double)i/(double)itot)-1.) * dxidxi;

  for(int i=itot/2+1; i<itot; i++)
    bmati[i] = bmati[itot-i];

  // allocate help variables for the matrix solver
  a  = new double[ktot];
  b  = new double[ktot];
  c  = new double[ktot];
  d  = new double[ktot];

  xin  = new double[ktot];
  xout = new double[ktot];

  // create vectors that go into the tridiagonal matrix solver
  for(int k=0; k<ktot; k++)
  {
    a[k] = grid->dz[k+kgc] * grid->dzhi[k+kgc  ];
    c[k] = grid->dz[k+kgc] * grid->dzhi[k+kgc+1];
  }

  return 0;
}

int cpres::pres_2nd_in(double * __restrict__ p, 
                       double * __restrict__ u , double * __restrict__ v , double * __restrict__ w , 
                       double * __restrict__ ut, double * __restrict__ vt, double * __restrict__ wt, 
                       double * __restrict__ dzi,
                       double dt)
{
  int    ijk,icells,ijcells,ii,jj,kk;
  double dxi, dyi;

  icells  = grid->icells;
  ijcells = grid->icells*grid->jcells;

  ii = 1;
  jj = 1*icells;
  kk = 1*ijcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*icells + k*ijcells;
        p[ijk] = ( (ut[ijk+ii] + u[ijk+ii] / dt) - (ut[ijk] + u[ijk] / dt) ) * dxi
               + ( (vt[ijk+jj] + v[ijk+jj] / dt) - (vt[ijk] + v[ijk] / dt) ) * dyi
               + ( (wt[ijk+kk] + w[ijk+kk] / dt) - (wt[ijk] + w[ijk] / dt) ) * dzi[k];
           
      }

  return 0;
}

int cpres::pres_2nd_solve(double * __restrict__ p, double * __restrict__ dz)
{
  int i,j,k,ii,jj,kk,ijk;
  int imax, jmax, kmax;
  int itot, jtot, ktot;
  int igc, jgc, kgc;
  int iindex, jindex;

  imax = grid->imax;
  jmax = grid->jmax;
  kmax = grid->kmax;
  itot = grid->itot;
  jtot = grid->jtot;
  ktot = grid->ktot;
  igc  = grid->igc;
  jgc  = grid->jgc;
  kgc  = grid->kgc;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  // do the first fourier transform
  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
    {
      for(int i=0;i<itot;i++)
      { 
        ijk = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        fftini[i] = p[ijk];
      }

      fftw_execute_r2r(iplanf, fftini, fftouti);

      for(int i=0;i<itot;i++)
      {
        ijk = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        p[ijk] = fftouti[i];
      }
    }

  // TRANSPOSE

  // do the second fourier transform
  for(int k=0; k<kmax; k++)
    for(int i=0; i<imax; i++)
    {
      for(int j=0;j<jtot;j++)
      { 
        ijk = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        fftinj[j] = p[ijk];
      }

      fftw_execute_r2r(jplanf, fftinj, fftoutj);

      for(int j=0;j<jtot;j++)
      {
        ijk = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        p[ijk] = fftoutj[j];
      }
    }

  // TRANSPOSE

  // solve the tridiagonal system

  for(j=0; j<jmax; j++)
    for(i=0; i<imax; i++)
    {
      // iindex = mpicoordx * imax + i
      // jindex = mpicoordy * jmax + j
      iindex = i;
      jindex = j;

      // create vectors that go into the tridiagonal matrix solver
      for(k=0; k<ktot; k++)
      {
        ijk = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        b  [k] = dz[k+kgc]*dz[k+kgc] * (bmati[iindex]+bmatj[jindex]) - (a[k]+c[k]);
        xin[k] = dz[k+kgc]*dz[k+kgc] * p[ijk];
      }

      // substitute BC's
      b[0] += a[0];

      // for wave number 0, which contains average, set pressure at top to zero
      if(iindex == 0 && jindex == 0)
        b[ktot] = b[ktot] - c[ktot];
      // set dp/dz at top to zero
      else
        b[ktot] = b[ktot] + c[ktot];

      // call tdma solver
      tdma(a, b, c, xin, xout, d, ktot);
        
      // update the pressure (in fourier space, still)
      for(int k=0;k<ktot;k++)
      {
        ijk = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        p[ijk] = xout[k];
      }
    }

  // TRANSPOSE
  
  // transform the second transform back
  for(int k=0; k<kmax; k++)
    for(int i=0; i<imax; i++)
    {
      for(int j=0;j<jtot;j++)
      { 
        ijk = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        fftinj[j] = p[ijk];
      }

      fftw_execute_r2r(jplanb, fftinj, fftoutj);

      for(int j=0;j<jtot;j++)
      {
        ijk = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        p[ijk] = fftoutj[j];
      }
    }

    // TRANSPOSE
    
    // transform the first transform back
    for(int k=0; k<kmax; k++)
      for(int j=0; j<jmax; j++)
      {
        for(int i=0;i<itot;i++)
        { 
          ijk = i+igc + (j+jgc)*jj + (k+kgc)*kk;
          fftini[i] = p[ijk];
        }

        fftw_execute_r2r(iplanb, fftini, fftouti);

        for(int i=0;i<itot;i++)
        {
          ijk = i+igc + (j+jgc)*jj + (k+kgc)*kk;
          p[ijk] = fftouti[i];
        }
      }

  return 0;
}

int cpres::pres_2nd_out(double * __restrict__ ut, double * __restrict__ vt, double * __restrict__ wt, 
                        double * __restrict__ p , double * __restrict__ dzhi)
{
  int    ijk,icells,ijcells,ii,jj,kk;
  double dxi, dyi;

  icells  = grid->icells;
  ijcells = grid->icells*grid->jcells;

  ii = 1;
  jj = 1*icells;
  kk = 1*ijcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*icells + k*ijcells;
        ut[ijk] = ut[ijk] - (p[ijk] - p[ijk-ii]) * dxi;
        vt[ijk] = vt[ijk] - (p[ijk] - p[ijk-jj]) * dyi;
        wt[ijk] = wt[ijk] - (p[ijk] - p[ijk-kk]) * dzhi[k];
      }

  return 0;
}

// tridiagonal matrix solver, taken from Numerical Recipes, Press
int cpres::tdma(double * __restrict__ a,   double * __restrict__ b,    double * __restrict__ c, 
                double * __restrict__ xin, double * __restrict__ xout, double * __restrict__ gam, 
                int size)
{
  int k;
  double tmp;

  tmp = b[0];

  xout[0] = xin[0] / tmp;

  for(k=1; k<size-1; k++)
  {
    gam[k]  = c[k-1] / tmp;
    tmp     = b[k] - a[k]*gam[k];
    xout[k] = (xin[k] - a[k] * xout[k-1]) / tmp;
  }

  gam[k]  = c[k-1] / tmp;
  tmp     = b[k] - a[k]*gam[k];
  xout[k] = (xin[k] - a[k]*xout[k-1]) / tmp;

  for(k=size-2; k>=0; k--)
    xout[k] = xout[k] - gam[k+1]*xout[k+1];

  return 0;
}

double cpres::calcdivergence(double * __restrict__ u, double * __restrict__ v, double * __restrict__ w, double * __restrict__ dzi)
{
  int    ijk,icells,ijcells,ii,jj,kk;
  double dxi, dyi;

  icells  = grid->icells;
  ijcells = grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  ii = 1;
  jj = 1*icells;
  kk = 1*ijcells;

  double div, divmax;
  div    = 0;
  divmax = 0;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*icells + k*ijcells;
        div = (u[ijk+ii]-u[ijk])*dxi + (v[ijk+jj]-v[ijk])*dyi + (w[ijk+kk]-w[ijk])*dzi[k];

        // if(k==1 && j==1)
        //  std::printf("%d, %24.14E\n", i, div);

        divmax = std::max(divmax, std::abs(div));
      }

  return divmax;
}

