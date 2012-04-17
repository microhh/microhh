#include <cstdio>
#include <cmath>
#include <algorithm>
#include <fftw3.h>
#include "grid.h"
#include "fields.h"
#include "pres.h"

#define restrict __restrict__

cpres::cpres(cgrid *gridin, cfields *fieldsin)
{
  std::printf("Creating instance of object pres\n");
  grid   = gridin;
  fields = fieldsin;

  allocated = false;
}

cpres::~cpres()
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

    delete[] a;
    delete[] b;
    delete[] c;
    delete[] work2d;
    delete[] work3d;

    delete[] bmati;
    delete[] bmatj;
  }

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
  pres_2nd_solve((*fields->p).data, grid->dz, fftini, fftouti, fftinj, fftoutj);

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

  allocated = true;

  return 0;
}

int cpres::divergence()
{
  double divmax;
  divmax = calcdivergence((*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi);

  std::printf("divmax = %24.14E\n", divmax);

  return 0;
}

int cpres::load()
{ 
  int itot, jtot;

  itot = grid->itot;
  jtot = grid->jtot;

  fftini  = fftw_alloc_real(itot);
  fftouti = fftw_alloc_real(itot);
  fftinj  = fftw_alloc_real(jtot);
  fftoutj = fftw_alloc_real(jtot);

  char filename[256];
  std::sprintf(filename, "%s.%06d", "fftwplan", 0);
  int n = fftw_import_wisdom_from_filename(filename);
  if(n == 0)
  {
    std::printf("ERROR \"%s\" does not exist\n", filename);
    return 1;
  }
  else
    std::printf("Loading \"%s\"\n", filename);

  iplanf = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_R2HC, FFTW_EXHAUSTIVE);
  iplanb = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_HC2R, FFTW_EXHAUSTIVE);
  jplanf = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_R2HC, FFTW_EXHAUSTIVE);
  jplanb = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_HC2R, FFTW_EXHAUSTIVE);

  fftw_forget_wisdom();

  return 0;
}

int cpres::save()
{
  int itot, jtot;

  itot = grid->itot;
  jtot = grid->jtot;

  fftini  = fftw_alloc_real(itot);
  fftouti = fftw_alloc_real(itot);
  fftinj  = fftw_alloc_real(jtot);
  fftoutj = fftw_alloc_real(jtot);

  iplanf = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_R2HC, FFTW_EXHAUSTIVE);
  iplanb = fftw_plan_r2r_1d(itot, fftini, fftouti, FFTW_HC2R, FFTW_EXHAUSTIVE);
  jplanf = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_R2HC, FFTW_EXHAUSTIVE);
  jplanb = fftw_plan_r2r_1d(jtot, fftinj, fftoutj, FFTW_HC2R, FFTW_EXHAUSTIVE);

  char filename[256];
  std::sprintf(filename, "%s.%06d", "fftwplan", 0);
  int n = fftw_export_wisdom_to_filename(filename);
  if(n == 0)
  {
    std::printf("ERROR \"%s\" cannot be saved\n", filename);
    return 1;
  }
  else
    std::printf("Saving \"%s\"\n", filename);

  return 0;
}

int cpres::pres_2nd_init()
{
  int itot, jtot, ktot, kgc;

  itot = grid->itot;
  jtot = grid->jtot;
  ktot = grid->ktot;
  kgc  = grid->kgc;

  bmati = new double[itot];
  bmatj = new double[jtot];
  
  // compute the modified wave numbers of the 2nd order scheme
  double dxidxi = 1./(grid->dx*grid->dx);
  double dyidyi = 1./(grid->dy*grid->dy);

  const double pi = std::acos(-1.);

  for(int j=0; j<jtot/2+1; j++)
    bmatj[j] = 2. * (std::cos(2.*pi*(double)j/(double)jtot)-1.) * dyidyi;

  for(int j=jtot/2+1; j<jtot; j++)
    bmatj[j] = bmatj[jtot-j];

  for(int i=0; i<itot/2+1; i++)
    bmati[i] = 2. * (std::cos(2.*pi*(double)i/(double)itot)-1.) * dxidxi;

  for(int i=itot/2+1; i<itot; i++)
    bmati[i] = bmati[itot-i];

  // allocate help variables for the matrix solver
  a = new double[ktot];
  b = new double[itot*jtot*ktot];
  c = new double[ktot];
  work2d = new double[itot*jtot];
  work3d = new double[itot*jtot*ktot];

  // create vectors that go into the tridiagonal matrix solver
  for(int k=0; k<ktot; k++)
  {
    a[k] = grid->dz[k+kgc] * grid->dzhi[k+kgc  ];
    c[k] = grid->dz[k+kgc] * grid->dzhi[k+kgc+1];
  }

  return 0;
}

int cpres::pres_2nd_in(double * restrict p, 
                       double * restrict u , double * restrict v , double * restrict w , 
                       double * restrict ut, double * restrict vt, double * restrict wt, 
                       double * restrict dzi,
                       double dt)
{
  int    ijk,ii,jj,kk;
  double dxi,dyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        p[ijk] = ( (ut[ijk+ii] + u[ijk+ii] / dt) - (ut[ijk] + u[ijk] / dt) ) * dxi
               + ( (vt[ijk+jj] + v[ijk+jj] / dt) - (vt[ijk] + v[ijk] / dt) ) * dyi
               + ( (wt[ijk+kk] + w[ijk+kk] / dt) - (wt[ijk] + w[ijk] / dt) ) * dzi[k];
      }

  return 0;
}

int cpres::pres_2nd_solve(double * restrict p, double * restrict dz,
                          double * restrict fftini, double * restrict fftouti, 
                          double * restrict fftinj, double * restrict fftoutj)

{
  int i,j,k,jj,kk,ijk,ijkb;
  int imax,jmax,kmax;
  int itot,jtot,ktot;
  int igc,jgc,kgc;
  int iindex,jindex;

  imax = grid->imax;
  jmax = grid->jmax;
  kmax = grid->kmax;
  itot = grid->itot;
  jtot = grid->jtot;
  ktot = grid->ktot;
  igc  = grid->igc;
  jgc  = grid->jgc;
  kgc  = grid->kgc;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  // do the first fourier transform
  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
    {
      for(int i=0; i<itot; i++)
      { 
        ijk = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        fftini[i] = p[ijk];
      }

      fftw_execute(iplanf);

      for(int i=0; i<itot; i++)
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
      for(int j=0; j<jtot; j++)
      { 
        ijk = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        fftinj[j] = p[ijk];
      }

      fftw_execute(jplanf);

      for(int j=0; j<jtot; j++)
      {
        ijk = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        p[ijk] = fftoutj[j];
      }
    }

  // TRANSPOSE

  // solve the tridiagonal system
  
  // create vectors that go into the tridiagonal matrix solver
  for(k=0; k<ktot; k++)
    for(j=0; j<jmax; j++)
      for(i=0; i<imax; i++)
      {
        // iindex = mpicoordx * imax + i
        // jindex = mpicoordy * jmax + j
        iindex = i;
        jindex = j;

        ijkb = i + j*itot + k*itot*jtot;
        ijk  = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        b[ijkb] = dz[k+kgc]*dz[k+kgc] * (bmati[iindex]+bmatj[jindex]) - (a[k]+c[k]);
        p[ijk]  = dz[k+kgc]*dz[k+kgc] * p[ijk];
      }

  for(j=0; j<jmax; j++)
    for(i=0; i<imax; i++)
    {
      iindex = i;
      jindex = j;

      // substitute BC's
      b[i+j*itot] += a[0];

      // for wave number 0, which contains average, set pressure at top to zero
      if(iindex == 0 && jindex == 0)
        b[i + j*itot + (ktot-1)*itot*jtot] -= c[ktot-1];
      // set dp/dz at top to zero
      else
        b[i + j*itot + (ktot-1)*itot*jtot] += c[ktot-1];
    }

  // call tdma solver
  tdma(a, b, c, p, work2d, work3d);
        
  // TRANSPOSE
  
  // transform the second transform back
  for(int k=0; k<kmax; k++)
    for(int i=0; i<imax; i++)
    {
      for(int j=0; j<jtot; j++)
      { 
        ijk = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        fftinj[j] = p[ijk];
      }

      fftw_execute(jplanb);

      for(int j=0; j<jtot; j++)
      {
        ijk = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        p[ijk] = fftoutj[j] / jtot;
      }
    }

  // TRANSPOSE
    
  // transform the first transform back
  for(int k=0; k<kmax; k++)
    for(int j=0; j<jmax; j++)
    {
      for(int i=0; i<itot; i++)
      { 
        ijk = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        fftini[i] = p[ijk];
      }

      fftw_execute(iplanb);

      for(int i=0; i<itot; i++)
      {
        ijk = i+igc + (j+jgc)*jj + (k+kgc)*kk;
        p[ijk] = fftouti[i] / itot;
      }
    }

  return 0;
}

int cpres::pres_2nd_out(double * restrict ut, double * restrict vt, double * restrict wt, 
                        double * restrict p , double * restrict dzhi)
{
  int    ijk,ii,jj,kk;
  double dxi,dyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ut[ijk] -= (p[ijk] - p[ijk-ii]) * dxi;
        vt[ijk] -= (p[ijk] - p[ijk-jj]) * dyi;
        wt[ijk] -= (p[ijk] - p[ijk-kk]) * dzhi[k];
      }

  return 0;
}

// tridiagonal matrix solver, taken from Numerical Recipes, Press
int cpres::tdma(double * restrict a, double * restrict b, double * restrict c, 
                double * restrict p, double * restrict work2d, double * restrict work3d)
                
{
  int i,j,k,jj,kk,ijk,ijkp,ij,jjp,kkp;
  int itot,jtot,ktot,igc,jgc,kgc;

  itot = grid->itot;
  jtot = grid->jtot;
  ktot = grid->ktot;

  igc = grid->igc;
  jgc = grid->jgc;
  kgc = grid->kgc;

  jj = itot;
  kk = itot*jtot;

  jjp = grid->icells;
  kkp = grid->icells*grid->jcells;

  for(j=0;j<jtot;j++)
    for(i=0;i<itot;i++)
    {
      ij = i + j*jj;
      work2d[ij] = b[ij];
    }

  for(j=0;j<jtot;j++)
    for(i=0;i<itot;i++)
    {
      ij   = i + j*jj;
      ijkp = i+igc + (j+jgc)*jjp + kgc*kkp;
      p[ijkp] /= work2d[ij];
    }

  for(k=1; k<ktot; k++)
  {
    for(j=0;j<jtot;j++)
      for(i=0;i<itot;i++)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + k*kk;
        work3d[ijk] = c[k-1] / work2d[ij];
      }
    for(j=0;j<jtot;j++)
      for(i=0;i<itot;i++)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + k*kk;
        work2d[ij] = b[ijk] - a[k]*work3d[ijk];
      }
    for(j=0;j<jtot;j++)
      for(i=0;i<itot;i++)
      {
        ij  = i + j*jj;
        ijk = i + j*jj + k*kk;
        ijkp = i+igc + (j+jgc)*jjp + (k+kgc)*kkp;
        p[ijkp] -= a[k]*p[ijkp-kkp];
        p[ijkp] /= work2d[ij];
      }
  }

  for(k=ktot-2; k>=0; k--)
    for(j=0;j<jtot;j++)
      for(i=0;i<itot;i++)
      {
        ijk  = i + j*jj + k*kk;
        ijkp = i+igc + (j+jgc)*jjp + (k+kgc)*kkp;
        p[ijkp] -= work3d[ijk+kk]*p[ijkp+kkp];
      }

  return 0;
}

double cpres::calcdivergence(double * restrict u, double * restrict v, double * restrict w, double * restrict dzi)
{
  int    ijk,ii,jj,kk;
  double dxi,dyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  double div, divmax;
  div    = 0;
  divmax = 0;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        div = (u[ijk+ii]-u[ijk])*dxi + (v[ijk+jj]-v[ijk])*dyi + (w[ijk+kk]-w[ijk])*dzi[k];

        divmax = std::max(divmax, std::abs(div));
      }

  return divmax;
}

