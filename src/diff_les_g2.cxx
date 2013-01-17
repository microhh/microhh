#include <cstdio>
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "diff_les_g2.h"
#include "defines.h"

cdiff_les_g2::cdiff_les_g2(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  // std::printf("Creating instance of object diff_les_g2\n");
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cdiff_les_g2::~cdiff_les_g2()
{
  // std::printf("Destroying instance of object diff_les_g2\n");
}

int cdiff_les_g2::diffu(double * restrict ut, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi, double * restrict dzhi, double visc)
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
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        ut[ijk] += visc * (
              // du/dx + du/dx
              + (  (u[ijk+ii]-u[ijk   ])*dxi
                 - (u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
              // du/dy + dv/dx
              + (  ((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                 - ((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
              // du/dz + dw/dx
              + (  ((u[ijk+kk]-u[ijk   ])* dzhi[k+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
                 - ((u[ijk   ]-u[ijk-kk])* dzhi[k  ] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) * dzi[k] );
      }

  return 0;
}

int cdiff_les_g2::diffv(double * restrict vt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi, double * restrict dzhi, double visc)
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
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        vt[ijk] += visc * (
              // dv/dx + du/dy
              + (  ((v[ijk+ii]-v[ijk   ])*dxi  + (u[ijk+jj]-u[ijk-ii+jj])*dyi)
                 - ((v[ijk   ]-v[ijk-ii])*dxi  + (u[ijk   ]-u[ijk-ii   ])*dyi) ) * dxi
              // dv/dy + dv/dy
              + (  (v[ijk+jj]-v[ijk   ])*dyi
                 - (v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
              // dv/dz + dw/dy
              + (  ((v[ijk+kk]-v[ijk   ])*dzhi[k+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
                 - ((v[ijk   ]-v[ijk-kk])*dzhi[k  ] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) * dzi[k] );
      }

  return 0;
}

int cdiff_les_g2::diffw(double * restrict wt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi, double * restrict dzhi, double visc)
{
  int    ijk,ii,jj,kk;
  double dxi,dyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        wt[ijk] += visc * (
              // dw/dx + du/dz
              + (  ((w[ijk+ii]-w[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-kk])*dzhi[k])
                 - ((w[ijk   ]-w[ijk-ii])*dxi + (u[ijk   ]-u[ijk+  -kk])*dzhi[k]) ) * dxi
              // dw/dy + dv/dz
              + (  ((w[ijk+jj]-w[ijk   ])*dyi + (v[ijk+jj]-v[ijk+jj-kk])*dzhi[k])
                 - ((w[ijk   ]-w[ijk-jj])*dyi + (v[ijk   ]-v[ijk+  -kk])*dzhi[k]) ) * dyi
              // dw/dz + dw/dz
              + (  (w[ijk+kk]-w[ijk   ])*dzi[k]
                 - (w[ijk   ]-w[ijk-kk])*dzi[k-1] ) * 2.* dzhi[k] );
      }

  return 0;
}

int cdiff_les_g2::diffc(double * restrict at, double * restrict a, double * restrict dzi, double * restrict dzhi, double visc)
{
  int    ijk,ii,jj,kk;
  double dxidxi,dyidyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dxidxi = 1./(grid->dx * grid->dx);
  dyidyi = 1./(grid->dy * grid->dy);

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        at[ijk] += visc * (
              + (  (a[ijk+ii]-a[ijk   ]) 
                 - (a[ijk   ]-a[ijk-ii]) )*dxidxi 
              + (  (a[ijk+jj]-a[ijk   ]) 
                 - (a[ijk   ]-a[ijk-jj]) )*dyidyi
              + (  (a[ijk+kk]-a[ijk   ])*dzhi[k+1]
                 - (a[ijk   ]-a[ijk-kk])*dzhi[k]  ) * dzi[k] );
      }

  return 0;
}

