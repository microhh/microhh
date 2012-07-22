#include <cstdio>
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "diff_g42.h"
#include "defines.h"

cdiff_g42::cdiff_g42(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  std::printf("Creating instance of object diff_g42\n");
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cdiff_g42::~cdiff_g42()
{
  std::printf("Destroying instance of object diff_g42\n");
}

int cdiff_g42::diffc(double * restrict at, double * restrict a, double * restrict dzi, double * restrict dzhi, double visc)
{
  int    ijk,kk;
  int    ii1,ii2,ii3,jj1,jj2,jj3;
  double dxidxi,dyidyi;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk = grid->icells*grid->jcells;

  dxidxi = 1./(grid->dx * grid->dx);
  dyidyi = 1./(grid->dy * grid->dy);

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk;
        at[ijk] += visc * (
              + divgrad4(a[ijk-ii3], a[ijk-ii2], a[ijk-ii1], a[ijk], a[ijk+ii1], a[ijk+ii2], a[ijk+ii3], dxidxi)
              + divgrad4(a[ijk-jj3], a[ijk-jj2], a[ijk-jj1], a[ijk], a[ijk+jj1], a[ijk+jj2], a[ijk+jj3], dyidyi) 
              + (  (a[ijk+kk] - a[ijk   ]) * dzhi[k+1]
                 - (a[ijk   ] - a[ijk-kk]) * dzhi[k]   ) * dzi[k] );
      }

  return 0;
}

int cdiff_g42::diffw(double * restrict wt, double * restrict w, double * restrict dzi, double * restrict dzhi, double visc)
{
  int    ijk,kk;
  int    ii1,ii2,ii3,jj1,jj2,jj3;
  double dxidxi,dyidyi;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk = grid->icells*grid->jcells;

  dxidxi = 1./(grid->dx*grid->dx);
  dyidyi = 1./(grid->dy*grid->dy);

  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk;
        wt[ijk] += visc * (
              + divgrad4(w[ijk-ii3], w[ijk-ii2], w[ijk-ii1], w[ijk], w[ijk+ii1], w[ijk+ii2], w[ijk+ii3], dxidxi) 
              + divgrad4(w[ijk-jj3], w[ijk-jj2], w[ijk-jj1], w[ijk], w[ijk+jj1], w[ijk+jj2], w[ijk+jj3], dyidyi) 
              + (  (w[ijk+kk] - w[ijk   ]) * dzi[k]
                 - (w[ijk   ] - w[ijk-kk]) * dzi[k-1] ) * dzhi[k] );
      }

  return 0;
}

inline double cdiff_g42::divgrad4(const double am3, const double am2, const double am1, const double a,
                                  const double ap1, const double ap2, const double ap3, const double dxidxi)
{
  return ( (1./576.)*(am3+ap3) - (54./576.)*(am2+ap2) + (783./576.)*(am1+ap1) - (1460./576.)*a ) * dxidxi;
}
