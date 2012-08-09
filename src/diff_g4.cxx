#include <cstdio>
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "diff_g4.h"
#include "defines.h"

cdiff_g4::cdiff_g4(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  std::printf("Creating instance of object diff_g4\n");
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cdiff_g4::~cdiff_g4()
{
  std::printf("Destroying instance of object diff_g4\n");
}

int cdiff_g4::diffc(double * restrict at, double * restrict a, double * restrict dzi4, double * restrict dzhi4, double visc)
{
  int    ijk,kstart,kend;
  int    ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  double dxidxi,dyidyi;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

  kstart = grid->kstart;
  kend   = grid->kend;

  dxidxi = 1./(grid->dx * grid->dx);
  dyidyi = 1./(grid->dy * grid->dy);

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + kstart*kk1;
      at[ijk] += visc * (
            + divgrad4(a[ijk-ii3], a[ijk-ii2], a[ijk-ii1], a[ijk], a[ijk+ii1], a[ijk+ii2], a[ijk+ii3], dxidxi)
            + divgrad4(a[ijk-jj3], a[ijk-jj2], a[ijk-jj1], a[ijk], a[ijk+jj1], a[ijk+jj2], a[ijk+jj3], dyidyi)
            + grad4x(grad4xbiasbot(a[ijk-kk2], a[ijk-kk1], a[ijk    ], a[ijk+kk1]) * dzhi4[kstart  ],
                     grad4x       (a[ijk-kk2], a[ijk-kk1], a[ijk    ], a[ijk+kk1]) * dzhi4[kstart+1],
                     grad4x       (a[ijk-kk1], a[ijk    ], a[ijk+kk1], a[ijk+kk2]) * dzhi4[kstart+2],
                     grad4x       (a[ijk    ], a[ijk+kk1], a[ijk+kk2], a[ijk+kk3]) * dzhi4[kstart+3])
              * dzi4[kstart] );
    }

  for(int k=grid->kstart+1; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        at[ijk] += visc * (
              + divgrad4(a[ijk-ii3], a[ijk-ii2], a[ijk-ii1], a[ijk], a[ijk+ii1], a[ijk+ii2], a[ijk+ii3], dxidxi)
              + divgrad4(a[ijk-jj3], a[ijk-jj2], a[ijk-jj1], a[ijk], a[ijk+jj1], a[ijk+jj2], a[ijk+jj3], dyidyi)
              + grad4x(grad4x(a[ijk-kk3], a[ijk-kk2], a[ijk-kk1], a[ijk    ]) * dzhi4[k-1],
                       grad4x(a[ijk-kk2], a[ijk-kk1], a[ijk    ], a[ijk+kk1]) * dzhi4[k  ],
                       grad4x(a[ijk-kk1], a[ijk    ], a[ijk+kk1], a[ijk+kk2]) * dzhi4[k+1],
                       grad4x(a[ijk    ], a[ijk+kk1], a[ijk+kk2], a[ijk+kk3]) * dzhi4[k+2])
                * dzi4[k] );
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      at[ijk] += visc * (
            + divgrad4(a[ijk-ii3], a[ijk-ii2], a[ijk-ii1], a[ijk], a[ijk+ii1], a[ijk+ii2], a[ijk+ii3], dxidxi)
            + divgrad4(a[ijk-jj3], a[ijk-jj2], a[ijk-jj1], a[ijk], a[ijk+jj1], a[ijk+jj2], a[ijk+jj3], dyidyi)
            + grad4x(grad4x       (a[ijk-kk3], a[ijk-kk2], a[ijk-kk1], a[ijk    ]) * dzhi4[kend-3],
                     grad4x       (a[ijk-kk2], a[ijk-kk1], a[ijk    ], a[ijk+kk1]) * dzhi4[kend-2],
                     grad4x       (a[ijk-kk1], a[ijk    ], a[ijk+kk1], a[ijk+kk2]) * dzhi4[kend-1],
                     grad4xbiastop(a[ijk-kk1], a[ijk    ], a[ijk+kk1], a[ijk+kk2]) * dzhi4[kend  ])
              * dzi4[kend-1] );
    }

  return 0;
}

int cdiff_g4::diffw(double * restrict at, double * restrict a, double * restrict dzi4, double * restrict dzhi4, double visc)
{
  int    ijk,kstart,kend;
  int    ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  double dxidxi,dyidyi;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

  kstart = grid->kstart;
  kend   = grid->kend;

  dxidxi = 1./(grid->dx * grid->dx);
  dyidyi = 1./(grid->dy * grid->dy);

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kstart+1)*kk1;
      at[ijk] += visc * (
            + divgrad4(a[ijk-ii3], a[ijk-ii2], a[ijk-ii1], a[ijk], a[ijk+ii1], a[ijk+ii2], a[ijk+ii3], dxidxi)
            + divgrad4(a[ijk-jj3], a[ijk-jj2], a[ijk-jj1], a[ijk], a[ijk+jj1], a[ijk+jj2], a[ijk+jj3], dyidyi)
            + grad4x(grad4xbiasbot(a[ijk-kk2], a[ijk-kk1], a[ijk    ], a[ijk+kk1]) * dzi4[kstart  ],
                     grad4x       (a[ijk-kk2], a[ijk-kk1], a[ijk    ], a[ijk+kk1]) * dzi4[kstart+1],
                     grad4x       (a[ijk-kk1], a[ijk    ], a[ijk+kk1], a[ijk+kk2]) * dzi4[kstart+2],
                     grad4x       (a[ijk    ], a[ijk+kk1], a[ijk+kk2], a[ijk+kk3]) * dzi4[kstart+3])
              * dzhi4[kstart+1] );
    }

  for(int k=grid->kstart+2; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        at[ijk] += visc * (
              + divgrad4(a[ijk-ii3], a[ijk-ii2], a[ijk-ii1], a[ijk], a[ijk+ii1], a[ijk+ii2], a[ijk+ii3], dxidxi)
              + divgrad4(a[ijk-jj3], a[ijk-jj2], a[ijk-jj1], a[ijk], a[ijk+jj1], a[ijk+jj2], a[ijk+jj3], dyidyi)
              + grad4x(grad4x(a[ijk-kk3], a[ijk-kk2], a[ijk-kk1], a[ijk    ]) * dzi4[k-2],
                       grad4x(a[ijk-kk2], a[ijk-kk1], a[ijk    ], a[ijk+kk1]) * dzi4[k-1],
                       grad4x(a[ijk-kk1], a[ijk    ], a[ijk+kk1], a[ijk+kk2]) * dzi4[k  ],
                       grad4x(a[ijk    ], a[ijk+kk1], a[ijk+kk2], a[ijk+kk3]) * dzi4[k+1])
                * dzhi4[k] );
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      at[ijk] += visc * (
            + divgrad4(a[ijk-ii3], a[ijk-ii2], a[ijk-ii1], a[ijk], a[ijk+ii1], a[ijk+ii2], a[ijk+ii3], dxidxi)
            + divgrad4(a[ijk-jj3], a[ijk-jj2], a[ijk-jj1], a[ijk], a[ijk+jj1], a[ijk+jj2], a[ijk+jj3], dyidyi)
            + grad4x(grad4x       (a[ijk-kk3], a[ijk-kk2], a[ijk-kk1], a[ijk    ]) * dzi4[kend-4],
                     grad4x       (a[ijk-kk2], a[ijk-kk1], a[ijk    ], a[ijk+kk1]) * dzi4[kend-3],
                     grad4x       (a[ijk-kk1], a[ijk    ], a[ijk+kk1], a[ijk+kk2]) * dzi4[kend-2],
                     grad4xbiastop(a[ijk-kk1], a[ijk    ], a[ijk+kk1], a[ijk+kk2]) * dzi4[kend-1])
              * dzhi4[kend-1] );
    }

  return 0;
}

inline double cdiff_g4::divgrad4(const double am3, const double am2, const double am1, const double a,
                                 const double ap1, const double ap2, const double ap3, const double dxidxi)
{
  return ( (1./576.)*(am3+ap3) - (54./576.)*(am2+ap2) + (783./576.)*(am1+ap1) - (1460./576.)*a ) * dxidxi;
}

inline double cdiff_g4::grad4x(const double a, const double b, const double c, const double d)
{
  return (-(d-a) + 27.*(c-b)); 
}

inline double cdiff_g4::grad4xbiasbot(const double a, const double b, const double c, const double d)
{
  return (-23.*a + 21.*b + 3.*c - d);
}

inline double cdiff_g4::grad4xbiastop(const double a, const double b, const double c, const double d)
{
  return ( 23.*d - 21.*c - 3.*b + a);
}

