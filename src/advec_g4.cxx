#include <cstdio>
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "advec_g4.h"
#include "defines.h"

cadvec_g4::cadvec_g4(cgrid *gridin, cfields *fieldsin, cmpi *mpiin)
{
  // std::printf("Creating instance of object advec_g4\n");
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cadvec_g4::~cadvec_g4()
{
  // std::printf("Destroying instance of object advec_g4\n");
}

double cadvec_g4::calccfl(double * restrict u, double * restrict v, double * restrict w, double * restrict dzi, double dt)
{
  int    ijk;
  int    ii1,ii2,jj1,jj2,kk1,kk2;
  double dxi,dyi;

  const double ci0 = -1./16.;
  const double ci1 =  9./16.;
  const double ci2 =  9./16.;
  const double ci3 = -1./16.;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  double cfl = 0;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj1 + k*kk1;
        cfl = std::max(cfl, std::abs(interp4(u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2]))*dxi 
                          + std::abs(interp4(v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2]))*dyi 
                          + std::abs(interp4(w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2]))*dzi[k] );
      }

  grid->getmax(&cfl);

  cfl = cfl*dt;

  return cfl;
}

int cadvec_g4::advecu(double * restrict ut, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi4)
{
  int    ijk,kstart,kend;
  int    ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  double dxi,dyi;

  const double cg0 =   1.;
  const double cg1 = -27.;
  const double cg2 =  27.;
  const double cg3 =  -1.;
  const double cgi =   1./24.;

  const double ci0 = -1./16.;
  const double ci1 =  9./16.;
  const double ci2 =  9./16.;
  const double ci3 = -1./16.;

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

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + kstart*kk1;
      ut[ijk] -= ( cg0*((ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]) * (ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]))
                 + cg1*((ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]) * (ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]))
                 + cg2*((ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]) * (ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]))
                 + cg3*((ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3]) * (ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3])) ) * cgi*dxi;

      ut[ijk] -= ( cg0*((ci0*v[ijk-ii2-jj1] + ci1*v[ijk-ii1-jj1] + ci2*v[ijk-jj1] + ci3*v[ijk+ii1-jj1]) * (ci0*u[ijk-jj3] + ci1*u[ijk-jj2] + ci2*u[ijk-jj1] + ci3*u[ijk    ]))
                 + cg1*((ci0*v[ijk-ii2    ] + ci1*v[ijk-ii1    ] + ci2*v[ijk    ] + ci3*v[ijk+ii1    ]) * (ci0*u[ijk-jj2] + ci1*u[ijk-jj1] + ci2*u[ijk    ] + ci3*u[ijk+jj1]))
                 + cg2*((ci0*v[ijk-ii2+jj1] + ci1*v[ijk-ii1+jj1] + ci2*v[ijk+jj1] + ci3*v[ijk+ii1+jj1]) * (ci0*u[ijk-jj1] + ci1*u[ijk    ] + ci2*u[ijk+jj1] + ci3*u[ijk+jj2]))
                 + cg3*((ci0*v[ijk-ii2+jj2] + ci1*v[ijk-ii1+jj2] + ci2*v[ijk+jj2] + ci3*v[ijk+ii1+jj2]) * (ci0*u[ijk    ] + ci1*u[ijk+jj1] + ci2*u[ijk+jj2] + ci3*u[ijk+jj3])) ) * cgi*dyi;

      ut[ijk] -= ( cg0*(interp4(w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1]) * interp4biasbot(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]))
                 + cg1*(interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp4       (u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]))
                 + cg2*(interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp4       (u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]))
                 + cg3*(interp4(w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2]) * interp4       (u[ijk    ], u[ijk+kk1], u[ijk+kk2], u[ijk+kk3])) )
                 * dzi4[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        ut[ijk] -= ( cg0*((ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]) * (ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]))
                   + cg1*((ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]) * (ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]))
                   + cg2*((ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]) * (ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]))
                   + cg3*((ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3]) * (ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3])) ) * cgi*dxi;

        ut[ijk] -= ( cg0*((ci0*v[ijk-ii2-jj1] + ci1*v[ijk-ii1-jj1] + ci2*v[ijk-jj1] + ci3*v[ijk+ii1-jj1]) * (ci0*u[ijk-jj3] + ci1*u[ijk-jj2] + ci2*u[ijk-jj1] + ci3*u[ijk    ]))
                   + cg1*((ci0*v[ijk-ii2    ] + ci1*v[ijk-ii1    ] + ci2*v[ijk    ] + ci3*v[ijk+ii1    ]) * (ci0*u[ijk-jj2] + ci1*u[ijk-jj1] + ci2*u[ijk    ] + ci3*u[ijk+jj1]))
                   + cg2*((ci0*v[ijk-ii2+jj1] + ci1*v[ijk-ii1+jj1] + ci2*v[ijk+jj1] + ci3*v[ijk+ii1+jj1]) * (ci0*u[ijk-jj1] + ci1*u[ijk    ] + ci2*u[ijk+jj1] + ci3*u[ijk+jj2]))
                   + cg3*((ci0*v[ijk-ii2+jj2] + ci1*v[ijk-ii1+jj2] + ci2*v[ijk+jj2] + ci3*v[ijk+ii1+jj2]) * (ci0*u[ijk    ] + ci1*u[ijk+jj1] + ci2*u[ijk+jj2] + ci3*u[ijk+jj3])) ) * cgi*dyi;

        ut[ijk] -= ( cg0*(interp4(w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1]) * interp4(u[ijk-kk3], u[ijk-kk2], u[ijk-kk1], u[ijk    ]))
                   + cg1*(interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp4(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]))
                   + cg2*(interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp4(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]))
                   + cg3*(interp4(w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2]) * interp4(u[ijk    ], u[ijk+kk1], u[ijk+kk2], u[ijk+kk3])) )
                   * dzi4[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      ut[ijk] -= ( cg0*((ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]) * (ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]))
                 + cg1*((ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]) * (ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]))
                 + cg2*((ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]) * (ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]))
                 + cg3*((ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3]) * (ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3])) ) * cgi*dxi;

      ut[ijk] -= ( cg0*((ci0*v[ijk-ii2-jj1] + ci1*v[ijk-ii1-jj1] + ci2*v[ijk-jj1] + ci3*v[ijk+ii1-jj1]) * (ci0*u[ijk-jj3] + ci1*u[ijk-jj2] + ci2*u[ijk-jj1] + ci3*u[ijk    ]))
                 + cg1*((ci0*v[ijk-ii2    ] + ci1*v[ijk-ii1    ] + ci2*v[ijk    ] + ci3*v[ijk+ii1    ]) * (ci0*u[ijk-jj2] + ci1*u[ijk-jj1] + ci2*u[ijk    ] + ci3*u[ijk+jj1]))
                 + cg2*((ci0*v[ijk-ii2+jj1] + ci1*v[ijk-ii1+jj1] + ci2*v[ijk+jj1] + ci3*v[ijk+ii1+jj1]) * (ci0*u[ijk-jj1] + ci1*u[ijk    ] + ci2*u[ijk+jj1] + ci3*u[ijk+jj2]))
                 + cg3*((ci0*v[ijk-ii2+jj2] + ci1*v[ijk-ii1+jj2] + ci2*v[ijk+jj2] + ci3*v[ijk+ii1+jj2]) * (ci0*u[ijk    ] + ci1*u[ijk+jj1] + ci2*u[ijk+jj2] + ci3*u[ijk+jj3])) ) * cgi*dyi;

      ut[ijk] -= ( cg0*(interp4(w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1]) * interp4       (u[ijk-kk3], u[ijk-kk2], u[ijk-kk1], u[ijk    ]))
                 + cg1*(interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp4       (u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]))
                 + cg2*(interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp4       (u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2]))
                 + cg3*(interp4(w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2]) * interp4biastop(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2])) )
                 * dzi4[kend-1];
    }

  return 0;
}

int cadvec_g4::advecv(double * restrict vt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi4)
{
  int    ijk,kstart,kend;
  int    ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  double dxi,dyi;

  const double cg0 =   1.;
  const double cg1 = -27.;
  const double cg2 =  27.;
  const double cg3 =  -1.;
  const double cgi =   1./24.;

  const double ci0 = -1./16.;
  const double ci1 =  9./16.;
  const double ci2 =  9./16.;
  const double ci3 = -1./16.;

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

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + kstart*kk1;
      vt[ijk] -= ( cg0*((ci0*u[ijk-ii1-jj2] + ci1*u[ijk-ii1-jj1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+jj1]) * (ci0*v[ijk-ii3] + ci1*v[ijk-ii2] + ci2*v[ijk-ii1] + ci3*v[ijk    ]))
                 + cg1*((ci0*u[ijk    -jj2] + ci1*u[ijk    -jj1] + ci2*u[ijk    ] + ci3*u[ijk    +jj1]) * (ci0*v[ijk-ii2] + ci1*v[ijk-ii1] + ci2*v[ijk    ] + ci3*v[ijk+ii1]))
                 + cg2*((ci0*u[ijk+ii1-jj2] + ci1*u[ijk+ii1-jj1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+jj1]) * (ci0*v[ijk-ii1] + ci1*v[ijk    ] + ci2*v[ijk+ii1] + ci3*v[ijk+ii2]))
                 + cg3*((ci0*u[ijk+ii2-jj2] + ci1*u[ijk+ii2-jj1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+jj1]) * (ci0*v[ijk    ] + ci1*v[ijk+ii1] + ci2*v[ijk+ii2] + ci3*v[ijk+ii3])) ) * cgi*dxi;

      vt[ijk] -= ( cg0*((ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]) * (ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]))
                 + cg1*((ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]) * (ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]))
                 + cg2*((ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]) * (ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]))
                 + cg3*((ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3]) * (ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3])) ) * cgi*dyi;

      vt[ijk] -= ( cg0*(interp4(w[ijk-jj2-kk1], w[ijk-jj1-kk1], w[ijk-kk1], w[ijk+jj1-kk1]) * interp4biasbot(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]))
                 + cg1*(interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp4       (v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]))
                 + cg2*(interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp4       (v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2]))
                 + cg3*(interp4(w[ijk-jj2+kk2], w[ijk-jj1+kk2], w[ijk+kk2], w[ijk+jj1+kk2]) * interp4       (v[ijk    ], v[ijk+kk1], v[ijk+kk2], v[ijk+kk3])) )
                 * dzi4[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        vt[ijk] -= ( cg0*((ci0*u[ijk-ii1-jj2] + ci1*u[ijk-ii1-jj1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+jj1]) * (ci0*v[ijk-ii3] + ci1*v[ijk-ii2] + ci2*v[ijk-ii1] + ci3*v[ijk    ]))
                   + cg1*((ci0*u[ijk    -jj2] + ci1*u[ijk    -jj1] + ci2*u[ijk    ] + ci3*u[ijk    +jj1]) * (ci0*v[ijk-ii2] + ci1*v[ijk-ii1] + ci2*v[ijk    ] + ci3*v[ijk+ii1]))
                   + cg2*((ci0*u[ijk+ii1-jj2] + ci1*u[ijk+ii1-jj1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+jj1]) * (ci0*v[ijk-ii1] + ci1*v[ijk    ] + ci2*v[ijk+ii1] + ci3*v[ijk+ii2]))
                   + cg3*((ci0*u[ijk+ii2-jj2] + ci1*u[ijk+ii2-jj1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+jj1]) * (ci0*v[ijk    ] + ci1*v[ijk+ii1] + ci2*v[ijk+ii2] + ci3*v[ijk+ii3])) ) * cgi*dxi;

        vt[ijk] -= ( cg0*((ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]) * (ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]))
                   + cg1*((ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]) * (ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]))
                   + cg2*((ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]) * (ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]))
                   + cg3*((ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3]) * (ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3])) ) * cgi*dyi;

        vt[ijk] -= ( cg0*(interp4(w[ijk-jj2-kk1], w[ijk-jj1-kk1], w[ijk-kk1], w[ijk+jj1-kk1]) * interp4(v[ijk-kk3], v[ijk-kk2], v[ijk-kk1], v[ijk    ]))
                   + cg1*(interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp4(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]))
                   + cg2*(interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp4(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2]))
                   + cg3*(interp4(w[ijk-jj2+kk2], w[ijk-jj1+kk2], w[ijk+kk2], w[ijk+jj1+kk2]) * interp4(v[ijk    ], v[ijk+kk1], v[ijk+kk2], v[ijk+kk3])) )
                   * dzi4[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      vt[ijk] -= ( cg0*((ci0*u[ijk-ii1-jj2] + ci1*u[ijk-ii1-jj1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+jj1]) * (ci0*v[ijk-ii3] + ci1*v[ijk-ii2] + ci2*v[ijk-ii1] + ci3*v[ijk    ]))
                 + cg1*((ci0*u[ijk    -jj2] + ci1*u[ijk    -jj1] + ci2*u[ijk    ] + ci3*u[ijk    +jj1]) * (ci0*v[ijk-ii2] + ci1*v[ijk-ii1] + ci2*v[ijk    ] + ci3*v[ijk+ii1]))
                 + cg2*((ci0*u[ijk+ii1-jj2] + ci1*u[ijk+ii1-jj1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+jj1]) * (ci0*v[ijk-ii1] + ci1*v[ijk    ] + ci2*v[ijk+ii1] + ci3*v[ijk+ii2]))
                 + cg3*((ci0*u[ijk+ii2-jj2] + ci1*u[ijk+ii2-jj1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+jj1]) * (ci0*v[ijk    ] + ci1*v[ijk+ii1] + ci2*v[ijk+ii2] + ci3*v[ijk+ii3])) ) * cgi*dxi;

      vt[ijk] -= ( cg0*((ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]) * (ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]))
                 + cg1*((ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]) * (ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]))
                 + cg2*((ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]) * (ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]))
                 + cg3*((ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3]) * (ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3])) ) * cgi*dyi;

      vt[ijk] -= ( cg0*(interp4(w[ijk-jj2-kk1], w[ijk-jj1-kk1], w[ijk-kk1], w[ijk+jj1-kk1]) * interp4       (v[ijk-kk3], v[ijk-kk2], v[ijk-kk1], v[ijk    ]))
                 + cg1*(interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp4       (v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]))
                 + cg2*(interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp4       (v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2]))
                 + cg3*(interp4(w[ijk-jj2+kk2], w[ijk-jj1+kk2], w[ijk+kk2], w[ijk+jj1+kk2]) * interp4biastop(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])) )
                 * dzi4[kend-1];
    }

  return 0;
}

int cadvec_g4::advecw(double * restrict wt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzhi4)
{
  int    ijk,kstart,kend;
  int    ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  double dxi,dyi;

  const double cg0 =   1.;
  const double cg1 = -27.;
  const double cg2 =  27.;
  const double cg3 =  -1.;
  const double cgi =   1./24.;

  const double ci0 = -1./16.;
  const double ci1 =  9./16.;
  const double ci2 =  9./16.;
  const double ci3 = -1./16.;

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

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kstart+1)*kk1;
      wt[ijk] -= ( cg0*((ci0*u[ijk-ii1-kk2] + ci1*u[ijk-ii1-kk1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+kk1]) * (ci0*w[ijk-ii3] + ci1*w[ijk-ii2] + ci2*w[ijk-ii1] + ci3*w[ijk    ]))
                 + cg1*((ci0*u[ijk    -kk2] + ci1*u[ijk    -kk1] + ci2*u[ijk    ] + ci3*u[ijk    +kk1]) * (ci0*w[ijk-ii2] + ci1*w[ijk-ii1] + ci2*w[ijk    ] + ci3*w[ijk+ii1]))
                 + cg2*((ci0*u[ijk+ii1-kk2] + ci1*u[ijk+ii1-kk1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+kk1]) * (ci0*w[ijk-ii1] + ci1*w[ijk    ] + ci2*w[ijk+ii1] + ci3*w[ijk+ii2]))
                 + cg3*((ci0*u[ijk+ii2-kk2] + ci1*u[ijk+ii2-kk1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+ii1] + ci2*w[ijk+ii2] + ci3*w[ijk+ii3])) ) * cgi*dxi;

      wt[ijk] -= ( cg0*((ci0*v[ijk-jj1-kk2] + ci1*v[ijk-jj1-kk1] + ci2*v[ijk-jj1] + ci3*v[ijk-jj1+kk1]) * (ci0*w[ijk-jj3] + ci1*w[ijk-jj2] + ci2*w[ijk-jj1] + ci3*w[ijk    ]))
                 + cg1*((ci0*v[ijk    -kk2] + ci1*v[ijk    -kk1] + ci2*v[ijk    ] + ci3*v[ijk    +kk1]) * (ci0*w[ijk-jj2] + ci1*w[ijk-jj1] + ci2*w[ijk    ] + ci3*w[ijk+jj1]))
                 + cg2*((ci0*v[ijk+jj1-kk2] + ci1*v[ijk+jj1-kk1] + ci2*v[ijk+jj1] + ci3*v[ijk+jj1+kk1]) * (ci0*w[ijk-jj1] + ci1*w[ijk    ] + ci2*w[ijk+jj1] + ci3*w[ijk+jj2]))
                 + cg3*((ci0*v[ijk+jj2-kk2] + ci1*v[ijk+jj2-kk1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+jj1] + ci2*w[ijk+jj2] + ci3*w[ijk+jj3])) ) * cgi*dyi;

      wt[ijk] -= ( cg0*(interp4biasbot(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) * interp4biasbot(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]))
                 + cg1*(interp4       (w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) * interp4       (w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]))
                 + cg2*(interp4       (w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) * interp4       (w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]))
                 + cg3*(interp4       (w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3]) * interp4       (w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3])) )
                 * dzhi4[kstart+1];
    }

  for(int k=grid->kstart+2; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        wt[ijk] -= ( cg0*((ci0*u[ijk-ii1-kk2] + ci1*u[ijk-ii1-kk1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+kk1]) * (ci0*w[ijk-ii3] + ci1*w[ijk-ii2] + ci2*w[ijk-ii1] + ci3*w[ijk    ]))
                   + cg1*((ci0*u[ijk    -kk2] + ci1*u[ijk    -kk1] + ci2*u[ijk    ] + ci3*u[ijk    +kk1]) * (ci0*w[ijk-ii2] + ci1*w[ijk-ii1] + ci2*w[ijk    ] + ci3*w[ijk+ii1]))
                   + cg2*((ci0*u[ijk+ii1-kk2] + ci1*u[ijk+ii1-kk1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+kk1]) * (ci0*w[ijk-ii1] + ci1*w[ijk    ] + ci2*w[ijk+ii1] + ci3*w[ijk+ii2]))
                   + cg3*((ci0*u[ijk+ii2-kk2] + ci1*u[ijk+ii2-kk1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+ii1] + ci2*w[ijk+ii2] + ci3*w[ijk+ii3])) ) * cgi*dxi;

        wt[ijk] -= ( cg0*((ci0*v[ijk-jj1-kk2] + ci1*v[ijk-jj1-kk1] + ci2*v[ijk-jj1] + ci3*v[ijk-jj1+kk1]) * (ci0*w[ijk-jj3] + ci1*w[ijk-jj2] + ci2*w[ijk-jj1] + ci3*w[ijk    ]))
                   + cg1*((ci0*v[ijk    -kk2] + ci1*v[ijk    -kk1] + ci2*v[ijk    ] + ci3*v[ijk    +kk1]) * (ci0*w[ijk-jj2] + ci1*w[ijk-jj1] + ci2*w[ijk    ] + ci3*w[ijk+jj1]))
                   + cg2*((ci0*v[ijk+jj1-kk2] + ci1*v[ijk+jj1-kk1] + ci2*v[ijk+jj1] + ci3*v[ijk+jj1+kk1]) * (ci0*w[ijk-jj1] + ci1*w[ijk    ] + ci2*w[ijk+jj1] + ci3*w[ijk+jj2]))
                   + cg3*((ci0*v[ijk+jj2-kk2] + ci1*v[ijk+jj2-kk1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+jj1] + ci2*w[ijk+jj2] + ci3*w[ijk+jj3])) ) * cgi*dyi;

        wt[ijk] -= ( cg0*(interp4(w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ]) * interp4(w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ]))
                   + cg1*(interp4(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) * interp4(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]))
                   + cg2*(interp4(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) * interp4(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]))
                   + cg3*(interp4(w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3]) * interp4(w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3])) )
                   * dzhi4[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      wt[ijk] -= ( cg0*((ci0*u[ijk-ii1-kk2] + ci1*u[ijk-ii1-kk1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+kk1]) * (ci0*w[ijk-ii3] + ci1*w[ijk-ii2] + ci2*w[ijk-ii1] + ci3*w[ijk    ]))
                 + cg1*((ci0*u[ijk    -kk2] + ci1*u[ijk    -kk1] + ci2*u[ijk    ] + ci3*u[ijk    +kk1]) * (ci0*w[ijk-ii2] + ci1*w[ijk-ii1] + ci2*w[ijk    ] + ci3*w[ijk+ii1]))
                 + cg2*((ci0*u[ijk+ii1-kk2] + ci1*u[ijk+ii1-kk1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+kk1]) * (ci0*w[ijk-ii1] + ci1*w[ijk    ] + ci2*w[ijk+ii1] + ci3*w[ijk+ii2]))
                 + cg3*((ci0*u[ijk+ii2-kk2] + ci1*u[ijk+ii2-kk1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+ii1] + ci2*w[ijk+ii2] + ci3*w[ijk+ii3])) ) * cgi*dxi;

      wt[ijk] -= ( cg0*((ci0*v[ijk-jj1-kk2] + ci1*v[ijk-jj1-kk1] + ci2*v[ijk-jj1] + ci3*v[ijk-jj1+kk1]) * (ci0*w[ijk-jj3] + ci1*w[ijk-jj2] + ci2*w[ijk-jj1] + ci3*w[ijk    ]))
                 + cg1*((ci0*v[ijk    -kk2] + ci1*v[ijk    -kk1] + ci2*v[ijk    ] + ci3*v[ijk    +kk1]) * (ci0*w[ijk-jj2] + ci1*w[ijk-jj1] + ci2*w[ijk    ] + ci3*w[ijk+jj1]))
                 + cg2*((ci0*v[ijk+jj1-kk2] + ci1*v[ijk+jj1-kk1] + ci2*v[ijk+jj1] + ci3*v[ijk+jj1+kk1]) * (ci0*w[ijk-jj1] + ci1*w[ijk    ] + ci2*w[ijk+jj1] + ci3*w[ijk+jj2]))
                 + cg3*((ci0*v[ijk+jj2-kk2] + ci1*v[ijk+jj2-kk1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+jj1] + ci2*w[ijk+jj2] + ci3*w[ijk+jj3])) ) * cgi*dyi;

      wt[ijk] -= ( cg0*(interp4       (w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ]) * interp4       (w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ]))
                 + cg1*(interp4       (w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) * interp4       (w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]))
                 + cg2*(interp4       (w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) * interp4       (w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]))
                 + cg3*(interp4biastop(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) * interp4biastop(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2])) )
                 * dzhi4[kend-1];
    }

  return 0;
}

int cadvec_g4::advecs(double * restrict st, double * restrict s, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi4)
{
  int    ijk,kstart,kend;
  int    ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  double dxi,dyi;

  const double cg0 =   1.;
  const double cg1 = -27.;
  const double cg2 =  27.;
  const double cg3 =  -1.;
  const double cgi =   1./24.;

  const double ci0 = -1./16.;
  const double ci1 =  9./16.;
  const double ci2 =  9./16.;
  const double ci3 = -1./16.;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk1 = 1*grid->icells*grid->jcells;
  kk2 = 2*grid->icells*grid->jcells;
  kk3 = 3*grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  kstart = grid->kstart;
  kend   = grid->kend;

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + kstart*kk1;
      st[ijk] -= ( cg0*(u[ijk-ii1] * (ci0*s[ijk-ii3] + ci1*s[ijk-ii2] + ci2*s[ijk-ii1] + ci3*s[ijk    ]))
                 + cg1*(u[ijk    ] * (ci0*s[ijk-ii2] + ci1*s[ijk-ii1] + ci2*s[ijk    ] + ci3*s[ijk+ii1]))
                 + cg2*(u[ijk+ii1] * (ci0*s[ijk-ii1] + ci1*s[ijk    ] + ci2*s[ijk+ii1] + ci3*s[ijk+ii2]))
                 + cg3*(u[ijk+ii2] * (ci0*s[ijk    ] + ci1*s[ijk+ii1] + ci2*s[ijk+ii2] + ci3*s[ijk+ii3])) ) * cgi*dxi;

      st[ijk] -= ( cg0*(v[ijk-jj1] * (ci0*s[ijk-jj3] + ci1*s[ijk-jj2] + ci2*s[ijk-jj1] + ci3*s[ijk    ]))
                 + cg1*(v[ijk    ] * (ci0*s[ijk-jj2] + ci1*s[ijk-jj1] + ci2*s[ijk    ] + ci3*s[ijk+jj1]))
                 + cg2*(v[ijk+jj1] * (ci0*s[ijk-jj1] + ci1*s[ijk    ] + ci2*s[ijk+jj1] + ci3*s[ijk+jj2]))
                 + cg3*(v[ijk+jj2] * (ci0*s[ijk    ] + ci1*s[ijk+jj1] + ci2*s[ijk+jj2] + ci3*s[ijk+jj3])) ) * cgi*dyi;

      st[ijk] -= ( cg0*(w[ijk-kk1] * interp4biasbot(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]))
                 + cg1*(w[ijk    ] * interp4       (s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]))
                 + cg2*(w[ijk+kk1] * interp4       (s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2]))
                 + cg3*(w[ijk+kk2] * interp4       (s[ijk    ], s[ijk+kk1], s[ijk+kk2], s[ijk+kk3])) )
                 * dzi4[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        st[ijk] -= ( cg0*(u[ijk-ii1] * (ci0*s[ijk-ii3] + ci1*s[ijk-ii2] + ci2*s[ijk-ii1] + ci3*s[ijk    ]))
                   + cg1*(u[ijk    ] * (ci0*s[ijk-ii2] + ci1*s[ijk-ii1] + ci2*s[ijk    ] + ci3*s[ijk+ii1]))
                   + cg2*(u[ijk+ii1] * (ci0*s[ijk-ii1] + ci1*s[ijk    ] + ci2*s[ijk+ii1] + ci3*s[ijk+ii2]))
                   + cg3*(u[ijk+ii2] * (ci0*s[ijk    ] + ci1*s[ijk+ii1] + ci2*s[ijk+ii2] + ci3*s[ijk+ii3])) ) * cgi*dxi;

        st[ijk] -= ( cg0*(v[ijk-jj1] * (ci0*s[ijk-jj3] + ci1*s[ijk-jj2] + ci2*s[ijk-jj1] + ci3*s[ijk    ]))
                   + cg1*(v[ijk    ] * (ci0*s[ijk-jj2] + ci1*s[ijk-jj1] + ci2*s[ijk    ] + ci3*s[ijk+jj1]))
                   + cg2*(v[ijk+jj1] * (ci0*s[ijk-jj1] + ci1*s[ijk    ] + ci2*s[ijk+jj1] + ci3*s[ijk+jj2]))
                   + cg3*(v[ijk+jj2] * (ci0*s[ijk    ] + ci1*s[ijk+jj1] + ci2*s[ijk+jj2] + ci3*s[ijk+jj3])) ) * cgi*dyi;

        st[ijk] -= ( cg0*(w[ijk-kk1] * interp4(s[ijk-kk3], s[ijk-kk2], s[ijk-kk1], s[ijk    ]))
                   + cg1*(w[ijk    ] * interp4(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]))
                   + cg2*(w[ijk+kk1] * interp4(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2]))
                   + cg3*(w[ijk+kk2] * interp4(s[ijk    ], s[ijk+kk1], s[ijk+kk2], s[ijk+kk3])) )
                   * dzi4[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ijk = i + j*jj1 + (kend-1)*kk1;
      st[ijk] -= ( cg0*(u[ijk-ii1] * (ci0*s[ijk-ii3] + ci1*s[ijk-ii2] + ci2*s[ijk-ii1] + ci3*s[ijk    ]))
                 + cg1*(u[ijk    ] * (ci0*s[ijk-ii2] + ci1*s[ijk-ii1] + ci2*s[ijk    ] + ci3*s[ijk+ii1]))
                 + cg2*(u[ijk+ii1] * (ci0*s[ijk-ii1] + ci1*s[ijk    ] + ci2*s[ijk+ii1] + ci3*s[ijk+ii2]))
                 + cg3*(u[ijk+ii2] * (ci0*s[ijk    ] + ci1*s[ijk+ii1] + ci2*s[ijk+ii2] + ci3*s[ijk+ii3])) ) * cgi*dxi;

      st[ijk] -= ( cg0*(v[ijk-jj1] * (ci0*s[ijk-jj3] + ci1*s[ijk-jj2] + ci2*s[ijk-jj1] + ci3*s[ijk    ]))
                 + cg1*(v[ijk    ] * (ci0*s[ijk-jj2] + ci1*s[ijk-jj1] + ci2*s[ijk    ] + ci3*s[ijk+jj1]))
                 + cg2*(v[ijk+jj1] * (ci0*s[ijk-jj1] + ci1*s[ijk    ] + ci2*s[ijk+jj1] + ci3*s[ijk+jj2]))
                 + cg3*(v[ijk+jj2] * (ci0*s[ijk    ] + ci1*s[ijk+jj1] + ci2*s[ijk+jj2] + ci3*s[ijk+jj3])) ) * cgi*dyi;

      st[ijk] -= ( cg0*(w[ijk-kk1] * interp4       (s[ijk-kk3], s[ijk-kk2], s[ijk-kk1], s[ijk    ]))
                 + cg1*(w[ijk    ] * interp4       (s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]))
                 + cg2*(w[ijk+kk1] * interp4       (s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2]))
                 + cg3*(w[ijk+kk2] * interp4biastop(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])) )
                 * dzi4[kend-1];
    }

  return 0;
}

inline double cadvec_g4::interp2(const double a, const double b)
{
  return 0.5*(a + b);
}

inline double cadvec_g4::interp4(const double a, const double b, const double c, const double d)
{
  return (-a + 9.*b + 9.*c - d) / 16.;
}

inline double cadvec_g4::interp4biasbot(const double a, const double b, const double c, const double d)
{
  return ((5./16.)*a + (15./16.)*b - (5./16.)*c + (1./16)*d);
}

inline double cadvec_g4::interp4biastop(const double a, const double b, const double c, const double d)
{
  return ((5./16.)*d + (15./16.)*c - (5./16.)*b + (1./16)*a);
}

inline double cadvec_g4::grad4(const double a, const double b, const double c, const double d, const double dxi)
{
  return ( -(1./24.)*(d-a) + (27./24.)*(c-b) ) * dxi;
}

inline double cadvec_g4::grad4x(const double a, const double b, const double c, const double d)
{
  return (-(d-a) + 27.*(c-b)); 
}

inline double cadvec_g4::grad4xbiasbot(const double a, const double b, const double c, const double d)
{
  return (-23.*a + 21.*b + 3.*c - d);
}

inline double cadvec_g4::grad4xbiastop(const double a, const double b, const double c, const double d)
{
  return ( 23.*d - 21.*c - 3.*b + a);
}

