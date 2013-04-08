#include <cstdio>
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "advec_g42.h"
#include "defines.h"

cadvec_g42::cadvec_g42(cgrid *gridin, cfields *fieldsin, cmpi *mpiin) : cadvec(gridin, fieldsin, mpiin)
{
  grid   = gridin;
  fields = fieldsin;
  mpi    = mpiin;
}

cadvec_g42::~cadvec_g42()
{
}

unsigned long cadvec_g42::gettimelim(unsigned long idt, double dt)
{
  unsigned long idtlim;

  idtlim = idt * cflmax / calccfl(fields->u->data, fields->v->data, fields->w->data, grid->dzi, dt);

  return idtlim;
}


double cadvec_g42::getcfl(double dt)
{
  double cfl;

  cfl = calccfl((*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi, dt);

  return cfl;
}

int cadvec_g42::exec()
{
  advecu((*fields->ut).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi );
  advecv((*fields->vt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi );
  advecw((*fields->wt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzhi);

  for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    advecs((*it->second).data, (*fields->s[it->first]).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi);

  return 0;
}

double cadvec_g42::calccfl(double * restrict u, double * restrict v, double * restrict w, double * restrict dzi, double dt)
{
  int    ijk,kk;
  int    ii1,ii2,jj1,jj2;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  kk  = grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;


  double cfl = 0;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj1 + k*kk;
        cfl = std::max(cfl, std::abs(interp4(u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2]))*dxi 
                          + std::abs(interp4(v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2]))*dyi 
                          + std::abs(interp2(w[ijk], w[ijk+kk]))*dzi[k]);
      }

  grid->getmax(&cfl);

  cfl = cfl*dt;

  return cfl;
}

int cadvec_g42::advecu(double * restrict ut, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi)
{
  int    ijk,kk;
  int    ii1,ii2,ii3,jj1,jj2,jj3;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk  = grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk;
        ut[ijk] += 
              // Morinishi
              // - grad4(interp4(u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ]) * interp2(u[ijk-ii3], u[ijk    ]),
              //         interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) * interp2(u[ijk-ii1], u[ijk    ]),
              //         interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]) * interp2(u[ijk    ], u[ijk+ii1]),
              //         interp4(u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]) * interp2(u[ijk    ], u[ijk+ii3]), dxi)

              - grad4(interp4(u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ]) * interp4(u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ]),
                      interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]),
                      interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]) * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]),
                      interp4(u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]) * interp4(u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]), dxi)

              - grad4(interp4(v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1]) * interp4(u[ijk-jj3], u[ijk-jj2], u[ijk-jj1], u[ijk    ]),
                      interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]),
                      interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2]),
                      interp4(v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2]) * interp4(u[ijk    ], u[ijk+jj1], u[ijk+jj2], u[ijk+jj3]), dyi)

              - (  interp4(w[ijk-ii2+kk], w[ijk-ii1+kk], w[ijk+kk], w[ijk+ii1+kk]) * interp2(u[ijk    ], u[ijk+kk ])
                 - interp4(w[ijk-ii2   ], w[ijk-ii1   ], w[ijk   ], w[ijk+ii1   ]) * interp2(u[ijk-kk ], u[ijk    ]) ) * dzi[k];
      }

  return 0;
}

int cadvec_g42::advecv(double * restrict vt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi)
{
  int    ijk,kk;
  int    ii1,ii2,ii3,jj1,jj2,jj3;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk  = grid->icells*grid->jcells;
  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk;
        vt[ijk] += 
              - grad4(interp4(u[ijk-ii1-jj2], u[ijk-ii1-jj1], u[ijk-ii1], u[ijk-ii1+jj1]) * interp4(v[ijk-ii3], v[ijk-ii2], v[ijk-ii1], v[ijk    ]),
                      interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]),
                      interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2]),
                      interp4(u[ijk+ii2-jj2], u[ijk+ii2-jj1], u[ijk+ii2], u[ijk+ii2+jj1]) * interp4(v[ijk    ], v[ijk+ii1], v[ijk+ii2], v[ijk+ii3]), dxi)

              - grad4(interp4(v[ijk-jj3], v[ijk-jj2], v[ijk-jj1], v[ijk    ]) * interp4(v[ijk-jj3], v[ijk-jj2], v[ijk-jj1], v[ijk    ]),
                      interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]),
                      interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]),
                      interp4(v[ijk    ], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]) * interp4(v[ijk    ], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]), dyi)

              - (  interp4(w[ijk-jj2+kk], w[ijk-jj1+kk], w[ijk+kk], w[ijk+jj1+kk]) * interp2(v[ijk   ], v[ijk+kk])
                 - interp4(w[ijk-jj2   ], w[ijk-jj1   ], w[ijk   ], w[ijk+jj1   ]) * interp2(v[ijk-kk], v[ijk   ]) ) * dzi[k];
      }

  return 0;
}

int cadvec_g42::advecw(double * restrict wt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzhi)
{
  int    ijk;
  int    ii1,ii2,ii3,jj1,jj2,jj3,kk1,kk2,kk3;
  double dxi,dyi;

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

  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk1;
        wt[ijk] +=
              // Morinishi
              // - grad4(interp4(u[ijk-ii1-kk2], u[ijk-ii1-kk1], u[ijk-ii1], u[ijk-ii1+kk1]) * interp2(w[ijk-ii3], w[ijk    ]),
              //         interp4(u[ijk    -kk2], u[ijk    -kk1], u[ijk    ], u[ijk    +kk1]) * interp2(w[ijk-ii1], w[ijk    ]),
              //         interp4(u[ijk+ii1-kk2], u[ijk+ii1-kk1], u[ijk+ii1], u[ijk+ii1+kk1]) * interp2(w[ijk    ], w[ijk+ii1]),
              //         interp4(u[ijk+ii2-kk2], u[ijk+ii2-kk1], u[ijk+ii2], u[ijk+ii2+kk1]) * interp2(w[ijk    ], w[ijk+ii3]), dxi)

              - grad4(interp2(u[ijk-ii1-kk1], u[ijk-ii1]) * interp4(w[ijk-ii3], w[ijk-ii2], w[ijk-ii1], w[ijk    ]),
                      interp2(u[ijk    -kk1], u[ijk    ]) * interp4(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]),
                      interp2(u[ijk+ii1-kk1], u[ijk+ii1]) * interp4(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2]),
                      interp2(u[ijk+ii2-kk1], u[ijk+ii2]) * interp4(w[ijk    ], w[ijk+ii1], w[ijk+ii2], w[ijk+ii3]), dxi)

              - grad4(interp2(v[ijk-jj1-kk1], v[ijk-jj1]) * interp4(w[ijk-jj3], w[ijk-jj2], w[ijk-jj1], w[ijk    ]),
                      interp2(v[ijk    -kk1], v[ijk    ]) * interp4(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]),
                      interp2(v[ijk+jj1-kk1], v[ijk+jj1]) * interp4(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2]),
                      interp2(v[ijk+jj2-kk1], v[ijk+jj2]) * interp4(w[ijk    ], w[ijk+jj1], w[ijk+jj2], w[ijk+jj3]), dyi)

              - (  interp2(w[ijk    ], w[ijk+kk1]) * interp2(w[ijk    ], w[ijk+kk1])
                 - interp2(w[ijk-kk1], w[ijk    ]) * interp2(w[ijk-kk1], w[ijk    ]) ) * dzhi[k];
      }

  return 0;
}

int cadvec_g42::advecs(double * restrict st, double * restrict s, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi)
{
  int    ijk,kk;
  int    ii1,ii2,ii3,jj1,jj2,jj3;
  double dxi,dyi;

  ii1 = 1;
  ii2 = 2;
  ii3 = 3;
  jj1 = 1*grid->icells;
  jj2 = 2*grid->icells;
  jj3 = 3*grid->icells;
  kk  = grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj1 + k*kk;
        st[ijk] += 
              // Morinishi
              // - grad4(u[ijk-ii1] * interp2(s[ijk-ii3], s[ijk    ]),
              //         u[ijk    ] * interp2(s[ijk-ii1], s[ijk    ]),
              //         u[ijk+ii1] * interp2(s[ijk    ], s[ijk+ii1]),
              //         u[ijk+ii2] * interp2(s[ijk    ], s[ijk+ii3]), dxi)

              // - grad4(v[ijk-jj1] * interp2(s[ijk-jj3], s[ijk    ]),
              //         v[ijk    ] * interp2(s[ijk-jj1], s[ijk    ]),
              //         v[ijk+jj1] * interp2(s[ijk    ], s[ijk+jj1]),
              //         v[ijk+jj2] * interp2(s[ijk    ], s[ijk+jj3]), dyi)

              - grad4(u[ijk-ii1] * interp4(s[ijk-ii3], s[ijk-ii2], s[ijk-ii1], s[ijk    ]),
                      u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]),
                      u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2]),
                      u[ijk+ii2] * interp4(s[ijk    ], s[ijk+ii1], s[ijk+ii2], s[ijk+ii3]), dxi)

              - grad4(v[ijk-jj1] * interp4(s[ijk-jj3], s[ijk-jj2], s[ijk-jj1], s[ijk    ]),
                      v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]),
                      v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2]),
                      v[ijk+jj2] * interp4(s[ijk    ], s[ijk+jj1], s[ijk+jj2], s[ijk+jj3]), dyi)

              - (  w[ijk+kk] * interp2(s[ijk   ], s[ijk+kk])
                 - w[ijk   ] * interp2(s[ijk-kk], s[ijk   ]) ) * dzi[k];
      }

  return 0;
}

inline double cadvec_g42::interp2(const double a, const double b)
{
  return 0.5*(a + b);
}

inline double cadvec_g42::interp4(const double a, const double b, const double c, const double d)
{
  return (-a + 9.*b + 9.*c - d) / 16.;
}

inline double cadvec_g42::grad4(const double a, const double b, const double c, const double d, const double dxi)
{
  return ( -(1./24.)*(d-a) + (27./24.)*(c-b) ) * dxi;
}
