#include <cstdio>
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "advec.h"

cadvec::cadvec(cgrid *gridin, cfields *fieldsin)
{
  std::printf("Creating instance of object advec\n");
  grid   = gridin;
  fields = fieldsin;
}

cadvec::~cadvec()
{
  std::printf("Destroying instance of object advec\n");
}

int cadvec::exec()
{
  // double c;
  // compute the courant number
  // c = courant((*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi, dt);

  // advect the flow
  advecu_2nd((*fields->ut).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi );
  advecv_2nd((*fields->vt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi );
  advecw_2nd((*fields->wt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzhi);

  return 0;
}

double cadvec::getcfl(double dt)
{
  double cfl;
  cfl = calccfl((*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzi, dt);

  return cfl;
}

double cadvec::calccfl(double * __restrict__ u, double * __restrict__ v, double * __restrict__ w, double * __restrict__ dzi, double dt)
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

  double cfl;

  cfl = 0;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*icells + k*ijcells;
        cfl = std::max(cfl, std::abs(interp2(u[ijk], u[ijk+ii]))*dxi + std::abs(interp2(v[ijk], v[ijk+jj]))*dyi + std::abs(interp2(w[ijk], w[ijk+kk]))*dzi[k]);
      }

  cfl = cfl*dt;

  return cfl;
}

int cadvec::advecu_2nd(double * __restrict__ ut, double * __restrict__ u, double * __restrict__ v, double * __restrict__ w, double * __restrict__ dzi)
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
        ut[ijk] += 
              - (  interp2(u[ijk   ], u[ijk+ii]) * interp2(u[ijk   ], u[ijk+ii])
                 - interp2(u[ijk-ii], u[ijk   ]) * interp2(u[ijk-ii], u[ijk   ]) ) * dxi

              - (  interp2(v[ijk-ii+jj], v[ijk+jj]) * interp2(u[ijk   ], u[ijk+jj])
                 - interp2(v[ijk-ii   ], v[ijk   ]) * interp2(u[ijk-jj], u[ijk   ]) ) * dyi 

              - (  interp2(w[ijk-ii+kk], w[ijk+kk]) * interp2(u[ijk   ], u[ijk+kk])
                 - interp2(w[ijk-ii   ], w[ijk   ]) * interp2(u[ijk-kk], u[ijk   ]) ) * dzi[k];
      }

  return 0;
}

int cadvec::advecv_2nd(double * __restrict__ vt, double * __restrict__ u, double * __restrict__ v, double * __restrict__ w, double * __restrict__ dzi)
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
        vt[ijk] += 
              - (  interp2(u[ijk+ii-jj], u[ijk+ii]) * interp2(v[ijk   ], v[ijk+ii])
                 - interp2(u[ijk   -jj], u[ijk   ]) * interp2(v[ijk-ii], v[ijk   ]) ) * dxi

              - (  interp2(v[ijk   ], v[ijk+jj]) * interp2(v[ijk   ], v[ijk+jj])
                 - interp2(v[ijk-jj], v[ijk   ]) * interp2(v[ijk-jj], v[ijk   ]) ) * dyi

              - (  interp2(w[ijk-jj+kk], w[ijk+kk]) * interp2(v[ijk   ], v[ijk+kk])
                 - interp2(w[ijk-jj   ], w[ijk   ]) * interp2(v[ijk-kk], v[ijk   ]) ) * dzi[k];
      }

  return 0;
}

int cadvec::advecw_2nd(double * __restrict__ wt, double * __restrict__ u, double * __restrict__ v, double * __restrict__ w, double * __restrict__ dzhi)
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

  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*icells + k*ijcells;
        wt[ijk] += 
              - (  interp2(u[ijk+ii-kk], u[ijk+ii]) * interp2(w[ijk   ], w[ijk+ii])
                 - interp2(u[ijk   -kk], u[ijk   ]) * interp2(w[ijk-ii], w[ijk   ]) ) * dxi

              - (  interp2(v[ijk   -kk], v[ijk+jj]) * interp2(w[ijk   ], w[ijk+jj])
                 - interp2(v[ijk-jj-kk], v[ijk   ]) * interp2(w[ijk-jj], w[ijk   ]) ) * dyi

              - (  interp2(w[ijk   ], w[ijk+kk]) * interp2(w[ijk   ], w[ijk+kk])
                 - interp2(w[ijk-kk], w[ijk   ]) * interp2(w[ijk-kk], w[ijk   ]) ) * dzhi[k];
      }

  return 0;
}

inline double cadvec::interp2(const double a, const double b)
{
  return 0.5*(a + b);
}
