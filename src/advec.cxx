#include <cstdio>
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
  // advect the flow
  advecu_2nd((*fields->ut).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dz );
  advecv_2nd((*fields->vt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dz );
  advecw_2nd((*fields->wt).data, (*fields->u).data, (*fields->v).data, (*fields->w).data, grid->dzh);
  return 0;
}

// high performance routine, restrict specifies that arrays are not aliased
int cadvec::advecu_2nd(double * __restrict__ ut, double * __restrict__ u, double * __restrict__ v, double * __restrict__ w, double * __restrict__ dz)
{
  int    ijk,icells,ijcells,ii,jj,kk;
  double dx, dy;

  icells  = grid->icells;
  ijcells = grid->icells*grid->jcells;

  ii = 1;
  jj = 1*icells;
  kk = 1*ijcells;

  dx = grid->dx;
  dy = grid->dy;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*icells + k*ijcells;
        ut[ijk] += 
              - (  interp2(u[ijk   ], u[ijk+ii]) * interp2(u[ijk   ], u[ijk+ii])
                 - interp2(u[ijk-ii], u[ijk   ]) * interp2(u[ijk-ii], u[ijk   ]) ) / dx

              - (  interp2(v[ijk-ii+jj], v[ijk+jj]) * interp2(u[ijk   ], u[ijk+jj])
                 - interp2(v[ijk-ii   ], v[ijk   ]) * interp2(u[ijk-jj], u[ijk   ]) ) / dy

              - (  interp2(w[ijk-ii+kk], w[ijk+kk]) * interp2(u[ijk   ], u[ijk+kk])
                 - interp2(w[ijk-ii   ], w[ijk   ]) * interp2(u[ijk-kk], u[ijk   ]) ) / dz[k];
      }

  return 0;
}

int cadvec::advecv_2nd(double * __restrict__ vt, double * __restrict__ u, double * __restrict__ v, double * __restrict__ w, double * __restrict__ dz)
{
  int    ijk,icells,ijcells,ii,jj,kk;
  double dx, dy;

  icells  = grid->icells;
  ijcells = grid->icells*grid->jcells;

  ii = 1;
  jj = 1*icells;
  kk = 1*ijcells;

  dx = grid->dx;
  dy = grid->dy;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*icells + k*ijcells;
        vt[ijk] += 
              - (  interp2(u[ijk+ii-jj], u[ijk+ii]) * interp2(v[ijk   ], v[ijk+ii])
                 - interp2(u[ijk   -jj], u[ijk   ]) * interp2(v[ijk-ii], v[ijk   ]) ) / dx

              - (  interp2(v[ijk   ], v[ijk+jj]) * interp2(v[ijk   ], v[ijk+jj])
                 - interp2(v[ijk-jj], v[ijk   ]) * interp2(v[ijk-jj], v[ijk   ]) ) / dy

              - (  interp2(w[ijk-jj+kk], w[ijk+kk]) * interp2(v[ijk   ], v[ijk+kk])
                 - interp2(w[ijk-jj   ], w[ijk   ]) * interp2(v[ijk-kk], v[ijk   ]) ) / dz[k];
      }

  return 0;
}

int cadvec::advecw_2nd(double * __restrict__ wt, double * __restrict__ u, double * __restrict__ v, double * __restrict__ w, double * __restrict__ dzh)
{
  int    ijk,icells,ijcells,ii,jj,kk;
  double dx, dy;

  icells  = grid->icells;
  ijcells = grid->icells*grid->jcells;

  ii = 1;
  jj = 1*icells;
  kk = 1*ijcells;

  dx  = grid->dx;
  dy  = grid->dy;

  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*icells + k*ijcells;
        wt[ijk] += 
              - (  interp2(u[ijk+ii-kk], u[ijk+ii]) * interp2(w[ijk   ], w[ijk+ii])
                 - interp2(u[ijk   -kk], u[ijk   ]) * interp2(w[ijk-ii], w[ijk   ]) ) / dx

              - (  interp2(v[ijk   -kk], v[ijk+jj]) * interp2(w[ijk   ], w[ijk+jj])
                 - interp2(v[ijk-jj-kk], v[ijk   ]) * interp2(w[ijk-jj], w[ijk   ]) ) / dy

              - (  interp2(w[ijk   ], w[ijk+kk]) * interp2(w[ijk   ], w[ijk+kk])
                 - interp2(w[ijk-kk], w[ijk   ]) * interp2(w[ijk-kk], w[ijk   ]) ) / dzh[k];
      }

  return 0;
}

inline double cadvec::interp2(const double a, const double b)
{
  return 0.5*(a + b);
}
