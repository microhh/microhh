#include <cstdio>
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

  pres_2nd_in((*fields->p ).data,
              (*fields->u ).data, (*fields->v ).data, (*fields->w ).data,
              (*fields->ut).data, (*fields->vt).data, (*fields->wt).data, 
              grid->dzi, dt);

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

