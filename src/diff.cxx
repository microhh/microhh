#include <cstdio>
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "diff.h"
#include "defines.h"

cdiff::cdiff(cgrid *gridin, cfields *fieldsin)
{
  std::printf("Creating instance of object diff\n");
  grid   = gridin;
  fields = fieldsin;
}

cdiff::~cdiff()
{
  std::printf("Destroying instance of object diff\n");
}

int cdiff::init()
{
  // get the maximum time step for diffusion

  return 0;
}

int cdiff::exec()
{
  // diffuse the flow
  diffc_2nd((*fields->ut).data, (*fields->u).data, grid->dzi, grid->dzhi, fields->visc);
  diffc_2nd((*fields->vt).data, (*fields->v).data, grid->dzi, grid->dzhi, fields->visc);
  diffw_2nd((*fields->wt).data, (*fields->w).data, grid->dzi, grid->dzhi, fields->visc);

  diffc_2nd((*fields->st).data, (*fields->s).data, grid->dzi, grid->dzhi, fields->viscs);
  return 0;
}

int cdiff::diffc_2nd(double * restrict at, double * restrict a, double * restrict dzi, double * restrict dzhi, double visc)
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
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        at[ijk] += visc * (
              + (  (a[ijk+ii] - a[ijk   ]) 
                 - (a[ijk   ] - a[ijk-ii]) ) * dxidxi 
              + (  (a[ijk+jj] - a[ijk   ]) 
                 - (a[ijk   ] - a[ijk-jj]) ) * dyidyi
              + (  (a[ijk+kk] - a[ijk   ]) * dzhi[k+1]
                 - (a[ijk   ] - a[ijk-kk]) * dzhi[k]   ) * dzi[k] );
      }

  return 0;
}

int cdiff::diffw_2nd(double * restrict wt, double * restrict w, double * restrict dzi, double * restrict dzhi, double visc)
{
  int    ijk,ii,jj,kk;
  double dxidxi,dyidyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dxidxi = 1./(grid->dx*grid->dx);
  dyidyi = 1./(grid->dy*grid->dy);

  for(int k=grid->kstart+1; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        wt[ijk] += visc * (
              + (  (w[ijk+ii] - w[ijk   ]) 
                 - (w[ijk   ] - w[ijk-ii]) ) * dxidxi 
              + (  (w[ijk+jj] - w[ijk   ]) 
                 - (w[ijk   ] - w[ijk-jj]) ) * dyidyi
              + (  (w[ijk+kk] - w[ijk   ]) * dzi[k]
                 - (w[ijk   ] - w[ijk-kk]) * dzi[k-1] ) * dzhi[k] );
      }

  return 0;
}

