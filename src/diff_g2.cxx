#include <cstdio>
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "diff_g2.h"
#include "defines.h"

cdiff_g2::cdiff_g2(cgrid *gridin, cfields *fieldsin, cmpi *mpiin) : cdiff(gridin, fieldsin, mpiin)
{
}

cdiff_g2::~cdiff_g2()
{
}

int cdiff_g2::setvalues()
{
  // get the maximum time step for diffusion
  double viscmax = fields->visc;
  for(fieldmap::iterator it = fields->sp.begin(); it!=fields->sp.end(); it++)
    viscmax = std::max(it->second->visc, viscmax);

  dnmul = 0;
  for(int k=grid->kstart; k<grid->kend; k++)
    dnmul = std::max(dnmul, std::abs(viscmax * (1./(grid->dx*grid->dx) + 1./(grid->dy*grid->dy) + 1./(grid->dz[k]*grid->dz[k]))));

  return 0;
}

double cdiff_g2::getdn(double dt)
{
  double dn;

  dn = dnmul*dt;

  return dn;
}

int cdiff_g2::exec()
{
  diffc(fields->ut->data, fields->u->data, grid->dzi, grid->dzhi, fields->visc);
  diffc(fields->vt->data, fields->v->data, grid->dzi, grid->dzhi, fields->visc);
  diffw(fields->wt->data, fields->w->data, grid->dzi, grid->dzhi, fields->visc);

  for(fieldmap::iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    diffc((*it->second).data, (*fields->s[it->first]).data, grid->dzi, grid->dzhi, fields->s[it->first]->visc);
  
  return 0;
}

int cdiff_g2::diffc(double * restrict at, double * restrict a, double * restrict dzi, double * restrict dzhi, double visc)
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
              + (  (a[ijk+ii] - a[ijk   ]) 
                 - (a[ijk   ] - a[ijk-ii]) ) * dxidxi 
              + (  (a[ijk+jj] - a[ijk   ]) 
                 - (a[ijk   ] - a[ijk-jj]) ) * dyidyi
              + (  (a[ijk+kk] - a[ijk   ]) * dzhi[k+1]
                 - (a[ijk   ] - a[ijk-kk]) * dzhi[k]   ) * dzi[k] );
      }

  return 0;
}

int cdiff_g2::diffw(double * restrict wt, double * restrict w, double * restrict dzi, double * restrict dzhi, double visc)
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
#pragma ivdep
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

