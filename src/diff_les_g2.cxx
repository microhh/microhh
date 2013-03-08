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

int cdiff_les_g2::evisc(double * restrict evisc, double * restrict u, double * restrict v, double * restrict w, double * restrict s, double * restrict z, double * restrict dz, double * restrict dzi, double * restrict dzhi, double tPr)
{
  int    ijk,ii,jj,kk;
  double dx,dy,dxi,dyi;

  // wall damping
  double mlen,mlen0,fac;
  const double cs = 0.23;
  const double z0 = 0.1;
  const double n  = 2.;

  double strain2, RitPrratio;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dx = grid->dx;
  dy = grid->dy;
  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart; k<grid->kend; k++)
  {
    // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
    mlen0 = cs*std::pow(dx*dy*dz[k], 1./3.);
    mlen  = std::pow(1./(1./std::pow(mlen0, n) + 1./(std::pow(kappa*(z[k]+z0), n))), 1./n);
    fac   = std::pow(mlen, 2.);

    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        strain2 = 2.*(
          // du/dx + du/dx
          + std::pow((u[ijk+ii]-u[ijk])*dxi, 2.)

          // dv/dy + dv/dy
          + std::pow((v[ijk+jj]-v[ijk])*dyi, 2.)

          // dw/dz + dw/dz
          + std::pow((w[ijk+kk]-w[ijk])*dzi[k], 2.)

          // du/dy + dv/dx
          + 0.125*std::pow((u[ijk      ]-u[ijk   -jj])*dyi  + (v[ijk      ]-v[ijk-ii   ])*dxi, 2.)
          + 0.125*std::pow((u[ijk+ii   ]-u[ijk+ii-jj])*dyi  + (v[ijk+ii   ]-v[ijk      ])*dxi, 2.)
          + 0.125*std::pow((u[ijk   +jj]-u[ijk      ])*dyi  + (v[ijk   +jj]-v[ijk-ii+jj])*dxi, 2.)
          + 0.125*std::pow((u[ijk+ii+jj]-u[ijk+ii   ])*dyi  + (v[ijk+ii+jj]-v[ijk   +jj])*dxi, 2.)

          // du/dz + dw/dx
          + 0.125*std::pow((u[ijk      ]-u[ijk   -kk])*dzhi[k  ] + (w[ijk      ]-w[ijk-ii   ])*dxi, 2.)
          + 0.125*std::pow((u[ijk+ii   ]-u[ijk+ii-kk])*dzhi[k  ] + (w[ijk+ii   ]-w[ijk      ])*dxi, 2.)
          + 0.125*std::pow((u[ijk   +kk]-u[ijk      ])*dzhi[k+1] + (w[ijk   +kk]-w[ijk-ii+kk])*dxi, 2.)
          + 0.125*std::pow((u[ijk+ii+kk]-u[ijk+ii   ])*dzhi[k+1] + (w[ijk+ii+kk]-w[ijk   +kk])*dxi, 2.)

          // dv/dz + dw/dy
          + 0.125*std::pow((v[ijk      ]-v[ijk   -kk])*dzhi[k  ] + (w[ijk      ]-w[ijk-jj   ])*dyi, 2.)
          + 0.125*std::pow((v[ijk+jj   ]-v[ijk+jj-kk])*dzhi[k  ] + (w[ijk+jj   ]-w[ijk      ])*dyi, 2.)
          + 0.125*std::pow((v[ijk   +kk]-v[ijk      ])*dzhi[k+1] + (w[ijk   +kk]-w[ijk-jj+kk])*dyi, 2.)
          + 0.125*std::pow((v[ijk+jj+kk]-v[ijk+jj   ])*dzhi[k+1] + (w[ijk+jj+kk]-w[ijk   +kk])*dyi, 2.) );

        // CvH use the thermal expansion coefficient from the input later, what to do if there is no buoyancy?
        // Add the buoyancy production to the TKE
        RitPrratio = (9.81/300.)*(s[ijk+kk]-s[ijk-kk])*0.5*dzi[k] / strain2 / tPr;
        RitPrratio = std::max(0., std::min(RitPrratio, 1.-dsmall));
        evisc[ijk] = fac * std::sqrt(strain2) * std::sqrt(1.-RitPrratio);
      }
  }

  grid->boundary_cyclic(evisc);

  return 0;
}

int cdiff_les_g2::diffu(double * restrict ut, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi, double * restrict dzhi, double * restrict evisc, double * restrict fluxbot, double * restrict fluxtop)
{
  int    ijk,ij,ii,jj,kk,kstart,kend;
  double dxi,dyi;
  double eviscn, eviscs, eviscb, evisct;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  kstart = grid->kstart;
  kend   = grid->kend;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      eviscn = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
      eviscs = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
      evisct = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
      eviscb = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);
      ut[ijk] +=
            // du/dx + du/dx
            + (  evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
               - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
            // du/dy + dv/dx
            + (  eviscn*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
               - eviscs*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
            // du/dz + dw/dx
            + (  evisct*((u[ijk+kk]-u[ijk   ])* dzhi[kstart+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
               + fluxbot[ij] ) * dzi[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        eviscn = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
        eviscs = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
        evisct = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
        eviscb = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);
        ut[ijk] +=
              // du/dx + du/dx
              + (  evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
                 - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
              // du/dy + dv/dx
              + (  eviscn*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                 - eviscs*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
              // du/dz + dw/dx
              + (  evisct*((u[ijk+kk]-u[ijk   ])* dzhi[k+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
                 - eviscb*((u[ijk   ]-u[ijk-kk])* dzhi[k  ] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) * dzi[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + (kend-1)*kk;
      eviscn = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
      eviscs = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
      evisct = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
      eviscb = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);
      ut[ijk] +=
            // du/dx + du/dx
            + (  evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
               - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
            // du/dy + dv/dx
            + (  eviscn*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
               - eviscs*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
            // du/dz + dw/dx
            + (- fluxtop[ij]
               - eviscb*((u[ijk   ]-u[ijk-kk])* dzhi[kend-1] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) * dzi[kend-1];
    }

  return 0;
}

int cdiff_les_g2::diffv(double * restrict vt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi, double * restrict dzhi, double * restrict evisc, double * restrict fluxbot, double * restrict fluxtop)
{
  int    ijk,ij,ii,jj,kk,kstart,kend;
  double dxi,dyi;
  double evisce,eviscw,eviscb,evisct;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  kstart = grid->kstart;
  kend   = grid->kend;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      evisce = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
      eviscw = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
      evisct = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
      eviscb = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);
      vt[ijk] +=
            // dv/dx + du/dy
            + (  evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
               - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
            // dv/dy + dv/dy
            + (  evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
               - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi;
            // dv/dz + dw/dy
            + (  evisct*((v[ijk+kk]-v[ijk   ])*dzhi[kstart+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
               + fluxbot[ij] ) * dzi[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        evisce = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
        eviscw = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
        evisct = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
        eviscb = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);
        vt[ijk] +=
              // dv/dx + du/dy
              + (  evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                 - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
              // dv/dy + dv/dy
              + (  evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
                 - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi;
              // dv/dz + dw/dy
              + (  evisct*((v[ijk+kk]-v[ijk   ])*dzhi[k+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
                 - eviscb*((v[ijk   ]-v[ijk-kk])*dzhi[k  ] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) * dzi[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + (kend-1)*kk;
      evisce = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
      eviscw = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
      evisct = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
      eviscb = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);
      vt[ijk] +=
            // dv/dx + du/dy
            + (  evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
               - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
            // dv/dy + dv/dy
            + (  evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
               - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi;
            // dv/dz + dw/dy
            + (- fluxtop[ij]
               - eviscb*((v[ijk   ]-v[ijk-kk])*dzhi[kend-1] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) * dzi[kend-1];
    }

  return 0;
}

int cdiff_les_g2::diffw(double * restrict wt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi, double * restrict dzhi, double * restrict evisc)
{
  int    ijk,ii,jj,kk;
  double dxi,dyi;
  double evisce, eviscw, eviscn, eviscs;

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
        evisce = 0.25*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+ii-kk] + evisc[ijk+ii]);
        eviscw = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-ii] + evisc[ijk   -kk] + evisc[ijk   ]);
        eviscn = 0.25*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+jj-kk] + evisc[ijk+jj]);
        eviscs = 0.25*(evisc[ijk-jj-kk] + evisc[ijk-jj] + evisc[ijk   -kk] + evisc[ijk   ]);
        wt[ijk] +=
              // dw/dx + du/dz
              + (  evisce*((w[ijk+ii]-w[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-kk])*dzhi[k])
                 - eviscw*((w[ijk   ]-w[ijk-ii])*dxi + (u[ijk   ]-u[ijk+  -kk])*dzhi[k]) ) * dxi
              // dw/dy + dv/dz
              + (  eviscn*((w[ijk+jj]-w[ijk   ])*dyi + (v[ijk+jj]-v[ijk+jj-kk])*dzhi[k])
                 - eviscs*((w[ijk   ]-w[ijk-jj])*dyi + (v[ijk   ]-v[ijk+  -kk])*dzhi[k]) ) * dyi
              // dw/dz + dw/dz
              + (  evisc[ijk   ]*(w[ijk+kk]-w[ijk   ])*dzi[k  ]
                 - evisc[ijk-kk]*(w[ijk   ]-w[ijk-kk])*dzi[k-1] ) * 2.* dzhi[k];
      }

  return 0;
}

int cdiff_les_g2::diffc(double * restrict at, double * restrict a, double * restrict dzi, double * restrict dzhi, double * restrict evisc, double * restrict fluxbot, double * restrict fluxtop, double tPr)
{
  int    ijk,ij,ii,jj,kk,kstart,kend;
  double dxidxi,dyidyi;
  double evisce,eviscw,eviscn,eviscs,evisct,eviscb;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  kstart = grid->kstart;
  kend   = grid->kend;

  dxidxi = 1./(grid->dx * grid->dx);
  dyidyi = 1./(grid->dy * grid->dy);

  // bottom boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])/tPr;
      eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])/tPr;
      eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])/tPr;
      eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])/tPr;
      evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])/tPr;
      eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])/tPr;

      at[ijk] +=
            + (  evisce*(a[ijk+ii]-a[ijk   ]) 
               - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi 
            + (  eviscn*(a[ijk+jj]-a[ijk   ]) 
               - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
            + (  evisct*(a[ijk+kk]-a[ijk   ])*dzhi[kstart+1]
               + fluxbot[ij] ) * dzi[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])/tPr;
        eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])/tPr;
        eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])/tPr;
        eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])/tPr;
        evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])/tPr;
        eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])/tPr;

        at[ijk] +=
              + (  evisce*(a[ijk+ii]-a[ijk   ]) 
                 - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi 
              + (  eviscn*(a[ijk+jj]-a[ijk   ]) 
                 - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
              + (  evisct*(a[ijk+kk]-a[ijk   ])*dzhi[k+1]
                 - eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[k]  ) * dzi[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; i++)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + (kend-1)*kk;
      evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])/tPr;
      eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])/tPr;
      eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])/tPr;
      eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])/tPr;
      evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])/tPr;
      eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])/tPr;

      at[ijk] +=
            + (  evisce*(a[ijk+ii]-a[ijk   ]) 
               - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi 
            + (  eviscn*(a[ijk+jj]-a[ijk   ]) 
               - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
            + (- fluxtop[ij]
               - eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[kend-1]  ) * dzi[kend-1];
    }

  return 0;
}

double cdiff_les_g2::getdn(double * restrict evisc, double * restrict dzi, double tPr)
{
  int    ijk,ij,ii,jj,kk;
  double dxidxi,dyidyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dxidxi = 1./(grid->dx * grid->dx);
  dyidyi = 1./(grid->dy * grid->dy);

  double tPrfac = std::min(1., tPr);
  double dnmul = 0;

  // get the maximum time step for diffusion
  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk = i + j*jj + k*kk;
        dnmul = std::max(dnmul, std::abs(tPrfac*evisc[ijk]*(dxidxi + dyidyi + dzi[k]*dzi[k])));
      }

  return dnmul;
}
