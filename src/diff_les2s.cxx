/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
 *
 * This file is part of MicroHH
 *
 * MicroHH is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * MicroHH is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with MicroHH.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdio>
#include <cmath>
#include <algorithm>
#include "grid.h"
#include "fields.h"
#include "master.h"
#include "diff_les2s.h"
#include "boundary_surface.h"
#include "defines.h"
#include "constants.h"
#include "thermo.h"
#include "model.h"

cdiff_les2s::cdiff_les2s(cmodel *modelin, cinput *inputin) : cdiff(modelin, inputin)
{
  swdiff = "les2s";

  int nerror = 0;
  nerror += inputin->getItem(&dnmax, "diff", "dnmax", "", 0.5  );
  nerror += inputin->getItem(&cs   , "diff", "cs"   , "", 0.23 );
  nerror += inputin->getItem(&tPr  , "diff", "tPr"  , "", 1./3.);

  nerror += fields->initdfld("evisc", "Eddy viscosity", "m2 s-1");

  if(nerror)
    throw 1;
}

cdiff_les2s::~cdiff_les2s()
{
}

unsigned long cdiff_les2s::gettimelim(unsigned long idt, double dt)
{
  unsigned long idtlim;
  double dnmul;

  dnmul = calcdnmul(fields->s["evisc"]->data, grid->dzi, this->tPr);
  // avoid zero division
  dnmul = std::max(constants::dsmall, dnmul);
  idtlim = idt * dnmax/(dnmul*dt);

  return idtlim;
}

//#ifndef USECUDA
int cdiff_les2s::execvisc()
{
  // do a cast because the base boundary class does not have the MOST related variables
  cboundary_surface *boundaryptr = static_cast<cboundary_surface *>(model->boundary);

  strain2(fields->s["evisc"]->data,
          fields->u->data, fields->v->data, fields->w->data,
          fields->u->datafluxbot, fields->v->datafluxbot,
          boundaryptr->ustar, boundaryptr->obuk,
          grid->z, grid->dzi, grid->dzhi);

  // start with retrieving the stability information
  if(model->thermo->getsw() == "0")
  {
    evisc_neutral(fields->s["evisc"]->data,
                  fields->u->data, fields->v->data, fields->w->data,
                  fields->u->datafluxbot, fields->v->datafluxbot,
                  grid->z, grid->dz, boundaryptr->z0m);
  }
  // assume buoyancy calculation is needed
  else
  {
    // store the buoyancyflux in tmp1
    model->thermo->getbuoyancyfluxbot(fields->sd["tmp1"]);
    // retrieve the full field in tmp1 and use tmp2 for temporary calculations
    model->thermo->getthermofield(fields->sd["tmp1"], fields->sd["tmp2"], "N2");
    // model->thermo->getthermofield(fields->sd["tmp1"], fields->sd["tmp2"], "b");

    evisc(fields->s["evisc"]->data,
          fields->u->data, fields->v->data, fields->w->data, fields->s["tmp1"]->data,
          fields->u->datafluxbot, fields->v->datafluxbot, fields->sd["tmp1"]->datafluxbot,
          boundaryptr->ustar, boundaryptr->obuk,
          grid->z, grid->dz, grid->dzi,
          boundaryptr->z0m);
  }

  return 0;
}
//#endif

double cdiff_les2s::getdn(double dt)
{
  double dnmul;

  // calculate eddy viscosity
  dnmul = calcdnmul(fields->s["evisc"]->data, grid->dzi, this->tPr);

  return dnmul*dt;
}

#ifndef USECUDA
int cdiff_les2s::exec()
{
  diffu(fields->ut->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi, grid->dzhi, fields->s["evisc"]->data, fields->u->datafluxbot, fields->u->datafluxtop, fields->rhoref, fields->rhorefh);
  diffv(fields->vt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi, grid->dzhi, fields->s["evisc"]->data, fields->v->datafluxbot, fields->v->datafluxtop, fields->rhoref, fields->rhorefh);
  diffw(fields->wt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi, grid->dzhi, fields->s["evisc"]->data, fields->rhoref, fields->rhorefh);

  for(fieldmap::const_iterator it = fields->st.begin(); it!=fields->st.end(); ++it)
    diffc(it->second->data, fields->s[it->first]->data, grid->dzi, grid->dzhi, fields->s["evisc"]->data, fields->s[it->first]->datafluxbot, fields->s[it->first]->datafluxtop, fields->rhoref, fields->rhorefh, this->tPr);

  return 0;
}
#endif

int cdiff_les2s::strain2(double * restrict strain2,
                          double * restrict u, double * restrict v, double * restrict w,
                          double * restrict ufluxbot, double * restrict vfluxbot,
                          double * restrict ustar, double * restrict obuk,
                          double * restrict z, double * restrict dzi, double * restrict dzhi)
{
  int    ij,ijk,ii,jj,kk,kstart;
  double dxi,dyi;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  kstart = grid->kstart;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      strain2[ijk] = 2.*(
        // du/dz
        + 0.5*std::pow(-0.5*(ufluxbot[ij]+ufluxbot[ij+ii])/(constants::kappa*z[kstart]*ustar[ij])*phim(z[kstart]/obuk[ij]), 2.)

        // dv/dz
        + 0.5*std::pow(-0.5*(vfluxbot[ij]+vfluxbot[ij+jj])/(constants::kappa*z[kstart]*ustar[ij])*phim(z[kstart]/obuk[ij]), 2.) );
      // add a small number to avoid zero divisions
      strain2[ijk] += constants::dsmall;
    }

  for(int k=grid->kstart+1; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        strain2[ijk] = 2.*(
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

        // add a small number to avoid zero divisions
        strain2[ijk] += constants::dsmall;
      }

  return 0;
}

int cdiff_les2s::evisc(double * restrict evisc,
                        double * restrict u, double * restrict v, double * restrict w,  double * restrict N2,
                        double * restrict ufluxbot, double * restrict vfluxbot, double * restrict bfluxbot,
                        double * restrict ustar, double * restrict obuk,
                        double * restrict z, double * restrict dz, double * restrict dzi,
                        double z0m)
{
  int    ij,ijk,jj,kk,kstart;
  double dx,dy;

  // wall damping
  double mlen,mlen0,fac;
  const double n = 2.;

  double RitPrratio;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  kstart = grid->kstart;

  dx = grid->dx;
  dy = grid->dy;

  // bottom boundary, here strain is fully parametrized using MO
  // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
  mlen0 = this->cs*std::pow(dx*dy*dz[kstart], 1./3.);
  mlen  = std::pow(1./(1./std::pow(mlen0, n) + 1./(std::pow(constants::kappa*(z[kstart]+z0m), n))), 1./n);
  fac   = std::pow(mlen, 2.);

  // local copies to aid vectorization
  double tPr = this->tPr;
  double cs  = this->cs;

  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
    {
      ij  = i + j*jj;
      ijk = i + j*jj + kstart*kk;
      // TODO use the thermal expansion coefficient from the input later, what to do if there is no buoyancy?
      // Add the buoyancy production to the TKE
      RitPrratio = -bfluxbot[ij]/(constants::kappa*z[kstart]*ustar[ij])*phih(z[kstart]/obuk[ij]) / evisc[ijk] / tPr;
      RitPrratio = std::min(RitPrratio, 1.-constants::dsmall);
      evisc[ijk] = fac * std::sqrt(evisc[ijk]) * std::sqrt(1.-RitPrratio);
    }

  for(int k=grid->kstart+1; k<grid->kend; ++k)
  {
    // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
    mlen0 = cs*std::pow(dx*dy*dz[k], 1./3.);
    mlen  = std::pow(1./(1./std::pow(mlen0, n) + 1./(std::pow(constants::kappa*(z[k]+z0m), n))), 1./n);
    fac   = std::pow(mlen, 2.);

    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        // Add the buoyancy production to the TKE
        RitPrratio = N2[ijk] / evisc[ijk] / tPr;
        RitPrratio = std::min(RitPrratio, 1.-constants::dsmall);
        evisc[ijk] = fac * std::sqrt(evisc[ijk]) * std::sqrt(1.-RitPrratio);
      }
  }

  grid->boundary_cyclic(evisc);

  return 0;
}

int cdiff_les2s::evisc_neutral(double * restrict evisc,
                               double * restrict u, double * restrict v, double * restrict w,
                               double * restrict ufluxbot, double * restrict vfluxbot,
                               double * restrict z, double * restrict dz, double z0m)
{
  int    ij,ijk,jj,kk,kstart;
  double dx,dy;

  // wall damping
  double mlen,mlen0,fac;
  const double n  = 2.;

  double RitPrratio;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;
  kstart = grid->kstart;

  dx = grid->dx;
  dy = grid->dy;

  // local copies to aid vectorization
  double tPr = this->tPr;
  double cs  = this->cs;

  for(int k=grid->kstart; k<grid->kend; ++k)
  {
    // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
    mlen0 = cs*std::pow(dx*dy*dz[k], 1./3.);
    mlen  = std::pow(1./(1./std::pow(mlen0, n) + 1./(std::pow(constants::kappa*(z[k]+z0m), n))), 1./n);
    fac   = std::pow(mlen, 2.);

    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        evisc[ijk] = fac * std::sqrt(evisc[ijk]);
      }
  }

  grid->boundary_cyclic(evisc);

  return 0;
}

int cdiff_les2s::diffu(double * restrict ut, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi, double * restrict dzhi, double * restrict evisc, double * restrict fluxbot, double * restrict fluxtop, double * restrict rhoref, double * restrict rhorefh)
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
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
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
            + (  rhorefh[kstart+1] * evisct*((u[ijk+kk]-u[ijk   ])* dzhi[kstart+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
               + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
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
              + (  rhorefh[k+1] * evisct*((u[ijk+kk]-u[ijk   ])* dzhi[k+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
                 - rhorefh[k  ] * eviscb*((u[ijk   ]-u[ijk-kk])* dzhi[k  ] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) / rhoref[k] * dzi[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
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
            + (- rhorefh[kend  ] * fluxtop[ij]
               - rhorefh[kend-1] * eviscb*((u[ijk   ]-u[ijk-kk])* dzhi[kend-1] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) / rhoref[kend-1] * dzi[kend-1];
    }

  return 0;
}

int cdiff_les2s::diffv(double * restrict vt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi, double * restrict dzhi, double * restrict evisc, double * restrict fluxbot, double * restrict fluxtop, double * restrict rhoref, double * restrict rhorefh)
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
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
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
               - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
            // dv/dz + dw/dy
            + (  rhorefh[kstart+1] * evisct*((v[ijk+kk]-v[ijk   ])*dzhi[kstart+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
               + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
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
                 - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
              // dv/dz + dw/dy
              + (  rhorefh[k+1] * evisct*((v[ijk+kk]-v[ijk   ])*dzhi[k+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
                 - rhorefh[k  ] * eviscb*((v[ijk   ]-v[ijk-kk])*dzhi[k  ] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) / rhoref[k] * dzi[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
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
               - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
            // dv/dz + dw/dy
            + (- rhorefh[kend  ] * fluxtop[ij]
               - rhorefh[kend-1] * eviscb*((v[ijk   ]-v[ijk-kk])*dzhi[kend-1] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) / rhoref[kend-1] * dzi[kend-1];
    }

  return 0;
}

int cdiff_les2s::diffw(double * restrict wt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi, double * restrict dzhi, double * restrict evisc, double * restrict rhoref, double * restrict rhorefh)
{
  int    ijk,ii,jj,kk;
  double dxi,dyi;
  double evisce, eviscw, eviscn, eviscs;

  ii = 1;
  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dxi = 1./grid->dx;
  dyi = 1./grid->dy;

  for(int k=grid->kstart+1; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
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
              + (  rhoref[k  ] * evisc[ijk   ]*(w[ijk+kk]-w[ijk   ])*dzi[k  ]
                 - rhoref[k-1] * evisc[ijk-kk]*(w[ijk   ]-w[ijk-kk])*dzi[k-1] ) / rhorefh[k] * 2.* dzhi[k];
      }

  return 0;
}

int cdiff_les2s::diffc(double * restrict at, double * restrict a, double * restrict dzi, double * restrict dzhi, double * restrict evisc, double * restrict fluxbot, double * restrict fluxtop, double * restrict rhoref, double * restrict rhorefh, double tPr)
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
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
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
            + (  rhorefh[kstart+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[kstart+1]
               + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
    }

  for(int k=grid->kstart+1; k<grid->kend-1; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
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
              + (  rhorefh[k+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[k+1]
                 - rhorefh[k  ] * eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[k]  ) / rhoref[k] * dzi[k];
      }

  // top boundary
  for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
    for(int i=grid->istart; i<grid->iend; ++i)
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
            + (- rhorefh[kend  ] * fluxtop[ij]
               - rhorefh[kend-1] * eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[kend-1] ) / rhoref[kend-1] * dzi[kend-1];
    }

  return 0;
}

double cdiff_les2s::calcdnmul(double * restrict evisc, double * restrict dzi, double tPr)
{
  int    ijk,jj,kk;
  double dxidxi,dyidyi;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  dxidxi = 1./(grid->dx * grid->dx);
  dyidyi = 1./(grid->dy * grid->dy);

  double tPrfac = std::min(1., tPr);
  double dnmul = 0;

  // get the maximum time step for diffusion
  for(int k=grid->kstart; k<grid->kend; ++k)
    for(int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
      for(int i=grid->istart; i<grid->iend; ++i)
      {
        ijk = i + j*jj + k*kk;
        dnmul = std::max(dnmul, std::abs(tPrfac*evisc[ijk]*(dxidxi + dyidyi + dzi[k]*dzi[k])));
      }

  grid->getmax(&dnmul);

  return dnmul;
}

inline double cdiff_les2s::phim(double zeta)
{
  double phim;
  if(zeta <= 0.)
  {
    // Businger-Dyer functions
    //x     = (1. - 16. * zeta) ** (0.25)
    //psim  = 3.14159265 / 2. - 2. * arctan(x) + log( (1.+x) ** 2. * (1. + x ** 2.) / 8.)
    // Wilson functions
    phim = std::pow(1. + 3.6*std::pow(std::abs(zeta), 2./3.), -1./2.);
  }
  else
    phim = 1. + 5.*zeta;

  return phim;
}

inline double cdiff_les2s::phih(double zeta)
{
  double phih;
  if(zeta <= 0.)
  {
    // Businger-Dyer functions
    // x     = (1. - 16. * zeta) ** (0.25)
    // psih  = 2. * log( (1. + x ** 2.) / 2. )
    // Wilson functions
    phih = std::pow(1. + 7.9*std::pow(std::abs(zeta), 2./3.), -1./2.);
  }
  else
    phih = 1. + 5.*zeta;

  return phih;
}
