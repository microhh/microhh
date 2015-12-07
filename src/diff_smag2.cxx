/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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
#include "diff_smag2.h"
#include "boundary_surface.h"
#include "defines.h"
#include "constants.h"
#include "thermo.h"
#include "model.h"
#include "monin_obukhov.h"

namespace
{
    namespace most = Monin_obukhov;
}

Diff_smag_2::Diff_smag_2(Model* modelin, Input* inputin) : Diff(modelin, inputin)
{
    swdiff = "smag2";

    #ifdef USECUDA
    mlen_g = 0;
    #endif

    fields->init_diagnostic_field("evisc", "Eddy viscosity", "m2 s-1");

    int nerror = 0;
    nerror += inputin->get_item(&dnmax, "diff", "dnmax", "", 0.5  );
    nerror += inputin->get_item(&cs   , "diff", "cs"   , "", 0.23 );
    nerror += inputin->get_item(&tPr  , "diff", "tPr"  , "", 1./3.);

    if (nerror)
        throw 1;
}

Diff_smag_2::~Diff_smag_2()
{
#ifdef USECUDA
    clear_device();
#endif
}

#ifndef USECUDA
unsigned long Diff_smag_2::get_time_limit(const unsigned long idt, const double dt)
{
    double dnmul = calc_dnmul(fields->sd["evisc"]->data, grid->dzi, this->tPr);
    // Avoid zero division.
    dnmul = std::max(Constants::dsmall, dnmul);

    return idt * dnmax/(dnmul*dt);
}
#endif

#ifndef USECUDA
double Diff_smag_2::get_dn(const double dt)
{
    // calculate eddy viscosity
    const double dnmul = calc_dnmul(fields->sd["evisc"]->data, grid->dzi, this->tPr);

    return dnmul*dt;
}
#endif

#ifndef USECUDA

void Diff_smag_2::exec_viscosity()
{
    // Do a cast because the base boundary class does not have the MOST related variables.
    Boundary_surface* boundaryptr = static_cast<Boundary_surface*>(model->boundary);

    calc_strain2(fields->sd["evisc"]->data,
                 fields->u->data, fields->v->data, fields->w->data,
                 fields->u->datafluxbot, fields->v->datafluxbot,
                 boundaryptr->ustar, boundaryptr->obuk,
                 grid->z, grid->dzi, grid->dzhi);

    // start with retrieving the stability information
    if (model->thermo->get_switch() == "0")
    {
        calc_evisc_neutral(fields->sd["evisc"]->data,
                           fields->u->data, fields->v->data, fields->w->data,
                           fields->u->datafluxbot, fields->v->datafluxbot,
                           grid->z, grid->dz, boundaryptr->z0m);
    }
    // assume buoyancy calculation is needed
    else
    {
        // store the buoyancyflux in tmp1
        model->thermo->get_buoyancy_fluxbot(fields->atmp["tmp1"]);
        // retrieve the full field in tmp1 and use tmp2 for temporary calculations
        model->thermo->get_thermo_field(fields->atmp["tmp1"], fields->atmp["tmp2"], "N2");
        // model->thermo->getThermoField(fields->sd["tmp1"], fields->sd["tmp2"], "b");

        calc_evisc(fields->sd["evisc"]->data,
                   fields->u->data, fields->v->data, fields->w->data, fields->atmp["tmp1"]->data,
                   fields->u->datafluxbot, fields->v->datafluxbot, fields->atmp["tmp1"]->datafluxbot,
                   boundaryptr->ustar, boundaryptr->obuk,
                   grid->z, grid->dz, grid->dzi,
                   boundaryptr->z0m);
    }
}
#endif

#ifndef USECUDA
void Diff_smag_2::exec()
{
    diff_u(fields->ut->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi, grid->dzhi, fields->sd["evisc"]->data,
           fields->u->datafluxbot, fields->u->datafluxtop, fields->rhoref, fields->rhorefh);
    diff_v(fields->vt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi, grid->dzhi, fields->sd["evisc"]->data,
           fields->v->datafluxbot, fields->v->datafluxtop, fields->rhoref, fields->rhorefh);
    diff_w(fields->wt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi, grid->dzhi, fields->sd["evisc"]->data,
           fields->rhoref, fields->rhorefh);

    for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); ++it)
        diff_c(it->second->data, fields->sp[it->first]->data, grid->dzi, grid->dzhi, fields->sd["evisc"]->data,
               fields->sp[it->first]->datafluxbot, fields->sp[it->first]->datafluxtop, fields->rhoref, fields->rhorefh, this->tPr);
}
#endif

void Diff_smag_2::calc_strain2(double* restrict strain2,
                               double* restrict u, double* restrict v, double* restrict w,
                               double* restrict ufluxbot, double* restrict vfluxbot,
                               double* restrict ustar, double* restrict obuk,
                               double* restrict z, double* restrict dzi, double* restrict dzhi)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            strain2[ijk] = 2.*(
                           // du/dz
                           + 0.5*std::pow(-0.5*(ufluxbot[ij]+ufluxbot[ij+ii])/(Constants::kappa*z[kstart]*ustar[ij])*most::phim(z[kstart]/obuk[ij]), 2)
                           // dv/dz
                           + 0.5*std::pow(-0.5*(vfluxbot[ij]+vfluxbot[ij+jj])/(Constants::kappa*z[kstart]*ustar[ij])*most::phim(z[kstart]/obuk[ij]), 2) );

            // add a small number to avoid zero divisions
            strain2[ijk] += Constants::dsmall;
        }

    for (int k=grid->kstart+1; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                strain2[ijk] = 2.*(
                               // du/dx + du/dx
                               + std::pow((u[ijk+ii]-u[ijk])*dxi, 2)

                               // dv/dy + dv/dy
                               + std::pow((v[ijk+jj]-v[ijk])*dyi, 2)

                               // dw/dz + dw/dz
                               + std::pow((w[ijk+kk]-w[ijk])*dzi[k], 2)

                               // du/dy + dv/dx
                               + 0.125*std::pow((u[ijk      ]-u[ijk   -jj])*dyi  + (v[ijk      ]-v[ijk-ii   ])*dxi, 2)
                               + 0.125*std::pow((u[ijk+ii   ]-u[ijk+ii-jj])*dyi  + (v[ijk+ii   ]-v[ijk      ])*dxi, 2)
                               + 0.125*std::pow((u[ijk   +jj]-u[ijk      ])*dyi  + (v[ijk   +jj]-v[ijk-ii+jj])*dxi, 2)
                               + 0.125*std::pow((u[ijk+ii+jj]-u[ijk+ii   ])*dyi  + (v[ijk+ii+jj]-v[ijk   +jj])*dxi, 2)

                               // du/dz + dw/dx
                               + 0.125*std::pow((u[ijk      ]-u[ijk   -kk])*dzhi[k  ] + (w[ijk      ]-w[ijk-ii   ])*dxi, 2)
                               + 0.125*std::pow((u[ijk+ii   ]-u[ijk+ii-kk])*dzhi[k  ] + (w[ijk+ii   ]-w[ijk      ])*dxi, 2)
                               + 0.125*std::pow((u[ijk   +kk]-u[ijk      ])*dzhi[k+1] + (w[ijk   +kk]-w[ijk-ii+kk])*dxi, 2)
                               + 0.125*std::pow((u[ijk+ii+kk]-u[ijk+ii   ])*dzhi[k+1] + (w[ijk+ii+kk]-w[ijk   +kk])*dxi, 2)

                               // dv/dz + dw/dy
                               + 0.125*std::pow((v[ijk      ]-v[ijk   -kk])*dzhi[k  ] + (w[ijk      ]-w[ijk-jj   ])*dyi, 2)
                               + 0.125*std::pow((v[ijk+jj   ]-v[ijk+jj-kk])*dzhi[k  ] + (w[ijk+jj   ]-w[ijk      ])*dyi, 2)
                               + 0.125*std::pow((v[ijk   +kk]-v[ijk      ])*dzhi[k+1] + (w[ijk   +kk]-w[ijk-jj+kk])*dyi, 2)
                               + 0.125*std::pow((v[ijk+jj+kk]-v[ijk+jj   ])*dzhi[k+1] + (w[ijk+jj+kk]-w[ijk   +kk])*dyi, 2) );

                       // Add a small number to avoid zero divisions.
                       strain2[ijk] += Constants::dsmall;
            }
}

void Diff_smag_2::calc_evisc(double* restrict evisc,
                             double* restrict u, double* restrict v, double* restrict w,  double* restrict N2,
                             double* restrict ufluxbot, double* restrict vfluxbot, double* restrict bfluxbot,
                             double* restrict ustar, double* restrict obuk,
                             double* restrict z, double* restrict dz, double* restrict dzi,
                             const double z0m)
{
    // Variables for the wall damping.
    double mlen,mlen0,fac;
    const double n = 2.;

    double RitPrratio;

    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;

    const double dx = grid->dx;
    const double dy = grid->dy;

    // bottom boundary, here strain is fully parametrized using MO
    // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
    mlen0 = this->cs*std::pow(dx*dy*dz[kstart], 1./3.);
    mlen  = std::pow(1./(1./std::pow(mlen0, n) + 1./(std::pow(Constants::kappa*(z[kstart]+z0m), n))), 1./n);
    fac   = std::pow(mlen, 2);

    // local copies to aid vectorization
    double tPr = this->tPr;
    double cs  = this->cs;

    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            // TODO use the thermal expansion coefficient from the input later, what to do if there is no buoyancy?
            // Add the buoyancy production to the TKE
            RitPrratio = -bfluxbot[ij]/(Constants::kappa*z[kstart]*ustar[ij])*most::phih(z[kstart]/obuk[ij]) / evisc[ijk] / tPr;
            RitPrratio = std::min(RitPrratio, 1.-Constants::dsmall);
            evisc[ijk] = fac * std::sqrt(evisc[ijk]) * std::sqrt(1.-RitPrratio);
        }

    for (int k=grid->kstart+1; k<grid->kend; ++k)
    {
        // calculate smagorinsky constant times filter width squared, use wall damping according to Mason
        mlen0 = cs*std::pow(dx*dy*dz[k], 1./3.);
        mlen  = std::pow(1./(1./std::pow(mlen0, n) + 1./(std::pow(Constants::kappa*(z[k]+z0m), n))), 1./n);
        fac   = std::pow(mlen, 2);

        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                // Add the buoyancy production to the TKE
                RitPrratio = N2[ijk] / evisc[ijk] / tPr;
                RitPrratio = std::min(RitPrratio, 1.-Constants::dsmall);
                evisc[ijk] = fac * std::sqrt(evisc[ijk]) * std::sqrt(1.-RitPrratio);
            }
    }

    grid->boundary_cyclic(evisc);
}

void Diff_smag_2::calc_evisc_neutral(double* restrict evisc,
                                     double* restrict u, double* restrict v, double* restrict w,
                                     double* restrict ufluxbot, double* restrict vfluxbot,
                                     double* restrict z, double* restrict dz, const double z0m)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    // Make local copies to aid vectorization.
    const double dx = grid->dx;
    const double dy = grid->dy;
    const double cs = this->cs;

    // Wall damping constant.
    const int n = 2;

    for (int k=grid->kstart; k<grid->kend; ++k)
    {
        // Calculate smagorinsky constant times filter width squared, use wall damping according to Mason's paper.
        const double mlen0 = cs*std::pow(dx*dy*dz[k], 1./3.);
        const double mlen  = std::pow(1./(1./std::pow(mlen0, n) + 1./(std::pow(Constants::kappa*(z[k]+z0m), n))), 1./n);
        const double fac   = std::pow(mlen, 2);

        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                evisc[ijk] = fac * std::sqrt(evisc[ijk]);
            }
    }

    grid->boundary_cyclic(evisc);
}

void Diff_smag_2::diff_u(double* restrict ut, double* restrict u, double* restrict v, double* restrict w,
                         double* restrict dzi, double* restrict dzhi, double* restrict evisc,
                         double* restrict fluxbot, double* restrict fluxtop,
                         double* restrict rhoref, double* restrict rhorefh)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    double eviscn, eviscs, eviscb, evisct;

    // bottom boundary
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            eviscn = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
            eviscs = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
            evisct = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
            eviscb = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);
            ut[ijk] +=
                     // du/dx + du/dx
                     + ( evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
                       - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
                     // du/dy + dv/dx
                     + ( eviscn*((u[ijk+jj]-u[ijk   ])*dyi + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                       - eviscs*((u[ijk   ]-u[ijk-jj])*dyi + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                     // du/dz + dw/dx
                     + ( rhorefh[kstart+1] * evisct*((u[ijk+kk]-u[ijk   ])* dzhi[kstart+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
                       + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
        }

    for (int k=grid->kstart+1; k<grid->kend-1; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                eviscn = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
                eviscs = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
                evisct = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
                eviscb = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);
                ut[ijk] +=
                         // du/dx + du/dx
                         + ( evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
                           - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
                         // du/dy + dv/dx
                         + ( eviscn*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                           - eviscs*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                         // du/dz + dw/dx
                         + ( rhorefh[k+1] * evisct*((u[ijk+kk]-u[ijk   ])* dzhi[k+1] + (w[ijk+kk]-w[ijk-ii+kk])*dxi)
                           - rhorefh[k  ] * eviscb*((u[ijk   ]-u[ijk-kk])* dzhi[k  ] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) / rhoref[k] * dzi[k];
            }

    // top boundary
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + (kend-1)*kk;
            eviscn = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+jj] + evisc[ijk+jj]);
            eviscs = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-jj] + evisc[ijk-ii   ] + evisc[ijk   ]);
            evisct = 0.25*(evisc[ijk-ii   ] + evisc[ijk   ] + evisc[ijk-ii+kk] + evisc[ijk+kk]);
            eviscb = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-kk] + evisc[ijk-ii   ] + evisc[ijk   ]);
            ut[ijk] +=
                     // du/dx + du/dx
                     + ( evisc[ijk   ]*(u[ijk+ii]-u[ijk   ])*dxi
                       - evisc[ijk-ii]*(u[ijk   ]-u[ijk-ii])*dxi ) * 2.* dxi
                     // du/dy + dv/dx
                     + ( eviscn*((u[ijk+jj]-u[ijk   ])*dyi  + (v[ijk+jj]-v[ijk-ii+jj])*dxi)
                       - eviscs*((u[ijk   ]-u[ijk-jj])*dyi  + (v[ijk   ]-v[ijk-ii   ])*dxi) ) * dyi
                     // du/dz + dw/dx
                     + (- rhorefh[kend  ] * fluxtop[ij]
                        - rhorefh[kend-1] * eviscb*((u[ijk   ]-u[ijk-kk])* dzhi[kend-1] + (w[ijk   ]-w[ijk-ii   ])*dxi) ) / rhoref[kend-1] * dzi[kend-1];
        }
}

void Diff_smag_2::diff_v(double* restrict vt, double* restrict u, double* restrict v, double* restrict w,
                         double* restrict dzi, double* restrict dzhi, double* restrict evisc,
                         double* restrict fluxbot, double* restrict fluxtop,
                         double* restrict rhoref, double* restrict rhorefh)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    double evisce, eviscw, eviscb, evisct;

    // bottom boundary
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            evisce = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
            eviscw = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
            evisct = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
            eviscb = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);
            vt[ijk] +=
                     // dv/dx + du/dy
                     + ( evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                       - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                     // dv/dy + dv/dy
                     + ( evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
                       - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
                     // dv/dz + dw/dy
                     + ( rhorefh[kstart+1] * evisct*((v[ijk+kk]-v[ijk   ])*dzhi[kstart+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
                       + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
        }

    for (int k=grid->kstart+1; k<grid->kend-1; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                evisce = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
                eviscw = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
                evisct = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
                eviscb = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);
                vt[ijk] +=
                         // dv/dx + du/dy
                         + ( evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                           - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                         // dv/dy + dv/dy
                         + ( evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
                           - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
                         // dv/dz + dw/dy
                         + ( rhorefh[k+1] * evisct*((v[ijk+kk]-v[ijk   ])*dzhi[k+1] + (w[ijk+kk]-w[ijk-jj+kk])*dyi)
                           - rhorefh[k  ] * eviscb*((v[ijk   ]-v[ijk-kk])*dzhi[k  ] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) / rhoref[k] * dzi[k];
            }

    // top boundary
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + (kend-1)*kk;
            evisce = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+ii-jj] + evisc[ijk+ii]);
            eviscw = 0.25*(evisc[ijk-ii-jj] + evisc[ijk-ii] + evisc[ijk   -jj] + evisc[ijk   ]);
            evisct = 0.25*(evisc[ijk   -jj] + evisc[ijk   ] + evisc[ijk+kk-jj] + evisc[ijk+kk]);
            eviscb = 0.25*(evisc[ijk-kk-jj] + evisc[ijk-kk] + evisc[ijk   -jj] + evisc[ijk   ]);
            vt[ijk] +=
                     // dv/dx + du/dy
                     + ( evisce*((v[ijk+ii]-v[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-jj])*dyi)
                       - eviscw*((v[ijk   ]-v[ijk-ii])*dxi + (u[ijk   ]-u[ijk   -jj])*dyi) ) * dxi
                     // dv/dy + dv/dy
                     + ( evisc[ijk   ]*(v[ijk+jj]-v[ijk   ])*dyi
                       - evisc[ijk-jj]*(v[ijk   ]-v[ijk-jj])*dyi ) * 2.* dyi
                     // dv/dz + dw/dy
                     + (- rhorefh[kend  ] * fluxtop[ij]
                        - rhorefh[kend-1] * eviscb*((v[ijk   ]-v[ijk-kk])*dzhi[kend-1] + (w[ijk   ]-w[ijk-jj   ])*dyi) ) / rhoref[kend-1] * dzi[kend-1];
        }
}

void Diff_smag_2::diff_w(double* restrict wt, double* restrict u, double* restrict v, double* restrict w,
                         double* restrict dzi, double* restrict dzhi, double* restrict evisc,
                         double* restrict rhoref, double* restrict rhorefh)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    double evisce, eviscw, eviscn, eviscs;

    for (int k=grid->kstart+1; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                evisce = 0.25*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+ii-kk] + evisc[ijk+ii]);
                eviscw = 0.25*(evisc[ijk-ii-kk] + evisc[ijk-ii] + evisc[ijk   -kk] + evisc[ijk   ]);
                eviscn = 0.25*(evisc[ijk   -kk] + evisc[ijk   ] + evisc[ijk+jj-kk] + evisc[ijk+jj]);
                eviscs = 0.25*(evisc[ijk-jj-kk] + evisc[ijk-jj] + evisc[ijk   -kk] + evisc[ijk   ]);
                wt[ijk] +=
                         // dw/dx + du/dz
                         + ( evisce*((w[ijk+ii]-w[ijk   ])*dxi + (u[ijk+ii]-u[ijk+ii-kk])*dzhi[k])
                           - eviscw*((w[ijk   ]-w[ijk-ii])*dxi + (u[ijk   ]-u[ijk+  -kk])*dzhi[k]) ) * dxi
                         // dw/dy + dv/dz
                         + ( eviscn*((w[ijk+jj]-w[ijk   ])*dyi + (v[ijk+jj]-v[ijk+jj-kk])*dzhi[k])
                           - eviscs*((w[ijk   ]-w[ijk-jj])*dyi + (v[ijk   ]-v[ijk+  -kk])*dzhi[k]) ) * dyi
                         // dw/dz + dw/dz
                         + ( rhoref[k  ] * evisc[ijk   ]*(w[ijk+kk]-w[ijk   ])*dzi[k  ]
                           - rhoref[k-1] * evisc[ijk-kk]*(w[ijk   ]-w[ijk-kk])*dzi[k-1] ) / rhorefh[k] * 2.* dzhi[k];
            }
}

void Diff_smag_2::diff_c(double* restrict at, double* restrict a,
                         double* restrict dzi, double* restrict dzhi, double* restrict evisc,
                         double* restrict fluxbot, double* restrict fluxtop, 
                         double* restrict rhoref, double* restrict rhorefh, double tPr)
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;
    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    const double dxidxi = 1./(grid->dx * grid->dx);
    const double dyidyi = 1./(grid->dy * grid->dy);

    double evisce,eviscw,eviscn,eviscs,evisct,eviscb;

    // bottom boundary
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + kstart*kk;
            evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])/tPr;
            eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])/tPr;
            eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])/tPr;
            eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])/tPr;
            evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])/tPr;
            eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])/tPr;

            at[ijk] +=
                     + ( evisce*(a[ijk+ii]-a[ijk   ]) 
                       - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi 
                     + ( eviscn*(a[ijk+jj]-a[ijk   ]) 
                       - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                     + ( rhorefh[kstart+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[kstart+1]
                       + rhorefh[kstart  ] * fluxbot[ij] ) / rhoref[kstart] * dzi[kstart];
        }

    for (int k=grid->kstart+1; k<grid->kend-1; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])/tPr;
                eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])/tPr;
                eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])/tPr;
                eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])/tPr;
                evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])/tPr;
                eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])/tPr;

                at[ijk] +=
                         + ( evisce*(a[ijk+ii]-a[ijk   ]) 
                           - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi 
                         + ( eviscn*(a[ijk+jj]-a[ijk   ]) 
                           - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                         + ( rhorefh[k+1] * evisct*(a[ijk+kk]-a[ijk   ])*dzhi[k+1]
                           - rhorefh[k  ] * eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[k]  ) / rhoref[k] * dzi[k];
            }

    // top boundary
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ij  = i + j*jj;
            const int ijk = i + j*jj + (kend-1)*kk;
            evisce = 0.5*(evisc[ijk   ]+evisc[ijk+ii])/tPr;
            eviscw = 0.5*(evisc[ijk-ii]+evisc[ijk   ])/tPr;
            eviscn = 0.5*(evisc[ijk   ]+evisc[ijk+jj])/tPr;
            eviscs = 0.5*(evisc[ijk-jj]+evisc[ijk   ])/tPr;
            evisct = 0.5*(evisc[ijk   ]+evisc[ijk+kk])/tPr;
            eviscb = 0.5*(evisc[ijk-kk]+evisc[ijk   ])/tPr;

            at[ijk] +=
                     + ( evisce*(a[ijk+ii]-a[ijk   ]) 
                       - eviscw*(a[ijk   ]-a[ijk-ii]) ) * dxidxi 
                     + ( eviscn*(a[ijk+jj]-a[ijk   ]) 
                       - eviscs*(a[ijk   ]-a[ijk-jj]) ) * dyidyi
                     + (-rhorefh[kend  ] * fluxtop[ij]
                       - rhorefh[kend-1] * eviscb*(a[ijk   ]-a[ijk-kk])*dzhi[kend-1] ) / rhoref[kend-1] * dzi[kend-1];
        }
}

double Diff_smag_2::calc_dnmul(double* restrict evisc, double* restrict dzi, double tPr)
{
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    const double dxidxi = 1./(grid->dx * grid->dx);
    const double dyidyi = 1./(grid->dy * grid->dy);

    const double tPrfac = std::min(1., tPr);
    double dnmul = 0;

    // get the maximum time step for diffusion
    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                dnmul = std::max(dnmul, std::abs(tPrfac*evisc[ijk]*(dxidxi + dyidyi + dzi[k]*dzi[k])));
            }

    grid->get_max(&dnmul);

    return dnmul;
}
