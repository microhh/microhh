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
#include "advec_2i4.h"
#include "defines.h"
#include "constants.h"
#include "fd.h"
#include "model.h"

using namespace fd::o4;
using namespace fd::o2;

Advec_2i4::Advec_2i4(Model* modelin, Input* inputin) : Advec(modelin, inputin)
{
    const int igc = 2;
    const int jgc = 2;
    const int kgc = 2;

    grid->set_minimum_ghost_cells(igc, jgc, kgc);
}

Advec_2i4::~Advec_2i4()
{
}

#ifndef USECUDA
unsigned long Advec_2i4::get_time_limit(unsigned long idt, double dt)
{
    double cfl = calc_cfl(fields->u->data, fields->v->data, fields->w->data, grid->dzi, dt);
    // Avoid zero divisons.
    cfl = std::max(cflmin, cfl);
    return idt * cflmax / cfl;
}
#endif

#ifndef USECUDA
double Advec_2i4::get_cfl(double dt)
{
    return calc_cfl(fields->u->data, fields->v->data, fields->w->data, grid->dzi, dt);
}
#endif

#ifndef USECUDA
void Advec_2i4::exec()
{
    advec_u(fields->ut->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi, 
            fields->rhoref, fields->rhorefh);
    advec_v(fields->vt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi,
            fields->rhoref, fields->rhorefh);
    advec_w(fields->wt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzhi,
            fields->rhoref, fields->rhorefh);

    for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
        advec_s(it->second->data, fields->sp[it->first]->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi,
                fields->rhoref, fields->rhorefh);

}
#endif

double Advec_2i4::calc_cfl(double* restrict u, double* restrict v, double* restrict w, double* restrict dzi, double dt)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int jj1 = 1*grid->icells;
    const int jj2 = 2*grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    double cfl = 0;

    int k = kstart;
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            cfl = std::max(cfl, std::abs(interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]))*dxi 
                              + std::abs(interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]))*dyi 
                              + std::abs(interp2(w[ijk    ], w[ijk+kk1]))*dzi[k]);
        }

    for (k=grid->kstart+1; k<grid->kend-1; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                cfl = std::max(cfl, std::abs(interp4(u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2]))*dxi 
                                  + std::abs(interp4(v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2]))*dyi 
                                  + std::abs(interp4(w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2]))*dzi[k]);
            }

    k = kend-1;
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk  = i + j*jj1 + k*kk1;
            cfl = std::max(cfl, std::abs(interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]))*dxi 
                              + std::abs(interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]))*dyi 
                              + std::abs(interp2(w[ijk    ], w[ijk+kk1]))*dzi[k]);
        }

    grid->getMax(&cfl);

    cfl = cfl*dt;

    return cfl;
}

void Advec_2i4::advec_u(double* restrict ut, double* restrict u, double* restrict v, double* restrict w, 
                        double* restrict dzi, double* restrict rhoref, double* restrict rhorefh)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int jj1 = 1*grid->icells;
    const int jj2 = 2*grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    int k = kstart; 

    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            ut[ijk] += 
                     // u*du/dx
                     - ( interp2(u[ijk        ], u[ijk+ii1]) * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                       - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

                     // v*du/dy
                     - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                       - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 

                     // w*du/dz -> second order interpolation for fluxtop, fluxbot = 0. as w=0
                     - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp2(u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k];
        }

    k = kstart + 1; 
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            ut[ijk] += 
                     // u*du/dx
                     - ( interp2(u[ijk        ], u[ijk+ii1]) * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                       - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

                     // v*du/dy
                     - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                       - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 

                     // w*du/dz -> second order interpolation for fluxbot
                     - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp4(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2])
                       - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp2(u[ijk-kk1], u[ijk    ]) ) / rhoref[k] * dzi[k];
        }

    for (k=grid->kstart+2; k<grid->kend-2; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                ut[ijk] += 
                         // u*du/dx
                         - ( interp2(u[ijk        ], u[ijk+ii1]) * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                           - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

                         // v*du/dy
                         - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                           - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 

                         // w*du/dz
                         - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp4(u[ijk-kk1], u[ijk    ], u[ijk+kk1], u[ijk+kk2])
                           - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp4(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }

    k = kend - 2; 
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            ut[ijk] += 
                     // u*du/dx
                     - ( interp2(u[ijk        ], u[ijk+ii1]) * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                       - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

                     // v*du/dy
                     - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                       - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 

                     // w*du/dz -> second order interpolation for fluxtop
                     - ( rhorefh[k+1] * interp2(w[ijk-ii1+kk1], w[ijk+kk1]) * interp2(u[ijk    ], u[ijk+kk1])
                       - rhorefh[k  ] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp4(u[ijk-kk2], u[ijk-kk1], u[ijk    ], u[ijk+kk1]) ) / rhoref[k] * dzi[k];
        }

    k = kend - 1; 
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            ut[ijk] += 
                     // u*du/dx
                     - ( interp2(u[ijk        ], u[ijk+ii1]) * interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2])
                       - interp2(u[ijk-ii1    ], u[ijk    ]) * interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) ) * dxi

                     // v*du/dy
                     - ( interp2(v[ijk-ii1+jj1], v[ijk+jj1]) * interp4(u[ijk-jj1], u[ijk    ], u[ijk+jj1], u[ijk+jj2])
                       - interp2(v[ijk-ii1    ], v[ijk    ]) * interp4(u[ijk-jj2], u[ijk-jj1], u[ijk    ], u[ijk+jj1]) ) * dyi 

                     // w*du/dz -> second order interpolation for fluxbot, fluxtop=0 as w=0
                     - ( -rhorefh[k] * interp2(w[ijk-ii1    ], w[ijk    ]) * interp2(u[ijk-kk1], u[ijk    ]) ) / rhoref[k] * dzi[k];
        }
}

void Advec_2i4::advec_v(double* restrict vt, double* restrict u, double* restrict v, double* restrict w,
                        double* restrict dzi, double* restrict rhoref, double* restrict rhorefh)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int jj1 = 1*grid->icells;
    const int jj2 = 2*grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    int k = kstart;
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            vt[ijk] += 
                     // u*dv/dx
                     - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                       - interp2(u[ijk    -jj1], u[ijk    ]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

                     // v*dv/dy
                     - ( interp2(v[ijk        ], v[ijk+jj1]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                       - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

                     // w*dv/dz
                     - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp2(v[ijk    ], v[ijk+kk1]) ) / rhoref[k] * dzi[k];
        }

    k = kstart+1;
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            vt[ijk] += 
                     // u*dv/dx
                     - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                       - interp2(u[ijk    -jj1], u[ijk    ]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

                     // v*dv/dy
                     - ( interp2(v[ijk        ], v[ijk+jj1]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                       - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

                     // w*dv/dz
                     - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp4(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])
                       - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp2(v[ijk-kk1], v[ijk    ]) ) / rhoref[k] * dzi[k];
        }

    for (k=grid->kstart+2; k<grid->kend-2; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                vt[ijk] += 
                         // u*dv/dx
                         - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                           - interp2(u[ijk    -jj1], u[ijk    ]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

                         // v*dv/dy
                         - ( interp2(v[ijk        ], v[ijk+jj1]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                           - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

                         // w*dv/dz
                         - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp4(v[ijk-kk1], v[ijk    ], v[ijk+kk1], v[ijk+kk2])
                           - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp4(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }

    k = kend-2;
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            vt[ijk] += 
                     // u*dv/dx
                     - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                       - interp2(u[ijk    -jj1], u[ijk    ]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

                     // v*dv/dy
                     - ( interp2(v[ijk        ], v[ijk+jj1]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                       - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

                     // w*dv/dz
                     - ( rhorefh[k+1] * interp2(w[ijk-jj1+kk1], w[ijk+kk1]) * interp2(v[ijk    ], v[ijk+kk1])
                       - rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp4(v[ijk-kk2], v[ijk-kk1], v[ijk    ], v[ijk+kk1]) ) / rhoref[k] * dzi[k];
        }

    k = kend-1;
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            vt[ijk] +=
                     // u*dv/dx
                     - ( interp2(u[ijk+ii1-jj1], u[ijk+ii1]) * interp4(v[ijk-ii1], v[ijk    ], v[ijk+ii1], v[ijk+ii2])
                       - interp2(u[ijk    -jj1], u[ijk    ]) * interp4(v[ijk-ii2], v[ijk-ii1], v[ijk    ], v[ijk+ii1]) ) * dxi

                     // v*dv/dy
                     - ( interp2(v[ijk        ], v[ijk+jj1]) * interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2])
                       - interp2(v[ijk-jj1    ], v[ijk    ]) * interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) ) * dyi

                     // w*dv/dz
                     - (- rhorefh[k  ] * interp2(w[ijk-jj1    ], w[ijk    ]) * interp2(v[ijk-kk1], v[ijk    ]) ) / rhoref[k] * dzi[k];
        }
}

void Advec_2i4::advec_w(double* restrict wt, double* restrict u, double* restrict v, double* restrict w,
                        double* restrict dzhi, double* restrict rhoref, double* restrict rhorefh)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int jj1 = 1*grid->icells;
    const int jj2 = 2*grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    int k = kstart+1;
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            wt[ijk] += 
                     // u*dw/dx 
                     - ( interp2(u[ijk+ii1-kk1], u[ijk+ii1]) * interp4(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
                       - interp2(u[ijk    -kk1], u[ijk    ]) * interp4(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi

                     // v*dw/dy 
                     - ( interp2(v[ijk+jj1-kk1], v[ijk+jj1]) * interp4(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
                       - interp2(v[ijk    -kk1], v[ijk    ]) * interp4(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi

                     // w*dw/dz 
                     - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp4(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2])
                       - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp2(w[ijk-kk1], w[ijk    ]) ) / rhorefh[k] * dzhi[k];
        }

    for (k=grid->kstart+2; k<grid->kend-1; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                wt[ijk] +=
                         // u*dw/dx 
                         - ( interp2(u[ijk+ii1-kk1], u[ijk+ii1]) * interp4(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
                           - interp2(u[ijk    -kk1], u[ijk    ]) * interp4(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi

                         // v*dw/dy 
                         - ( interp2(v[ijk+jj1-kk1], v[ijk+jj1]) * interp4(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
                           - interp2(v[ijk    -kk1], v[ijk    ]) * interp4(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi

                         // w*dw/dz 
                         - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp4(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2])
                           - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp4(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) / rhorefh[k] * dzhi[k];
            }

    k = kend-1;
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            wt[ijk] += 
                     // u*dw/dx 
                     - ( interp2(u[ijk+ii1-kk1], u[ijk+ii1]) * interp4(w[ijk-ii1], w[ijk    ], w[ijk+ii1], w[ijk+ii2])
                       - interp2(u[ijk    -kk1], u[ijk    ]) * interp4(w[ijk-ii2], w[ijk-ii1], w[ijk    ], w[ijk+ii1]) ) * dxi

                     // v*dw/dy 
                     - ( interp2(v[ijk+jj1-kk1], v[ijk+jj1]) * interp4(w[ijk-jj1], w[ijk    ], w[ijk+jj1], w[ijk+jj2])
                       - interp2(v[ijk    -kk1], v[ijk    ]) * interp4(w[ijk-jj2], w[ijk-jj1], w[ijk    ], w[ijk+jj1]) ) * dyi

                     // w*dw/dz 
                     - ( rhoref[k  ] * interp2(w[ijk        ], w[ijk+kk1]) * interp2(w[ijk    ], w[ijk+kk1])
                       - rhoref[k-1] * interp2(w[ijk-kk1    ], w[ijk    ]) * interp4(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) ) / rhorefh[k] * dzhi[k];
        }
}

void Advec_2i4::advec_s(double* restrict st, double* restrict s, double* restrict u, double* restrict v, double* restrict w,
                        double* restrict dzi, double* restrict rhoref, double* restrict rhorefh)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int jj1 = 1*grid->icells;
    const int jj2 = 2*grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    // assume that w at the boundary equals zero...
    int k = kstart;
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            st[ijk] += 
                     - ( u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                       - u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                     - ( v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                       - v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 

                     - ( rhorefh[k+1] * w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
        }

    k = kstart+1;
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            st[ijk] += 
                     - ( u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                       - u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                     - ( v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                       - v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 

                     - ( rhorefh[k+1] * w[ijk+kk1] * interp4(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                       - rhorefh[k  ] * w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]) ) / rhoref[k] * dzi[k];
        }

    for (k=grid->kstart+2; k<grid->kend-2; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                st[ijk] += 
                         - ( u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                           - u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                         - ( v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                           - v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 

                         - ( rhorefh[k+1] * w[ijk+kk1] * interp4(s[ijk-kk1], s[ijk    ], s[ijk+kk1], s[ijk+kk2])
                           - rhorefh[k  ] * w[ijk    ] * interp4(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
            }

    k = kend-2;
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            st[ijk] += 
                     - ( u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                       - u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                     - ( v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                       - v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 

                     - ( rhorefh[k+1] * w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1])
                       - rhorefh[k  ] * w[ijk    ] * interp4(s[ijk-kk2], s[ijk-kk1], s[ijk    ], s[ijk+kk1]) ) / rhoref[k] * dzi[k];
        }

    // assume that w at the boundary equals zero...
    k = kend-1;
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + k*kk1;
            st[ijk] += 
                     - ( u[ijk+ii1] * interp4(s[ijk-ii1], s[ijk    ], s[ijk+ii1], s[ijk+ii2])
                       - u[ijk    ] * interp4(s[ijk-ii2], s[ijk-ii1], s[ijk    ], s[ijk+ii1]) ) * dxi

                     - ( v[ijk+jj1] * interp4(s[ijk-jj1], s[ijk    ], s[ijk+jj1], s[ijk+jj2])
                       - v[ijk    ] * interp4(s[ijk-jj2], s[ijk-jj1], s[ijk    ], s[ijk+jj1]) ) * dyi 

                     - (- rhorefh[k  ] * w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]) ) / rhoref[k] * dzi[k];
        }
}
