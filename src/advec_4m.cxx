/*
 * MicroHH
 * Copyright (c) 2011-2017 Chiel van Heerwaarden
 * Copyright (c) 2011-2017 Thijs Heus
 * Copyright (c) 2014-2017 Bart van Stratum
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
#include "advec_4m.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "model.h"

using Finite_difference::O2::interp2;
using Finite_difference::O4::interp4;
using Finite_difference::O4::grad4;
using Finite_difference::O4::grad4x;

Advec_4m::Advec_4m(Model* modelin, Input* inputin) : Advec(modelin, inputin)
{
    swadvec = "4m";
}

Advec_4m::~Advec_4m()
{
}

#ifndef USECUDA
unsigned long Advec_4m::get_time_limit(unsigned long idt, double dt)
{
    // Calculate cfl and prevent zero divisons.
    double cfl = calc_cfl(fields->u->data, fields->v->data, fields->w->data, grid->dzi, dt);
    cfl = std::max(cflmin, cfl);
    return idt * cflmax / cfl;
}

double Advec_4m::get_cfl(double dt)
{
    return calc_cfl(fields->u->data, fields->v->data, fields->w->data, grid->dzi, dt);
}

void Advec_4m::exec()
{
    advec_u(fields->ut->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi4 );
    advec_v(fields->vt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi4 );
    advec_w(fields->wt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzhi4);

    for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); ++it)
        advec_s(it->second->data, fields->sp[it->first]->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi4);
}
#endif

double Advec_4m::calc_cfl(double * restrict u, double * restrict v, double * restrict w, double * restrict dzi, double dt)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int jj1 = 1*grid->icells;
    const int jj2 = 2*grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    double cfl = 0;

    for (int k=grid->kstart; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                cfl = std::max(cfl, std::abs(interp4(u[ijk-ii1], u[ijk], u[ijk+ii1], u[ijk+ii2]))*dxi 
                                  + std::abs(interp4(v[ijk-jj1], v[ijk], v[ijk+jj1], v[ijk+jj2]))*dyi 
                                  + std::abs(interp4(w[ijk-kk1], w[ijk], w[ijk+kk1], w[ijk+kk2]))*dzi[k] );
            }

    grid->get_max(&cfl);

    cfl = cfl*dt;

    return cfl;
}

void Advec_4m::advec_u(double * restrict ut, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi4)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int ii3 = 3;
    const int jj1 = 1*grid->icells;
    const int jj2 = 2*grid->icells;
    const int jj3 = 3*grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;
    const int kk3 = 3*grid->ijcells;

    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    // bottom boundary
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + kstart*kk1;
            ut[ijk] +=
                     - grad4(interp4(u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ]) * interp2(u[ijk-ii3], u[ijk    ]),
                             interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) * interp2(u[ijk-ii1], u[ijk    ]),
                             interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]) * interp2(u[ijk    ], u[ijk+ii1]),
                             interp4(u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]) * interp2(u[ijk    ], u[ijk+ii3]), dxi)

                     - grad4(interp4(v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1]) * interp2(u[ijk-jj3], u[ijk    ]),
                             interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp2(u[ijk-jj1], u[ijk    ]),
                             interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp2(u[ijk    ], u[ijk+jj1]),
                             interp4(v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2]) * interp2(u[ijk    ], u[ijk+jj3]), dyi)

                     // boundary condition
                     - grad4x(-interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp2(u[ijk-kk1], u[ijk+kk2]),
                               interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp2(u[ijk-kk1], u[ijk    ]),
                               interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp2(u[ijk    ], u[ijk+kk1]),
                               interp4(w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2]) * interp2(u[ijk    ], u[ijk+kk3]))
                       * dzi4[kstart];
        }

    for (int k=grid->kstart+1; k<grid->kend-1; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                ut[ijk] +=
                         - grad4(interp4(u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ]) * interp2(u[ijk-ii3], u[ijk    ]),
                                 interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) * interp2(u[ijk-ii1], u[ijk    ]),
                                 interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]) * interp2(u[ijk    ], u[ijk+ii1]),
                                 interp4(u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]) * interp2(u[ijk    ], u[ijk+ii3]), dxi)

                         - grad4(interp4(v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1]) * interp2(u[ijk-jj3], u[ijk    ]),
                                 interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp2(u[ijk-jj1], u[ijk    ]),
                                 interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp2(u[ijk    ], u[ijk+jj1]),
                                 interp4(v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2]) * interp2(u[ijk    ], u[ijk+jj3]), dyi)

                         - grad4x(interp4(w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1]) * interp2(u[ijk-kk3], u[ijk    ]),
                                  interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp2(u[ijk-kk1], u[ijk    ]),
                                  interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp2(u[ijk    ], u[ijk+kk1]),
                                  interp4(w[ijk-ii2+kk2], w[ijk-ii1+kk2], w[ijk+kk2], w[ijk+ii1+kk2]) * interp2(u[ijk    ], u[ijk+kk3]))
                           * dzi4[k];
            }

    // top boundary
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + (kend-1)*kk1;
            ut[ijk] +=
                     - grad4(interp4(u[ijk-ii3], u[ijk-ii2], u[ijk-ii1], u[ijk    ]) * interp2(u[ijk-ii3], u[ijk    ]),
                             interp4(u[ijk-ii2], u[ijk-ii1], u[ijk    ], u[ijk+ii1]) * interp2(u[ijk-ii1], u[ijk    ]),
                             interp4(u[ijk-ii1], u[ijk    ], u[ijk+ii1], u[ijk+ii2]) * interp2(u[ijk    ], u[ijk+ii1]),
                             interp4(u[ijk    ], u[ijk+ii1], u[ijk+ii2], u[ijk+ii3]) * interp2(u[ijk    ], u[ijk+ii3]), dxi)

                     - grad4(interp4(v[ijk-ii2-jj1], v[ijk-ii1-jj1], v[ijk-jj1], v[ijk+ii1-jj1]) * interp2(u[ijk-jj3], u[ijk    ]),
                             interp4(v[ijk-ii2    ], v[ijk-ii1    ], v[ijk    ], v[ijk+ii1    ]) * interp2(u[ijk-jj1], u[ijk    ]),
                             interp4(v[ijk-ii2+jj1], v[ijk-ii1+jj1], v[ijk+jj1], v[ijk+ii1+jj1]) * interp2(u[ijk    ], u[ijk+jj1]),
                             interp4(v[ijk-ii2+jj2], v[ijk-ii1+jj2], v[ijk+jj2], v[ijk+ii1+jj2]) * interp2(u[ijk    ], u[ijk+jj3]), dyi)

                     - grad4x( interp4(w[ijk-ii2-kk1], w[ijk-ii1-kk1], w[ijk-kk1], w[ijk+ii1-kk1]) * interp2(u[ijk-kk3], u[ijk    ]),
                               interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp2(u[ijk-kk1], u[ijk    ]),
                               interp4(w[ijk-ii2+kk1], w[ijk-ii1+kk1], w[ijk+kk1], w[ijk+ii1+kk1]) * interp2(u[ijk    ], u[ijk+kk1]),
                              -interp4(w[ijk-ii2    ], w[ijk-ii1    ], w[ijk    ], w[ijk+ii1    ]) * interp2(u[ijk-kk2], u[ijk+kk1]))
                       * dzi4[kend-1];
        }
}

void Advec_4m::advec_v(double * restrict vt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi4)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int ii3 = 3;
    const int jj1 = 1*grid->icells;
    const int jj2 = 2*grid->icells;
    const int jj3 = 3*grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;
    const int kk3 = 3*grid->ijcells;

    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    // bottom boundary
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + kstart*kk1;
            vt[ijk] +=
                     - grad4(interp4(u[ijk-ii1-jj2], u[ijk-ii1-jj1], u[ijk-ii1], u[ijk-ii1+jj1]) * interp2(v[ijk-ii3], v[ijk    ]),
                             interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp2(v[ijk-ii1], v[ijk    ]),
                             interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp2(v[ijk    ], v[ijk+ii1]),
                             interp4(u[ijk+ii2-jj2], u[ijk+ii2-jj1], u[ijk+ii2], u[ijk+ii2+jj1]) * interp2(v[ijk    ], v[ijk+ii3]), dxi)

                     - grad4(interp4(v[ijk-jj3], v[ijk-jj2], v[ijk-jj1], v[ijk    ]) * interp2(v[ijk-jj3], v[ijk    ]),
                             interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) * interp2(v[ijk-jj1], v[ijk    ]),
                             interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]) * interp2(v[ijk    ], v[ijk+jj1]),
                             interp4(v[ijk    ], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]) * interp2(v[ijk    ], v[ijk+jj3]), dyi)

                     - grad4x(-interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp2(v[ijk-kk1], v[ijk+kk2]),
                               interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp2(v[ijk-kk1], v[ijk    ]),
                               interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp2(v[ijk    ], v[ijk+kk1]),
                               interp4(w[ijk-jj2+kk2], w[ijk-jj1+kk2], w[ijk+kk2], w[ijk+jj1+kk2]) * interp2(v[ijk    ], v[ijk+kk3]))
                       * dzi4[kstart];
        }

    for (int k=grid->kstart+1; k<grid->kend-1; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                vt[ijk] +=
                         - grad4(interp4(u[ijk-ii1-jj2], u[ijk-ii1-jj1], u[ijk-ii1], u[ijk-ii1+jj1]) * interp2(v[ijk-ii3], v[ijk    ]),
                                 interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp2(v[ijk-ii1], v[ijk    ]),
                                 interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp2(v[ijk    ], v[ijk+ii1]),
                                 interp4(u[ijk+ii2-jj2], u[ijk+ii2-jj1], u[ijk+ii2], u[ijk+ii2+jj1]) * interp2(v[ijk    ], v[ijk+ii3]), dxi)

                         - grad4(interp4(v[ijk-jj3], v[ijk-jj2], v[ijk-jj1], v[ijk    ]) * interp2(v[ijk-jj3], v[ijk    ]),
                                 interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) * interp2(v[ijk-jj1], v[ijk    ]),
                                 interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]) * interp2(v[ijk    ], v[ijk+jj1]),
                                 interp4(v[ijk    ], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]) * interp2(v[ijk    ], v[ijk+jj3]), dyi)

                         - grad4x(interp4(w[ijk-jj2-kk1], w[ijk-jj1-kk1], w[ijk-kk1], w[ijk+jj1-kk1]) * interp2(v[ijk-kk3], v[ijk    ]),
                                  interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp2(v[ijk-kk1], v[ijk    ]),
                                  interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp2(v[ijk    ], v[ijk+kk1]),
                                  interp4(w[ijk-jj2+kk2], w[ijk-jj1+kk2], w[ijk+kk2], w[ijk+jj1+kk2]) * interp2(v[ijk    ], v[ijk+kk3]))
                           * dzi4[k];
            }

    // top boundary
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + (kend-1)*kk1;
            vt[ijk] +=
                - grad4(interp4(u[ijk-ii1-jj2], u[ijk-ii1-jj1], u[ijk-ii1], u[ijk-ii1+jj1]) * interp2(v[ijk-ii3], v[ijk    ]),
                        interp4(u[ijk    -jj2], u[ijk    -jj1], u[ijk    ], u[ijk    +jj1]) * interp2(v[ijk-ii1], v[ijk    ]),
                        interp4(u[ijk+ii1-jj2], u[ijk+ii1-jj1], u[ijk+ii1], u[ijk+ii1+jj1]) * interp2(v[ijk    ], v[ijk+ii1]),
                        interp4(u[ijk+ii2-jj2], u[ijk+ii2-jj1], u[ijk+ii2], u[ijk+ii2+jj1]) * interp2(v[ijk    ], v[ijk+ii3]), dxi)

                - grad4(interp4(v[ijk-jj3], v[ijk-jj2], v[ijk-jj1], v[ijk    ]) * interp2(v[ijk-jj3], v[ijk    ]),
                        interp4(v[ijk-jj2], v[ijk-jj1], v[ijk    ], v[ijk+jj1]) * interp2(v[ijk-jj1], v[ijk    ]),
                        interp4(v[ijk-jj1], v[ijk    ], v[ijk+jj1], v[ijk+jj2]) * interp2(v[ijk    ], v[ijk+jj1]),
                        interp4(v[ijk    ], v[ijk+jj1], v[ijk+jj2], v[ijk+jj3]) * interp2(v[ijk    ], v[ijk+jj3]), dyi)

                - grad4x( interp4(w[ijk-jj2-kk1], w[ijk-jj1-kk1], w[ijk-kk1], w[ijk+jj1-kk1]) * interp2(v[ijk-kk3], v[ijk    ]),
                          interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp2(v[ijk-kk1], v[ijk    ]),
                          interp4(w[ijk-jj2+kk1], w[ijk-jj1+kk1], w[ijk+kk1], w[ijk+jj1+kk1]) * interp2(v[ijk    ], v[ijk+kk1]),
                         -interp4(w[ijk-jj2    ], w[ijk-jj1    ], w[ijk    ], w[ijk+jj1    ]) * interp2(v[ijk-kk2], v[ijk+kk1]))
                  * dzi4[kend-1];
        }
}

void Advec_4m::advec_w(double * restrict wt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzhi4)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int ii3 = 3;
    const int jj1 = 1*grid->icells;
    const int jj2 = 2*grid->icells;
    const int jj3 = 3*grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;
    const int kk3 = 3*grid->ijcells;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    /*
    // bottom boundary 
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
for (int i=grid->istart; i<grid->iend; ++i)
{
ijk = i + j*jj1 + (kstart+1)*kk1;
wt[ijk] +=
- grad4(interp4(u[ijk-ii1-kk2], u[ijk-ii1-kk1], u[ijk-ii1], u[ijk-ii1+kk1]) * interp2(w[ijk-ii3], w[ijk    ]),
interp4(u[ijk    -kk2], u[ijk    -kk1], u[ijk    ], u[ijk    +kk1]) * interp2(w[ijk-ii1], w[ijk    ]),
interp4(u[ijk+ii1-kk2], u[ijk+ii1-kk1], u[ijk+ii1], u[ijk+ii1+kk1]) * interp2(w[ijk    ], w[ijk+ii1]),
interp4(u[ijk+ii2-kk2], u[ijk+ii2-kk1], u[ijk+ii2], u[ijk+ii2+kk1]) * interp2(w[ijk    ], w[ijk+ii3]), dxi)

- grad4(interp4(v[ijk-jj1-kk2], v[ijk-jj1-kk1], v[ijk-jj1], v[ijk-jj1+kk1]) * interp2(w[ijk-jj3], w[ijk    ]),
interp4(v[ijk    -kk2], v[ijk    -kk1], v[ijk    ], v[ijk    +kk1]) * interp2(w[ijk-jj1], w[ijk    ]),
interp4(v[ijk+jj1-kk2], v[ijk+jj1-kk1], v[ijk+jj1], v[ijk+jj1+kk1]) * interp2(w[ijk    ], w[ijk+jj1]),
interp4(v[ijk+jj2-kk2], v[ijk+jj2-kk1], v[ijk+jj2], v[ijk+jj2+kk1]) * interp2(w[ijk    ], w[ijk+jj3]), dyi)

- grad4x(interp4bot(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) * interp4bot(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]),
interp4       (w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) * interp4       (w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]),
interp4       (w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) * interp4       (w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]),
interp4       (w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3]) * interp4       (w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3]))
     * dzhi4[kstart+1];
     }

*/
    for (int k=grid->kstart+1; k<grid->kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                wt[ijk] +=
                         - grad4(interp4(u[ijk-ii1-kk2], u[ijk-ii1-kk1], u[ijk-ii1], u[ijk-ii1+kk1]) * interp2(w[ijk-ii3], w[ijk    ]),
                                 interp4(u[ijk    -kk2], u[ijk    -kk1], u[ijk    ], u[ijk    +kk1]) * interp2(w[ijk-ii1], w[ijk    ]),
                                 interp4(u[ijk+ii1-kk2], u[ijk+ii1-kk1], u[ijk+ii1], u[ijk+ii1+kk1]) * interp2(w[ijk    ], w[ijk+ii1]),
                                 interp4(u[ijk+ii2-kk2], u[ijk+ii2-kk1], u[ijk+ii2], u[ijk+ii2+kk1]) * interp2(w[ijk    ], w[ijk+ii3]), dxi)
            
                         - grad4(interp4(v[ijk-jj1-kk2], v[ijk-jj1-kk1], v[ijk-jj1], v[ijk-jj1+kk1]) * interp2(w[ijk-jj3], w[ijk    ]),
                                 interp4(v[ijk    -kk2], v[ijk    -kk1], v[ijk    ], v[ijk    +kk1]) * interp2(w[ijk-jj1], w[ijk    ]),
                                 interp4(v[ijk+jj1-kk2], v[ijk+jj1-kk1], v[ijk+jj1], v[ijk+jj1+kk1]) * interp2(w[ijk    ], w[ijk+jj1]),
                                 interp4(v[ijk+jj2-kk2], v[ijk+jj2-kk1], v[ijk+jj2], v[ijk+jj2+kk1]) * interp2(w[ijk    ], w[ijk+jj3]), dyi)
            
                         - grad4x(interp4(w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ]) * interp2(w[ijk-kk3], w[ijk    ]),
                                  interp4(w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) * interp2(w[ijk-kk1], w[ijk    ]),
                                  interp4(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) * interp2(w[ijk    ], w[ijk+kk1]),
                                  interp4(w[ijk    ], w[ijk+kk1], w[ijk+kk2], w[ijk+kk3]) * interp2(w[ijk    ], w[ijk+kk3]))
                           * dzhi4[k];
            }

/*
// top boundary
for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
for (int i=grid->istart; i<grid->iend; ++i)
{
ijk = i + j*jj1 + (kend-1)*kk1;
wt[ijk] +=
- grad4(interp4(u[ijk-ii1-kk2], u[ijk-ii1-kk1], u[ijk-ii1], u[ijk-ii1+kk1]) * interp2(w[ijk-ii3], w[ijk    ]),
interp4(u[ijk    -kk2], u[ijk    -kk1], u[ijk    ], u[ijk    +kk1]) * interp2(w[ijk-ii1], w[ijk    ]),
interp4(u[ijk+ii1-kk2], u[ijk+ii1-kk1], u[ijk+ii1], u[ijk+ii1+kk1]) * interp2(w[ijk    ], w[ijk+ii1]),
interp4(u[ijk+ii2-kk2], u[ijk+ii2-kk1], u[ijk+ii2], u[ijk+ii2+kk1]) * interp2(w[ijk    ], w[ijk+ii3]), dxi)

- grad4(interp4(v[ijk-jj1-kk2], v[ijk-jj1-kk1], v[ijk-jj1], v[ijk-jj1+kk1]) * interp2(w[ijk-jj3], w[ijk    ]),
interp4(v[ijk    -kk2], v[ijk    -kk1], v[ijk    ], v[ijk    +kk1]) * interp2(w[ijk-jj1], w[ijk    ]),
interp4(v[ijk+jj1-kk2], v[ijk+jj1-kk1], v[ijk+jj1], v[ijk+jj1+kk1]) * interp2(w[ijk    ], w[ijk+jj1]),
interp4(v[ijk+jj2-kk2], v[ijk+jj2-kk1], v[ijk+jj2], v[ijk+jj2+kk1]) * interp2(w[ijk    ], w[ijk+jj3]), dyi)

- grad4x(interp4       (w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ]) * interp4       (w[ijk-kk3], w[ijk-kk2], w[ijk-kk1], w[ijk    ]),
interp4       (w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]) * interp4       (w[ijk-kk2], w[ijk-kk1], w[ijk    ], w[ijk+kk1]),
interp4       (w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) * interp4       (w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]),
interp4top(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]) * interp4top(w[ijk-kk1], w[ijk    ], w[ijk+kk1], w[ijk+kk2]))
 * dzhi4[kend-1];
 }
 */
}

void Advec_4m::advec_s(double * restrict st, double * restrict s, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi4)
{
    const int ii1 = 1;
    const int ii2 = 2;
    const int ii3 = 3;
    const int jj1 = 1*grid->icells;
    const int jj2 = 2*grid->icells;
    const int jj3 = 3*grid->icells;
    const int kk1 = 1*grid->ijcells;
    const int kk2 = 2*grid->ijcells;
    const int kk3 = 3*grid->ijcells;

    const double dxi = 1./grid->dx;
    const double dyi = 1./grid->dy;

    const int kstart = grid->kstart;
    const int kend   = grid->kend;

    // bottom boundary
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + kstart*kk1;
            st[ijk] +=
                     - grad4(u[ijk-ii1] * interp2(s[ijk-ii3], s[ijk    ]),
                             u[ijk    ] * interp2(s[ijk-ii1], s[ijk    ]),
                             u[ijk+ii1] * interp2(s[ijk    ], s[ijk+ii1]),
                             u[ijk+ii2] * interp2(s[ijk    ], s[ijk+ii3]), dxi)

                     - grad4(v[ijk-jj1] * interp2(s[ijk-jj3], s[ijk    ]),
                             v[ijk    ] * interp2(s[ijk-jj1], s[ijk    ]),
                             v[ijk+jj1] * interp2(s[ijk    ], s[ijk+jj1]),
                             v[ijk+jj2] * interp2(s[ijk    ], s[ijk+jj3]), dyi)

                     - grad4x(-w[ijk+kk1] * interp2(s[ijk-kk1], s[ijk+kk2]),
                               w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]),
                               w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1]),
                               w[ijk+kk2] * interp2(s[ijk    ], s[ijk+kk3])) 
                       * dzi4[kstart];
        }

    for (int k=grid->kstart+1; k<grid->kend-1; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; ++i)
            {
                const int ijk = i + j*jj1 + k*kk1;
                st[ijk] +=
                         - grad4(u[ijk-ii1] * interp2(s[ijk-ii3], s[ijk    ]),
                                 u[ijk    ] * interp2(s[ijk-ii1], s[ijk    ]),
                                 u[ijk+ii1] * interp2(s[ijk    ], s[ijk+ii1]),
                                 u[ijk+ii2] * interp2(s[ijk    ], s[ijk+ii3]), dxi)

                         - grad4(v[ijk-jj1] * interp2(s[ijk-jj3], s[ijk    ]),
                                 v[ijk    ] * interp2(s[ijk-jj1], s[ijk    ]),
                                 v[ijk+jj1] * interp2(s[ijk    ], s[ijk+jj1]),
                                 v[ijk+jj2] * interp2(s[ijk    ], s[ijk+jj3]), dyi)

                         - grad4x(w[ijk-kk1] * interp2(s[ijk-kk3], s[ijk    ]),
                                  w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]),
                                  w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1]),
                                  w[ijk+kk2] * interp2(s[ijk    ], s[ijk+kk3])) 
                           * dzi4[k];
            }

    // top boundary
    for (int j=grid->jstart; j<grid->jend; ++j)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; ++i)
        {
            const int ijk = i + j*jj1 + (kend-1)*kk1;
            st[ijk] +=
                     - grad4(u[ijk-ii1] * interp2(s[ijk-ii3], s[ijk    ]),
                             u[ijk    ] * interp2(s[ijk-ii1], s[ijk    ]),
                             u[ijk+ii1] * interp2(s[ijk    ], s[ijk+ii1]),
                             u[ijk+ii2] * interp2(s[ijk    ], s[ijk+ii3]), dxi)

                     - grad4(v[ijk-jj1] * interp2(s[ijk-jj3], s[ijk    ]),
                             v[ijk    ] * interp2(s[ijk-jj1], s[ijk    ]),
                             v[ijk+jj1] * interp2(s[ijk    ], s[ijk+jj1]),
                             v[ijk+jj2] * interp2(s[ijk    ], s[ijk+jj3]), dyi)

                     - grad4x( w[ijk-kk1] * interp2(s[ijk-kk3], s[ijk    ]),
                               w[ijk    ] * interp2(s[ijk-kk1], s[ijk    ]),
                               w[ijk+kk1] * interp2(s[ijk    ], s[ijk+kk1]),
                              -w[ijk    ] * interp2(s[ijk-kk2], s[ijk+kk1])) 
                       * dzi4[kend-1];
        }
}
