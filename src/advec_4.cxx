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
#include "advec_4.h"
#include "defines.h"
#include "constants.h"
#include "finite_difference.h"
#include "model.h"

using namespace Finite_difference::O4;

Advec_4::Advec_4(Model* modelin, Input* inputin) : Advec(modelin, inputin)
{
    swadvec = "4";
}

Advec_4::~Advec_4()
{
}

#ifndef USECUDA
unsigned long Advec_4::get_time_limit(unsigned long idt, double dt)
{
    // Calculate cfl and prevent zero divisons.
    double cfl = calc_cfl(fields->u->data, fields->v->data, fields->w->data, grid->dzi, dt);
    cfl = std::max(cflmin, cfl);
    return idt * cflmax / cfl;
}

double Advec_4::get_cfl(double dt)
{
    return calc_cfl(fields->u->data, fields->v->data, fields->w->data, grid->dzi, dt);
}

void Advec_4::exec()
{
    // In case of a two-dimensional run, strip v component out of all kernels and do 
    // not calculate v-advection tendency.
    if (grid->jtot == 1)
    {
        advec_u<false>(fields->ut->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi4 );
        advec_w<false>(fields->wt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzhi4);

        for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
            advec_s<false>(it->second->data, fields->sp[it->first]->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi4);
    }
    else
    {
        advec_u<true>(fields->ut->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi4 );
        advec_v<true>(fields->vt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi4 );
        advec_w<true>(fields->wt->data, fields->u->data, fields->v->data, fields->w->data, grid->dzhi4);

        for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
            advec_s<true>(it->second->data, fields->sp[it->first]->data, fields->u->data, fields->v->data, fields->w->data, grid->dzi4);
    }
}
#endif

double Advec_4::calc_cfl(double * restrict u, double * restrict v, double * restrict w, double * restrict dzi, double dt)
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

    for (int k=grid->kstart; k<grid->kend; k++)
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj1 + k*kk1;
                cfl = std::max(cfl, std::abs(ci0*u[ijk-ii1] + ci1*u[ijk] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2])*dxi 
                                  + std::abs(ci0*v[ijk-jj1] + ci1*v[ijk] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2])*dyi 
                                  + std::abs(ci0*w[ijk-kk1] + ci1*w[ijk] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2])*dzi[k] );
            }

    grid->get_max(&cfl);

    cfl = cfl*dt;

    return cfl;
}

    template<bool dim3>
void Advec_4::advec_u(double * restrict ut, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi4)
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
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ijk = i + j*jj1 + kstart*kk1;
            ut[ijk] -= ( cg0*((ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]) * (ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]))
                       + cg1*((ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]) * (ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]))
                       + cg2*((ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]) * (ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]))
                       + cg3*((ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3]) * (ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3])) ) * cgi*dxi;

            if (dim3)
            {
                ut[ijk] -= ( cg0*((ci0*v[ijk-ii2-jj1] + ci1*v[ijk-ii1-jj1] + ci2*v[ijk-jj1] + ci3*v[ijk+ii1-jj1]) * (ci0*u[ijk-jj3] + ci1*u[ijk-jj2] + ci2*u[ijk-jj1] + ci3*u[ijk    ]))
                           + cg1*((ci0*v[ijk-ii2    ] + ci1*v[ijk-ii1    ] + ci2*v[ijk    ] + ci3*v[ijk+ii1    ]) * (ci0*u[ijk-jj2] + ci1*u[ijk-jj1] + ci2*u[ijk    ] + ci3*u[ijk+jj1]))
                           + cg2*((ci0*v[ijk-ii2+jj1] + ci1*v[ijk-ii1+jj1] + ci2*v[ijk+jj1] + ci3*v[ijk+ii1+jj1]) * (ci0*u[ijk-jj1] + ci1*u[ijk    ] + ci2*u[ijk+jj1] + ci3*u[ijk+jj2]))
                           + cg3*((ci0*v[ijk-ii2+jj2] + ci1*v[ijk-ii1+jj2] + ci2*v[ijk+jj2] + ci3*v[ijk+ii1+jj2]) * (ci0*u[ijk    ] + ci1*u[ijk+jj1] + ci2*u[ijk+jj2] + ci3*u[ijk+jj3])) ) * cgi*dyi;
            }

            ut[ijk] -= ( cg0*((ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk-kk1] + ci3*w[ijk+ii1-kk1]) * (bi0*u[ijk-kk2] + bi1*u[ijk-kk1] + bi2*u[ijk    ] + bi3*u[ijk+kk1]))
                       + cg1*((ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk    ] + ci3*w[ijk+ii1    ]) * (ci0*u[ijk-kk2] + ci1*u[ijk-kk1] + ci2*u[ijk    ] + ci3*u[ijk+kk1]))
                       + cg2*((ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk+kk1] + ci3*w[ijk+ii1+kk1]) * (ci0*u[ijk-kk1] + ci1*u[ijk    ] + ci2*u[ijk+kk1] + ci3*u[ijk+kk2]))
                       + cg3*((ci0*w[ijk-ii2+kk2] + ci1*w[ijk-ii1+kk2] + ci2*w[ijk+kk2] + ci3*w[ijk+ii1+kk2]) * (ci0*u[ijk    ] + ci1*u[ijk+kk1] + ci2*u[ijk+kk2] + ci3*u[ijk+kk3])) )
                     * dzi4[kstart];
        }

    for (int k=grid->kstart+1; k<grid->kend-1; k++)
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj1 + k*kk1;
                ut[ijk] -= ( cg0*((ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]) * (ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]))
                           + cg1*((ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]) * (ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]))
                           + cg2*((ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]) * (ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]))
                           + cg3*((ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3]) * (ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3])) ) * cgi*dxi;

                if (dim3)
                {
                    ut[ijk] -= ( cg0*((ci0*v[ijk-ii2-jj1] + ci1*v[ijk-ii1-jj1] + ci2*v[ijk-jj1] + ci3*v[ijk+ii1-jj1]) * (ci0*u[ijk-jj3] + ci1*u[ijk-jj2] + ci2*u[ijk-jj1] + ci3*u[ijk    ]))
                               + cg1*((ci0*v[ijk-ii2    ] + ci1*v[ijk-ii1    ] + ci2*v[ijk    ] + ci3*v[ijk+ii1    ]) * (ci0*u[ijk-jj2] + ci1*u[ijk-jj1] + ci2*u[ijk    ] + ci3*u[ijk+jj1]))
                               + cg2*((ci0*v[ijk-ii2+jj1] + ci1*v[ijk-ii1+jj1] + ci2*v[ijk+jj1] + ci3*v[ijk+ii1+jj1]) * (ci0*u[ijk-jj1] + ci1*u[ijk    ] + ci2*u[ijk+jj1] + ci3*u[ijk+jj2]))
                               + cg3*((ci0*v[ijk-ii2+jj2] + ci1*v[ijk-ii1+jj2] + ci2*v[ijk+jj2] + ci3*v[ijk+ii1+jj2]) * (ci0*u[ijk    ] + ci1*u[ijk+jj1] + ci2*u[ijk+jj2] + ci3*u[ijk+jj3])) ) * cgi*dyi;
                }

                ut[ijk] -= ( cg0*((ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk-kk1] + ci3*w[ijk+ii1-kk1]) * (ci0*u[ijk-kk3] + ci1*u[ijk-kk2] + ci2*u[ijk-kk1] + ci3*u[ijk    ]))
                           + cg1*((ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk    ] + ci3*w[ijk+ii1    ]) * (ci0*u[ijk-kk2] + ci1*u[ijk-kk1] + ci2*u[ijk    ] + ci3*u[ijk+kk1]))
                           + cg2*((ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk+kk1] + ci3*w[ijk+ii1+kk1]) * (ci0*u[ijk-kk1] + ci1*u[ijk    ] + ci2*u[ijk+kk1] + ci3*u[ijk+kk2]))
                           + cg3*((ci0*w[ijk-ii2+kk2] + ci1*w[ijk-ii1+kk2] + ci2*w[ijk+kk2] + ci3*w[ijk+ii1+kk2]) * (ci0*u[ijk    ] + ci1*u[ijk+kk1] + ci2*u[ijk+kk2] + ci3*u[ijk+kk3])) )
                         * dzi4[k];
            }

    // top boundary
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ijk = i + j*jj1 + (kend-1)*kk1;
            ut[ijk] -= ( cg0*((ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]) * (ci0*u[ijk-ii3] + ci1*u[ijk-ii2] + ci2*u[ijk-ii1] + ci3*u[ijk    ]))
                       + cg1*((ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]) * (ci0*u[ijk-ii2] + ci1*u[ijk-ii1] + ci2*u[ijk    ] + ci3*u[ijk+ii1]))
                       + cg2*((ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]) * (ci0*u[ijk-ii1] + ci1*u[ijk    ] + ci2*u[ijk+ii1] + ci3*u[ijk+ii2]))
                       + cg3*((ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3]) * (ci0*u[ijk    ] + ci1*u[ijk+ii1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii3])) ) * cgi*dxi;

            if (dim3)
            {
                ut[ijk] -= ( cg0*((ci0*v[ijk-ii2-jj1] + ci1*v[ijk-ii1-jj1] + ci2*v[ijk-jj1] + ci3*v[ijk+ii1-jj1]) * (ci0*u[ijk-jj3] + ci1*u[ijk-jj2] + ci2*u[ijk-jj1] + ci3*u[ijk    ]))
                           + cg1*((ci0*v[ijk-ii2    ] + ci1*v[ijk-ii1    ] + ci2*v[ijk    ] + ci3*v[ijk+ii1    ]) * (ci0*u[ijk-jj2] + ci1*u[ijk-jj1] + ci2*u[ijk    ] + ci3*u[ijk+jj1]))
                           + cg2*((ci0*v[ijk-ii2+jj1] + ci1*v[ijk-ii1+jj1] + ci2*v[ijk+jj1] + ci3*v[ijk+ii1+jj1]) * (ci0*u[ijk-jj1] + ci1*u[ijk    ] + ci2*u[ijk+jj1] + ci3*u[ijk+jj2]))
                           + cg3*((ci0*v[ijk-ii2+jj2] + ci1*v[ijk-ii1+jj2] + ci2*v[ijk+jj2] + ci3*v[ijk+ii1+jj2]) * (ci0*u[ijk    ] + ci1*u[ijk+jj1] + ci2*u[ijk+jj2] + ci3*u[ijk+jj3])) ) * cgi*dyi;
            }

            ut[ijk] -= ( cg0*((ci0*w[ijk-ii2-kk1] + ci1*w[ijk-ii1-kk1] + ci2*w[ijk-kk1] + ci3*w[ijk+ii1-kk1]) * (ci0*u[ijk-kk3] + ci1*u[ijk-kk2] + ci2*u[ijk-kk1] + ci3*u[ijk    ]))
                       + cg1*((ci0*w[ijk-ii2    ] + ci1*w[ijk-ii1    ] + ci2*w[ijk    ] + ci3*w[ijk+ii1    ]) * (ci0*u[ijk-kk2] + ci1*u[ijk-kk1] + ci2*u[ijk    ] + ci3*u[ijk+kk1]))
                       + cg2*((ci0*w[ijk-ii2+kk1] + ci1*w[ijk-ii1+kk1] + ci2*w[ijk+kk1] + ci3*w[ijk+ii1+kk1]) * (ci0*u[ijk-kk1] + ci1*u[ijk    ] + ci2*u[ijk+kk1] + ci3*u[ijk+kk2]))
                       + cg3*((ci0*w[ijk-ii2+kk2] + ci1*w[ijk-ii1+kk2] + ci2*w[ijk+kk2] + ci3*w[ijk+ii1+kk2]) * (ti0*u[ijk-kk1] + ti1*u[ijk    ] + ti2*u[ijk+kk1] + ti3*u[ijk+kk2])) )
                     * dzi4[kend-1];
        }
}

    template<bool dim3>
void Advec_4::advec_v(double * restrict vt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi4)
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
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ijk = i + j*jj1 + kstart*kk1;
            vt[ijk] -= ( cg0*((ci0*u[ijk-ii1-jj2] + ci1*u[ijk-ii1-jj1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+jj1]) * (ci0*v[ijk-ii3] + ci1*v[ijk-ii2] + ci2*v[ijk-ii1] + ci3*v[ijk    ]))
                       + cg1*((ci0*u[ijk    -jj2] + ci1*u[ijk    -jj1] + ci2*u[ijk    ] + ci3*u[ijk    +jj1]) * (ci0*v[ijk-ii2] + ci1*v[ijk-ii1] + ci2*v[ijk    ] + ci3*v[ijk+ii1]))
                       + cg2*((ci0*u[ijk+ii1-jj2] + ci1*u[ijk+ii1-jj1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+jj1]) * (ci0*v[ijk-ii1] + ci1*v[ijk    ] + ci2*v[ijk+ii1] + ci3*v[ijk+ii2]))
                       + cg3*((ci0*u[ijk+ii2-jj2] + ci1*u[ijk+ii2-jj1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+jj1]) * (ci0*v[ijk    ] + ci1*v[ijk+ii1] + ci2*v[ijk+ii2] + ci3*v[ijk+ii3])) ) * cgi*dxi;

            if (dim3)
            {
                vt[ijk] -= ( cg0*((ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]) * (ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]))
                           + cg1*((ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]) * (ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]))
                           + cg2*((ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]) * (ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]))
                           + cg3*((ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3]) * (ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3])) ) * cgi*dyi;
            }

            vt[ijk] -= ( cg0*((ci0*w[ijk-jj2-kk1] + ci1*w[ijk-jj1-kk1] + ci2*w[ijk-kk1] + ci3*w[ijk+jj1-kk1]) * (bi0*v[ijk-kk2] + bi1*v[ijk-kk1] + bi2*v[ijk    ] + bi3*v[ijk+kk1]))
                       + cg1*((ci0*w[ijk-jj2    ] + ci1*w[ijk-jj1    ] + ci2*w[ijk    ] + ci3*w[ijk+jj1    ]) * (ci0*v[ijk-kk2] + ci1*v[ijk-kk1] + ci2*v[ijk    ] + ci3*v[ijk+kk1]))
                       + cg2*((ci0*w[ijk-jj2+kk1] + ci1*w[ijk-jj1+kk1] + ci2*w[ijk+kk1] + ci3*w[ijk+jj1+kk1]) * (ci0*v[ijk-kk1] + ci1*v[ijk    ] + ci2*v[ijk+kk1] + ci3*v[ijk+kk2]))
                       + cg3*((ci0*w[ijk-jj2+kk2] + ci1*w[ijk-jj1+kk2] + ci2*w[ijk+kk2] + ci3*w[ijk+jj1+kk2]) * (ci0*v[ijk    ] + ci1*v[ijk+kk1] + ci2*v[ijk+kk2] + ci3*v[ijk+kk3])) )
                     * dzi4[kstart];
        }

    for (int k=grid->kstart+1; k<grid->kend-1; k++)
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj1 + k*kk1;
                vt[ijk] -= ( cg0*((ci0*u[ijk-ii1-jj2] + ci1*u[ijk-ii1-jj1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+jj1]) * (ci0*v[ijk-ii3] + ci1*v[ijk-ii2] + ci2*v[ijk-ii1] + ci3*v[ijk    ]))
                           + cg1*((ci0*u[ijk    -jj2] + ci1*u[ijk    -jj1] + ci2*u[ijk    ] + ci3*u[ijk    +jj1]) * (ci0*v[ijk-ii2] + ci1*v[ijk-ii1] + ci2*v[ijk    ] + ci3*v[ijk+ii1]))
                           + cg2*((ci0*u[ijk+ii1-jj2] + ci1*u[ijk+ii1-jj1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+jj1]) * (ci0*v[ijk-ii1] + ci1*v[ijk    ] + ci2*v[ijk+ii1] + ci3*v[ijk+ii2]))
                           + cg3*((ci0*u[ijk+ii2-jj2] + ci1*u[ijk+ii2-jj1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+jj1]) * (ci0*v[ijk    ] + ci1*v[ijk+ii1] + ci2*v[ijk+ii2] + ci3*v[ijk+ii3])) ) * cgi*dxi;

                if (dim3)
                {
                    vt[ijk] -= ( cg0*((ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]) * (ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]))
                               + cg1*((ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]) * (ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]))
                               + cg2*((ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]) * (ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]))
                               + cg3*((ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3]) * (ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3])) ) * cgi*dyi;
                }

                vt[ijk] -= ( cg0*((ci0*w[ijk-jj2-kk1] + ci1*w[ijk-jj1-kk1] + ci2*w[ijk-kk1] + ci3*w[ijk+jj1-kk1]) * (ci0*v[ijk-kk3] + ci1*v[ijk-kk2] + ci2*v[ijk-kk1] + ci3*v[ijk    ]))
                           + cg1*((ci0*w[ijk-jj2    ] + ci1*w[ijk-jj1    ] + ci2*w[ijk    ] + ci3*w[ijk+jj1    ]) * (ci0*v[ijk-kk2] + ci1*v[ijk-kk1] + ci2*v[ijk    ] + ci3*v[ijk+kk1]))
                           + cg2*((ci0*w[ijk-jj2+kk1] + ci1*w[ijk-jj1+kk1] + ci2*w[ijk+kk1] + ci3*w[ijk+jj1+kk1]) * (ci0*v[ijk-kk1] + ci1*v[ijk    ] + ci2*v[ijk+kk1] + ci3*v[ijk+kk2]))
                           + cg3*((ci0*w[ijk-jj2+kk2] + ci1*w[ijk-jj1+kk2] + ci2*w[ijk+kk2] + ci3*w[ijk+jj1+kk2]) * (ci0*v[ijk    ] + ci1*v[ijk+kk1] + ci2*v[ijk+kk2] + ci3*v[ijk+kk3])) )
                         * dzi4[k];
            }

    // top boundary
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ijk = i + j*jj1 + (kend-1)*kk1;
            vt[ijk] -= ( cg0*((ci0*u[ijk-ii1-jj2] + ci1*u[ijk-ii1-jj1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+jj1]) * (ci0*v[ijk-ii3] + ci1*v[ijk-ii2] + ci2*v[ijk-ii1] + ci3*v[ijk    ]))
                       + cg1*((ci0*u[ijk    -jj2] + ci1*u[ijk    -jj1] + ci2*u[ijk    ] + ci3*u[ijk    +jj1]) * (ci0*v[ijk-ii2] + ci1*v[ijk-ii1] + ci2*v[ijk    ] + ci3*v[ijk+ii1]))
                       + cg2*((ci0*u[ijk+ii1-jj2] + ci1*u[ijk+ii1-jj1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+jj1]) * (ci0*v[ijk-ii1] + ci1*v[ijk    ] + ci2*v[ijk+ii1] + ci3*v[ijk+ii2]))
                       + cg3*((ci0*u[ijk+ii2-jj2] + ci1*u[ijk+ii2-jj1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+jj1]) * (ci0*v[ijk    ] + ci1*v[ijk+ii1] + ci2*v[ijk+ii2] + ci3*v[ijk+ii3])) ) * cgi*dxi;

            if (dim3)
            {
                vt[ijk] -= ( cg0*((ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]) * (ci0*v[ijk-jj3] + ci1*v[ijk-jj2] + ci2*v[ijk-jj1] + ci3*v[ijk    ]))
                           + cg1*((ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]) * (ci0*v[ijk-jj2] + ci1*v[ijk-jj1] + ci2*v[ijk    ] + ci3*v[ijk+jj1]))
                           + cg2*((ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]) * (ci0*v[ijk-jj1] + ci1*v[ijk    ] + ci2*v[ijk+jj1] + ci3*v[ijk+jj2]))
                           + cg3*((ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3]) * (ci0*v[ijk    ] + ci1*v[ijk+jj1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj3])) ) * cgi*dyi;
            }

            vt[ijk] -= ( cg0*((ci0*w[ijk-jj2-kk1] + ci1*w[ijk-jj1-kk1] + ci2*w[ijk-kk1] + ci3*w[ijk+jj1-kk1]) * (ci0*v[ijk-kk3] + ci1*v[ijk-kk2] + ci2*v[ijk-kk1] + ci3*v[ijk    ]))
                       + cg1*((ci0*w[ijk-jj2    ] + ci1*w[ijk-jj1    ] + ci2*w[ijk    ] + ci3*w[ijk+jj1    ]) * (ci0*v[ijk-kk2] + ci1*v[ijk-kk1] + ci2*v[ijk    ] + ci3*v[ijk+kk1]))
                       + cg2*((ci0*w[ijk-jj2+kk1] + ci1*w[ijk-jj1+kk1] + ci2*w[ijk+kk1] + ci3*w[ijk+jj1+kk1]) * (ci0*v[ijk-kk1] + ci1*v[ijk    ] + ci2*v[ijk+kk1] + ci3*v[ijk+kk2]))
                       + cg3*((ci0*w[ijk-jj2+kk2] + ci1*w[ijk-jj1+kk2] + ci2*w[ijk+kk2] + ci3*w[ijk+jj1+kk2]) * (ti0*v[ijk-kk1] + ti1*v[ijk    ] + ti2*v[ijk+kk1] + ti3*v[ijk+kk2])) )
                     * dzi4[kend-1];
        }
}

    template<bool dim3>
void Advec_4::advec_w(double * restrict wt, double * restrict u, double * restrict v, double * restrict w, double * restrict dzhi4)
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
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ijk = i + j*jj1 + (kstart+1)*kk1;
            wt[ijk] -= ( cg0*((ci0*u[ijk-ii1-kk2] + ci1*u[ijk-ii1-kk1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+kk1]) * (ci0*w[ijk-ii3] + ci1*w[ijk-ii2] + ci2*w[ijk-ii1] + ci3*w[ijk    ]))
                    + cg1*((ci0*u[ijk    -kk2] + ci1*u[ijk    -kk1] + ci2*u[ijk    ] + ci3*u[ijk    +kk1]) * (ci0*w[ijk-ii2] + ci1*w[ijk-ii1] + ci2*w[ijk    ] + ci3*w[ijk+ii1]))
                    + cg2*((ci0*u[ijk+ii1-kk2] + ci1*u[ijk+ii1-kk1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+kk1]) * (ci0*w[ijk-ii1] + ci1*w[ijk    ] + ci2*w[ijk+ii1] + ci3*w[ijk+ii2]))
                    + cg3*((ci0*u[ijk+ii2-kk2] + ci1*u[ijk+ii2-kk1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+ii1] + ci2*w[ijk+ii2] + ci3*w[ijk+ii3])) ) * cgi*dxi;

            if (dim3)
            {
                wt[ijk] -= ( cg0*((ci0*v[ijk-jj1-kk2] + ci1*v[ijk-jj1-kk1] + ci2*v[ijk-jj1] + ci3*v[ijk-jj1+kk1]) * (ci0*w[ijk-jj3] + ci1*w[ijk-jj2] + ci2*w[ijk-jj1] + ci3*w[ijk    ]))
                        + cg1*((ci0*v[ijk    -kk2] + ci1*v[ijk    -kk1] + ci2*v[ijk    ] + ci3*v[ijk    +kk1]) * (ci0*w[ijk-jj2] + ci1*w[ijk-jj1] + ci2*w[ijk    ] + ci3*w[ijk+jj1]))
                        + cg2*((ci0*v[ijk+jj1-kk2] + ci1*v[ijk+jj1-kk1] + ci2*v[ijk+jj1] + ci3*v[ijk+jj1+kk1]) * (ci0*w[ijk-jj1] + ci1*w[ijk    ] + ci2*w[ijk+jj1] + ci3*w[ijk+jj2]))
                        + cg3*((ci0*v[ijk+jj2-kk2] + ci1*v[ijk+jj2-kk1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+jj1] + ci2*w[ijk+jj2] + ci3*w[ijk+jj3])) ) * cgi*dyi;
            }

            wt[ijk] -= ( cg0*((bi0*w[ijk-kk2] + bi1*w[ijk-kk1] + bi2*w[ijk    ] + bi3*w[ijk+kk1]) * (bi0*w[ijk-kk2] + bi1*w[ijk-kk1] + bi2*w[ijk    ] + bi3*w[ijk+kk1]))
                    + cg1*((ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1]) * (ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1]))
                    + cg2*((ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2]) * (ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2]))
                    + cg3*((ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3]) * (ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3])) )
                * dzhi4[kstart+1];
        }

    for (int k=grid->kstart+2; k<grid->kend-1; k++)
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj1 + k*kk1;
                wt[ijk] -= ( cg0*((ci0*u[ijk-ii1-kk2] + ci1*u[ijk-ii1-kk1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+kk1]) * (ci0*w[ijk-ii3] + ci1*w[ijk-ii2] + ci2*w[ijk-ii1] + ci3*w[ijk    ]))
                        + cg1*((ci0*u[ijk    -kk2] + ci1*u[ijk    -kk1] + ci2*u[ijk    ] + ci3*u[ijk    +kk1]) * (ci0*w[ijk-ii2] + ci1*w[ijk-ii1] + ci2*w[ijk    ] + ci3*w[ijk+ii1]))
                        + cg2*((ci0*u[ijk+ii1-kk2] + ci1*u[ijk+ii1-kk1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+kk1]) * (ci0*w[ijk-ii1] + ci1*w[ijk    ] + ci2*w[ijk+ii1] + ci3*w[ijk+ii2]))
                        + cg3*((ci0*u[ijk+ii2-kk2] + ci1*u[ijk+ii2-kk1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+ii1] + ci2*w[ijk+ii2] + ci3*w[ijk+ii3])) ) * cgi*dxi;

                if (dim3)
                {
                    wt[ijk] -= ( cg0*((ci0*v[ijk-jj1-kk2] + ci1*v[ijk-jj1-kk1] + ci2*v[ijk-jj1] + ci3*v[ijk-jj1+kk1]) * (ci0*w[ijk-jj3] + ci1*w[ijk-jj2] + ci2*w[ijk-jj1] + ci3*w[ijk    ]))
                            + cg1*((ci0*v[ijk    -kk2] + ci1*v[ijk    -kk1] + ci2*v[ijk    ] + ci3*v[ijk    +kk1]) * (ci0*w[ijk-jj2] + ci1*w[ijk-jj1] + ci2*w[ijk    ] + ci3*w[ijk+jj1]))
                            + cg2*((ci0*v[ijk+jj1-kk2] + ci1*v[ijk+jj1-kk1] + ci2*v[ijk+jj1] + ci3*v[ijk+jj1+kk1]) * (ci0*w[ijk-jj1] + ci1*w[ijk    ] + ci2*w[ijk+jj1] + ci3*w[ijk+jj2]))
                            + cg3*((ci0*v[ijk+jj2-kk2] + ci1*v[ijk+jj2-kk1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+jj1] + ci2*w[ijk+jj2] + ci3*w[ijk+jj3])) ) * cgi*dyi;
                }

                wt[ijk] -= ( cg0*((ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ]) * (ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ]))
                        + cg1*((ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1]) * (ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1]))
                        + cg2*((ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2]) * (ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2]))
                        + cg3*((ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3]) * (ci0*w[ijk    ] + ci1*w[ijk+kk1] + ci2*w[ijk+kk2] + ci3*w[ijk+kk3])) )
                    * dzhi4[k];
            }

    // top boundary
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ijk = i + j*jj1 + (kend-1)*kk1;
            wt[ijk] -= ( cg0*((ci0*u[ijk-ii1-kk2] + ci1*u[ijk-ii1-kk1] + ci2*u[ijk-ii1] + ci3*u[ijk-ii1+kk1]) * (ci0*w[ijk-ii3] + ci1*w[ijk-ii2] + ci2*w[ijk-ii1] + ci3*w[ijk    ]))
                    + cg1*((ci0*u[ijk    -kk2] + ci1*u[ijk    -kk1] + ci2*u[ijk    ] + ci3*u[ijk    +kk1]) * (ci0*w[ijk-ii2] + ci1*w[ijk-ii1] + ci2*w[ijk    ] + ci3*w[ijk+ii1]))
                    + cg2*((ci0*u[ijk+ii1-kk2] + ci1*u[ijk+ii1-kk1] + ci2*u[ijk+ii1] + ci3*u[ijk+ii1+kk1]) * (ci0*w[ijk-ii1] + ci1*w[ijk    ] + ci2*w[ijk+ii1] + ci3*w[ijk+ii2]))
                    + cg3*((ci0*u[ijk+ii2-kk2] + ci1*u[ijk+ii2-kk1] + ci2*u[ijk+ii2] + ci3*u[ijk+ii2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+ii1] + ci2*w[ijk+ii2] + ci3*w[ijk+ii3])) ) * cgi*dxi;

            if (dim3)
            {
                wt[ijk] -= ( cg0*((ci0*v[ijk-jj1-kk2] + ci1*v[ijk-jj1-kk1] + ci2*v[ijk-jj1] + ci3*v[ijk-jj1+kk1]) * (ci0*w[ijk-jj3] + ci1*w[ijk-jj2] + ci2*w[ijk-jj1] + ci3*w[ijk    ]))
                        + cg1*((ci0*v[ijk    -kk2] + ci1*v[ijk    -kk1] + ci2*v[ijk    ] + ci3*v[ijk    +kk1]) * (ci0*w[ijk-jj2] + ci1*w[ijk-jj1] + ci2*w[ijk    ] + ci3*w[ijk+jj1]))
                        + cg2*((ci0*v[ijk+jj1-kk2] + ci1*v[ijk+jj1-kk1] + ci2*v[ijk+jj1] + ci3*v[ijk+jj1+kk1]) * (ci0*w[ijk-jj1] + ci1*w[ijk    ] + ci2*w[ijk+jj1] + ci3*w[ijk+jj2]))
                        + cg3*((ci0*v[ijk+jj2-kk2] + ci1*v[ijk+jj2-kk1] + ci2*v[ijk+jj2] + ci3*v[ijk+jj2+kk1]) * (ci0*w[ijk    ] + ci1*w[ijk+jj1] + ci2*w[ijk+jj2] + ci3*w[ijk+jj3])) ) * cgi*dyi;
            }

            wt[ijk] -= ( cg0*((ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ]) * (ci0*w[ijk-kk3] + ci1*w[ijk-kk2] + ci2*w[ijk-kk1] + ci3*w[ijk    ]))
                    + cg1*((ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1]) * (ci0*w[ijk-kk2] + ci1*w[ijk-kk1] + ci2*w[ijk    ] + ci3*w[ijk+kk1]))
                    + cg2*((ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2]) * (ci0*w[ijk-kk1] + ci1*w[ijk    ] + ci2*w[ijk+kk1] + ci3*w[ijk+kk2]))
                    + cg3*((ti0*w[ijk-kk1] + ti1*w[ijk    ] + ti2*w[ijk+kk1] + ti3*w[ijk+kk2]) * (ti0*w[ijk-kk1] + ti1*w[ijk    ] + ti2*w[ijk+kk1] + ti3*w[ijk+kk2])) )
                * dzhi4[kend-1];
        }
}

    template<bool dim3>
void Advec_4::advec_s(double * restrict st, double * restrict s, double * restrict u, double * restrict v, double * restrict w, double * restrict dzi4)
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
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ijk = i + j*jj1 + kstart*kk1;
            st[ijk] -= ( cg0*(u[ijk-ii1] * (ci0*s[ijk-ii3] + ci1*s[ijk-ii2] + ci2*s[ijk-ii1] + ci3*s[ijk    ]))
                       + cg1*(u[ijk    ] * (ci0*s[ijk-ii2] + ci1*s[ijk-ii1] + ci2*s[ijk    ] + ci3*s[ijk+ii1]))
                       + cg2*(u[ijk+ii1] * (ci0*s[ijk-ii1] + ci1*s[ijk    ] + ci2*s[ijk+ii1] + ci3*s[ijk+ii2]))
                       + cg3*(u[ijk+ii2] * (ci0*s[ijk    ] + ci1*s[ijk+ii1] + ci2*s[ijk+ii2] + ci3*s[ijk+ii3])) ) * cgi*dxi;

            if (dim3)
            {
                st[ijk] -= ( cg0*(v[ijk-jj1] * (ci0*s[ijk-jj3] + ci1*s[ijk-jj2] + ci2*s[ijk-jj1] + ci3*s[ijk    ]))
                           + cg1*(v[ijk    ] * (ci0*s[ijk-jj2] + ci1*s[ijk-jj1] + ci2*s[ijk    ] + ci3*s[ijk+jj1]))
                           + cg2*(v[ijk+jj1] * (ci0*s[ijk-jj1] + ci1*s[ijk    ] + ci2*s[ijk+jj1] + ci3*s[ijk+jj2]))
                           + cg3*(v[ijk+jj2] * (ci0*s[ijk    ] + ci1*s[ijk+jj1] + ci2*s[ijk+jj2] + ci3*s[ijk+jj3])) ) * cgi*dyi;
            }

            st[ijk] -= ( cg0*(w[ijk-kk1] * (bi0*s[ijk-kk2] + bi1*s[ijk-kk1] + bi2*s[ijk    ] + bi3*s[ijk+kk1]))
                       + cg1*(w[ijk    ] * (ci0*s[ijk-kk2] + ci1*s[ijk-kk1] + ci2*s[ijk    ] + ci3*s[ijk+kk1]))
                       + cg2*(w[ijk+kk1] * (ci0*s[ijk-kk1] + ci1*s[ijk    ] + ci2*s[ijk+kk1] + ci3*s[ijk+kk2]))
                       + cg3*(w[ijk+kk2] * (ci0*s[ijk    ] + ci1*s[ijk+kk1] + ci2*s[ijk+kk2] + ci3*s[ijk+kk3])) )
                     * dzi4[kstart];
        }

    for (int k=grid->kstart+1; k<grid->kend-1; k++)
        for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
            for (int i=grid->istart; i<grid->iend; i++)
            {
                const int ijk = i + j*jj1 + k*kk1;
                st[ijk] -= ( cg0*(u[ijk-ii1] * (ci0*s[ijk-ii3] + ci1*s[ijk-ii2] + ci2*s[ijk-ii1] + ci3*s[ijk    ]))
                           + cg1*(u[ijk    ] * (ci0*s[ijk-ii2] + ci1*s[ijk-ii1] + ci2*s[ijk    ] + ci3*s[ijk+ii1]))
                           + cg2*(u[ijk+ii1] * (ci0*s[ijk-ii1] + ci1*s[ijk    ] + ci2*s[ijk+ii1] + ci3*s[ijk+ii2]))
                           + cg3*(u[ijk+ii2] * (ci0*s[ijk    ] + ci1*s[ijk+ii1] + ci2*s[ijk+ii2] + ci3*s[ijk+ii3])) ) * cgi*dxi;

                if (dim3)
                {
                    st[ijk] -= ( cg0*(v[ijk-jj1] * (ci0*s[ijk-jj3] + ci1*s[ijk-jj2] + ci2*s[ijk-jj1] + ci3*s[ijk    ]))
                               + cg1*(v[ijk    ] * (ci0*s[ijk-jj2] + ci1*s[ijk-jj1] + ci2*s[ijk    ] + ci3*s[ijk+jj1]))
                               + cg2*(v[ijk+jj1] * (ci0*s[ijk-jj1] + ci1*s[ijk    ] + ci2*s[ijk+jj1] + ci3*s[ijk+jj2]))
                               + cg3*(v[ijk+jj2] * (ci0*s[ijk    ] + ci1*s[ijk+jj1] + ci2*s[ijk+jj2] + ci3*s[ijk+jj3])) ) * cgi*dyi;
                }

                st[ijk] -= ( cg0*(w[ijk-kk1] * (ci0*s[ijk-kk3] + ci1*s[ijk-kk2] + ci2*s[ijk-kk1] + ci3*s[ijk    ]))
                           + cg1*(w[ijk    ] * (ci0*s[ijk-kk2] + ci1*s[ijk-kk1] + ci2*s[ijk    ] + ci3*s[ijk+kk1]))
                           + cg2*(w[ijk+kk1] * (ci0*s[ijk-kk1] + ci1*s[ijk    ] + ci2*s[ijk+kk1] + ci3*s[ijk+kk2]))
                           + cg3*(w[ijk+kk2] * (ci0*s[ijk    ] + ci1*s[ijk+kk1] + ci2*s[ijk+kk2] + ci3*s[ijk+kk3])) )
                         * dzi4[k];
            }

    // top boundary
    for (int j=grid->jstart; j<grid->jend; j++)
#pragma ivdep
        for (int i=grid->istart; i<grid->iend; i++)
        {
            const int ijk = i + j*jj1 + (kend-1)*kk1;
            st[ijk] -= ( cg0*(u[ijk-ii1] * (ci0*s[ijk-ii3] + ci1*s[ijk-ii2] + ci2*s[ijk-ii1] + ci3*s[ijk    ]))
                       + cg1*(u[ijk    ] * (ci0*s[ijk-ii2] + ci1*s[ijk-ii1] + ci2*s[ijk    ] + ci3*s[ijk+ii1]))
                       + cg2*(u[ijk+ii1] * (ci0*s[ijk-ii1] + ci1*s[ijk    ] + ci2*s[ijk+ii1] + ci3*s[ijk+ii2]))
                       + cg3*(u[ijk+ii2] * (ci0*s[ijk    ] + ci1*s[ijk+ii1] + ci2*s[ijk+ii2] + ci3*s[ijk+ii3])) ) * cgi*dxi;

            if (dim3)
            {
                st[ijk] -= ( cg0*(v[ijk-jj1] * (ci0*s[ijk-jj3] + ci1*s[ijk-jj2] + ci2*s[ijk-jj1] + ci3*s[ijk    ]))
                           + cg1*(v[ijk    ] * (ci0*s[ijk-jj2] + ci1*s[ijk-jj1] + ci2*s[ijk    ] + ci3*s[ijk+jj1]))
                           + cg2*(v[ijk+jj1] * (ci0*s[ijk-jj1] + ci1*s[ijk    ] + ci2*s[ijk+jj1] + ci3*s[ijk+jj2]))
                           + cg3*(v[ijk+jj2] * (ci0*s[ijk    ] + ci1*s[ijk+jj1] + ci2*s[ijk+jj2] + ci3*s[ijk+jj3])) ) * cgi*dyi;
            }

            st[ijk] -= ( cg0*(w[ijk-kk1] * (ci0*s[ijk-kk3] + ci1*s[ijk-kk2] + ci2*s[ijk-kk1] + ci3*s[ijk    ]))
                       + cg1*(w[ijk    ] * (ci0*s[ijk-kk2] + ci1*s[ijk-kk1] + ci2*s[ijk    ] + ci3*s[ijk+kk1]))
                       + cg2*(w[ijk+kk1] * (ci0*s[ijk-kk1] + ci1*s[ijk    ] + ci2*s[ijk+kk1] + ci3*s[ijk+kk2]))
                       + cg3*(w[ijk+kk2] * (ti0*s[ijk-kk1] + ti1*s[ijk    ] + ti2*s[ijk+kk1] + ti3*s[ijk+kk2])) )
                     * dzi4[kend-1];
        }
}
