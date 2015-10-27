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
#include <iostream>
#include "master.h"
#include "model.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "finite_difference.h"
#include "immersed_boundary.h"

Immersed_boundary::Immersed_boundary(Model* modelin, Input* inputin)
{
    model  = modelin;
    fields = model->fields;
    grid   = model->grid;
}

Immersed_boundary::~Immersed_boundary()
{
}

void Immersed_boundary::init()
{
    const int ii = 1;
    const int jj = grid->icells;
    const int kk = grid->ijcells;

    // Put a solid block at 1/2 of the horizontal dimension in the
    // middle of the channel.
    const int ibc_istart = grid->istart +  4 * (grid->iend-grid->istart) / 16;
    const int ibc_iend   = grid->istart + 12 * (grid->iend-grid->istart) / 16;
    const int ibc_kstart = grid->kstart +  3 * (grid->kend-grid->kstart) / 8;
    const int ibc_kend   = grid->kstart +  5 * (grid->kend-grid->kstart) / 8;

    IB_cell tmp;

    for (int k=ibc_kstart; k<ibc_kend; ++k)
        for (int j=grid->jstart; j<grid->jend; ++j)
            for (int i=ibc_istart; i<ibc_iend; ++i)
            {
                tmp.sw_outer_wall.reset();  // Set all walls to false
                tmp.i = i;
                tmp.j = j;
                tmp.k = k;

                if(i == ibc_istart)
                    tmp.sw_outer_wall.set(West_edge);
                if(i == ibc_iend-1)
                    tmp.sw_outer_wall.set(East_edge);
                if(k == ibc_kstart)
                    tmp.sw_outer_wall.set(Bottom_edge);
                if(k == ibc_kend-1)
                    tmp.sw_outer_wall.set(Top_edge);

                IB_cells.push_back(tmp);
            }
}

namespace
{
    void set_no_penetration_new(double* const restrict ut, double* const restrict wt,
                                double* const restrict u,  double* const restrict w,
                                std::vector<IB_cell> IB_cells,
                                const int icells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        for (std::vector<IB_cell>::iterator it=IB_cells.begin(); it<IB_cells.end(); ++it) // c++11 :(
        {
            const int ijk = it->i + it->j*jj + it->k*kk;

            if (it->sw_outer_wall[West_edge])
            {
                u [ijk] = 0;
                ut[ijk] = 0;
            }

            if (it->sw_outer_wall[East_edge])
            {
                u [ijk+ii] = 0;
                ut[ijk+ii] = 0;
            }

            if (it->sw_outer_wall[Bottom_edge])
            {
                w [ijk] = 0;
                wt[ijk] = 0;
            }

            if (it->sw_outer_wall[Top_edge])
            {
                w [ijk+kk] = 0;
                wt[ijk+kk] = 0;
            }
        }
    }

    void set_no_penetration(double* const restrict ut, double* const restrict wt,
                            double* const restrict u, double* const restrict w,
                            const int istart, const int iend,
                            const int jstart, const int jend,
                            const int kstart, const int kend,
                            const int icells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        // Put an solid block at 1/2 of the horizontal dimension in the
        // middle of the channel.
        const int ibc_istart = istart +  4 * (iend-istart) / 16;
        const int ibc_iend   = istart + 12 * (iend-istart) / 16;
        const int ibc_kstart = kstart +  3 * (kend-kstart) / 8;
        const int ibc_kend   = kstart +  5 * (kend-kstart) / 8;

        // Set the u ghost cells, no flow in the block.
        for (int k=ibc_kstart; k<ibc_kend; ++k)
            for (int j=jstart; j<jend; ++j)
            {
                const int ijk_istart  = ibc_istart + j*jj + k*kk;
                const int ijk_iend    = ibc_iend   + j*jj + k*kk;
                u [ijk_istart] = 0.;
                u [ijk_iend  ] = 0.;
                ut[ijk_istart] = 0.;
                ut[ijk_iend  ] = 0.;
            }

        // Set the w ghost cells, no flow in the block.
        for (int j=jstart; j<jend; ++j)
            for (int i=ibc_istart; i<ibc_iend; ++i)
            {
                const int ijk_kstart = i + j*jj + ibc_kstart*kk;
                const int ijk_kend   = i + j*jj + ibc_kend  *kk;
                w [ijk_kstart] = 0.;
                w [ijk_kend  ] = 0.;
                wt[ijk_kstart] = 0.;
                wt[ijk_kend  ] = 0.;
            }
    }

    void set_no_slip(double* const restrict ut, double* const restrict wt,
                     const double* const restrict u, const double* const restrict w,
                     const double* const restrict rhoref, const double* const restrict rhorefh,
                     const double* const restrict dzi, const double* const restrict dzhi,
                     const double dxi,
                     const double visc,
                     const int istart, const int iend,
                     const int jstart, const int jend,
                     const int kstart, const int kend,
                     const int icells, const int ijcells)
    {
        using namespace Finite_difference::O2;

        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        const double dxidxi = dxi*dxi;

        // Put an solid block at 1/2 of the horizontal dimension in the
        // middle of the channel.
        const int ibc_istart = istart +  4 * (iend-istart) / 16;
        const int ibc_iend   = istart + 12 * (iend-istart) / 16;
        const int ibc_kstart = kstart +  3 * (kend-kstart) / 8;
        const int ibc_kend   = kstart +  5 * (kend-kstart) / 8;

        // Set the w no slip at the vertical walls, by reverting the advection 
        // and diffusion towards the wall and adding the proper diffusion
        for (int k=ibc_kstart; k<ibc_kend+1; ++k)
            for (int j=jstart; j<jend; ++j)
            {
                int ijk = ibc_istart-1 + j*jj + k*kk;
                wt[ijk] +=
                        + ( interp2(u[ijk+ii-kk], u[ijk+ii]) * interp2(w[ijk], w[ijk+ii]) ) * dxi
                        - visc * ( (w[ijk+ii] - w[ijk]) ) * dxidxi
                        + visc * ( -2.*w[ijk] ) * dxidxi;

                ijk = ibc_iend + j*jj + k*kk;
                wt[ijk] +=
                        - ( interp2(u[ijk-kk], u[ijk]) * interp2(w[ijk-ii], w[ijk]) ) * dxi
                        + visc * ( (w[ijk] - w[ijk-ii]) ) * dxidxi
                        - visc * ( 2.*w[ijk] ) * dxidxi;
            }

        // Set the u no slip at the horizontal walls
        for (int j=jstart; j<jend; ++j)
            for (int i=ibc_istart; i<ibc_iend+1; ++i)
            {
                int k = ibc_kstart-1;
                int ijk = i + j*jj + k*kk;
                ut[ijk] +=
                        + ( rhorefh[k+1] * interp2(w[ijk-ii+kk], w[ijk+kk]) * interp2(u[ijk], u[ijk+kk]) ) / rhoref[k] * dzi[k];
                        - visc * ( (u[ijk+kk] - u[ijk]) * dzhi[k+1]) * dzi[k];
                        + visc * ( -2.*u[ijk] * dzhi[k+1] ) * dzi[k];

                k = ibc_kend;
                ijk = i + j*jj + k*kk;
                ut[ijk] +=
                        - ( rhorefh[k] * interp2(w[ijk-ii], w[ijk]) * interp2(u[ijk-kk], u[ijk]) ) / rhoref[k] * dzi[k];
                        + visc * ( (u[ijk] - u[ijk-kk]) * dzhi[k] ) * dzi[k];
                        - visc * ( 2.*u[ijk] * dzhi[k] ) * dzi[k];
            }
    }

    void set_scalar(double* const restrict st, double* const restrict s,
                    const double* const restrict u, const double* const restrict w,
                    const double* const restrict rhoref, const double* const restrict rhorefh,
                    const double* const restrict dzi, const double* const restrict dzhi,
                    const double dxi,
                    const double visc,
                    const int istart, const int iend,
                    const int jstart, const int jend,
                    const int kstart, const int kend,
                    const int icells, const int ijcells)
    {
        using namespace Finite_difference::O2;

        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        const double dxidxi = dxi*dxi;

        // Put an solid block at 1/2 of the horizontal dimension in the
        // middle of the channel.
        const int ibc_istart = istart +  4 * (iend-istart) / 16;
        const int ibc_iend   = istart + 12 * (iend-istart) / 16;
        const int ibc_kstart = kstart +  3 * (kend-kstart) / 8;
        const int ibc_kend   = kstart +  5 * (kend-kstart) / 8;

        // Set no flow through the object at the vertical wall and a neumann BC.
        for (int k=ibc_kstart; k<ibc_kend; ++k)
            for (int j=jstart; j<jend; ++j)
            {
                int ijk = ibc_istart-1 + j*jj + k*kk;
                st[ijk] +=
                         // + ( u[ijk+ii] * interp2(s[ijk], s[ijk+ii]) ) * dxi
                         - visc * ( s[ijk+ii] - s[ijk] ) * dxidxi;

                ijk = ibc_iend + j*jj + k*kk;
                st[ijk] +=
                         // - ( u[ijk] * interp2(s[ijk-ii], s[ijk]) ) * dxi
                         + visc * ( (s[ijk] - s[ijk-ii]) ) * dxidxi;
            }

        // Set no flow through the object at the horizontal wall
        for (int j=jstart; j<jend; ++j)
            for (int i=ibc_istart; i<ibc_iend; ++i)
            {
                int k = ibc_kstart-1;
                int ijk = i + j*jj + k*kk;
                st[ijk] +=
                         // + ( rhorefh[k+1] * w[ijk+kk] * interp2(s[ijk], s[ijk+kk]) ) / rhoref[k] * dzi[k]
                         - visc * (s[ijk+kk] - s[ijk]) * dzhi[k+1] * dzi[k];

                k = ibc_kend;
                ijk = i + j*jj + k*kk;
                st[ijk] +=
                         // - ( rhorefh[k] * w[ijk] * interp2(s[ijk-kk], s[ijk]) ) / rhoref[k] * dzi[k];
                         + visc * (s[ijk] - s[ijk-kk]) * dzhi[k] * dzi[k];
            }
    }

}

void Immersed_boundary::exec()
{
    //set_no_penetration(fields->ut->data, fields->wt->data,
    //                   fields->u->data, fields->w->data,
    //                   grid->istart, grid->iend,
    //                   grid->jstart, grid->jend,
    //                   grid->kstart, grid->kend,
    //                   grid->icells, grid->ijcells);

    set_no_penetration_new(fields->ut->data, fields->wt->data,
                           fields->u->data, fields->w->data, IB_cells,
                           grid->icells, grid->ijcells);

    set_no_slip(fields->ut->data, fields->wt->data,
                fields->u->data, fields->w->data,
                fields->rhoref, fields->rhorefh,
                grid->dzi, grid->dzhi,
                grid->dxi,
                fields->visc,
                grid->istart, grid->iend,
                grid->jstart, grid->jend,
                grid->kstart, grid->kend,
                grid->icells, grid->ijcells);

    for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
        set_scalar(it->second->data, fields->sp[it->first]->data,
                   fields->u->data, fields->w->data,
                   fields->rhoref, fields->rhorefh,
                   grid->dzi, grid->dzhi,
                   grid->dxi,
                   fields->sp[it->first]->visc,
                   grid->istart, grid->iend,
                   grid->jstart, grid->jend,
                   grid->kstart, grid->kend,
                   grid->icells, grid->ijcells);
}
