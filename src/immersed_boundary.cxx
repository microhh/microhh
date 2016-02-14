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
#include <cmath>
#include "master.h"
#include "model.h"
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "stats.h"
#include "finite_difference.h"
#include "immersed_boundary.h"

Immersed_boundary::Immersed_boundary(Model* modelin, Input* inputin)
{
    model  = modelin;
    fields = model->fields;
    grid   = model->grid;

    inputin->get_item(&sw_ib, "boundary", "sw_ib", "", "0");
}

Immersed_boundary::~Immersed_boundary()
{
}

namespace
{
    void add_IB_block(std::vector<IB_cell> &IB_cells,
                      const int istart, const int iend,
                      const int jstart, const int jend,
                      const int kstart, const int kend)
    {
        IB_cell tmp;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    tmp.sw_outer_wall.reset();  // Set all walls to false

                    tmp.i = i;
                    tmp.j = j;
                    tmp.k = k;

                    if(i == istart)
                        tmp.sw_outer_wall.set(West_edge);
                    if(i == iend-1)
                        tmp.sw_outer_wall.set(East_edge);
                    if(j == jstart)
                        tmp.sw_outer_wall.set(South_edge);
                    if(j == jend-1)
                        tmp.sw_outer_wall.set(North_edge);
                    if(k == kstart)
                        tmp.sw_outer_wall.set(Bottom_edge);
                    if(k == kend-1)
                        tmp.sw_outer_wall.set(Top_edge);

                    IB_cells.push_back(tmp);
                }
        }

    void set_no_penetration_new(double* const restrict ut, double* const restrict vt, double* const restrict wt,
                                double* const restrict u,  double* const restrict v, double* const restrict w,
                                std::vector<IB_cell> IB_cells,
                                const int icells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        for (std::vector<IB_cell>::iterator it=IB_cells.begin(); it<IB_cells.end(); ++it) // c++11 :(
        {
            const int ijk = it->i + it->j*jj + it->k*kk;

            u [ijk] = 0;
            //v [ijk] = 0;
            w [ijk] = 0;
            ut[ijk] = 0;
            //vt[ijk] = 0;
            wt[ijk] = 0;

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

            //if (it->sw_outer_wall[South_edge])
            //{
            //    v [ijk] = 0;
            //    vt[ijk] = 0;
            //}

            //if (it->sw_outer_wall[North_edge])
            //{
            //    v [ijk+jj] = 0;
            //    vt[ijk+jj] = 0;
            //}

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
                            double* const restrict u,  double* const restrict w,
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

    void set_no_slip_new(double* const restrict ut, double* const restrict wt,
                         const double* const restrict u, const double* const restrict w,
                         const double* const restrict rhoref, const double* const restrict rhorefh,
                         const double* const restrict dzi, const double* const restrict dzhi,
                         const double dxi,
                         const double visc,
                         std::vector<IB_cell> IB_cells,
                         const int icells, const int ijcells)
    {
        using namespace Finite_difference::O2;

        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        const double dxidxi = dxi*dxi;

        for (std::vector<IB_cell>::iterator it=IB_cells.begin(); it<IB_cells.end(); ++it) // c++11 :(
        {
            // Adjust u-component directly below cell
            if (it->sw_outer_wall[Bottom_edge])
            {
                const int k   = it->k-1;
                const int ijk = it->i + it->j*jj + k*kk; 

                ut[ijk] +=
                        + ( rhorefh[k+1] * interp2(w[ijk-ii+kk], w[ijk+kk]) * interp2(u[ijk], u[ijk+kk]) ) / rhoref[k] * dzi[k];
                        - visc * ( (u[ijk+kk] - u[ijk]) * dzhi[k+1]) * dzi[k];
                        + visc * ( -2.*u[ijk] * dzhi[k+1] ) * dzi[k];

                // For the bottom-right corner, also adjust u-component
                if(it->sw_outer_wall[East_edge])
                {
                    const int ijk = it->i+1 + it->j*jj + k*kk;

                    ut[ijk] +=
                            + ( rhorefh[k+1] * interp2(w[ijk-ii+kk], w[ijk+kk]) * interp2(u[ijk], u[ijk+kk]) ) / rhoref[k] * dzi[k];
                            - visc * ( (u[ijk+kk] - u[ijk]) * dzhi[k+1]) * dzi[k];
                            + visc * ( -2.*u[ijk] * dzhi[k+1] ) * dzi[k];
                }
            }

            // Adjust u-component directly above cell
            if (it->sw_outer_wall[Top_edge])
            {
                const int k   = it->k+1;
                const int ijk = it->i + it->j*jj + k*kk; 

                ut[ijk] +=
                        - ( rhorefh[k] * interp2(w[ijk-ii], w[ijk]) * interp2(u[ijk-kk], u[ijk]) ) / rhoref[k] * dzi[k];
                        + visc * ( (u[ijk] - u[ijk-kk]) * dzhi[k] ) * dzi[k];
                        - visc * ( 2.*u[ijk] * dzhi[k] ) * dzi[k];

                // For the top-right corner, also adjust u-component
                if(it->sw_outer_wall[East_edge])
                {
                    const int ijk = it->i+1 + it->j*jj + k*kk; 

                    ut[ijk] +=
                            - ( rhorefh[k] * interp2(w[ijk-ii], w[ijk]) * interp2(u[ijk-kk], u[ijk]) ) / rhoref[k] * dzi[k];
                            + visc * ( (u[ijk] - u[ijk-kk]) * dzhi[k] ) * dzi[k];
                            - visc * ( 2.*u[ijk] * dzhi[k] ) * dzi[k];
                }
            }

            // Adjust w-component west of cell  
            if (it->sw_outer_wall[West_edge])
            {
                const int k   = it->k;
                const int ijk = it->i-1 + it->j*jj + k*kk; 

                wt[ijk] +=
                        + ( interp2(u[ijk+ii-kk], u[ijk+ii]) * interp2(w[ijk], w[ijk+ii]) ) * dxi
                        - visc * ( (w[ijk+ii] - w[ijk]) ) * dxidxi
                        + visc * ( -2.*w[ijk] ) * dxidxi;
                
                // For the top-left corner, also adjust w-component
                if (it->sw_outer_wall[Top_edge])
                {
                    const int k   = it->k+1;
                    const int ijk = it->i-1 + it->j*jj + k*kk; 

                    wt[ijk] +=
                            + ( interp2(u[ijk+ii-kk], u[ijk+ii]) * interp2(w[ijk], w[ijk+ii]) ) * dxi
                            - visc * ( (w[ijk+ii] - w[ijk]) ) * dxidxi
                            + visc * ( -2.*w[ijk] ) * dxidxi;
                }
            }

            // Adjust w-component east of cell  
            if (it->sw_outer_wall[East_edge])
            {
                const int k   = it->k;
                const int ijk = it->i+1 + it->j*jj + k*kk; 

                wt[ijk] +=
                        - ( interp2(u[ijk-kk], u[ijk]) * interp2(w[ijk-ii], w[ijk]) ) * dxi
                        + visc * ( (w[ijk] - w[ijk-ii]) ) * dxidxi
                        - visc * ( 2.*w[ijk] ) * dxidxi;
                
                // For the top-right corner, also adjust w-component
                if (it->sw_outer_wall[Top_edge])
                {
                    const int k   = it->k+1;
                    const int ijk = it->i+1 + it->j*jj + k*kk; 

                    wt[ijk] +=
                            - ( interp2(u[ijk-kk], u[ijk]) * interp2(w[ijk-ii], w[ijk]) ) * dxi
                            + visc * ( (w[ijk] - w[ijk-ii]) ) * dxidxi
                            - visc * ( 2.*w[ijk] ) * dxidxi;
                }
            }
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

    void set_scalar_new(double* const restrict st, double* const restrict s,
                        const double* const restrict u, const double* const restrict w,
                        const double* const restrict rhoref, const double* const restrict rhorefh,
                        const double* const restrict dzi, const double* const restrict dzhi,
                        const double dxi,
                        const double visc,
                        std::vector<IB_cell> IB_cells,
                        const int icells, const int ijcells)
    {
        using namespace Finite_difference::O2;

        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        const double dxidxi = dxi*dxi;

        for (std::vector<IB_cell>::iterator it=IB_cells.begin(); it<IB_cells.end(); ++it) // c++11 :(
        {
            if (it->sw_outer_wall[West_edge])
            {
                const int ijk = it->i-1 + it->j*jj + it->k*kk;
                st[ijk] -= visc * ( s[ijk+ii] - s[ijk] ) * dxidxi;
            }
       
            if (it->sw_outer_wall[East_edge])
            {
                const int ijk = it->i+1 + it->j*jj + it->k*kk;
                st[ijk] += visc * ( (s[ijk] - s[ijk-ii]) ) * dxidxi;
            }
      
            if(it->sw_outer_wall[Bottom_edge])
            {
                const int k   = it->k-1;
                const int ijk = it->i + it->j*jj + k*kk;
                st[ijk] -= visc * (s[ijk+kk] - s[ijk]) * dzhi[k+1] * dzi[k];
            }

            if(it->sw_outer_wall[Top_edge])
            {
                const int k   = it->k+1;
                const int ijk = it->i + it->j*jj + k*kk;
                st[ijk] += visc * (s[ijk] - s[ijk-kk]) * dzhi[k] * dzi[k];
            }
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

    void check_no_penetration(double* u_max, double* w_max, 
                              const double* const restrict u, const double* const restrict w, 
                              std::vector<IB_cell> IB_cells,
                              const int icells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        *u_max = 0;
        *w_max = 0;

        for (std::vector<IB_cell>::iterator it=IB_cells.begin(); it<IB_cells.end(); ++it) // c++11 :(
        {
            const int ijk = it->i + it->j*jj + it->k*kk;

            if (it->sw_outer_wall[West_edge])
                if(std::abs(u[ijk]) > *u_max)
                    *u_max = std::abs(u[ijk]);

            if (it->sw_outer_wall[East_edge])
                if(std::abs(u[ijk+ii]) > *u_max)
                    *u_max = std::abs(u[ijk+ii]);

            if (it->sw_outer_wall[Bottom_edge])
                if(std::abs(w[ijk]) > *w_max)
                    *w_max = std::abs(w[ijk]);

            if (it->sw_outer_wall[Top_edge])
                if(std::abs(w[ijk+kk]) > *w_max)
                    *w_max = std::abs(w[ijk+kk]);
        }
    }
}

void Immersed_boundary::init()
{
    stats = model->stats;

    if (sw_ib != "1")
        return;

    int ibc_istart;
    int ibc_iend;
    int ibc_jstart;
    int ibc_jend;
    int ibc_kstart;
    int ibc_kend;

    // Inflow grid
    //const int nvert=16;
    //for (int k=0; k<grid->ktot/nvert; ++k)
    //{
    //    if(k%2==0)
    //    {
    //        ibc_kstart = grid->kstart + k*nvert+nvert/2;
    //        ibc_kend   = grid->kstart + (k+1)*nvert+nvert/2;
    //        add_IB_block(IB_cells, 64, 70, grid->jstart, grid->jend+1, ibc_kstart, ibc_kend);
    //    }
    //}

    // Set of 2d building...
    //const int w = grid->itot/15;

    //ibc_istart = grid->istart + 0.3*grid->itot - w/2;
    //ibc_iend   = grid->istart + 0.3*grid->itot + w/2;
    //ibc_kstart = grid->kstart;
    //ibc_kend   = grid->kstart + grid->ktot/5;
    //add_IB_block(IB_cells, ibc_istart, ibc_iend, grid->jstart, grid->jend+1, ibc_kstart, ibc_kend);

    //ibc_istart = grid->istart + 0.5*grid->itot - w/2;
    //ibc_iend   = grid->istart + 0.5*grid->itot + w/2;
    //ibc_kstart = grid->kstart;
    //ibc_kend   = grid->kstart + grid->ktot/3;
    //add_IB_block(IB_cells, ibc_istart, ibc_iend, grid->jstart, grid->jend+1, ibc_kstart, ibc_kend);

    //ibc_istart = grid->istart + 0.7*grid->itot - w/2;
    //ibc_iend   = grid->istart + 0.7*grid->itot + w/2;
    //ibc_kstart = grid->kstart;
    //ibc_kend   = grid->kstart + grid->ktot/5;
    //add_IB_block(IB_cells, ibc_istart, ibc_iend, grid->jstart, grid->jend+1, ibc_kstart, ibc_kend);

    // Block in channel
    ibc_istart = grid->istart +  6 * (grid->iend-grid->istart) / 16;
    ibc_iend   = grid->istart + 10 * (grid->iend-grid->istart) / 16;
    ibc_jstart = grid->jstart;
    ibc_jend   = grid->jend+1;
    ibc_kstart = grid->kstart;
    ibc_kend   = grid->kstart+ 0.2 * (grid->kend-grid->kstart);
    add_IB_block(IB_cells, ibc_istart, ibc_iend, ibc_jstart, ibc_jend, ibc_kstart, ibc_kend);
}

void Immersed_boundary::create()
{
    if (sw_ib != "1")
        return;

    if (stats->getSwitch() == "1")
    {
        stats->add_time_series("IB_u_max", "Maximum u-component flow through IB walls", "m s-1");
        stats->add_time_series("IB_w_max", "Maximum w-component flow through IB walls", "m s-1");
    }
}

void Immersed_boundary::exec_stats(Mask *m)
{
    if (sw_ib != "1")
        return;
    
    // Temporary hack (BvS) to get wall velocities in statistics, doesn't account for mask!!
    m->tseries["IB_u_max"].data = boundary_u_max;
    m->tseries["IB_w_max"].data = boundary_w_max;
}

void Immersed_boundary::exec()
{
    if (sw_ib != "1")
        return;

    if (stats->getSwitch() == "1")
        check_no_penetration(&boundary_u_max, &boundary_w_max, fields->u->data, fields->w->data, IB_cells, grid->icells, grid->ijcells);

    set_no_penetration_new(fields->ut->data, fields->vt->data, fields->wt->data,
                           fields->u->data, fields->v->data, fields->w->data, 
                           IB_cells, grid->icells, grid->ijcells);

    set_no_slip_new(fields->ut->data, fields->wt->data,
                    fields->u->data, fields->w->data,
                    fields->rhoref, fields->rhorefh,
                    grid->dzi, grid->dzhi,
                    grid->dxi,
                    fields->visc,
                    IB_cells,
                    grid->icells, grid->ijcells);

    for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
        set_scalar_new(it->second->data, fields->sp[it->first]->data,
                       fields->u->data, fields->w->data,
                       fields->rhoref, fields->rhorefh,
                       grid->dzi, grid->dzhi,
                       grid->dxi,
                       fields->sp[it->first]->visc,
                       IB_cells,
                       grid->icells, grid->ijcells);

    //set_no_penetration(fields->ut->data, fields->wt->data,
    //                   fields->u->data, fields->w->data,
    //                   grid->istart, grid->iend,
    //                   grid->jstart, grid->jend,
    //                   grid->kstart, grid->kend,
    //                   grid->icells, grid->ijcells);

    //set_no_slip(fields->ut->data, fields->wt->data,
    //            fields->u->data, fields->w->data,
    //            fields->rhoref, fields->rhorefh,
    //            grid->dzi, grid->dzhi,
    //            grid->dxi,
    //            fields->visc,
    //            grid->istart, grid->iend,
    //            grid->jstart, grid->jend,
    //            grid->kstart, grid->kend,
    //            grid->icells, grid->ijcells);

    //for (FieldMap::const_iterator it = fields->st.begin(); it!=fields->st.end(); it++)
    //    set_scalar(it->second->data, fields->sp[it->first]->data,
    //               fields->u->data, fields->w->data,
    //               fields->rhoref, fields->rhorefh,
    //               grid->dzi, grid->dzhi,
    //               grid->dxi,
    //               fields->sp[it->first]->visc,
    //               grid->istart, grid->iend,
    //               grid->jstart, grid->jend,
    //               grid->kstart, grid->kend,
    //               grid->icells, grid->ijcells);
}


