/*
 * MicroHH
 * Copyright (c) 2011-2020 Chiel van Heerwaarden
 * Copyright (c) 2011-2020 Thijs Heus
 * Copyright (c) 2014-2020 Bart van Stratum
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
#include <sstream>
#include <algorithm>

#include "boundary_lateral.h"
#include "netcdf_interface.h"
#include "grid.h"
#include "fields.h"
#include "input.h"
#include "master.h"
#include "timeloop.h"

namespace
{
    template<typename TF>
    bool in_list(const TF value, const std::vector<TF>& list)
    {
        if(std::find(list.begin(), list.end(), value) != list.end())
            return true;
        else
            return false;
    }


    //template<typename TF, Lbc_location location>
    //void set_ghost_cell_kernel_u(
    //        TF* const restrict u,
    //        const TF* const restrict lbc_u,
    //        const int istart, const int iend, const int igc,
    //        const int jstart, const int jend,
    //        const int kstart, const int kend,
    //        const int icells, const int jcells,
    //        const int ijcells)
    //{
    //    for (int k=kstart; k<kend; ++k)
    //        for (int j=jstart; j<jend; ++j)
    //        {
    //            const int jk = j+k*jcells;

    //            // Set boundary values directly.
    //            if (location == Lbc_location::West)
    //            {
    //                const int ijk_b = istart + j*icells + k*ijcells;
    //                const int ijk_d = istart+1 + j*icells +k*ijcells;

    //                u[ijk_b] = lbc_u[jk];

    //                for (int i=0; i<igc; ++i)
    //                {
    //                    const int ijk_gc = ijk_b - (i+1);
    //                    u[ijk_gc] = lbc_u[jk] - (i+1)*(u[ijk_d]-lbc_u[jk]);
    //                }
    //            }
    //            else if (location == Lbc_location::East)
    //            {
    //                const int ijk_b = iend + j*icells + k*ijcells;
    //                const int ijk_d = iend-1 + j*icells +k*ijcells;

    //                u[ijk_b] = lbc_u[jk];

    //                for (int i=0; i<igc-1; ++i)
    //                {
    //                    const int ijk_gc = ijk_b + (i+1);
    //                    u[ijk_gc] = lbc_u[jk] + (i+1)*(lbc_u[jk]-u[ijk_d]);
    //                }
    //            }
    //        }
    //}

    //template<typename TF, Lbc_location location>
    //void set_ghost_cell_kernel_v(
    //        TF* const restrict v,
    //        const TF* const restrict lbc_v,
    //        const int istart, const int iend,
    //        const int jstart, const int jend, const int jgc,
    //        const int kstart, const int kend,
    //        const int icells, const int jcells,
    //        const int ijcells)
    //{
    //    for (int k=kstart; k<kend; ++k)
    //        for (int i=istart; i<iend; ++i)
    //        {
    //            const int ik = i+k*icells;

    //            // Set boundary values directly.
    //            if (location == Lbc_location::South)
    //            {
    //                const int ijk_b = i + jstart*icells + k*ijcells;
    //                const int ijk_d = i + (jstart+1)*icells +k*ijcells;

    //                v[ijk_b] = lbc_v[ik];

    //                for (int j=0; j<jgc; ++j)
    //                {
    //                    const int ijk_gc = ijk_b - (j+1)*icells;
    //                    v[ijk_gc] = lbc_v[ik] - (j+1)*(v[ijk_d]-lbc_v[ik]);
    //                }
    //            }
    //            else if (location == Lbc_location::North)
    //            {
    //                const int ijk_b = i + jend*icells + k*ijcells;
    //                const int ijk_d = i + (jend-1)*icells +k*ijcells;

    //                v[ijk_b] = lbc_v[ik];

    //                for (int j=0; j<jgc-1; ++j)
    //                {
    //                    const int ijk_gc = ijk_b + (j+1)*icells;
    //                    v[ijk_gc] = lbc_v[ik] + (j+1)*(lbc_v[ik]-v[ijk_d]);
    //                }
    //            }
    //        }
    //}

    // This kernel enforces a Neumann BC of 0 on w.
    template<typename TF, Lbc_location location>
    void set_ghost_cell_kernel_w(
            TF* const restrict a,
            const int istart, const int iend, const int igc,
            const int jstart, const int jend, const int jgc,
            const int kstart, const int kend,
            const int icells, const int jcells, const int kcells,
            const int ijcells)
    {
        int ijk;
        int ijk_gc;
        int ijk_d;

        // Set the ghost cells using extrapolation.
        if (location == Lbc_location::West || location == Lbc_location::East)
        {
            for (int k=kstart; k<kend; ++k)
                for (int j=jstart; j<jend; ++j)
                    for (int i=0; i<igc; ++i)
                    {
                        if (location == Lbc_location::West)
                        {
                            ijk_d  = (istart    ) + j*icells + k*ijcells;
                            ijk_gc = (istart-1-i) + j*icells + k*ijcells;
                        }
                        else if (location == Lbc_location::East)
                        {
                            ijk_d  = (iend-1  ) + j*icells + k*ijcells;
                            ijk_gc = (iend+i  ) + j*icells + k*ijcells;
                        }

                        a[ijk_gc] = a[ijk_d];
                    }
        }
        else if (location == Lbc_location::North || location == Lbc_location::South)
        {
            for (int k=kstart; k<kend; ++k)
                for (int i=istart; i<iend; ++i)
                    for (int j=0; j<jgc; ++j)
                    {
                        if (location == Lbc_location::South)
                        {
                            ijk_d  = i + (jstart    )*icells + k*ijcells;
                            ijk_gc = i + (jstart-1-j)*icells + k*ijcells;
                        }
                        else if (location == Lbc_location::North)
                        {
                            ijk_d  = i + (jend-1  )*icells + k*ijcells;
                            ijk_gc = i + (jend+j  )*icells + k*ijcells;
                        }

                        a[ijk_gc] = a[ijk_d];
                    }
        }
    }


    //template<typename TF, Lbc_location location>
    //void set_ghost_cell_kernel_s(
    //        TF* const restrict a,
    //        const TF* const restrict lbc,
    //        const int istart, const int iend, const int igc,
    //        const int jstart, const int jend, const int jgc,
    //        const int kstart, const int kend,
    //        const int icells, const int jcells, const int kcells,
    //        const int ijcells)
    //{
    //    int ijk;
    //    int ijk_gc;
    //    int ijk_d;

    //    // Set the ghost cells using extrapolation.
    //    if (location == Lbc_location::West || location == Lbc_location::East)
    //    {
    //        for (int k=kstart; k<kend; ++k)
    //            for (int j=jstart; j<jend; ++j)
    //            {
    //                const int jk = j+k*jcells;

    //                for (int i=0; i<igc; ++i)
    //                {
    //                    if (location == Lbc_location::West)
    //                    {
    //                        ijk_d  = (istart    ) + j*icells + k*ijcells;
    //                        ijk_gc = (istart-1-i) + j*icells + k*ijcells;
    //                    }
    //                    else if (location == Lbc_location::East)
    //                    {
    //                        ijk_d  = (iend-1  ) + j*icells + k*ijcells;
    //                        ijk_gc = (iend+i  ) + j*icells + k*ijcells;
    //                    }

    //                    a[ijk_gc] = a[ijk_d] - (i+1)*TF(2)*(a[ijk_d] - lbc[jk]);
    //                }
    //            }

    //    }
    //    else if (location == Lbc_location::North || location == Lbc_location::South)
    //    {
    //        for (int k=kstart; k<kend; ++k)
    //            for (int i=istart; i<iend; ++i)
    //            {
    //                const int ik = i+k*icells;

    //                for (int j=0; j<jgc; ++j)
    //                {
    //                    if (location == Lbc_location::South)
    //                    {
    //                        ijk_d  = i + (jstart    )*icells + k*ijcells;
    //                        ijk_gc = i + (jstart-1-j)*icells + k*ijcells;
    //                    }
    //                    else if (location == Lbc_location::North)
    //                    {
    //                        ijk_d  = i + (jend-1  )*icells + k*ijcells;
    //                        ijk_gc = i + (jend+j  )*icells + k*ijcells;
    //                    }

    //                    const TF lbc_val = lbc[ik];
    //                    a[ijk_gc] = a[ijk_d] - (j+1)*TF(2)*(a[ijk_d] - lbc[ik]);
    //                }
    //            }
    //    }
    //}

    template<typename TF>
    TF diffusion_3x3x3(
        const TF* const restrict fld,
        const int ijk,
        const int icells,
        const int ijcells)
    {
        auto index = [&](
                const int i3, const int j3, const int k3)
        {
            return ijk + i3-1 + (j3-1)*icells + (k3-1)*ijcells;
        };

        const TF fld_diff =
                - TF(1) * fld[index(0,0,0)] + TF(2) * fld[index(0,1,0)] - TF(1) * fld[index(0,2,0)]
                + TF(2) * fld[index(1,0,0)] - TF(4) * fld[index(1,1,0)] + TF(2) * fld[index(1,2,0)]
                - TF(1) * fld[index(2,0,0)] + TF(2) * fld[index(2,1,0)] - TF(1) * fld[index(2,2,0)]
                + TF(2) * fld[index(0,0,1)] - TF(4) * fld[index(0,1,1)] + TF(2) * fld[index(0,2,1)]
                - TF(4) * fld[index(1,0,1)] + TF(8) * fld[index(1,1,1)] - TF(4) * fld[index(1,2,1)]
                + TF(2) * fld[index(2,0,1)] - TF(4) * fld[index(2,1,1)] + TF(2) * fld[index(2,2,1)]
                - TF(1) * fld[index(0,0,2)] + TF(2) * fld[index(0,1,2)] - TF(1) * fld[index(0,2,2)]
                + TF(2) * fld[index(1,0,2)] - TF(4) * fld[index(1,1,2)] + TF(2) * fld[index(1,2,2)]
                - TF(1) * fld[index(2,0,2)] + TF(2) * fld[index(2,1,2)] - TF(1) * fld[index(2,2,2)];

        return fld_diff;
    }

    template<typename TF, Lbc_location location>
    void lateral_sponge_kernel_u(
            TF* const restrict ut,
            const TF* const restrict u,
            const TF* const restrict lbc_u,
            const TF tau_nudge,
            const TF w_diff,
            const int N_sponge,
            const int npy,
            const int mpiidy,
            const int igc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells,
            const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        const TF w_dt = TF(1) / tau_nudge;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
            {
                const int ilbc = (location==Lbc_location::West) ? igc : 0;
                const int jk = ilbc + j*igc + k*igc*jcells;

                for (int n=2; n<=N_sponge; ++n)
                {
                    const int i = (location==Lbc_location::West) ? istart+(n-1) : iend-(n-1);
                    const int ijk = i + j*icells + k*ijcells;

                    // Calculate diffusion term over 3x3x3 stencil.
                    // Offset block near lateral boundaries to avoid using ghost cells.
                    // No offset needed for `i`, as stencil center starts at `istart+1` or `iend-1`,
                    // and `istart` and `iend` contain the correct boundary values.
                    //const int jo =
                    //        (mpiidy == 0 && j == jstart) ? 1 :
                    //        (mpiidy == npy-1 && j == jend-1) ? -1 : 0;

                    const int ko =
                            (k == kstart) ? 1 :
                            (k == kend-1) ? -1 : 0;

                    //const int ijkc = i + (j+jo)*icells + (k+ko)*(ijcells);
                    const int ijkc = i + j*icells + (k+ko)*(ijcells);

                    const TF u_diff = diffusion_3x3x3(
                            u, ijkc, icells, ijcells);

                    // Nudge coefficient.
                    const TF w1 = w_dt * (TF(1)+N_sponge-n) / N_sponge;

                    ut[ijk] += w1*(lbc_u[jk]-u[ijk]);
                    ut[ijk] -= w_diff*u_diff;
                }
            }
    }

    template<typename TF, Lbc_location location>
    void lateral_sponge_kernel_v(
            TF* const restrict vt,
            const TF* const restrict v,
            const TF* const restrict lbc_v,
            const TF tau_nudge,
            const TF w_diff,
            const int N_sponge,
            const int npx,
            const int mpiidx,
            const int jgc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells,
            const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        const TF w_dt = TF(1) / tau_nudge;

        for (int k=kstart; k<kend; ++k)
            for (int i=istart; i<iend; ++i)
            {
                const int jlbc = (location==Lbc_location::South) ? jgc : 0;
                const int ik = i + jlbc*icells + k*icells*jgc;

                for (int n=2; n<=N_sponge; ++n)
                {
                    const int j = (location==Lbc_location::South) ? jstart+(n-1) : jend-(n-1);
                    const int ijk = i + j*icells + k*ijcells;

                    // Calculate diffusion term over 3x3x3 stencil.
                    // Offset block near lateral boundaries to avoid using ghost cells.
                    // No offset needed for `j`, as stencil center starts at `jstart+1` or `jend-1`,
                    // and `jstart` and `jend` contain the correct boundary values.
                    //const int io =
                    //        (mpiidx == 0 && i == istart) ? 1 :
                    //        (mpiidx == npx-1 && i == iend-1) ? -1 : 0;

                    const int ko =
                            (k == kstart) ? 1 :
                            (k == kend-1) ? -1 : 0;

                    //const int ijkc = i+io + j*icells + (k+ko)*(ijcells);
                    const int ijkc = i + j*icells + (k+ko)*(ijcells);

                    const TF v_diff = diffusion_3x3x3(
                            v, ijkc, icells, ijcells);

                    // Nudge coefficient.
                    const TF w1 = w_dt * (TF(1)+N_sponge-n) / N_sponge;

                    vt[ijk] += w1*(lbc_v[ik]-v[ijk]);
                    vt[ijk] -= w_diff*v_diff;
                }
            }
    }

    template<typename TF, Lbc_location location, bool sw_recycle>
    void lateral_sponge_kernel_s(
            TF* const restrict at,
            const TF* const restrict a,
            const TF* const restrict lbc,
            const TF tau_nudge,
            const TF w_diff,
            const int N_sponge,
            const TF tau_recycle,
            const TF recycle_offset,
            const int npx, const int npy,
            const int mpiidx, const int mpiidy,
            const int igc, const int jgc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells,
            const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        const TF w_dt = TF(1) / tau_nudge;
        const TF r_dt = TF(1) / tau_recycle;

        if (location == Lbc_location::West || location == Lbc_location::East)
        {
            for (int k=kstart; k<kend; ++k)
                for (int n=1; n<=N_sponge; ++n)
                {
                    // Offset in y-direction for domain corners.
                    const int jstart_loc = mpiidy == 0     ? jstart + (n-1) : jstart;
                    const int jend_loc   = mpiidy == npy-1 ? jend   - (n-1) : jend;

                    for (int j=jstart_loc; j<jend_loc; ++j)
                    {
                        // Index in LBC:
                        const int ilbc = (location==Lbc_location::West) ? igc-1 : 0;
                        const int jk = ilbc + j*igc + k*igc*jcells;

                        const int i = (location==Lbc_location::West) ? istart+(n-1) : iend-n;
                        const int ijk = i + j*icells + k*ijcells;

                        // Calculate diffusion term over 3x3x3 stencil.
                        // Offset block near lateral boundaries to avoid using ghost cells.
                        //const int io =
                        //        (location == Lbc_location::West && n==1) ? 1 :
                        //        (location == Lbc_location::East && n==1) ? -1 : 0;

                        //const int jo =
                        //        (mpiidy == 0 && j == jstart) ? 1 :
                        //        (mpiidy == npy-1 && j == jend-1) ? -1 : 0;

                        const int ko =
                                (k == kstart) ? 1 :
                                (k == kend-1) ? -1 : 0;

                        const int ijkc = i + j*icells + (k+ko)*(ijcells);
                        //const int ijkc = (i+io) + (j+jo)*icells + (k+ko)*(ijcells);

                        const TF a_diff = diffusion_3x3x3(
                                a, ijk, icells, ijcells);

                        // Nudge coefficient.
                        const TF f_sponge = (TF(1)+N_sponge-(n+TF(0.5))) / N_sponge;
                        const TF w1 = w_dt * f_sponge;

                        const TF lbc_val = lbc[jk];
                        const TF fld_val = a[ijk];

                        // Apply nudge and sponge tendencies.
                        at[ijk] += w1*(lbc[jk]-a[ijk]);
                        at[ijk] -= w_diff*a_diff;

                        // Turbulence recycling.
                        //if (sw_recycle)
                        //{
                        //    // Recycle strength; 0 at boundary, 1 at edge nudging zone.
                        //    const TF f_recycle = TF(1) - f_sponge;

                        //    // Source of recycling.
                        //    const int offset = (location == Lbc_location::West) ? recycle_offset : -recycle_offset;
                        //    const int ijko = (i+offset) + j*icells + k*ijcells;

                        //    TF a_mean = 0;
                        //    for (int jc=-3; jc<4; ++jc)
                        //        for (int ic=-3; ic<4; ++ic)
                        //            a_mean += a[ijko + ic + jc*icells];
                        //    a_mean /= TF(49.);

                        //    at[ijk] += f_recycle * r_dt * ((a[ijko] - a_mean) - (a[ijk] - lbc[jk]));
                        //}
                    }
                }
        }
        else if (location == Lbc_location::South || location == Lbc_location::North)
        {
            for (int k=kstart; k<kend; ++k)
                for (int n=1; n<=N_sponge; ++n)
                {
                    const int istart_loc = (mpiidx == 0)     ? istart + n : istart;
                    const int iend_loc   = (mpiidx == npx-1) ? iend   - n : iend;

                    for (int i=istart_loc; i<iend_loc; ++i)
                    {
                        const int jlbc = (location==Lbc_location::South) ? jgc-1 : 0;
                        const int ik = i + jlbc*icells + k*icells*jgc;

                        const int j = (location==Lbc_location::South) ? jstart+(n-1) : jend-n;
                        const int ijk = i + j*icells + k*ijcells;

                        // Calculate diffusion term over 3x3x3 stencil.
                        // Offset block near lateral boundaries to avoid using ghost cells.
                        //const int io =
                        //        (mpiidx == 0 && i == istart) ? 1 :
                        //        (mpiidx == npx-1 && i == iend-1) ? -1 : 0;

                        //const int jo =
                        //        (location == Lbc_location::South && n==1) ? 1 :
                        //        (location == Lbc_location::North && n==1) ? -1 : 0;

                        const int ko =
                                (k == kstart) ? 1 :
                                (k == kend-1) ? -1 : 0;

                        //const int ijkc = (i+io) + (j+jo)*icells + (k+ko)*(ijcells);
                        const int ijkc = i + j*icells + (k+ko)*(ijcells);

                        const TF a_diff = diffusion_3x3x3(
                                a, ijkc, icells, ijcells);

                        // Nudge coefficient.
                        const TF f_sponge = (TF(1)+N_sponge-(n+TF(0.5))) / N_sponge;
                        const TF w1 = w_dt * f_sponge;

                        at[ijk] += w1*(lbc[ik]-a[ijk]);
                        at[ijk] -= w_diff*a_diff;

                        //if (sw_recycle)
                        //{
                        //    // Recycle strength; 0 at boundary, 1 at edge nudging zone.
                        //    const TF f_recycle = TF(1) - f_sponge;

                        //    // Source of recycling.
                        //    const int offset = (location == Lbc_location::South) ? recycle_offset : -recycle_offset;
                        //    const int ijko = i + (j+offset)*icells + k*ijcells;

                        //    TF a_mean = 0;
                        //    for (int jc=-3; jc<4; ++jc)
                        //        for (int ic=-3; ic<4; ++ic)
                        //            a_mean += a[ijko + ic + jc*icells];
                        //    a_mean /= TF(49.);

                        //    at[ijk] += f_recycle * r_dt * ((a[ijko] - a_mean) - (a[ijk] - lbc[ik]));
                        //}
                    }
                }
        }
    }


    template<typename TF>
    void set_corner_ghost_cell_kernel(
            TF* const restrict fld,
            const int mpiidx, const int mpiidy,
            const int npx, const int npy,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells,
            const int kcells, const int ijcells)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        const int igc = istart;
        const int jgc = jstart;

        if (mpiidx == 0 && mpiidy == 0)
        {
            for (int k=kstart; k<kend; ++k)
            {
                const int ijk0 = istart + jstart*jj + k*kk;
                const TF dfdi = fld[ijk0] - fld[ijk0-ii];
                const TF dfdj = fld[ijk0] - fld[ijk0-jj];

                for (int dj=1; dj<jgc+1; ++dj)
                    for (int di=1; di<igc+1; ++di)
                    {
                        const int i = istart-di;
                        const int j = jstart-dj;
                        const int ijk = i + j*jj + k*kk;

                        fld[ijk] = fld[ijk0] - di*dfdi - dj*dfdj;
                    }
            }

            for (int k=0; k<kstart; ++k)
                for (int j=jstart-1; j>=0; --j)
                    for (int i=istart-1; i>=0; --i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        const int ijks = i+j*jj + kstart*kk;
                        fld[ijk] = fld[ijks];
                    }

            for (int k=kend; k<kcells; ++k)
                for (int j=jstart-1; j>=0; --j)
                    for (int i=istart-1; i>=0; --i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        const int ijke = i+j*jj + (kend-1)*kk;
                        fld[ijk] = fld[ijke];
                    }
        }

        if (mpiidx == npx-1 && mpiidy == 0)
        {
            // South-east corner
            for (int k=kstart; k<kend; ++k)
            {
                const int ijk0 = (iend-1) + jstart*jj + k*kk;
                const TF dfdi = fld[ijk0+ii] - fld[ijk0];
                const TF dfdj = fld[ijk0] - fld[ijk0-jj];

                for (int dj=1; dj<jgc+1; ++dj)
                    for (int di=1; di<igc+1; ++di)
                    {
                        const int i = (iend-1)+di;
                        const int j = jstart-dj;
                        const int ijk = i + j*jj + k*kk;

                        fld[ijk] = fld[ijk0] + di*dfdi - dj*dfdj;
                    }
            }

            for (int k=0; k<kstart; ++k)
                for (int j=jstart-1; j>=0; --j)
                    for (int i=iend; i<icells; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        const int ijks = i+j*jj + kstart*kk;
                        fld[ijk] = fld[ijks];
                    }

            for (int k=kend; k<kcells; ++k)
                for (int j=jstart-1; j>=0; --j)
                    for (int i=iend; i<icells; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        const int ijke = i+j*jj + (kend-1)*kk;
                        fld[ijk] = fld[ijke];
                    }
        }

        if (mpiidx == 0 && mpiidy == npy-1)
        {
            // North-west corner
            for (int k=kstart; k<kend; ++k)
            {
                const int ijk0 = istart + (jend-1)*jj + k*kk;
                const TF dfdi = fld[ijk0] - fld[ijk0-ii];
                const TF dfdj = fld[ijk0+jj] - fld[ijk0];

                for (int dj=1; dj<jgc+1; ++dj)
                    for (int di=1; di<igc+1; ++di)
                    {
                        const int i = istart-di;
                        const int j = (jend-1)+dj;
                        const int ijk = i + j*jj + k*kk;

                        fld[ijk] = fld[ijk0] - di*dfdi + dj*dfdj;
                    }
            }

            for (int k=0; k<kstart; ++k)
                for (int j=jend; j<jcells; ++j)
                    for (int i=istart-1; i>=0; --i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        const int ijks = i+j*jj + kstart*kk;
                        fld[ijk] = fld[ijks];
                    }

            for (int k=kend; k<kcells; ++k)
                for (int j=jend; j<jcells; ++j)
                    for (int i=istart-1; i>=0; --i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        const int ijke = i+j*jj + (kend-1)*kk;
                        fld[ijk] = fld[ijke];
                    }
        }

        if (mpiidx == npx-1 && mpiidy == npy-1)
        {
            // North-east corner
            for (int k=kstart; k<kend; ++k)
            {
                const int ijk0 = (iend-1) + (jend-1)*jj + k*kk;
                const TF dfdi = fld[ijk0+ii] - fld[ijk0];
                const TF dfdj = fld[ijk0+jj] - fld[ijk0];

                for (int dj=1; dj<jgc+1; ++dj)
                    for (int di=1; di<igc+1; ++di)
                    {
                        const int i = (iend-1)+di;
                        const int j = (jend-1)+dj;
                        const int ijk = i + j*jj + k*kk;

                        fld[ijk] = fld[ijk0] + di*dfdi + dj*dfdj;
                    }
            }

            for (int k=0; k<kstart; ++k)
                for (int j=jend; j<jcells; ++j)
                    for (int i=iend; i<icells; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        const int ijks = i+j*jj + kstart*kk;
                        fld[ijk] = fld[ijks];
                    }

            for (int k=kend; k<kcells; ++k)
                for (int j=jend; j<jcells; ++j)
                    for (int i=iend; i<icells; ++i)
                    {
                        const int ijk = i + j*jj + k*kk;
                        const int ijke = i+j*jj + (kend-1)*kk;
                        fld[ijk] = fld[ijke];
                    }
        }
    }

    template<typename TF>
    void interpolate_lbc_kernel(
            TF* const restrict fld,
            const TF* const restrict fld_in,
            const int stride,
            const int t0,
            const TF f0)
    {
        for (int n=0; n<stride; ++n)
        {
            const int n0 = n + t0*stride;
            const int n1 = n + (t0+1)*stride;

            fld[n] = f0 * fld_in[n0] + (TF(1)-f0) * fld_in[n1];
        }
    }


    //template<typename TF>
    //void calc_div_h(
    //        TF* const restrict div,
    //        const TF* const restrict lbc_u,
    //        const TF* const restrict lbc_d,
    //        const TF* const restrict rhoref,
    //        const TF* const restrict dz,
    //        const TF dx_or_dy,
    //        const int ntime,
    //        const int ntot, const int ktot,
    //        const int kgc)
    //{
    //    for (int t=0; t<ntime; ++t)
    //    {
    //        div[t] = TF(0);
    //        for (int k=0; k<ktot; ++k)
    //            for (int n=0; n<ntot; ++n)
    //            {
    //                const int nk = n + k*ntot + t*ntot*ktot;

    //                div[t] += rhoref[k+kgc] * dx_or_dy * dz[k+kgc] * (lbc_u[nk] - lbc_d[nk]);
    //            }
    //    }
    //}


    template<typename TF, Lbc_location location>
    void perturb_boundary_kernel(
            TF* const restrict tend,
            const TF* const restrict fld,
            const TF* const restrict lbc,
            const TF amplitude, const TF dt,
            const int width, const int block_size,
            const int mpiidx, const int npx,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells,
            const int ijcells)
    {
        const int n_blocks = width / block_size;

        const TF rand_max_i = TF(1) / RAND_MAX;
        const TF dt_i = TF(1) / dt;

        if (location == Lbc_location::West || location == Lbc_location::East)
        {
            const int i0 = (location == Lbc_location::West) ? istart : iend-width;
            const int n_blocks_ew = (jend-jstart)  / block_size;

            for (int k=kstart; k<kend; ++k)
                for (int bj=0; bj<n_blocks_ew; ++bj)
                    for (int bi=0; bi<n_blocks; ++bi)
                    {
                        const TF pert = (TF(std::rand()) * rand_max_i - TF(0.5)) * amplitude;

                        for (int dj=0; dj<block_size; ++dj)
                            for (int di=0; di<block_size; ++di)
                            {
                                const int i = i0 + bi*block_size + di;
                                const int j = jstart + bj*block_size + dj;
                                const int ijk = i + j*icells + k*ijcells;
                                const int jk = j + k*jcells;

                                tend[ijk] -= ((fld[ijk] - (lbc[jk] + pert)) * dt_i + tend[ijk]);
                            }
                    }
        }

        if (location == Lbc_location::South || location == Lbc_location::North)
        {
            const int j0 = (location == Lbc_location::South) ? jstart : jend-width;
            const int istart_loc = (mpiidx == 0) ? istart + width : istart;
            const int iend_loc = (mpiidx == npx-1) ? iend - width : iend;
            const int n_blocks_ns = (iend_loc-istart_loc)  / block_size;

            for (int k=kstart; k<kend; ++k)
                for (int bj=0; bj<n_blocks; ++bj)
                    for (int bi=0; bi<n_blocks_ns; ++bi)
                    {
                        const TF pert = (TF(std::rand()) * rand_max_i - TF(0.5)) * amplitude;

                        for (int dj=0; dj<block_size; ++dj)
                            for (int di=0; di<block_size; ++di)
                            {
                                const int i = istart_loc + bi*block_size + di;
                                const int j = j0 + bj*block_size + dj;
                                const int ijk = i + j*icells + k*ijcells;
                                const int ik = i + k*icells;

                                tend[ijk] -= ((fld[ijk] - (lbc[ik] + pert)) * dt_i + tend[ijk]);
                            }
                    }
        }
    }
}



namespace
{
    template<typename TF>
    void set_lbc_gcs(
            TF* const restrict fld,
            const TF* const restrict lbc,
            const int ngc,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells,
            const int kcells,
            Lbc_location location)
    {
        const int jstride_out = icells;
        const int kstride_out = icells * jcells;

        if (location == Lbc_location::West)
        {
            const int jstride_w = ngc;
            const int kstride_w = jstride_w * jcells;

            for (int k=kstart; k<kend; k++)
                for (int j=0; j<jcells; j++)
                    for (int i=0; i<ngc; i++)
                    {
                        const int ijk_in = i + j*jstride_w + k*kstride_w;
                        const int ijk_out = i + j*jstride_out + k*kstride_out;
                        fld[ijk_out] = lbc[ijk_in];
                    }
        }

        if (location == Lbc_location::East)
        {
            const int jstride_e = ngc;
            const int kstride_e = jstride_e * jcells;

            for (int k=kstart; k<kend; k++)
                for (int j=0; j<jcells; j++)
                    for (int i=0; i<ngc; i++)
                    {
                        const int ijk_in = i + j*jstride_e + k*kstride_e;
                        const int ijk_out = (i+iend) + j*jstride_out + k*kstride_out;
                        fld[ijk_out] = lbc[ijk_in];
                    }
        }

        if (location == Lbc_location::South)
        {
            const int jstride_s = icells;
            const int kstride_s = jstride_s * ngc;

            for (int k=kstart; k<kend; k++)
                for (int j=0; j<ngc; j++)
                    for (int i=istart; i<iend; i++)
                    {
                        const int ijk_in = i + j*jstride_s + k*kstride_s;
                        const int ijk_out = i + j*jstride_out + k*kstride_out;
                        fld[ijk_out] = lbc[ijk_in];
                    }
        }

        if (location == Lbc_location::North)
        {
            const int jstride_n= icells;
            const int kstride_n = jstride_n * ngc;

            for (int k=kstart; k<kend; k++)
                for (int j=0; j<ngc; j++)
                    for (int i=istart; i<iend; i++)
                    {
                        const int ijk_in = i + j*jstride_n + k*kstride_n;
                        const int ijk_out = i + (j+jend)*jstride_out + k*kstride_out;
                        fld[ijk_out] = lbc[ijk_in];
                    }
        }
    };


    template<typename TF>
    void calc_div_x(
            TF* const restrict div,
            const TF* const restrict lbc_u_west,
            const TF* const restrict lbc_u_east,
            const TF* const restrict rhoref,
            const TF* const restrict dz,
            const TF dy,
            const int ntime,
            const int ngc, const int kgc,
            const int jstart, const int jend,
            const int ktot,
            const int jcells)
    {
        const int jstride_w = ngc+1;
        const int kstride_w = jstride_w * jcells;
        const int tstride_w = kstride_w * ktot;

        const int jstride_e = ngc;
        const int kstride_e = jstride_e * jcells;
        const int tstride_e = kstride_e * ktot;

        const int iw = ngc;
        const int ie = 0;

        for (int t=0; t<ntime; ++t)
        {
            div[t] = 0;
            for (int k=0; k<ktot; ++k)
                for (int j=jstart; j<jend; ++j)
                {
                    const int ijk_w = iw + j*jstride_w + k*kstride_w + t*tstride_w;
                    const int ijk_e = ie + j*jstride_e + k*kstride_e + t*tstride_e;

                    const TF v1 = lbc_u_east[ijk_e];
                    const TF v2 = lbc_u_west[ijk_w];

                    div[t] += rhoref[k+kgc] * dy * dz[k+kgc] * (lbc_u_east[ijk_e] - lbc_u_west[ijk_w]);
                }
        }
    }

    template<typename TF>
    void calc_div_y(
            TF* const restrict div,
            const TF* const restrict lbc_v_south,
            const TF* const restrict lbc_v_north,
            const TF* const restrict rhoref,
            const TF* const restrict dz,
            const TF dx,
            const int ntime,
            const int ngc, const int kgc,
            const int istart, const int iend,
            const int ktot,
            const int icells)
    {
        const int jstride_s = icells;
        const int kstride_s = jstride_s * (ngc+1);
        const int tstride_s = kstride_s * ktot;

        const int jstride_n = icells;
        const int kstride_n = jstride_n * ngc;
        const int tstride_n = kstride_n * ktot;

        const int js = ngc;
        const int jn = 0;

        for (int t=0; t<ntime; ++t)
        {
            div[t] = 0;
            for (int k=0; k<ktot; ++k)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk_s = i + js*jstride_s + k*kstride_s + t*tstride_s;
                    const int ijk_n = i + jn*jstride_n + k*kstride_n + t*tstride_n;

                    const TF v1 = lbc_v_north[ijk_n];
                    const TF v2 = lbc_v_south[ijk_s];

                    div[t] += rhoref[k+kgc] * dx * dz[k+kgc] * (lbc_v_north[ijk_n] - lbc_v_south[ijk_s]);
                }
        }
    }
}


template<typename TF>
Boundary_lateral<TF>::Boundary_lateral(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
        master(masterin), grid(gridin), fields(fieldsin), field3d_io(masterin, gridin)
{
    sw_inoutflow = inputin.get_item<bool>("boundary", "sw_inoutflow", "", false);

    if (sw_inoutflow)
    {
        sw_inoutflow_uv = inputin.get_item<bool>("boundary", "sw_inoutflow_uv", "", true);
        sw_inoutflow_w = inputin.get_item<bool>("boundary", "sw_inoutflow_w", "", false);
        sw_neumann_w = inputin.get_item<bool>("boundary", "sw_neumann_w", "", true);
        sw_wtop_2d = inputin.get_item<bool>("boundary", "sw_wtop_2d", "", false);
        inoutflow_s = inputin.get_list<std::string>("boundary", "inoutflow_slist", "", std::vector<std::string>());

        // Check....
        if (sw_inoutflow_w && sw_neumann_w)
            throw std::runtime_error("Cant have both \"sw_inoutflow_w\" and \"sw_neumann_w\" = true!");

        sw_timedep = inputin.get_item<bool>("boundary", "sw_timedep", "", false);
        if (sw_wtop_2d && sw_timedep)
            wtop_2d_loadtime = inputin.get_item<int>("boundary", "wtop_2d_loadtime", "");

        // Lateral sponge / diffusion layer.
        sw_sponge = inputin.get_item<bool>("boundary", "sw_sponge", "", false);
        if (sw_sponge)
        {
            n_sponge = inputin.get_item<int>("boundary", "n_sponge", "", 5);
            tau_nudge = inputin.get_item<TF>("boundary", "tau_nudge", "", 60);
            w_diff = inputin.get_item<TF>("boundary", "w_diff", "", 0.0033);
        }

        //// Inflow perturbations
        //sw_perturb = inputin.get_item<bool>("boundary", "sw_perturb", "", false);
        //if (sw_perturb)
        //{
        //    perturb_list = inputin.get_list<std::string>("boundary", "perturb_list", "", std::vector<std::string>());
        //    perturb_width = inputin.get_item<int>("boundary", "perturb_width", "", 4);
        //    perturb_block = inputin.get_item<int>("boundary", "perturb_block", "", 2);
        //    perturb_seed = inputin.get_item<int>("boundary", "perturb_seed", "", 0);

        //    for (auto& fld : perturb_list)
        //    {
        //        const TF ampl = inputin.get_item<TF>("boundary", "perturb_ampl", fld);
        //        perturb_ampl.emplace(fld, ampl);
        //    }
        //}

        //// Turbulence recycling.
        //recycle_list = inputin.get_list<std::string>("boundary", "recycle_list", "", std::vector<std::string>());
        //if (recycle_list.size() > 0)
        //{
        //    tau_recycle = inputin.get_item<TF>("boundary", "tau_recycle", "");
        //    recycle_offset = inputin.get_item<int>("boundary", "recycle_offset", "");
        //}
    }
}

template <typename TF>
Boundary_lateral<TF>::~Boundary_lateral()
{
}

template <typename TF>
void Boundary_lateral<TF>::init()
{
    if (!sw_inoutflow)
        return;

    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    // Checks!
    //if (recycle_list.size() > 0)
    //{
    //    if (!sw_sponge)
    //        throw std::runtime_error("Turbulence recycling only works combined with sw_sponge=1");

    //    if (n_sponge+recycle_offset > gd.imax+gd.igc)
    //        throw std::runtime_error("Turbulence recycling offset too large for domain decomposition in x-direction");
    //    if (n_sponge+recycle_offset > gd.jmax+gd.jgc)
    //        throw std::runtime_error("Turbulence recycling offset too large for domain decomposition in y-direction");
    //}

    auto add_lbc = [&](const std::string& name)
    {
        const int igc_pad = (name == "u") ? gd.igc+1 : gd.igc;
        const int jgc_pad = (name == "v") ? gd.jgc+1 : gd.jgc;

        if (md.mpicoordx == 0)
            lbc_w.emplace(name, std::vector<TF>(igc_pad*gd.kcells*gd.jcells));
        if (md.mpicoordx == md.npx-1)
            lbc_e.emplace(name, std::vector<TF>(gd.igc*gd.kcells*gd.jcells));
        if (md.mpicoordy == 0)
            lbc_s.emplace(name, std::vector<TF>(gd.icells*jgc_pad*gd.kcells));
        if (md.mpicoordy == md.npy-1)
            lbc_n.emplace(name, std::vector<TF>(gd.icells*gd.jgc*gd.kcells));
    };

    if (sw_inoutflow_uv)
    {
        add_lbc("u");
        add_lbc("v");
    }

    if (sw_inoutflow_w)
        add_lbc("w");

    for (auto& fld : inoutflow_s)
        add_lbc(fld);

    //// Make sure every MPI task has different seed.
    //if (sw_perturb)
    //{
    //    perturb_seed += master.get_mpiid();
    //    std::srand(perturb_seed);
    //}
}

template <typename TF>
void Boundary_lateral<TF>::create(
        Input& inputin,
        Timeloop<TF>& timeloop,
        const std::string& sim_name)
{
    if (!sw_inoutflow)
        return;

    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    // Read NetCDF file with boundary data.
    Netcdf_file input_nc = Netcdf_file(master, sim_name + "_lbc_input.nc", Netcdf_mode::Read);

    const int ntime = input_nc.get_dimension_size("time");
    time_in = input_nc.get_variable<TF>("time", {ntime});

    TF* rhoref = fields.rhoref.data();

    // Domain total divergence in u and v direction.
    if (sw_inoutflow_uv)
    {
        div_u.resize(ntime);
        div_v.resize(ntime);
        w_top.resize(gd.ijcells);
    }

    auto dump_vector = [&](
            std::vector<TF>& fld,
            const std::string& name)
    {
        std::string name_out = name + "." + std::to_string(md.mpicoordx) + "." + std::to_string(md.mpicoordy) + ".bin";

        FILE *pFile;
        pFile = fopen(name_out.c_str(), "wb");

        if (pFile == NULL)
            throw std::runtime_error("Opening raw dump field failed.");

        fwrite(fld.data(), sizeof(TF), fld.size(), pFile);
        fclose(pFile);
    };

    auto copy_boundary = [&](
            std::map<std::string, std::vector<TF>>& map_out,
            const std::vector<TF>& fld_in,
            const int isize_in, const int jsize_in,
            const int isize_out, const int jsize_out,
            const int istart_in, const int jstart_in,
            const std::string& name, const std::string& loc)
    {
        const int size_in = ntime * gd.ktot * jsize_in * isize_in;
        const int size_out = ntime * gd.kcells * jsize_out * isize_out;

        std::vector<TF> fld_out = std::vector<TF>(size_out);

        const int jstride_in = isize_in;
        const int kstride_in = jstride_in * jsize_in;
        const int tstride_in = kstride_in * gd.ktot;

        const int jstride_out = isize_out;
        const int kstride_out = jstride_out * jsize_out;
        const int tstride_out = kstride_out * gd.kcells;

        for (int t=0; t<ntime; t++)
            for (int k=0; k<gd.ktot; k++)
                for (int j=0; j<jsize_out; j++)
                    for (int i=0; i<isize_out; i++)
                    {
                        const int ijk_in = i+istart_in + (j+jstart_in)*jstride_in + k*kstride_in + t*tstride_in;
                        const int ijk_out = i + j*jstride_out + (k+gd.kstart)*kstride_out + t*tstride_out;

                        fld_out[ijk_out] = fld_in[ijk_in];
                    }

        map_out.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(name),
                std::forward_as_tuple(std::move(fld_out)));
    };

    auto copy_boundaries = [&](const std::string& name)
    {
        // `u` at west boundary and `v` at south boundary also contain `u` at `istart`
        // and `v` at `jstart`, and are therefore larger.
        const int igc_pad = (name == "u") ? gd.igc+1 : gd.igc;
        const int jgc_pad = (name == "v") ? gd.jgc+1 : gd.jgc;

        // Read full boundaries for entire domain.
        std::vector<TF> lbc_w_full = input_nc.get_variable<TF>(name + "_west",  {ntime, gd.ktot, gd.jtot+2*gd.jgc, igc_pad});
        std::vector<TF> lbc_e_full = input_nc.get_variable<TF>(name + "_east",  {ntime, gd.ktot, gd.jtot+2*gd.jgc, gd.igc});
        std::vector<TF> lbc_s_full = input_nc.get_variable<TF>(name + "_south", {ntime, gd.ktot, jgc_pad, gd.itot+2*gd.igc});
        std::vector<TF> lbc_n_full = input_nc.get_variable<TF>(name + "_north", {ntime, gd.ktot, gd.jgc, gd.itot+2*gd.igc});

        if (md.mpicoordx == 0)
            copy_boundary(
                    lbc_w_in, lbc_w_full,
                    igc_pad, gd.jtot+2*gd.jgc,
                    igc_pad, gd.jcells,
                    0, md.mpicoordy*gd.jmax,
                    name, "west");

        if (md.mpicoordx == md.npx-1)
            copy_boundary(
                    lbc_e_in, lbc_e_full,
                    gd.igc, gd.jtot+2*gd.jgc,
                    gd.igc, gd.jcells,
                    0, md.mpicoordy*gd.jmax,
                    name, "east");

        if (md.mpicoordy == 0)
            copy_boundary(
                    lbc_s_in, lbc_s_full,
                    gd.itot+2*gd.igc, jgc_pad,
                    gd.icells, jgc_pad,
                    md.mpicoordx*gd.imax, 0,
                    name, "south");

        if (md.mpicoordy == md.npy-1)
            copy_boundary(
                    lbc_n_in, lbc_n_full,
                    gd.itot+2*gd.igc, gd.jgc,
                    gd.icells, gd.jgc,
                    md.mpicoordx*gd.imax, 0,
                    name, "north");

        // Calculate domain total mass imbalance in kg s-1.
        if (name == "u")
            calc_div_x(
                    div_u.data(),
                    lbc_w_full.data(),
                    lbc_e_full.data(),
                    fields.rhoref.data(),
                    gd.dz.data(),
                    gd.dy,
                    ntime,
                    gd.igc, gd.kgc,
                    gd.jgc, gd.jtot+gd.jgc,
                    gd.ktot, gd.jtot+(2*gd.jgc));
        else if (name == "v")
            calc_div_y(
                    div_v.data(),
                    lbc_s_full.data(),
                    lbc_n_full.data(),
                    fields.rhoref.data(),
                    gd.dz.data(),
                    gd.dx,
                    ntime,
                    gd.jgc, gd.kgc,
                    gd.igc, gd.itot+gd.igc,
                    gd.ktot, gd.itot+(2*gd.igc));

        if (!sw_timedep)
        {
            if (md.mpicoordx == 0)
            {
                for (int n=0; n<igc_pad*gd.jcells*gd.kcells; ++n)
                    lbc_w.at(name)[n] = lbc_w_in.at(name)[n];
            }

            if (md.mpicoordx == md.npx-1)
            {
                for (int n=0; n<gd.igc*gd.jcells*gd.kcells; ++n)
                    lbc_e.at(name)[n] = lbc_e_in.at(name)[n];
            }

            if (md.mpicoordy == 0)
            {
                for (int n=0; n<gd.icells*jgc_pad*gd.kcells; ++n)
                    lbc_s.at(name)[n] = lbc_s_in.at(name)[n];
            }

            if (md.mpicoordy == md.npy-1)
            {
                for (int n=0; n<gd.icells*gd.jgc*gd.kcells; ++n)
                    lbc_n.at(name)[n] = lbc_n_in.at(name)[n];
            }
        }
    };

    if (sw_inoutflow_uv)
    {
        copy_boundaries("u");
        copy_boundaries("v");
    }

    if (sw_inoutflow_w)
        copy_boundaries("w");

    for (auto& fld : inoutflow_s)
        copy_boundaries(fld);

    // Calculate domain mean vertical velocity.
    if (sw_inoutflow_uv)
    {
        if (sw_wtop_2d)
        {
            // Read constant or time varying w_top fields.
            if (sw_timedep)
            {
                // Find previous and next times.
                const double time = timeloop.get_time();
                const double ifactor = timeloop.get_ifactor();
                unsigned long iiotimeprec = timeloop.get_iiotimeprec();

                // Read first two input times
                itime_w_top_prev = ifactor * int(time/wtop_2d_loadtime) * wtop_2d_loadtime;
                itime_w_top_next = itime_w_top_prev + wtop_2d_loadtime*ifactor;

                // IO time accounting for iotimeprec
                const unsigned long iotime0 = int(itime_w_top_prev / iiotimeprec);
                const unsigned long iotime1 = int(itime_w_top_next / iiotimeprec);

                w_top_prev.resize(gd.ijcells);
                w_top_next.resize(gd.ijcells);

                // Read the first two w_top fields.
                read_xy_slice(w_top_prev, "w_top", iotime0);
                read_xy_slice(w_top_next, "w_top", iotime1);
            }
            else
                read_xy_slice(w_top, "w_top", 0);
        }
        else
        {
            // Calculate domain mean `w_top`.
            w_top_in.resize(ntime);

            for (int t=0; t<ntime; ++t)
            {
                // w_top is the total mass in/outflow at the top divided by the total area and the local density.
                w_top_in[t] = -(div_u[t] + div_v[t]) / (fields.rhorefh[gd.kend] * gd.xsize * gd.ysize);

                std::string message=
                        "- div(u) = " + std::to_string(div_u[t])
                      + ", div(v) = " + std::to_string(div_v[t])
                      + ", w_top = " + std::to_string(w_top_in[t])
                      + " m/s @ t=" + std::to_string(time_in[t]) + " sec.";
                master.print_message(message);
            }

            if (!sw_timedep)
                std::fill(w_top.begin(), w_top.end(), w_top_in[0]);
        }
    }

    //if (sw_perturb)
    //{
    //    const TF perturb_zmax = inputin.get_item<TF>("boundary", "perturb_zmax", "");

    //    for (int k=gd.kstart; k<gd.kend; ++k)
    //        if (gd.z[k] < perturb_zmax && gd.z[k+1] >= perturb_zmax)
    //        {
    //            perturb_kend = k+1;
    //            break;
    //        }
    //}
}



template <typename TF>
void Boundary_lateral<TF>::set_ghost_cells(Timeloop<TF>& timeloop)
{
    if (!sw_inoutflow)
        return;

    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    auto dump_vector = [&](
            std::vector<TF>& fld,
            const std::string& name)
    {
        std::string name_out = name + "." + std::to_string(md.mpicoordx) + "." + std::to_string(md.mpicoordy) + ".bin";

        FILE *pFile;
        pFile = fopen(name_out.c_str(), "wb");

        if (pFile == NULL)
            throw std::runtime_error("Opening raw dump field failed.");

        fwrite(fld.data(), sizeof(TF), fld.size(), pFile);
        fclose(pFile);
    };


    auto set_lbc_gcs_wrapper = [&](
            const std::string& fld,
            std::vector<TF>& lbc,
            const Lbc_location location)
    {
        int ngc = gd.igc;
        if (fld == "u" && location == Lbc_location::West)
            ngc += 1;
        if (fld == "v" && location == Lbc_location::South)
            ngc += 1;

        set_lbc_gcs(
                fields.ap.at(fld)->fld.data(),
                lbc.data(),
                ngc,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.jcells,
                gd.kcells,
                location);
    };

    auto set_gcs = [&](const std::string& fld)
    {
        if (md.mpicoordx == 0)
            set_lbc_gcs_wrapper(fld, lbc_w.at(fld), Lbc_location::West);
        if (md.mpicoordx == md.npx-1)
            set_lbc_gcs_wrapper(fld, lbc_e.at(fld), Lbc_location::East);
        if (md.mpicoordy == 0)
            set_lbc_gcs_wrapper(fld, lbc_s.at(fld), Lbc_location::South);
        if (md.mpicoordy == md.npy-1)
            set_lbc_gcs_wrapper(fld, lbc_n.at(fld), Lbc_location::North);
    };

    if (sw_inoutflow_uv)
    {
        set_gcs("u");
        set_gcs("v");
    }

    if (sw_inoutflow_w)
        set_gcs("w");

    for (auto& fld : inoutflow_s)
        set_gcs(fld);


    // Set vertical velocity at domain top.
    if (sw_inoutflow_uv)
    {
        const int k = gd.kend;
        for (int j=gd.jstart; j<gd.jend; ++j)
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*gd.icells + k*gd.ijcells;
                const int ij = i + j*gd.icells;
                fields.mp.at("w")->fld[ijk] = w_top[ij];
            }
    }

    if (sw_neumann_w)
    {
        // Here, we enfore a Neumann BC of 0 over the boundaries. This works if the large scale w is approximately
        // constant in the horizontal plane. Note that w must be derived from u and v if the large-scale field is to
        // be divergence free. If there is a horizontal gradient in w_top, then it is probably better to extrapolate that
        // gradient into the ghost cells.
        auto set_ghost_cell_w_wrapper = [&]<Lbc_location location>()
        {
            set_ghost_cell_kernel_w<TF, location>(
                    fields.mp.at("w")->fld.data(),
                    gd.istart, gd.iend, gd.igc,
                    gd.jstart, gd.jend, gd.jgc,
                    gd.kstart, gd.kend+1,
                    gd.icells, gd.jcells, gd.kcells,
                    gd.ijcells);
        };

        auto set_corner_ghost_cell_wrapper = [&](
                std::vector<TF>& fld,
                const int kend)
        {
            set_corner_ghost_cell_kernel(
                    fld.data(),
                    md.mpicoordx, md.mpicoordy,
                    md.npx, md.npy,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, kend,
                    gd.icells, gd.jcells,
                    gd.kcells, gd.ijcells);
        };

        if (md.mpicoordx == 0)
            set_ghost_cell_w_wrapper.template operator()<Lbc_location::West>();
        if (md.mpicoordx == md.npx-1)
            set_ghost_cell_w_wrapper.template operator()<Lbc_location::East>();
        if (md.mpicoordy == 0)
            set_ghost_cell_w_wrapper.template operator()<Lbc_location::South>();
        if (md.mpicoordy == md.npy-1)
            set_ghost_cell_w_wrapper.template operator()<Lbc_location::North>();

        set_corner_ghost_cell_wrapper(fields.mp.at("w")->fld, gd.kend+1);
    }

    //auto set_ghost_cell_s_wrapper = [&]<Lbc_location location>(
    //        std::map<std::string, std::vector<TF>>& lbc_map,
    //        const std::string& name)
    //{
    //    set_ghost_cell_kernel_s<TF, location>(
    //            fields.ap.at(name)->fld.data(),
    //            lbc_map.at(name).data(),
    //            gd.istart, gd.iend, gd.igc,
    //            gd.jstart, gd.jend, gd.jgc,
    //            gd.kstart, gd.kend,
    //            gd.icells, gd.jcells, gd.kcells,
    //            gd.ijcells);
    //};

    auto sponge_layer_wrapper = [&]<Lbc_location location, bool sw_recycle>(
            std::map<std::string, std::vector<TF>>& lbc_map,
            const std::string& name)
    {
        if (!sw_sponge)
            return;

        lateral_sponge_kernel_s<TF, location, sw_recycle>(
                fields.at.at(name)->fld.data(),
                fields.ap.at(name)->fld.data(),
                lbc_map.at(name).data(),
                tau_nudge,
                w_diff,
                n_sponge,
                tau_recycle,
                recycle_offset,
                md.npx, md.npy,
                md.mpicoordx, md.mpicoordy,
                gd.igc, gd.jgc,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.jcells,
                gd.ijcells);
    };

    //auto set_ghost_cell_u_wrapper = [&]<Lbc_location location>(
    //        std::map<std::string, std::vector<TF>>& lbc_map)
    //{
    //    set_ghost_cell_kernel_u<TF, location>(
    //            fields.mp.at("u")->fld.data(),
    //            lbc_map.at("u").data(),
    //            gd.istart, gd.iend, gd.igc,
    //            gd.jstart, gd.jend,
    //            gd.kstart, gd.kend,
    //            gd.icells, gd.jcells,
    //            gd.ijcells);
    //};

    auto sponge_layer_u_wrapper = [&]<Lbc_location location>(
            std::map<std::string, std::vector<TF>>& lbc_map)
    {
        if (!sw_sponge)
            return;

        lateral_sponge_kernel_u<TF, location>(
                fields.mt.at("u")->fld.data(),
                fields.mp.at("u")->fld.data(),
                lbc_map.at("u").data(),
                tau_nudge,
                w_diff,
                n_sponge,
                md.npy,
                md.mpicoordy,
                gd.igc,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.jcells,
                gd.ijcells);
    };

    //auto set_ghost_cell_v_wrapper = [&]<Lbc_location location>(
    //        std::map<std::string, std::vector<TF>>& lbc_map)
    //{
    //    set_ghost_cell_kernel_v<TF, location>(
    //            fields.mp.at("v")->fld.data(),
    //            lbc_map.at("v").data(),
    //            gd.istart, gd.iend,
    //            gd.jstart, gd.jend, gd.jgc,
    //            gd.kstart, gd.kend,
    //            gd.icells, gd.jcells,
    //            gd.ijcells);
    //};

    auto sponge_layer_v_wrapper = [&]<Lbc_location location>(
            std::map<std::string, std::vector<TF>>& lbc_map)
    {
        if (!sw_sponge)
            return;

        lateral_sponge_kernel_v<TF, location>(
                fields.mt.at("v")->fld.data(),
                fields.mp.at("v")->fld.data(),
                lbc_map.at("v").data(),
                tau_nudge,
                w_diff,
                n_sponge,
                md.npx,
                md.mpicoordx,
                gd.jgc,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.jcells,
                gd.ijcells);
    };

    if (sw_inoutflow_uv)
    {
        if (md.mpicoordy == 0)
        {
            sponge_layer_wrapper.template operator()<Lbc_location::South, false>(lbc_s, "u");
            sponge_layer_v_wrapper.template operator()<Lbc_location::South>(lbc_s);
        }
        if (md.mpicoordy == md.npy-1)
        {
            sponge_layer_wrapper.template operator()<Lbc_location::North, false>(lbc_n, "u");
            sponge_layer_v_wrapper.template operator()<Lbc_location::North>(lbc_n);
        }
        if (md.mpicoordx == 0)
        {
            sponge_layer_u_wrapper.template operator()<Lbc_location::West>(lbc_w);
            sponge_layer_wrapper.template operator()<Lbc_location::West, false>(lbc_w, "v");
        }
        if (md.mpicoordx == md.npx-1)
        {
            sponge_layer_u_wrapper.template operator()<Lbc_location::East>(lbc_e);
            sponge_layer_wrapper.template operator()<Lbc_location::East, false>(lbc_e, "v");
        }
    }

    //if (sw_inoutflow_u)
    //{
    //    if (md.mpicoordy == 0)
    //    {
    //        set_ghost_cell_s_wrapper.template operator()<Lbc_location::South>(lbc_s, "u");
    //        sponge_layer_wrapper.template operator()<Lbc_location::South, false>(lbc_s, "u");
    //    }
    //    if (md.mpicoordy == md.npy-1)
    //    {
    //        set_ghost_cell_s_wrapper.template operator()<Lbc_location::North>(lbc_n, "u");
    //        sponge_layer_wrapper.template operator()<Lbc_location::North, false>(lbc_n, "u");
    //    }
    //    if (md.mpicoordx == 0)
    //    {
    //        set_ghost_cell_u_wrapper.template operator()<Lbc_location::West>(lbc_w);
    //        sponge_layer_u_wrapper.template operator()<Lbc_location::West>(lbc_w);
    //    }
    //    if (md.mpicoordx == md.npx-1)
    //    {
    //        set_ghost_cell_u_wrapper.template operator()<Lbc_location::East>(lbc_e);
    //        sponge_layer_u_wrapper.template operator()<Lbc_location::East>(lbc_e);
    //    }

    //    set_corner_ghost_cell_wrapper(fields.mp.at("u")->fld, gd.kend);
    //}

    //if (sw_inoutflow_v)
    //{
    //    if (md.mpicoordx == 0)
    //    {
    //        set_ghost_cell_s_wrapper.template operator()<Lbc_location::West>(lbc_w, "v");
    //        sponge_layer_wrapper.template operator()<Lbc_location::West, false>(lbc_w, "v");
    //    }
    //    if (md.mpicoordx == md.npx-1)
    //    {
    //        set_ghost_cell_s_wrapper.template operator()<Lbc_location::East>(lbc_e, "v");
    //        sponge_layer_wrapper.template operator()<Lbc_location::East, false>(lbc_e, "v");
    //    }
    //    if (md.mpicoordy == 0)
    //    {
    //        set_ghost_cell_v_wrapper.template operator()<Lbc_location::South>(lbc_s);
    //        sponge_layer_v_wrapper.template operator()<Lbc_location::South>(lbc_s);
    //    }
    //    if (md.mpicoordy == md.npy-1)
    //    {
    //        set_ghost_cell_v_wrapper.template operator()<Lbc_location::North>(lbc_n);
    //        sponge_layer_v_wrapper.template operator()<Lbc_location::North>(lbc_n);
    //    }

    //    set_corner_ghost_cell_wrapper(fields.mp.at("v")->fld, gd.kend);
    //}

    //// Here, we enfore a Neumann BC of 0 over the boundaries. This works if the large scale w is approximately
    //// constant in the horizontal plane. Note that w must be derived from u and v if the large-scale field is to
    //// be divergence free. If there is a horizontal gradient in w_top, then it is probably better to extrapolate that
    //// gradient into the ghost cells.
    //if (sw_inoutflow_w)
    //{
    //    if (md.mpicoordx == 0)
    //    {
    //        set_ghost_cell_w_wrapper.template operator()<Lbc_location::West>();
    //        // sponge_layer_wrapper.template operator()<Lbc_location::West>(lbc_w, "v");
    //    }
    //    if (md.mpicoordx == md.npx-1)
    //    {
    //        set_ghost_cell_w_wrapper.template operator()<Lbc_location::East>();
    //        // sponge_layer_wrapper.template operator()<Lbc_location::East>(lbc_e, "v");
    //    }
    //    if (md.mpicoordy == 0)
    //    {
    //        set_ghost_cell_w_wrapper.template operator()<Lbc_location::South>();
    //        // sponge_layer_v_wrapper.template operator()<Lbc_location::South>(lbc_s);
    //    }
    //    if (md.mpicoordy == md.npy-1)
    //    {
    //        set_ghost_cell_w_wrapper.template operator()<Lbc_location::North>();
    //        // sponge_layer_v_wrapper.template operator()<Lbc_location::North>(lbc_n);
    //    }

    //    set_corner_ghost_cell_wrapper(fields.mp.at("w")->fld, gd.kend+1);
    //}

    for (auto& fld : inoutflow_s)
    {
        const bool sw_recycle = in_list<std::string>(fld, recycle_list);

        if (md.mpicoordx == 0)
        {
            //set_ghost_cell_s_wrapper.template operator()<Lbc_location::West>(lbc_w, fld);
            if (sw_recycle)
                sponge_layer_wrapper.template operator()<Lbc_location::West, true>(lbc_w, fld);
            else
                sponge_layer_wrapper.template operator()<Lbc_location::West, false>(lbc_w, fld);
        }
        if (md.mpicoordx == md.npx-1)
        {
            //set_ghost_cell_s_wrapper.template operator()<Lbc_location::East>(lbc_e, fld);
            if (sw_recycle)
                sponge_layer_wrapper.template operator()<Lbc_location::East, true>(lbc_e, fld);
            else
                sponge_layer_wrapper.template operator()<Lbc_location::East, false>(lbc_e, fld);
        }
        if (md.mpicoordy == 0)
        {
            //set_ghost_cell_s_wrapper.template operator()<Lbc_location::South>(lbc_s, fld);
            if (sw_recycle)
                sponge_layer_wrapper.template operator()<Lbc_location::South, true>(lbc_s, fld);
            else
                sponge_layer_wrapper.template operator()<Lbc_location::South, false>(lbc_s, fld);
        }
        if (md.mpicoordy == md.npy-1)
        {
            //set_ghost_cell_s_wrapper.template operator()<Lbc_location::North>(lbc_n, fld);
            if (sw_recycle)
                sponge_layer_wrapper.template operator()<Lbc_location::North, true>(lbc_n, fld);
            else
                sponge_layer_wrapper.template operator()<Lbc_location::North, false>(lbc_n, fld);
        }

        //set_corner_ghost_cell_wrapper(fields.ap.at(fld)->fld, gd.kend);
    }

    //if (sw_perturb)
    //{
    //    auto perturb_boundary_wrapper = [&]<Lbc_location location>(
    //            std::map<std::string, std::vector<TF>>& lbc_map,
    //            const std::string& fld)
    //    {
    //        perturb_boundary_kernel<TF, location>(
    //                fields.at.at(fld)->fld.data(),
    //                fields.ap.at(fld)->fld.data(),
    //                lbc_map.at(fld).data(),
    //                perturb_ampl.at(fld),
    //                timeloop.get_sub_time_step(),
    //                perturb_width, perturb_block,
    //                md.mpicoordx, md.npx,
    //                gd.istart, gd.iend,
    //                gd.jstart, gd.jend,
    //                gd.kstart, perturb_kend,
    //                gd.icells, gd.jcells,
    //                gd.ijcells);
    //    };

    //    // Add random perturbation in a certain block size to the fields near the lateral boundaries.
    //    for (auto& fld : perturb_list)
    //    {
    //        if (md.mpicoordx == 0)
    //            perturb_boundary_wrapper.template operator()<Lbc_location::West>(lbc_w, fld);
    //        //if (md.mpicoordx == md.npx-1)
    //        //    perturb_boundary_wrapper.template operator()<Lbc_location::East>(lbc_e, fld);
    //        //if (md.mpicoordy == 0)
    //        //    perturb_boundary_wrapper.template operator()<Lbc_location::South>(lbc_s, fld);
    //        //if (md.mpicoordy == md.npy-1)
    //        //    perturb_boundary_wrapper.template operator()<Lbc_location::North>(lbc_n, fld);
    //    }
    //}
}


template <typename TF>
void Boundary_lateral<TF>::update_time_dependent(
        Timeloop<TF>& timeloop,
        const bool pres_fix)
{
    if (!sw_inoutflow || !sw_timedep)
        return;

    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    // Find index in time array.
    double time = timeloop.get_time();

    // CvH: this is an UGLY hack, because it only works for RK4.
    // We need to know from the time whether we are in the last iter.
    // Also, it will fail miserably if we perturb the velocities at the walls
    if (pres_fix && timeloop.get_substep() == 4)
        time += timeloop.get_dt();

    int t0;
    for (int i=0; i<time_in.size()-1; ++i)
        if (time_in[i] <= time and time_in[i+1] > time)
        {
            t0=i;
            break;
        }

    // Interpolation factor.
    const TF f0 = TF(1) - ((time - time_in[t0]) / (time_in[t0+1] - time_in[t0]));
    const TF f1 = TF(1) - f0;

    // Interpolate mean domain top velocity
    if (sw_wtop_2d)
    {
        unsigned long itime = timeloop.get_itime();

        if (itime > itime_w_top_next)
        {
            // Read new w_top field
            const double ifactor = timeloop.get_ifactor();
            unsigned long iiotimeprec = timeloop.get_iiotimeprec();

            itime_w_top_prev = itime_w_top_next;
            itime_w_top_next = itime_w_top_prev + wtop_2d_loadtime*ifactor;
            const int iotime1 = int(itime_w_top_next / iiotimeprec);

            // Copy of data from next to prev. time
            w_top_prev = w_top_next;

            // Read new w_top slice.
            read_xy_slice(w_top_next, "w_top", iotime1);
        }

        // Interpolate `w_top` field in time.
        for (int n=0; n<gd.ijcells; ++n)
            w_top[n] = f0 * w_top_prev[n] + f1 * w_top_next[n];
    }
    else
    {
        const TF w_top_int = f0 * w_top_in[t0] + f1 * w_top_in[t0+1];
        std::fill(w_top.begin(), w_top.end(), w_top_int);
    }

    // Interpolate boundaries in time.
    if (md.mpicoordx == 0)
    {
        for (auto& it : lbc_w)
        {
            const int ngc = (it.first == "u") ? gd.igc+1 : gd.igc;

            interpolate_lbc_kernel(
                    lbc_w.at(it.first).data(),
                    lbc_w_in.at(it.first).data(),
                    gd.kcells*gd.jcells*ngc,
                    t0, f0);
        }
    }

    if (md.mpicoordx == md.npx-1)
    {
        for (auto& it : lbc_e)
            interpolate_lbc_kernel(
                    lbc_e.at(it.first).data(),
                    lbc_e_in.at(it.first).data(),
                    gd.kcells*gd.jcells*gd.igc,
                    t0, f0);
    }

    if (md.mpicoordy == 0)
    {
        for (auto& it : lbc_s)
        {
            const int ngc = (it.first == "v") ? gd.jgc+1 : gd.jgc;
            interpolate_lbc_kernel(
                    lbc_s.at(it.first).data(),
                    lbc_s_in.at(it.first).data(),
                    gd.kcells*ngc*gd.icells,
                    t0, f0);
        }
    }

    if (md.mpicoordy == md.npy-1)
    {
        for (auto& it : lbc_n)
            interpolate_lbc_kernel(
                    lbc_n.at(it.first).data(),
                    lbc_n_in.at(it.first).data(),
                    gd.kcells*gd.jgc*gd.icells,
                    t0, f0);
    }
}


template <typename TF>
void Boundary_lateral<TF>::read_xy_slice(
        std::vector<TF>& field,
        const std::string& name,
        const int time)
{
    char filename[256];
    std::sprintf(filename, "%s.%07d", name.c_str(), time);
    master.print_message("Loading \"%s\" ... ", filename);

    auto tmp  = fields.get_tmp();

    if (field3d_io.load_xy_slice(field.data(), tmp->fld.data(), filename))
        master.print_message("FAILED\n");
    else
        master.print_message("OK\n");

    fields.release_tmp(tmp);
}

#ifdef FLOAT_SINGLE
template class Boundary_lateral<float>;
#else
template class Boundary_lateral<double>;
#endif
