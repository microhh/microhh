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

#include "boundary_lateral.h"
#include "netcdf_interface.h"
#include "grid.h"
#include "fields.h"
#include "input.h"
#include "master.h"
#include "timeloop.h"

namespace
{
    template<typename TF, Lbc_location location>
    void set_ghost_cell_kernel_u(
            TF* const restrict u,
            const TF* const restrict lbc_u,
            const int istart, const int iend, const int igc,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int jcells,
            const int ijcells)
    {
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
            {
                const int jk = j+k*jcells;

                // Set boundary values directly.
                if (location == Lbc_location::West)
                {
                    const int ijk_b = istart + j*icells + k*ijcells;
                    const int ijk_d = istart+1 + j*icells +k*ijcells;

                    u[ijk_b] = lbc_u[jk];

                    for (int i=0; i<igc; ++i)
                    {
                        const int ijk_gc = ijk_b - (i+1);
                        u[ijk_gc] = lbc_u[jk] - (i+1)*(u[ijk_d]-lbc_u[jk]);
                    }
                }
                else if (location == Lbc_location::East)
                {
                    const int ijk_b = iend + j*icells + k*ijcells;
                    const int ijk_d = iend-1 + j*icells +k*ijcells;

                    u[ijk_b] = lbc_u[jk];

                    for (int i=0; i<igc-1; ++i)
                    {
                        const int ijk_gc = ijk_b + (i+1);
                        u[ijk_gc] = lbc_u[jk] + (i+1)*(lbc_u[jk]-u[ijk_d]);
                    }
                }
            }
    }

    template<typename TF, Lbc_location location>
    void set_ghost_cell_kernel_v(
            TF* const restrict v,
            const TF* const restrict lbc_v,
            const int istart, const int iend,
            const int jstart, const int jend, const int jgc,
            const int kstart, const int kend,
            const int icells, const int jcells,
            const int ijcells)
    {
        for (int k=kstart; k<kend; ++k)
            for (int i=istart; i<iend; ++i)
            {
                const int ik = i+k*icells;

                // Set boundary values directly.
                if (location == Lbc_location::South)
                {
                    const int ijk_b = i + jstart*icells + k*ijcells;
                    const int ijk_d = i + (jstart+1)*icells +k*ijcells;

                    v[ijk_b] = lbc_v[ik];

                    for (int j=0; j<jgc; ++j)
                    {
                        const int ijk_gc = ijk_b - (j+1)*icells;
                        v[ijk_gc] = lbc_v[ik] - (j+1)*(v[ijk_d]-lbc_v[ik]);
                    }
                }
                else if (location == Lbc_location::North)
                {
                    const int ijk_b = i + jend*icells + k*ijcells;
                    const int ijk_d = i + (jend-1)*icells +k*ijcells;

                    v[ijk_b] = lbc_v[ik];

                    for (int j=0; j<jgc-1; ++j)
                    {
                        const int ijk_gc = ijk_b + (j+1)*icells;
                        v[ijk_gc] = lbc_v[ik] + (j+1)*(lbc_v[ik]-v[ijk_d]);
                    }
                }
            }
    }

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


    template<typename TF, Lbc_location location>
    void set_ghost_cell_kernel_s(
            TF* const restrict a,
            const TF* const restrict lbc,
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
                {
                    const int jk = j+k*jcells;

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

                        a[ijk_gc] = a[ijk_d] - (i+1)*TF(2)*(a[ijk_d] - lbc[jk]);
                    }
                }

        }
        else if (location == Lbc_location::North || location == Lbc_location::South)
        {
            for (int k=kstart; k<kend; ++k)
                for (int i=istart; i<iend; ++i)
                {
                    const int ik = i+k*icells;

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

                        const TF lbc_val = lbc[ik];
                        a[ijk_gc] = a[ijk_d] - (j+1)*TF(2)*(a[ijk_d] - lbc[ik]);
                    }
                }
        }
    }

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
                const int jk = j + k*jcells;

                for (int n=2; n<=N_sponge; ++n)
                {
                    const int i = (location==Lbc_location::West) ? istart+(n-1) : iend-(n-1);
                    const int ijk = i + j*icells + k*ijcells;

                    // Calculate diffusion term over 3x3x3 stencil.
                    // Offset block near lateral boundaries to avoid using ghost cells.
                    // No offset needed for `i`, as stencil center starts at `istart+1` or `iend-1`,
                    // and `istart` and `iend` contain the correct boundary values.
                    const int jo =
                            (mpiidy == 0 && j == jstart) ? 1 :
                            (mpiidy == npy-1 && j == jend-1) ? -1 : 0;

                    const int ko =
                            (k == kstart) ? 1 :
                            (k == kend-1) ? -1 : 0;

                    const int ijkc = i + (j+jo)*icells + (k+ko)*(ijcells);

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
                const int ik = i + k*icells;

                for (int n=2; n<=N_sponge; ++n)
                {
                    const int j = (location==Lbc_location::South) ? jstart+(n-1) : jend-(n-1);
                    const int ijk = i + j*icells + k*ijcells;

                    // Calculate diffusion term over 3x3x3 stencil.
                    // Offset block near lateral boundaries to avoid using ghost cells.
                    // No offset needed for `j`, as stencil center starts at `jstart+1` or `jend-1`,
                    // and `jstart` and `jend` contain the correct boundary values.
                    const int io =
                            (mpiidx == 0 && i == istart) ? 1 :
                            (mpiidx == npx-1 && i == iend-1) ? -1 : 0;

                    const int ko =
                            (k == kstart) ? 1 :
                            (k == kend-1) ? -1 : 0;

                    const int ijkc = i+io + j*icells + (k+ko)*(ijcells);

                    const TF v_diff = diffusion_3x3x3(
                            v, ijkc, icells, ijcells);

                    // Nudge coefficient.
                    const TF w1 = w_dt * (TF(1)+N_sponge-n) / N_sponge;

                    vt[ijk] += w1*(lbc_v[ik]-v[ijk]);
                    vt[ijk] -= w_diff*v_diff;
                }
            }
    }

    template<typename TF, Lbc_location location>
    void lateral_sponge_kernel_s(
            TF* const restrict at,
            const TF* const restrict a,
            const TF* const restrict lbc,
            const TF tau_nudge,
            const TF w_diff,
            const int N_sponge,
            const int npx, const int npy,
            const int mpiidx, const int mpiidy,
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

        if (location == Lbc_location::West || location == Lbc_location::East)
        {
            for (int k=kstart; k<kend; ++k)
                for (int n=1; n<=N_sponge; ++n)
                {
                    const int jstart_loc = mpiidy == 0     ? jstart + (n-1) : jstart;
                    const int jend_loc   = mpiidy == npy-1 ? jend   - (n-1) : jend;

                    for (int j=jstart_loc; j<jend_loc; ++j)
                    {
                        const int jk = j + k*jcells;

                        const int i = (location==Lbc_location::West) ? istart+(n-1) : iend-n;
                        const int ijk = i + j*icells + k*ijcells;

                        // Calculate diffusion term over 3x3x3 stencil.
                        // Offset block near lateral boundaries to avoid using ghost cells.
                        const int io =
                                (location == Lbc_location::West && n==1) ? 1 :
                                (location == Lbc_location::East && n==1) ? -1 : 0;

                        const int jo =
                                (mpiidy == 0 && j == jstart) ? 1 :
                                (mpiidy == npy-1 && j == jend-1) ? -1 : 0;

                        const int ko =
                                (k == kstart) ? 1 :
                                (k == kend-1) ? -1 : 0;

                        const int ijkc = (i+io) + (j+jo)*icells + (k+ko)*(ijcells);

                        const TF a_diff = diffusion_3x3x3(
                                a, ijkc, icells, ijcells);

                        // Nudge coefficient.
                        const TF w1 = w_dt * (TF(1)+N_sponge-(n+TF(0.5))) / N_sponge;

                        at[ijk] += w1*(lbc[jk]-a[ijk]);
                        at[ijk] -= w_diff*a_diff;
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
                        const int ik = i + k*icells;

                        const int j = (location==Lbc_location::South) ? jstart+(n-1) : jend-n;
                        const int ijk = i + j*icells + k*ijcells;

                        // Calculate diffusion term over 3x3x3 stencil.
                        // Offset block near lateral boundaries to avoid using ghost cells.
                        const int io =
                                (mpiidx == 0 && i == istart) ? 1 :
                                (mpiidx == npx-1 && i == iend-1) ? -1 : 0;

                        const int jo =
                                (location == Lbc_location::South && n==1) ? 1 :
                                (location == Lbc_location::North && n==1) ? -1 : 0;

                        const int ko =
                                (k == kstart) ? 1 :
                                (k == kend-1) ? -1 : 0;

                        const int ijkc = (i+io) + (j+jo)*icells + (k+ko)*(ijcells);

                        const TF a_diff = diffusion_3x3x3(
                                a, ijkc, icells, ijcells);

                        // Nudge coefficient.
                        const TF w1 = w_dt * (TF(1)+N_sponge-(n+TF(0.5))) / N_sponge;

                        at[ijk] += w1*(lbc[ik]-a[ijk]);
                        at[ijk] -= w_diff*a_diff;
                    }
                }
        }
    }


    template<typename TF>
    void blend_corners_kernel(
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


    template<typename TF>
    void calc_div_h(
            TF* const restrict div,
            const TF* const restrict lbc_u,
            const TF* const restrict lbc_d,
            const TF* const restrict rhoref,
            const TF* const restrict dz,
            const TF dx_or_dy,
            const int ntime,
            const int ntot, const int ktot,
            const int kgc)
    {
        for (int t=0; t<ntime; ++t)
        {
            div[t] = TF(0);
            for (int k=0; k<ktot; ++k)
                for (int n=0; n<ntot; ++n)
                {
                    const int nk = n + k*ntot + t*ntot*ktot;

                    div[t] += rhoref[k+kgc] * dx_or_dy * dz[k+kgc] * (lbc_u[nk] - lbc_d[nk]);
                }
        }
    }

    template<typename TF>
    void check_div(
        const TF* const restrict u,
        const TF* const restrict v,
        const TF* const restrict w,
        const TF* const restrict dzi,
        const TF* const restrict rhoref,
        const TF* const restrict rhorefh,
        const TF dx, const TF dy,
        const int istart, const int iend,
        const int jstart, const int jend,
        const int kstart, const int kend,
        const int icells, const int ijcells,
        Master& master)
    {
        const int ii = 1;
        const int jj = icells;
        const int kk = ijcells;

        const TF dxi = TF(1.)/dx;
        const TF dyi = TF(1.)/dy;

        TF divmax = 0.;
        int imax = 0;
        int jmax = 0;
        int kmax = 0;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {

                    const int ijk = i + j*jj + k*kk;
                    const TF div = rhoref[k]*((u[ijk+ii]-u[ijk])*dxi + (v[ijk+jj]-v[ijk])*dyi)
                          + (rhorefh[k+1]*w[ijk+kk]-rhorefh[k]*w[ijk])*dzi[k];

                    if (div > divmax)
                    {
                        divmax = div;
                        imax = i;
                        jmax = j;
                        kmax = k;
                    }
                }

        master.max(&divmax, 1);

        std::string message = "Max div. = " + std::to_string(divmax) + " @ i,j,k = "
                + std::to_string(imax) + ", "
                + std::to_string(jmax) + ", "
                + std::to_string(kmax);

        master.print_message(message);
    }


    template<typename TF, Lbc_location location>
    void perturb_boundary_kernel(
            TF* const restrict tend,
            const TF amplitude, const TF dt,
            const int width, const int block_size,
            const int mpiidx, const int npx,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
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
                        const TF tend_value = (TF(std::rand()) * rand_max_i - TF(0.5)) * amplitude * dt_i;

                        for (int dj=0; dj<block_size; ++dj)
                            for (int di=0; di<block_size; ++di)
                            {
                                const int i = i0 + bi*block_size + di;
                                const int j = jstart + bj*block_size + dj;
                                const int ijk = i + j*icells + k*ijcells;

                                tend[ijk] += tend_value;
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
                        const TF tend_value = (TF(std::rand()) * rand_max_i - TF(0.5)) * amplitude * dt_i;

                        for (int dj=0; dj<block_size; ++dj)
                            for (int di=0; di<block_size; ++di)
                            {
                                const int i = istart_loc + bi*block_size + di;
                                const int j = j0 + bj*block_size + dj;
                                const int ijk = i + j*icells + k*ijcells;

                                tend[ijk] += tend_value;
                            }
                    }
        }
    }
}


template<typename TF>
Boundary_lateral<TF>::Boundary_lateral(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
        master(masterin), grid(gridin), fields(fieldsin)
{
    sw_inoutflow = inputin.get_item<bool>("boundary", "sw_inoutflow", "", false);

    if (sw_inoutflow)
    {
        sw_timedep = inputin.get_item<bool>("boundary", "sw_timedep", "", false);
        sw_inoutflow_u = inputin.get_item<bool>("boundary", "sw_inoutflow_u", "", true);
        sw_inoutflow_v = inputin.get_item<bool>("boundary", "sw_inoutflow_v", "", true);
        sw_inoutflow_w = inputin.get_item<bool>("boundary", "sw_inoutflow_w", "", true);
        inoutflow_s = inputin.get_list<std::string>("boundary", "inoutflow_slist", "", std::vector<std::string>());

        // Lateral sponge / diffusion layer.
        sw_sponge = inputin.get_item<bool>("boundary", "sw_sponge", "", true);
        if (sw_sponge)
        {
            n_sponge = inputin.get_item<int>("boundary", "n_sponge", "", 5);
            tau_nudge = inputin.get_item<TF>("boundary", "tau_nudge", "", 60);
            w_diff = inputin.get_item<TF>("boundary", "w_diff", "", 0.0033);
        }

        // Inflow perturbations
        sw_perturb = inputin.get_item<bool>("boundary", "sw_perturb", "", false);
        if (sw_perturb)
        {
            perturb_list = inputin.get_list<std::string>("boundary", "perturb_list", "", std::vector<std::string>());
            perturb_width = inputin.get_item<int>("boundary", "perturb_width", "", 4);
            perturb_block = inputin.get_item<int>("boundary", "perturb_block", "", 2);
            perturb_seed = inputin.get_item<int>("boundary", "perturb_seed", "", 0);

            for (auto& fld : perturb_list)
            {
                const TF ampl = inputin.get_item<TF>("boundary", "perturb_ampl", fld);
                perturb_ampl.emplace(fld, ampl);
            }
        }
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

    auto add_lbc = [&](const std::string& name)
    {
        if (md.mpicoordx == 0)
            lbc_w.emplace(name, std::vector<TF>(gd.kcells*gd.jcells));
        if (md.mpicoordx == md.npx-1)
            lbc_e.emplace(name, std::vector<TF>(gd.kcells*gd.jcells));
        if (md.mpicoordy == 0)
            lbc_s.emplace(name, std::vector<TF>(gd.kcells*gd.icells));
        if (md.mpicoordy == md.npy-1)
            lbc_n.emplace(name, std::vector<TF>(gd.kcells*gd.icells));
    };

    if (sw_inoutflow_u)
        add_lbc("u");
    if (sw_inoutflow_v)
        add_lbc("v");

    for (auto& fld : inoutflow_s)
        add_lbc(fld);

    // Make sure every MPI task has different seed.
    if (sw_perturb)
    {
        perturb_seed += master.get_mpiid();
        std::srand(perturb_seed);
    }
}

template <typename TF>
void Boundary_lateral<TF>::create(Input& inputin, const std::string& sim_name)
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
    if (sw_inoutflow_u && sw_inoutflow_v)
    {
        div_u.resize(ntime);
        div_v.resize(ntime);
        w_top_in.resize(ntime);
    }

    auto copy_boundary = [&](
            std::map<std::string, std::vector<TF>>& map_out,
            std::vector<TF>& fld_in,
            const int istart, const int iend,
            const int igc, const int imax,
            const int istride_in, const int istride_out,
            const int mpicoord, const std::string& name)
    {
        std::vector<TF> fld_out = std::vector<TF>(ntime * gd.kcells * istride_out);

        for (int t=0; t<ntime; t++)
            for (int k=gd.kstart; k<gd.kend; k++)
                for (int i=istart; i<iend; i++)
                {
                    const int kk = k-gd.kgc;
                    const int ii = i-igc;

                    const int ikt_out = i + k*istride_out + t*istride_out*gd.kcells;
                    const int ikt_in = (ii+mpicoord*imax) + kk*istride_in + t*istride_in*gd.ktot;

                    fld_out[ikt_out] = fld_in[ikt_in];
                }

        map_out.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(name),
                std::forward_as_tuple(std::move(fld_out)));
    };

    auto copy_boundaries = [&](const std::string& name)
    {
        std::vector<TF> lbc_w_full = input_nc.get_variable<TF>(name + "_west", {ntime, gd.ktot, gd.jtot});
        std::vector<TF> lbc_e_full = input_nc.get_variable<TF>(name + "_east", {ntime, gd.ktot, gd.jtot});
        std::vector<TF> lbc_s_full = input_nc.get_variable<TF>(name + "_south", {ntime, gd.ktot, gd.itot});
        std::vector<TF> lbc_n_full = input_nc.get_variable<TF>(name + "_north", {ntime, gd.ktot, gd.itot});

        if (md.mpicoordx == 0)
            copy_boundary(
                    lbc_w_in, lbc_w_full,
                    gd.jstart, gd.jend, gd.jgc, gd.jmax,
                    gd.jtot, gd.jcells, md.mpicoordy, name);

        if (md.mpicoordx == md.npx-1)
            copy_boundary(
                    lbc_e_in, lbc_e_full,
                    gd.jstart, gd.jend, gd.jgc, gd.jmax,
                    gd.jtot, gd.jcells, md.mpicoordy, name);

        if (md.mpicoordy == 0)
            copy_boundary(
                    lbc_s_in, lbc_s_full,
                    gd.istart, gd.iend, gd.igc, gd.imax,
                    gd.itot, gd.icells, md.mpicoordx, name);

        if (md.mpicoordy == md.npy-1)
            copy_boundary(
                    lbc_n_in, lbc_n_full,
                    gd.istart, gd.iend, gd.igc, gd.imax,
                    gd.itot, gd.icells, md.mpicoordx, name);

        // Calculate domain total mass imbalance in kg s-1.
        if (name == "u")
            calc_div_h(
                    div_u.data(),
                    lbc_e_full.data(),
                    lbc_w_full.data(),
                    fields.rhoref.data(),
                    gd.dz.data(),
                    gd.dy,
                    ntime,
                    gd.jtot, gd.ktot, gd.kgc);
        else if (name == "v")
            calc_div_h(
                    div_v.data(),
                    lbc_n_full.data(),
                    lbc_s_full.data(),
                    fields.rhoref.data(),
                    gd.dz.data(),
                    gd.dx,
                    ntime,
                    gd.itot, gd.ktot, gd.kgc);

        //if (!sw_timedep)
        //{
            if (md.mpicoordx == 0)
            {
                for (int n=0; n<gd.jcells*gd.kcells; ++n)
                    lbc_w.at(name)[n] = lbc_w_in.at(name)[n];
            }

            if (md.mpicoordx == md.npx-1)
            {
                for (int n=0; n<gd.jcells*gd.kcells; ++n)
                    lbc_e.at(name)[n] = lbc_e_in.at(name)[n];
            }

            if (md.mpicoordy == 0)
            {
                for (int n=0; n<gd.icells*gd.kcells; ++n)
                    lbc_s.at(name)[n] = lbc_s_in.at(name)[n];
            }

            if (md.mpicoordy == md.npy-1)
            {
                for (int n=0; n<gd.icells*gd.kcells; ++n)
                    lbc_n.at(name)[n] = lbc_n_in.at(name)[n];
            }
        //}
    };

    if (sw_inoutflow_u)
        copy_boundaries("u");

    if (sw_inoutflow_v)
        copy_boundaries("v");

    for (auto& fld : inoutflow_s)
        copy_boundaries(fld);

    // Calculate domain mean vertical velocity.
    if (sw_inoutflow_u && sw_inoutflow_v)
    {
        for (int t=0; t<ntime; ++t)
        {
            // w_top is the total mass in/outflow at the top divided by the total area and the local density.
            w_top_in[t] = -(div_u[t] + div_v[t]) / (fields.rhorefh[gd.kend]*gd.xsize*gd.ysize);

            std::string message = "<w_top> =" + std::to_string(w_top_in[t]) + " m/s @ t=" + std::to_string(time_in[t]);
            master.print_message(message);
        }

        if (!sw_timedep)
            w_top = w_top_in[0];
    }
}

template <typename TF>
void Boundary_lateral<TF>::set_ghost_cells(Timeloop<TF>& timeloop)
{
    if (!sw_inoutflow)
        return;

    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    auto dump_fld3d = [&](
            std::vector<TF>& fld,
            const std::string& name)
    {
        if (md.npx > 1 || md.npy > 1)
            throw std::runtime_error("Raw dump function does not support npx,npy>1");

        FILE *pFile;
        pFile = fopen(name.c_str(), "wbx");

        if (pFile == NULL)
            throw std::runtime_error("Opening raw dump field failed.");

        const int jj = gd.icells;
        const int kk = gd.icells*gd.jcells;

        for (int k=0; k<gd.kcells; ++k)
            for (int j=0; j<gd.jcells; ++j)
            {
                const int ijk = j*jj + k*kk;
                fwrite(&fld.data()[ijk], sizeof(TF), gd.icells, pFile);
            }
        fclose(pFile);
    };

    auto set_ghost_cell_s_wrapper = [&]<Lbc_location location>(
            std::map<std::string, std::vector<TF>>& lbc_map,
            const std::string& name)
    {
        set_ghost_cell_kernel_s<TF, location>(
                fields.ap.at(name)->fld.data(),
                lbc_map.at(name).data(),
                gd.istart, gd.iend, gd.igc,
                gd.jstart, gd.jend, gd.jgc,
                gd.kstart, gd.kend,
                gd.icells, gd.jcells, gd.kcells,
                gd.ijcells);
    };

    auto sponge_layer_wrapper = [&]<Lbc_location location>(
            std::map<std::string, std::vector<TF>>& lbc_map,
            const std::string& name)
    {
        if (!sw_sponge)
            return;

        lateral_sponge_kernel_s<TF, location>(
                fields.at.at(name)->fld.data(),
                fields.ap.at(name)->fld.data(),
                lbc_map.at(name).data(),
                tau_nudge,
                w_diff,
                n_sponge,
                md.npx, md.npy,
                md.mpicoordx, md.mpicoordy,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.jcells,
                gd.ijcells);
    };

    auto set_ghost_cell_u_wrapper = [&]<Lbc_location location>(
            std::map<std::string, std::vector<TF>>& lbc_map)
    {
        set_ghost_cell_kernel_u<TF, location>(
                fields.mp.at("u")->fld.data(),
                lbc_map.at("u").data(),
                gd.istart, gd.iend, gd.igc,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.jcells,
                gd.ijcells);
    };

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
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.jcells,
                gd.ijcells);
    };

    auto set_ghost_cell_v_wrapper = [&]<Lbc_location location>(
            std::map<std::string, std::vector<TF>>& lbc_map)
    {
        set_ghost_cell_kernel_v<TF, location>(
                fields.mp.at("v")->fld.data(),
                lbc_map.at("v").data(),
                gd.istart, gd.iend,
                gd.jstart, gd.jend, gd.jgc,
                gd.kstart, gd.kend,
                gd.icells, gd.jcells,
                gd.ijcells);
    };

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
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, gd.kend,
                gd.icells, gd.jcells,
                gd.ijcells);
    };

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

    auto blend_corners_wrapper = [&](
            std::vector<TF>& fld,
            const int kend)
    {
        blend_corners_kernel(
                fld.data(),
                md.mpicoordx, md.mpicoordy,
                md.npx, md.npy,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, kend,
                gd.icells, gd.jcells,
                gd.kcells, gd.ijcells);
    };


    if (sw_inoutflow_u)
    {
        if (md.mpicoordy == 0)
        {
            set_ghost_cell_s_wrapper.template operator()<Lbc_location::South>(lbc_s, "u");
            sponge_layer_wrapper.template operator()<Lbc_location::South>(lbc_s, "u");
        }
        if (md.mpicoordy == md.npy-1)
        {
            set_ghost_cell_s_wrapper.template operator()<Lbc_location::North>(lbc_n, "u");
            sponge_layer_wrapper.template operator()<Lbc_location::North>(lbc_n, "u");
        }
        if (md.mpicoordx == 0)
        {
            set_ghost_cell_u_wrapper.template operator()<Lbc_location::West>(lbc_w);
            sponge_layer_u_wrapper.template operator()<Lbc_location::West>(lbc_w);
        }
        if (md.mpicoordx == md.npx-1)
        {
            set_ghost_cell_u_wrapper.template operator()<Lbc_location::East>(lbc_e);
            sponge_layer_u_wrapper.template operator()<Lbc_location::East>(lbc_e);
        }

        blend_corners_wrapper(fields.mp.at("u")->fld, gd.kend);
    }

    if (sw_inoutflow_v)
    {
        if (md.mpicoordx == 0)
        {
            set_ghost_cell_s_wrapper.template operator()<Lbc_location::West>(lbc_w, "v");
            sponge_layer_wrapper.template operator()<Lbc_location::West>(lbc_w, "v");
        }
        if (md.mpicoordx == md.npx-1)
        {
            set_ghost_cell_s_wrapper.template operator()<Lbc_location::East>(lbc_e, "v");
            sponge_layer_wrapper.template operator()<Lbc_location::East>(lbc_e, "v");
        }
        if (md.mpicoordy == 0)
        {
            set_ghost_cell_v_wrapper.template operator()<Lbc_location::South>(lbc_s);
            sponge_layer_v_wrapper.template operator()<Lbc_location::South>(lbc_s);
        }
        if (md.mpicoordy == md.npy-1)
        {
            set_ghost_cell_v_wrapper.template operator()<Lbc_location::North>(lbc_n);
            sponge_layer_v_wrapper.template operator()<Lbc_location::North>(lbc_n);
        }

        blend_corners_wrapper(fields.mp.at("v")->fld, gd.kend);
    }

    if (sw_inoutflow_u || sw_inoutflow_v)
    {
        //const int jstride = gd.icells;
        //const int kstride = gd.icells*gd.jcells;

        //// Correct the vertical velocity top BC for the lateral BCs to ensure divergence free field.
        //// CvH THIS IS A FIRST ATTEMPT THAT ASSUMES UNIFORM W.
        //TF w_top = TF(0);

        //for (int k=gd.kstart; k<gd.kend; ++k)
        //{
        //    const TF hor_div = fields.rhoref[k]*(TF(0.) - TF(2.)) / gd.xsize  // * gd.z[k]/gd.zsize
        //                     + fields.rhoref[k]*(TF(2.) - TF(0.)) / gd.ysize; // * gd.z[k]/gd.zsize;

        //    w_top = (fields.rhorefh[k]*w_top - gd.dz[k]*hor_div) / fields.rhorefh[k+1];
        //}

        //for (int j=gd.jstart; j<gd.jend; ++j)
        //    for (int i=gd.istart; i<gd.iend; ++i)
        //    {
        //        const int ijk = i + j*jstride + gd.kend*kstride;
        //        fields.mp.at("w")->fld[ijk] = w_top;
        //    }


        //const int ii = 1;
        //const int jj = gd.icells;
        //const int kk = gd.ijcells;

        //TF* u = fields.mp.at("u")->fld.data();
        //TF* v = fields.mp.at("v")->fld.data();
        //TF* w = fields.mp.at("w")->fld.data();

        //TF* rho  = fields.rhoref.data();
        //TF* rhoh = fields.rhorefh.data();

        //const TF dxi = TF(1)/gd.dx;
        //const TF dyi = TF(1)/gd.dy;

        //const int k = gd.kend-1;
        //for (int j=gd.jstart; j<gd.jend; ++j)
        //    for (int i=gd.istart; i<gd.iend; ++i)
        //    {
        //        const int ijk = i + j*jj + k*kk;
        //        w[ijk+kk] = -(rho[k] * ((u[ijk+ii] - u[ijk]) * dxi + (v[ijk+jj] - v[ijk]) * dyi) * gd.dz[k] - rhoh[k] * w[ijk]) / rhoh[k+1];
        //    }

        const int k = gd.kend;
        for (int j=gd.jstart; j<gd.jend; ++j)
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*gd.icells + k*gd.ijcells;
                fields.mp.at("w")->fld[ijk] = w_top;
            }

        //check_div(
        //    fields.mp.at("u")->fld.data(),
        //    fields.mp.at("v")->fld.data(),
        //    fields.mp.at("w")->fld.data(),
        //    gd.dzi.data(),
        //    fields.rhoref.data(),
        //    fields.rhorefh.data(),
        //    gd.dx, gd.dy,
        //    gd.istart, gd.iend,
        //    gd.jstart, gd.jend,
        //    gd.kstart, gd.kend,
        //    gd.icells, gd.ijcells,
        //    master);
    }
    // END

    // Here, we enfore a Neumann BC of 0 over the boundaries. This works if the large scale w is approximately
    // constant in the horizontal plane. Note that w must be derived from u and v if the large-scale field is to
    // be divergence free. If there is a horizontal gradient in w_top, then it is probably better to extrapolate that
    // gradient into the ghost cells.
    if (sw_inoutflow_w)
    {
        if (md.mpicoordx == 0)
        {
            set_ghost_cell_w_wrapper.template operator()<Lbc_location::West>();
            // sponge_layer_wrapper.template operator()<Lbc_location::West>(lbc_w, "v");
        }
        if (md.mpicoordx == md.npx-1)
        {
            set_ghost_cell_w_wrapper.template operator()<Lbc_location::East>();
            // sponge_layer_wrapper.template operator()<Lbc_location::East>(lbc_e, "v");
        }
        if (md.mpicoordy == 0)
        {
            set_ghost_cell_w_wrapper.template operator()<Lbc_location::South>();
            // sponge_layer_v_wrapper.template operator()<Lbc_location::South>(lbc_s);
        }
        if (md.mpicoordy == md.npy-1)
        {
            set_ghost_cell_w_wrapper.template operator()<Lbc_location::North>();
            // sponge_layer_v_wrapper.template operator()<Lbc_location::North>(lbc_n);
        }

        blend_corners_wrapper(fields.mp.at("w")->fld, gd.kend+1);
    }

    for (auto& fld : inoutflow_s)
    {
        if (md.mpicoordx == 0)
        {
            set_ghost_cell_s_wrapper.template operator()<Lbc_location::West>(lbc_w, fld);
            sponge_layer_wrapper.template operator()<Lbc_location::West>(lbc_w, fld);
        }
        if (md.mpicoordx == md.npx-1)
        {
            set_ghost_cell_s_wrapper.template operator()<Lbc_location::East>(lbc_e, fld);
            sponge_layer_wrapper.template operator()<Lbc_location::East>(lbc_e, fld);
        }
        if (md.mpicoordy == 0)
        {
            set_ghost_cell_s_wrapper.template operator()<Lbc_location::South>(lbc_s, fld);
            sponge_layer_wrapper.template operator()<Lbc_location::South>(lbc_s, fld);
        }
        if (md.mpicoordy == md.npy-1)
        {
            set_ghost_cell_s_wrapper.template operator()<Lbc_location::North>(lbc_n, fld);
            sponge_layer_wrapper.template operator()<Lbc_location::North>(lbc_n, fld);
        }

        blend_corners_wrapper(fields.ap.at(fld)->fld, gd.kend);
    }


    if (sw_perturb)
    {
        auto perturb_boundary_wrapper = [&]<Lbc_location location>(
                const std::string& fld)
        {
            perturb_boundary_kernel<TF, location>(
                    fields.at.at(fld)->fld.data(),
                    perturb_ampl.at(fld),
                    timeloop.get_sub_time_step(),
                    perturb_width, perturb_block,
                    md.mpicoordx, md.npx,
                    gd.istart, gd.iend,
                    gd.jstart, gd.jend,
                    gd.kstart, gd.kend,
                    gd.icells, gd.ijcells);
        };

        // Add random perturbation in a certain block size to the fields near the lateral boundaries.
        for (auto& fld : perturb_list)
        {
            if (md.mpicoordx == 0)
                perturb_boundary_wrapper.template operator()<Lbc_location::West>(fld);
            if (md.mpicoordx == md.npx-1)
                perturb_boundary_wrapper.template operator()<Lbc_location::East>(fld);
            if (md.mpicoordy == 0)
                perturb_boundary_wrapper.template operator()<Lbc_location::South>(fld);
            if (md.mpicoordy == md.npy-1)
                perturb_boundary_wrapper.template operator()<Lbc_location::North>(fld);
        }
    }

    //dump_fld3d(fields.ap.at("th")->fld, "th0");
    //dump_fld3d(fields.mp.at("u")->fld, "u0");
    //dump_fld3d(fields.mp.at("v")->fld, "v0");
    //dump_fld3d(fields.mp.at("w")->fld, "w0");
    //throw 1;
}


template <typename TF>
void Boundary_lateral<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (!sw_inoutflow || !sw_timedep)
        return;

    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    // Find index in time array.
    const double time = timeloop.get_time();

    int t0;
    for (int i=0; i<time_in.size()-1; ++i)
        if (time_in[i] <= time and time_in[i+1] > time)
        {
            t0=i;
            break;
        }

    // Interpolation factor.
    const TF f0 = TF(1) - ((time - time_in[t0]) / (time_in[t0+1] - time_in[t0]));

    // Interpolate mean domain top velocity
    w_top = f0 * w_top_in[t0] + (TF(1) - f0) * w_top_in[t0+1];

    // Interpolate boundaries in time.
    if (md.mpicoordx == 0)
    {
        for (auto& it : lbc_w)
            interpolate_lbc_kernel(
                lbc_w.at(it.first).data(),
                lbc_w_in.at(it.first).data(),
                gd.kcells*gd.jcells,
                t0, f0);
    }

    if (md.mpicoordx == md.npx-1)
    {
        for (auto& it : lbc_e)
            interpolate_lbc_kernel(
                    lbc_e.at(it.first).data(),
                    lbc_e_in.at(it.first).data(),
                    gd.kcells*gd.jcells,
                    t0, f0);
    }

    if (md.mpicoordy == 0)
    {
        for (auto& it : lbc_s)
            interpolate_lbc_kernel(
                    lbc_s.at(it.first).data(),
                    lbc_s_in.at(it.first).data(),
                    gd.kcells*gd.icells,
                    t0, f0);
    }

    if (md.mpicoordy == md.npy-1)
    {
        for (auto& it : lbc_n)
            interpolate_lbc_kernel(
                    lbc_n.at(it.first).data(),
                    lbc_n_in.at(it.first).data(),
                    gd.kcells*gd.icells,
                    t0, f0);
    }
}

template class Boundary_lateral<double>;
template class Boundary_lateral<float>;
