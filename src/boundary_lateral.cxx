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


    template<typename TF>
    TF diffusion_3x3x3(
        const TF* const restrict fld,
        const TF lbc_val,
        const int ijk,
        const int icells,
        const int ijcells)
    {
        auto index = [&](
                const int i3, const int j3, const int k3)
        {
            return ijk + i3-1 + (j3-1)*icells + (k3-1)*ijcells;
        };

        //const TF fld_diff =
        //        - TF(1) * fld[index(0,0,0)] + TF(2) * fld[index(0,1,0)] - TF(1) * fld[index(0,2,0)]
        //        + TF(2) * fld[index(1,0,0)] - TF(4) * fld[index(1,1,0)] + TF(2) * fld[index(1,2,0)]
        //        - TF(1) * fld[index(2,0,0)] + TF(2) * fld[index(2,1,0)] - TF(1) * fld[index(2,2,0)]
        //        + TF(2) * fld[index(0,0,1)] - TF(4) * fld[index(0,1,1)] + TF(2) * fld[index(0,2,1)]
        //        - TF(4) * fld[index(1,0,1)] + TF(8) * fld[index(1,1,1)] - TF(4) * fld[index(1,2,1)]
        //        + TF(2) * fld[index(2,0,1)] - TF(4) * fld[index(2,1,1)] + TF(2) * fld[index(2,2,1)]
        //        - TF(1) * fld[index(0,0,2)] + TF(2) * fld[index(0,1,2)] - TF(1) * fld[index(0,2,2)]
        //        + TF(2) * fld[index(1,0,2)] - TF(4) * fld[index(1,1,2)] + TF(2) * fld[index(1,2,2)]
        //        - TF(1) * fld[index(2,0,2)] + TF(2) * fld[index(2,1,2)] - TF(1) * fld[index(2,2,2)];

        const TF vc = lbc_val - fld[ijk];
        const TF v1 = lbc_val - fld[ijk-1];
        const TF v2 = lbc_val - fld[ijk+1];
        const TF v3 = lbc_val - fld[ijk-icells];
        const TF v4 = lbc_val - fld[ijk+icells];

        const TF fld_diff = v1 + v2 + v3 + v3 - TF(4) * vc;

        return fld_diff;
    }


    template<typename TF, Lbc_location location>
    void lateral_sponge_kernel_u(
            TF* const restrict ut,
            const TF* const restrict u,
            const TF* const restrict lbc_u,
            const TF tau_sponge,
            const TF w_diff,
            const int nsponge,
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

        const int igc_pad = (location==Lbc_location::West) ? igc+1 : igc;
        const int nstart = (location==Lbc_location::West) ? 2 : 1;
        const int jstride_lbc = igc_pad+nsponge;

        const TF w_dt = TF(1) / tau_sponge;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int n=nstart; n<=nsponge; ++n)
                {
                    const int ilbc = (location==Lbc_location::West) ? igc+n-1 : nsponge-n;
                    const int ijk_lbc = ilbc + j*jstride_lbc + k*jstride_lbc*jcells;

                    const int i = (location==Lbc_location::West) ? istart+(n-1) : iend-n;
                    const int ijk = i + j*icells + k*ijcells;

                    const TF u_diff = diffusion_3x3x3(
                            u, lbc_u[ijk_lbc], ijk, icells, ijcells);

                    // Nudge coefficient.
                    const TF f_sponge = (TF(1)+nsponge-n) / nsponge;
                    const TF w1n = w_dt * f_sponge;
                    const TF w2n = w_diff * f_sponge;

                    ut[ijk] += w1n * (lbc_u[ijk_lbc]-u[ijk]);
                    ut[ijk] -= w2n * u_diff;
                }
    }

    template<typename TF, Lbc_location location>
    void lateral_sponge_kernel_v(
            TF* const restrict vt,
            const TF* const restrict v,
            const TF* const restrict lbc_v,
            const TF tau_sponge,
            const TF w_diff,
            const int nsponge,
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

        const int jgc_pad = (location==Lbc_location::South) ? jgc+1 : jgc;

        const TF w_dt = TF(1) / tau_sponge;

        for (int k=kstart; k<kend; ++k)
            for (int i=istart; i<iend; ++i)
                for (int n=2; n<=nsponge; ++n)
                {
                    const int kstride_lbc = jgc + nsponge;
                    const int jlbc = (location==Lbc_location::South) ? jgc+n-1 : nsponge-n;
                    const int ijk_lbc = i + jlbc*icells + k*icells*(jgc_pad+nsponge);

                    const int j = (location==Lbc_location::South) ? jstart+(n-1) : jend-(n-1);
                    const int ijk = i + j*icells + k*ijcells;

                    const TF v_diff = diffusion_3x3x3(
                            v, lbc_v[ijk_lbc], ijk, icells, ijcells);

                    // Nudge coefficient.
                    const TF f_sponge = (TF(1)+nsponge-n) / nsponge;
                    const TF w1n = w_dt * f_sponge;
                    const TF w2n = w_diff * f_sponge;

                    vt[ijk] += w1n * (lbc_v[ijk_lbc]-v[ijk]);
                    vt[ijk] -= w2n * v_diff;
                }
    }

    template<typename TF, Lbc_location location, bool sw_recycle>
    void lateral_sponge_kernel_s(
            TF* const restrict at,
            const TF* const restrict a,
            const TF* const restrict lbc,
            const TF tau_sponge,
            const TF w_diff,
            const int nsponge,
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

        const TF w_dt = TF(1) / tau_sponge;
        const TF r_dt = TF(1) / tau_recycle;

        if (location == Lbc_location::West || location == Lbc_location::East)
        {
            for (int k=kstart; k<kend; ++k)
                for (int n=1; n<=nsponge; ++n)
                {
                    // Offset in y-direction for domain corners.
                    const int jstart_loc = mpiidy == 0     ? jstart + (n-1) : jstart;
                    const int jend_loc   = mpiidy == npy-1 ? jend   - (n-1) : jend;

                    for (int j=jstart_loc; j<jend_loc; ++j)
                    {
                        // Index in LBC:
                        const int jstride_lbc = igc+nsponge;
                        const int ilbc = (location==Lbc_location::West) ? igc+n-1 : nsponge-n;
                        const int ijk_lbc = ilbc + j*jstride_lbc + k*jstride_lbc*jcells;

                        const int i = (location==Lbc_location::West) ? istart+n-1 : iend-n;
                        const int ijk = i + j*icells + k*ijcells;

                        const TF a_diff = diffusion_3x3x3(
                                a, lbc[ijk_lbc], ijk, icells, ijcells);

                        // Nudge coefficient.
                        const TF f_sponge = (TF(1)+nsponge-(n+TF(0.5))) / nsponge;
                        const TF w1n = w_dt * f_sponge;
                        const TF w2n = w_diff * f_sponge;

                        // Apply nudge and sponge tendencies.
                        at[ijk] += w1n * (lbc[ijk_lbc]-a[ijk]);
                        at[ijk] -= w2n * a_diff;

                        // Turbulence recycling.
                        if (sw_recycle)
                        {
                            // Recycle strength; 0 at boundary, 1 at edge nudging zone.
                            const TF f_recycle = TF(1) - f_sponge;

                            // Source of recycling.
                            const int offset = (location == Lbc_location::West) ? recycle_offset : -recycle_offset;
                            const int ijko = (i+offset) + j*icells + k*ijcells;

                            TF a_mean = 0;
                            for (int jc=-3; jc<4; ++jc)
                                for (int ic=-3; ic<4; ++ic)
                                    a_mean += a[ijko + ic + jc*icells];
                            a_mean /= TF(49.);

                            at[ijk] += f_recycle * r_dt * ((a[ijko] - a_mean) - (a[ijk] - lbc[ijk_lbc]));
                        }
                    }
                }
        }
        if (location == Lbc_location::South || location == Lbc_location::North)
        {
            for (int k=kstart; k<kend; ++k)
                for (int n=1; n<=nsponge; ++n)
                {
                    const int istart_loc = (mpiidx == 0)     ? istart + n : istart;
                    const int iend_loc   = (mpiidx == npx-1) ? iend   - n : iend;

                    for (int i=istart_loc; i<iend_loc; ++i)
                    {
                        const int kstride_lbc = jgc + nsponge;
                        const int jlbc = (location==Lbc_location::South) ? jgc+n-1 : nsponge-n;
                        const int ijk_lbc = i + jlbc*icells + k*icells*kstride_lbc;

                        const int j = (location==Lbc_location::South) ? jstart+n-1 : jend-n;
                        const int ijk = i + j*icells + k*ijcells;

                        const TF a_diff = diffusion_3x3x3(
                                a, lbc[ijk_lbc], ijk, icells, ijcells);

                        // Nudge coefficient.
                        const TF f_sponge = (TF(1)+nsponge-(n+TF(0.5))) / nsponge;
                        const TF w1n = w_dt * f_sponge;
                        const TF w2n = w_diff * f_sponge;

                        at[ijk] += w1n*(lbc[ijk_lbc]-a[ijk]);
                        at[ijk] -= w2n*a_diff;

                        if (sw_recycle)
                        {
                            // Recycle strength; 0 at boundary, 1 at edge nudging zone.
                            const TF f_recycle = TF(1) - f_sponge;

                            // Source of recycling.
                            const int offset = (location == Lbc_location::South) ? recycle_offset : -recycle_offset;
                            const int ijko = i + (j+offset)*icells + k*ijcells;

                            TF a_mean = 0;
                            for (int jc=-3; jc<4; ++jc)
                                for (int ic=-3; ic<4; ++ic)
                                    a_mean += a[ijko + ic + jc*icells];
                            a_mean /= TF(49.);

                            at[ijk] += f_recycle * r_dt * ((a[ijko] - a_mean) - (a[ijk] - lbc[ijk_lbc]));
                        }
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
            const TF* const restrict fld_prev,
            const TF* const restrict fld_next,
            const int size,
            const TF f0)
    {
        const TF f1 = TF(1) - f0;
        for (int n=0; n<size; ++n)
            fld[n] = f0 * fld_prev[n] + f1 * fld_next[n];
    }


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


    template<typename TF>
    void set_lbc_gcs(
            TF* const restrict fld,
            const TF* const restrict lbc,
            const int ngc, const int nsponge,
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
            const int jstride_w = ngc + nsponge;
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
            const int jstride_e = ngc + nsponge;
            const int kstride_e = jstride_e * jcells;

            for (int k=kstart; k<kend; k++)
                for (int j=0; j<jcells; j++)
                    for (int i=0; i<ngc; i++)
                    {
                        const int ijk_in = (i+nsponge) + j*jstride_e + k*kstride_e;
                        const int ijk_out = (i+iend) + j*jstride_out + k*kstride_out;
                        fld[ijk_out] = lbc[ijk_in];
                    }
        }

        if (location == Lbc_location::South)
        {
            const int jstride_s = icells;
            const int kstride_s = jstride_s * (ngc + nsponge);

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
            const int jstride_n = icells;
            const int kstride_n = jstride_n * (ngc + nsponge);

            for (int k=kstart; k<kend; k++)
                for (int j=0; j<ngc; j++)
                    for (int i=istart; i<iend; i++)
                    {
                        const int ijk_in = i + (j+nsponge)*jstride_n + k*kstride_n;
                        const int ijk_out = i + (j+jend)*jstride_out + k*kstride_out;
                        fld[ijk_out] = lbc[ijk_in];
                    }
        }
    };


    template<typename TF, Lbc_location lbc_location>
    void calc_div_x(
            TF& div,
            const TF* const restrict lbc_u,
            const TF* const restrict rhoref,
            const TF* const restrict dz,
            const TF dy,
            const int nsponge,
            const int t,
            const int ngc, const int kgc,
            const int jstart, const int jend,
            const int ktot,
            const int jcells)
    {
        const int jstride_w = ngc + nsponge + 1;
        const int kstride_w = jstride_w * jcells;

        const int jstride_e = ngc + nsponge;
        const int kstride_e = jstride_e * jcells;

        const int iw = ngc;
        const int ie = nsponge;

        for (int k=0; k<ktot; ++k)
            for (int j=jstart; j<jend; ++j)
            {
                const int ijk_w = iw + j*jstride_w + k*kstride_w;
                const int ijk_e = ie + j*jstride_e + k*kstride_e;

                // Div = east-west.
                if (lbc_location == Lbc_location::East)
                    div += rhoref[k+kgc] * dy * dz[k+kgc] * lbc_u[ijk_e];
                else
                    div -= rhoref[k+kgc] * dy * dz[k+kgc] * lbc_u[ijk_w];
            }
    }

    template<typename TF, Lbc_location lbc_location>
    void calc_div_y(
            TF& div,
            const TF* const restrict lbc_v,
            const TF* const restrict rhoref,
            const TF* const restrict dz,
            const TF dx,
            const int nsponge,
            const int t,
            const int ngc, const int kgc,
            const int istart, const int iend,
            const int ktot,
            const int icells)
    {
        const int jstride_s = icells;
        const int kstride_s = jstride_s * (ngc + nsponge + 1);
        const int tstride_s = kstride_s * ktot;

        const int jstride_n = icells;
        const int kstride_n = jstride_n * (ngc + nsponge);
        const int tstride_n = kstride_n * ktot;

        const int js = ngc;
        const int jn = nsponge;

        for (int k=0; k<ktot; ++k)
            for (int i=istart; i<iend; ++i)
            {
                const int ijk_s = i + js*jstride_s + k*kstride_s;
                const int ijk_n = i + jn*jstride_n + k*kstride_n;

                // Div = north-south.
                if (lbc_location == Lbc_location::North)
                    div += rhoref[k+kgc] * dx * dz[k+kgc] * lbc_v[ijk_n];
                else
                    div -= rhoref[k+kgc] * dx * dz[k+kgc] * lbc_v[ijk_s];
            }
    }
}


template<typename TF>
Boundary_lateral<TF>::Boundary_lateral(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
        master(masterin), grid(gridin), fields(fieldsin), field3d_io(masterin, gridin)
{
    sw_openbc = inputin.get_item<bool>("boundary_lateral", "sw_openbc", "", false);

    if (sw_openbc)
    {
        sw_openbc_uv = inputin.get_item<bool>("boundary_lateral", "sw_openbc_uv", "", true);
        sw_openbc_w = inputin.get_item<bool>("boundary_lateral", "sw_openbc_w", "", false);
        sw_neumann_w = inputin.get_item<bool>("boundary_lateral", "sw_neumann_w", "", true);
        sw_wtop_2d = inputin.get_item<bool>("boundary_lateral", "sw_wtop_2d", "", false);
        slist = inputin.get_list<std::string>("boundary_lateral", "slist", "", std::vector<std::string>());

        // Check....
        if (sw_openbc_w && sw_neumann_w)
            throw std::runtime_error("Cant have both \"sw_openbc_w\" and \"sw_neumann_w\" = true!");

        sw_timedep = inputin.get_item<bool>("boundary_lateral", "sw_timedep", "", false);
        loadfreq = inputin.get_item<int>("boundary_lateral", "loadfreq", "");

        // Lateral sponge / diffusion layer.
        sw_sponge = inputin.get_item<bool>("boundary_lateral", "sw_sponge", "", false);
        if (sw_sponge)
        {
            n_sponge = inputin.get_item<int>("boundary_lateral", "n_sponge", "", 5);
            tau_sponge = inputin.get_item<TF>("boundary_lateral", "tau_sponge", "", 60);
            w_diff = inputin.get_item<TF>("boundary_lateral", "w_diff", "", 0.0033);
        }

        // Inflow perturbations
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

        // Turbulence recycling.
        sw_recycle = inputin.get_item<bool>("boundary_lateral", "sw_recycle", "", false);
        if (sw_recycle)
        {
            recycle_list = inputin.get_list<std::string>(
                "boundary_lateral", "recycle_list", "", std::vector<std::string>());
            tau_recycle = inputin.get_item<TF>("boundary_lateral", "tau_recycle", "");
            recycle_offset = inputin.get_item<int>("boundary_lateral", "recycle_offset", "");
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
    if (!sw_openbc)
        return;

    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    // Checks!
    if (sw_recycle)
    {
        if (!sw_sponge)
            throw std::runtime_error("Turbulence recycling only works combined with sw_sponge=1");

        if (n_sponge+recycle_offset > gd.imax+gd.igc)
            throw std::runtime_error("Turbulence recycling offset too large for domain decomposition in x-direction");
        if (n_sponge+recycle_offset > gd.jmax+gd.jgc)
            throw std::runtime_error("Turbulence recycling offset too large for domain decomposition in y-direction");
    }

    auto add_lbc_var = [&](
            Lbc_map<TF>& lbc_w_in,
            Lbc_map<TF>& lbc_e_in,
            Lbc_map<TF>& lbc_s_in,
            Lbc_map<TF>& lbc_n_in,
            const std::string& name)
    {
        const int igc_pad = (name == "u") ? gd.igc+1 : gd.igc;
        const int jgc_pad = (name == "v") ? gd.jgc+1 : gd.jgc;

        const int nlbc_w = igc_pad + n_sponge;
        const int nlbc_e = gd.igc + n_sponge;
        const int nlbc_s = jgc_pad + n_sponge;
        const int nlbc_n = gd.jgc + n_sponge;

        if (md.mpicoordx == 0)
            lbc_w_in.emplace(name, std::vector<TF>(nlbc_w * gd.kcells * gd.jcells));
        if (md.mpicoordx == md.npx-1)
            lbc_e_in.emplace(name, std::vector<TF>(nlbc_e * gd.kcells * gd.jcells));
        if (md.mpicoordy == 0)
            lbc_s_in.emplace(name, std::vector<TF>(gd.icells * nlbc_s * gd.kcells));
        if (md.mpicoordy == md.npy-1)
            lbc_n_in.emplace(name, std::vector<TF>(gd.icells * nlbc_n * gd.kcells));
    };

    auto add_lbcs = [&](
            Lbc_map<TF>& lbc_w_in,
            Lbc_map<TF>& lbc_e_in,
            Lbc_map<TF>& lbc_s_in,
            Lbc_map<TF>& lbc_n_in)
    {
        if (sw_openbc_uv)
        {
            add_lbc_var(lbc_w_in, lbc_e_in, lbc_s_in, lbc_n_in, "u");
            add_lbc_var(lbc_w_in, lbc_e_in, lbc_s_in, lbc_n_in, "v");
        }

        if (sw_openbc_w)
            add_lbc_var(lbc_w_in, lbc_e_in, lbc_s_in, lbc_n_in, "w");

        for (auto& fld : slist)
            add_lbc_var(lbc_w_in, lbc_e_in, lbc_s_in, lbc_n_in, fld);
    };

    // Create fixed/constant (or time interpolated in case of `sw_timedep`) LBC arrays.
    add_lbcs(lbc_w, lbc_e, lbc_s, lbc_n);

    if (sw_timedep)
    {
        // Add LBC arrays at previous and next time steps.
        add_lbcs(lbc_w_prev, lbc_e_prev, lbc_s_prev, lbc_n_prev);
        add_lbcs(lbc_w_next, lbc_e_next, lbc_s_next, lbc_n_next);
    }

    //// Make sure every MPI task has different seed.
    //if (sw_perturb)
    //{
    //    perturb_seed += master.get_mpiid();
    //    std::srand(perturb_seed);
    //}
}

template <typename TF>
void Boundary_lateral<TF>::read_lbc(
        TF& div_u, TF& div_v,
        Lbc_map<TF>& lbc_w_in,
        Lbc_map<TF>& lbc_e_in,
        Lbc_map<TF>& lbc_s_in,
        Lbc_map<TF>& lbc_n_in,
        const int time_index)
{
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


    auto print_minmax = [&](const std::vector<TF>& vec, const std::string& name)
    {
        TF min_val = 1e9;
        TF max_val = -1e9;

        for (auto & val : vec)
        {
            min_val = std::min(min_val, val);
            max_val = std::max(max_val, val);
        }

        printf(" - %s @ %d : min=%f, max=%f\n", name.c_str(), md.mpiid, min_val, max_val);
    };


    auto read_binary = [&](
            std::vector<TF>& vec,
            const std::string file_name,
            const unsigned long size,
            const unsigned long time_index)
    {
        FILE *pFile;
        pFile = fopen(file_name.c_str(), "rb");

        bool success = true;
        if (pFile == NULL)
            success = false;

        if (success)
        {
            // Jump to offset & read requested chunk.
            const size_t offset = time_index * size * sizeof(TF);
            fseek(pFile, offset, SEEK_SET);

            if (fread(vec.data(), sizeof(TF), size, pFile) != (unsigned)size)
                success = false;
        }

        if (!success)
        {
            #ifdef USEMPI
            std::cout << "SINGLE PROCESS EXCEPTION: reading binary failed." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            #else
            throw std::runtime_error("ERROR: reading binary failed");
            #endif
        }

        fclose(pFile);
    };


    auto copy_boundary = [&](
            std::vector<TF>& fld_out,
            const std::vector<TF>& fld_in,
            const int isize_in, const int jsize_in,
            const int isize_out, const int jsize_out,
            const int istart_in, const int jstart_in,
            const std::string& name, const std::string& loc)
    {
        // Copy LBC data from vector `fld_in`, which holds an entire
        // domain edge with data for all MPI tasks, to the `Lbc_map`
        // containing only data needed for the current MPI subdomain.

        const int size_in = gd.ktot * jsize_in * isize_in;
        const int size_out = gd.kcells * jsize_out * isize_out;

        const int jstride_in = isize_in;
        const int kstride_in = jstride_in * jsize_in;

        const int jstride_out = isize_out;
        const int kstride_out = jstride_out * jsize_out;

        for (int k=0; k<gd.ktot; k++)
            for (int j=0; j<jsize_out; j++)
                for (int i=0; i<isize_out; i++)
                {
                    const int ijk_in = i+istart_in + (j+jstart_in)*jstride_in + k*kstride_in;
                    const int ijk_out = i + j*jstride_out + (k+gd.kstart)*kstride_out;

                    fld_out[ijk_out] = fld_in[ijk_in];
                }
    };


    auto copy_boundaries = [&](const std::string& name)
    {
        // Read LBC data from binary files, and copy out the
        // part needed on the current MPI subdomain.

        // `u` at west boundary and `v` at south boundary also contain
        // `u` at `istart` and `v` at `jstart`, and are therefore larger.
        const int igc_pad = (name == "u") ? gd.igc+1 : gd.igc;
        const int jgc_pad = (name == "v") ? gd.jgc+1 : gd.jgc;

        // Number of ghost + sponge cells.
        const int nlbc_w = igc_pad + n_sponge;
        const int nlbc_e = gd.igc + n_sponge;
        const int nlbc_s = jgc_pad + n_sponge;
        const int nlbc_n = gd.jgc + n_sponge;

        const int ncells_w = gd.ktot * (gd.jtot+2*gd.jgc) * nlbc_w;
        const int ncells_e = gd.ktot * (gd.jtot+2*gd.jgc) * nlbc_e;
        const int ncells_s = gd.ktot * nlbc_s * (gd.itot+2*gd.igc);
        const int ncells_n = gd.ktot * nlbc_n * (gd.itot+2*gd.igc);

        // Arrays which hold data for the full domain edge,
        // i.e. for all MPI tasks. The `copy_boundary()` function
        // later copies out the local data needed on each MPI subdomain.
        std::vector<TF> lbc_w_full;
        std::vector<TF> lbc_e_full;
        std::vector<TF> lbc_s_full;
        std::vector<TF> lbc_n_full;

        if (md.mpicoordx == 0)
            lbc_w_full.resize(ncells_w);
        if (md.mpicoordx == md.npx-1)
            lbc_e_full.resize(ncells_e);
        if (md.mpicoordy == 0)
            lbc_s_full.resize(ncells_s);
        if (md.mpicoordy == md.npy-1)
            lbc_n_full.resize(ncells_n);


        if (md.mpicoordx == 0)
        {
            if (md.mpicoordy == 0)
                read_binary(lbc_w_full, "lbc_" + name + "_west.0000000", ncells_w, time_index);
            master.broadcast_y(lbc_w_full.data(), ncells_w, 0);

            copy_boundary(
                    lbc_w_in.at(name), lbc_w_full,
                    nlbc_w, gd.jtot+2*gd.jgc,
                    nlbc_w, gd.jcells,
                    0, md.mpicoordy*gd.jmax,
                    name, "west");

            // Calculate total inflow over west boundary.
            if (name == "u" && md.mpicoordy == 0)
                calc_div_x<TF, Lbc_location::West>(
                        div_u,
                        lbc_w_full.data(),
                        fields.rhoref.data(),
                        gd.dz.data(),
                        gd.dy,
                        n_sponge,
                        time_index,
                        gd.igc, gd.kgc,
                        gd.jgc, gd.jtot+gd.jgc,
                        gd.ktot, gd.jtot+(2*gd.jgc));
        }

        if (md.mpicoordx == md.npx-1)
        {
            if (md.mpicoordy == md.npy-1)
                read_binary(lbc_e_full, "lbc_" + name + "_east.0000000", ncells_e, time_index);
            master.broadcast_y(lbc_e_full.data(), ncells_e, md.npy-1);

            copy_boundary(
                    lbc_e_in.at(name), lbc_e_full,
                    nlbc_e, gd.jtot+2*gd.jgc,
                    nlbc_e, gd.jcells,
                    0, md.mpicoordy*gd.jmax,
                    name, "east");

            // Calculate total outflow over east boundary.
            if (name == "u" && md.mpicoordy == md.npy-1)
                calc_div_x<TF, Lbc_location::East>(
                        div_u,
                        lbc_e_full.data(),
                        fields.rhoref.data(),
                        gd.dz.data(),
                        gd.dy,
                        n_sponge,
                        time_index,
                        gd.igc, gd.kgc,
                        gd.jgc, gd.jtot+gd.jgc,
                        gd.ktot, gd.jtot+(2*gd.jgc));
	    }

        if (md.mpicoordy == 0)
	    {
            if (md.mpicoordx == md.npx-1)
                read_binary(lbc_s_full, "lbc_" + name + "_south.0000000", ncells_s, time_index);
            master.broadcast_x(lbc_s_full.data(), ncells_s, md.npx-1);

            copy_boundary(
                    lbc_s_in.at(name), lbc_s_full,
                    gd.itot+2*gd.igc, nlbc_s,
                    gd.icells, nlbc_s,
                    md.mpicoordx*gd.imax, 0,
                    name, "south");

            if (name == "v" && md.mpicoordx == md.npx-1)
                calc_div_y<TF, Lbc_location::South>(
                        div_v,
                        lbc_s_full.data(),
                        fields.rhoref.data(),
                        gd.dz.data(),
                        gd.dx,
                        n_sponge,
                        time_index,
                        gd.jgc, gd.kgc,
                        gd.igc, gd.itot+gd.igc,
                        gd.ktot, gd.itot+(2*gd.igc));
	    }

        if (md.mpicoordy == md.npy-1)
	    {
            if (md.mpicoordx == 0)
                read_binary(lbc_n_full, "lbc_" + name + "_north.0000000", ncells_n, time_index);
            master.broadcast_x(lbc_n_full.data(), ncells_n, 0);

            copy_boundary(
                    lbc_n_in.at(name), lbc_n_full,
                    gd.itot+2*gd.igc, nlbc_n,
                    gd.icells, nlbc_n,
                    md.mpicoordx*gd.imax, 0,
                    name, "north");

            if (name == "v" && md.mpicoordx == 0)
                calc_div_y<TF, Lbc_location::North>(
                        div_v,
                        lbc_n_full.data(),
                        fields.rhoref.data(),
                        gd.dz.data(),
                        gd.dx,
                        n_sponge,
                        time_index,
                        gd.jgc, gd.kgc,
                        gd.igc, gd.itot+gd.igc,
                        gd.ktot, gd.itot+(2*gd.igc));
	    }
    };

    if (sw_openbc_uv)
    {
        copy_boundaries("u");
        copy_boundaries("v");

        // Inflow is calculated at south+west edges, outflow at north+east edges.
        // Take sum, to get the net inflow in both directions.
        master.sum(&div_u, 1);
        master.sum(&div_v, 1);
    }

    if (sw_openbc_w)
        copy_boundaries("w");

    for (auto& fld : slist)
        copy_boundaries(fld);
}

template <typename TF>
void Boundary_lateral<TF>::create(
        Input& inputin,
        Timeloop<TF>& timeloop,
        const std::string& sim_name)
{
    if (!sw_openbc)
        return;

    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();
    TF* rhoref = fields.rhoref.data();

    // Domain total divergence in u and v direction.
    if (sw_openbc_uv)
    {
        w_top_2d.resize(gd.ijcells);

        if (sw_timedep)
        {
            w_top_2d_prev.resize(gd.ijcells);
            w_top_2d_next.resize(gd.ijcells);
        }
    }

    if (!sw_timedep)
    {
        // Read LBC data directly into `lbc_{w/e/n/s}`.
        TF div_u = 0;
        TF div_v = 0;
        const int time_index = 0;
        read_lbc(div_u, div_v, lbc_w, lbc_e, lbc_s, lbc_n, time_index);

        if (sw_wtop_2d)
            read_xy_slice(w_top_2d, "w_top", 0);
        else
        {
            const TF w_top_mean = -(div_u + div_v) / (fields.rhorefh[gd.kend] * gd.xsize * gd.ysize);
            std::fill(w_top_2d.begin(), w_top_2d.end(), w_top_mean);

            std::string message =
                    "- div(u) = " + std::to_string(div_u)
                  + ", div(v) = " + std::to_string(div_v)
                  + ", w_top = " + std::to_string(w_top_mean*100) + " cm/s";
            master.print_message(message);
        }
    }
    else
    {
        // Read previous and next input times.
        const double time = timeloop.get_time();
        const unsigned long itime = timeloop.get_itime();
        unsigned long iiotimeprec = timeloop.get_iiotimeprec();
        unsigned long iloadtime = convert_to_itime(loadfreq);

        // Determine time index for LBCs.
        const int prev_index = itime / iloadtime;
        const int next_index = prev_index + 1;

        // Determine `iotime` for `w_top` fields.
        prev_itime = prev_index * iloadtime;
        next_itime = next_index * iloadtime;

        unsigned long prev_iotime = prev_index * iloadtime / iiotimeprec;
        unsigned long next_iotime = next_index * iloadtime / iiotimeprec;

        // Read previous and next LBC values.
        TF div_u_prev = 0;
        TF div_v_prev = 0;

        TF div_u_next = 0;
        TF div_v_next = 0;

        read_lbc(div_u_prev, div_v_prev, lbc_w_prev, lbc_e_prev, lbc_s_prev, lbc_n_prev, prev_index);
        read_lbc(div_u_next, div_v_next, lbc_w_next, lbc_e_next, lbc_s_next, lbc_n_next, next_index);

        if (sw_wtop_2d)
        {
            read_xy_slice(w_top_2d_prev, "w_top", prev_iotime);
            read_xy_slice(w_top_2d_next, "w_top", next_iotime);
        }
        else
        {
            w_top_prev = -(div_u_prev + div_v_prev) / (fields.rhorefh[gd.kend] * gd.xsize * gd.ysize);
            w_top_next = -(div_u_next + div_v_next) / (fields.rhorefh[gd.kend] * gd.xsize * gd.ysize);

            std::string message1 =
                    "- div(u) = " + std::to_string(div_u_prev)
                    + ", div(v) = " + std::to_string(div_v_prev)
                    + ", w_top = " + std::to_string(w_top_prev*100)
                    + " cm/s @ t= " + std::to_string(prev_index*loadfreq) + " sec.";
            master.print_message(message1);

            std::string message2 =
                    "- div(u) = " + std::to_string(div_u_next)
                    + ", div(v) = " + std::to_string(div_v_next)
                    + ", w_top = " + std::to_string(w_top_next*100)
                    + " cm/s @ t= " + std::to_string(next_index*loadfreq) + " sec.";
            master.print_message(message2);
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
    if (!sw_openbc)
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
                ngc, n_sponge,
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

    if (sw_openbc_uv)
    {
        set_gcs("u");
        set_gcs("v");
    }

    if (sw_openbc_w)
        set_gcs("w");

    for (auto& fld : slist)
        set_gcs(fld);


    // Set vertical velocity at domain top.
    if (sw_openbc_uv)
    {
        const int k = gd.kend;
        for (int j=gd.jstart; j<gd.jend; ++j)
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*gd.icells + k*gd.ijcells;
                const int ij = i + j*gd.icells;
                fields.mp.at("w")->fld[ijk] = w_top_2d[ij];
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


    auto sponge_layer_wrapper = [&]<Lbc_location location, bool sw_recycle>(
            std::map<std::string, std::vector<TF>>& lbc_map,
            const std::string& name)
    {
        if (!sw_sponge)
            return;

        const int kstart = (name == "w") ? gd.kstart+1 : gd.kstart;

        lateral_sponge_kernel_s<TF, location, sw_recycle>(
                fields.at.at(name)->fld.data(),
                fields.ap.at(name)->fld.data(),
                lbc_map.at(name).data(),
                tau_sponge,
                w_diff,
                n_sponge,
                tau_recycle,
                recycle_offset,
                md.npx, md.npy,
                md.mpicoordx, md.mpicoordy,
                gd.igc, gd.jgc,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                kstart, gd.kend,
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
                tau_sponge,
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


    auto sponge_layer_v_wrapper = [&]<Lbc_location location>(
            std::map<std::string, std::vector<TF>>& lbc_map)
    {
        if (!sw_sponge)
            return;

        lateral_sponge_kernel_v<TF, location>(
                fields.mt.at("v")->fld.data(),
                fields.mp.at("v")->fld.data(),
                lbc_map.at("v").data(),
                tau_sponge,
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

    if (sw_openbc_uv)
    {
        // NOTE: order of calls here is important; the `_u` and `_v` wrappers don't take
        //       the offset in corners into account, so always first call the `_u` and `_v`
        //       wrappers, followed by the generic `sponge_layer_wrapper()`.

        if (md.mpicoordx == 0)
            sponge_layer_u_wrapper.template operator()<Lbc_location::West>(lbc_w);
        if (md.mpicoordx == md.npx-1)
            sponge_layer_u_wrapper.template operator()<Lbc_location::East>(lbc_e);
        if (md.mpicoordy == 0)
            sponge_layer_wrapper.template operator()<Lbc_location::South, false>(lbc_s, "u");
        if (md.mpicoordy == md.npy-1)
            sponge_layer_wrapper.template operator()<Lbc_location::North, false>(lbc_n, "u");

        if (md.mpicoordy == 0)
            sponge_layer_v_wrapper.template operator()<Lbc_location::South>(lbc_s);
        if (md.mpicoordy == md.npy-1)
            sponge_layer_v_wrapper.template operator()<Lbc_location::North>(lbc_n);
        if (md.mpicoordx == 0)
            sponge_layer_wrapper.template operator()<Lbc_location::West, false>(lbc_w, "v");
        if (md.mpicoordx == md.npx-1)
            sponge_layer_wrapper.template operator()<Lbc_location::East, false>(lbc_e, "v");
    }

    if (sw_openbc_w)
    {
        if (md.mpicoordx == 0)
            sponge_layer_wrapper.template operator()<Lbc_location::West, false>(lbc_w, "w");
        if (md.mpicoordx == md.npx-1)
            sponge_layer_wrapper.template operator()<Lbc_location::East, false>(lbc_e, "w");
        if (md.mpicoordy == 0)
            sponge_layer_wrapper.template operator()<Lbc_location::South, false>(lbc_s, "w");
        if (md.mpicoordy == md.npy-1)
            sponge_layer_wrapper.template operator()<Lbc_location::North, false>(lbc_n, "w");
    }

    for (auto& fld : slist)
    {
        const bool sw_recycle = in_list<std::string>(fld, recycle_list);

        if (md.mpicoordx == 0)
        {
            if (sw_recycle)
                sponge_layer_wrapper.template operator()<Lbc_location::West, true>(lbc_w, fld);
            else
                sponge_layer_wrapper.template operator()<Lbc_location::West, false>(lbc_w, fld);
        }
        if (md.mpicoordx == md.npx-1)
        {
            if (sw_recycle)
                sponge_layer_wrapper.template operator()<Lbc_location::East, true>(lbc_e, fld);
            else
                sponge_layer_wrapper.template operator()<Lbc_location::East, false>(lbc_e, fld);
        }
        if (md.mpicoordy == 0)
        {
            if (sw_recycle)
                sponge_layer_wrapper.template operator()<Lbc_location::South, true>(lbc_s, fld);
            else
                sponge_layer_wrapper.template operator()<Lbc_location::South, false>(lbc_s, fld);
        }
        if (md.mpicoordy == md.npy-1)
        {
            if (sw_recycle)
                sponge_layer_wrapper.template operator()<Lbc_location::North, true>(lbc_n, fld);
            else
                sponge_layer_wrapper.template operator()<Lbc_location::North, false>(lbc_n, fld);
        }
    }

    //dump_vector(fields.ap.at("u")->fld, "u");
    //dump_vector(fields.ap.at("v")->fld, "v");
    //dump_vector(fields.ap.at("w")->fld, "w");
    //dump_vector(fields.ap.at("th")->fld, "th");

    //dump_vector(fields.at.at("u")->fld, "ut");
    //dump_vector(fields.at.at("v")->fld, "vt");
    //dump_vector(fields.at.at("w")->fld, "wt");
    //dump_vector(fields.at.at("th")->fld, "tht");

    //throw 1;

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
    if (!sw_openbc || !sw_timedep)
        return;

    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    // Find index in time array.
    double time = timeloop.get_time();
    double endtime = timeloop.get_endtime();
    unsigned long itime = timeloop.get_itime();
    unsigned long iiotimeprec = timeloop.get_iiotimeprec();
    unsigned long iloadtime = convert_to_itime(loadfreq);
    unsigned long iendtime = convert_to_itime(endtime);

    // CvH: this is an UGLY hack, because it only works for RK4.
    // We need to know from the time whether we are in the last iter.
    // Also, it will fail miserably if we perturb the velocities at the walls
    if (pres_fix && timeloop.get_substep() == 4)
    {
        time += timeloop.get_dt();
        itime += timeloop.get_idt();
    }

    if (itime >= next_itime)
    {
        // BvS: We are setting `w_top` for the next time step;
        //      skip if next time is beyond the endtime.
        if (next_itime + iloadtime > iendtime)
        {
            master.print_warning("Timedep boundary_lateral out-of-bounds at t+dt, skipping update.\n");
            return;
        }

        // Advance time and read new files.
        prev_itime = next_itime;
        next_itime = prev_itime + iloadtime;

        unsigned long next_iotime = next_itime / iiotimeprec;

        // Move LBCs from next to previous values.
        for (auto& it : lbc_w_next)
            lbc_w_prev.at(it.first) = it.second;
        for (auto& it : lbc_e_next)
            lbc_e_prev.at(it.first) = it.second;
        for (auto& it : lbc_s_next)
            lbc_s_prev.at(it.first) = it.second;
        for (auto& it : lbc_n_next)
            lbc_n_prev.at(it.first) = it.second;

        // Read new LBC values.
        const int next_index = next_itime / iloadtime;

        TF div_u_next = 0;
        TF div_v_next = 0;
        read_lbc(div_u_next, div_v_next, lbc_w_next, lbc_e_next, lbc_s_next, lbc_n_next, next_index);

        // Read in or calculate new `w_top`.
        if (sw_wtop_2d)
        {
            w_top_2d_prev = w_top_2d_next;
            read_xy_slice(w_top_2d_next, "w_top", next_iotime);
        }
        else
        {
            w_top_prev = w_top_next;
            w_top_next = -(div_u_next + div_v_next) / (fields.rhorefh[gd.kend] * gd.xsize * gd.ysize);

            std::string message2 =
                    "- div(u) = " + std::to_string(div_u_next)
                    + ", div(v) = " + std::to_string(div_v_next)
                    + ", w_top = " + std::to_string(w_top_next*100)
                    + " cm/s @ t= " + std::to_string(next_index*loadfreq) + " sec.";
            master.print_message(message2);
        }
    }

    // Interpolate LBCs and w_top to current time.
    const TF f0 = TF(1) - ((itime - prev_itime) / TF(iloadtime));
    const TF f1 = TF(1) - f0;

    // Interpolate mean domain top velocity
    if (sw_wtop_2d)
    {
        // Interpolate `w_top` field in time.
        for (int n=0; n<gd.ijcells; ++n)
            w_top_2d[n] = f0 * w_top_2d_prev[n] + f1 * w_top_2d_next[n];
    }
    else
    {
        const TF w_top = f0 * w_top_prev + f1 * w_top_next;
        std::fill(w_top_2d.begin(), w_top_2d.end(), w_top);
    }

    // Interpolate boundaries in time.
    if (md.mpicoordx == 0)
    {
        for (auto& it : lbc_w)
        {
            const int ngc = (it.first == "u") ? gd.igc+1 : gd.igc;
            interpolate_lbc_kernel(
                    lbc_w.at(it.first).data(),
                    lbc_w_prev.at(it.first).data(),
                    lbc_w_next.at(it.first).data(),
                    gd.kcells * gd.jcells * (ngc+n_sponge),
                    f0);
        }
    }

    if (md.mpicoordx == md.npx-1)
    {
        for (auto& it : lbc_e)
            interpolate_lbc_kernel(
                    lbc_e.at(it.first).data(),
                    lbc_e_prev.at(it.first).data(),
                    lbc_e_next.at(it.first).data(),
                    gd.kcells * gd.jcells * (gd.igc+n_sponge),
                    f0);
    }

    if (md.mpicoordy == 0)
    {
        for (auto& it : lbc_s)
        {
            const int ngc = (it.first == "v") ? gd.jgc+1 : gd.jgc;
            interpolate_lbc_kernel(
                    lbc_s.at(it.first).data(),
                    lbc_s_prev.at(it.first).data(),
                    lbc_s_next.at(it.first).data(),
                    gd.kcells * (ngc+n_sponge) * gd.icells,
                    f0);
        }
    }

    if (md.mpicoordy == md.npy-1)
    {
        for (auto& it : lbc_n)
            interpolate_lbc_kernel(
                    lbc_n.at(it.first).data(),
                    lbc_n_prev.at(it.first).data(),
                    lbc_n_next.at(it.first).data(),
                    gd.kcells * (gd.jgc+n_sponge) * gd.icells,
                    f0);
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
