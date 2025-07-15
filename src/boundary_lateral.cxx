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
#include <cmath>

#include "boundary_lateral.h"
#include "boundary_lateral_kernels.h"

#include "netcdf_interface.h"
#include "grid.h"
#include "fields.h"
#include "input.h"
#include "master.h"
#include "timeloop.h"
#include "stats.h"
#include "constants.h"

namespace blk = boundary_lateral_kernels;

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

    bool any_true(std::map<Lbc_location, bool>& map_in)
    {
        return std::any_of(map_in.begin(), map_in.end(), [](const auto& p) { return p.second; });
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

    template<typename TF>
    bool is_equal(const TF a, const TF b)
    {
        const TF epsilon = std::max(TF(10) * std::numeric_limits<TF>::epsilon(), std::max(std::abs(a), std::abs(b)) * TF(10) * std::numeric_limits<TF>::epsilon());
        return std::abs(a - b) <= epsilon;
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

        // Turbulence recycling.
        sw_recycle[Lbc_location::North] = inputin.get_item<bool>("boundary_lateral", "sw_recycle", "north", false);
        sw_recycle[Lbc_location::East]  = inputin.get_item<bool>("boundary_lateral", "sw_recycle", "east", false);
        sw_recycle[Lbc_location::South] = inputin.get_item<bool>("boundary_lateral", "sw_recycle", "south", false);
        sw_recycle[Lbc_location::West]  = inputin.get_item<bool>("boundary_lateral", "sw_recycle", "west", false);

        if (any_true(sw_recycle))
        {
            recycle_list = inputin.get_list<std::string>(
                "boundary_lateral", "recycle_list", "", std::vector<std::string>());
            tau_recycle = inputin.get_item<TF>("boundary_lateral", "tau_recycle", "");
            recycle_offset = inputin.get_item<int>("boundary_lateral", "recycle_offset", "");
        }

        //sw_recycle = inputin.get_item<bool>("boundary_lateral", "sw_recycle", "", false);
        //if (sw_recycle)
        //{
        //    recycle_list = inputin.get_list<std::string>(
        //        "boundary_lateral", "recycle_list", "", std::vector<std::string>());
        //    tau_recycle = inputin.get_item<TF>("boundary_lateral", "tau_recycle", "");
        //    recycle_offset = inputin.get_item<int>("boundary_lateral", "recycle_offset", "");
        //}
    }

    // Output for sub-domain.
    // Keep out of the `if (sw_openbc)` block; domains can have periodic
    // boundaries, but still need to output sub-domain data.
    sw_subdomain = inputin.get_item<bool>("subdomain", "sw_subdomain", "", false);

    if (sw_subdomain)
    {
        xstart_sub = inputin.get_item<TF>("subdomain", "xstart", "");
        xend_sub   = inputin.get_item<TF>("subdomain", "xend", "");
        ystart_sub = inputin.get_item<TF>("subdomain", "ystart", "");
        yend_sub   = inputin.get_item<TF>("subdomain", "yend", "");

        grid_ratio_sub = inputin.get_item<int>("subdomain", "grid_ratio", "");
        n_ghost_sub = inputin.get_item<int>("subdomain", "n_ghost", "");
        n_sponge_sub = inputin.get_item<int>("subdomain", "n_sponge", "");
        savetime_sub = inputin.get_item<int>("subdomain", "savetime", "");
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
    if (any_true(sw_recycle))
    {
        if (!sw_sponge)
            throw std::runtime_error("Turbulence recycling only works combined with sw_sponge=1");

        if ((sw_recycle[Lbc_location::West] || sw_recycle[Lbc_location::East])  && recycle_offset > gd.imax+gd.igc)
            throw std::runtime_error("Turbulence recycling offset too large for domain decomposition in x-direction");
        if ((sw_recycle[Lbc_location::North] || sw_recycle[Lbc_location::South])  && recycle_offset > gd.jmax+gd.jgc)
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
}


template <typename TF>
void Boundary_lateral<TF>::read_lbc(
        TF& div_u, TF& div_v,
        Lbc_map<TF>& lbc_w_in,
        Lbc_map<TF>& lbc_e_in,
        Lbc_map<TF>& lbc_s_in,
        Lbc_map<TF>& lbc_n_in,
        const int iotime)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    //auto dump_vector = [&](
    //        std::vector<TF>& fld,
    //        const std::string& name)
    //{
    //    std::string name_out = name + "." + std::to_string(md.mpicoordx) + "." + std::to_string(md.mpicoordy) + ".bin";

    //    FILE *pFile;
    //    pFile = fopen(name_out.c_str(), "wb");

    //    if (pFile == NULL)
    //        throw std::runtime_error("Opening raw dump field failed.");

    //    fwrite(fld.data(), sizeof(TF), fld.size(), pFile);
    //    fclose(pFile);
    //};


    auto read_binary = [&](
            std::vector<TF>& vec,
            const std::string name,
            const unsigned long size)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", name.c_str(), iotime);

        master.print_message("Loading \"%s\" ... ", filename);

        FILE *pFile;
        pFile = fopen(filename, "rb");

        bool success = true;
        if (pFile == NULL)
            success = false;

        if (success)
        {
            // Jump to offset & read requested chunk.
            //const size_t offset = time_index * size * sizeof(TF);
            //fseek(pFile, offset, SEEK_SET);

            if (fread(vec.data(), sizeof(TF), size, pFile) != (unsigned)size)
                success = false;
        }

        if (success)
            master.print_message("OK\n");
        else
        {
            master.print_message("FAILED\n");

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
                read_binary(lbc_w_full, "lbc_" + name + "_west", ncells_w);
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
                        gd.igc, gd.kgc,
                        gd.jgc, gd.jtot+gd.jgc,
                        gd.ktot, gd.jtot+(2*gd.jgc));
        }

        if (md.mpicoordx == md.npx-1)
        {
            if (md.mpicoordy == md.npy-1)
                read_binary(lbc_e_full, "lbc_" + name + "_east", ncells_e);
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
                        gd.igc, gd.kgc,
                        gd.jgc, gd.jtot+gd.jgc,
                        gd.ktot, gd.jtot+(2*gd.jgc));
	    }

        if (md.mpicoordy == 0)
	    {
            if (md.mpicoordx == md.npx-1)
                read_binary(lbc_s_full, "lbc_" + name + "_south", ncells_s);
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
                        gd.jgc, gd.kgc,
                        gd.igc, gd.itot+gd.igc,
                        gd.ktot, gd.itot+(2*gd.igc));
	    }

        if (md.mpicoordy == md.npy-1)
	    {
            if (md.mpicoordx == 0)
                read_binary(lbc_n_full, "lbc_" + name + "_north", ncells_n);
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
        Stats<TF>& stats,
        const std::string& sim_name)
{
    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    if (sw_subdomain)
    {
        bool error = false;

        auto check_subdomain = [&](const TF value, const TF spacing, const std::string& name)
        {
            if (!is_equal(std::fmod(value, spacing), TF(0)))
            {
                error = true;
                master.print_message(
                    "ERROR: " + name + " must be an integer multiple of " + std::to_string(spacing) + ".");
            }
        };

        // Sub-domain has perfectly align with parent grid at x/y half levels.
        check_subdomain(xstart_sub, gd.dx, "xstart_sub");
        check_subdomain(xend_sub, gd.dx, "xend_sub");
        check_subdomain(ystart_sub, gd.dy, "ystart_sub");
        check_subdomain(yend_sub, gd.dy, "yend_sub");

        if (error)
            throw std::runtime_error("Sub-domain boundaries not aligned with parent grid.");

        // Initialise LBCs instance.
        xsize_sub = xend_sub - xstart_sub;
        ysize_sub = yend_sub - ystart_sub;

        itot_sub = static_cast<int>(xsize_sub / gd.dx * grid_ratio_sub + 0.5);
        jtot_sub = static_cast<int>(ysize_sub / gd.dy * grid_ratio_sub + 0.5);

        dx_sub = xsize_sub / itot_sub;
        dy_sub = ysize_sub / jtot_sub;
    }

    // Only proceed if open boundary conditions are enabled.
    if (!sw_openbc)
        return;

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
        const int iotime = 0;
        read_lbc(div_u, div_v, lbc_w, lbc_e, lbc_s, lbc_n, iotime);

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

        read_lbc(div_u_prev, div_v_prev, lbc_w_prev, lbc_e_prev, lbc_s_prev, lbc_n_prev, prev_iotime);
        read_lbc(div_u_next, div_v_next, lbc_w_next, lbc_e_next, lbc_s_next, lbc_n_next, next_iotime);

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

    if (sw_sponge)
    {
        stats.add_tendency(*fields.mt.at("u"), "z", tend_name, tend_longname);
        stats.add_tendency(*fields.mt.at("v"), "z", tend_name, tend_longname);
        stats.add_tendency(*fields.mt.at("w"), "zh", tend_name, tend_longname);

        for (auto& fld : slist)
            stats.add_tendency(*fields.at.at(fld), "zh", tend_name, tend_longname);
    }
}


template<typename TF>
unsigned long Boundary_lateral<TF>::get_time_limit(unsigned long itime)
{
    unsigned long idtlim = Constants::ulhuge;

    if (sw_openbc && sw_timedep)
    {
        const unsigned long ifreq = convert_to_itime(loadfreq);
        idtlim = std::min(idtlim, ifreq - itime % ifreq);
    }

    if (sw_subdomain)
    {
        const unsigned long ifreq = convert_to_itime(savetime_sub);
        idtlim = std::min(idtlim, ifreq - itime % ifreq);
    }

    return idtlim;
}


template <typename TF>
void Boundary_lateral<TF>::set_ghost_cells(
        Timeloop<TF>& timeloop)
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

    //dump_vector(fields.ap.at("u")->fld, "u");
    //dump_vector(fields.ap.at("v")->fld, "v");
    //dump_vector(fields.ap.at("w")->fld, "w");
    //dump_vector(fields.ap.at("th")->fld, "th");

    //dump_vector(fields.at.at("u")->fld, "ut");
    //dump_vector(fields.at.at("v")->fld, "vt");
    //dump_vector(fields.at.at("w")->fld, "wt");
    //dump_vector(fields.at.at("th")->fld, "tht");

    //throw 1;
}


template <typename TF>
void Boundary_lateral<TF>::exec_lateral_sponge(
        Stats<TF>& stats)
{
    if (!sw_openbc or !sw_sponge)
        return;

    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    auto sponge_layer_wrapper = [&]<Lbc_location location, bool sw_recycle>(
            std::map<std::string, std::vector<TF>>& lbc_map,
            const std::string& name)
    {
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

        stats.calc_tend(*fields.mt.at("u"), tend_name);
        stats.calc_tend(*fields.mt.at("v"), tend_name);
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

        stats.calc_tend(*fields.mt.at("w"), tend_name);
    }

    for (auto& fld : slist)
    {
        const bool sw_recycle_fld = in_list<std::string>(fld, recycle_list);

        if (md.mpicoordx == 0)
        {
            if (sw_recycle_fld && sw_recycle[Lbc_location::West])
                sponge_layer_wrapper.template operator()<Lbc_location::West, true>(lbc_w, fld);
            else
                sponge_layer_wrapper.template operator()<Lbc_location::West, false>(lbc_w, fld);
        }
        if (md.mpicoordx == md.npx-1 && sw_recycle[Lbc_location::East])
        {
            if (sw_recycle_fld && sw_recycle[Lbc_location::East])
                sponge_layer_wrapper.template operator()<Lbc_location::East, true>(lbc_e, fld);
            else
                sponge_layer_wrapper.template operator()<Lbc_location::East, false>(lbc_e, fld);
        }
        if (md.mpicoordy == 0 && sw_recycle[Lbc_location::South])
        {
            if (sw_recycle_fld && sw_recycle[Lbc_location::South])
                sponge_layer_wrapper.template operator()<Lbc_location::South, true>(lbc_s, fld);
            else
                sponge_layer_wrapper.template operator()<Lbc_location::South, false>(lbc_s, fld);
        }
        if (md.mpicoordy == md.npy-1 && sw_recycle[Lbc_location::North])
        {
            if (sw_recycle_fld && sw_recycle[Lbc_location::North])
                sponge_layer_wrapper.template operator()<Lbc_location::North, true>(lbc_n, fld);
            else
                sponge_layer_wrapper.template operator()<Lbc_location::North, false>(lbc_n, fld);
        }

        stats.calc_tend(*fields.at.at(fld), tend_name);
    }
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
            master.print_warning("Timedep boundary_lateral out-of-bounds at t+dt, skipping update.\n");
        else
        {
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
            read_lbc(div_u_next, div_v_next, lbc_w_next, lbc_e_next, lbc_s_next, lbc_n_next, next_iotime);

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


namespace
{
    template <typename TF>
    void fetch_lbcs_scalar(
        TF* const restrict lbc,
        const TF* const restrict fld,
        const int iref_g,       // Base index in global domain
        const int jref_g,
        const int istart_s,     // Start index in sub-domain
        const int jstart_s,
        const int kstart_s,
        const int iend_s,       // End index in sub-domain
        const int jend_s,
        const int kend_s,
        const int jstride_g,    // Stride in global domain
        const int kstride_g,
        const int jstride_lbc,  // Stride in lbc field.
        const int kstride_lbc,
        const int n_ghost,
        const int n_sponge,
        const int grid_ratio,
        const int kgc)
    {
        for (int k=kstart_s; k<kend_s; ++k)
            for (int j=jstart_s; j<jend_s; ++j)
            {
                const int j_lbc = j + n_ghost;
                const int j_g = jref_g + int(std::floor(TF(j) / grid_ratio));     // NN-index

                for (int i=istart_s; i<iend_s; ++i)
                {
                    const int i_lbc = i + n_ghost;
                    const int i_g = iref_g + int(std::floor(TF(i) / grid_ratio));    // NN-index

                    //std::cout << "i=" << i << " i_lbc=" << i_lbc << " ig=" << i_g << std::endl;

                    const int ijk_lbc = i_lbc + j_lbc * jstride_lbc + (k-kgc) * kstride_lbc;
                    const int ijk_g = i_g + j_g * jstride_g + k * kstride_g;

                    lbc[ijk_lbc] = fld[ijk_g];
                }
            }
    }
}



template <typename TF>
void Boundary_lateral<TF>::save_lbcs(
        Timeloop<TF>& timeloop)
{
    if (!sw_subdomain || timeloop.in_substep() || timeloop.get_itime() % convert_to_itime(savetime_sub) != 0)
        return;

    std::string msg = "Saving sub-domain LBCs for time " + std::to_string(timeloop.get_time()) + " ...";
    master.print_message(msg);

    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

    auto save_binary = [&](
            std::vector<TF>& fld,
            const std::string& name,
            const int time)
    {
        char filename[256];
        std::sprintf(filename, "%s.%07d", name.c_str(), time);
        master.print_message("Saving \"%s\" ... ", filename);

        FILE *pFile;
        pFile = fopen(filename, "wb");

        if (pFile == NULL)
            master.print_message("FAILED\n");
        else
            master.print_message("OK\n");

        fwrite(fld.data(), sizeof(TF), fld.size(), pFile);
        fclose(pFile);
    };
}

#ifdef FLOAT_SINGLE
template class Boundary_lateral<float>;
#else
template class Boundary_lateral<double>;
#endif
