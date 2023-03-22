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

        sw_sponge = inputin.get_item<bool>("boundary", "sw_sponge", "", true);
        n_sponge = inputin.get_item<int>("boundary", "n_sponge", "", 5);

        tau_nudge = inputin.get_item<TF>("boundary", "tau_sponge", "", 60);
        w_diff = inputin.get_item<TF>("boundary", "w_diff", "", 0.0033);
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

    auto copy_boundaries = [&](
                const std::string& name)
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

        if (!sw_timedep)
        {
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
        }
    };

    if (sw_inoutflow_u)
        copy_boundaries("u");
    if (sw_inoutflow_v)
        copy_boundaries("v");

    for (auto& fld : inoutflow_s)
        copy_boundaries(fld);
}

template <typename TF>
void Boundary_lateral<TF>::set_ghost_cells()
{
    if (!sw_inoutflow)
        return;

    auto& gd = grid.get_grid_data();
    auto& md = master.get_MPI_data();

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

    if (sw_inoutflow_u)
    {
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
        //    const TF hor_div = fields.rhoref[k]*(TF(0.) - TF(2.)) / gd.xsize // * gd.z[k]/gd.zsize
        //                     + fields.rhoref[k]*(TF(2.) - TF(0.)) / gd.ysize; // * gd.z[k]/gd.zsize;

        //    w_top = (fields.rhorefh[k]*w_top - gd.dz[k]*hor_div) / fields.rhorefh[k+1];
        //}

        //for (int j=gd.jstart; j<gd.jend; ++j)
        //    for (int i=gd.istart; i<gd.iend; ++i)
        //    {
        //        const int ijk = i + j*jstride + gd.kend*kstride;
        //        fields.mp.at("w")->fld[ijk] = w_top;
        //    }

        const int ii = 1;
        const int jj = gd.icells;
        const int kk = gd.ijcells;

        TF* u = fields.mp.at("u")->fld.data();
        TF* v = fields.mp.at("v")->fld.data();
        TF* w = fields.mp.at("w")->fld.data();

        TF* rho  = fields.rhoref.data();
        TF* rhoh = fields.rhorefh.data();

        const TF dxi = TF(1)/gd.dx;
        const TF dyi = TF(1)/gd.dy;

        const int k = gd.kend-1;
        for (int j=gd.jstart; j<gd.jend; ++j)
            for (int i=gd.istart; i<gd.iend; ++i)
            {
                const int ijk = i + j*jj + k*kk;
                w[ijk+kk] = -(rho[k] * ((u[ijk+ii] - u[ijk]) * dxi + (v[ijk+jj] - v[ijk]) * dyi) * gd.dz[k] - rhoh[k] * w[ijk]) / rhoh[k+1];
            }
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
    }
}

template <typename TF>
void Boundary_lateral<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (!sw_inoutflow)
        return;
}

template class Boundary_lateral<double>;
template class Boundary_lateral<float>;
