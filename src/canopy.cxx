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

#include <algorithm>
#include <cmath>
#include  <stdexcept>

#include "canopy.h"
#include "input.h"
#include "field3d_operators.h"
#include "grid.h"
#include "netcdf_interface.h"
#include "fast_math.h"
#include "fields.h"
#include "master.h"

namespace
{
    namespace fm = Fast_math;

    template<typename TF, bool sw_pad_3d>
    void canopy_drag_u(
            TF* const restrict ut,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict pad,
            const TF utrans,
            const TF vtrans,
            const TF cd,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        const int ii = 1;
        const int jj = jstride;
        const int kk = kstride;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // Interpolate `pad` to `u` locations.
                    const TF pad_u = sw_pad_3d ? TF(0.5) * (pad[ijk-ii] + pad[ijk]) : pad[k];

                    // Interpolate `v` and `w` to `u` locations.
                    const TF u_on_u = u[ijk] + utrans;
                    const TF v_on_u = TF(0.25) * (v[ijk] + v[ijk-ii] + v[ijk-ii+jj] + v[ijk+jj]) + vtrans;
                    const TF w_on_u = TF(0.25) * (w[ijk] + w[ijk-ii] + w[ijk-ii+kk] + w[ijk+kk]);

                    const TF ftau = -cd * pad_u *
                        std::pow( fm::pow2(u_on_u) +
                                  fm::pow2(v_on_u) +
                                  fm::pow2(w_on_u), TF(0.5) );

                    ut[ijk] += ftau * u_on_u;
                }
    }

    template<typename TF, bool sw_pad_3d>
    void canopy_drag_v(
            TF* const restrict vt,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict pad,
            const TF utrans,
            const TF vtrans,
            const TF cd,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        const int ii = 1;
        const int jj = jstride;
        const int kk = kstride;

        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // Interpolate `pad` to `v` locations.
                    const TF pad_v = sw_pad_3d ? TF(0.5) * (pad[ijk-jj] + pad[ijk]) : pad[k];

                    // Interpolate `u` and `w` to `v` locations.
                    const TF u_on_v = TF(0.25) * (u[ijk] + u[ijk+ii] + u[ijk+ii-jj] + u[ijk-jj]) + utrans;
                    const TF v_on_v = v[ijk] + vtrans;
                    const TF w_on_v = TF(0.25) * (w[ijk] + w[ijk+kk] + w[ijk-jj] + w[ijk-jj+kk]);

                    const TF ftau = -cd * pad_v *
                        std::pow( fm::pow2(u_on_v) +
                                  fm::pow2(v_on_v) +
                                  fm::pow2(w_on_v), TF(0.5) );

                    vt[ijk] += ftau * v_on_v;
                }
    }

    template<typename TF, bool sw_pad_3d>
    void canopy_drag_w(
            TF* const restrict wt,
            const TF* const restrict u,
            const TF* const restrict v,
            const TF* const restrict w,
            const TF* const restrict pad,
            const TF utrans,
            const TF vtrans,
            const TF cd,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride)
    {
        const int ii = 1;
        const int jj = jstride;
        const int kk = kstride;

        for (int k=kstart+1; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jj + k*kk;

                    // Interpolate `pad` to `w` locations.
                    const TF pad_w = sw_pad_3d ? TF(0.5) * (pad[ijk-kk] + pad[ijk]) : TF(0.5) * (pad[k-1] + pad[k]);

                    // Interpolate `u` and `v` to `w` locations.
                    const TF u_on_w = TF(0.25) * (u[ijk] + u[ijk+ii] + u[ijk+ii-kk] + u[ijk-kk]) + utrans;
                    const TF v_on_w = TF(0.25) * (v[ijk] + v[ijk+jj] + v[ijk+jj-kk] + v[ijk-kk]) + vtrans;
                    const TF w_on_w = w[ijk];

                    const TF ftau = -cd * pad_w *
                        std::pow( fm::pow2(u_on_w) +
                                  fm::pow2(v_on_w) +
                                  fm::pow2(w_on_w), TF(0.5) );

                    wt[ijk] += ftau * w_on_w;
                }
    }
}

template<typename TF>
Canopy<TF>::Canopy(
        Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
        master(masterin), grid(gridin), fields(fieldsin),
        field3d_operators(masterin, gridin, fieldsin), field3d_io(masterin, gridin),
        boundary_cyclic(masterin, gridin)
{
    sw_canopy = inputin.get_item<bool>("canopy", "sw_canopy", "", false);

    if (sw_canopy)
    {
        sw_3d_pad = inputin.get_item<bool>("canopy", "sw_3d_pad", "", false);
        cd = inputin.get_item<TF>("canopy", "cd", "");

        if (sw_3d_pad)
            ktot_canopy = inputin.get_item<int>("canopy", "ktot_canopy", "");
    }
}

template <typename TF>
Canopy<TF>::~Canopy()
{
}

template <typename TF>
void Canopy<TF>::init()
{
    if (!sw_canopy)
        return;

    auto& gd = grid.get_grid_data();

    if (sw_3d_pad)
    {
        if (ktot_canopy < 1 || ktot_canopy > gd.ktot)
            throw std::runtime_error("Invalid ktot_canopy, should be between 1 and ktot.");

        kend_canopy = gd.kstart + ktot_canopy;

        // One extra (zero) level on top, required for the interpolation
        // of `pad` to the `w` location at the top of the canopy.
        pad.resize((kend_canopy+1) * gd.ijcells, TF(0));

        boundary_cyclic.init();
    }
    else
        pad.resize(gd.kcells, TF(0));
}

template <typename TF>
void Canopy<TF>::create(
        Input& inputin, Netcdf_handle& input_nc, Stats<TF>& stats)
{
    if (!sw_canopy)
        return;

    auto& gd = grid.get_grid_data();

    if (sw_3d_pad)
    {
        // Read plant area density at full levels from binary
        auto tmp1 = fields.get_tmp();
        auto tmp2 = fields.get_tmp();

        char filename[256];
        std::sprintf(filename, "%s.%07d", "pad", 0);
        master.print_message("Loading \"%s\" ... ", filename);

        const TF no_offset = TF(0);
        if (field3d_io.load_field3d(
                pad.data(),
                tmp1->fld.data(),
                tmp2->fld.data(),
                filename,
                no_offset,
                gd.kstart,
                kend_canopy))
        {
            master.print_message("FAILED\n");
            throw std::runtime_error("Loading canopy area density failed!");
        }
        else
            master.print_message("OK\n");

        fields.release_tmp(tmp1);
        fields.release_tmp(tmp2);

        // Fill the horizontal ghost cells, required for the interpolation to the `u` and `v` locations.
        for (int k=gd.kstart; k<kend_canopy; ++k)
            boundary_cyclic.exec_2d(&pad[k*gd.ijcells]);
    }
    else
    {
        // Read plant area density at full levels.
        Netcdf_group& group_nc = input_nc.get_group("init");

        if (!group_nc.variable_exists("pad") && group_nc.variable_exists("padh"))
            throw std::runtime_error("Canopy: \"padh\" at half levels is no longer supported, provide \"pad\" at full levels.");

        group_nc.get_variable(pad, "pad", {0}, {gd.ktot});
        std::rotate(pad.rbegin(), pad.rbegin() + gd.kstart, pad.rend());

        // Determine end of canopy, as the first full level above the highest non-zero `pad`.
        kend_canopy = gd.kstart;
        for (int k=gd.kend-1; k>=gd.kstart; --k)
            if (pad[k] > TF(0))
            {
                kend_canopy = k+1;
                break;
            }
    }
}

#ifndef USECUDA
template <typename TF>
void Canopy<TF>::exec()
{
    if (!sw_canopy)
        return;

    auto& gd = grid.get_grid_data();

    auto canopy_drag_wrapper = [&]<const bool sw_pad_3d>()
    {
        canopy_drag_u<TF, sw_pad_3d>(
                fields.mt.at("u")->fld.data(),
                fields.mp.at("u")->fld.data(),
                fields.mp.at("v")->fld.data(),
                fields.mp.at("w")->fld.data(),
                pad.data(),
                gd.utrans,
                gd.vtrans,
                cd,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, kend_canopy,
                gd.icells, gd.ijcells);

        canopy_drag_v<TF, sw_pad_3d>(
                fields.mt.at("v")->fld.data(),
                fields.mp.at("u")->fld.data(),
                fields.mp.at("v")->fld.data(),
                fields.mp.at("w")->fld.data(),
                pad.data(),
                gd.utrans,
                gd.vtrans,
                cd,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, kend_canopy,
                gd.icells, gd.ijcells);

        // Include the half level at the top of the canopy, where `pad` interpolated to `w` is non-zero.
        const int kend_canopy_w = std::min(kend_canopy+1, gd.kend);

        canopy_drag_w<TF, sw_pad_3d>(
                fields.mt.at("w")->fld.data(),
                fields.mp.at("u")->fld.data(),
                fields.mp.at("v")->fld.data(),
                fields.mp.at("w")->fld.data(),
                pad.data(),
                gd.utrans,
                gd.vtrans,
                cd,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                gd.kstart, kend_canopy_w,
                gd.icells, gd.ijcells);
    };

    if (sw_3d_pad)
        canopy_drag_wrapper.template operator()<true>();
    else
        canopy_drag_wrapper.template operator()<false>();
}
#endif

#ifdef FLOAT_SINGLE
template class Canopy<float>;
#else
template class Canopy<double>;
#endif
