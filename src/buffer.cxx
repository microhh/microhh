/*
 * MicroHH
 * Copyright (c) 2011-2023 Chiel van Heerwaarden
 * Copyright (c) 2011-2023 Thijs Heus
 * Copyright (c) 2014-2023 Bart van Stratum
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

#include <cmath>
#include <algorithm>

#include "master.h"
#include "input.h"
#include "grid.h"
#include "fields.h"
#include "stats.h"
#include "buffer.h"
#include "defines.h"
#include "netcdf_interface.h"

namespace
{
    template<typename TF>
    void calc_buffer(
            TF* const restrict at,
            const TF* const restrict a,
            const TF* const restrict abuf,
            const TF* const restrict z,
            const TF zstart, const TF zsize,
            const TF beta, const TF sigma,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart_buffer, const int kend,
            const int jstride, const int kstride)
    {
        const TF zsizebuf = zsize - zstart;

        for (int k=kstart_buffer; k<kend; ++k)
        {
            const TF sigmaz = sigma*std::pow((z[k]-zstart)/zsizebuf, beta);

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j * jstride + k * kstride;
                    at[ijk] -= sigmaz * (a[ijk]-abuf[k]);
                }
        }
    }


    template<typename TF>
    void calc_buffer_local(
            TF* const restrict at,
            const TF* const restrict a,
            TF* const restrict abuf,
            const TF* const restrict z,
            const TF zstart, const TF zsize,
            const TF beta, const TF sigma,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart_buffer, const int kend,
            const int jstride, const int kstride)
    {
        const TF zsizebuf = zsize - zstart;

        for (int k=kstart_buffer; k<kend; ++k)
        {
            const TF sigmaz = sigma*std::pow((z[k]-zstart)/zsizebuf, beta);

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jstride;
                    const int ijk = i + j*jstride + k*kstride;

                    abuf[ij] = TF(0.);

                    for (int jc=-3; jc<4; ++jc)
                        for (int ic=-3; ic<4; ++ic)
                            abuf[ij] += a[ijk + ic + jc*jstride];

                    abuf[ij] /= TF(49.);
                }

            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij = i + j*jstride;
                    const int ijk = i + j*jstride + k*kstride;

                    at[ijk] -= sigmaz*(a[ijk]-abuf[ij]);
                }
        }
    }
}


template<typename TF>
Buffer<TF>::Buffer(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin)
{
    swbuffer = inputin.get_item<bool>("buffer", "swbuffer", "", false);

    if (swbuffer)
    {
        swupdate = inputin.get_item<bool>("buffer", "swupdate", "", false);
        swupdate_local = inputin.get_item<bool>("buffer", "swupdate_local", "", false);

        zstart = inputin.get_item<TF>("buffer", "zstart", "");
        sigma  = inputin.get_item<TF>("buffer", "sigma", "", 2.);
        beta   = inputin.get_item<TF>("buffer", "beta", "", 2.);
    }

    if (swbuffer && swupdate)
        fields.set_calc_mean_profs(true);
}

template<typename TF>
Buffer<TF>::~Buffer()
{
    #ifdef USECUDA
    clear_device();
    #endif
}

template<typename TF>
void Buffer<TF>::init()
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    if (swbuffer)
    {
        // Create vectors of zero for buffer.
        for (auto& it : fields.ap)
            bufferprofs.emplace(it.first, std::vector<TF>(gd.kcells));
    }
}

template<typename TF>
void Buffer<TF>::create(Input& inputin, Netcdf_handle& input_nc, Stats<TF>& stats)
{
    if (swbuffer)
    {
        const Grid_data<TF>& gd = grid.get_grid_data();

        // Find the starting points.
        bufferkstart  = gd.kstart;
        bufferkstarth = gd.kstart;

        for (int k=gd.kstart; k<gd.kend; ++k)
        {
            // Check if the cell center is in the buffer zone.
            if (gd.z[k] < zstart)
                ++bufferkstart;
            // Check if the cell face is in the buffer zone.
            if (gd.zh[k] < zstart)
                ++bufferkstarth;
        }

        // Check whether the lowest of the two levels is contained in the buffer layer.
        if (bufferkstarth == gd.kend)
        {
            std::string msg = "Buffer is too close to the model top";
            throw std::runtime_error(msg);
        }

        if (!swupdate)
        {
            // Set the buffers according to the initial profiles of the variables.
            const std::vector<int> start = {0};
            const std::vector<int> count = {gd.ktot};

            Netcdf_group& group_nc = input_nc.get_group("init");
            group_nc.get_variable(bufferprofs.at("u"), "u", start, count);
            group_nc.get_variable(bufferprofs.at("v"), "v", start, count);
            std::rotate(bufferprofs.at("u").rbegin(), bufferprofs.at("u").rbegin() + gd.kstart, bufferprofs.at("u").rend());
            std::rotate(bufferprofs.at("v").rbegin(), bufferprofs.at("v").rbegin() + gd.kstart, bufferprofs.at("v").rend());

            // In case of u and v, subtract the grid velocity.
            for (int k=gd.kstart; k<gd.kend; ++k)
            {
                bufferprofs.at("u")[k] -= gd.utrans;
                bufferprofs.at("v")[k] -= gd.vtrans;
            }

            for (auto& it : fields.sp)
            {
                group_nc.get_variable(bufferprofs.at(it.first), it.first, start, count, fields.required_read.at(it.first));
                std::rotate(bufferprofs.at(it.first).rbegin(), bufferprofs.at(it.first).rbegin() + gd.kstart, bufferprofs.at(it.first).rend());
            }
        }

        stats.add_tendency(*fields.mt.at("u"), "z", tend_name, tend_longname);
        stats.add_tendency(*fields.mt.at("v"), "z", tend_name, tend_longname);
        stats.add_tendency(*fields.mt.at("w"), "zh", tend_name, tend_longname);

        for (auto it : fields.st)
            stats.add_tendency(*it.second, "z", tend_name, tend_longname);
    }
}

#ifndef USECUDA
template<typename TF>
void Buffer<TF>::exec(Stats<TF>& stats)
{
    if (!swbuffer)
        return;

    const Grid_data<TF>& gd = grid.get_grid_data();
    auto tmp = fields.get_tmp();

    auto buffer_local_wrapper = [&](
            TF* const restrict tend,
            const TF* const restrict field,
            const int kstart)
    {
        calc_buffer_local(
                tend,
                field,
                tmp->fld.data(),
                gd.z.data(),
                zstart, gd.zsize,
                beta, sigma,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                kstart, gd.kend,
                gd.icells, gd.ijcells);
    };

    auto buffer_wrapper = [&](
            TF* const restrict tend,
            const TF* const restrict field,
            const TF* const restrict buffer_prof,
            const int kstart)
    {
        calc_buffer(
                tend,
                field,
                buffer_prof,
                gd.z.data(), zstart, gd.zsize,
                beta, sigma,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                kstart, gd.kend,
                gd.icells, gd.ijcells);
    };

    if (swupdate)
    {
        if (swupdate_local)
        {
            // Buffer local to a 49 grid point mean.
            buffer_local_wrapper(
                    fields.mt.at("u")->fld.data(),
                    fields.mp.at("u")->fld.data(),
                    bufferkstart);

            buffer_local_wrapper(
                    fields.mt.at("v")->fld.data(),
                    fields.mp.at("v")->fld.data(),
                    bufferkstart);

            buffer_local_wrapper(
                    fields.mt.at("w")->fld.data(),
                    fields.mp.at("w")->fld.data(),
                    bufferkstarth);

            for (auto& it : fields.sp)
                buffer_local_wrapper(
                        fields.st.at(it.first)->fld.data(),
                        fields.sp.at(it.first)->fld.data(),
                        bufferkstart);
        }
        else
        {
            // Buffer to domain mean time updated profiles.
            buffer_wrapper(
                    fields.mt.at("u")->fld.data(),
                    fields.mp.at("u")->fld.data(),
                    fields.mp.at("u")->fld_mean.data(),
                    bufferkstart);

            buffer_wrapper(
                    fields.mt.at("v")->fld.data(),
                    fields.mp.at("v")->fld.data(),
                    fields.mp.at("v")->fld_mean.data(),
                    bufferkstart);

            buffer_wrapper(
                    fields.mt.at("w")->fld.data(),
                    fields.mp.at("w")->fld.data(),
                    fields.mp.at("w")->fld_mean.data(),
                    bufferkstarth);

            for (auto& it : fields.sp)
                buffer_wrapper(
                        fields.st.at(it.first)->fld.data(),
                        fields.sp.at(it.first)->fld.data(),
                        fields.sp.at(it.first)->fld_mean.data(),
                        bufferkstart);
        }
    }
    else
    {
        // Buffer to initial profiles.
        buffer_wrapper(
                fields.mt.at("u")->fld.data(),
                fields.mp.at("u")->fld.data(),
                bufferprofs.at("u").data(),
                bufferkstart);

        buffer_wrapper(
                fields.mt.at("v")->fld.data(),
                fields.mp.at("v")->fld.data(),
                bufferprofs.at("v").data(),
                bufferkstart);

        buffer_wrapper(
                fields.mt.at("w")->fld.data(),
                fields.mp.at("w")->fld.data(),
                bufferprofs.at("w").data(),
                bufferkstarth);

        for (auto& it : fields.sp)
            buffer_wrapper(
                    fields.st.at(it.first)->fld.data(),
                    fields.sp.at(it.first)->fld.data(),
                    bufferprofs.at(it.first).data(),
                    bufferkstart);
    }

    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);

    for (auto it : fields.st)
        stats.calc_tend(*it.second, tend_name);

    fields.release_tmp(tmp);
}
#endif


#ifdef FLOAT_SINGLE
template class Buffer<float>;
#else
template class Buffer<double>;
#endif
