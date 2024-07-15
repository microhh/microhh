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
#include <iostream>

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
    void calc_buffer_3d(
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
                    const int ijk2 = i + j * jstride + (k-kstart_buffer) * kstride;

                    const TF aa = a[ijk];
                    const TF bb = abuf[ijk2];

                    at[ijk] -= sigmaz * (a[ijk]-abuf[ijk2]);
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


    template<typename TF>
    void interpolate_buffer(
            TF* const restrict buf_out,
            const TF* const restrict buf_prev,
            const TF* const restrict buf_next,
            const TF f0, const TF f1,
            const int ncells)
    {
        for (int n=0; n<ncells; ++n)
            buf_out[n] = f0 * buf_prev[n] + f1 * buf_next[n];
    }
}


template<typename TF>
Buffer<TF>::Buffer(Master& masterin, Grid<TF>& gridin, Fields<TF>& fieldsin, Input& inputin) :
    master(masterin), grid(gridin), fields(fieldsin), field3d_io(masterin, gridin)
{
    swbuffer = inputin.get_item<bool>("buffer", "swbuffer", "", false);

    if (swbuffer)
    {
        swupdate = inputin.get_item<bool>("buffer", "swupdate", "", false);
        swupdate_local = inputin.get_item<bool>("buffer", "swupdate_local", "", false);
        swbuffer_3d = inputin.get_item<bool>("buffer", "swbuffer_3d", "", false);

        zstart = inputin.get_item<TF>("buffer", "zstart", "");
        sigma  = inputin.get_item<TF>("buffer", "sigma", "", 2.);
        beta   = inputin.get_item<TF>("buffer", "beta", "", 2.);

        if (swbuffer_3d)
        {
            buffer3d_list = inputin.get_list<std::string>("buffer", "buffer3d_list", "", std::vector<std::string>());
            swtimedep_buffer_3d = inputin.get_item<bool>("buffer", "swtimedep_buffer_3d", "", false);

            if (swtimedep_buffer_3d)
                loadfreq = inputin.get_item<int>("buffer", "loadfreq", "");
        }
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
    if (swbuffer && !swbuffer_3d)
    {
        const Grid_data<TF>& gd = grid.get_grid_data();

        // Create vectors of zero for buffer.
        for (auto& it : fields.ap)
            bufferprofs.emplace(it.first, std::vector<TF>(gd.kcells));
    }
}

template<typename TF>
void Buffer<TF>::create(
        Input& inputin,
        Netcdf_handle& input_nc,
        Stats<TF>& stats,
        Timeloop<TF>& timeloop)
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

        if (swbuffer_3d)
        {
            auto tmp1 = fields.get_tmp();
            auto tmp2 = fields.get_tmp();

            const TF no_offset = TF(0);
            int nerror = 0;

            // Read 3D buffers from binary files.
            auto load_3d_field = [&](
                    TF* const restrict field,
                    const std::string& name,
                    const int itime,
                    const int kstart,
                    const int kend)
            {
                char filename[256];
                std::sprintf(filename, "%s_buffer.%07d", name.c_str(), itime);
                master.print_message("Loading \"%s\" ... ", filename);

                if (field3d_io.load_field3d(
                        field,
                        tmp1->fld.data(), tmp2->fld.data(),
                        filename, no_offset,
                        kstart, kend))
                {
                    master.print_message("FAILED\n");
                    nerror += 1;
                }
                else
                    master.print_message("OK\n");
            };

            if (swtimedep_buffer_3d)
            {
                std::pair<unsigned long, unsigned long> iotimes = timeloop.get_prev_and_next_iotime(loadfreq);

                for (auto& fld : buffer3d_list)
                {
                    const int kstart = 0;
                    const int ksize = fld == "w" ? gd.kend-bufferkstarth : gd.kend-bufferkstart;

                    // Allocate arrays in buffer map.
                    bufferprofs.emplace(fld, std::vector<TF>(gd.ijcells*ksize));
                    buffer_data_prev.emplace(fld, std::vector<TF>(gd.ijcells * ksize));
                    buffer_data_next.emplace(fld, std::vector<TF>(gd.ijcells * ksize));

                    // Read data from binary.
                    load_3d_field(
                            buffer_data_prev.at(fld).data(),
                            fld, iotimes.first, kstart, ksize);

                    load_3d_field(
                            buffer_data_next.at(fld).data(),
                            fld, iotimes.second, kstart, ksize);
                }
            }
            else
            {
                // Read fixed 3D buffers from binary files.
                for (auto& fld : buffer3d_list)
                {
                    const int kstart = 0;
                    const int ksize = fld == "w" ? gd.kend-bufferkstarth : gd.kend-bufferkstart;
                    const int itime = 0;

                    // Allocate array in std::map.
                    bufferprofs.emplace(fld, std::vector<TF>(gd.ijcells*ksize));

                    // Read data from binary.
                    load_3d_field(
                            bufferprofs.at(fld).data(),
                            fld, itime, kstart, ksize);
                }
            }

            if (nerror > 0)
                throw std::runtime_error("Error loading buffer fields.");

            fields.release_tmp(tmp1);
            fields.release_tmp(tmp2);
        }

        if (!swupdate && !swbuffer_3d)
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


    auto buffer_3d_wrapper = [&](
            TF* const restrict tend,
            const TF* const restrict field,
            const TF* const restrict buffer_data,
            const int kstart)
    {
        calc_buffer_3d(
                tend,
                field,
                buffer_data,
                gd.z.data(), zstart, gd.zsize,
                beta, sigma,
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                kstart, gd.kend,
                gd.icells, gd.ijcells);
    };


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

    if (swbuffer_3d)
    {
        for (auto& fld: buffer3d_list)
        {
            const int kstart = fld == "w" ? bufferkstarth : bufferkstart;

            buffer_3d_wrapper(
                    fields.at.at(fld)->fld.data(),
                    fields.ap.at(fld)->fld.data(),
                    bufferprofs.at(fld).data(),
                    kstart);
        }
    }
    else
    {
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

                for (auto& it: fields.sp)
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

                for (auto& it: fields.sp)
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

            for (auto& it: fields.sp)
                buffer_wrapper(
                        fields.st.at(it.first)->fld.data(),
                        fields.sp.at(it.first)->fld.data(),
                        bufferprofs.at(it.first).data(),
                        bufferkstart);
        }
    }

    stats.calc_tend(*fields.mt.at("u"), tend_name);
    stats.calc_tend(*fields.mt.at("v"), tend_name);
    stats.calc_tend(*fields.mt.at("w"), tend_name);

    for (auto it : fields.st)
        stats.calc_tend(*it.second, tend_name);

    fields.release_tmp(tmp);
}


template <typename TF>
void Buffer<TF>::update_time_dependent(
        Timeloop<TF>& timeloop)
{
    if (!swtimedep_buffer_3d)
        return;

    const Grid_data<TF>& gd = grid.get_grid_data();

    double time = timeloop.get_time();
    unsigned long itime = timeloop.get_itime();
    unsigned long iloadtime = convert_to_itime(loadfreq);
    unsigned long iiotimeprec = timeloop.get_iiotimeprec();

    if (itime > next_itime)
    {
        // Advance time and read new files.
        prev_itime = next_itime;
        next_itime = prev_itime + iloadtime;

        unsigned long next_iotime = next_itime / iiotimeprec;

        auto tmp1 = fields.get_tmp();
        auto tmp2 = fields.get_tmp();

        const TF no_offset = TF(0);
        int nerror = 0;

        // Read 3D buffers from binary files.
        auto load_3d_field = [&](
                TF* const restrict field,
                const std::string& name,
                const int itime,
                const int kstart,
                const int kend)
        {
            char filename[256];
            std::sprintf(filename, "%s_buffer.%07d", name.c_str(), itime);
            master.print_message("Loading \"%s\" ... ", filename);

            if (field3d_io.load_field3d(
                    field,
                    tmp1->fld.data(), tmp2->fld.data(),
                    filename, no_offset,
                    kstart, kend))
            {
                master.print_message("FAILED\n");
                nerror += 1;
            }
            else
                master.print_message("OK\n");
        };

        for (auto& fld : buffer3d_list)
        {
            // Copy old next to new prev buffer.
            buffer_data_prev.at(fld) = buffer_data_next.at(fld);

            // Read new 3D data.
            const int kstart = 0;
            const int ksize = fld == "w" ? gd.kend-bufferkstarth : gd.kend-bufferkstart;

            load_3d_field(
                    buffer_data_next.at(fld).data(),
                    fld, next_iotime, kstart, ksize);
        }
    }

    // Interpolate in time.
    const TF f0 = TF(1) - ((itime - prev_itime) / TF(iloadtime));
    const TF f1 = TF(1) - f0;

    for (auto& fld : buffer3d_list)
    {
        const int ksize = fld == "w" ? gd.kend-bufferkstarth : gd.kend-bufferkstart;
        const int ncells = ksize * gd.ijcells;

        interpolate_buffer(
            bufferprofs.at(fld).data(),
            buffer_data_prev.at(fld).data(),
            buffer_data_next.at(fld).data(),
            f0, f1, ncells);
    }
}
#endif


#ifdef FLOAT_SINGLE
template class Buffer<float>;
#else
template class Buffer<double>;
#endif
