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

#include <cstdio>
#include <iostream>
#include <cmath>
#include <stdlib.h>
#include "grid.h"
#include "fields.h"
#include "buffer.h"
#include "constants.h"
#include "stats.h"
#include "tools.h"

#include "buffer_kernels.cuh"
#include "cuda_launcher.h"
#include "cuda_tiling.h"

template<typename TF>
void Buffer<TF>::prepare_device()
{
    const Grid_data<TF>& gd = grid.get_grid_data();

    if (swbuffer)
    {
        const int nmemsize = gd.kcells*sizeof(TF);

        // Allocate the buffer arrays at GPU.
        for (auto& it : fields.ap)
        {
            bufferprofs_g.emplace(it.first, cuda_vector<TF>(gd.kcells));
            cuda_safe_call(cudaMemcpy(bufferprofs_g.at(it.first), bufferprofs.at(it.first).data(), nmemsize, cudaMemcpyHostToDevice));
        }

        // Pre-calculate buffer factor.
        sigma_z.allocate(gd.kcells);
        sigma_zh.allocate(gd.kcells);

        auto tmp = fields.get_tmp();
        const TF zsizebufi = 1./(gd.zsize-zstart);

        // Calculate & copy to device.
        for (int k=bufferkstart; k<gd.kend; k++)
            tmp->fld_mean[k] = sigma * pow((gd.z[k]-zstart)*zsizebufi, beta);
        cuda_safe_call(cudaMemcpy(sigma_z, tmp->fld_mean.data(), nmemsize, cudaMemcpyHostToDevice));

        for (int k=bufferkstarth; k<gd.kend; k++)
            tmp->fld_mean[k] = sigma * pow((gd.zh[k]-zstart)*zsizebufi, beta);
        cuda_safe_call(cudaMemcpy(sigma_zh, tmp->fld_mean.data(), nmemsize, cudaMemcpyHostToDevice));
    }
}

template<typename TF>
void Buffer<TF>::clear_device()
{
}

#ifdef USECUDA
template<typename TF>
void Buffer<TF>::exec(Stats<TF>& stats)
{
    if (swbuffer)
    {
        const Grid_data<TF>& gd = grid.get_grid_data();

        Grid_layout grid_layout_full = {
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                bufferkstart, gd.kend,
                gd.istride,
                gd.jstride,
                gd.kstride};

        Grid_layout grid_layout_half = {
                gd.istart, gd.iend,
                gd.jstart, gd.jend,
                bufferkstarth, gd.kend,
                gd.istride,
                gd.jstride,
                gd.kstride};

        if (swupdate)
        {
            launch_grid_kernel<Buffer_kernels::buffer_g<TF>>(
                    grid_layout_full,
                    fields.mt.at("u")->fld_g.view(),
                    fields.mp.at("u")->fld_g,
                    fields.mp.at("u")->fld_mean_g,
                    sigma_z);

            launch_grid_kernel<Buffer_kernels::buffer_g<TF>>(
                    grid_layout_full,
                    fields.mt.at("v")->fld_g.view(),
                    fields.mp.at("v")->fld_g,
                    fields.mp.at("v")->fld_mean_g,
                    sigma_z);

            launch_grid_kernel<Buffer_kernels::buffer_g<TF>>(
                    grid_layout_half,
                    fields.mt.at("w")->fld_g.view(),
                    fields.mp.at("w")->fld_g,
                    fields.mp.at("w")->fld_mean_g,
                    sigma_zh);

            for (auto& it : fields.sp)
                launch_grid_kernel<Buffer_kernels::buffer_g<TF>>(
                        grid_layout_full,
                        fields.st.at(it.first)->fld_g.view(),
                        fields.sp.at(it.first)->fld_g,
                        fields.sp.at(it.first)->fld_mean_g,
                        sigma_z);
        }
        else
        {
            launch_grid_kernel<Buffer_kernels::buffer_g<TF>>(
                    grid_layout_full,
                    fields.mt.at("u")->fld_g.view(),
                    fields.mp.at("u")->fld_g,
                    bufferprofs_g.at("u"),
                    sigma_z);

            launch_grid_kernel<Buffer_kernels::buffer_g<TF>>(
                    grid_layout_full,
                    fields.mt.at("v")->fld_g.view(),
                    fields.mp.at("v")->fld_g,
                    bufferprofs_g.at("v"),
                    sigma_z);

            launch_grid_kernel<Buffer_kernels::buffer_g<TF>>(
                    grid_layout_half,
                    fields.mt.at("w")->fld_g.view(),
                    fields.mp.at("w")->fld_g,
                    bufferprofs_g.at("w"),
                    sigma_zh);

            for (auto& it : fields.sp)
                launch_grid_kernel<Buffer_kernels::buffer_g<TF>>(
                        grid_layout_full,
                        fields.st.at(it.first)->fld_g.view(),
                        fields.sp.at(it.first)->fld_g,
                        bufferprofs_g.at(it.first),
                        sigma_z);
        }

        cudaDeviceSynchronize();
        stats.calc_tend(*fields.mt.at("u"), tend_name);
        stats.calc_tend(*fields.mt.at("v"), tend_name);
        stats.calc_tend(*fields.mt.at("w"), tend_name);
        for (auto it : fields.st)
            stats.calc_tend(*it.second, tend_name);
    }
}
#endif


#ifdef FLOAT_SINGLE
template class Buffer<float>;
#else
template class Buffer<double>;
#endif
