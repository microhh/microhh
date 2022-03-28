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

#ifndef LAND_SURFACE_KERNELS_GPU_H
#define LAND_SURFACE_KERNELS_GPU_H

#include "tools.h"

namespace Land_surface_kernels_g
{
    template<typename TF>
    void init_tile(Surface_tile<TF>& tile, const int ijcells)
    {
        const int memsize_tf  = ijcells * sizeof(TF);
        const int memsize_int = ijcells * sizeof(int);

        cuda_safe_call(cudaMalloc(&tile.fraction_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.thl_bot_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.qt_bot_g, memsize_tf));

        cuda_safe_call(cudaMalloc(&tile.obuk_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.ustar_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.bfluxbot_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.ra_g, memsize_tf));

        cuda_safe_call(cudaMalloc(&tile.nobuk_g, memsize_int));

        cuda_safe_call(cudaMalloc(&tile.rs_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.H_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.LE_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.G_g, memsize_tf));
        cuda_safe_call(cudaMalloc(&tile.S_g, memsize_tf));
    }

    template<typename TF>
    void forward_device_tile(Surface_tile<TF>& tile, const int ijcells)
    {
        const int memsize_tf  = ijcells * sizeof(TF);
        const int memsize_int = ijcells * sizeof(int);

        cuda_safe_call(cudaMemcpy(tile.fraction_g, tile.fraction.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.thl_bot_g, tile.thl_bot.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.qt_bot_g, tile.qt_bot.data(), memsize_tf, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(tile.obuk_g, tile.obuk.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.ustar_g, tile.ustar.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.bfluxbot_g, tile.bfluxbot.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.ra_g, tile.ra.data(), memsize_tf, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(tile.nobuk_g, tile.nobuk.data(), memsize_int, cudaMemcpyHostToDevice));

        cuda_safe_call(cudaMemcpy(tile.rs_g, tile.rs.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.H_g, tile.H.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.LE_g, tile.LE.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.G_g, tile.G.data(), memsize_tf, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(tile.S_g, tile.S.data(), memsize_tf, cudaMemcpyHostToDevice));
    }
}
#endif
