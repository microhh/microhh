/*
 * MicroHH
 * Copyright (c) 2011-2024 Chiel van Heerwaarden
 * Copyright (c) 2011-2024 Thijs Heus
 * Copyright (c) 2014-2024 Bart van Stratum
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

#include "tools.h"
#include "grid.h"
#include "fields.h"

#include "source_3d.h"
#include "source_3d_kernels.cuh"

namespace s3k = Source_3d_kernels_g;


#ifdef USECUDA
template<typename TF>
void Source_3d<TF>::exec()
{
    auto& gd = grid.get_grid_data();

    const int blocki = gd.ithread_block;
    const int blockj = gd.jthread_block;
    const int gridi  = gd.imax/blocki + (gd.imax%blocki > 0);
    const int gridj  = gd.jmax/blockj + (gd.jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, this->ktot);
    dim3 blockGPU(blocki, blockj, 1);

    for (auto& specie : sourcelist)
        s3k::add_source_tend_g<TF><<<gridGPU, blockGPU>>>(
            fields.st.at(specie)->fld_g,
            emission_g.at(specie),
            gd.istart, gd.iend,
            gd.jstart, gd.jend,
            gd.kstart, gd.kstart + this->ktot,
            gd.icells, gd.ijcells);
    cuda_check_error();
}


template<typename TF>
void Source_3d<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (!sw_timedep)
        return;

    throw std::runtime_error("Time dependent 3D emissions not (yet) implemented.");
}

    
template<typename TF>
void Source_3d<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();

    const int ncells = gd.ijcells * this->ktot;
    const int memsize = ncells * sizeof(TF);

    for (auto& specie : sourcelist)
    {
        emission_g[specie] = cuda_vector<TF>(ncells);
        cuda_safe_call(cudaMemcpy(emission_g.at(specie), emission.at(specie).data(), memsize, cudaMemcpyHostToDevice));
    }
}
#endif


#ifdef FLOAT_SINGLE
template class Source_3d<float>;
#else
template class Source_3d<double>;
#endif
