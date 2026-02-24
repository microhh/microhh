
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

#include <iostream>

#include "chemistry.h"
#include "tools.h"

#ifdef USECUDA
template <typename TF>
void Chemistry<TF>::exec(Thermo<TF>& thermo, const double sdt, const double dt)
{



}

template<typename TF>
void Chemistry<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();

    const int time_size = time.size();

    // Allocate GPU arrays.
    jval_g.allocate(n_jval);
    jo31d_g.allocate(time_size);
    jh2o2_g.allocate(time_size);
    jno2_g.allocate(time_size);
    jno3_g.allocate(time_size);
    jn2o5_g.allocate(time_size);
    jch2or_g.allocate(time_size);
    jch2om_g.allocate(time_size);
    jch3o2h_g.allocate(time_size);

    vdo3_g.allocate(gd.ijcells);
    vdno_g.allocate(gd.ijcells);
    vdno2_g.allocate(gd.ijcells);
    vdhno3_g.allocate(gd.ijcells);
    vdh2o2_g.allocate(gd.ijcells);
    vdrooh_g.allocate(gd.ijcells);
    vdhcho_g.allocate(gd.ijcells);

    // Copy data from host to device.
    cuda_safe_call(cudaMemcpy(jval_g,    jval.data(),    n_jval*sizeof(TF),     cudaMemcpyHostToDevice));

    const int memsize_time = time_size * sizeof(TF);
    cuda_safe_call(cudaMemcpy(jo31d_g,   jo31d.data(),   memsize_time, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(jh2o2_g,   jh2o2.data(),   memsize_time, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(jno2_g,    jno2.data(),    memsize_time, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(jno3_g,    jno3.data(),    memsize_time, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(jn2o5_g,   jn2o5.data(),   memsize_time, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(jch2or_g,  jch2or.data(),  memsize_time, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(jch2om_g,  jch2om.data(),  memsize_time, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(jch3o2h_g, jch3o2h.data(), memsize_time, cudaMemcpyHostToDevice));

    const int memsize_ij = gd.ijcells * sizeof(TF);
    cuda_safe_call(cudaMemcpy(vdo3_g,    vdo3.data(),    memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdno_g,    vdno.data(),    memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdno2_g,   vdno2.data(),   memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdhno3_g,  vdhno3.data(),  memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdh2o2_g,  vdh2o2.data(),  memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdrooh_g,  vdrooh.data(),  memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(vdhcho_g,  vdhcho.data(),  memsize_ij, cudaMemcpyHostToDevice));
}

template<typename TF>
void Chemistry<TF>::clear_device()
{
}
#endif

#ifdef FLOAT_SINGLE
template class Chemistry<float>;
#else
template class Chemistry<double>;
#endif