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

#include "boundary_surface_lsm.h"
#include "boundary.h"
#include "tools.h"


#ifdef USECUDA
template<typename TF>
void Boundary_surface_lsm<TF>::exec(
        Thermo<TF>& thermo, Radiation<TF>& radiation,
        Microphys<TF>& microphys, Timeloop<TF>& timeloop)
{
}

template<typename TF>
void Boundary_surface_lsm<TF>::exec_column(Column<TF>& column)
{
}

template<typename TF>
void Boundary_surface_lsm<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();

    const int tf_memsize_ij  = gd.ijcells*sizeof(TF);
    const int int_memsize_ij = gd.ijcells*sizeof(int);
    const int float_memsize_lut = nzL_lut*sizeof(float);

    // Monin-Obukhov stuff:
    cuda_safe_call(cudaMalloc(&obuk_g,  tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&ustar_g, tf_memsize_ij));

    cuda_safe_call(cudaMalloc(&z0m_g,   tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&z0h_g,   tf_memsize_ij));

    cuda_safe_call(cudaMalloc(&dudz_mo_g, tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&dvdz_mo_g, tf_memsize_ij));
    cuda_safe_call(cudaMalloc(&dbdz_mo_g, tf_memsize_ij));

    cuda_safe_call(cudaMalloc(&nobuk_g, int_memsize_ij));
    cuda_safe_call(cudaMalloc(&zL_sl_g, float_memsize_lut));
    cuda_safe_call(cudaMalloc(&f_sl_g,  float_memsize_lut));

    // Land-surface stuff:

    // Copy data from host to device
    forward_device();
}

template<typename TF>
void Boundary_surface_lsm<TF>::forward_device()
{
    auto& gd = grid.get_grid_data();

    const int tf_memsize_ij  = gd.ijcells*sizeof(TF);
    const int int_memsize_ij = gd.ijcells*sizeof(int);
    const int float_memsize_lut = nzL_lut*sizeof(float);

    cuda_safe_call(cudaMemcpy(obuk_g,  obuk.data(),  tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(ustar_g, ustar.data(), tf_memsize_ij, cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy(z0m_g, z0m.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(z0h_g, z0h.data(), tf_memsize_ij, cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy(dudz_mo_g, dudz_mo.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(dvdz_mo_g, dvdz_mo.data(), tf_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(dbdz_mo_g, dbdz_mo.data(), tf_memsize_ij, cudaMemcpyHostToDevice));

    cuda_safe_call(cudaMemcpy(nobuk_g, nobuk.data(), int_memsize_ij, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(zL_sl_g, zL_sl.data(), float_memsize_lut, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(f_sl_g,  f_sl.data(),  float_memsize_lut, cudaMemcpyHostToDevice));
}

template<typename TF>
void Boundary_surface_lsm<TF>::backward_device()
{
    auto& gd = grid.get_grid_data();

    const int tf_memsize_ij  = gd.ijcells*sizeof(TF);
    const int int_memsize_ij = gd.ijcells*sizeof(int);

    // NOTE: only copy back the required/useful data...

    cuda_safe_call(cudaMemcpy(obuk.data(),  obuk_g,  tf_memsize_ij, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(ustar.data(), ustar_g, tf_memsize_ij, cudaMemcpyDeviceToHost));

    cuda_safe_call(cudaMemcpy(dudz_mo.data(), dudz_mo_g, tf_memsize_ij, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(dvdz_mo.data(), dvdz_mo_g, tf_memsize_ij, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(dbdz_mo.data(), dbdz_mo_g, tf_memsize_ij, cudaMemcpyDeviceToHost));

    cuda_safe_call(cudaMemcpy(nobuk.data(), nobuk_g, int_memsize_ij, cudaMemcpyDeviceToHost));
}

template<typename TF>
void Boundary_surface_lsm<TF>::clear_device()
{
    //
    // De-llocate fields on GPU
    //
    // Monin-Obukhov stuff:
    cuda_safe_call(cudaFree(obuk_g));
    cuda_safe_call(cudaFree(ustar_g));

    cuda_safe_call(cudaFree(z0m_g));
    cuda_safe_call(cudaFree(z0h_g));

    cuda_safe_call(cudaFree(dudz_mo_g));
    cuda_safe_call(cudaFree(dvdz_mo_g));
    cuda_safe_call(cudaFree(dbdz_mo_g));

    cuda_safe_call(cudaFree(nobuk_g));
    cuda_safe_call(cudaFree(zL_sl_g));
    cuda_safe_call(cudaFree(f_sl_g));

    // Land-surface stuff:
}
#endif

template class Boundary_surface_lsm<double>;
template class Boundary_surface_lsm<float>;
