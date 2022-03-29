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

#include "radiation_prescribed.h"
#include "grid.h"
#include "tools.h"

#ifdef USECUDA
template<typename TF>
void Radiation_prescribed<TF>::exec(
        Thermo<TF>& thermo, const double time, Timeloop<TF>& timeloop, Stats<TF>& stats)
{
    if (swtimedep_prescribed)
    {
        auto& gd = grid.get_grid_data();
        const int memsize = gd.ijcells*sizeof(TF);

        cuda_safe_call(cudaMemset(sw_flux_dn_g, sw_flux_dn_value, memsize));
        cuda_safe_call(cudaMemset(sw_flux_up_g, sw_flux_up_value, memsize));
        cuda_safe_call(cudaMemset(lw_flux_dn_g, lw_flux_dn_value, memsize));
        cuda_safe_call(cudaMemset(lw_flux_up_g, lw_flux_up_value, memsize));
    }
}

template<typename TF>
TF* Radiation_prescribed<TF>::get_surface_radiation_g(std::string name)
{
    if (name == "sw_down")
        return sw_flux_dn_g;
    else if (name == "sw_up")
        return sw_flux_up_g;
    else if (name == "lw_down")
        return lw_flux_dn_g;
    else if (name == "lw_up")
        return lw_flux_up_g;
    else
    {
        std::string error = "Variable \"" + name + "\" is not a valid surface radiation field";
        throw std::runtime_error(error);
    }
}

template<typename TF>
void Radiation_prescribed<TF>::prepare_device()
{
    auto& gd = grid.get_grid_data();
    const int memsize = gd.ijcells*sizeof(TF);

    // Allocate surface radiation fields.
    cuda_safe_call(cudaMalloc(&lw_flux_dn_g, memsize));
    cuda_safe_call(cudaMalloc(&lw_flux_up_g, memsize));
    cuda_safe_call(cudaMalloc(&sw_flux_dn_g, memsize));
    cuda_safe_call(cudaMalloc(&sw_flux_up_g, memsize));

    // Send data to GPU, in case timedep is disabled.
    forward_device();
}

template<typename TF>
void Radiation_prescribed<TF>::clear_device()
{
    cuda_safe_call(cudaFree(lw_flux_dn_g));
    cuda_safe_call(cudaFree(lw_flux_up_g));
    cuda_safe_call(cudaFree(sw_flux_dn_g));
    cuda_safe_call(cudaFree(sw_flux_up_g));
}

template<typename TF>
void Radiation_prescribed<TF>::forward_device()
{
    auto& gd = grid.get_grid_data();
    const int memsize = gd.ijcells*sizeof(TF);

    cuda_safe_call(cudaMemcpy(lw_flux_dn_g, lw_flux_dn.data(), memsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(lw_flux_up_g, lw_flux_up.data(), memsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(sw_flux_dn_g, sw_flux_dn.data(), memsize, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(sw_flux_up_g, sw_flux_up.data(), memsize, cudaMemcpyHostToDevice));
}

template<typename TF>
void Radiation_prescribed<TF>::backward_device()
{
    auto& gd = grid.get_grid_data();
    const int memsize = gd.ijcells*sizeof(TF);

    cuda_safe_call(cudaMemcpy(lw_flux_dn.data(), lw_flux_dn_g, memsize, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(lw_flux_up.data(), lw_flux_up_g, memsize, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(sw_flux_dn.data(), sw_flux_dn_g, memsize, cudaMemcpyDeviceToHost));
    cuda_safe_call(cudaMemcpy(sw_flux_up.data(), sw_flux_up_g, memsize, cudaMemcpyDeviceToHost));
}
#endif

template class Radiation_prescribed<double>;
template class Radiation_prescribed<float>;
