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

#include "radiation_prescribed.h"
#include "grid.h"
#include "tools.h"

namespace
{
    template<typename TF> __global__
    void set_surface_fluxes(
            TF* const __restrict__ sw_flux_dn_sfc,
            TF* const __restrict__ sw_flux_up_sfc,
            TF* const __restrict__ lw_flux_dn_sfc,
            TF* const __restrict__ lw_flux_up_sfc,
            const TF sw_flux_dn_in,
            const TF sw_flux_up_in,
            const TF lw_flux_dn_in,
            const TF lw_flux_up_in,
            const int icells, const int jcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij  = i + j*icells;

            sw_flux_dn_sfc[ij] = sw_flux_dn_in;
            sw_flux_up_sfc[ij] = sw_flux_up_in;
            lw_flux_dn_sfc[ij] = lw_flux_dn_in;
            lw_flux_up_sfc[ij] = lw_flux_up_in;
        }
    }
}

#ifdef USECUDA
template<typename TF>
void Radiation_prescribed<TF>::exec(
        Thermo<TF>& thermo, const double time, Timeloop<TF>& timeloop, Stats<TF>& stats)
{
    if (swtimedep_prescribed)
    {
        auto& gd = grid.get_grid_data();

        const int blocki = gd.ithread_block;
        const int blockj = gd.jthread_block;
        const int gridi = gd.icells/blocki + (gd.icells%blocki > 0);
        const int gridj = gd.jcells/blockj + (gd.jcells%blockj > 0);
        dim3 gridGPU (gridi,  gridj,  1);
        dim3 blockGPU(blocki, blockj, 1);

        set_surface_fluxes<<<gridGPU, blockGPU>>>(
            sw_flux_dn_sfc_g,
            sw_flux_up_sfc_g,
            lw_flux_dn_sfc_g,
            lw_flux_up_sfc_g,
            sw_flux_dn,
            sw_flux_up,
            lw_flux_dn,
            lw_flux_up,
            gd.icells, gd.jcells);
        cuda_check_error();
    }
}
#endif

template class Radiation_prescribed<double>;
template class Radiation_prescribed<float>;
