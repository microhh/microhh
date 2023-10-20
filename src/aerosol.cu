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

#include <algorithm>
#include "grid.h"
#include "tools.h"
#include "timedep.h"
#include "aerosol.h"
#include "Array.h"

template<typename TF>
void Aerosol<TF>::prepare_device()
{
    if (sw_aerosol && sw_timedep)
    {
        auto &gd = grid.get_grid_data();
        const int nmemsize = gd.kcells * sizeof(TF);

        cuda_safe_call(cudaMalloc(&aermr01_g, nmemsize));
        cuda_safe_call(cudaMalloc(&aermr02_g, nmemsize));
        cuda_safe_call(cudaMalloc(&aermr03_g, nmemsize));
        cuda_safe_call(cudaMalloc(&aermr04_g, nmemsize));
        cuda_safe_call(cudaMalloc(&aermr05_g, nmemsize));
        cuda_safe_call(cudaMalloc(&aermr06_g, nmemsize));
        cuda_safe_call(cudaMalloc(&aermr07_g, nmemsize));
        cuda_safe_call(cudaMalloc(&aermr08_g, nmemsize));
        cuda_safe_call(cudaMalloc(&aermr09_g, nmemsize));
        cuda_safe_call(cudaMalloc(&aermr10_g, nmemsize));
        cuda_safe_call(cudaMalloc(&aermr11_g, nmemsize));

        cuda_safe_call(cudaMemcpy(aermr01_g, aermr01.data(), nmemsize, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(aermr02_g, aermr02.data(), nmemsize, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(aermr03_g, aermr03.data(), nmemsize, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(aermr04_g, aermr04.data(), nmemsize, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(aermr05_g, aermr05.data(), nmemsize, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(aermr06_g, aermr06.data(), nmemsize, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(aermr07_g, aermr07.data(), nmemsize, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(aermr08_g, aermr08.data(), nmemsize, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(aermr09_g, aermr09.data(), nmemsize, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(aermr10_g, aermr10.data(), nmemsize, cudaMemcpyHostToDevice));
        cuda_safe_call(cudaMemcpy(aermr11_g, aermr11.data(), nmemsize, cudaMemcpyHostToDevice));
    }
}

template<typename TF>
void Aerosol<TF>::clear_device()
{
    if (sw_aerosol && sw_timedep)
    {
        cuda_safe_call(cudaFree(aermr01_g));
        cuda_safe_call(cudaFree(aermr02_g));
        cuda_safe_call(cudaFree(aermr03_g));
        cuda_safe_call(cudaFree(aermr04_g));
        cuda_safe_call(cudaFree(aermr05_g));
        cuda_safe_call(cudaFree(aermr06_g));
        cuda_safe_call(cudaFree(aermr07_g));
        cuda_safe_call(cudaFree(aermr08_g));
        cuda_safe_call(cudaFree(aermr09_g));
        cuda_safe_call(cudaFree(aermr10_g));
        cuda_safe_call(cudaFree(aermr11_g));

        tdep_aermr01->clear_device();
        tdep_aermr02->clear_device();
        tdep_aermr03->clear_device();
        tdep_aermr04->clear_device();
        tdep_aermr05->clear_device();
        tdep_aermr06->clear_device();
        tdep_aermr07->clear_device();
        tdep_aermr08->clear_device();
        tdep_aermr09->clear_device();
        tdep_aermr10->clear_device();
        tdep_aermr11->clear_device();
    }
}

#ifdef USECUDA
template <typename TF>
void Aerosol<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (!sw_aerosol)
        return;

    if (sw_timedep)
    {
        tdep_aermr01 ->update_time_dependent_prof_g(aermr01_g, timeloop);
        tdep_aermr02 ->update_time_dependent_prof_g(aermr02_g, timeloop);
        tdep_aermr03 ->update_time_dependent_prof_g(aermr03_g, timeloop);
        tdep_aermr04 ->update_time_dependent_prof_g(aermr04_g, timeloop);
        tdep_aermr05 ->update_time_dependent_prof_g(aermr05_g, timeloop);
        tdep_aermr06 ->update_time_dependent_prof_g(aermr06_g, timeloop);
        tdep_aermr07 ->update_time_dependent_prof_g(aermr07_g, timeloop);
        tdep_aermr08 ->update_time_dependent_prof_g(aermr08_g, timeloop);
        tdep_aermr09 ->update_time_dependent_prof_g(aermr09_g, timeloop);
        tdep_aermr10 ->update_time_dependent_prof_g(aermr10_g, timeloop);
        tdep_aermr11 ->update_time_dependent_prof_g(aermr11_g, timeloop);
    }
}

template<typename TF>
void Aerosol<TF>::get_radiation_fields(std::unique_ptr<Aerosol_concs_gpu>& aerosol_concs_gpu)
{
    auto& gd = grid.get_grid_data();
    const int ncol = 1;

    Array_gpu<Float,2> aermr01_a(aermr01_g, {1, int(gd.kcells)});
    Array_gpu<Float,2> aermr02_a(aermr02_g, {1, int(gd.kcells)});
    Array_gpu<Float,2> aermr03_a(aermr03_g, {1, int(gd.kcells)});
    Array_gpu<Float,2> aermr04_a(aermr04_g, {1, int(gd.kcells)});
    Array_gpu<Float,2> aermr05_a(aermr05_g, {1, int(gd.kcells)});
    Array_gpu<Float,2> aermr06_a(aermr06_g, {1, int(gd.kcells)});
    Array_gpu<Float,2> aermr07_a(aermr07_g, {1, int(gd.kcells)});
    Array_gpu<Float,2> aermr08_a(aermr08_g, {1, int(gd.kcells)});
    Array_gpu<Float,2> aermr09_a(aermr09_g, {1, int(gd.kcells)});
    Array_gpu<Float,2> aermr10_a(aermr10_g, {1, int(gd.kcells)});
    Array_gpu<Float,2> aermr11_a(aermr11_g, {1, int(gd.kcells)});

    const int kstart = gd.kgc+1;
    const int kend = gd.ktot + gd.kgc;

    aerosol_concs_gpu->set_vmr("aermr01", aermr01_a.subset({ {{1, ncol}, {kstart, kend}}} ));
    aerosol_concs_gpu->set_vmr("aermr02", aermr02_a.subset({ {{1, ncol}, {kstart, kend}}} ));
    aerosol_concs_gpu->set_vmr("aermr03", aermr03_a.subset({ {{1, ncol}, {kstart, kend}}} ));
    aerosol_concs_gpu->set_vmr("aermr04", aermr04_a.subset({ {{1, ncol}, {kstart, kend}}} ));
    aerosol_concs_gpu->set_vmr("aermr05", aermr05_a.subset({ {{1, ncol}, {kstart, kend}}} ));
    aerosol_concs_gpu->set_vmr("aermr06", aermr06_a.subset({ {{1, ncol}, {kstart, kend}}} ));
    aerosol_concs_gpu->set_vmr("aermr07", aermr07_a.subset({ {{1, ncol}, {kstart, kend}}} ));
    aerosol_concs_gpu->set_vmr("aermr08", aermr08_a.subset({ {{1, ncol}, {kstart, kend}}} ));
    aerosol_concs_gpu->set_vmr("aermr09", aermr09_a.subset({ {{1, ncol}, {kstart, kend}}} ));
    aerosol_concs_gpu->set_vmr("aermr10", aermr10_a.subset({ {{1, ncol}, {kstart, kend}}} ));
    aerosol_concs_gpu->set_vmr("aermr11", aermr11_a.subset({ {{1, ncol}, {kstart, kend}}} ));

}
#endif


#ifdef FLOAT_SINGLE
template class Aerosol<float>;
#else
template class Aerosol<double>;
#endif
