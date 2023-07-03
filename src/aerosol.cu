#include <algorithm>
#include "grid.h"
#include "tools.h"
#include "timedep.h"
#include "aerosol.h"

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
        auto& gd = grid.get_grid_data();
        const int nmemsize  = gd.kcells*sizeof(TF);

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

        cuda_safe_call(cudaMemcpy(aermr01.data(), aermr01_g, nmemsize, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(aermr02.data(), aermr02_g, nmemsize, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(aermr03.data(), aermr03_g, nmemsize, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(aermr04.data(), aermr04_g, nmemsize, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(aermr05.data(), aermr05_g, nmemsize, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(aermr06.data(), aermr06_g, nmemsize, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(aermr07.data(), aermr07_g, nmemsize, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(aermr08.data(), aermr08_g, nmemsize, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(aermr09.data(), aermr09_g, nmemsize, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(aermr10.data(), aermr10_g, nmemsize, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(aermr11.data(), aermr11_g, nmemsize, cudaMemcpyDeviceToHost));
    }
}
#endif

template class Aerosol<double>;
template class Aerosol<float>;