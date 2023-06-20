#include <algorithm>
#include "tools.h"
#include "timedep.h"
#include "background_profs.h"

template<typename TF>
void Background<TF>::prepare_device()
{
    if (!sw_update_background)
        return;

    const int nmemsize_lay  = n_era_layers*sizeof(TF);
    const int nmemsize_lev  = n_era_levels*sizeof(TF);

    cuda_safe_call(cudaMalloc(&t_lay_g, nmemsize_lay));
    cuda_safe_call(cudaMalloc(&t_lev_g, nmemsize_lev));
    cuda_safe_call(cudaMalloc(&p_lay_g, nmemsize_lay));
    cuda_safe_call(cudaMalloc(&p_lev_g, nmemsize_lev));
    cuda_safe_call(cudaMalloc(&h2o_g, nmemsize_lay));
    cuda_safe_call(cudaMemcpy(t_lay_g, t_lay.data(), nmemsize_lay, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(t_lev_g, t_lev.data(), nmemsize_lev, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(p_lay_g, p_lay.data(), nmemsize_lay, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(p_lev_g, p_lev.data(), nmemsize_lev, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(h2o_g, h2o.data(), nmemsize_lay, cudaMemcpyHostToDevice));

    for (auto& it : gaslist)
    {
        gasprofs_g.emplace(it, nullptr);
        cuda_safe_call(cudaMalloc(&gasprofs_g.at(it), nmemsize_lay));
        cuda_safe_call(cudaMemcpy(gasprofs_g.at(it), gasprofs.at(it).data(), nmemsize_lay, cudaMemcpyHostToDevice));
    }

    cuda_safe_call(cudaMalloc(&aermr01_g, nmemsize_lay));
    cuda_safe_call(cudaMalloc(&aermr02_g, nmemsize_lay));
    cuda_safe_call(cudaMalloc(&aermr03_g, nmemsize_lay));
    cuda_safe_call(cudaMalloc(&aermr04_g, nmemsize_lay));
    cuda_safe_call(cudaMalloc(&aermr05_g, nmemsize_lay));
    cuda_safe_call(cudaMalloc(&aermr06_g, nmemsize_lay));
    cuda_safe_call(cudaMalloc(&aermr07_g, nmemsize_lay));
    cuda_safe_call(cudaMalloc(&aermr08_g, nmemsize_lay));
    cuda_safe_call(cudaMalloc(&aermr09_g, nmemsize_lay));
    cuda_safe_call(cudaMalloc(&aermr10_g, nmemsize_lay));
    cuda_safe_call(cudaMalloc(&aermr11_g, nmemsize_lay));
    cuda_safe_call(cudaMemcpy(aermr01_g, aermr01.data(), nmemsize_lay, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(aermr02_g, aermr02.data(), nmemsize_lay, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(aermr03_g, aermr03.data(), nmemsize_lay, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(aermr04_g, aermr04.data(), nmemsize_lay, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(aermr05_g, aermr05.data(), nmemsize_lay, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(aermr06_g, aermr06.data(), nmemsize_lay, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(aermr07_g, aermr07.data(), nmemsize_lay, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(aermr08_g, aermr08.data(), nmemsize_lay, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(aermr09_g, aermr09.data(), nmemsize_lay, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(aermr10_g, aermr10.data(), nmemsize_lay, cudaMemcpyHostToDevice));
    cuda_safe_call(cudaMemcpy(aermr11_g, aermr11.data(), nmemsize_lay, cudaMemcpyHostToDevice));
}

template<typename TF>
void Background<TF>::clear_device()
{
    if (!sw_update_background)
        return;

    cuda_safe_call(cudaFree(t_lay_g));
    cuda_safe_call(cudaFree(t_lev_g));
    cuda_safe_call(cudaFree(p_lay_g));
    cuda_safe_call(cudaFree(p_lev_g));
    cuda_safe_call(cudaFree(h2o_g));
    tdep_t_lay->clear_device();
    tdep_t_lev->clear_device();
    tdep_p_lay->clear_device();
    tdep_p_lev->clear_device();
    tdep_h2o->clear_device();

    for (auto& it : gasprofs_g)
        cuda_safe_call(cudaFree(it.second));

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

#ifdef USECUDA
template <typename TF>
void Background<TF>::update_time_dependent(Timeloop<TF>& timeloop)
{
    if (!sw_update_background)
        return;

    const bool do_radiation = ((timeloop.get_itime() % idt_rad == 0) && !timeloop.in_substep()) ;
    const int nmemsize_lay  = n_era_layers*sizeof(TF);
    const int nmemsize_lev = n_era_levels*sizeof(TF);

    if (do_radiation)
    {
        // temperature, pressure and moisture
        tdep_t_lay   ->update_time_dependent_background_prof_g(t_lay_g, timeloop, n_era_layers);
        tdep_t_lev   ->update_time_dependent_background_prof_g(t_lev_g, timeloop, n_era_levels);
        tdep_p_lay   ->update_time_dependent_background_prof_g(p_lay_g, timeloop, n_era_layers);
        tdep_p_lev   ->update_time_dependent_background_prof_g(p_lev_g, timeloop, n_era_levels);
        tdep_h2o     ->update_time_dependent_background_prof_g(h2o_g, timeloop, n_era_layers);
        cuda_safe_call(cudaMemcpy(t_lay.data(), t_lay_g, nmemsize_lay, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(t_lev.data(), t_lev_g, nmemsize_lev, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(p_lay.data(), p_lay_g, nmemsize_lay, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(p_lev.data(), p_lev_g, nmemsize_lev, cudaMemcpyDeviceToHost));
        cuda_safe_call(cudaMemcpy(h2o.data(), h2o_g, nmemsize_lay, cudaMemcpyDeviceToHost));

        // gasses
        for (auto& it : tdep_gases)
        {
            it.second->update_time_dependent_background_prof_g(gasprofs_g.at(it.first), timeloop, n_era_layers);
            cuda_safe_call(cudaMemcpy(gasprofs.at(it.first).data(), gasprofs_g.at(it.first), nmemsize_lay, cudaMemcpyDeviceToHost));
        }

        // aerosols
        if (sw_aerosol && sw_aerosol_timedep)
        {
            tdep_aermr01 ->update_time_dependent_background_prof_g(aermr01_g, timeloop, n_era_layers);
            tdep_aermr02 ->update_time_dependent_background_prof_g(aermr02_g, timeloop, n_era_layers);
            tdep_aermr03 ->update_time_dependent_background_prof_g(aermr03_g, timeloop, n_era_layers);
            tdep_aermr04 ->update_time_dependent_background_prof_g(aermr04_g, timeloop, n_era_layers);
            tdep_aermr05 ->update_time_dependent_background_prof_g(aermr05_g, timeloop, n_era_layers);
            tdep_aermr06 ->update_time_dependent_background_prof_g(aermr06_g, timeloop, n_era_layers);
            tdep_aermr07 ->update_time_dependent_background_prof_g(aermr07_g, timeloop, n_era_layers);
            tdep_aermr08 ->update_time_dependent_background_prof_g(aermr08_g, timeloop, n_era_layers);
            tdep_aermr09 ->update_time_dependent_background_prof_g(aermr09_g, timeloop, n_era_layers);
            tdep_aermr10 ->update_time_dependent_background_prof_g(aermr10_g, timeloop, n_era_layers);
            tdep_aermr11 ->update_time_dependent_background_prof_g(aermr11_g, timeloop, n_era_layers);

            cuda_safe_call(cudaMemcpy(aermr01.data(), aermr01_g, nmemsize_lay, cudaMemcpyDeviceToHost));
            cuda_safe_call(cudaMemcpy(aermr02.data(), aermr02_g, nmemsize_lay, cudaMemcpyDeviceToHost));
            cuda_safe_call(cudaMemcpy(aermr03.data(), aermr03_g, nmemsize_lay, cudaMemcpyDeviceToHost));
            cuda_safe_call(cudaMemcpy(aermr04.data(), aermr04_g, nmemsize_lay, cudaMemcpyDeviceToHost));
            cuda_safe_call(cudaMemcpy(aermr05.data(), aermr05_g, nmemsize_lay, cudaMemcpyDeviceToHost));
            cuda_safe_call(cudaMemcpy(aermr06.data(), aermr06_g, nmemsize_lay, cudaMemcpyDeviceToHost));
            cuda_safe_call(cudaMemcpy(aermr07.data(), aermr07_g, nmemsize_lay, cudaMemcpyDeviceToHost));
            cuda_safe_call(cudaMemcpy(aermr08.data(), aermr08_g, nmemsize_lay, cudaMemcpyDeviceToHost));
            cuda_safe_call(cudaMemcpy(aermr09.data(), aermr09_g, nmemsize_lay, cudaMemcpyDeviceToHost));
            cuda_safe_call(cudaMemcpy(aermr10.data(), aermr10_g, nmemsize_lay, cudaMemcpyDeviceToHost));
            cuda_safe_call(cudaMemcpy(aermr11.data(), aermr11_g, nmemsize_lay, cudaMemcpyDeviceToHost));
        }
    }

}
#endif

template class Background<double>;
template class Background<float>;
