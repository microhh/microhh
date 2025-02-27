#pragma once


class Master;


namespace Grid_kernels
{
    template<typename TF, bool use_gpu = true>
    TF calc_mean_kernel(
            const TF* const __restrict__ a,
            const TF* const __restrict__ dz, const TF zsize,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride, const int ijtot,
            const Master& master);


    template<typename TF, bool use_gpu = true>
    void calc_mean_prof_kernel(
            TF* const __restrict__ a_mean,
            const TF* const __restrict__ a,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride,
            const int ijtot,
            const Master& master);


    template<typename TF, bool use_gpu = true>
    TF calc_max_kernel(
            const TF* const __restrict__ a,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride,
            const Master& master);
}
