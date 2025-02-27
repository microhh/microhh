#include "grid_kernels.h"
#include "master.h"


namespace Grid_kernels
{
    template<typename TF, bool use_gpu>
    TF calc_mean_kernel(
            const TF* const __restrict__ a,
            const TF* const __restrict__ dz, const TF zsize,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride, const int ijtot,
            const Master& master)
    {
        TF a_mean = static_cast<TF>(0);

        #pragma acc parallel loop gang deviceptr(a, dz) if (use_gpu)
        for (int k=kstart; k<kend; ++k)
        {
            // Use always double precision for accuracy reasons, otherwise noisy stats.
            double a_mean_tmp = static_cast<TF>(0);

            #pragma acc loop collapse(2) reduction(+:a_mean_tmp)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jstride + k*kstride;
                    a_mean_tmp += a[ijk];
                }

            #pragma acc atomic
            a_mean += static_cast<TF>(a_mean_tmp * dz[k] / zsize / static_cast<double>(ijtot));
        }

        // Sum over the MPI processes.
        master.sum(&a_mean, 1);

        return a_mean;
    }


    template<typename TF, bool use_gpu>
    void calc_mean_prof_kernel(
            TF* const __restrict__ a_mean,
            const TF* const __restrict__ a,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart_mean, const int kend_mean,
            const int jstride, const int kstride,
            const int ijtot,
            const Master& master)
    {
        #pragma acc parallel loop gang deviceptr(a, a_mean) if (use_gpu)
        for (int k=kstart_mean; k<kend_mean; ++k)
        {
            // Use always double precision for accuracy reasons, otherwise noisy stats.
            double a_mean_tmp = static_cast<TF>(0);

            #pragma acc loop collapse(2) reduction(+:a_mean_tmp)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jstride + k*kstride;
                    a_mean_tmp += a[ijk];
                }

            a_mean[k] = static_cast<TF>(a_mean_tmp / static_cast<double>(ijtot));
        }

        // Average over the MPI processes.
        master.sum(&a_mean[kstart_mean], kend_mean - kstart_mean);
    }


    /*
    template<typename TF, bool use_gpu = true>
    inline void calc_mean_masked_prof_kernel(
            TF* const __restrict__ a_mean,
            const TF* const __restrict__ a,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart_mean, const int kend_mean,
            const int jstride, const int kstride,
            const int ijtot,
            const Master& master)
    {
        #pragma acc parallel loop gang deviceptr(a, a_mean) if (use_gpu)
        for (int k=kstart_mean; k<kend_mean; ++k)
        {
            // Use always double precision for accuracy reasons, otherwise noisy stats.
            double a_mean_tmp = static_cast<TF>(0);

            #pragma acc loop collapse(2) reduction(+:a_mean_tmp)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jstride + k*kstride;
                    a_mean_tmp += a[ijk];
                }

            a_mean[k] = static_cast<TF>(a_mean_tmp / static_cast<double>(ijtot));
        }

        // Average over the MPI processes.
        master.sum(&a_mean[kstart_mean], kend_mean - kstart_mean);
    }*/


    template<typename TF, bool use_gpu>
    TF calc_max_kernel(
            const TF* const __restrict__ a,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride,
            const Master& master)
    {
        TF a_max = static_cast<TF>(0);

        #pragma acc parallel loop collapse(3) deviceptr(a) reduction(max:a_max) if (use_gpu)
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*jstride + k*kstride;
                    a_max = max(a_max, a[ijk]);
                }

        // Max over the MPI processes.
        master.max(&a_max, 1);

        return a_max;
    }
}


#ifdef FLOAT_SINGLE
template void Grid_kernels::calc_mean_prof_kernel<float, true>(
            float* const __restrict__ a_mean,
            const float* const __restrict__ a,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride,
            const int ijtot,
            const Master& master);

template void Grid_kernels::calc_mean_prof_kernel<float, false>(
            float* const __restrict__ a_mean,
            const float* const __restrict__ a,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride,
            const int ijtot,
            const Master& master);

template float Grid_kernels::calc_max_kernel<float, true>(
            const float* const __restrict__ a,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride,
            const Master& master);

template float Grid_kernels::calc_max_kernel<float, false>(
            const float* const __restrict__ a,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride,
            const Master& master);
#else
template void Grid_kernels::calc_mean_prof_kernel<double, true>(
            double* const __restrict__ a_mean,
            const double* const __restrict__ a,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride,
            const int ijtot,
            const Master& master);

template void Grid_kernels::calc_mean_prof_kernel<double, false>(
            double* const __restrict__ a_mean,
            const double* const __restrict__ a,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride,
            const int ijtot,
            const Master& master);

template double Grid_kernels::calc_max_kernel<double, true>(
            const double* const __restrict__ a,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride,
            const Master& master);

template double Grid_kernels::calc_max_kernel<double, false>(
            const double* const __restrict__ a,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int jstride, const int kstride,
            const Master& master);
#endif
