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

#include "constants.h"

#ifdef __CUDACC__
#  define CUDA_MACRO __host__ __device__
#else
#  define CUDA_MACRO
#endif


#ifndef STATS_FUNCTIONS_H
#define STATS_FUNCTIONS_H

namespace Stats_functions
{
    using namespace Constants;

    // // Help function(s) to switch between the different NetCDF data types
    // template<typename TF> TF netcdf_fp_fillvalue();
    // template<> double netcdf_fp_fillvalue<double>() { return NC_FILL_DOUBLE; }
    // template<> float  netcdf_fp_fillvalue<float>()  { return NC_FILL_FLOAT; }

    template<typename TF>
    TF netcdf_fp_fillvalue()
    {
        if (sizeof(TF) == 4)
            return NC_FILL_FLOAT;
        else
            return NC_FILL_DOUBLE;
    }

    template<typename TF, Stats_mask_type mode> CUDA_MACRO
    TF is_false(const TF value, const TF threshold)
    {
        if (mode == Stats_mask_type::Plus)
            return (value <= threshold);
        else if (mode == Stats_mask_type::Min)
            return (value > threshold);
    }

    template<typename TF>
    TF in_mask(const unsigned int mask, const unsigned int flag)
    {
        return static_cast<TF>( (mask & flag) != 0 );
    }

    template<typename TF>
    void set_flag(unsigned int& flag, const int*& restrict nmask, const Mask<TF>& m, const int loc)
    {
        if (loc == 0)
        {
            flag = m.flag;
            nmask = m.nmask.data();
        }
        else
        {
            flag = m.flagh;
            nmask = m.nmaskh.data();
        }
    }

    template<typename TF>
    void set_fillvalue_prof(TF* const restrict data, const int* const restrict nmask, const int kstart, const int kcells)
    {
        for (int k=0; k<kcells; ++k)
        {
            if (nmask[k] == 0)
                data[k] = 77777;
        }
    }

    template<typename TF, Stats_mask_type mode>
    void calc_mask_thres(
            unsigned int* const restrict mfield, unsigned int* const restrict mfield_bot,
            const unsigned int flag, const unsigned int flagh,
            const TF* const restrict fld, const TF* const restrict fldh,
            const TF* const restrict fld_bot, const TF threshold,
            const int istart, const int jstart, const int kstart,
            const int iend, const int jend, const int kend,
            const int icells, const int ijcells, const int kcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    mfield[ijk] -= (mfield[ijk] & flag ) * is_false<TF, mode>(fld [ijk], threshold);
                    mfield[ijk] -= (mfield[ijk] & flagh) * is_false<TF, mode>(fldh[ijk], threshold);
                }

        // Set the top value for the flux level.
        #pragma omp parallel for
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ijk = i + j*icells + kend*ijcells;
                mfield[ijk] -= (mfield[ijk] & flagh) * is_false<TF, mode>(fldh[ijk], threshold);
            }

        // Set the ghost cells equal to the first model level.
        #pragma omp parallel for
        for (int k=0; k<kstart; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    const int ijk_ref = i + j*icells + kstart*ijcells;
                    mfield[ijk] -= (mfield[ijk] & flag ) * is_false<TF, mode>(fld [ijk_ref], threshold);
                    mfield[ijk] -= (mfield[ijk] & flagh) * is_false<TF, mode>(fldh[ijk_ref], threshold);
                }

        // Set the ghost cells for the full level equal to kend-1.
        #pragma omp parallel for
        for (int k=kend; k<kcells; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    const int ijk_ref = i + j*icells + (kend-1)*ijcells;
                    mfield[ijk] -= (mfield[ijk] & flag) * is_false<TF, mode>(fld [ijk_ref], threshold);
                }

        // Set the ghost cells for the flux level equal to kend.
        #pragma omp parallel for
        for (int k=kend+1; k<kcells; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    const int ijk_ref = i + j*icells + kend*ijcells;
                    mfield[ijk] -= (mfield[ijk] & flagh) * is_false<TF, mode>(fldh[ijk_ref], threshold);
                }

        // Set the mask for surface projected quantities
        #pragma omp parallel for
        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij = i + j*icells;
                mfield_bot[ij] -= (mfield_bot[ij] & flag) * is_false<TF, mode>(fld_bot[ij], threshold);
            }
    }

#ifdef __CUDACC__
    template<typename TF, Stats_mask_type mode> __global__
    void calc_mask_thres_g(
            unsigned int* const restrict mfield, unsigned int* const restrict mfield_bot,
            const unsigned int flag, const unsigned int flagh,
            const TF* const restrict fld, const TF* const restrict fldh,
            const TF* const restrict fld_bot, const TF threshold,
            const int istart, const int jstart, const int kstart,
            const int iend, const int jend, const int kend,
            const int icells, const int jcells, const int ijcells, const int kcells)
    {
        const int i  = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j  = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k  = blockIdx.z + kstart;

        if (i < icells && j < jcells && k < kcells)
        {
            const int ijk = i + j*icells + k*ijcells;
            mfield[ijk] -= (mfield[ijk] & flag ) * is_false<TF, mode>(fld [ijk], threshold);
            mfield[ijk] -= (mfield[ijk] & flagh) * is_false<TF, mode>(fldh[ijk], threshold);
            // Set the top value for the flux level.
            if (k == kend)
            {
                mfield[ijk] -= (mfield[ijk] & flagh) * is_false<TF, mode>(fldh[ijk], threshold);
                const int ijk_ref = i + j*icells + (kend-1)*ijcells;
                mfield[ijk] -= (mfield[ijk] & flag) * is_false<TF, mode>(fld [ijk_ref], threshold);
            }
            // Set the ghost cells equal to the first model level.
            if (k < kstart)
            {
                const int ijk_ref = i + j*icells + kstart*ijcells;
                mfield[ijk] -= (mfield[ijk] & flag ) * is_false<TF, mode>(fld [ijk_ref], threshold);
                mfield[ijk] -= (mfield[ijk] & flagh) * is_false<TF, mode>(fldh[ijk_ref], threshold);
            }
            // Set the ghost cells for the full level equal to kend-1.
            // Set the ghost cells for the flux level equal to kend.
            if (k > kend)
            {
                const int ijk_ref = i + j*icells + (kend-1)*ijcells;
                const int ijk_refh = i + j*icells + kend*ijcells;
                mfield[ijk] -= (mfield[ijk] & flag) * is_false<TF, mode>(fld [ijk_ref], threshold);
                mfield[ijk] -= (mfield[ijk] & flagh) * is_false<TF, mode>(fldh[ijk_refh], threshold);
            }
            // Set the mask for surface projected quantities
            if (k == kstart)
            {
                const int ij = i + j*icells;
                mfield_bot[ij] -= (mfield_bot[ij] & flag) * is_false<TF, mode>(fld_bot[ij], threshold);
            }
        }
    }
#endif

    template<typename TF>
    void calc_area(
            TF* const restrict area, const int loc[3], const int* const restrict nmask,
            const int kstart, const int kend, const int ijtot)
    {
        for (int k=kstart; k<kend+loc[2]; ++k)
        {
            if (nmask[k])
                area[k] = static_cast<TF>(nmask[k]) / static_cast<TF>(ijtot);
            else
                area[k] = 0.;
        }
    }

    // Calculate the number of points contained in the mask.
    template<typename TF>
    void calc_nmask(
            int* restrict nmask_full, int* restrict nmask_half, int& nmask_bottom,
            const unsigned int* const mfield, const unsigned int* const mfield_bot,
            const unsigned int flag, const unsigned int flagh,
            const int istart, const int iend, const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells, const int kcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
        {
            nmask_full[k] = 0;
            nmask_half[k] = 0;

            // #pragma omp parallel for reduction (+:nmask_full[k], nmask_half[k]) collapse(2)
            for (int j=jstart; j<jend; ++j)
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    nmask_full[k] += in_mask<int>(mfield[ijk], flag );
                    nmask_half[k] += in_mask<int>(mfield[ijk], flagh);
                }
        }

        nmask_bottom     = 0;
        // nmask_half[kend] = 0;
        // #pragma omp parallel for reduction (+:nmask_bottom) collapse(2)
        // #pragma omp parallel for reduction (+:nmask_bottom, nmask_half[kend]) collapse(2)
        for (int j=jstart; j<jend; ++j)
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                const int ijk = i + j*icells + kend*ijcells;
                nmask_bottom += in_mask<int>(mfield_bot[ij], flag);
                // nmask_half[kend] += in_mask<int>(mfield[ijk], flagh);
            }
    }

#ifdef __CUDACC__
    // Calculate the number of points contained in the mask.
    template<typename TF> __global__
    void calc_nmask_g(
            int* restrict nmask_full, int* restrict nmask_half, int& nmask_bottom,
            const unsigned int* const mfield, const unsigned int* const mfield_bot,
            const unsigned int flag, const unsigned int flagh,
            const int istart, const int iend, const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells, const int kcells)
    {
    }
#endif
    template<typename TF>
    void calc_mean(
            TF* const restrict prof, const TF* const restrict fld,
            const unsigned int* const mask, const unsigned int flag, const int* const nmask,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        #pragma omp parallel for
        for (int k=kstart-1; k<kend+1; ++k)
        {
            if (nmask[k])
            {
                double tmp = 0.;
                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk  = i + j*icells + k*ijcells;
                        tmp += in_mask<double>(mask[ijk], flag) * fld[ijk];
                    }

                prof[k] = tmp / nmask[k];
            }
        }
    }

    template<typename TF>
    void calc_mean_2d(
            TF& out, const TF* const restrict fld,
            const int istart, const int iend, const int jstart, const int jend,
            const int icells, const int itot, const int jtot)
    {
        double tmp = 0.;

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                const int ij  = i + j*icells;
                tmp += fld[ij];
            }

        out = tmp / (itot*jtot);
    }


    template<typename TF>
    void calc_mean_projected_mask(
            TF* const restrict prof, const TF* const restrict fld,
            const unsigned int* const mask_bot, const int nmask_bot,
            const int flag,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        double tmp;

        if (nmask_bot == 0)
            return;

        for (int k=kstart; k<kend; ++k)
        {
            tmp = 0;
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ij  = i  + j*icells;
                    const int ijk = ij + k*ijcells;

                    tmp += in_mask<double>(mask_bot[ij], flag) * fld[ijk];
                }
            prof[k] = tmp / nmask_bot;
        }
    }


    template<typename TF>
    void calc_moment(
            TF* const restrict prof, const TF* const restrict fld, const TF* const restrict fld_mean, const TF offset,
            const unsigned int* const mask, const unsigned int flag, const int* const nmask, const int power,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend+1; ++k)
        {
            if (nmask[k])
            {
                double tmp = 0.;
                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk  = i + j*icells + k*ijcells;
                        tmp += in_mask<double>(mask[ijk], flag)*std::pow(fld[ijk] - fld_mean[k] + offset, power);
                    }

                prof[k] = tmp / nmask[k];
            }
        }
    }

    template<typename TF>
    void calc_cov(
            TF* const restrict prof, const TF* const restrict fld1, const TF* const restrict fld1_mean,
            const TF offset1, const int pow1,
            const TF* const restrict fld2, const TF* const restrict fld2_mean, const TF offset2, const int pow2,
            const unsigned int* const mask, const unsigned int flag, const int* const nmask,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend+1; ++k)
        {
            if (nmask[k])
            {
                double tmp = 0.;
                if ((fld1_mean[k] != 77777) && (fld2_mean[k] != 77777))
                {
                    for (int j=jstart; j<jend; ++j)
                        #pragma ivdep
                        for (int i=istart; i<iend; ++i)
                        {
                            const int ijk  = i + j*icells + k*ijcells;
                            tmp += in_mask<double>(mask[ijk], flag)
                                * std::pow(fld1[ijk] - fld1_mean[k] + offset1, pow1)
                                * std::pow(fld2[ijk] - fld2_mean[k] + offset2, pow2);
                        }

                    prof[k] = tmp / nmask[k];
                }
            }
        }
    }

    template<typename TF>
    void add_fluxes(
            TF* const restrict flux, const TF* const restrict turb, const TF* const restrict diff,
            const int kstart, const int kend)
    {
        for (int k=kstart; k<kend+1; ++k)
            flux[k] = turb[k] + diff[k];
    }


    template<typename TF>
    void calc_frac(
            TF* const restrict prof, const TF* const restrict fld, const TF offset, const TF threshold,
            const unsigned int* const mask, const unsigned int flag, const int* const nmask,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend+1; ++k)
        {
            if (nmask[k])
            {
                double tmp = 0.;
                for (int j=jstart; j<jend; ++j)
                    #pragma ivdep
                    for (int i=istart; i<iend; ++i)
                    {
                        const int ijk  = i + j*icells + k*ijcells;
                        tmp += in_mask<double>(mask[ijk], flag)*((fld[ijk] + offset) > threshold);
                    }
                prof[k] = tmp / nmask[k];
            }
        }
    }

    template<typename TF>
    std::pair<TF, int> calc_path(
            const TF* const restrict data, const TF* const restrict dz, const TF* const restrict rho,
            const unsigned int* const mask, const unsigned int flag, const int* const nmask,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int jj, const int kk)
    {
        int nmask_proj = 0;
        TF path = TF(0.);

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                for (int k=kstart; k<kend; ++k)
                {
                    const int ijk = i + j*jj + k*kk;
                    if (in_mask<bool>(mask[ijk], flag))
                    {
                        ++nmask_proj;
                        break;
                    }
                }

                if (nmask_proj > 0)
                {
                    for (int k=kstart; k<kend; ++k)
                    {
                        const int ijk = i + j*jj + k*kk;
                        if (in_mask<bool>(mask[ijk], flag))
                        {
                            path += data[ijk]*rho[k]*dz[k];
                        }
                    }
                }
            }

        return std::make_pair(path, nmask_proj);
    }

    template<typename TF>
    std::pair<int, int> calc_cover(
            const TF* const restrict fld, const TF offset, const TF threshold,
            const unsigned int* const mask, const unsigned int flag, const int* const nmask,
            const int istart, const int iend, const int jstart, const int jend, const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        int cover = 0.;
        int nmaskcover = 0;

        for (int j=jstart; j<jend; ++j)
            #pragma ivdep
            for (int i=istart; i<iend; ++i)
            {
                int maskincolumn = 0;
                for (int k=kstart; k<kend; ++k)
                {
                    const int ijk  = i + j*icells + k*ijcells;
                    if (in_mask<bool>(mask[ijk], flag))
                    {
                        maskincolumn = 1;
                        if ((fld[ijk] + offset) > threshold)
                        {
                            ++cover;
                            break;
                        }
                    }
                }
                nmaskcover += maskincolumn;
            }

        return std::make_pair(cover, nmaskcover);
    }

    template<typename TF>
    bool has_only_digits(const std::string& s)
    {
        return s.find_first_not_of( "23456789" ) == std::string::npos;
    }

    template<typename TF>
    void subtract_mean(
            TF* const restrict fld_prime,
            const TF* const restrict fld,
            const TF* const restrict fld_mean,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        #pragma omp parallel for
        for (int k=kstart; k<kend; ++k)
            for (int j=jstart; j<jend; ++j)
                #pragma ivdep
                for (int i=istart; i<iend; ++i)
                {
                    const int ijk = i + j*icells + k*ijcells;
                    fld_prime[ijk] = fld[ijk] - fld_mean[k];
                }
    }
#ifdef __CUDACC__
    template<typename TF> __global__
    void apply_mask_g(
            TF* slice, const TF* const __restrict__ fld,
            const unsigned int* const __restrict__ mask, const unsigned int flag,
            const int icells, const int jcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x;
        const int j = blockIdx.y*blockDim.y + threadIdx.y;

        if (i < icells && j < jcells)
        {
            const int ij = i + j*icells;
            slice[ij] = 1.;//fld[ij] * ( mask[ij] & flag );
        }
    }
#endif
}
#endif