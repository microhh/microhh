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

#ifndef MICROPHYS_SEDI_KERNELS_H
#define MICROPHYS_SEDI_KERNELS_H

namespace Micro_sedimentation_kernels
{
    template<typename TF> __global__
    void set_bc_g(
            TF* const __restrict__ w,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ijkb = i + j*icells + kstart*ijcells;
            const int ijkt = i + j*icells + kend  *ijcells;
            const int kk = ijcells;

            // Constant fall speed from lowest grid point to surface (?)
            w[ijkb-kk] = w[ijkb];

            // Zero velocity in top ghost cell
            w[ijkt] = TF(0.);
        }
    }

    template<typename TF> __global__
    void calc_cfl_g(
            TF* const __restrict__ cfl,
            const TF* const __restrict__ w,
            const TF* const __restrict__ dzi,
            const TF dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;
        const int kk = ijcells;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*icells + k*ijcells;

            cfl[ijk] = TF(0.25) * (w[ijk-kk] + w[ijk] + w[ijk] + w[ijk+kk]) * dzi[k] * dt;
        }
    }

    template<typename TF> __device__
    inline TF minmod(const TF a, const TF b)
    {
        return copysign(TF(1.), a) * fmax(TF(0.), fmin(fabs(a), TF(copysign(TF(1.), a))*b));
    }

    template<typename TF> __global__
    void calc_slope_g(
            TF* const __restrict__ slope,
            const TF* const __restrict__ fld,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;
        const int kk = ijcells;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*icells + k*ijcells;

            slope[ijk] = minmod(fld[ijk]-fld[ijk-kk], fld[ijk+kk]-fld[ijk]);
        }
    }

    template<typename TF> __global__
    void calc_flux_g(
            TF* const __restrict__ flux,
            const TF* const __restrict__ fld,
            const TF* const __restrict__ slope,
            const TF* const __restrict__ cfl,
            const TF* const __restrict__ dz,
            const TF* const __restrict__ dzi,
            const TF* const __restrict__ rho,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        if (i < iend && j < jend && k <= kend)
        {
            const int ijk = i + j*icells + k*ijcells;

            if (k == kend || cfl[ijk] < TF(1e-3))
                flux[ijk] = TF(0.);
            else
            {
                int kk;
                TF dzz, cc;

                kk         = k;      // current grid level
                flux[ijk]  = TF(0);  // cumulative 'flux'
                dzz        = TF(0);  // distance from zh[k]
                cc         = fmin(TF(1), cfl[ijk]);
                while (cc > 0 && kk < kend)
                {
                    const int ijk2 = i + j*icells + kk*ijcells;

                    flux[ijk] += rho[kk] * (fld[ijk2] + TF(0.5) * slope[ijk2] * (TF(1.)-cc)) * cc * dz[kk];

                    dzz += dz[kk];
                    kk  += 1;
                    cc   = min(TF(1.), cfl[ijk2] - dzz*dzi[kk]);
                }
            }
        }
    }

    template<typename TF> __global__
    void limit_flux_g(
            TF* const __restrict__ flux,
            const TF* const __restrict__ fld,
            const TF* const __restrict__ dz,
            const TF* const __restrict__ rho,
            const TF dt,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            for (int k=kend-1; k>kstart-1; k--)
            {
                const int ijk = i + j*icells + k*ijcells;

                const TF ftot = fmin(flux[ijk], rho[k]*dz[k]*fld[ijk] - flux[ijk+ijcells]*dt);
                flux[ijk] = -ftot / dt;
            }
        }
    }

    template<typename TF> __global__
    void calc_flux_div_g(
            TF* const __restrict__ fldt,
            const TF* const __restrict__ flux,
            const TF* const __restrict__ dzi,
            const TF* const __restrict__ rho,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart, const int kend,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;
        const int kk = ijcells;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*icells + k*ijcells;

            fldt[ijk] += -(flux[ijk+kk] - flux[ijk]) / rho[k] * dzi[k];
        }
    }

    template<typename TF> __global__
    void copy_precip_rate_bot_g(
            TF* const __restrict__ rate_bot,
            const TF* const __restrict__ flux,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*icells;
            const int ijk = ij + kstart*ijcells;

            rate_bot[ij] = -flux[ijk];
        }
    }
}
#endif
