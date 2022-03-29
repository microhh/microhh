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

#ifndef BOUNDARY_SURFACE_KERNELS_GPU_H
#define BOUNDARY_SURFACE_KERNELS_GPU_H

#include "monin_obukhov.h"
#include "boundary.h"

namespace Boundary_surface_kernels_g
{
    namespace most = Monin_obukhov;

    /* Calculate absolute wind speed */
    template<typename TF> __global__
    void calc_dutot_g(
            TF* const __restrict__ dutot,
            const TF* const __restrict__ u,
            const TF* const __restrict__ v,
            const TF* const __restrict__ ubot,
            const TF* const __restrict__ vbot,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ii  = 1;
            const int ii2 = 2;
            const int jj  = icells;
            const int jj2 = 2*icells;

            const int ij  = i + j*icells;
            const int ijk = i + j*icells + kstart*ijcells;
            const TF minval = 1.e-1;

            const TF u_filtered = TF(1./9) *
                ( TF(0.5)*u[ijk-ii-jj] + u[ijk-jj] + u[ijk+ii-jj] + TF(0.5)*u[ijk+ii2-jj]
                + TF(0.5)*u[ijk-ii   ] + u[ijk   ] + u[ijk+ii   ] + TF(0.5)*u[ijk+ii2   ]
                + TF(0.5)*u[ijk-ii+jj] + u[ijk+jj] + u[ijk+ii+jj] + TF(0.5)*u[ijk+ii2+jj] );

            const TF v_filtered = TF(1./9) *
                ( TF(0.5)*v[ijk-ii-jj] + v[ijk-ii] + v[ijk-ii+jj] + TF(0.5)*v[ijk-ii+jj2]
                + TF(0.5)*v[ijk   -jj] + v[ijk   ] + v[ijk   +jj] + TF(0.5)*v[ijk   +jj2]
                + TF(0.5)*v[ijk+ii-jj] + v[ijk+ii] + v[ijk+ii+jj] + TF(0.5)*v[ijk+ii+jj2] );

            const TF du2 = fm::pow2(u_filtered - TF(0.5)*(ubot[ij] + ubot[ij+ii]))
                         + fm::pow2(v_filtered - TF(0.5)*(vbot[ij] + vbot[ij+jj]));

            // Prevent the absolute wind gradient from reaching values less than 0.01 m/s,
            // otherwise evisc at k = kstart blows up
            dutot[ij] = fmax(std::pow(du2, TF(0.5)), minval);
        }
    }

    template<typename TF> __global__
    void calc_duvdz_mo_g(
            TF* const __restrict__ dudz,
            TF* const __restrict__ dvdz,
            const TF* const __restrict__ u,
            const TF* const __restrict__ v,
            const TF* const __restrict__ ubot,
            const TF* const __restrict__ vbot,
            const TF* const __restrict__ ufluxbot,
            const TF* const __restrict__ vfluxbot,
            const TF* const __restrict__ ustar,
            const TF* const __restrict__ obuk,
            const TF* const __restrict__ z0m,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int kstart,
            const int icells, const int ijcells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        const int ii = 1;
        const int jj = icells;

        if (i < iend && j < jend)
        {
            const int ij  = i + j*icells;
            const int ijk = i + j*icells + kstart*ijcells;

            const TF du_c = TF(0.5)*((u[ijk] - ubot[ij]) + (u[ijk+ii] - ubot[ij+ii]));
            const TF dv_c = TF(0.5)*((v[ijk] - vbot[ij]) + (v[ijk+jj] - vbot[ij+jj]));

            const TF ufluxbot = -du_c * ustar[ij] * most::fm(zsl, z0m[ij], obuk[ij]);
            const TF vfluxbot = -dv_c * ustar[ij] * most::fm(zsl, z0m[ij], obuk[ij]);

            dudz[ij] = -ufluxbot / (Constants::kappa<TF> * zsl * ustar[ij]) * most::phim(zsl/obuk[ij]);
            dvdz[ij] = -vfluxbot / (Constants::kappa<TF> * zsl * ustar[ij]) * most::phim(zsl/obuk[ij]);
        }
    }

    template<typename TF> __global__
    void calc_dbdz_mo_g(
            TF* const __restrict__ dbdz,
            const TF* const __restrict__ bfluxbot,
            const TF* const __restrict__ ustar,
            const TF* const __restrict__ obuk,
            const TF zsl,
            const int istart, const int iend,
            const int jstart, const int jend,
            const int icells)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;

        if (i < iend && j < jend)
        {
            const int ij = i + j*icells;
            dbdz[ij] = -bfluxbot[ij] / (Constants::kappa<TF> * zsl * ustar[ij]) * most::phih(zsl/obuk[ij]);
        }
    }

    template<typename TF> __device__
    TF find_obuk_g(
            const float* const __restrict__ zL,
            const float* const __restrict__ f,
            int &n,
            const TF Ri,
            const TF zsl)
    {
        // Determine search direction.
        if ((f[n]-Ri) > 0.f)
            while ( (f[n-1]-Ri) > 0.f && n > 0) { --n; }
        else
            while ( (f[n]-Ri) < 0.f && n < (nzL_lut-1) ) { ++n; }

        const TF zL0 = (n == 0 || n == nzL_lut-1) ? zL[n] : zL[n-1] + (Ri-f[n-1]) / (f[n]-f[n-1]) * (zL[n]-zL[n-1]);

        return zsl/zL0;
    }

    template<typename TF> __device__
    TF calc_obuk_noslip_dirichlet_lookup_g(
            const float* const __restrict__ zL,
            const float* const __restrict__ f,
            int& n,
            const TF du,
            const TF db,
            const TF zsl)
    {
        // Calculate the appropriate Richardson number.
        const TF Ri = Constants::kappa<TF> * db * zsl / fm::pow2(du);
        return find_obuk_g(zL, f, n, Ri, zsl);
    }
}
#endif
