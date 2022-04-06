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

#include "timeloop.h"
#include "grid.h"
#include "soil_grid.h"
#include "master.h"
#include "fields.h"
#include "soil_field3d.h"
#include "constants.h"
#include "tools.h"


namespace
{
    template<typename TF, int substep> __global__
    void rk3_g(TF* __restrict__ a, TF* __restrict__ at, double dt,
               const int jj, const int kk,
               const int istart, const int jstart, const int kstart,
               const int iend,   const int jend,   const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        constexpr TF cA1 = -5./9.;
        constexpr TF cA2 = -153./128.;

        constexpr TF cB0 =  1./ 3.;
        constexpr TF cB1 = 15./16.;
        constexpr TF cB2 =  8./15.;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;

            switch (substep)
            {
                case 0:
                    a [ijk] = a[ijk] + cB0*dt*at[ijk];
                    at[ijk] = cA1*at[ijk];
                    break;
                case 1:
                    a [ijk] = a[ijk] + cB1*dt*at[ijk];
                    at[ijk] = cA2*at[ijk];
                    break;
                case 2:
                    a [ijk] = a[ijk] + cB2*dt*at[ijk];
                    at[ijk] = TF(0.);
                    break;
            }
        }
    }

    template<typename TF, int substep> __global__
    void rk4_g(TF* __restrict__ a, TF* __restrict__ at, double dt,
               const int jj, const int kk,
               const int istart, const int jstart, const int kstart,
               const int iend,   const int jend,   const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        constexpr TF cA1 = - 567301805773./1357537059087.;
        constexpr TF cA2 = -2404267990393./2016746695238.;
        constexpr TF cA3 = -3550918686646./2091501179385.;
        constexpr TF cA4 = -1275806237668./ 842570457699.;

        constexpr TF cB0 = 1432997174477./ 9575080441755.;
        constexpr TF cB1 = 5161836677717./13612068292357.;
        constexpr TF cB2 = 1720146321549./ 2090206949498.;
        constexpr TF cB3 = 3134564353537./ 4481467310338.;
        constexpr TF cB4 = 2277821191437./14882151754819.;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;

            switch (substep)
            {
                case 0:
                    a [ijk] = a[ijk] + cB0*dt*at[ijk];
                    at[ijk] = cA1*at[ijk];
                    break;
                case 1:
                    a [ijk] = a[ijk] + cB1*dt*at[ijk];
                    at[ijk] = cA2*at[ijk];
                    break;
                case 2:
                    a [ijk] = a[ijk] + cB2*dt*at[ijk];
                    at[ijk] = cA3*at[ijk];
                    break;
                case 3:
                    a [ijk] = a[ijk] + cB3*dt*at[ijk];
                    at[ijk] = cA4*at[ijk];
                    break;
                case 4:
                    a [ijk] = a[ijk] + cB4*dt*at[ijk];
                    at[ijk] = 0;
                    break;
            }
        }
    }
}

#ifdef USECUDA
template<typename TF>
void Timeloop<TF>::exec()
{
    auto& gd = grid.get_grid_data();
    auto& sgd = soil_grid.get_grid_data();

    // Time integration 2D fields with 3D routines:
    const int kstart_2d = 0;
    const int kend_2d = 1;

    if (rkorder == 3)
    {
        auto rk3_substep = [&](
            TF* const __restrict__ tend,
            TF* const __restrict__ fld,
            const int kstart, const int kend)
        {
            const int kmax = kend-kstart;

            const int blocki = gd.ithread_block;
            const int blockj = gd.jthread_block;
            const int gridi = gd.imax/blocki + (gd.imax%blocki > 0);
            const int gridj = gd.jmax/blockj + (gd.jmax%blockj > 0);

            dim3 gridGPU (gridi, gridj, kmax);
            dim3 blockGPU(blocki, blockj, 1);

            if (substep == 0)
                rk3_g<TF, 0><<<gridGPU, blockGPU>>>(
                    tend, fld, dt,
                    gd.icells, gd.ijcells,
                    gd.istart,  gd.jstart, kstart,
                    gd.iend,    gd.jend,   kend);
            else if (substep == 1)
                rk3_g<TF, 1><<<gridGPU, blockGPU>>>(
                    tend, fld, dt,
                    gd.icells, gd.ijcells,
                    gd.istart,  gd.jstart, kstart,
                    gd.iend,    gd.jend,   kend);
            else if (substep == 2)
                rk3_g<TF, 2><<<gridGPU, blockGPU>>>(
                    tend, fld, dt,
                    gd.icells, gd.ijcells,
                    gd.istart,  gd.jstart, kstart,
                    gd.iend,    gd.jend,   kend);
        };

        // Atmospheric fields
        for (auto& f : fields.at)
            rk3_substep(fields.ap.at(f.first)->fld_g, f.second->fld_g, gd.kstart, gd.kend);

        // Soil fields
        for (auto& f : fields.sts)
            rk3_substep(fields.sps.at(f.first)->fld_g, f.second->fld_g, sgd.kstart, sgd.kend);

        // 2D fields
        for (auto& f : fields.at2d)
            rk3_substep(fields.ap2d.at(f.first)->fld_g, f.second->fld_g, kstart_2d, kend_2d);

        substep = (substep+1) % 3;

        /*
           rk3_kernel<<<gridGPU, blockGPU>>>(a, at, dt,
           substep, gd.icells, gd.ijcells,
           gd.istart, gd.jstart, gd.kstart,
           gd.iend, gd.jend, gd.kend);
         */
    }

    //else if (rkorder == 4)
    //{
    //    for (auto& it : fields.at)
    //    {
    //        if (substep==0)
    //            rk4_g<TF,0><<<gridGPU, blockGPU>>>(
    //                fields.ap.at(it.first)->fld_g, it.second->fld_g, dt,
    //                gd.icells, gd.ijcells,
    //                gd.istart,  gd.jstart, gd.kstart,
    //                gd.iend,    gd.jend,   gd.kend);
    //        else if (substep==1)
    //            rk4_g<TF,1><<<gridGPU, blockGPU>>>(
    //                fields.ap.at(it.first)->fld_g, it.second->fld_g, dt,
    //                gd.icells, gd.ijcells,
    //                gd.istart,  gd.jstart, gd.kstart,
    //                gd.iend,    gd.jend,   gd.kend);
    //        else if (substep==2)
    //            rk4_g<TF,2><<<gridGPU, blockGPU>>>(
    //                fields.ap.at(it.first)->fld_g, it.second->fld_g, dt,
    //                gd.icells, gd.ijcells,
    //                gd.istart,  gd.jstart, gd.kstart,
    //                gd.iend,    gd.jend,   gd.kend);
    //        else if (substep==3)
    //            rk4_g<TF,3><<<gridGPU, blockGPU>>>(
    //                fields.ap.at(it.first)->fld_g, it.second->fld_g, dt,
    //                gd.icells, gd.ijcells,
    //                gd.istart,  gd.jstart, gd.kstart,
    //                gd.iend,    gd.jend,   gd.kend);
    //        else if (substep==4)
    //            rk4_g<TF,4><<<gridGPU, blockGPU>>>(
    //                fields.ap.at(it.first)->fld_g, it.second->fld_g, dt,
    //                gd.icells, gd.ijcells,
    //                gd.istart,  gd.jstart, gd.kstart,
    //                gd.iend,    gd.jend,   gd.kend);
    //    }

    //    substep = (substep+1) % 5;

    //    /*
    //       rk4_kernel<<<gridGPU, blockGPU>>>(a, at, dt,
    //       substep, gd.icells, gd.ijcells,
    //       gd.istart, gd.jstart, gd.kstart,
    //       gd.iend, gd.jend, gd.kend);
    //     */
    //}

    cuda_check_error();
}
#endif

template class Timeloop<double>;
template class Timeloop<float>;
