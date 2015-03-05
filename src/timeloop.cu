/*
 * MicroHH
 * Copyright (c) 2011-2015 Chiel van Heerwaarden
 * Copyright (c) 2011-2015 Thijs Heus
 * Copyright (c) 2014-2015 Bart van Stratum
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
#include "master.h"
#include "fields.h"
#include "constants.h"
#include "tools.h"

namespace
{
    /*
       __global__ void rk3_kernel(double * __restrict__ a, double * __restrict__ at, double dt,
       const int substep, const int jj, const int kk,
       const int istart, const int jstart, const int kstart,
       const int iend, const int jend, const int kend)
       {
       const double cA[] = {0., -5./9., -153./128.};
       const double cB[] = {1./3., 15./16., 8./15.};

       const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
       const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
       const int k = blockIdx.z + kstart;

       if(i < iend && j < jend && k < kend)
       {
       const int ijk = i + j*jj + k*kk;
       a[ijk] = a[ijk] + cB[substep]*dt*at[ijk];

       const int substepn = (substep+1) % 3;
    // substep 0 resets the tendencies, because cA[0] == 0
    at[ijk] = cA[substepn]*at[ijk];
    }
    }
     */

    template<int substep> __global__ 
    void rk3_g(double* __restrict__ a, double* __restrict__ at, double dt,
               const int jj, const int kk,
               const int istart, const int jstart, const int kstart,
               const int iend,   const int jend,   const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        // const double cA0 =  0.;
        const double cA1 = -5./9.;
        const double cA2 = -153./128.;

        const double cB0 =  1./ 3.;
        const double cB1 = 15./16.;
        const double cB2 =  8./15.;

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
                    at[ijk] = 0.; 
                    break;
            }
        }
    }

    template<int substep> __global__ 
    void rk4_g(double* __restrict__ a, double* __restrict__ at, double dt,
               const int jj, const int kk,
               const int istart, const int jstart, const int kstart,
               const int iend,   const int jend,   const int kend)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart;
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart;
        const int k = blockIdx.z + kstart;

        // const double cA0 =   0.;
        const double cA1 = - 567301805773./1357537059087.;
        const double cA2 = -2404267990393./2016746695238.;
        const double cA3 = -3550918686646./2091501179385.;
        const double cA4 = -1275806237668./ 842570457699.;

        const double cB0 = 1432997174477./ 9575080441755.;
        const double cB1 = 5161836677717./13612068292357.;
        const double cB2 = 1720146321549./ 2090206949498.;
        const double cB3 = 3134564353537./ 4481467310338.;
        const double cB4 = 2277821191437./14882151754819.;

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
void Timeloop::exec()
{
    const int blocki = grid->ithread_block;
    const int blockj = grid->jthread_block;
    const int gridi = grid->imax/blocki + (grid->imax%blocki > 0);
    const int gridj = grid->jmax/blockj + (grid->jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, grid->kmax);
    dim3 blockGPU(blocki, blockj, 1);

    const int offs = grid->memoffset;

    if (rkorder == 3)
    {
        for (FieldMap::const_iterator it = fields->at.begin(); it!=fields->at.end(); ++it)
        {
            if (substep == 0)
                rk3_g<0><<<gridGPU, blockGPU>>>(
                    &fields->ap[it->first]->data_g[offs], &it->second->data_g[offs], dt,
                    grid->icellsp, grid->ijcellsp,
                    grid->istart,  grid->jstart, grid->kstart,
                    grid->iend,    grid->jend,   grid->kend);
            else if (substep == 1)
                rk3_g<1><<<gridGPU, blockGPU>>>(
                    &fields->ap[it->first]->data_g[offs], &it->second->data_g[offs], dt,
                    grid->icellsp, grid->ijcellsp,
                    grid->istart,  grid->jstart, grid->kstart,
                    grid->iend,    grid->jend,   grid->kend);
            else if (substep == 2)
                rk3_g<2><<<gridGPU, blockGPU>>>(
                    &fields->ap[it->first]->data_g[offs], &it->second->data_g[offs], dt,
                    grid->icellsp, grid->ijcellsp,
                    grid->istart,  grid->jstart, grid->kstart,
                    grid->iend,    grid->jend,   grid->kend);
        }

        substep = (substep+1) % 3;

        /*
           rk3_kernel<<<gridGPU, blockGPU>>>(a, at, dt,
           substep, grid->icells, grid->ijcells,
           grid->istart, grid->jstart, grid->kstart,
           grid->iend, grid->jend, grid->kend);
         */
    }

    else if (rkorder == 4)
    {
        for (FieldMap::const_iterator it = fields->at.begin(); it!=fields->at.end(); ++it)
        {
            if (substep==0)
                rk4_g<0><<<gridGPU, blockGPU>>>(
                    &fields->ap[it->first]->data_g[offs], &it->second->data_g[offs], dt,
                    grid->icellsp, grid->ijcellsp,
                    grid->istart,  grid->jstart, grid->kstart,
                    grid->iend,    grid->jend,   grid->kend);
            else if (substep==1)
                rk4_g<1><<<gridGPU, blockGPU>>>(
                    &fields->ap[it->first]->data_g[offs], &it->second->data_g[offs], dt,
                    grid->icellsp, grid->ijcellsp,
                    grid->istart,  grid->jstart, grid->kstart,
                    grid->iend,    grid->jend,   grid->kend);
            else if (substep==2)
                rk4_g<2><<<gridGPU, blockGPU>>>(
                    &fields->ap[it->first]->data_g[offs], &it->second->data_g[offs], dt,
                    grid->icellsp, grid->ijcellsp,
                    grid->istart,  grid->jstart, grid->kstart,
                    grid->iend,    grid->jend,   grid->kend);
            else if (substep==3)
                rk4_g<3><<<gridGPU, blockGPU>>>(
                    &fields->ap[it->first]->data_g[offs], &it->second->data_g[offs], dt,
                    grid->icellsp, grid->ijcellsp,
                    grid->istart,  grid->jstart, grid->kstart,
                    grid->iend,    grid->jend,   grid->kend);
            else if (substep==4)
                rk4_g<4><<<gridGPU, blockGPU>>>(
                    &fields->ap[it->first]->data_g[offs], &it->second->data_g[offs], dt,
                    grid->icellsp, grid->ijcellsp,
                    grid->istart,  grid->jstart, grid->kstart,
                    grid->iend,    grid->jend,   grid->kend);
        }

        substep = (substep+1) % 5;

        /*
           rk4_kernel<<<gridGPU, blockGPU>>>(a, at, dt,
           substep, grid->icells, grid->ijcells,
           grid->istart, grid->jstart, grid->kstart,
           grid->iend, grid->jend, grid->kend);
         */
    }

    cuda_check_error();
}
#endif
