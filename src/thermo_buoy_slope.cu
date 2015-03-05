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

#include <cstdio>
#include "grid.h"
#include "fields.h"
#include "thermo_buoy_slope.h"
#include "master.h"
#include "finite_difference.h"
#include "tools.h"

namespace
{
    __global__ 
    void calc_buoyancy_tend_u_4th_g(double* const __restrict__ ut, const double* const __restrict__ b,
                                    const double sinalpha,
                                    const int istart, const int jstart, const int kstart,
                                    const int iend,   const int jend,   const int kend,
                                    const int jj,     const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
        const int k = blockIdx.z + kstart;

        const int ii1 = 1;
        const int ii2 = 2;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            ut[ijk] += sinalpha * (fd::o4::ci0*b[ijk-ii2] + fd::o4::ci1*b[ijk-ii1] + fd::o4::ci2*b[ijk] + fd::o4::ci3*b[ijk+ii1]);
        }
    }

    __global__ 
    void calc_buoyancy_tend_w_4th_g(double* __restrict__ wt, const double* const __restrict__ b,
                                    const double cosalpha,
                                    const int istart, const int jstart, const int kstart,
                                    const int iend, const int jend, const int kend,
                                    const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
        const int k = blockIdx.z + kstart;

        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            wt[ijk] += cosalpha * (fd::o4::ci0*b[ijk-kk2] + fd::o4::ci1*b[ijk-kk1] + fd::o4::ci2*b[ijk] + fd::o4::ci3*b[ijk+kk1]);
        }
    }

    __global__ 
    void calc_buoyancy_tend_b_4th_g(double* const __restrict__ bt,
                                    const double* const __restrict__ u, const double* const __restrict__ w,
                                    const double utrans, const double n2, const double sinalpha, const double cosalpha,
                                    const int istart, const int jstart, const int kstart,
                                    const int iend, const int jend, const int kend,
                                    const int jj, const int kk)
    {
        const int i = blockIdx.x*blockDim.x + threadIdx.x + istart; 
        const int j = blockIdx.y*blockDim.y + threadIdx.y + jstart; 
        const int k = blockIdx.z + kstart;

        const int ii1 = 1;
        const int ii2 = 2;

        const int kk1 = 1*kk;
        const int kk2 = 2*kk;

        if (i < iend && j < jend && k < kend)
        {
            const int ijk = i + j*jj + k*kk;
            bt[ijk] -= n2 * ( sinalpha * ( (fd::o4::ci0*u[ijk-ii1] + fd::o4::ci1*u[ijk] + fd::o4::ci2*u[ijk+ii1] + fd::o4::ci3*u[ijk+ii2]) + utrans )
                            + cosalpha * (  fd::o4::ci0*w[ijk-kk1] + fd::o4::ci1*w[ijk] + fd::o4::ci2*w[ijk+kk1] + fd::o4::ci3*w[ijk+kk2]) );
        }
    }
} // End namespace.

#ifdef USECUDA
void Thermo_buoy_slope::exec()
{
    const int blocki = 128;
    const int blockj = 2;
    const int gridi  = grid->imax/blocki + (grid->imax%blocki > 0);
    const int gridj  = grid->jmax/blockj + (grid->jmax%blockj > 0);

    dim3 gridGPU (gridi, gridj, grid->kmax-1);
    dim3 blockGPU(blocki, blockj, 1);

    const int offs = grid->memoffset;

    if (grid->swspatialorder== "2")
    {
        master->print_error("Second order accuracy not implemented for slope flow thermodynamics\n");
        throw 1;
    }
    else if (grid->swspatialorder== "4")
    {
        const double sinalpha = std::sin(this->alpha);
        const double cosalpha = std::cos(this->alpha);

        calc_buoyancy_tend_u_4th_g<<<gridGPU, blockGPU>>>(
            &fields->ut->data_g[offs], &fields->sp["b"]->data_g[offs],
            sinalpha,
            grid->istart,  grid->jstart, grid->kstart,
            grid->iend,    grid->jend,   grid->kend,
            grid->icellsp, grid->ijcellsp);
        cuda_check_error(); 

        calc_buoyancy_tend_w_4th_g<<<gridGPU, blockGPU>>>(
            &fields->wt->data_g[offs], &fields->sp["b"]->data_g[offs],
            cosalpha,
            grid->istart,  grid->jstart, grid->kstart+1,
            grid->iend,    grid->jend,   grid->kend,
            grid->icellsp, grid->ijcellsp);
        cuda_check_error(); 

        calc_buoyancy_tend_b_4th_g<<<gridGPU, blockGPU>>>(
            &fields->st["b"]->data_g[offs],
            &fields->u->data_g[offs], &fields->w->data_g[offs],
            grid->utrans, n2, sinalpha, cosalpha,
            grid->istart,  grid->jstart, grid->kstart,
            grid->iend,    grid->jend,   grid->kend,
            grid->icellsp, grid->ijcellsp);
        cuda_check_error();
    }
}
#endif
