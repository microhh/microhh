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

#include "field3d.h"
#include "grid.h"
#include "master.h"
#include "tools.h"

#ifdef USECUDA
Field3d::~Field3d()
{
    cuda_safe_call(cudaFreeHost(data));
    cuda_safe_call(cudaFreeHost(databot));
    cuda_safe_call(cudaFreeHost(datatop));
    cuda_safe_call(cudaFreeHost(datagradbot));
    cuda_safe_call(cudaFreeHost(datagradtop));
    cuda_safe_call(cudaFreeHost(datafluxbot));
    cuda_safe_call(cudaFreeHost(datafluxtop));
    cuda_safe_call(cudaFreeHost(datamean));
}

int Field3d::init()
{
    const int ijksize = grid->ncells *sizeof(double);
    const int ijsize  = grid->ijcells*sizeof(double);
    const int ksize   = grid->kcells *sizeof(double);

    // Allocate the 3d field.
    cuda_safe_call(cudaMallocHost(&data, ijksize));

    // Allocate the boundary cells.
    cuda_safe_call(cudaMallocHost(&databot, ijsize));
    cuda_safe_call(cudaMallocHost(&datatop, ijsize));
    cuda_safe_call(cudaMallocHost(&datagradbot, ijsize));
    cuda_safe_call(cudaMallocHost(&datagradtop, ijsize));
    cuda_safe_call(cudaMallocHost(&datafluxbot, ijsize));
    cuda_safe_call(cudaMallocHost(&datafluxtop, ijsize));
    cuda_safe_call(cudaMallocHost(&datamean, ksize));

    // Set all values to zero
    for (int n=0; n<grid->ncells; n++)
        data[n] = 0.;

    for (int n=0; n<grid->kcells; n++)
        datamean[n] = 0.;

    for (int n=0; n<grid->icells*grid->jcells; n++)
    {
        databot    [n] = 0.;
        datatop    [n] = 0.;
        datagradbot[n] = 0.;
        datagradtop[n] = 0.;
        datafluxbot[n] = 0.;
        datafluxtop[n] = 0.;
    }

    return 0;
}
#endif

void Field3d::init_device()
{
    const int nmemsize   = grid->ncellsp*sizeof(double);
    const int nmemsize1d = grid->kcells *sizeof(double);
    const int nmemsize2d = (grid->ijcellsp+grid->memoffset)*sizeof(double);

    cuda_safe_call(cudaMalloc(&data_g,        nmemsize  ));
    cuda_safe_call(cudaMalloc(&databot_g,     nmemsize2d));
    cuda_safe_call(cudaMalloc(&datatop_g,     nmemsize2d));
    cuda_safe_call(cudaMalloc(&datagradbot_g, nmemsize2d));
    cuda_safe_call(cudaMalloc(&datagradtop_g, nmemsize2d));
    cuda_safe_call(cudaMalloc(&datafluxbot_g, nmemsize2d));
    cuda_safe_call(cudaMalloc(&datafluxtop_g, nmemsize2d));
    cuda_safe_call(cudaMalloc(&datamean_g,    nmemsize1d));
}

void Field3d::clear_device()
{
    cuda_safe_call(cudaFree(data_g));
    cuda_safe_call(cudaFree(databot_g));
    cuda_safe_call(cudaFree(datatop_g));
    cuda_safe_call(cudaFree(datagradbot_g));
    cuda_safe_call(cudaFree(datagradtop_g));
    cuda_safe_call(cudaFree(datafluxbot_g));
    cuda_safe_call(cudaFree(datafluxtop_g));
    cuda_safe_call(cudaFree(datamean_g));
}
