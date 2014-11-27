/*
 * MicroHH
 * Copyright (c) 2011-2014 Chiel van Heerwaarden
 * Copyright (c) 2011-2014 Thijs Heus
 * Copyright (c)      2014 Bart van Stratum
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
#include <cufft.h>
#include "grid.h"
#include "fields.h"
#include "pres.h"
#include "pres_2.h"
#include "model.h"
#include "tools.h"


namespace Pres_g
{
  const int TILE_DIM = 16; // Size of shared memory array used for transpose

  __global__ void transpose(double *fieldOut, const double *fieldIn, const int itot, const int jtot, const int ktot)
  {
    __shared__ double tile[TILE_DIM][TILE_DIM+1];
    
    int i,j,k,ijk;
   
    // Index in fieldIn 
    i = blockIdx.x * TILE_DIM + threadIdx.x;
    j = blockIdx.y * TILE_DIM + threadIdx.y;
    k = blockIdx.z;
    ijk = i + j*itot + k*itot*jtot;
  
    // Read to shared memory
    if(i < itot && j < jtot)
      tile[threadIdx.y][threadIdx.x] = fieldIn[ijk];
   
    __syncthreads();
    
    // Transposed index
    i = blockIdx.y * TILE_DIM + threadIdx.x;
    j = blockIdx.x * TILE_DIM + threadIdx.y;
    ijk = i + j*jtot + k*itot*jtot;
  
    // Write transposed field back from shared to global memory 
    if(i < jtot && j < itot) 
      fieldOut[ijk] = tile[threadIdx.x][threadIdx.y];
  }

  __global__ void complex_double_x(cufftDoubleComplex * __restrict__ cdata, double * __restrict__ ddata, const unsigned int itot, const unsigned int jtot, unsigned int kk, unsigned int kki, bool forward)
  {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;
    int k = blockIdx.z;

    int ij   = i + j*itot + k*kk;         // index real part in ddata
    int ij2  = (itot-i) + j*itot + k*kk;  // index complex part in ddata
    int imax = itot/2+1;
    int ijc  = i + j*imax + k*kki;        // index in cdata

    if((j < jtot) && (i < imax))
    {
      if(forward) // complex -> double
      {
        ddata[ij]  = cdata[ijc].x;
        if(i>0 && i<imax-1)
          ddata[ij2] = cdata[ijc].y;
      }
      else // double -> complex
      {
        cdata[ijc].x = ddata[ij];
        if(i>0 && i<imax-1)
          cdata[ijc].y = ddata[ij2];
      }
    }
  }

  __global__ void complex_double_y(cufftDoubleComplex * __restrict__ cdata, double * __restrict__ ddata, const unsigned int itot, const unsigned int jtot, unsigned int kk, unsigned int kkj, bool forward)
  {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;
    int k = blockIdx.z;

    int ij   = i + j*itot + k*kk;        // index real part in ddata
    int ij2  = i + (jtot-j)*itot + k*kk;    // index complex part in ddata
    int jmax = jtot/2+1;
    int ijc  = i + j*itot + k*kkj;

    if((i < itot) && (j < jmax))
    {
      if(forward) // complex -> double
      {
        ddata[ij] = cdata[ijc].x;
        if(j>0 && j<jmax-1)
          ddata[ij2] = cdata[ijc].y;
      }
      else // double -> complex
      {
        cdata[ijc].x = ddata[ij];
        if(j>0 && j<jmax-1)
          cdata[ijc].y = ddata[ij2];
      }
    }
  }

   __global__ void normalize(double * const __restrict__ data, const int itot, const int jtot, const int ktot, const double in)
  {
    int i = blockIdx.x*blockDim.x + threadIdx.x;
    int j = blockIdx.y*blockDim.y + threadIdx.y;
    int k = blockIdx.z;

    int ijk = i + j*itot + k*itot*jtot;
    if((i < itot) && (j < jtot) && (k < ktot))
      data[ijk] = data[ijk] * in;
  }
}

#ifdef USECUDA
void Pres::fftForward(double * __restrict__ p, double * __restrict__ tmp1, double * __restrict__ tmp2)
{
  const int blocki = grid->iThreadBlock;
  const int blockj = grid->jThreadBlock;
  int gridi = grid->imax/blocki + (grid->imax%blocki > 0);
  int gridj = grid->jmax/blockj + (grid->jmax%blockj > 0);

  // 3D grid
  dim3 gridGPU (gridi,  gridj,  grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  // Square grid for transposes 
  const int gridiT = grid->imax/Pres_g::TILE_DIM + (grid->imax%Pres_g::TILE_DIM > 0);
  const int gridjT = grid->jmax/Pres_g::TILE_DIM + (grid->jmax%Pres_g::TILE_DIM > 0);
  dim3 gridGPUTf(gridiT, gridjT, grid->ktot);
  dim3 gridGPUTb(gridjT, gridiT, grid->ktot);
  dim3 blockGPUT(Pres_g::TILE_DIM, Pres_g::TILE_DIM, 1);

  // Transposed grid
  gridi = grid->jmax/blocki + (grid->jmax%blocki > 0);
  gridj = grid->imax/blockj + (grid->imax%blockj > 0);
  dim3 gridGPUji (gridi,  gridj,  grid->kmax);

  const int kk = grid->itot*grid->jtot;
  const int kki = (grid->itot/2+1)*grid->jtot;
  const int kkj = (grid->jtot/2+1)*grid->itot;

  // Forward FFT in the x-direction.
  if(FFTPerSlice) // Batched FFT per horizontal slice
  {
    for (int k=0; k<grid->ktot; ++k)
    {
      int ijk  = k*kk;
      int ijk2 = 2*k*kki;
      cufftExecD2Z(iplanf, (cufftDoubleReal*)&p[ijk], (cufftDoubleComplex*)&tmp1[ijk2]);
    }
  }
  else // Single batched FFT over entire 3D field
  {
    
    cufftExecD2Z(iplanf, (cufftDoubleReal*)p, (cufftDoubleComplex*)tmp1);
    cudaThreadSynchronize();
  }

  // Transform complex to double output. Allows for creating parallel cuda version at a later stage
  Pres_g::complex_double_x<<<gridGPU,blockGPU>>>((cufftDoubleComplex*)tmp1, p, grid->itot, grid->jtot, kk, kki,  true);
  cudaCheckError();

  // Forward FFT in the y-direction.
  if(grid->jtot > 1)
  {
    if(FFTPerSlice) // Batched FFT per horizontal slice
    {
      for (int k=0; k<grid->ktot; ++k)
      {
        int ijk  = k*kk;
        int ijk2 = 2*k*kkj;
        cufftExecD2Z(jplanf, (cufftDoubleReal*)&p[ijk], (cufftDoubleComplex*)&tmp1[ijk2]);
      }

      cudaThreadSynchronize();
      cudaCheckError();

      Pres_g::complex_double_y<<<gridGPU,blockGPU>>>((cufftDoubleComplex*)tmp1, p, grid->itot, grid->jtot, kk, kkj, true);
      cudaCheckError();
    }
    else // Single batched FFT over entire 3D field. Y-direction FFT requires transpose of field
    {
      Pres_g::transpose<<<gridGPUTf, blockGPUT>>>(tmp2, p, grid->itot, grid->jtot, grid->ktot); 
      cudaCheckError();

      cufftExecD2Z(jplanf, (cufftDoubleReal*)tmp2, (cufftDoubleComplex*)tmp1);
      cudaThreadSynchronize();

      Pres_g::complex_double_x<<<gridGPUji,blockGPU>>>((cufftDoubleComplex*)tmp1, p, grid->jtot, grid->itot, kk, kkj,  true);
      cudaCheckError();

      Pres_g::transpose<<<gridGPUTb, blockGPUT>>>(tmp1, p, grid->jtot, grid->itot, grid->ktot); 
      cudaSafeCall(cudaMemcpy(p, tmp1, grid->ncellsp*sizeof(double), cudaMemcpyDeviceToDevice));
      cudaCheckError();
    }
  }
}

void Pres::fftBackward(double * __restrict__ p, double * __restrict__ tmp1, double * __restrict__ tmp2)
{
  const int blocki = grid->iThreadBlock;
  const int blockj = grid->jThreadBlock;
  int gridi = grid->imax/blocki + (grid->imax%blocki > 0);
  int gridj = grid->jmax/blockj + (grid->jmax%blockj > 0);

  // 3D grid
  dim3 gridGPU (gridi,  gridj,  grid->kmax);
  dim3 blockGPU(blocki, blockj, 1);

  // Square grid for transposes 
  const int gridiT = grid->imax/Pres_g::TILE_DIM + (grid->imax%Pres_g::TILE_DIM > 0);
  const int gridjT = grid->jmax/Pres_g::TILE_DIM + (grid->jmax%Pres_g::TILE_DIM > 0);
  dim3 gridGPUTf(gridiT, gridjT, grid->ktot); 
  dim3 gridGPUTb(gridjT, gridiT, grid->ktot); 
  dim3 blockGPUT(Pres_g::TILE_DIM, Pres_g::TILE_DIM, 1);

  // Transposed grid
  gridi = grid->jmax/blocki + (grid->jmax%blocki > 0);
  gridj = grid->imax/blockj + (grid->imax%blockj > 0);
  dim3 gridGPUji (gridi,  gridj,  grid->kmax);

  const int kk = grid->itot*grid->jtot;
  const int kki = (grid->itot/2+1)*grid->jtot;
  const int kkj = (grid->jtot/2+1)*grid->itot;

  // Backward FFT in the y-direction.
  if(grid->jtot > 1)
  {
    if(FFTPerSlice) // Batched FFT per horizontal slice
    {
      Pres_g::complex_double_y<<<gridGPU,blockGPU>>>((cufftDoubleComplex*)tmp1, p, grid->itot, grid->jtot, kk, kkj, false);
      cudaCheckError();
      for (int k=0; k<grid->ktot; ++k)
      {
        int ijk = k*kk;
        int ijk2 = 2*k*kkj;
        cufftExecZ2D(jplanb, (cufftDoubleComplex*)&tmp1[ijk2], (cufftDoubleReal*)&p[ijk]);
      }
      cudaThreadSynchronize();
      cudaCheckError();
    }
    else // Single batched FFT over entire 3D field. Y-direction FFT requires transpose of field
    {
      Pres_g::transpose<<<gridGPUTf, blockGPUT>>>(tmp2, p, grid->itot, grid->jtot, grid->ktot); 
      cudaCheckError();

      Pres_g::complex_double_x<<<gridGPUji,blockGPU>>>((cufftDoubleComplex*)tmp1, tmp2, grid->jtot, grid->itot, kk, kkj, false);
      cudaCheckError();

      cufftExecZ2D(jplanb, (cufftDoubleComplex*)tmp1, (cufftDoubleReal*)p);
      cudaThreadSynchronize();
      cudaCheckError();

      Pres_g::transpose<<<gridGPUTb, blockGPUT>>>(tmp1, p, grid->jtot, grid->itot, grid->ktot); 
      cudaCheckError();
      cudaSafeCall(cudaMemcpy(p, tmp1, grid->ncellsp*sizeof(double), cudaMemcpyDeviceToDevice));
      cudaCheckError();
    }
 }

  // Backward FFT in the x-direction
  Pres_g::complex_double_x<<<gridGPU,blockGPU>>>((cufftDoubleComplex*)tmp1, p, grid->itot, grid->jtot, kk, kki,  false);
  cudaCheckError();

  if(FFTPerSlice) // Batched FFT per horizontal slice
  {
    for (int k=0; k<grid->ktot; ++k)
    {
      int ijk = k*kk;
      int ijk2 = 2*k*kki;
      cufftExecZ2D(iplanb, (cufftDoubleComplex*)&tmp1[ijk2], (cufftDoubleReal*)&p[ijk]);
    }
    cudaThreadSynchronize();
    cudaCheckError();
  }
  else // Batch FFT over entire domain
  {
    cufftExecZ2D(iplanb, (cufftDoubleComplex*)tmp1, (cufftDoubleReal*)p);
    cudaThreadSynchronize();
    cudaCheckError();
  }

  // Normalize output
  Pres_g::normalize<<<gridGPU,blockGPU>>>(p, grid->itot, grid->jtot, grid->ktot, 1./(grid->itot*grid->jtot));
  cudaCheckError();
}
#endif
