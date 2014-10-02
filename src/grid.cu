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

#include "grid.h"
#include "tools.h"

__global__ void grid_cyclic_x(double * const __restrict__ data,
                              const int icells, const int jcells, const int kcells,
                              const int icellsp,
                              const int istart, const int jstart,
                              const int iend,   const int jend, 
                              const int igc,    const int jgc)
{
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  const int j = blockIdx.y*blockDim.y + threadIdx.y;
  const int k = blockIdx.z;

  const int jj = icellsp;
  const int kk = icellsp*jcells;

  // East-west
  if(k < kcells && j < jcells && i < igc)
  {
    const int ijk0 = i          + j*jj + k*kk;
    const int ijk1 = iend-igc+i + j*jj + k*kk;
    const int ijk2 = i+iend     + j*jj + k*kk;
    const int ijk3 = i+istart   + j*jj + k*kk;

    data[ijk0] = data[ijk1];
    data[ijk2] = data[ijk3];
  }
}

__global__ void grid_cyclic_y(double * const __restrict__ data,
                              const int icells, const int jcells, const int kcells,
                              const int icellsp,
                              const int istart, const int jstart,
                              const int iend,   const int jend, 
                              const int igc,    const int jgc)
{
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  const int j = blockIdx.y*blockDim.y + threadIdx.y;
  const int k = blockIdx.z;

  const int jj = icellsp;
  const int kk = icellsp*jcells;

  // North-south
  if(jend-jstart == 1)
  {
    if(k < kcells && j < jgc && i < icells)
    {
      const int ijkref   = i + jstart*jj   + k*kk;
      const int ijknorth = i + j*jj        + k*kk;
      const int ijksouth = i + (jend+j)*jj + k*kk;
      data[ijknorth] = data[ijkref];
      data[ijksouth] = data[ijkref];
    }
  }
  else
  {
    if(k < kcells && j < jgc && i < icells)
    {
      const int ijk0 = i + j           *jj + k*kk;
      const int ijk1 = i + (jend-jgc+j)*jj + k*kk;
      const int ijk2 = i + (j+jend  )*jj + k*kk;
      const int ijk3 = i + (j+jstart)*jj + k*kk;

      data[ijk0] = data[ijk1];
      data[ijk2] = data[ijk3];
    }
  }
}

int cgrid::prepareDevice()
{
  /* Align the interior of the grid (i.e. excluding ghost cells) with 
     the 128 byte memory blocks of the GPU's global memory */
  memoffset = 16 - igc;           // Padding at start of array 
  int padl  = 16-(int)imax%16;    // Elements left in last 128 byte block
  icellsp   = imax + padl + (padl < 2*igc) * 16;
  ijcellsp  = icellsp * jcells;  
  ncellsp   = ijcellsp * kcells + memoffset;

  const int kmemsize = kcells*sizeof(double);

  cudaSafeCall(cudaMalloc((void**)&z_g    , kmemsize));
  cudaSafeCall(cudaMalloc((void**)&zh_g   , kmemsize));
  cudaSafeCall(cudaMalloc((void**)&dz_g   , kmemsize));
  cudaSafeCall(cudaMalloc((void**)&dzh_g  , kmemsize));
  cudaSafeCall(cudaMalloc((void**)&dzi_g  , kmemsize));
  cudaSafeCall(cudaMalloc((void**)&dzhi_g , kmemsize));
  cudaSafeCall(cudaMalloc((void**)&dzi4_g , kmemsize));
  cudaSafeCall(cudaMalloc((void**)&dzhi4_g, kmemsize));

  cudaSafeCall(cudaMemcpy(z_g    , z    , kmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(zh_g   , zh   , kmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(dz_g   , dz   , kmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(dzh_g  , dzh  , kmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(dzi_g  , dzi  , kmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(dzhi_g , dzhi , kmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(dzi4_g , dzi4 , kmemsize, cudaMemcpyHostToDevice));
  cudaSafeCall(cudaMemcpy(dzhi4_g, dzhi4, kmemsize, cudaMemcpyHostToDevice));

  return 0;
}

int cgrid::clearDevice()
{
  cudaSafeCall(cudaFree(z_g    ));
  cudaSafeCall(cudaFree(zh_g   ));
  cudaSafeCall(cudaFree(dz_g   ));
  cudaSafeCall(cudaFree(dzh_g  ));
  cudaSafeCall(cudaFree(dzi_g  ));
  cudaSafeCall(cudaFree(dzhi_g ));
  cudaSafeCall(cudaFree(dzi4_g ));
  cudaSafeCall(cudaFree(dzhi4_g));

  return 0;
}

int cgrid::boundary_cyclic_g(double * data)
{
  const int blocki_x = igc;
  const int blockj_x = 256 / igc + (256%igc > 0);
  const int gridi_x  = 1;
  const int gridj_x  = jcells/blockj_x + (jcells%blockj_x > 0);

  const int blocki_y = 256 / jgc + (256%jgc > 0);
  const int blockj_y = jgc;
  const int gridi_y  = icells/blocki_y + (icells%blocki_y > 0);
  const int gridj_y  = 1;

  dim3 gridGPUx (gridi_x, gridj_x, kcells);
  dim3 blockGPUx(blocki_x, blockj_x, 1);

  dim3 gridGPUy (gridi_y, gridj_y, kcells);
  dim3 blockGPUy(blocki_y, blockj_y, 1);

  grid_cyclic_x<<<gridGPUx,blockGPUx>>>(data, icells, jcells, kcells, icellsp,
                                        istart, jstart,
                                        iend,   jend,
                                        igc,    jgc);

  grid_cyclic_y<<<gridGPUy,blockGPUy>>>(data, icells, jcells, kcells, icellsp,
                                        istart, jstart,
                                        iend,   jend,
                                        igc,    jgc);

  cudaCheckError();

  return 0;
}

int cgrid::boundary_cyclic2d_g(double * data)
{
  const int blocki_x = igc;
  const int blockj_x = 256 / igc + (256%igc > 0);
  const int gridi_x  = 1;
  const int gridj_x  = jcells/blockj_x + (jcells%blockj_x > 0);

  const int blocki_y = 256 / jgc + (256%jgc > 0);
  const int blockj_y = jgc;
  const int gridi_y  = icells/blocki_y + (icells%blocki_y > 0);
  const int gridj_y  = 1;

  dim3 gridGPUx (gridi_x, gridj_x, 1);
  dim3 blockGPUx(blocki_x, blockj_x, 1);

  dim3 gridGPUy (gridi_y, gridj_y, 1);
  dim3 blockGPUy(blocki_y, blockj_y, 1);

  grid_cyclic_x<<<gridGPUx,blockGPUx>>>(data, icells, jcells, kcells, icellsp,
                                        istart, jstart,
                                        iend,   jend,
                                        igc,    jgc);

  grid_cyclic_y<<<gridGPUy,blockGPUy>>>(data, icells, jcells, kcells, icellsp,
                                        istart, jstart,
                                        iend,   jend,
                                        igc,    jgc);

  cudaCheckError();

  return 0;
}


double cgrid::getmax_g(double *data, double *tmp)
{
  const unsigned int max = 1;
  const double scalefac = 1.;
  double maxvalue;

  // Reduce 3D field excluding ghost cells and padding to jtot*ktot values
  reduceInterior(data, tmp, itot, istart, iend, jtot, jstart, jend, ktot, kstart, icellsp, ijcellsp, max);
  // Reduce jtot*ktot to ktot values
  reduceAll     (tmp, &tmp[jtot*ktot], jtot*ktot, ktot, jtot, max, scalefac);
  // Reduce ktot values to a single value
  reduceAll     (&tmp[jtot*ktot], tmp, ktot, 1, ktot, max, scalefac);
  // Copy back result from GPU
  cudaSafeCall(cudaMemcpy(&maxvalue, &tmp[0], sizeof(double), cudaMemcpyDeviceToHost));
  
  return maxvalue;
}

double cgrid::getsum_g(double *data, double *tmp)
{
  const unsigned int sum = 0;
  const double scalefac = 1.;
  double sumvalue;

  reduceInterior(data, tmp, itot, istart, iend, jtot, jstart, jend, ktot, kstart, icellsp, ijcellsp, sum);
  reduceAll     (tmp, &tmp[jtot*ktot], jtot*ktot, ktot, jtot, sum, scalefac);
  reduceAll     (&tmp[jtot*ktot], tmp, ktot, 1, ktot, sum, scalefac);

  cudaSafeCall(cudaMemcpy(&sumvalue, &tmp[0], sizeof(double), cudaMemcpyDeviceToHost));
  
  return sumvalue;
}

int cgrid::calcmean_g(double *prof, double *data, double *tmp)
{
  const unsigned int sum = 0;
  const double scalefac = 1./(itot*jtot);

  // Reduce 3D field excluding ghost cells and padding to jtot*kcells values
  reduceInterior(data, tmp, itot, istart, iend, jtot, jstart, jend, kcells, 0, icellsp, ijcellsp, sum);
  // Reduce jtot*kcells to kcells values
  reduceAll     (tmp, prof, jtot*kcells, kcells, jtot, sum, scalefac);

  return 0;
} 

