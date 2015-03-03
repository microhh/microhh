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
#include <cmath>
#include "grid.h"
#include "fields.h"
#include "mpiinterface.h"
#include "mpicheck.h"

cmpicheck::cmpicheck(cgrid *gridin, cmpi *mpiin)
{
  std::printf("Creating instance of object mpicheck\n");
  grid   = gridin;
  mpi    = mpiin;
}

cmpicheck::~cmpicheck()
{
  std::printf("Destroying instance of object mpicheck\n");
}

int cmpicheck::checkLayout()
{
  
  std::printf("MPI id, mpicoordx, mpicoordy, neast, nwest, nnorth, nsouth, nprocs: %2d, %2d, %2d, %2d, %2d, %2d, %2d, %2d\n",
    mpi->mpiid, mpi->mpicoordx, mpi->mpicoordy, mpi->neast, mpi->nwest, mpi->nnorth, mpi->nsouth, mpi->nprocs);

  return 0;
}

int cmpicheck::create()
{
  s     = new cfield3d(grid, mpi, "s");
  temp1 = new cfield3d(grid, mpi, "temp1");
  temp2 = new cfield3d(grid, mpi, "temp2");

  s->init();
  temp1->init();
  temp2->init();

  int n;

  for(n=0; n<grid->ncells; n++)
    s->data[n] = (double)mpi->mpiid;
  
  return 0;
}

int cmpicheck::checkBoundary()
{
  grid->boundary_cyclic(s->data);

  int i,j,k;
  int ijk,jj,kk;

  jj = grid->icells;
  kk = grid->icells*grid->jcells;

  k = grid->kstart;
  j = grid->jstart;
  for(i=0; i<grid->icells; i++)
  {
    ijk = i + j*jj + k*kk;
    std::printf("MPI i-line id %d, s(%d,%d,%d) = %4.0f\n", mpi->mpiid, i, j, k, s->data[ijk]);
  }

  k = grid->kstart;
  i = grid->istart;
  for(j=0; j<grid->jcells; j++)
  {
    ijk = i + j*jj + k*kk;
    std::printf("MPI j-line id %d, s(%d,%d,%d) = %4.0f\n", mpi->mpiid, i, j, k, s->data[ijk]);
  }

  return 0;
}

int cmpicheck::checkTranspose()
{
  int ijk,ijkw;
  int jj,kk,jjw,kkw;
  int igc,jgc,kgc;

  jj  = grid->icells;
  kk  = grid->icells*grid->jcells;
  jjw = grid->imax;
  kkw = grid->imax*grid->jmax;
  igc = grid->igc;
  jgc = grid->jgc;
  kgc = grid->kgc;

  for(int k=grid->kstart; k<grid->kend; k++)
    for(int j=grid->jstart; j<grid->jend; j++)
      for(int i=grid->istart; i<grid->iend; i++)
      {
        ijk  = i + j*jj + k*kk;
        ijkw = (i-igc) + (j-jgc)*jjw + (k-kgc)*kkw;

        temp1->data[ijkw] = mpi->mpiid; //k - kgc + 100*mpi->mpiid;
      }

  grid->transposezx(temp2->data, temp1->data);

  jj = grid->itot;
  kk = grid->itot*grid->jmax;

  for(int k=0; k<grid->kblock; k++)
    for(int j=0; j<grid->jmax; j++)
      for(int i=0; i<grid->itot; i++)
      {
        ijk = i + j*jj + k*kk;
        std::printf("MPI transzx id %d, (%d,%d,%d) = %4.0f\n", mpi->mpiid, i, j, k, temp2->data[ijk]);
      }

  grid->transposexz(temp1->data, temp2->data);

  jj = grid->imax;
  kk = grid->imax*grid->jmax;

  for(int k=0; k<grid->kmax; k++)
    for(int j=0; j<grid->jmax; j++)
      for(int i=0; i<grid->imax; i++)
      {
        ijk = i + j*jj + k*kk;
        std::printf("MPI transxz id %d, (%d,%d,%d) = %4.0f\n", mpi->mpiid, i, j, k, temp1->data[ijk]);
      }

  grid->transposezx(temp2->data, temp1->data);
  grid->transposexy(temp1->data, temp2->data);

  jj = grid->iblock;
  kk = grid->iblock*grid->jtot;

  for(int k=0; k<grid->kblock; k++)
    for(int j=0; j<grid->jtot; j++)
      for(int i=0; i<grid->iblock; i++)
      {
        ijk = i + j*jj + k*kk;
        std::printf("MPI transxy id %d, (%d,%d,%d) = %4.0f\n", mpi->mpiid, i, j, k, temp1->data[ijk]);
      }

  grid->transposeyx(temp2->data, temp1->data);

  jj = grid->itot;
  kk = grid->itot*grid->jmax;

  for(int k=0; k<grid->kblock; k++)
    for(int j=0; j<grid->jmax; j++)
      for(int i=0; i<grid->itot; i++)
      {
        ijk = i + j*jj + k*kk;
        std::printf("MPI transyx id %d, (%d,%d,%d) = %4.0f\n", mpi->mpiid, i, j, k, temp2->data[ijk]);
      }

  grid->transposexz(temp1->data, temp2->data);

  jj = grid->imax;
  kk = grid->imax*grid->jmax;

  for(int k=0; k<grid->kmax; k++)
    for(int j=0; j<grid->jmax; j++)
      for(int i=0; i<grid->imax; i++)
      {
        ijk = i + j*jj + k*kk;
        std::printf("MPI transxz id %d, (%d,%d,%d) = %4.0f\n", mpi->mpiid, i, j, k, temp1->data[ijk]);
      }

  grid->transposezx(temp2->data, temp1->data);
  grid->transposexy(temp1->data, temp2->data);
  grid->transposeyz(temp2->data, temp1->data);

  jj = grid->iblock;
  kk = grid->iblock*grid->jblock;

  for(int k=0; k<grid->kmax; k++)
    for(int j=0; j<grid->jblock; j++)
      for(int i=0; i<grid->iblock; i++)
      {
        ijk = i + j*jj + k*kk;
        std::printf("MPI transyz id %d, (%d,%d,%d) = %4.0f\n", mpi->mpiid, i, j, k, temp2->data[ijk]);
      }

  jj = grid->imax;
  kk = grid->imax*grid->jmax;

  grid->transposezy(temp1->data, temp2->data);
  grid->transposeyx(temp2->data, temp1->data);
  grid->transposexz(temp1->data, temp2->data);

  for(int k=0; k<grid->kmax; k++)
    for(int j=0; j<grid->jmax; j++)
      for(int i=0; i<grid->imax; i++)
      {
        ijk = i + j*jj + k*kk;
        std::printf("MPI transzy id %d, (%d,%d,%d) = %4.0f\n", mpi->mpiid, i, j, k, temp1->data[ijk]);
      }

  return 0;
}

