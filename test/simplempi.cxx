/*
 * MicroHH
 * Copyright (c) 2011-2013 Chiel van Heerwaarden
 * Copyright (c) 2011-2013 Thijs Heus
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

#include <mpi.h>
#include <cstdio>

int main()
{
  // hard set processors
  const int npx = 8;
  const int npy = 5;

  // variables from header
  MPI_Comm commxy;
  MPI_Comm commx;
  MPI_Comm commy;

  int nprocs;
  int mpiid;
  int mpicoordx;
  int mpicoordy;

  // initialization function
  char err_buffer[MPI_MAX_ERROR_STRING];
  int  n, resultlen;

  n = MPI_Init(NULL, NULL);
  if(n != MPI_SUCCESS)
  {
    MPI_Error_string(n, err_buffer, &resultlen);
    std::printf("MPI ERROR: %s\n", err_buffer);
    return 1;
  }

  n = MPI_Comm_rank(MPI_COMM_WORLD, &mpiid);
  if(n != MPI_SUCCESS)
  {
    MPI_Error_string(n, err_buffer, &resultlen);
    std::printf("MPI ERROR: %s\n", err_buffer);
    return 1;
  }

  n = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if(n != MPI_SUCCESS)
  {
    MPI_Error_string(n, err_buffer, &resultlen);
    std::printf("MPI ERROR: %s\n", err_buffer);
    return 1;
  }

  int dims    [2] = {npy, npx};
  int periodic[2] = {true, true};

  n = MPI_Dims_create(nprocs, 2, dims);
  if(n != MPI_SUCCESS)
  {
    MPI_Error_string(n, err_buffer, &resultlen);
    std::printf("MPI ERROR: %s\n", err_buffer);
    return 1;
  }

  n = MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, true, &commxy);
  if(n != MPI_SUCCESS)
  {
    MPI_Error_string(n, err_buffer, &resultlen);
    std::printf("MPI ERROR: %s\n", err_buffer);
    return 1;
  }

  n = MPI_Comm_rank(commxy, &mpiid);
  if(n != MPI_SUCCESS)
  {
    MPI_Error_string(n, err_buffer, &resultlen);
    std::printf("MPI ERROR: %s\n", err_buffer);
    return 1;
  }

  int mpicoords[2];
  n = MPI_Cart_coords(commxy, mpiid, 2, mpicoords);
  if(n != MPI_SUCCESS)
  {
    MPI_Error_string(n, err_buffer, &resultlen);
    std::printf("MPI ERROR: %s\n", err_buffer);
    return 1;
  }

  mpicoordx = mpicoords[1];
  mpicoordy = mpicoords[0];

  int dimx[2] = {false, true };
  int dimy[2] = {true , false};

  n = MPI_Cart_sub(commxy, dimx, &commx);
  if(n != MPI_SUCCESS)
  {
    MPI_Error_string(n, err_buffer, &resultlen);
    std::printf("MPI ERROR: %s\n", err_buffer);
    return 1;
  }

  n = MPI_Cart_sub(commxy, dimy, &commy);
  if(n != MPI_SUCCESS)
  {
    MPI_Error_string(n, err_buffer, &resultlen);
    std::printf("MPI ERROR: %s\n", err_buffer);
    return 1;
  }

  MPI_Finalize();

  return 0;
}

