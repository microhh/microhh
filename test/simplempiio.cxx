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

int checkerror(int);

int main()
{
  // MPI vars
  int mpiid;
  int nprocs;
  MPI_Comm commxy;
  MPI_Comm commx;
  MPI_Comm commy;
  int mpicoordx;
  int mpicoordy;

  // grid vars
  int npx  = 16;
  int npy  = 32;
  int kmax = 1024;
  int itot = 1024;
  int jtot = 1024;

  int imax = itot / npx;
  int jmax = jtot / npy;

  int count = imax*jmax*kmax;

  double *u;

  // CREATE THE GRID
  u = new double[count];

  // START THE MPI
  int n;

  // initialize the MPI
  n = MPI_Init(NULL, NULL);
  if(checkerror(n))
    return 1;

  // get the rank of the current process
  n = MPI_Comm_rank(MPI_COMM_WORLD, &mpiid);
  if(checkerror(n))
    return 1;

  // get the total number of processors
  n = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  if(checkerror(n))
    return 1;

  if(mpiid == 0) std::printf("Starting run on %d processes\n", nprocs);

  // CREATE 2D GRID
  int dims    [2] = {npy, npx};
  int periodic[2] = {true, true};

  // define the dimensions of the 2-D grid layout
  n = MPI_Dims_create(nprocs, 2, dims);
  if(checkerror(n))
    return 1;

  // create a 2-D grid communicator that is optimized for grid to grid transfer
  n = MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, true , &commxy);
  if(checkerror(n))
    return 1;

  // retrieve the x- and y-coordinates in the 2-D grid for each process
  int mpicoords[2];
  n = MPI_Cart_coords(commxy, mpiid, 2, mpicoords);
  if(checkerror(n))
    return 1;

  mpicoordx = mpicoords[1];
  mpicoordy = mpicoords[0];

  int dimx[2] = {false, true };
  int dimy[2] = {true , false};

  n = MPI_Cart_sub(commxy, dimx, &commx);
  if(checkerror(n))
    return 1;
  n = MPI_Cart_sub(commxy, dimy, &commy);
  if(checkerror(n))
    return 1;

  // CREATE THE SUBARRAY
  MPI_Datatype subarray;
  int totsize [3] = {kmax, jtot, itot};
  int subsize [3] = {kmax, jmax, imax};
  int substart[3] = {0, mpicoordy*jmax, mpicoordx*imax};
  MPI_Type_create_subarray(3, totsize, subsize, substart, MPI_ORDER_C, MPI_DOUBLE, &subarray);
  MPI_Type_commit(&subarray);

  // SET THE VALUE OF THE GRID EQUAL TO THE PROCESS ID FOR CHECKING
  if(mpiid == 0) std::printf("Setting the value of the array\n");
  for(int i=0; i<count; i++)
    u[i] = (double)mpiid;

  // WRITE THE FULL GRID USING MPI-IO
  if(mpiid == 0) std::printf("Write the full array to disk\n");
  char filename[] = "u.dump";
  MPI_File fh;
  n = MPI_File_open(commxy, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, MPI_INFO_NULL, &fh);
  if(checkerror(n))
    return 1;

  // select noncontiguous part of 3d array to store the selected data
  MPI_Offset fileoff = 0; // the offset within the file (header size)
  char name[] = "native";

  n = MPI_File_set_view(fh, fileoff, MPI_DOUBLE, subarray, name, MPI_INFO_NULL);
  if(checkerror(n))
    return 1;

  n = MPI_File_write_all(fh, u, count, MPI_DOUBLE, MPI_STATUS_IGNORE);
  if(checkerror(n))
    return 1;

  MPI_File_close(&fh);
  if(checkerror(n))
    return 1;

  if(mpiid == 0) std::printf("Writing done\n");

  // CLOSE THE MPI
  MPI_Finalize();
  if(mpiid == 0) std::printf("Finished run on %d processes\n", nprocs);

  delete[] u;

  return 0;
}

int checkerror(int n)
{
  char errbuffer[MPI_MAX_ERROR_STRING];
  int errlen;

  if(n != MPI_SUCCESS)
  {
    MPI_Error_string(n, errbuffer, &errlen);
    std::printf("ERROR MPI %s\n", errbuffer);
    return 1;
  }

  return 0;
}

