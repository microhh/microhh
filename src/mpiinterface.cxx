#include <cstdio>
#include <mpi.h>
#include "grid.h"
#include "fields.h"
#include "defines.h"
#include "mpiinterface.h"

// cmpi::cmpi(cgrid *gridin, cfields *fieldsin)
cmpi::cmpi()
{
  std::printf("Creating instance of object mpi\n");
  // grid   = gridin;
  // fields = fieldsin;
}

cmpi::~cmpi()
{
  MPI_Finalize();
  std::printf("Destroying instance of object mpi\n");
}

int cmpi::readinifile(cinput *inputin)
{
  int n = 0;

  n += inputin->getItem(&npx, "grid", "npx");
  n += inputin->getItem(&npy, "grid", "npy");

  if(n > 0)
    return 1;
  
  return 0;
}

int cmpi::init()
{
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  int dims    [2] = {npy, npx};
  int periodic[2] = {true, true};

  MPI_Dims_create(nprocs, 2, dims);
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, true, &commxy);
  MPI_Comm_rank(commxy, &mpiid);

  int mpicoords[2];
  MPI_Cart_coords(commxy, mpiid, 2, mpicoords);

  mpicoordx = mpicoords[1];
  mpicoordy = mpicoords[0];

  int dimx[2] = {false, true };
  int dimy[2] = {true , false};

  // MPI_Cart_sub(commxy, dimx, &commx);
  // MPI_Cart_sub(commxy, dimy, &commy);

  MPI_Cart_shift(commxy, 1, 1, &nwest , &neast );
  MPI_Cart_shift(commxy, 0, 1, &nsouth, &nnorth);

  std::printf("MPI id, mpicoordx, mpicoordy, neast, nwest, nnorth, nsouth, nprocs: %2d, %2d, %2d, %2d, %2d, %2d, %2d, %2d\n", mpiid, mpicoordx, mpicoordy, neast, nwest, nnorth, nsouth, nprocs);

  return 0;
}
