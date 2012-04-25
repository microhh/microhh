#include <cstdio>
#include <mpi.h>
#include "grid.h"
#include "defines.h"
#include "mpiinterface.h"

cmpi::cmpi(cgrid *gridin)
{
  std::printf("Creating instance of object mpi\n");
  grid = gridin;
}

cmpi::~cmpi()
{
  MPI_Finalize();
  std::printf("Destroying instance of object mpi\n");
}

int cmpi::readinifile(cinput *inputin)
{
  int n = 0;

  n += inputin->getItem(&npx, "grid", "npx", 1);
  n += inputin->getItem(&npy, "grid", "npy", 1);

  if(n > 0)
    return 1;
  
  return 0;
}

int cmpi::init()
{
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if(nprocs != npx*npy)
  {
    std::printf("ERROR nprocs = %d does not equal npx*npy = %d*%d\n", nprocs, npx, npy);
    return 1;
  }

  int dims    [2] = {npy, npx};
  int periodic[2] = {true, true};

  int n;

  if(MPI_Dims_create(nprocs, 2, dims))
    return 1;

  if(MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periodic, true, &commxy))
    return 1;
  if(MPI_Comm_rank(commxy, &mpiid))
    return 1;

  int mpicoords[2];
  n = MPI_Cart_coords(commxy, mpiid, 2, mpicoords);

  mpicoordx = mpicoords[1];
  mpicoordy = mpicoords[0];

  int dimx[2] = {false, true };
  int dimy[2] = {true , false};

  // MPI_Cart_sub(commxy, dimx, &commx);
  // MPI_Cart_sub(commxy, dimy, &commy);

  if(MPI_Cart_shift(commxy, 1, 1, &nwest , &neast ))
    return 1;
  if(MPI_Cart_shift(commxy, 0, 1, &nsouth, &nnorth))
    return 1;

  // create the MPI types for the cyclic boundary conditions
  int datacount, datablock, datastride;

  // east west
  datacount  = grid->jcells*grid->kcells;
  datablock  = grid->igc;
  datastride = grid->icells;
    
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE_PRECISION, &eastwestedge);
  MPI_Type_commit(&eastwestedge);

  // north south
  datacount  = grid->kcells;
  datablock  = grid->icells*grid->jgc;
  datastride = grid->icells*grid->jcells;
    
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE_PRECISION, &northsouthedge);
  MPI_Type_commit(&northsouthedge);

  return 0;
} 

int cmpi::boundary_cyclic(double * restrict data)
{
  int n;
  int ncount = 1;

  // communicate east-west edges
  int eastout = grid->iend-grid->igc;
  int westin  = 0;
  int westout = grid->istart;
  int eastin  = grid->iend;

  n = MPI_Sendrecv(&data[eastout], ncount, eastwestedge, neast, 1,
                   &data[westin ], ncount, eastwestedge, nwest, 1,
                   commxy, MPI_STATUS_IGNORE);

  n = MPI_Sendrecv(&data[westout], ncount, eastwestedge, nwest, 2,
                   &data[eastin ], ncount, eastwestedge, neast, 2,
                   commxy, MPI_STATUS_IGNORE);

  // communicate north-south edges
  int northout = (grid->jend-grid->jgc)*grid->icells;
  int southin  = 0;
  int southout = grid->jstart*grid->icells;
  int northin  = grid->jend  *grid->icells;

  n = MPI_Sendrecv(&data[northout], ncount, northsouthedge, nnorth, 1,
                   &data[southin ], ncount, northsouthedge, nsouth, 1,
                   commxy, MPI_STATUS_IGNORE);

  n = MPI_Sendrecv(&data[southout], ncount, northsouthedge, nsouth, 2,
                   &data[northin ], ncount, northsouthedge, nnorth, 2,
                   commxy, MPI_STATUS_IGNORE);

  return 0;
}
