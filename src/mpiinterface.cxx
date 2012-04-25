#include <cstdio>
#include <mpi.h>
#include "grid.h"
#include "defines.h"
#include "mpiinterface.h"

cmpi::cmpi()
{
  std::printf("Creating instance of object mpi\n");
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

int cmpi::init(cgrid *gridin)
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
  datacount  = gridin->jcells*gridin->kcells;
  datablock  = gridin->igc;
  datastride = gridin->icells;
    
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE_PRECISION, &eastwestedge);
  MPI_Type_commit(&eastwestedge);

  // north south
  datacount  = gridin->kcells;
  datablock  = gridin->icells*gridin->jgc;
  datastride = gridin->icells*gridin->jcells;
    
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE_PRECISION, &northsouthedge);
  MPI_Type_commit(&northsouthedge);

  return 0;
} 

int cmpi::boundary_cyclic(double * restrict data, cgrid *gridin)
{
  int ncount = 1;
  MPI_Status status;

  int eastout = gridin->iend - gridin->igc;
  int westin  = 0;
  int westout = gridin->istart;
  int eastin  = gridin->iend;

  // communicate east-west edges
  MPI_Sendrecv(&data[eastout], ncount, eastwestedge, neast, 1,
               &data[westin ], ncount, eastwestedge, nwest, 1,
               commxy, &status);

  MPI_Sendrecv(&data[westout], ncount, eastwestedge, nwest, 2,
               &data[eastin ], ncount, eastwestedge, neast, 2,
               commxy, &status);
/*
  // communicate north-south edges
  call MPI_SENDRECV(var(1-kgc, 1-igc, 1     ), ncount, northsouthedge, nsouth, 3, &
                    var(1-kgc, 1-igc, jmax+1), ncount, northsouthedge, nnorth, 3, &
                    comm2d, mpistatus, mpierr)

  call MPI_SENDRECV(var(1-kgc, 1-igc, jmax-jgc+1), ncount, northsouthedge, nnorth, 4, &
                    var(1-kgc, 1-igc, 1-jgc     ), ncount, northsouthedge, nsouth, 4, &
                    comm2d, mpistatus, mpierr)*/

  return 0;
}
