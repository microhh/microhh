#include <cstdio>
#include <mpi.h>
#include "grid.h"
#include "defines.h"
#include "mpiinterface.h"

cmpi::cmpi(cgrid *gridin)
{
  std::printf("Creating instance of object mpi\n");
  grid = gridin;

  initialized = false;
}

cmpi::~cmpi()
{
  if(initialized)
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
  initialized = true;

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

  // transposez
  datacount = grid->imax*grid->jmax*grid->kblock;
  MPI_Type_contiguous(datacount, MPI_DOUBLE_PRECISION, &transposez);
  MPI_Type_commit(&transposez);
  std::printf("CvH: %d, %d, %d\n", datacount, datablock, datastride);

  // transposex
  datacount  = grid->jmax*grid->kblock;
  datablock  = grid->imax;
  datastride = grid->itot;
  MPI_Type_vector(datacount, datablock, datastride, MPI_DOUBLE_PRECISION, &transposex);
  MPI_Type_commit(&transposex);


  // create the requests arrays for the nonblocking sends
  reqsx = new MPI_Request[npx*2];
  reqsy = new MPI_Request[npy*2];

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

int cmpi::transposezx(double * restrict ar, double * restrict as)
{
  int startk;
  int nblock;
  int ncount = 1;

  int jj = grid->imax;
  int kk = grid->imax*grid->jmax;

  int kblock = grid->kblock;

  int reqidx = 0;

  for(int k=0; k<npx; k++)
  {
    // determine where to send it to
    nblock = mpiid - mpiid % npx + k;

    // determine where to fetch the send data
    int ijks = k*kblock*kk;

    // determine where to store the receive data
    int ijkr = k*jj;

    // send the block, tag it with the height (in kblocks) where it should come
    int sendtag = k;
    MPI_Isend(&as[ijks], ncount, transposez, nblock, sendtag, commxy, &reqsx[reqidx]);
    reqidx++;
    // and determine what has to be delivered at height k (in kblocks)
    int recvtag = mpiid % npx;
    MPI_Irecv(&ar[ijkr], ncount, transposex, nblock, recvtag, commxy, &reqsx[reqidx]);
    reqidx++;
    std::printf("MPI isend id %d, %d, send: %d, recv: %d\n", mpiid, nblock, sendtag, recvtag);
  }

  MPI_Waitall(reqidx, reqsx, MPI_STATUSES_IGNORE);

  return 0;
}

int cmpi::transposexz(double * restrict ar, double * restrict as)
{
  int starti, nblock;
  int ncount = 1;

  int jj = grid->imax;
  int kk = grid->imax*grid->jmax;

  for(int i=0; i<npx; i++)
  {
    // determine where to send it to
    nblock = mpiid - mpiid % npx + i;

    // determine what to send
    int ijks = i*jj;

    // determine what to receive 
    int ijkr = i*kk;

    MPI_Sendrecv(&as[ijks], ncount, transposex, nblock, 1,
                 &ar[ijkr], ncount, transposez, nblock, 1,
                 commxy, MPI_STATUS_IGNORE);
  }

  return 0;
}
